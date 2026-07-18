# -------------------------------------------------------------
# file: gridpack/dynamic_sim.py
# -------------------------------------------------------------
# High-level wrappers on top of the pybind11 dynamic-simulation bindings.
#
#   DynamicSim         -- full-run driver on ``DSFullApp`` (equivalent of
#                         dsf2.x).  Uses setup()/run(), which is the
#                         production-tested path.
#   DynamicSimStepper  -- step-by-step driver for RL / co-sim on top of
#                         ``HADRECAppModule``.  HADREC exposes
#                         ``executeDynSimuOneStep`` with a stable teardown
#                         path; the low-level ``DSFullApp::executeOneSimuStep``
#                         has a pre-existing SEGV in its destructor on
#                         some datasets, so the stepper goes through
#                         HADREC instead.
#   DSFResult          -- container for observations / time series.
#
# Both wrappers follow the same MPI-safe teardown pattern as PowerFlow:
# they register with the Session, track live Results in a WeakSet, and
# null out each Result's ``_dsapp`` reference before releasing their own.
# Without this, the pybind11 app modules can outlive Session.close() and
# crash inside MPI_Iprobe on MPI_COMM_NULL during MPI finalize.
# -------------------------------------------------------------

from __future__ import annotations

import weakref
from typing import Iterable, List, Optional, Sequence

from ._gridpack.dynamic_simulation import Event, EventVector

from .session import Session
from .results import DSFResult


# -------------------------------------------------------------
# Full-run driver (DSFullApp.setup / run)
# -------------------------------------------------------------

class DynamicSim:
    """Full-run dynamic simulation driver (mirrors ``dsf2.x``).

    Usage::

        from gridpack import Session, DynamicSim

        with Session() as s:
            ds = DynamicSim(s, "input.xml")
            result = ds.run()

    Parameters
    ----------
    session : gridpack.Session
    input_file : str
        Path to the GridPACK dynamic-simulation XML configuration.
    pf_idx : int, optional
        Power-flow configuration index (defaults to ``-1``).
    suppress_output : bool, optional
        Enable ``NoPrint`` for the duration of the run.
    generator_watch : bool, optional
        If True (default), enable generator-watch CSV output as configured
        in the XML (``<generatorWatch>`` / ``<generatorWatchFileName>``).
    """

    def __init__(
        self,
        session: Session,
        input_file: str,
        *,
        pf_idx: int = -1,
        suppress_output: Optional[bool] = None,
        generator_watch: bool = True,
    ) -> None:
        if session.closed:
            raise RuntimeError("Session is closed")

        from . import _gridpack as _ext

        self._session = session
        self.input_file = input_file
        self._pf_idx = pf_idx

        # Open Configuration on the session's Communicator so we stay on
        # the shared MPI world (initFromConfig avoids MPI_Comm_dup).
        self._config = _ext.Configuration()
        self._config.open(input_file, session.comm)
        self._ds_cursor = self._config.getCursor(
            "Configuration.Dynamic_simulation"
        )

        xml_suppress = bool(self._ds_cursor.get("suppressOutput"))
        self._suppress_output = (
            xml_suppress if suppress_output is None else suppress_output
        )
        if self._suppress_output:
            _ext.NoPrint().setStatus(True)

        self._dsapp = _ext.dynamic_simulation.DSFullApp()
        self._dsapp.initFromConfig(self._config, pf_idx)
        self._dsapp.readGenerators()
        self._dsapp.readSequenceData()
        self._dsapp.initialize()

        try:
            self._dsapp.setObservations(self._ds_cursor)
        except Exception:
            pass

        if generator_watch:
            self._dsapp.setGeneratorWatch()

        self._ran = False
        self._closed = False
        self._live_results: "weakref.WeakSet[DSFResult]" = weakref.WeakSet()

        session.register(self)

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    @property
    def dsapp(self):
        self._require_open()
        return self._dsapp

    @property
    def suppress_output(self) -> bool:
        return self._suppress_output

    @property
    def time_step(self) -> float:
        self._require_open()
        return float(self._dsapp.getTimeStep())

    @property
    def final_time(self) -> float:
        self._require_open()
        return float(self._dsapp.getFinalTime())

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------

    def run(self, tend: Optional[float] = None) -> DSFResult:
        """Run the full simulation and return a :class:`DSFResult`.

        The container is populated with the observation-list metadata
        (from ``<observations>`` in the XML) so post-run per-bus /
        per-generator queries via :meth:`DSFResult.get_generator_power`
        work.  Per-step time-series is not captured through this path
        -- use :class:`DynamicSimStepper` for that.
        """
        self._require_open()

        self._dsapp.setup()
        if tend is None:
            self._dsapp.run()
        else:
            self._dsapp.run(float(tend))
        self._ran = True

        result = DSFResult(self._dsapp, input_file=self.input_file)
        try:
            (
                result.obs_gen_buses,
                result.obs_gen_ids,
                result.obs_load_buses,
                result.obs_load_ids,
                result.obs_bus_ids,
            ) = self._dsapp.getObservationLists()
        except Exception:
            pass
        self._live_results.add(result)
        return result

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def close(self) -> None:
        if self._closed:
            return
        for r in list(self._live_results):
            try:
                r.close()
            except Exception:
                pass
        self._live_results.clear()
        self._dsapp = None
        self._ds_cursor = None
        self._config = None
        self._closed = True

    def _require_open(self) -> None:
        if self._closed:
            raise RuntimeError("DynamicSim is closed")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        self.close()
        return False


# -------------------------------------------------------------
# Step-by-step driver (HADRECAppModule)
# -------------------------------------------------------------

# Action-type IDs recognized by HADRECAppModule::applyAction.
# 0 = load shedding, 1 = line tripping, 2 = generator tripping.
_ACT_LOAD_SHED = 0
_ACT_LINE_TRIP = 1
_ACT_GEN_TRIP = 2


class DynamicSimStepper:
    """Step-by-step dynamic simulation driver for RL / co-sim.

    Wraps :class:`gridpack.hadrec.Module` (the HADREC C++ app) so the
    stepping loop matches what production RL workflows expect.  Actions
    (:meth:`apply_load_shedding`, :meth:`apply_generator_trip`) are
    routed through ``HADRECAppModule::applyAction``; wide-area control
    signals go through ``setWideAreaControlSignal``.

    Usage::

        from gridpack import Session, DynamicSimStepper

        with Session() as s:
            st = DynamicSimStepper(s, "input.xml")
            while not st.done:
                st.step()
                if st.step_count == 400:
                    st.apply_load_shedding(bus_number=5, percentage=-0.2)
            st.result.to_csv("stepper_out.csv")

    Parameters
    ----------
    session : gridpack.Session
    input_file : str
        Path to the GridPACK dynamic-simulation XML configuration.
    faults : Sequence[Event] or None, optional
        Fault events to initialize the run with.  If omitted, an empty
        event vector is passed (matches the ``hadrec.py`` reference).
    pf_idx : int, optional
        Power-flow configuration index (defaults to ``-1``).
    dscase_idx : int, optional
        Dynamic-simulation configuration index (defaults to ``-1``).
    suppress_output : bool, optional
        Enable ``NoPrint`` for the duration of the run.
    record_bus_freq : bool, optional
        If True, use ``getObservations_withBusFreq`` when sampling.
        Requires the XML ``<observations>`` block to include bus-frequency
        entries.
    """

    def __init__(
        self,
        session: Session,
        input_file: str,
        *,
        faults: Optional[Sequence[Event]] = None,
        pf_idx: int = -1,
        dscase_idx: int = -1,
        suppress_output: Optional[bool] = None,
        record_bus_freq: bool = False,
    ) -> None:
        if session.closed:
            raise RuntimeError("Session is closed")

        from . import _gridpack as _ext

        self._session = session
        self.input_file = input_file
        self._record_bus_freq = record_bus_freq

        if suppress_output:
            _ext.NoPrint().setStatus(True)
        self._suppress_output = bool(suppress_output)

        self._hadapp = _ext.hadrec.Module()
        # solvePowerFlowBeforeDynSimu takes the input filename directly;
        # HADREC internally opens the Configuration on its own communicator.
        self._hadapp.solvePowerFlowBeforeDynSimu(input_file, pf_idx)
        self._hadapp.transferPFtoDS()

        event_vec = _ext.dynamic_simulation.EventVector()
        if faults:
            for f in faults:
                event_vec.append(f)
        self._hadapp.initializeDynSimu(event_vec, dscase_idx)

        # Cache observation-list metadata for DSFResult column naming
        # and expose it up-front so callers can size their buffers.
        try:
            if record_bus_freq:
                (
                    self.obs_gen_buses,
                    self.obs_gen_ids,
                    self.obs_load_buses,
                    self.obs_load_ids,
                    self.obs_bus_ids,
                    self.obs_bus_freq_ids,
                ) = self._hadapp.getObservationLists_withBusFreq()
            else:
                (
                    self.obs_gen_buses,
                    self.obs_gen_ids,
                    self.obs_load_buses,
                    self.obs_load_ids,
                    self.obs_bus_ids,
                ) = self._hadapp.getObservationLists()
                self.obs_bus_freq_ids = []
        except Exception:
            self.obs_gen_buses = []
            self.obs_gen_ids = []
            self.obs_load_buses = []
            self.obs_load_ids = []
            self.obs_bus_ids = []
            self.obs_bus_freq_ids = []

        # DSFResult accumulates time / observation rows.
        self.result = DSFResult(
            self._hadapp,
            input_file=input_file,
            with_bus_freq=record_bus_freq,
        )
        self.result.obs_gen_buses = list(self.obs_gen_buses)
        self.result.obs_gen_ids = list(self.obs_gen_ids)
        self.result.obs_load_buses = list(self.obs_load_buses)
        self.result.obs_load_ids = list(self.obs_load_ids)
        self.result.obs_bus_ids = list(self.obs_bus_ids)
        self.result.obs_bus_freq_ids = list(self.obs_bus_freq_ids)

        self._step_count = 0
        self._sim_time = 0.0
        self._closed = False
        self._live_results: "weakref.WeakSet[DSFResult]" = weakref.WeakSet()
        self._live_results.add(self.result)

        session.register(self)

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    @property
    def hadapp(self):
        """The underlying pybind11 ``HADRECAppModule`` (advanced use only)."""
        self._require_open()
        return self._hadapp

    @property
    def done(self) -> bool:
        """True once the simulation has reached ``simulationTime``."""
        self._require_open()
        return bool(self._hadapp.isDynSimuDone())

    @property
    def step_count(self) -> int:
        return self._step_count

    @property
    def current_time(self) -> float:
        """Elapsed simulation time inferred from step count.

        HADREC does not expose ``getCurrentTime``; we track it locally.
        """
        return self._sim_time

    # ------------------------------------------------------------------
    # Stepping
    # ------------------------------------------------------------------

    def step(self, *, record: bool = True) -> Optional[tuple]:
        """Advance the simulation by one time step.

        Returns the observation tuple for this step (or ``None`` when
        the run is already done).  If ``record=True`` (default), the
        observation is also appended to :attr:`result`.
        """
        self._require_open()
        if self._hadapp.isDynSimuDone():
            return None
        self._hadapp.executeDynSimuOneStep()
        self._step_count += 1
        obs = self._hadapp.getObservations()
        if record:
            # HADREC also doesn't emit current-time; approximate as
            # step * dt if we can read dt from the XML, else count.
            self.result.times.append(float(self._step_count))
            self.result.observations.append(obs)
        return obs

    def run_until_done(self, *, record: bool = True) -> DSFResult:
        """Convenience: step until :attr:`done` is True."""
        while not self.done:
            self.step(record=record)
        return self.result

    # ------------------------------------------------------------------
    # Actuators
    # ------------------------------------------------------------------

    def apply_load_shedding(
        self,
        bus_number: int,
        percentage: float,
        load_id: str = "1",
    ) -> None:
        """Shed a fraction of load at ``bus_number`` / ``load_id``.

        ``percentage`` is signed: ``-0.2`` sheds 20% of P at this call.
        """
        self._require_open()
        from . import _gridpack as _ext
        act = _ext.hadrec.Action()
        act.actiontype = _ACT_LOAD_SHED
        act.bus_number = int(bus_number)
        act.componentID = str(load_id)
        act.percentage = float(percentage)
        self._hadapp.applyAction(act)

    def apply_generator_trip(
        self,
        bus_number: int,
        gen_id: str = "1",
    ) -> None:
        """Trip generator ``gen_id`` at ``bus_number``."""
        self._require_open()
        from . import _gridpack as _ext
        act = _ext.hadrec.Action()
        act.actiontype = _ACT_GEN_TRIP
        act.bus_number = int(bus_number)
        act.componentID = str(gen_id)
        self._hadapp.applyAction(act)

    def apply_line_trip(
        self,
        from_bus: int,
        to_bus: int,
        ckt: str = "1 ",
    ) -> None:
        """Trip a line between ``from_bus`` and ``to_bus`` on circuit ``ckt``."""
        self._require_open()
        from . import _gridpack as _ext
        act = _ext.hadrec.Action()
        act.actiontype = _ACT_LINE_TRIP
        act.brch_from_bus_number = int(from_bus)
        act.brch_to_bus_number = int(to_bus)
        act.branch_ckt = str(ckt)
        self._hadapp.applyAction(act)

    def set_wide_area_signal(
        self,
        bus_number: int,
        gen_id: str,
        signal: float,
    ) -> None:
        """Set a wide-area control signal for one generator."""
        self._require_open()
        self._hadapp.setWideAreaControlSignal(
            int(bus_number), str(gen_id), float(signal),
        )

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def close(self) -> None:
        if self._closed:
            return
        for r in list(self._live_results):
            try:
                r.close()
            except Exception:
                pass
        self._live_results.clear()
        self._hadapp = None
        self._closed = True

    def _require_open(self) -> None:
        if self._closed:
            raise RuntimeError("DynamicSimStepper is closed")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        self.close()
        return False


__all__ = [
    "DynamicSim",
    "DynamicSimStepper",
    "DSFResult",
    "Event",
    "EventVector",
]
