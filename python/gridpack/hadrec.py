# -------------------------------------------------------------
# file: gridpack/hadrec.py
# -------------------------------------------------------------
# Re-exports the compiled pybind11 hadrec submodule and adds the
# high-level :class:`Hadrec` driver on top.
#
# The :class:`Hadrec` class exposes the HADRECAppModule surface
# (PF-then-DS chain, direct Action objects, RTPR-style step-by-step
# observation queries).  If you just want a full run, use
# :class:`gridpack.DynamicSim`; if you want a Gym-style step API that
# hides the underlying app, use :class:`gridpack.DynamicSimStepper`.
# Both of those already exist in :mod:`gridpack.dynamic_sim` and are
# implemented in terms of HADRECAppModule.
#
# See also :mod:`gridpack.powerflow` for the wrapper pattern
# (register with Session, WeakSet of live results, null-out on close).
# -------------------------------------------------------------

from __future__ import annotations

import weakref
from typing import List, Optional, Sequence, Tuple

from ._gridpack.hadrec import *  # noqa: F401,F403
from ._gridpack import hadrec as _mod  # noqa: F401
from ._gridpack.dynamic_simulation import Event, EventVector

from .session import Session
from .results import DSFResult


# Action-type IDs recognized by HADRECAppModule::applyAction.
# 0 = load shedding, 1 = line tripping, 2 = generator tripping.
ACT_LOAD_SHED = 0
ACT_LINE_TRIP = 1
ACT_GEN_TRIP = 2


def __getattr__(name):
    # Fallback for any pybind11 attribute not explicitly re-exported.
    return getattr(_mod, name)


def make_load_shed_action(bus_number: int, percentage: float,
                          load_id: str = "1"):
    """Return a :class:`Action` configured for load shedding.

    ``percentage`` is signed: ``-0.2`` sheds 20% of P.
    """
    from . import _gridpack as _ext
    act = _ext.hadrec.Action()
    act.actiontype = ACT_LOAD_SHED
    act.bus_number = int(bus_number)
    act.componentID = str(load_id)
    act.percentage = float(percentage)
    return act


def make_line_trip_action(from_bus: int, to_bus: int, ckt: str = "1 "):
    """Return a :class:`Action` configured for line tripping."""
    from . import _gridpack as _ext
    act = _ext.hadrec.Action()
    act.actiontype = ACT_LINE_TRIP
    act.brch_from_bus_number = int(from_bus)
    act.brch_to_bus_number = int(to_bus)
    act.branch_ckt = str(ckt)
    return act


def make_generator_trip_action(bus_number: int, gen_id: str = "1"):
    """Return a :class:`Action` configured for generator tripping."""
    from . import _gridpack as _ext
    act = _ext.hadrec.Action()
    act.actiontype = ACT_GEN_TRIP
    act.bus_number = int(bus_number)
    act.componentID = str(gen_id)
    return act


class Hadrec:
    """High-level HADREC driver.

    Wraps :class:`gridpack.hadrec.Module` for direct HADREC / RTPR use.
    Runs the power-flow-then-dynamic-simulation chain, optionally
    accepts fault events, and exposes ``executeDynSimuOneStep`` +
    ``getObservations`` for tight step loops.

    Usage::

        from gridpack import Session, Hadrec
        from gridpack.hadrec import make_line_trip_action

        with Session() as s:
            h = Hadrec(s, "input.xml")
            h.initialize_dyn_simu()               # optional faults=[...]
            trip = make_line_trip_action(6, 7)

            while not h.done:
                if h.step_count == 400:
                    h.apply_action(trip)
                h.step()
            # per-step time series is in h.result

    Parameters
    ----------
    session : gridpack.Session
    input_file : str
        GridPACK dynamic-simulation XML.
    pf_idx : int, optional
        Power-flow configuration index (default ``-1``).
    suppress_output : bool, optional
        Enable ``NoPrint`` for the run.
    record_bus_freq : bool, optional
        Use ``getObservations_withBusFreq`` and
        ``getObservationLists_withBusFreq`` (requires bus-frequency
        entries in the XML ``<observations>`` block).
    """

    def __init__(
        self,
        session: Session,
        input_file: str,
        *,
        pf_idx: int = -1,
        suppress_output: bool = False,
        record_bus_freq: bool = False,
    ) -> None:
        if session.closed:
            raise RuntimeError("Session is closed")

        from . import _gridpack as _ext

        self._session = session
        self.input_file = input_file
        self._pf_idx = pf_idx
        self._record_bus_freq = record_bus_freq

        if suppress_output:
            _ext.NoPrint().setStatus(True)
        self._suppress_output = bool(suppress_output)

        # HADREC solves PF as part of its own initialization; it opens
        # its own Configuration internally from the file path.
        self._hadapp = _ext.hadrec.Module()
        self._hadapp.solvePowerFlowBeforeDynSimu(input_file, pf_idx)
        self._pf_solved = True

        self._initialized = False
        self._step_count = 0
        self._closed = False
        self._result: Optional[DSFResult] = None
        self._live_results: "weakref.WeakSet[DSFResult]" = weakref.WeakSet()

        session.register(self)

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    @property
    def hadapp(self):
        """The underlying :class:`gridpack.hadrec.Module` (advanced use)."""
        self._require_open()
        return self._hadapp

    @property
    def step_count(self) -> int:
        return self._step_count

    @property
    def done(self) -> bool:
        """True once the dynamic simulation has completed."""
        self._require_open()
        return bool(self._hadapp.isDynSimuDone())

    @property
    def result(self) -> Optional[DSFResult]:
        """The :class:`DSFResult` accumulating this run's observations.

        ``None`` until :meth:`initialize_dyn_simu` is called.
        """
        return self._result

    # ------------------------------------------------------------------
    # PF-to-DS transfer + initialization
    # ------------------------------------------------------------------

    def transfer_pf_to_ds(self) -> None:
        """Copy the post-PF state into the DS initial condition.

        Called automatically by :meth:`initialize_dyn_simu` on first
        entry; expose it here so RTPR-style loops that re-run PF between
        events can also invoke it.
        """
        self._require_open()
        self._hadapp.transferPFtoDS()

    def initialize_dyn_simu(
        self,
        faults: Optional[Sequence[Event]] = None,
        *,
        dscase_idx: int = -1,
    ) -> DSFResult:
        """Initialize the dynamic simulation and return a fresh
        :class:`DSFResult` that accumulates per-step observations.

        Parameters
        ----------
        faults : Sequence[Event] or None
            Fault events to seed the simulation with.  If omitted, an
            empty :class:`EventVector` is passed.
        dscase_idx : int
            Dynamic-simulation configuration index (default ``-1``).
        """
        self._require_open()
        if not self._pf_solved:
            raise RuntimeError(
                "Hadrec: solvePowerFlowBeforeDynSimu has not been called"
            )
        # Do the PF -> DS transfer here (idempotent from HADREC's side).
        self._hadapp.transferPFtoDS()

        from . import _gridpack as _ext
        event_vec = _ext.dynamic_simulation.EventVector()
        if faults:
            for f in faults:
                event_vec.append(f)
        self._hadapp.initializeDynSimu(event_vec, dscase_idx)
        self._initialized = True
        self._step_count = 0

        # Fresh DSFResult for this initialization pass.  Populate
        # observation-list metadata up front so downstream tooling
        # (channel names, dataframes) works.
        result = DSFResult(
            self._hadapp,
            input_file=self.input_file,
            with_bus_freq=self._record_bus_freq,
        )
        try:
            if self._record_bus_freq:
                (
                    result.obs_gen_buses,
                    result.obs_gen_ids,
                    result.obs_load_buses,
                    result.obs_load_ids,
                    result.obs_bus_ids,
                    result.obs_bus_freq_ids,
                ) = self._hadapp.getObservationLists_withBusFreq()
            else:
                (
                    result.obs_gen_buses,
                    result.obs_gen_ids,
                    result.obs_load_buses,
                    result.obs_load_ids,
                    result.obs_bus_ids,
                ) = self._hadapp.getObservationLists()
        except Exception:
            pass

        self._result = result
        self._live_results.add(result)
        return result

    # ------------------------------------------------------------------
    # Stepping
    # ------------------------------------------------------------------

    def step(self, *, record: bool = True):
        """Advance the simulation one time step.

        Returns the observation tuple for this step (or ``None`` if the
        run is already done).  If ``record=True`` (default), the
        observation is appended to :attr:`result`.
        """
        self._require_open()
        if not self._initialized:
            raise RuntimeError(
                "Hadrec.step called before initialize_dyn_simu"
            )
        if self._hadapp.isDynSimuDone():
            return None
        self._hadapp.executeDynSimuOneStep()
        self._step_count += 1
        if self._record_bus_freq:
            obs = self._hadapp.getObservations_withBusFreq() \
                if hasattr(self._hadapp, "getObservations_withBusFreq") \
                else self._hadapp.getObservations()
        else:
            obs = self._hadapp.getObservations()
        if record and self._result is not None:
            self._result.times.append(float(self._step_count))
            self._result.observations.append(obs)
        return obs

    def run_until_done(self, *, record: bool = True) -> DSFResult:
        """Convenience: step until :attr:`done` is True."""
        self._require_open()
        while not self.done:
            self.step(record=record)
        return self._result

    # ------------------------------------------------------------------
    # Actions / signals
    # ------------------------------------------------------------------

    def apply_action(self, action) -> None:
        """Apply a pre-built :class:`Action` (see :func:`make_*_action`)."""
        self._require_open()
        self._hadapp.applyAction(action)

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
    # Observation and query passthroughs
    # ------------------------------------------------------------------

    def get_observations(self):
        """Return the current-step observation tuple.

        Handy for callers that want to sample without stepping (e.g.
        immediately after :meth:`initialize_dyn_simu`).
        """
        self._require_open()
        return self._hadapp.getObservations()

    def get_observation_lists(self):
        """Return ``(gen_buses, gen_ids, load_buses, load_ids, bus_ids)``.

        With ``record_bus_freq=True`` a 6-tuple is returned instead,
        with an appended ``bus_freq_ids`` list.
        """
        self._require_open()
        if self._record_bus_freq:
            return self._hadapp.getObservationLists_withBusFreq()
        return self._hadapp.getObservationLists()

    def get_generator_power(self, bus_id: int, gen_id: str = "1"):
        """Return ``(P, Q)`` for a generator, or ``None`` if absent."""
        self._require_open()
        return self._hadapp.getGeneratorPower(bus_id, gen_id)

    def get_bus_total_load_power(self, bus_id: int):
        """Return ``(P, Q)`` for the aggregate load at a bus, or ``None``."""
        self._require_open()
        return self._hadapp.getBusTotalLoadPower(bus_id)

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def close(self) -> None:
        """Release the underlying pybind11 module.

        Idempotent; called automatically at Session close.
        """
        if self._closed:
            return
        for r in list(self._live_results):
            try:
                r.close()
            except Exception:
                pass
        self._live_results.clear()
        self._hadapp = None
        self._result = None
        self._closed = True

    def _require_open(self) -> None:
        if self._closed:
            raise RuntimeError("Hadrec is closed")

    def __enter__(self) -> "Hadrec":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        self.close()
        return False


__all__ = [
    "Hadrec",
    "make_load_shed_action",
    "make_line_trip_action",
    "make_generator_trip_action",
    "ACT_LOAD_SHED",
    "ACT_LINE_TRIP",
    "ACT_GEN_TRIP",
]
