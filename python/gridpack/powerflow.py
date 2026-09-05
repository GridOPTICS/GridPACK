# -------------------------------------------------------------
# file: gridpack/powerflow.py
# -------------------------------------------------------------
# Re-exports the compiled pybind11 powerflow submodule (so existing
# code like `from gridpack.powerflow import Powerflow` keeps working)
# and adds the high-level PowerFlow driver class on top.
#
# See also :mod:`gridpack.results` for the returned PowerFlowResult.
# -------------------------------------------------------------

from __future__ import annotations

from typing import Optional

# Re-export the pybind11 submodule for backwards compatibility.
from ._gridpack.powerflow import *  # noqa: F401,F403
from ._gridpack import powerflow as _pfmod  # noqa: F401

from .session import Session
from .results import PowerFlowResult
from .exceptions import PowerFlowDiverged


def __getattr__(name):
    # Fallback for any pybind11 attribute not explicitly re-exported.
    return getattr(_pfmod, name)


def _xml_bool(raw, default: bool = False) -> bool:
    """Parse an XML scalar as a bool.  cursor.get() returns a string, so
    bool(raw) would be True for "false"."""
    if raw is None:
        return default
    if isinstance(raw, bool):
        return raw
    text = str(raw).strip().lower()
    if text in ("true", "yes", "on", "1"):
        return True
    if text in ("false", "no", "off", "0", ""):
        return False
    return default


def _xml_number(raw, cast, default=None):
    """Parse an XML scalar as ``cast``, or ``default`` if unusable."""
    if raw is None:
        return default
    try:
        return cast(str(raw).strip())
    except (TypeError, ValueError):
        return default


class PowerFlow:
    """High-level power flow driver.

    Usage::

        from gridpack import Session, PowerFlow

        with Session() as s:
            pf = PowerFlow(s, "input_14.xml")
            result = pf.solve()
            result.write()                    # dump to stdout
            df = result.to_dataframe()        # all buses, all columns
            bad = result.violations()         # voltage / overload report

    Parameters
    ----------
    session : gridpack.Session
        Active session owning the MPI environment.
    input_file : str
        Path to the GridPACK powerflow XML configuration.
    idx : int, optional
        Index of the power flow problem within a multi-network XML.
        Defaults to -1 (single-network XML).
    suppress_output : bool, optional
        If True, calls ``pfapp.suppressOutput(True)`` before solving.
        Reads ``Configuration.Powerflow.suppressOutput`` from XML if not
        overridden here.
    comm : gridpack.Communicator, optional
        Build the network on this communicator instead of the session's.
        ContingencyAnalysis passes a task communicator so each task owns a
        private copy; results are then task-local and never gathered.
    """

    def __init__(
        self,
        session: Session,
        input_file: str,
        *,
        idx: int = -1,
        suppress_output: Optional[bool] = None,
        comm: Optional[object] = None,
    ) -> None:
        if session.closed:
            raise RuntimeError("Session is closed")

        from . import _gridpack as _ext

        self._session = session
        self.input_file = input_file
        self._idx = idx

        # Read configuration.
        self._config = _ext.Configuration()
        self._config.open(input_file, session.comm)
        cursor = self._config.getCursor("Configuration.Powerflow")
        self._cursor = cursor

        # Precompute XML-controlled options so solve() can respect them.
        self._xml_nonlinear = _xml_bool(cursor.get("UseNonLinear"))
        xml_suppress = _xml_bool(cursor.get("suppressOutput"))
        self._suppress_output = (
            xml_suppress if suppress_output is None else suppress_output
        )

        # Reported in PowerFlowDiverged.
        self._tolerance = _xml_number(cursor.get("tolerance"), float)
        self._max_iteration = _xml_number(cursor.get("maxIteration"), int)

        # PSS/E exports requested by the XML, applied by the caller.
        self._psse_exports = [
            (v, cursor.get("exportPSSE_v%d" % v))
            for v in (23, 33, 34)
        ]
        self._psse_exports = [(v, f) for v, f in self._psse_exports if f]

        # Create the pybind11 Powerflow application.
        self._pfapp = _ext.powerflow.Powerflow()
        if self._suppress_output:
            self._pfapp.suppressOutput(True)

        # A network on a sub-communicator is private to this task, so
        # gathering its results over the session's world comm would splice
        # together unrelated copies.
        self._task_local = comm is not None
        if comm is None:
            self._pfapp.readNetwork(self._config, idx)
        else:
            self._pfapp.readNetwork(self._config, idx, comm)
        self._pfapp.initialize()

        self._output_open = False

        self._solved = False
        self._closed = False

        # Track live results so we can null out their pfapp reference
        # when we close (otherwise the pybind11 PFAppModule outlives
        # this wrapper and crashes at MPI finalize).
        import weakref
        self._live_results: "weakref.WeakSet[PowerFlowResult]" = weakref.WeakSet()

        # Ensure the session tears us down before it releases the
        # Communicator/Environment.
        session.register(self)

    # ------------------------------------------------------------------
    # Configuration accessors (read from XML at construction time)
    # ------------------------------------------------------------------

    @property
    def use_nonlinear_from_xml(self) -> bool:
        """Value of ``UseNonLinear`` in the XML config."""
        return self._xml_nonlinear

    @property
    def suppress_output(self) -> bool:
        return self._suppress_output

    @property
    def network_is_task_local(self) -> bool:
        """True if the network was built on a caller-supplied communicator."""
        return self._task_local

    @property
    def psse_exports_from_xml(self):
        """``[(version, path)]`` from the ``exportPSSE_v*`` XML keys."""
        return list(self._psse_exports)

    # ------------------------------------------------------------------
    # Solve
    # ------------------------------------------------------------------

    def solve(
        self,
        nonlinear: Optional[bool] = None,
        *,
        strict: bool = True,
    ) -> PowerFlowResult:
        """Run the solver and return a :class:`PowerFlowResult`.

        ``nonlinear`` overrides the XML ``UseNonLinear`` key: True picks
        ``nl_solve``, False the Newton-Raphson loop.

        Raises :class:`gridpack.PowerFlowDiverged` if the tolerance is not
        reached; every rank raises, so it is safe to catch under MPI.  Pass
        ``strict=False`` to get the unconverged result instead --
        ``result.mismatch`` names the worst bus.

        ``nl_solve`` populates no convergence summary, so on that path the
        verdict is its return value and ``result.convergence`` is None.
        """
        use_nl = self._xml_nonlinear if nonlinear is None else nonlinear
        if use_nl:
            ok = self._pfapp.nl_solve()
        else:
            ok = self._pfapp.solve()
        self._pfapp.saveData()
        self._solved = True

        result = PowerFlowResult(
            self._pfapp,
            nonlinear=use_nl,
            input_file=self.input_file,
            solver_converged=bool(ok),
            mpi_comm=None if self._task_local else self._session.mpi_comm,
        )
        self._live_results.add(result)

        if strict and result.converged is False:
            raise PowerFlowDiverged(
                convergence=result.convergence,
                nonlinear=use_nl,
                tolerance=self._tolerance,
                max_iteration=self._max_iteration,
                input_file=self.input_file,
            )
        return result

    # ------------------------------------------------------------------
    # Contingency primitives
    # ------------------------------------------------------------------
    # Thin wrappers over PFAppModule, driven by ContingencyAnalysis; also
    # the supported way to script an ad-hoc contingency.

    def set_voltage_limits(self, vmin: float, vmax: float) -> None:
        """Set the voltage band used by :meth:`check_voltage_violations`."""
        self._require_open()
        self._pfapp.setVoltageLimits(float(vmin), float(vmax))

    def set_contingency(self, contingency) -> bool:
        """Trip the lines or generators named by ``contingency``.

        False if no element matched -- the next solve would just repeat
        the base case.
        """
        self._require_open()
        return bool(self._pfapp.setContingency(contingency))

    def unset_contingency(self, contingency) -> None:
        """Restore the elements tripped by :meth:`set_contingency`."""
        self._require_open()
        self._pfapp.unSetContingency(contingency)

    def reset_voltages(self) -> None:
        """Restore bus voltages to their initial values.

        Required between contingencies, or convergence depends on order.
        """
        self._require_open()
        self._pfapp.resetVoltages()

    def ignore_voltage_violations(self) -> None:
        """Exempt currently violating buses from later voltage checks.

        Called on the base case so pre-existing violations are not
        re-reported against every contingency.
        """
        self._require_open()
        self._pfapp.ignoreVoltageViolations()

    def check_qlim_violations(self) -> bool:
        """True if no generator is outside its reactive limits.

        False means a bus was converted PV -> PQ, so re-solve.
        """
        self._require_open()
        return bool(self._pfapp.checkQlimViolations())

    def clear_qlim_violations(self) -> None:
        """Undo the PV -> PQ conversions from :meth:`check_qlim_violations`."""
        self._require_open()
        self._pfapp.clearQlimViolations()

    def check_voltage_violations(self) -> bool:
        """True if every bus is inside the band, ignoring exempted buses."""
        self._require_open()
        return bool(self._pfapp.checkVoltageViolations())

    def check_line_overload_violations(self) -> bool:
        """True if no branch exceeds its rating.

        Tests one end per branch, unlike
        :meth:`PowerFlowResult.violations`, which uses ``max(from, to)``.
        """
        self._require_open()
        return bool(self._pfapp.checkLineOverloadViolations())

    # ------------------------------------------------------------------
    # Output redirection
    # ------------------------------------------------------------------

    def open_output(self, path: str) -> None:
        """Redirect this application's output to ``path``."""
        self._require_open()
        self._pfapp.open(path)
        self._output_open = True

    def close_output(self) -> None:
        """Stop redirecting output.  Idempotent."""
        if self._closed or not self._output_open:
            return
        self._pfapp.close()
        self._output_open = False

    def print_output(self, text: str) -> None:
        """Write ``text`` through the application's output stream."""
        self._require_open()
        self._pfapp.print(text)

    def _require_open(self) -> None:
        if self._closed:
            raise RuntimeError("PowerFlow is closed")

    def close(self) -> None:
        """Release the underlying pybind11 objects in the correct order.

        Idempotent; called automatically at Session close.
        """
        if self._closed:
            return
        # Invalidate any outstanding Result objects so their strong
        # reference to _pfapp is dropped BEFORE we release it here.
        for r in list(self._live_results):
            try:
                r.close()
            except Exception:
                pass
        self._live_results.clear()
        self.close_output()
        self._pfapp = None
        self._cursor = None
        self._config = None
        self._closed = True

    def __enter__(self) -> "PowerFlow":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        self.close()
        return False


__all__ = ["PowerFlow"]
