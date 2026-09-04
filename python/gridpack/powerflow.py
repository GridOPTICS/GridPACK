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
    """

    def __init__(
        self,
        session: Session,
        input_file: str,
        *,
        idx: int = -1,
        suppress_output: Optional[bool] = None,
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

        # Create the pybind11 Powerflow application.
        self._pfapp = _ext.powerflow.Powerflow()
        if self._suppress_output:
            self._pfapp.suppressOutput(True)

        self._pfapp.readNetwork(self._config, idx)
        self._pfapp.initialize()

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
            mpi_comm=self._session.mpi_comm,
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
