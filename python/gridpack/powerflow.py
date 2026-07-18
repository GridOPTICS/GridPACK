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


def __getattr__(name):
    # Fallback for any pybind11 attribute not explicitly re-exported.
    return getattr(_pfmod, name)


class PowerFlow:
    """High-level power flow driver.

    Usage::

        from gridpack import Session, PowerFlow

        with Session() as s:
            pf = PowerFlow(s, "input_14.xml")
            result = pf.solve()
            result.write()                          # dump to stdout
            df = result.to_dataframe([1, 2, 3])      # per-bus DataFrame

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
        self._xml_nonlinear = bool(cursor.get("UseNonLinear"))
        xml_suppress = bool(cursor.get("suppressOutput"))
        self._suppress_output = (
            xml_suppress if suppress_output is None else suppress_output
        )

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

    def solve(self, nonlinear: Optional[bool] = None) -> PowerFlowResult:
        """Run the power flow solver and return a :class:`PowerFlowResult`.

        Parameters
        ----------
        nonlinear : bool, optional
            If given, override the XML ``UseNonLinear`` setting.
            ``True`` uses the math library non-linear solver
            (``pfapp.nl_solve``); ``False`` uses the custom Newton-Raphson
            loop (``pfapp.solve``).
        """
        use_nl = self._xml_nonlinear if nonlinear is None else nonlinear
        if use_nl:
            self._pfapp.nl_solve()
        else:
            self._pfapp.solve()
        self._pfapp.saveData()
        self._solved = True

        result = PowerFlowResult(
            self._pfapp,
            nonlinear=use_nl,
            input_file=self.input_file,
        )
        self._live_results.add(result)
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
