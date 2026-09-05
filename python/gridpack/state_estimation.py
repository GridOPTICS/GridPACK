# -------------------------------------------------------------
# file: gridpack/state_estimation.py
# -------------------------------------------------------------
# Re-exports the compiled pybind11 state_estimation submodule (so
# existing code like `import gridpack.state_estimation` /
# `gridpack.state_estimation.SEApp()` keeps working) and adds the
# high-level StateEstimation driver class on top.
#
# See the PowerFlow wrapper (powerflow.py) for the pattern followed
# here: register with Session, track live Results in a WeakSet, and
# null out the pybind11 app reference at close() to keep MPI finalize
# happy.
# -------------------------------------------------------------

from __future__ import annotations

import weakref
from typing import List, Optional, Sequence

from ._gridpack.state_estimation import *  # noqa: F401,F403
from ._gridpack import state_estimation as _mod  # noqa: F401

from .session import Session


def __getattr__(name):
    return getattr(_mod, name)


class StateEstimation:
    """High-level state estimation driver.

    Usage::

        from gridpack import Session, StateEstimation

        with Session() as s:
            se = StateEstimation(s, "input.xml")
            result = se.solve()
            if not result.converged:
                print("did not converge")
            result.write()

    Parameters
    ----------
    session : gridpack.Session
        Active session owning the MPI environment.
    input_file : str
        Path to the GridPACK state estimation XML configuration.
    measurement_file : str, optional
        Override the ``measurementList`` value in the XML.  If ``None``
        (default), the wrapper reads the filename from
        ``Configuration.State_estimation.measurementList``.
    suppress_output : bool, optional
        If True, enable ``NoPrint`` for the duration of the run.
    """

    def __init__(
        self,
        session: Session,
        input_file: str,
        *,
        measurement_file: Optional[str] = None,
        suppress_output: bool = False,
    ) -> None:
        if session.closed:
            raise RuntimeError("Session is closed")

        from . import _gridpack as _ext

        self._session = session
        self.input_file = input_file

        if suppress_output:
            _ext.NoPrint().setStatus(True)
        self._suppress_output = bool(suppress_output)

        # Open the main configuration on the session's Communicator.
        self._config = _ext.Configuration()
        self._config.open(input_file, session.comm)

        # Resolve the measurement file: CLI override wins, otherwise
        # read from the XML.
        if measurement_file is None:
            measurement_file = self._config.get(
                "Configuration.State_estimation.measurementList"
            )
            if not measurement_file:
                raise ValueError(
                    "State estimation XML does not specify "
                    "'Configuration.State_estimation.measurementList' and "
                    "no measurement_file override was provided."
                )
        self.measurement_file = measurement_file

        # Load measurements through a second Configuration.  Keep it
        # alive on `self` because SEApp caches pointers into cursors.
        self._mconfig = _ext.Configuration()
        self._mconfig.open(measurement_file, session.comm)

        self._seapp = _ext.state_estimation.SEApp()
        self._measurements = self._seapp.getMeasurements(self._mconfig)
        self._seapp.readNetwork(self._config)
        self._seapp.initialize()
        self._seapp.setMeasurements(self._measurements)

        self._solved = False
        self._closed = False

        # Track live results so we null out their _seapp reference
        # before we release ours (same pattern as PowerFlow).
        self._live_results: "weakref.WeakSet[StateEstimationResult]" = (
            weakref.WeakSet()
        )

        session.register(self)

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    @property
    def measurements(self):
        """The list of :class:`Measurement` objects being solved against."""
        return self._measurements

    @property
    def suppress_output(self) -> bool:
        return self._suppress_output

    # ------------------------------------------------------------------
    # Solve
    # ------------------------------------------------------------------

    def solve(self) -> "StateEstimationResult":
        """Run the state estimation solver and return a
        :class:`StateEstimationResult`.
        """
        if self._closed:
            raise RuntimeError("StateEstimation is closed")
        self._seapp.solve()
        self._seapp.saveData()
        converged = bool(self._seapp.hasConverged())
        self._solved = True

        result = StateEstimationResult(
            self._seapp,
            converged=converged,
            measurements=self._measurements,
            input_file=self.input_file,
        )
        self._live_results.add(result)
        return result

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def close(self) -> None:
        """Release the pybind11 app in the correct order.

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
        self._seapp = None
        self._mconfig = None
        self._config = None
        self._closed = True

    def __enter__(self) -> "StateEstimation":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        self.close()
        return False


class StateEstimationResult:
    """Solution container for :class:`StateEstimation`.

    Holds a reference to the underlying pybind11 SEApp so
    :meth:`write` continues to work until the parent wrapper closes.
    """

    def __init__(
        self,
        seapp,
        *,
        converged: bool,
        measurements,
        input_file: Optional[str] = None,
    ) -> None:
        self._seapp = seapp
        self.converged = bool(converged)
        self._measurements = measurements
        self.input_file = input_file

    def close(self) -> None:
        """Drop the reference to the pybind11 SEApp.

        Called automatically when the parent :class:`StateEstimation`
        closes.  Subsequent calls to :meth:`write` raise ``RuntimeError``.
        """
        self._seapp = None

    def _check(self):
        if self._seapp is None:
            raise RuntimeError("StateEstimationResult is closed")

    @property
    def measurements(self):
        """The list of :class:`Measurement` objects used for the solve."""
        return self._measurements

    def has_converged(self) -> bool:
        return self.converged

    def write(self, path: Optional[str] = None) -> None:
        """Write final results of state estimation to standard output.

        ``path`` is rejected: SEApp has no ``open``/``close`` binding, so
        there is nothing to redirect through.  Raising here beats the
        AttributeError that reaching for ``seapp.open`` used to produce
        *after* the solve had already printed to stdout.
        """
        self._check()
        if path is not None:
            raise NotImplementedError(
                "state estimation output cannot be redirected to a file; "
                "the SEApp binding exposes no open()/close(). Redirect the "
                "process's stdout instead."
            )
        self._seapp.write()

    def __repr__(self) -> str:
        state = "closed" if self._seapp is None else "open"
        return (
            f"<StateEstimationResult {state} "
            f"converged={self.converged} input={self.input_file!r}>"
        )


__all__ = ["StateEstimation", "StateEstimationResult"]
