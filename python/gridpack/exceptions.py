# -------------------------------------------------------------
# file: gridpack/exceptions.py
# -------------------------------------------------------------
# Exceptions raised by the high-level GridPACK wrappers.
#
# Created September, 2026 by Yousu Chen
# -------------------------------------------------------------

from __future__ import annotations

from typing import Optional


class GridPACKError(Exception):
    """Base class for errors raised by the high-level wrappers."""


class PowerFlowDiverged(GridPACKError):
    """The solver did not reach the configured tolerance.

    Raised by :meth:`gridpack.PowerFlow.solve` unless ``strict=False``.
    Every rank raises: the verdict is globally reduced before it reaches
    Python, so this is safe to catch in MPI code.

    ``convergence``, ``mismatch``, ``iterations`` and ``final_tolerance``
    are None on the non-linear path, which never populates the summary.
    """

    def __init__(
        self,
        *,
        convergence=None,
        nonlinear: bool = False,
        tolerance: Optional[float] = None,
        max_iteration: Optional[int] = None,
        input_file: Optional[str] = None,
    ) -> None:
        self.convergence = convergence
        self.nonlinear = nonlinear
        self.tolerance = tolerance
        self.max_iteration = max_iteration
        self.input_file = input_file

        self.mismatch = getattr(convergence, "finalMismatch", None)
        self.iterations = getattr(convergence, "iterations", None)
        self.final_tolerance = getattr(convergence, "finalTolerance", None)

        super().__init__(self._build_message())

    def _build_message(self) -> str:
        solver = "nl_solve" if self.nonlinear else "solve"
        parts = ["power flow did not converge (%s)" % solver]

        if self.iterations is not None:
            it = "%d iteration%s" % (
                self.iterations, "" if self.iterations == 1 else "s"
            )
            if self.max_iteration is not None:
                it += " of %d allowed" % self.max_iteration
            parts.append(it)

        if self.final_tolerance is not None:
            tol = "mismatch %.6g" % self.final_tolerance
            if self.tolerance is not None:
                tol += " vs tolerance %.6g" % self.tolerance
            parts.append(tol)

        if self.mismatch is not None:
            parts.append(
                "worst dP %.6g MW at bus %d, dQ %.6g MVAr at bus %d" % (
                    self.mismatch.maxPMismatch, self.mismatch.maxPBus,
                    self.mismatch.maxQMismatch, self.mismatch.maxQBus,
                )
            )

        msg = "; ".join(parts)
        if self.input_file:
            msg += " [%s]" % self.input_file
        return msg


__all__ = ["GridPACKError", "PowerFlowDiverged"]
