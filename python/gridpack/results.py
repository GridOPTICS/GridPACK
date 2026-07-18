# -------------------------------------------------------------
# file: gridpack/results.py
# -------------------------------------------------------------
# Results containers for high-level GridPACK wrappers.
#
# Pandas / matplotlib are OPTIONAL. Import them lazily so a plain
# `pip install gridpack` (without the [results] or [plot] extras) still
# works for users who only need the raw solver.
# -------------------------------------------------------------

from __future__ import annotations

from typing import Iterable, List, Optional, Sequence


def _try_import_pandas():
    try:
        import pandas as pd
        return pd
    except ImportError:  # pragma: no cover - only exercised without extras
        return None


class PowerFlowResult:
    """Solution container for :class:`gridpack.PowerFlow`.

    Holds a reference to the underlying pybind11 application module so
    per-bus queries stay live after :meth:`gridpack.PowerFlow.solve`
    returns.  Bulk extraction is opt-in via :meth:`to_dataframe` and
    :meth:`to_csv`, which take the bus IDs to include.
    """

    def __init__(
        self,
        pfapp,
        *,
        nonlinear: bool,
        input_file: Optional[str] = None,
    ) -> None:
        self._pfapp = pfapp
        self.nonlinear = nonlinear
        self.input_file = input_file

    def close(self) -> None:
        """Drop the reference to the pybind11 application.

        Called automatically when the parent :class:`PowerFlow` closes.
        Subsequent queries raise ``RuntimeError``.
        """
        self._pfapp = None

    def _check(self):
        if self._pfapp is None:
            raise RuntimeError("PowerFlowResult is closed")

    # ------------------------------------------------------------------
    # Per-bus query
    # ------------------------------------------------------------------

    def get_bus_solution(self, bus_id: int):
        """Return (vmag, vangle_deg) for the given original bus number.

        Returns ``None`` if the bus is not present on this rank.
        """
        self._check()
        return self._pfapp.getPFSolutionSingleBus(bus_id)

    # ------------------------------------------------------------------
    # Bulk extraction
    # ------------------------------------------------------------------

    def to_records(self, bus_ids: Sequence[int]) -> List[dict]:
        """Return a list of ``{'bus': id, 'vmag': ..., 'vangle': ...}``.

        Ranks other than the owner return ``None`` for a given bus; those
        rows are dropped.
        """
        self._check()
        rows = []
        for bid in bus_ids:
            sol = self._pfapp.getPFSolutionSingleBus(bid)
            if sol is None:
                continue
            vmag, vangle = sol
            rows.append({"bus": int(bid), "vmag": float(vmag),
                         "vangle": float(vangle)})
        return rows

    def to_dataframe(self, bus_ids: Sequence[int]):
        """Return per-bus results as a pandas DataFrame.

        Requires the ``[results]`` extra (``pip install gridpack[results]``).
        """
        pd = _try_import_pandas()
        if pd is None:
            raise ImportError(
                "pandas is not installed. Install with "
                "'pip install gridpack[results]' or 'pip install pandas'."
            )
        return pd.DataFrame(self.to_records(bus_ids))

    def to_csv(self, path: str, bus_ids: Sequence[int]) -> None:
        """Write per-bus solutions to a CSV file.

        Uses pandas if available, otherwise the stdlib ``csv`` module.
        """
        pd = _try_import_pandas()
        if pd is not None:
            self.to_dataframe(bus_ids).to_csv(path, index=False)
            return

        import csv
        rows = self.to_records(bus_ids)
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["bus", "vmag", "vangle"])
            w.writeheader()
            w.writerows(rows)

    # ------------------------------------------------------------------
    # Passthrough
    # ------------------------------------------------------------------

    def write(self) -> None:
        """Write the full solver output to stdout (mimics ``pf.x``)."""
        self._check()
        self._pfapp.write()

    def export_psse(self, path: str, version: int = 33) -> None:
        """Export the solved network to PSS/E format.

        Parameters
        ----------
        path : str
        version : {23, 33, 34}
        """
        self._check()
        if version == 23:
            self._pfapp.exportPSSE23(path)
        elif version == 33:
            self._pfapp.exportPSSE33(path)
        elif version == 34:
            self._pfapp.exportPSSE34(path)
        else:
            raise ValueError(f"Unsupported PSS/E version: {version}")

    def __repr__(self) -> str:
        return (f"<PowerFlowResult nonlinear={self.nonlinear} "
                f"input={self.input_file!r}>")


__all__ = ["PowerFlowResult"]
