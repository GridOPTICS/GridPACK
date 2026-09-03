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

from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


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
        solver_converged: Optional[bool] = None,
    ) -> None:
        self._pfapp = pfapp
        self.nonlinear = nonlinear
        self.input_file = input_file

        # Snapshot now: re-solving overwrites it.  nl_solve never writes
        # p_convergence, so that path falls back to the solver's return.
        self._solver_converged = solver_converged
        self._convergence = None if nonlinear else pfapp.getConvergence()

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
    # Convergence
    # ------------------------------------------------------------------

    @property
    def convergence(self):
        """Iteration history, or None on the non-linear path."""
        return self._convergence

    @property
    def converged(self) -> Optional[bool]:
        """Whether the tolerance was reached.  Rank-uniform."""
        if self._convergence is not None:
            return self._convergence.converged
        return self._solver_converged

    @property
    def iterations(self) -> Optional[int]:
        """Newton iterations, or None on the non-linear path."""
        return None if self._convergence is None else self._convergence.iterations

    @property
    def mismatch(self):
        """Largest final P/Q mismatch in MW / MVAr, or None."""
        if self._convergence is None:
            return None
        if not self._convergence.perIteration:
            return None
        return self._convergence.finalMismatch

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
        conv = self.converged
        state = {True: "converged", False: "DIVERGED", None: "unknown"}[conv]
        iters = "" if self.iterations is None else f" iterations={self.iterations}"
        return (f"<PowerFlowResult {state}{iters} "
                f"nonlinear={self.nonlinear} input={self.input_file!r}>")


def _try_import_matplotlib():
    try:
        import matplotlib.pyplot as plt
        return plt
    except ImportError:  # pragma: no cover - only exercised without extras
        return None


# The seven parallel arrays returned by DSFullApp::getObservations.
# Order comes from the pybind11 tuple in gridpack_ds.cpp::getObservations.
_DSF_OBS_LABELS = (
    "vmag", "vangle", "rspeed", "rangle", "genP", "genQ", "fOnline",
)


class DSFResult:
    """Container for :class:`gridpack.DynamicSim` / :class:`DynamicSimStepper` output.

    Holds:

    * per-step time series of observations (list of tuples, one per step),
    * the observation-list metadata (which generators / loads / buses were
      configured under ``<observations>`` in the XML),
    * a live reference to the underlying pybind11 ``DSFullApp`` so on-demand
      queries (``get_generator_power``, ``get_bus_total_load_power``, ...)
      continue to work until the parent wrapper closes.

    The container itself is time-series first: each row is one simulation
    step.  Use :meth:`to_dataframe` for a flattened DataFrame (requires
    the ``[results]`` extra) or :meth:`plot` for a quick visualization
    (requires the ``[plot]`` extra).
    """

    def __init__(
        self,
        dsapp,
        *,
        input_file: Optional[str] = None,
        with_bus_freq: bool = False,
    ) -> None:
        self._dsapp = dsapp
        self.input_file = input_file
        self.with_bus_freq = with_bus_freq

        # (time, obs-tuple) rows.  Populated by DynamicSim / DynamicSimStepper.
        self.times: List[float] = []
        self.observations: List[Tuple] = []

        # Metadata from getObservationLists (populated lazily by the wrapper).
        self.obs_gen_buses: List[int] = []
        self.obs_gen_ids: List[str] = []
        self.obs_load_buses: List[int] = []
        self.obs_load_ids: List[str] = []
        self.obs_bus_ids: List[int] = []
        self.obs_bus_freq_ids: List[int] = []

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def close(self) -> None:
        """Drop the reference to the pybind11 application.

        Called automatically when the parent wrapper closes.  Subsequent
        queries that need the C++ app raise ``RuntimeError``; the
        already-collected time series remains accessible.
        """
        self._dsapp = None

    def _check(self):
        if self._dsapp is None:
            raise RuntimeError("DSFResult is closed")

    # ------------------------------------------------------------------
    # Time-series accessors
    # ------------------------------------------------------------------

    @property
    def n_steps(self) -> int:
        return len(self.observations)

    def channel_names(self) -> List[str]:
        """Flattened list of column names for :meth:`to_dataframe`."""
        cols: List[str] = ["time"]
        for bus, gid in zip(self.obs_gen_buses, self.obs_gen_ids):
            gid_s = str(gid).strip() or "1"
            cols.append(f"gen_{bus}_{gid_s}_rspeed")
            cols.append(f"gen_{bus}_{gid_s}_rangle")
            cols.append(f"gen_{bus}_{gid_s}_P")
            cols.append(f"gen_{bus}_{gid_s}_Q")
            cols.append(f"gen_{bus}_{gid_s}_online")
        for bus in self.obs_bus_ids:
            cols.append(f"bus_{bus}_vmag")
            cols.append(f"bus_{bus}_vangle")
        if self.with_bus_freq:
            for bus in self.obs_bus_freq_ids:
                cols.append(f"bus_{bus}_freq")
        return cols

    def _flatten_row(self, t: float, obs: Tuple) -> List[float]:
        # obs = (vMag, vAng, rSpd, rAng, genP, genQ, fOnline[, busfreq])
        # HADREC returns an empty tuple when the XML has no <observations>.
        empty: list = []
        padded = tuple(obs) + (empty,) * max(0, 7 - len(obs))
        vMag, vAng, rSpd, rAng, genP, genQ, fOnline = padded[:7]
        busfreq = obs[7] if self.with_bus_freq and len(obs) > 7 else []

        row: List[float] = [t]
        n_gen = len(self.obs_gen_buses)
        for i in range(n_gen):
            row.append(float(rSpd[i]) if i < len(rSpd) else float("nan"))
            row.append(float(rAng[i]) if i < len(rAng) else float("nan"))
            row.append(float(genP[i]) if i < len(genP) else float("nan"))
            row.append(float(genQ[i]) if i < len(genQ) else float("nan"))
            row.append(float(fOnline[i]) if i < len(fOnline) else float("nan"))
        n_bus = len(self.obs_bus_ids)
        for i in range(n_bus):
            row.append(float(vMag[i]) if i < len(vMag) else float("nan"))
            row.append(float(vAng[i]) if i < len(vAng) else float("nan"))
        if self.with_bus_freq:
            for i in range(len(self.obs_bus_freq_ids)):
                row.append(float(busfreq[i]) if i < len(busfreq) else float("nan"))
        return row

    def to_records(self) -> List[Dict[str, float]]:
        """Return the time series as a list of ``{column: value}`` dicts."""
        cols = self.channel_names()
        rows = []
        for t, obs in zip(self.times, self.observations):
            values = self._flatten_row(t, obs)
            rows.append(dict(zip(cols, values)))
        return rows

    def to_dataframe(self):
        """Return the time series as a pandas DataFrame.

        Requires the ``[results]`` extra (``pip install gridpack[results]``).
        """
        pd = _try_import_pandas()
        if pd is None:
            raise ImportError(
                "pandas is not installed. Install with "
                "'pip install gridpack[results]' or 'pip install pandas'."
            )
        cols = self.channel_names()
        data = [self._flatten_row(t, obs)
                for t, obs in zip(self.times, self.observations)]
        return pd.DataFrame(data, columns=cols)

    def to_csv(self, path: str) -> None:
        """Write the collected time series to CSV.

        Uses pandas if available, otherwise the stdlib ``csv`` module.
        """
        cols = self.channel_names()
        pd = _try_import_pandas()
        if pd is not None:
            self.to_dataframe().to_csv(path, index=False)
            return
        import csv
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(cols)
            for t, obs in zip(self.times, self.observations):
                w.writerow(self._flatten_row(t, obs))

    # ------------------------------------------------------------------
    # Live query passthroughs (require the app to still be open)
    # ------------------------------------------------------------------

    def get_generator_power(self, bus_id: int, gen_id: str = "1"):
        """Return ``(P, Q)`` for a generator, or ``None`` if not present."""
        self._check()
        return self._dsapp.getGeneratorPower(bus_id, gen_id)

    def get_bus_total_load_power(self, bus_id: int):
        """Return ``(P, Q)`` for a bus load, or ``None`` if not present."""
        self._check()
        return self._dsapp.getBusTotalLoadPower(bus_id)

    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------

    def plot(self, channels: Optional[Sequence[str]] = None, *, ax=None):
        """Plot one or more channels against time.

        Requires the ``[plot]`` extra (``pip install gridpack[plot]``).

        If ``channels`` is omitted, every generator's rotor speed
        (``gen_<bus>_<id>_rspeed``) is plotted.
        """
        plt = _try_import_matplotlib()
        if plt is None:
            raise ImportError(
                "matplotlib is not installed. Install with "
                "'pip install gridpack[plot]' or 'pip install matplotlib'."
            )
        cols = self.channel_names()
        if channels is None:
            channels = [c for c in cols if c.endswith("_rspeed")]

        if ax is None:
            _, ax = plt.subplots()
        rows = [self._flatten_row(t, obs)
                for t, obs in zip(self.times, self.observations)]
        col_idx = {c: i for i, c in enumerate(cols)}
        t_idx = col_idx["time"]
        times = [r[t_idx] for r in rows]
        for ch in channels:
            if ch not in col_idx:
                raise KeyError(f"Unknown channel {ch!r}. "
                               f"Available: {cols}")
            j = col_idx[ch]
            ax.plot(times, [r[j] for r in rows], label=ch)
        ax.set_xlabel("time [s]")
        ax.legend(loc="best", fontsize="small")
        return ax

    # ------------------------------------------------------------------
    # Repr
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        state = "closed" if self._dsapp is None else "open"
        return (f"<DSFResult {state} steps={self.n_steps} "
                f"gens={len(self.obs_gen_buses)} buses={len(self.obs_bus_ids)}>")


__all__ = ["PowerFlowResult", "DSFResult"]
