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


# Record fields, in C++ struct order.  Dicts rather than the pybind11
# records: records do not pickle, so they cannot cross MPI.
_BUS_FIELDS = (
    "busId", "type", "area", "zone", "baseKV", "voltage", "angle",
    "pInjection", "qInjection", "pLoad", "qLoad", "pGen", "qGen",
    "shuntMvar",
)
_BRANCH_FIELDS = (
    "fromBus", "toBus", "circuitId", "pFrom", "qFrom", "pTo", "qTo",
    "pLoss", "qLoss", "mvaFrom", "mvaTo", "rateA", "loadingPercent",
)
_GEN_FIELDS = (
    "busId", "genId", "pGen", "qGen", "qMax", "qMin", "voltageSetpoint",
    "status",
)

# Sort keys make row order independent of the rank count.
_BUS_KEY = lambda r: r["busId"]
_BRANCH_KEY = lambda r: (r["fromBus"], r["toBus"], r["circuitId"])
_GEN_KEY = lambda r: (r["busId"], r["genId"])

_TABLES = {
    "buses": (_BUS_FIELDS, _BUS_KEY),
    "branches": (_BRANCH_FIELDS, _BRANCH_KEY),
    "generators": (_GEN_FIELDS, _GEN_KEY),
}


def _records_to_dicts(records, fields):
    """Copy pybind11 records into dicts."""
    return [{f: getattr(rec, f) for f in fields} for rec in records]


class PowerFlowResult:
    """Solution container for :class:`gridpack.PowerFlow`.

    Holds a reference to the underlying pybind11 application module so
    per-bus queries stay live after :meth:`gridpack.PowerFlow.solve`
    returns.

    :meth:`buses`, :meth:`branches`, :meth:`generators` and
    :meth:`violations` return the whole network, gathered across ranks and
    sorted, so their output does not depend on the rank count.  They are
    COLLECTIVE -- every rank must call them.
    """

    def __init__(
        self,
        pfapp,
        *,
        nonlinear: bool,
        input_file: Optional[str] = None,
        solver_converged: Optional[bool] = None,
        mpi_comm: Optional[object] = None,
    ) -> None:
        self._pfapp = pfapp
        self.nonlinear = nonlinear
        self.input_file = input_file

        # Snapshot now: re-solving overwrites it.  nl_solve never writes
        # p_convergence, so that path falls back to the solver's return.
        self._solver_converged = solver_converged
        self._convergence = None if nonlinear else pfapp.getConvergence()

        # Gathering is collective, so cache each table after the first call.
        self._mpi_comm = mpi_comm
        self._tables: Dict[str, List[dict]] = {}

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
    # Gathered network tables
    # ------------------------------------------------------------------
    # Per-rank lists are a disjoint partition: concatenate, no dedup.
    # Sorting fixes row order across rank counts; values still move ~1e-14.

    def _table(self, name: str) -> List[dict]:
        """Gather one table across ranks; cached.

        COLLECTIVE: every rank must call, or the allgather deadlocks.
        """
        if name in self._tables:
            return self._tables[name]
        self._check()
        fields, key = _TABLES[name]

        results = self._pfapp.collectResults()
        rows = _records_to_dicts(getattr(results, name), fields)

        comm = self._mpi_comm
        if comm is not None and comm.Get_size() > 1:
            rows = [r for chunk in comm.allgather(rows) for r in chunk]
        rows.sort(key=key)

        self._tables[name] = rows
        return rows

    def buses(self) -> List[dict]:
        """All buses, sorted by bus number.  Collective."""
        return self._table("buses")

    def branches(self) -> List[dict]:
        """All branches, sorted by (from, to, circuit).  Collective."""
        return self._table("branches")

    def generators(self) -> List[dict]:
        """All generators, sorted by (bus, id).  Collective."""
        return self._table("generators")

    # ------------------------------------------------------------------
    # Violations
    # ------------------------------------------------------------------

    def violations(
        self,
        *,
        min_voltage: float = 0.9,
        max_voltage: float = 1.1,
        overload_threshold: float = 100.0,
    ) -> Dict[str, Any]:
        """Report buses outside the voltage band and overloaded branches.

        Defaults match ca_driver's minVoltage/maxVoltage.  Collective.

        ``unrated_branches`` counts branches with ``rateA <= 0``: without
        it, an unrated network looks like one with no overloads.

        Differs from the C++ checks: overload uses ``loadingPercent``
        (``max(mvaFrom, mvaTo)``) where ``checkLineOverloadViolations()``
        tests one end, so this can flag more; and limits are uniform,
        since per-bus ``BUS_VOLTAGE_MIN/MAX`` are not in ``BusResult``.
        """
        volt = [
            {"busId": b["busId"], "voltage": b["voltage"],
             "limit": max_voltage if b["voltage"] > max_voltage else min_voltage,
             "kind": "high" if b["voltage"] > max_voltage else "low"}
            for b in self.buses()
            if b["voltage"] > max_voltage or b["voltage"] < min_voltage
        ]

        branches = self.branches()
        rated = [br for br in branches if br["rateA"] > 0.0]
        over = [
            {"fromBus": br["fromBus"], "toBus": br["toBus"],
             "circuitId": br["circuitId"], "rateA": br["rateA"],
             "mva": max(br["mvaFrom"], br["mvaTo"]),
             "loadingPercent": br["loadingPercent"]}
            for br in rated
            if br["loadingPercent"] > overload_threshold
        ]

        return {
            "voltage": volt,
            "overload": over,
            "unrated_branches": len(branches) - len(rated),
            "n_branches": len(branches),
            "limits": {"min_voltage": min_voltage, "max_voltage": max_voltage,
                       "overload_threshold": overload_threshold},
        }

    # ------------------------------------------------------------------
    # Bulk extraction
    # ------------------------------------------------------------------

    def to_records(
        self,
        bus_ids: Optional[Sequence[int]] = None,
        *,
        table: str = "buses",
    ) -> List[dict]:
        """Return results as a list of dicts.

        With ``bus_ids``: the narrow rank-local view, dropping buses this
        rank does not own.  Without: the full gathered ``table``
        ("buses", "branches", "generators"), which is collective.
        """
        if bus_ids is None:
            if table not in _TABLES:
                raise ValueError(
                    f"unknown table {table!r}; expected one of "
                    f"{sorted(_TABLES)}"
                )
            return list(self._table(table))

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

    def to_dataframe(
        self,
        bus_ids: Optional[Sequence[int]] = None,
        *,
        table: str = "buses",
    ):
        """Return results as a pandas DataFrame.  See :meth:`to_records`.

        Requires the ``[results]`` extra (``pip install gridpack[results]``).
        """
        pd = _try_import_pandas()
        if pd is None:
            raise ImportError(
                "pandas is not installed. Install with "
                "'pip install gridpack[results]' or 'pip install pandas'."
            )
        rows = self.to_records(bus_ids, table=table)
        if bus_ids is None:
            # Fix column order even when the table is empty.
            return pd.DataFrame(rows, columns=list(_TABLES[table][0]))
        return pd.DataFrame(rows, columns=["bus", "vmag", "vangle"])

    def to_csv(
        self,
        path: str,
        bus_ids: Optional[Sequence[int]] = None,
        *,
        table: str = "buses",
    ) -> None:
        """Write results to a CSV file.  See :meth:`to_records`.

        Uses pandas if available, otherwise the stdlib ``csv`` module.
        """
        rows = self.to_records(bus_ids, table=table)
        cols = (list(_TABLES[table][0]) if bus_ids is None
                else ["bus", "vmag", "vangle"])

        pd = _try_import_pandas()
        if pd is not None:
            pd.DataFrame(rows, columns=cols).to_csv(path, index=False)
            return

        import csv
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)

    # ------------------------------------------------------------------
    # Passthrough
    # ------------------------------------------------------------------

    def write(self, path: Optional[str] = None) -> None:
        """Write the full solver output (mimics ``pf.x``).

        Goes to stdout unless ``path`` is given, in which case the C++ side
        redirects to that file and the redirect is closed again on return.
        """
        self._check()
        if path is None:
            self._pfapp.write()
            return
        self._pfapp.open(path)
        try:
            self._pfapp.write()
        finally:
            self._pfapp.close()

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
