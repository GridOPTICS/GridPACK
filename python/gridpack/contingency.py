# -------------------------------------------------------------
# file: gridpack/contingency.py
# -------------------------------------------------------------
# N-1/N-k contingency analysis driver, following the C++ ca_driver.
#
# groupSize is fixed at 1: a contingency spanning two task groups needs
# cross-group bookkeeping this driver does not do.  Parallelism comes from
# the MPI rank count instead.
# -------------------------------------------------------------

from __future__ import annotations

import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import List, Optional

from .session import Session
from .powerflow import PowerFlow, _xml_bool, _xml_number


# -------------------------------------------------------------
# Contingency definitions
# -------------------------------------------------------------

@dataclass
class Contingency:
    """One contingency from a contingency-list XML.

    Line cases use ``from_buses``/``to_buses``/``circuit_ids``, generator
    cases ``buses``/``generator_ids``.  More than one element is N-k.
    """

    name: str
    kind: str
    from_buses: List[int] = None
    to_buses: List[int] = None
    circuit_ids: List[str] = None
    buses: List[int] = None
    generator_ids: List[str] = None

    @property
    def is_line(self) -> bool:
        return self.kind == "Line"

    def to_pybind(self):
        """Build the pybind11 ``powerflow.Contingency`` this describes."""
        from . import _gridpack as _ext

        ctg = _ext.powerflow.Contingency()
        ctg.p_name = self.name
        if self.is_line:
            # p_type: 1 = Branch, 0 = Generator.
            ctg.p_type = 1
            ctg.p_from = self.from_buses
            ctg.p_to = self.to_buses
            ctg.p_ckt = self.circuit_ids
            ctg.p_saveLineStatus = [1] * len(self.circuit_ids)
        else:
            ctg.p_type = 0
            ctg.p_busid = self.buses
            ctg.p_genid = self.generator_ids
            ctg.p_saveGenStatus = [1] * len(self.generator_ids)
        return ctg


@dataclass
class ContingencyResult:
    """Outcome of one contingency.

    ``found`` is False when no element matched the network: a no-op case,
    not a genuine divergence.
    """

    name: str
    found: bool
    converged: bool
    voltage_ok: bool = True
    overload_ok: bool = True

    @property
    def ok(self) -> bool:
        return self.found and self.converged and self.voltage_ok and self.overload_ok

    @property
    def status(self) -> str:
        if not self.found:
            return "NOT FOUND"
        if not self.converged:
            return "DIVERGENT"
        if not self.voltage_ok and not self.overload_ok:
            return "BUS+BRANCH VIOLATION"
        if not self.voltage_ok:
            return "BUS VIOLATION"
        if not self.overload_ok:
            return "BRANCH VIOLATION"
        return "OK"

    @property
    def report_lines(self) -> List[str]:
        """Lines for the ``.out`` file, worded as ``ca.x``.

        A list, not :attr:`status`: a doubly-violating case emits both.
        """
        if not self.found:
            return ["Not found for contingency %s" % self.name]
        if not self.converged:
            return ["Divergent for contingency %s" % self.name]
        out = []
        if not self.voltage_ok:
            out.append("Bus Violation for contingency %s" % self.name)
        if not self.overload_ok:
            out.append("Branch Violation for contingency %s" % self.name)
        return out or ["No violation for contingency %s" % self.name]


# -------------------------------------------------------------
# Contingency list parsing
# -------------------------------------------------------------

def _fixed_width_ids(raw: str) -> List[str]:
    """Pad/truncate each whitespace-separated id to the 2 chars PSS/E uses."""
    return [tok.ljust(2)[:2] for tok in raw.split()]


def parse_contingency_list(path: str) -> List[Contingency]:
    """Parse a contingency-list XML into :class:`Contingency` objects.

    ElementTree, not the GridPACK config reader: Boost property_tree does
    not validate closing-tag names, so typos parse silently.
    """
    root = ET.parse(path).getroot()

    out: List[Contingency] = []
    for elem in root.iter("Contingency"):
        kind = elem.findtext("contingencyType", "").strip()
        name = elem.findtext("contingencyName", "")

        if kind == "Line":
            bus_ids = [int(x) for x in
                       elem.findtext("contingencyLineBuses", "").split()]
            ckts = _fixed_width_ids(elem.findtext("contingencyLineNames", ""))
            # Buses are listed as flat from/to pairs, one pair per circuit.
            if len(bus_ids) < 2 * len(ckts):
                raise ValueError(
                    "contingency %r lists %d line names but only %d buses; "
                    "expected 2 buses per name" % (name, len(ckts), len(bus_ids))
                )
            out.append(Contingency(
                name=name,
                kind="Line",
                from_buses=[bus_ids[2 * i] for i in range(len(ckts))],
                to_buses=[bus_ids[2 * i + 1] for i in range(len(ckts))],
                circuit_ids=ckts,
            ))
        elif kind == "Generator":
            out.append(Contingency(
                name=name,
                kind="Generator",
                buses=[int(x) for x in
                       elem.findtext("contingencyBuses", "").split()],
                generator_ids=_fixed_width_ids(
                    elem.findtext("contingencyGenerators", "")),
            ))
        elif kind:
            raise ValueError(
                "contingency %r has unknown contingencyType %r "
                "(expected 'Line' or 'Generator')" % (name, kind)
            )

    return out


# -------------------------------------------------------------
# Driver
# -------------------------------------------------------------

class ContingencyAnalysis:
    """N-1/N-k contingency analysis driver (mirrors ``ca.x``).

    Usage::

        with Session() as s:
            ca = ContingencyAnalysis(s, "input_14.xml")
            for r in ca.run():
                print(r.name, r.status)

    Each rank gets a private network copy and pulls from a shared task
    queue, so ``mpiexec -np N`` gives N-way parallelism.  :meth:`run`
    returns this rank's share, :meth:`gather` all of them.

    ``voltage_limits`` and ``print_calc_files`` override the XML
    ``minVoltage``/``maxVoltage`` and ``printCalcFiles``;
    ``suppress_output`` also silences the per-contingency ``.out`` files.
    """

    def __init__(
        self,
        session: Session,
        input_file: str,
        *,
        voltage_limits: Optional[tuple] = None,
        print_calc_files: Optional[bool] = None,
        suppress_output: Optional[bool] = None,
    ) -> None:
        if session.closed:
            raise RuntimeError("Session is closed")

        from . import _gridpack as _ext

        self._session = session
        self.input_file = input_file

        # Own Configuration; PowerFlow opens the file again for its block.
        self._config = _ext.Configuration()
        self._config.open(input_file, session.comm)
        cursor = self._config.getCursor("Configuration.Contingency_analysis")

        xml_group_size = _xml_number(cursor.get("groupSize"), int, 1)
        if xml_group_size != 1:
            sys.stderr.write(
                "Warning: groupSize=%d is ignored; a contingency spanning two "
                "task groups is not supported. Using groupSize=1 -- get "
                "parallelism from the MPI rank count instead.\n" % xml_group_size
            )
        self.group_size = 1

        if voltage_limits is None:
            self.min_voltage = _xml_number(cursor.get("minVoltage"), float, 0.9)
            self.max_voltage = _xml_number(cursor.get("maxVoltage"), float, 1.1)
        else:
            self.min_voltage, self.max_voltage = (
                float(voltage_limits[0]), float(voltage_limits[1]))

        self.print_calc_files = (
            _xml_bool(cursor.get("printCalcFiles"), True)
            if print_calc_files is None else bool(print_calc_files)
        )
        # The Powerflow block owns qlim for the solve; the CA block mirrors it.
        self.check_qlim = _xml_bool(cursor.get("qlim"), True)

        self._contingency_file = cursor.get("contingencyList")
        self._auto_n1 = (_xml_bool(cursor.get("FullBranchN1"))
                         or _xml_bool(cursor.get("FullGeneratorN1")))

        # Each task needs its own network, or ranks running different
        # contingencies deadlock inside the collective solve.
        self._task_comm = session.comm.divide(self.group_size)
        self._pf = PowerFlow(
            session, input_file,
            suppress_output=suppress_output,
            comm=self._task_comm,
        )

        self._contingencies: Optional[List[Contingency]] = None
        self._results: List[ContingencyResult] = []
        self._base_result = None
        self._closed = False

        # Registered after PowerFlow so the LIFO drain closes us first.
        session.register(self)

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    @property
    def powerflow(self) -> PowerFlow:
        """The :class:`~gridpack.PowerFlow` driving each contingency."""
        self._require_open()
        return self._pf

    @property
    def base_result(self):
        """:class:`~gridpack.PowerFlowResult` for the base case, or None."""
        return self._base_result

    @property
    def contingencies(self) -> List[Contingency]:
        """The contingency list, parsed on first access.

        Empty if the XML names no ``contingencyList``; ``FullBranchN1`` /
        ``FullGeneratorN1`` auto-generation is not implemented.
        """
        if self._contingencies is None:
            self._contingencies = (
                parse_contingency_list(self._contingency_file)
                if self._contingency_file else []
            )
        return self._contingencies

    @property
    def requests_auto_n1(self) -> bool:
        """True if the XML asks for unimplemented N-1 auto-generation."""
        return self._auto_n1

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------

    def solve_base_case(self):
        """Solve the base case and exempt its pre-existing violations.

        Returns the :class:`~gridpack.PowerFlowResult`; check ``.converged``
        -- contingencies against a diverged base are meaningless.
        """
        self._require_open()
        self._pf.set_voltage_limits(self.min_voltage, self.max_voltage)
        result = self._pf.solve(strict=False)
        if result.converged:
            if self.check_qlim and not self._pf.check_qlim_violations():
                # PV -> PQ conversion changed the network; re-solve.
                self._pf.solve(strict=False)
            self._pf.ignore_voltage_violations()
        self._base_result = result
        return result

    def run(self, progress=None) -> List[ContingencyResult]:
        """Evaluate this rank's share of the contingencies.

        Solves the base case first if not already done, and raises if it
        diverged.  ``progress(task_id, contingency)`` runs before each
        solve; the driver prints nothing itself.
        """
        self._require_open()

        if self._base_result is None:
            self.solve_base_case()
        if not self._base_result.converged:
            raise RuntimeError(
                "contingency analysis base case did not converge; "
                "contingencies would be measured against an invalid solution"
            )

        contingencies = self.contingencies
        if not contingencies:
            return []

        from . import _gridpack as _ext

        taskmgr = _ext.TaskManager(self._session.comm)
        taskmgr.set(len(contingencies))
        task = _ext.TaskCounter()

        results: List[ContingencyResult] = []
        # nextTask returns the same id on every rank of task_comm, so a
        # multi-rank group evaluates one contingency together.
        while taskmgr.nextTask(self._task_comm, task):
            if progress is not None:
                progress(task.task_id, contingencies[task.task_id])
            results.append(self._run_one(contingencies[task.task_id]))

        taskmgr.printStats()
        self._results = results
        return results

    def _run_one(self, contingency: Contingency) -> ContingencyResult:
        pf = self._pf
        ctg = contingency.to_pybind()
        out_file = contingency.name.rstrip() + ".out"

        if self.print_calc_files:
            pf.open_output(out_file)
        try:
            pf.reset_voltages()
            found = pf.set_contingency(ctg)

            converged = False
            voltage_ok = overload_ok = True
            if found:
                result = pf.solve(strict=False)
                converged = bool(result.converged)
                if converged:
                    if self.check_qlim and not pf.check_qlim_violations():
                        converged = bool(pf.solve(strict=False).converged)
                    if self.print_calc_files:
                        result.write()
                    voltage_ok = pf.check_voltage_violations()
                    overload_ok = pf.check_line_overload_violations()

            outcome = ContingencyResult(
                name=contingency.name,
                found=found,
                converged=converged,
                voltage_ok=voltage_ok,
                overload_ok=overload_ok,
            )
            if self.print_calc_files:
                for line in outcome.report_lines:
                    pf.print_output("\n%s\n" % line)

            # Restore the network before the next task reuses it.
            pf.unset_contingency(ctg)
            if self.check_qlim:
                pf.clear_qlim_violations()
            return outcome
        finally:
            if self.print_calc_files:
                pf.close_output()

    def gather(self) -> List[ContingencyResult]:
        """All ranks' results, sorted by name.  Collective.

        Falls back to this rank's own results if mpi4py is unavailable.
        """
        comm = self._session.mpi_comm
        if comm is None or comm.Get_size() == 1:
            return sorted(self._results, key=lambda r: r.name)
        merged = [r for chunk in comm.allgather(self._results) for r in chunk]
        merged.sort(key=lambda r: r.name)
        return merged

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def close(self) -> None:
        """Idempotent; called automatically at Session close."""
        if self._closed:
            return
        try:
            self._pf.close()
        except Exception:
            pass
        self._pf = None
        self._base_result = None
        self._task_comm = None
        self._config = None
        self._closed = True

    def _require_open(self) -> None:
        if self._closed:
            raise RuntimeError("ContingencyAnalysis is closed")

    def __enter__(self) -> "ContingencyAnalysis":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        self.close()
        return False

    def __repr__(self) -> str:
        if self._closed:
            return "<gridpack.ContingencyAnalysis closed>"
        return ("<gridpack.ContingencyAnalysis input=%r n_contingencies=%d>"
                % (self.input_file, len(self.contingencies)))


__all__ = [
    "Contingency",
    "ContingencyResult",
    "ContingencyAnalysis",
    "parse_contingency_list",
]
