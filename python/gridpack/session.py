# -------------------------------------------------------------
# file: gridpack/session.py
# -------------------------------------------------------------
# GridPACK Session: MPI-safe context manager for Environment and
# Communicator lifecycle.
#
# The GridPACK C++ Environment must be constructed first and destroyed
# last; Communicator, Configuration, and any App modules must all be
# released before the Environment. Getting this wrong is the single
# most common cause of Python-side segfaults on interpreter shutdown.
# Session enforces the correct ordering.
# -------------------------------------------------------------

from __future__ import annotations

from typing import Optional


# Module-level guard: GridPACK's C++ Environment (MPI + PETSc + GA)
# cannot be initialized twice within the same process.  We enforce
# process-lifetime singleton semantics.
_SESSION_CREATED = False


class Session:
    """MPI-safe context manager for a GridPACK run.

    Usage::

        from gridpack import Session

        with Session() as s:
            print("rank", s.rank, "of", s.size)
            # build wrappers that take the session, e.g.
            #   result = PowerFlow(s, "input.xml").solve()

    Or imperative form::

        s = Session()
        try:
            ...
        finally:
            s.close()

    Parameters
    ----------
    mpi_comm : mpi4py.MPI.Comm, optional
        An externally created MPI communicator (from ``mpi4py``).  If
        given, GridPACK is initialized over this communicator instead
        of ``MPI_COMM_WORLD``.  Required for multi-physics / co-simulation
        scenarios that share an MPI world.
    suppress_output : bool, optional
        If True, enable ``gridpack.NoPrint`` immediately so the C++ side
        stays quiet.  Equivalent to ``NoPrint().setStatus(True)``.
    """

    def __init__(
        self,
        mpi_comm: Optional[object] = None,
        *,
        suppress_output: bool = False,
    ) -> None:
        # Local import so `import gridpack` doesn't recurse through this
        # module before the compiled extension is ready.
        from . import _gridpack as _ext

        global _SESSION_CREATED
        if _SESSION_CREATED:
            raise RuntimeError(
                "gridpack.Session has already been created in this process. "
                "GridPACK's MPI/PETSc/GA environment can only be initialized "
                "once per process; instantiate Session at most once."
            )
        _SESSION_CREATED = True

        # Environment MUST be constructed first.
        if mpi_comm is None:
            self._env = _ext.Environment()
        else:
            self._env = _ext.Environment(mpi_comm)

        self._comm = _ext.Communicator()

        if suppress_output:
            _ext.NoPrint().setStatus(True)

        self._closed = False

        # Registry of wrapper objects to close before we release the
        # Environment.  Wrappers with pybind11-owned resources (PowerFlow,
        # DynamicSim, etc.) call `session.register(self)` at construction
        # and their `.close()` runs here in LIFO order.  This makes the
        # teardown order explicit -- app modules first, then Communicator,
        # then Environment -- matching the deletion order in pf.py.
        self._resources: list = []

    # ------------------------------------------------------------------
    # Resource registry
    # ------------------------------------------------------------------

    def register(self, resource) -> None:
        """Register a wrapper to be ``close()``'d when this Session closes.

        The registry is drained in LIFO order.  A resource that has
        already been closed by the user is skipped safely as long as
        its ``close()`` is idempotent.
        """
        self._require_open()
        self._resources.append(resource)

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    @property
    def comm(self):
        """The GridPACK ``Communicator`` bound to this session."""
        self._require_open()
        return self._comm

    @property
    def env(self):
        """The GridPACK ``Environment`` bound to this session.

        Usually you do not need this directly; wrapper classes take the
        session and pull what they need.
        """
        self._require_open()
        return self._env

    @property
    def rank(self) -> int:
        """This process's rank within the session communicator."""
        return self.comm.rank()

    @property
    def size(self) -> int:
        """Number of processes in the session communicator."""
        return self.comm.size()

    @property
    def closed(self) -> bool:
        return self._closed

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def barrier(self) -> None:
        """Synchronize all ranks in the session communicator."""
        self.comm.barrier()

    def close(self) -> None:
        """Release Communicator, then Environment, in that order.

        Idempotent.  Safe to call from ``__exit__`` and from ``__del__``.
        After close(), most Session methods raise ``RuntimeError``.

        Runs a full ``gc.collect()`` before releasing MPI-owning objects
        so that any pybind11 app modules with cyclic references (e.g.
        PowerFlow -> PowerFlowResult -> pfapp) are torn down while the
        Communicator and Environment are still alive.  Without this,
        Python's shutdown sequence can release the app module after the
        Environment is gone, causing MPI_Iprobe on MPI_COMM_NULL crashes.
        """
        if self._closed:
            return
        # Drain registered resources in LIFO order.
        while self._resources:
            r = self._resources.pop()
            try:
                r.close()
            except Exception:
                pass
        # Release Communicator first, Environment last.  Matching the
        # explicit ordering used in the standalone scripts (pf.py etc.)
        # avoids MPI finalize / Boost.MPI teardown crashes.
        self._comm = None
        self._env = None
        self._closed = True

    def _require_open(self) -> None:
        if self._closed:
            raise RuntimeError("Session is closed")

    # ------------------------------------------------------------------
    # Context-manager protocol
    # ------------------------------------------------------------------

    def __enter__(self) -> "Session":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        self.close()
        return False

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    def __repr__(self) -> str:
        if self._closed:
            return "<gridpack.Session closed>"
        return f"<gridpack.Session rank={self.rank}/{self.size}>"


__all__ = ["Session"]
