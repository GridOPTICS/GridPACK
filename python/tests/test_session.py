# -------------------------------------------------------------
# file: python/tests/test_session.py
# -------------------------------------------------------------
# Session-level unit tests that need no MPI environment.
# -------------------------------------------------------------

# -------------------------------------------------------------
# mpi_comm resolution
# -------------------------------------------------------------
# Tested against a stub rather than a live Session: Session is a process
# singleton, and the guard is pure logic over two sizes.


class _FakeGPComm:
    def __init__(self, n):
        self._n = n

    def size(self):
        return self._n


class _FakeMPIComm:
    def __init__(self, n):
        self._n = n

    def Get_size(self):
        return self._n


def _resolve(gp_size, mpi_comm):
    from gridpack.session import Session

    stub = _FakeStub(gp_size, mpi_comm)
    return Session._resolve_mpi_comm(stub)


class _FakeStub:
    def __init__(self, gp_size, mpi_comm):
        self._comm = _FakeGPComm(gp_size)
        self._mpi_comm_arg = mpi_comm


def test_mpi_comm_accepts_matching_custom_comm():
    comm = _FakeMPIComm(4)
    assert _resolve(4, comm) is comm


def test_mpi_comm_rejects_size_mismatch():
    """Gathering over a larger comm would mix in ranks from another job."""
    import pytest

    with pytest.raises(RuntimeError, match="does not match"):
        _resolve(2, _FakeMPIComm(4))
