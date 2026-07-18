# -------------------------------------------------------------
# file: python/tests/test_hadrec.py
# -------------------------------------------------------------
# Integration tests for the high-level gridpack.Hadrec wrapper.
# Needs os._exit(0) to bypass the upstream DSFullApp destructor SEGV.
# -------------------------------------------------------------

import pytest

from .conftest import run_inline


@pytest.mark.integration
def test_hadrec_full_run(dsf_build_dir):
    r = run_inline(
        """
        import os, sys
        from gridpack import Session, Hadrec
        from gridpack.hadrec import make_line_trip_action

        s = Session(suppress_output=True)
        try:
            h = Hadrec(s, "input_9b3g.xml", suppress_output=True)
            result = h.initialize_dyn_simu()
            trip = make_line_trip_action(from_bus=6, to_bus=7, ckt="1 ")
            while not h.done:
                if h.step_count == 400:
                    h.apply_action(trip)
                h.step()
            if s.rank == 0:
                print(f"STEPS={h.step_count}")
                print(f"ROWS={result.n_steps}")
                gp = h.get_generator_power(1, "1")
                assert gp is not None, "generator (1,'1') missing"
                print(f"GENP={gp[0]}")
                print("OK")
            sys.stdout.flush()
        finally:
            os._exit(0)
        """,
        cwd=dsf_build_dir,
        timeout=180,
    )
    assert r.returncode == 0, (
        f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    )
    assert "OK" in r.stdout


@pytest.mark.integration
def test_hadrec_run_until_done(dsf_build_dir):
    r = run_inline(
        """
        import os, sys
        from gridpack import Session, Hadrec
        s = Session(suppress_output=True)
        try:
            h = Hadrec(s, "input_9b3g.xml", suppress_output=True)
            r = h.initialize_dyn_simu()
            result = h.run_until_done()
            if s.rank == 0:
                print(f"STEPS={h.step_count}")
                print(f"ROWS={result.n_steps}")
                print("OK")
            sys.stdout.flush()
        finally:
            os._exit(0)
        """,
        cwd=dsf_build_dir,
        timeout=180,
    )
    assert r.returncode == 0, (
        f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    )
    assert "OK" in r.stdout
