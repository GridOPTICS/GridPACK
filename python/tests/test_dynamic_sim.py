# -------------------------------------------------------------
# file: python/tests/test_dynamic_sim.py
# -------------------------------------------------------------
# Dynamic simulation integration tests: DynamicSim.run() (clean
# teardown) and DynamicSimStepper (needs os._exit(0) escape hatch
# to bypass upstream DSFullApp destructor SEGV).
# -------------------------------------------------------------

import pytest

from .conftest import run_inline


@pytest.mark.integration
def test_dynamic_sim_full_run(dsf_build_dir):
    r = run_inline(
        """
        from gridpack import Session, DynamicSim
        with Session() as s:
            ds = DynamicSim(s, "input_9b3g.xml", suppress_output=True)
            result = ds.run()
            if s.rank == 0:
                print(f"OK steps={result.n_steps}")
        """,
        cwd=dsf_build_dir,
        timeout=180,
    )
    assert r.returncode == 0, (
        f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    )
    assert "OK " in r.stdout


@pytest.mark.integration
def test_dynamic_sim_stepper(dsf_build_dir):
    r = run_inline(
        """
        import os, sys
        from gridpack import Session, DynamicSimStepper

        s = Session()
        try:
            st = DynamicSimStepper(s, "input_9b3g.xml", suppress_output=True)
            while not st.done:
                st.step()
            if s.rank == 0:
                print(f"STEPS={st.step_count}")
                print("OK")
            sys.stdout.flush()
        finally:
            # Bypass upstream DSFullApp destructor SEGV.
            os._exit(0)
        """,
        cwd=dsf_build_dir,
        timeout=180,
    )
    assert r.returncode == 0, (
        f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    )
    assert "OK" in r.stdout
    steps = _parse_int(r.stdout, "STEPS=")
    assert steps > 0


@pytest.mark.integration
def test_dynamic_sim_stepper_actuators(dsf_build_dir):
    r = run_inline(
        """
        import os, sys
        from gridpack import Session, DynamicSimStepper

        s = Session()
        try:
            st = DynamicSimStepper(s, "input_9b3g.xml", suppress_output=True)
            while not st.done:
                st.step()
                if st.step_count == 20:
                    st.apply_load_shedding(bus_number=5, percentage=-0.2)
                elif st.step_count == 40:
                    st.apply_line_trip(from_bus=6, to_bus=7, ckt="1 ")
                elif st.step_count == 60:
                    st.set_wide_area_signal(bus_number=2, gen_id="1", signal=0.05)
                elif st.step_count == 80:
                    st.apply_generator_trip(bus_number=3, gen_id="1")
            if s.rank == 0:
                print(f"STEPS={st.step_count}")
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


def _parse_int(text: str, prefix: str) -> int:
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(prefix):
            return int(line[len(prefix):])
    raise AssertionError(f"'{prefix}' not found in output:\n{text}")
