# -------------------------------------------------------------
# file: python/tests/test_dynamic_sim.py
# -------------------------------------------------------------
# Dynamic simulation integration tests: DynamicSim.run() (clean
# teardown) and DynamicSimStepper (needs os._exit(0) escape hatch
# to bypass upstream DSFullApp destructor SEGV).
# -------------------------------------------------------------

import shutil

import pytest

from .conftest import run_inline


@pytest.mark.integration
def test_dynamic_sim_full_run(dsf_data_dir):
    r = run_inline(
        """
        from gridpack import Session, DynamicSim
        with Session() as s:
            ds = DynamicSim(s, "input_9b3g.xml", suppress_output=True)
            result = ds.run()
            if s.rank == 0:
                print(f"OK steps={result.n_steps}")
        """,
        cwd=dsf_data_dir,
        timeout=180,
    )
    assert r.returncode == 0, (
        f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    )
    assert "OK " in r.stdout


@pytest.mark.integration
def test_dynamic_sim_stepper(dsf_data_dir):
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
        cwd=dsf_data_dir,
        timeout=180,
    )
    assert r.returncode == 0, (
        f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    )
    assert "OK" in r.stdout
    steps = _parse_int(r.stdout, "STEPS=")
    assert steps > 0


@pytest.mark.integration
def test_dynamic_sim_stepper_actuators(dsf_data_dir):
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
        cwd=dsf_data_dir,
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


def _parse_float(text: str, prefix: str) -> float:
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(prefix):
            return float(line[len(prefix):])
    raise AssertionError(f"'{prefix}' not found in output:\n{text}")


@pytest.mark.integration
@pytest.mark.parametrize("xml_value,expected", [("false", False), ("true", True)])
def test_dynamic_sim_reads_xml_suppress_output(dsf_data_dir, tmp_path,
                                               xml_value, expected):
    """cursor.get() returns a string, so bool("false") was True."""
    for name in ("9b3g.raw", "9b3g.dyr"):
        shutil.copy(dsf_data_dir / name, tmp_path / name)
    xml = (dsf_data_dir / "input_9b3g.xml").read_text()
    (tmp_path / "input_9b3g.xml").write_text(xml.replace(
        "<Dynamic_simulation>",
        "<Dynamic_simulation>\n<suppressOutput>%s</suppressOutput>" % xml_value,
        1))

    r = run_inline("""
        from gridpack import Session, DynamicSim
        with Session() as s:
            print("SUPPRESS", DynamicSim(s, "input_9b3g.xml").suppress_output)
    """, cwd=tmp_path, timeout=180)
    assert r.returncode == 0, r.stderr[-2000:]
    assert ("SUPPRESS %s" % expected) in r.stdout


@pytest.mark.integration
@pytest.mark.parametrize("xml_value,nonlinear", [("false", False), ("true", True)])
def test_dynamic_sim_honors_xml_nonlinear(dsf_data_dir, tmp_path,
                                          xml_value, nonlinear):
    """initFromConfig always called solve(), so UseNonLinear was ignored.

    PETSc's SNES banner, which -ksp_view enables, is the only externally
    visible difference between the two power flow solvers.
    """
    for name in ("9b3g.raw", "9b3g.dyr"):
        shutil.copy(dsf_data_dir / name, tmp_path / name)
    xml = (dsf_data_dir / "input_9b3g.xml").read_text()
    xml = xml.replace("        -ksp_type richardson",
                      "        -ksp_view\n        -ksp_type richardson", 1)
    xml = xml.replace("<UseNonLinear>false</UseNonLinear>",
                      "<UseNonLinear>%s</UseNonLinear>" % xml_value, 1)
    (tmp_path / "input_9b3g.xml").write_text(xml)

    r = run_inline("""
        from gridpack import Session, DynamicSim
        with Session() as s:
            DynamicSim(s, "input_9b3g.xml")
    """, cwd=tmp_path, timeout=180)
    assert r.returncode == 0, r.stderr[-2000:]
    assert ("SNES" in r.stdout + r.stderr) is nonlinear


@pytest.mark.integration
def test_dynamic_sim_writes_power_flow_report(dsf_data_dir, tmp_path):
    """initFromConfig skipped pf_app.write(), dropping the PF tables."""
    for name in ("9b3g.raw", "9b3g.dyr", "input_9b3g.xml"):
        shutil.copy(dsf_data_dir / name, tmp_path / name)

    r = run_inline("""
        from gridpack import Session, DynamicSim
        with Session() as s:
            DynamicSim(s, "input_9b3g.xml")
    """, cwd=tmp_path, timeout=180)
    assert r.returncode == 0, r.stderr[-2000:]
    assert "Branch Power Flow" in r.stdout


# -------------------------------------------------------------
# Observation-vector layout
# -------------------------------------------------------------
# getObservations returns ONE FLAT vector grouped by quantity, not one
# group per generator.  The wrapper used to unpack it as a 7-tuple of
# lists, which raised TypeError on every run with a generator
# observation -- to_records, to_dataframe, to_csv and plot were all
# broken and nothing tested the flattening.  These run without MPI: the
# layout is the contract, and a container is enough to pin it.

def _result(gens=(("1", 1), ("2", 4)), loads=(("1", 6),), buses=(6, 10),
            freqs=(6,)):
    from gridpack import DSFResult
    r = DSFResult(None)
    r.obs_gen_ids = [g for g, _ in gens]
    r.obs_gen_buses = [b for _, b in gens]
    r.obs_load_ids = [l for l, _ in loads]
    r.obs_load_buses = [b for _, b in loads]
    r.obs_bus_ids = list(buses)
    r.obs_bus_freq_ids = list(freqs)
    return r


def test_channel_names_are_quantity_major():
    assert _result().channel_names() == [
        "time",
        "gen_1_1_rspeed", "gen_4_2_rspeed",
        "gen_1_1_rangle", "gen_4_2_rangle",
        "gen_1_1_P", "gen_4_2_P",
        "gen_1_1_Q", "gen_4_2_Q",
        "bus_6_vmag", "bus_10_vmag",
        "bus_6_vangle", "bus_10_vangle",
        "load_6_1_online",
        "bus_6_freq",
    ]


def test_online_channel_is_per_load_not_per_generator():
    """fOnline comes from getOnlineLoadFraction, one entry per observed load."""
    cols = _result().channel_names()
    assert [c for c in cols if c.endswith("_online")] == ["load_6_1_online"]


def test_bus_freq_columns_are_named_without_the_flag():
    """getObservations appends frequencies regardless of with_bus_freq."""
    r = _result(freqs=(6, 10))
    assert r.with_bus_freq is False
    assert [c for c in r.channel_names() if c.endswith("_freq")] == [
        "bus_6_freq", "bus_10_freq"]


def test_flatten_row_maps_a_flat_observation_vector():
    r = _result()
    cols = r.channel_names()
    obs = [float(i) for i in range(len(cols) - 1)]
    row = dict(zip(cols, r._flatten_row(0.5, obs)))
    assert row["time"] == 0.5
    # Gen-major unpacking would put 1.0 here.
    assert row["gen_4_2_rspeed"] == 1.0
    assert row["gen_1_1_P"] == 4.0
    assert row["bus_10_vmag"] == 9.0
    assert row["load_6_1_online"] == 12.0
    assert row["bus_6_freq"] == 13.0


def test_flatten_row_pads_short_vectors_with_nan():
    import math
    r = _result()
    row = r._flatten_row(0.0, [1.0, 2.0])
    assert len(row) == len(r.channel_names())
    assert math.isnan(row[-1])


def test_flatten_row_drops_values_past_the_named_columns():
    r = _result(freqs=())
    n = len(r.channel_names())
    row = r._flatten_row(0.0, [1.0] * (n + 5))
    assert len(row) == n


def test_to_records_accepts_a_flat_observation_vector():
    """The TypeError surfaced here first: 'float' has no len()."""
    r = _result()
    r.times = [0.0, 0.01]
    r.observations = [[1.0] * 14, [2.0] * 14]
    recs = r.to_records()
    assert [rec["time"] for rec in recs] == [0.0, 0.01]
    assert recs[1]["gen_1_1_rspeed"] == 2.0


# -------------------------------------------------------------
# Time axis
# -------------------------------------------------------------
# current_time was initialized to 0.0 and never advanced, so it always
# read 0.0 and result.times held step indices under a column named
# "time".  HADREC exposes no getTimeStep, so the wrapper reads
# <timeStep> from the XML.

@pytest.mark.integration
def test_stepper_time_axis_is_seconds(dsf_data_dir):
    r = run_inline(
        """
        import os, sys
        from gridpack import Session, DynamicSimStepper

        s = Session()
        try:
            st = DynamicSimStepper(s, "input_9b3g.xml", suppress_output=True)
            for _ in range(5):
                st.step()
            if s.rank == 0:
                print("DT=%r" % st.time_step)
                print("NOW=%r" % st.current_time)
                print("LAST=%r" % st.result.times[-1])
            sys.stdout.flush()
        finally:
            os._exit(0)
        """,
        cwd=dsf_data_dir,
        timeout=180,
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    # <timeStep> in input_9b3g.xml; step indices would give 1.0 and 5.0.
    assert _parse_float(r.stdout, "DT=") == pytest.approx(0.01)
    assert _parse_float(r.stdout, "NOW=") == pytest.approx(0.05)
    assert _parse_float(r.stdout, "LAST=") == pytest.approx(0.05)
