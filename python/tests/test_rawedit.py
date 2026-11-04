# -------------------------------------------------------------
# file: python/tests/test_rawedit.py
# -------------------------------------------------------------
# gridpack.rawedit: typed edits to RAW records, every version.
# -------------------------------------------------------------

import shutil

import pytest

from gridpack.rawedit import (EditError, edit_raw, raw_blocks, raw_generators,
                              raw_version, _tokens)

from .conftest import _data_sets

V23, V33 = "IEEE14.raw", "kundur-twoarea_v33.raw"
V34 = "240busWECC_2018_PSS_fixedshunt.raw"
V36 = "240busWECC_2018_PSS_PQLoad_fixedshunt_noremotebus_yuan_v36.raw"


def _raw(name):
    return _data_sets() / "raw" / name


def _edit(tmp_path, name, edits):
    dst = tmp_path / name
    report = edit_raw(_raw(name), dst, edits)
    return report, _raw(name).read_text().splitlines(), dst.read_text().splitlines()


def test_tokens_keep_quotes_comments_and_spacing():
    fields, delim, tail = _tokens("  1,'A, B',  2.5 / note\n")
    assert (fields, delim, tail) == (["  1", "'A, B'", "  2.5 "], ",", "/ note\n")
    assert _tokens("1 2 3\r\n")[:2] == (["1", "2", "3"], " ")


@pytest.mark.parametrize("name, version, blocks", [
    (V23, 23, {"bus": 14, "generator": 5, "branch": 20}),
    (V33, 33, {"bus": 11, "load": 3, "fixed_shunt": 2, "generator": 4,
               "branch": 8, "transformer": 16}),
    (V34, 34, {"bus": 243, "load": 139, "fixed_shunt": 7, "generator": 146,
               "branch": 329}),
    (V36, 36, {"bus": 243, "load": 139, "generator": 146, "branch": 329}),
])
def test_blocks_by_label_or_position(name, version, blocks):
    assert raw_version(_raw(name)) == version
    got = raw_blocks(_raw(name))
    assert {k: len(got[k]) for k in blocks} == blocks
    assert len(raw_generators(_raw(name))) == blocks["generator"]


def test_v23_edits_touch_only_their_records(tmp_path):
    report, a, b = _edit(tmp_path, V23, [
        {"op": "set_load_active_power", "bus": 9, "mw": 60.0},
        {"op": "scale_load", "bus": 14, "factor": 1.2},
        {"op": "set_shunt", "bus": 9, "mvar": 25.0},
        {"op": "set_generator_active_power", "bus": 2, "mw": 60.0},
        {"op": "set_generator_voltage_setpoint", "bus": 6, "vm_pu": 1.05},
        {"op": "set_generator_in_service", "bus": 3, "in_service": False},
        {"op": "set_branch_in_service", "from": 5, "to": 1, "ckt": "BL",
         "in_service": False},
        {"op": "set_param", "kind": "generator", "bus": 3, "field": "QT",
         "value": 30.0}])
    assert [(r["was"], r["now"]) for r in report] == [
        ({"PL": 29.5}, {"PL": 60.0}),
        ({"PL": 14.9, "QL": 5.0}, {"PL": 17.88, "QL": 6.0}),
        ({"BL": 19.0}, {"BL": 25.0}),
        ({"PG": 40.0}, {"PG": 60.0}),
        ({"VS": 1.07}, {"VS": 1.05}),
        ({"STAT": 1}, {"STAT": 0}),
        ({"STAT": 1}, {"STAT": 0}),
        ({"QT": 40.0}, {"QT": 30.0})]
    changed = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
    assert len(a) == len(b) and len(changed) == 6      # bus 9 and gen 3 twice
    # Field widths are kept, so columns still line up.
    assert all(len(a[i]) == len(b[i]) for i in changed)
    assert "    60.000," in b[changed[0]] and "1.05000" in b[changed[4]]


def test_set_param_writes_integer_fields_as_integers(tmp_path):
    report, a, b = _edit(tmp_path, V23, [
        {"op": "set_param", "kind": "branch", "from": 1, "to": 2, "ckt": "BL",
         "field": "STAT", "value": 0.0}])
    assert report[0]["now"] == {"STAT": 0} and type(report[0]["now"]["STAT"]) is int
    changed = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
    assert len(changed) == 1 and b[changed[0]].rstrip().endswith(", 0")


def test_v33_edits(tmp_path):
    report, a, b = _edit(tmp_path, V33, [
        {"op": "set_load_active_power", "bus": 7, "mw": 900.0},
        {"op": "set_load_in_service", "bus": 8, "in_service": True},
        {"op": "set_shunt", "bus": 9, "mvar": 300.0, "mw": 1.5},
        {"op": "set_shunt_in_service", "bus": 7, "in_service": False},
        {"op": "set_branch_in_service", "from": 8, "to": 7, "ckt": "2",
         "in_service": False},
        {"op": "set_generator_in_service", "bus": 4, "in_service": False}])
    assert report[0]["was"] == {"PL": 967.0} and report[1]["was"] == {"STAT": 0}
    assert report[2]["now"] == {"BL": 300.0, "GL": 1.5}
    assert sum(x != y for x, y in zip(a, b)) == 6 and len(a) == len(b)


@pytest.mark.parametrize("name", [V34, V36])
def test_v34_and_v36_edits(tmp_path, name):
    report, a, b = _edit(tmp_path, name, [
        {"op": "set_generator_in_service", "bus": 1032, "id": "C",
         "in_service": False},
        {"op": "set_branch_in_service", "from": 1001, "to": 1201, "ckt": "1",
         "in_service": False},
        {"op": "set_load_active_power", "bus": 1101, "mw": 100.0}])
    assert [r["was"] for r in report] == [{"STAT": 1}, {"STAT": 1},
                                          {"PL": report[2]["was"]["PL"]}]
    assert [r["now"] for r in report] == [{"STAT": 0}, {"STAT": 0}, {"PL": 100.0}]
    assert sum(x != y for x, y in zip(a, b)) == 3 and len(a) == len(b)


@pytest.mark.parametrize("name, edits, msg", [
    (V23, [], "non-empty list"),
    (V23, [{"bus": 9}], "must be an object with an \"op\""),
    (V23, [{"op": "trip", "bus": 9}], "unknown op"),
    (V23, [{"op": "set_load_active_power", "bus": 9}], "missing mw"),
    (V23, [{"op": "set_load_active_power", "bus": 9, "mw": 1, "x": 2}],
     "unexpected x"),
    (V23, [{"op": "set_load_active_power", "bus": "9", "mw": 1}],
     "bus must be a bus number"),
    (V23, [{"op": "set_load_active_power", "bus": 9, "mw": "60"}],
     "mw must be a number"),
    (V23, [{"op": "set_branch_in_service", "from": 1, "to": 5, "ckt": "BL",
            "in_service": 0}], "in_service must be true or false"),
    (V23, [{"op": "set_load_active_power", "bus": 99, "mw": 1}],
     "no bus 99 in IEEE14.raw"),
    (V23, [{"op": "set_generator_active_power", "bus": 9, "mw": 1}],
     "no generator 9 id 1"),
    (V23, [{"op": "set_load_active_power", "bus": 9, "id": "2", "mw": 1}],
     "one load per bus"),
    (V23, [{"op": "set_load_in_service", "bus": 9, "in_service": False}],
     "no load status"),
    (V23, [{"op": "set_branch_in_service", "from": 1, "to": 5, "ckt": "1",
            "in_service": False}], "no branch 1-5 circuit 1"),
    (V23, [{"op": "set_param", "kind": "generator", "bus": 2, "field": "NREG",
            "value": 0}], "field NREG is not in"),
    (V23, [{"op": "set_param", "kind": "line", "bus": 2, "field": "R",
            "value": 0}], "kind must be one of"),
    (V23, [{"op": "set_shunt", "bus": 9}], "give mvar and/or mw"),
    (V33, [{"op": "set_branch_in_service", "from": 5, "to": 1, "ckt": "1",
            "in_service": False}], "transformers are in their own block"),
    (V33, [{"op": "set_shunt", "bus": 8, "mvar": 1.0}],
     "no fixed_shunt 8 id 1"),
])
def test_rejected_edits_write_nothing(tmp_path, name, edits, msg):
    # One good edit first: the bad one must still stop everything.
    good = {"op": "set_generator_active_power", "bus": 1, "mw": 50.0}
    with pytest.raises(EditError, match=msg):
        edit_raw(_raw(name), tmp_path / name, [good] + edits if edits else edits)
    assert not (tmp_path / name).exists()


def test_line_endings_and_untouched_bytes_survive(tmp_path):
    src = tmp_path / "crlf.raw"
    src.write_bytes(_raw(V23).read_bytes().replace(b"\n", b"\r\n"))
    edit_raw(src, tmp_path / "out.raw", [
        {"op": "set_generator_active_power", "bus": 2, "mw": 60.0}])
    out = (tmp_path / "out.raw").read_bytes()
    assert out.count(b"\r\n") == src.read_bytes().count(b"\r\n")
    a, b = src.read_bytes().split(b"\r\n"), out.split(b"\r\n")
    assert sum(x != y for x, y in zip(a, b)) == 1


def test_scale_everything_and_edits_see_each_other(tmp_path):
    report, a, b = _edit(tmp_path, V23, [
        {"op": "scale_all_loads", "factor": 1.1},
        {"op": "scale_all_generation", "factor": 0.9},
        {"op": "scale_load", "bus": 9, "factor": 2.0}])
    assert report[0]["was"] == {"PL": 259.0, "QL": 73.5}
    assert report[0]["now"] == {"PL": 284.9, "QL": 80.85} and report[0]["records"] == 14
    assert report[1]["was"] == {"PG": 272.4} and report[1]["records"] == 5
    # The later edit starts from the scaled value, not the file's.
    assert report[2]["was"] == {"PL": 32.45, "QL": 18.26} and "records" not in report[2]
    # Buses without load (1, 7, 8) and generators at 0 MW (3, 6, 8) are
    # rewritten unchanged.
    assert sum(x != y for x, y in zip(a, b)) == 11 + 2
    report, _, _ = _edit(tmp_path, V33, [{"op": "scale_all_loads", "factor": 1.1}])
    assert report[0]["records"] == 3 and report[0]["was"]["PL"] == 2734.0
