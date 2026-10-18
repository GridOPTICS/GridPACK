# -------------------------------------------------------------
# file: python/tests/test_dyredit.py
# -------------------------------------------------------------
# gridpack.dyredit: DYR records, synthesis and edits, no GridPACK run.
# -------------------------------------------------------------

import pytest

from gridpack.dyredit import (FIELDS, MODEL_SETS, EditError, dyr_records,
                              edit_dyr, synthesize_dyr)
from .conftest import _data_sets


def _dyr(name):
    return _data_sets() / "dyr" / name


def test_records_of_every_bundled_style():
    recs = dyr_records(_dyr("9b3g.dyr"))                  # commas, quoted ids
    assert [(r.bus, r.model, r.id) for r in recs] == [(1, "GENCLS", "1"),
                                                      (2, "GENCLS", "1"),
                                                      (3, "GENCLS", "1")]
    assert recs[0].get("H") == 8.0 and recs[0].get("XDP") == 0.05
    kundur = dyr_records(_dyr("kundur-twoarea_esst1a.dyr"))   # blanks, multi-line
    assert len(kundur) == 12 and {r.model for r in kundur} == {
        "GENROU", "ESST1A", "SEXS", "TGOV1"}
    g1 = next(r for r in kundur if r.model == "GENROU" and r.bus == 1)
    assert (g1.get("TDOP"), g1.get("H"), g1.get("XD"), g1.get("XL")) == (
        8.0, 6.5, 1.8, 0.2)
    wecc = dyr_records(_dyr("240busWECC_2018_PSS_mod.dyr"))
    assert len(wecc) == 438 and wecc[0].id == "C"             # lettered ids
    # Unknown models still expose their values by position.
    e = next(r for r in kundur if r.model == "ESST1A")
    assert e.names()[:2] == ["p3", "p4"] and isinstance(e.get("p3"), float)


def test_synthesized_models_parse_back():
    text = synthesize_dyr([(1, "1"), (7, "G2")], "detailed")
    # No comment line (the parser aborts on one); ids padded to 2 chars.
    assert text.startswith("1, 'GENROU', '1 ', 6, 0.05, 1, 0.05, 4, 0, 1.8")
    import tempfile, pathlib
    p = pathlib.Path(tempfile.mkdtemp()) / "s.dyr"
    p.write_text(text)
    recs = dyr_records(p)
    assert [(r.bus, r.model, r.id) for r in recs] == [
        (1, "GENROU", "1"), (1, "SEXS", "1"), (1, "TGOV1", "1"),
        (7, "GENROU", "G2"), (7, "SEXS", "G2"), (7, "TGOV1", "G2")]
    for r in recs:
        names = FIELDS[r.model]
        params = dict(MODEL_SETS["detailed"])[r.model]
        assert [r.get(n) for n in names] == params
    cls = dyr_records(_write(synthesize_dyr([(3, "1")], "classical")))
    assert [(cls[0].model, cls[0].get("H"), cls[0].get("D"), cls[0].get("XDP"))] == [
        ("GENCLS", 5.0, 2.0, 0.3)]
    with pytest.raises(EditError, match="model_set is one of"):
        synthesize_dyr([(1, "1")], "fancy")


def _write(text):
    import tempfile, pathlib
    p = pathlib.Path(tempfile.mkdtemp()) / "t.dyr"
    p.write_text(text)
    return p


def test_edits_touch_only_their_fields(tmp_path):
    report = edit_dyr(_dyr("kundur-twoarea_esst1a.dyr"), tmp_path / "k.dyr", [
        {"bus": 1, "model": "GENROU", "field": "H", "value": 9.0},
        {"bus": 1, "model": "sexs", "field": "K", "value": 150},
        {"bus": 3, "id": "1", "model": "TGOV1", "field": "R", "value": 0.03}])
    assert [(r["was"], r["now"]) for r in report] == [(6.5, 9.0), (200.0, 150.0),
                                                      (0.04, 0.03)]
    a = _dyr("kundur-twoarea_esst1a.dyr").read_text().splitlines()
    b = (tmp_path / "k.dyr").read_text().splitlines()
    assert len(a) == len(b) and sum(x != y for x, y in zip(a, b)) == 3
    recs = dyr_records(tmp_path / "k.dyr")
    assert next(r for r in recs if r.model == "GENROU" and r.bus == 1).get("H") == 9.0
    assert next(r for r in recs if r.model == "GENROU" and r.bus == 2).get("H") == 6.5
    # Two edits to one record.
    edit_dyr(_dyr("9b3g.dyr"), tmp_path / "n.dyr", [
        {"bus": 2, "model": "GENCLS", "field": "H", "value": 12.5},
        {"bus": 2, "model": "GENCLS", "field": "XDP", "value": 0.07}])
    r2 = dyr_records(tmp_path / "n.dyr")[1]
    assert (r2.get("H"), r2.get("D"), r2.get("XDP")) == (12.5, 5.0, 0.07)


@pytest.mark.parametrize("edits, msg", [
    ([], "non-empty list"),
    ([{"bus": 1, "model": "GENCLS", "field": "H"}], "missing value"),
    ([{"bus": 1, "model": "GENCLS", "field": "H", "value": "9"}], "must be a number"),
    ([{"bus": 9, "model": "GENCLS", "field": "H", "value": 9}],
     "no GENCLS record for generator 9 id 1"),
    ([{"bus": 1, "model": "GENROU", "field": "H", "value": 9}],
     "no GENROU record for generator 1 id 1; it has GENCLS"),
    ([{"bus": 1, "model": "GENCLS", "field": "XD", "value": 9}],
     "GENCLS has no field XD; fields are H, D, XDP"),
    ([{"bus": 1, "model": "GENCLS", "field": "p9", "value": 9}], "has no field p9"),
])
def test_rejected_edits_write_nothing(tmp_path, edits, msg):
    good = {"bus": 1, "model": "GENCLS", "field": "H", "value": 7.0}
    with pytest.raises(EditError, match=msg):
        edit_dyr(_dyr("9b3g.dyr"), tmp_path / "x.dyr", [good] + edits if edits else edits)
    assert not (tmp_path / "x.dyr").exists()
