# -------------------------------------------------------------
# file: python/tests/test_report.py
# -------------------------------------------------------------
# gridpack.mcp.report rendering and gridpack.mcp.plots paper figures.
# -------------------------------------------------------------

import math

import pytest

from gridpack.mcp import plots, report


def _dyn(n=51, step=0.01):
    cols = ["t", "1_1_speed", "1_1_angle", "2_1_speed", "2_1_angle"]
    rows = [[round(i * step, 6), 1.0 + 0.002 * math.sin(i / 5.0), 0.5 + 0.01 * i,
             1.0, 0.3] for i in range(n)]
    return {"watch": {"columns": cols, "rows": rows}}


def _pf():
    return {"converged": True, "iterations": 2,
            "buses": [{"busId": b, "name": "B%d" % b, "type": 1, "voltage": 1.0 + 0.01 * b,
                       "angle": -b, "pGen": 0.0, "qGen": 0.0, "pLoad": 1.0,
                       "qLoad": 0.0, "shuntMvar": 0.0} for b in range(1, 6)],
            "branches": [], "generators": []}


SECTIONS = [("Run", [("kv", {"case": "x.xml", "ranks": 3, "note": None}),
                     ("text", "a <b> text")]),
            ("Verdict", [("verdict", "line one\nline two")]),
            ("Tables", [("table", ["bus", "v"], [[1, 1.0195138], [2, None]]),
                        ("pre", "{\"edits\": []}")]),
            ("Figures", [("figure", "speeds", b"\x89PNG\r\n\x1a\nfake"),
                         ("figure_file", "angles", "r-angles.png")])]


def test_html_report_is_self_contained_and_escaped():
    h = report.render("html", "Run <1>", SECTIONS)
    assert h.startswith("<!doctype html>") and "<title>Run &lt;1&gt;</title>" in h
    assert "a &lt;b&gt; text" in h and "<style>" in h
    assert "line one<br>line two" in h
    assert "<td>1.01951</td>" in h and "<td></td>" in h          # 6 significant digits, None blank
    assert "data:image/png;base64,iVBORw0KGg" in h               # inline PNG
    assert "src='r-angles.png'" in h
    assert "<th>case</th><td>x.xml</td>" in h


def test_markdown_report_links_figures():
    m = report.render("md", "Run 1", SECTIONS)
    assert m.startswith("# Run 1\n") and "## Verdict\n\n> line one\n> line two" in m
    assert "| bus | v |\n|---|---|\n| 1 | 1.01951 |\n| 2 |  |" in m
    assert "![angles](r-angles.png)" in m and "```\n{\"edits\": []}\n```" in m
    assert "![speeds]((embedded figure))" in m


def test_long_tables_are_truncated_with_a_note():
    rows = [[i, 1.0] for i in range(1000)]
    h = report.render("html", "T", [("S", [("table", ["i", "v"], rows)])])
    assert h.count("<tr>") == 401 and "400 of 1000 rows shown" in h
    with pytest.raises(ValueError, match="format is one of"):
        report.render("pdf", "T", [])


def test_paper_figures_in_every_format(tmp_path):
    pytest.importorskip("matplotlib")
    pdf = plots.save(tmp_path / "f.pdf", [_dyn()], ["dsf"], channels=["*_speed", "*_angle"],
                     events_list=[[{"start": 0.1, "end": 0.15}]], width="single")
    assert pdf[0] == 3.5 and (tmp_path / "f.pdf").read_bytes()[:5] == b"%PDF-"
    svg = plots.save(tmp_path / "f.svg", [_dyn(), _dyn()], ["dsf", "dsf"],
                     channels=["1_1_speed"], labels=["base", "modified"], fmt="svg",
                     width="double", height=2.5)
    assert svg == (7.16, 2.5) and b"<svg" in (tmp_path / "f.svg").read_bytes()[:400]
    png = plots.save(tmp_path / "f.png", [_pf()], ["powerflow"], fmt="png", width=4.0,
                     style="screen", font="sans")
    assert png[0] == 4.0 and (tmp_path / "f.png").read_bytes()[:4] == b"\x89PNG"


def test_paper_style_dashes_and_no_title():
    pytest.importorskip("matplotlib")
    fig, plt = plots.build([_dyn()], ["dsf"], channels=["*_speed"], style="paper",
                           width="single", title=None)
    ax = fig.axes[0]
    assert ax.get_title(loc="left") == "" and len(ax.lines) == 2
    assert ax.lines[0].get_linestyle() == "-" and ax.lines[1].get_linestyle() == "--"
    assert ax.get_ylabel() == "rotor speed (p.u.)" and ax.get_xlabel() == "time (s)"
    assert plt.rcParams["font.size"] == 8 and plt.rcParams["pdf.fonttype"] == 42
    plt.close(fig)
    # Overlay: a color per job, a dash per channel, legend names both.
    fig, plt = plots.build([_dyn(), _dyn()], ["dsf", "dsf"],
                           channels=["1_1_speed", "2_1_speed"], labels=["a", "b"],
                           style="paper")
    ax = fig.axes[0]
    colors = [l.get_color() for l in ax.lines]
    assert colors[0] == colors[1] != colors[2] == colors[3]
    assert [l.get_linestyle() for l in ax.lines] == ["-", "--", "-", "--"]
    plt.close(fig)
    # One channel, several jobs: the dash follows the job instead.
    fig, plt = plots.build([_dyn(), _dyn()], ["dsf", "dsf"], channels=["1_1_speed"],
                           labels=["a", "b"], style="paper")
    assert [l.get_linestyle() for l in fig.axes[0].lines] == ["-", "--"]
    assert [l.get_label() for l in ax.lines] == ["a 1_1_speed", "a 2_1_speed",
                                                 "b 1_1_speed", "b 2_1_speed"]
    plt.close(fig)
    fig, plt = plots.build([_dyn()], ["dsf"], channels=["*_speed"], style="screen")
    assert fig.axes[0].get_title(loc="left") == "Dynamic simulation run, 2 channels"
    plt.close(fig)


@pytest.mark.parametrize("kwargs, msg", [
    ({"results": [_dyn()], "apps": ["dsf"], "style": "glossy"}, "style is one of"),
    ({"results": [_dyn(), _pf()], "apps": ["dsf", "powerflow"]}, "same kind"),
    ({"results": [_dyn()], "apps": ["dsf"], "width": -1}, "width is single or double"),
    ({"results": [_dyn()], "apps": ["dsf"], "labels": ["a", "b"]}, "one label per job"),
    ({"results": [_dyn(), {"watch": {"columns": ["t", "9_1_speed"], "rows": [[0, 1]]}}],
      "apps": ["dsf", "dsf"], "channels": ["1_1_speed"]}, "not every job has channel"),
    ({"results": [], "apps": []}, "one result per job"),
])
def test_figure_rejections(kwargs, msg):
    pytest.importorskip("matplotlib")
    with pytest.raises(ValueError, match=msg):
        plots.build(**kwargs)
    with pytest.raises(ValueError, match="format is one of"):
        plots.save("/tmp/x.gif", [_dyn()], ["dsf"], fmt="gif")
