# -------------------------------------------------------------
# file: gridpack/mcp/report.py
# -------------------------------------------------------------
# Self-contained HTML and Markdown reports of jobs and studies.
# -------------------------------------------------------------

"""Everything the server already knows about a run, in one file a reader
can open without the server: case, how it was derived, events, settings,
the verdict, tables and figures (inline in HTML, beside the file in
Markdown).
"""

import base64
import datetime
import html
import json

__all__ = ["FORMATS", "render"]

FORMATS = ("html", "md")
_ROWS = 400                   # table rows shown before truncating

_CSS = """
body{font:14px/1.5 -apple-system,Segoe UI,Helvetica,Arial,sans-serif;color:#15181c;
background:#fff;max-width:980px;margin:32px auto;padding:0 16px}
h1{font-size:22px;margin:0 0 4px}h2{font-size:16px;margin:28px 0 8px;border-bottom:1px solid #e3e6ea;padding-bottom:4px}
.sub{color:#5b6470;margin:0 0 18px}
table{border-collapse:collapse;font-size:13px;font-variant-numeric:tabular-nums;margin:8px 0}
th,td{padding:3px 10px;border-bottom:1px solid #e3e6ea;text-align:right;white-space:nowrap}
th{color:#5b6470;font-weight:500}td:first-child,th:first-child{text-align:left}
pre{background:#f5f6f7;padding:10px 12px;border-radius:4px;font-size:12.5px;overflow-x:auto;white-space:pre-wrap}
.verdict{background:#f5f6f7;border-left:3px solid #2a78d6;padding:10px 14px;margin:8px 0}
.note{color:#5b6470;font-size:13px}img{max-width:100%;display:block;margin:10px 0}
.kv td{text-align:left}
"""


def render(fmt, title, sections):
    """A document from ``sections``: ``[(heading, [block, ...])]`` where a
    block is ``("text", str)``, ``("verdict", str)``, ``("kv", dict)``,
    ``("table", columns, rows)``, ``("pre", str)``, ``("figure", caption,
    png bytes)`` or ``("figure_file", caption, filename)``."""
    if fmt not in FORMATS:
        raise ValueError("format is one of %s, not %r" % (", ".join(FORMATS), fmt))
    stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    return (_html if fmt == "html" else _md)(title, stamp, sections)


def _fmt(v):
    if isinstance(v, float):
        return "%.6g" % v
    if isinstance(v, (dict, list)):
        return json.dumps(v)
    return "" if v is None else str(v)


def _html(title, stamp, sections):
    out = ["<!doctype html><html><head><meta charset='utf-8'><title>%s</title>"
           "<style>%s</style></head><body>" % (html.escape(title), _CSS),
           "<h1>%s</h1><p class='sub'>GridPACK report, %s</p>" % (html.escape(title), stamp)]
    for heading, blocks in sections:
        out.append("<h2>%s</h2>" % html.escape(heading))
        for b in blocks:
            kind = b[0]
            if kind == "text":
                out.append("<p>%s</p>" % html.escape(b[1]))
            elif kind == "note":
                out.append("<p class='note'>%s</p>" % html.escape(b[1]))
            elif kind == "verdict":
                out.append("<div class='verdict'>%s</div>" % "<br>".join(
                    html.escape(l) for l in b[1].splitlines()))
            elif kind == "kv":
                out.append("<table class='kv'>" + "".join(
                    "<tr><th>%s</th><td>%s</td></tr>" % (html.escape(str(k)), html.escape(_fmt(v)))
                    for k, v in b[1].items()) + "</table>")
            elif kind == "table":
                cols, rows = b[1], b[2]
                out.append("<table><tr>" + "".join("<th>%s</th>" % html.escape(str(c)) for c in cols)
                           + "</tr>" + "".join(
                    "<tr>" + "".join("<td>%s</td>" % html.escape(_fmt(v)) for v in r) + "</tr>"
                    for r in rows[:_ROWS]) + "</table>")
                if len(rows) > _ROWS:
                    out.append("<p class='note'>%d of %d rows shown.</p>" % (_ROWS, len(rows)))
            elif kind == "pre":
                out.append("<pre>%s</pre>" % html.escape(b[1]))
            elif kind == "figure":
                out.append("<figure><img alt='%s' src='data:image/png;base64,%s'>"
                           "<figcaption class='note'>%s</figcaption></figure>" % (
                               html.escape(b[1]), base64.b64encode(b[2]).decode(),
                               html.escape(b[1])))
            elif kind == "figure_file":
                out.append("<figure><img alt='%s' src='%s'><figcaption class='note'>%s"
                           "</figcaption></figure>" % (html.escape(b[1]), html.escape(b[2]),
                                                       html.escape(b[1])))
    out.append("</body></html>\n")
    return "".join(out)


def _md(title, stamp, sections):
    out = ["# %s" % title, "", "GridPACK report, %s" % stamp, ""]
    for heading, blocks in sections:
        out += ["## %s" % heading, ""]
        for b in blocks:
            kind = b[0]
            if kind in ("text", "note"):
                out += [b[1], ""]
            elif kind == "verdict":
                out += ["> " + l for l in b[1].splitlines()] + [""]
            elif kind == "kv":
                out += ["| | |", "|---|---|"] + ["| %s | %s |" % (k, _fmt(v))
                                                 for k, v in b[1].items()] + [""]
            elif kind == "table":
                cols, rows = b[1], b[2]
                out += ["| " + " | ".join(str(c) for c in cols) + " |",
                        "|" + "---|" * len(cols)]
                out += ["| " + " | ".join(_fmt(v) for v in r) + " |" for r in rows[:_ROWS]]
                if len(rows) > _ROWS:
                    out.append("\n%d of %d rows shown." % (_ROWS, len(rows)))
                out.append("")
            elif kind == "pre":
                out += ["```", b[1], "```", ""]
            elif kind in ("figure", "figure_file"):
                name = b[2] if kind == "figure_file" else "(embedded figure)"
                out += ["![%s](%s)" % (b[1], name), "", "*%s*" % b[1], ""]
    return "\n".join(out)
