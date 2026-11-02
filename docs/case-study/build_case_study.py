#!/usr/bin/env python3
"""Build the one-page case study from CASE_STUDY.md: figure (SVG), HTML, then PDF through headless Chrome.

Usage: python3 docs/case-study/build_case_study.py
The figure's numbers are read from the frozen real-data summary, which is never written to.
"""
import csv, html, pathlib, re, subprocess, sys, tempfile

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
SUMMARY = ROOT / "validation/real-data-study/summary.csv"
BLUE, ORANGE, MUTED, INK, INK2, HAIR = "#2a78d6", "#eb6834", "#898781", "#14242f", "#46555f", "#d9dcd9"
FONT = "'Helvetica Neue', Helvetica, Arial, sans-serif"


def rmse(family, method):
    rows = [r for r in csv.DictReader(SUMMARY.open())
            if r["family"] == family and r["scenario"] == "audit_25pct" and r["method"] == method]
    assert len(rows) == 1
    return float(rows[0]["rmse"])


def figure_svg():
    """Horizontal bars on one shared scale, directly labeled. Width is the page column at 96 px per inch."""
    rows = [("Discard recorded value", "reference_normal", MUTED, False), ("New model", "contamination", ORANGE, True),
            ("Simple linear MI", "cole_mi", BLUE, True), ("All references", "oracle", MUTED, False)]
    panels = [("ClinicalTrials.gov", "ctgov"), ("NHANES", "nhanes")]
    width, label_w, panel_w, top, step, bar_h, r = 352, 112, 120, 15, 15, 8, 4
    scale = 84 / 0.08
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {top + 4 * step}" width="{width}" height="{top + 4 * step}" '
           f'role="img" aria-label="Typical error at a 25% audit, by method and dataset" font-family="{FONT}" font-size="10">']
    for i, (label, method, color, strong) in enumerate(rows):
        y = top + i * step
        weight = ' font-weight="700"' if strong else ""
        out.append(f'<text x="{label_w - 6}" y="{y + 9}" text-anchor="end" fill="{INK if strong else INK2}"{weight}>{label}</text>')
    for p, (title, family) in enumerate(panels):
        x0 = label_w + p * panel_w
        out.append(f'<text x="{x0}" y="10" fill="{INK}" font-weight="700">{title}</text>')
        out.append(f'<line x1="{x0}" x2="{x0}" y1="{top - 1}" y2="{top + 4 * step - 3}" stroke="{HAIR}"/>')
        for i, (_, method, color, strong) in enumerate(rows):
            value, y = rmse(family, method), top + i * step
            length = value * scale
            out.append(f'<path d="M{x0} {y + 1}h{length - r:.1f}a{r} {r} 0 0 1 {r} {r}a{r} {r} 0 0 1 {-r} {r}h{-(length - r):.1f}z" fill="{color}"/>')
            out.append(f'<text x="{x0 + length + 4:.1f}" y="{y + 9}" fill="{INK}"{" font-weight=\"700\"" if strong else ""}>{value:.4f}</text>')
    out.append("</svg>")
    return "\n".join(out)


def inline(text):
    """The small Markdown subset CASE_STUDY.md uses. Link targets are dropped, because print has no repository URL."""
    text = html.escape(text, quote=False)
    text = re.sub(r"\[`([^`]+)`\]\([^)]+\)", r"<code>\1</code>", text)
    text = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", text)
    text = re.sub(r"`([^`]+)`", r"<code>\1</code>", text)
    text = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", text)
    return text


def parse(md):
    title = re.match(r"# (.+)\n", md).group(1)
    intro = re.search(r"\n\n(\*\*A case study.+?)\n\n", md, re.S).group(1)
    table = [[c.strip() for c in line.strip("|").split("|")] for line in re.findall(r"^\|.+\|$", md, re.M)][2:]
    sections = {m.group(1): m.group(2).strip() for m in re.finditer(r"^## (.+?)\n\n(.+?)(?=^## |\Z)", md, re.M | re.S)}
    return title, intro, table, sections


def blocks(body):
    out, items = [], []
    for block in body.split("\n\n"):
        block = block.strip()
        if block.startswith("!["):
            out.append("@@FIGURE@@")
        elif block.startswith("- "):
            items = [inline(line[2:]) for line in block.split("\n")]
            out.append("<ul>" + "".join(f"<li>{i}</li>" for i in items) + "</ul>")
        elif block:
            out.append(f"<p>{inline(block)}</p>")
    return out


CSS = """
@page { size: 8.5in 11in; margin: 0; }
:root { --ink: #14242f; --ink2: #46555f; --navy: #10222f; --card: #1b3244; --line: #2d4a5f; --on-navy: #f2f6f8; --on-navy2: #b7c7d2;
        --sheet: #ffffff; --tint: #eef2f5; --hair: #d9dcd9; --accent: #2a78d6; --warn: #fab219; --crit: #d03b3b; }
* { box-sizing: border-box; -webkit-print-color-adjust: exact; print-color-adjust: exact; }
html, body { margin: 0; padding: 0; background: var(--sheet); color: var(--ink); font: 8.75pt/1.36 'Helvetica Neue', Helvetica, Arial, sans-serif; }
.page { width: 8.5in; }
header { background: var(--navy); color: var(--on-navy); padding: .27in .45in .2in; }
h1 { font: 600 20.5pt/1.12 'Avenir Next', 'Helvetica Neue', sans-serif; margin: .03in 0 .05in; letter-spacing: -.005em; }
.intro { color: var(--on-navy2); font-size: 8.8pt; line-height: 1.34; margin: 0; max-width: 7.3in; }
.intro strong { color: var(--on-navy); font-weight: 600; }
.strip { display: grid; grid-template-columns: repeat(3, 1fr); gap: .08in; margin-top: .13in; }
.tile { background: var(--card); border: 1px solid var(--line); border-radius: 5px; padding: .075in .09in .085in; }
.tile b { display: block; font-weight: 500; font-size: 7.7pt; line-height: 1.22; color: var(--on-navy); }
.chip { display: inline-block; font: 700 7.5pt/1 'Avenir Next', 'Helvetica Neue', sans-serif; padding: .028in .06in; border-radius: 3px; margin: .035in 0 .045in; }
.chip.good { background: var(--accent); color: #fff; } .chip.mixed { background: var(--warn); color: #1c1402; } .chip.stop { background: var(--crit); color: #fff; }
.tile span.note { display: block; font-size: 7.4pt; line-height: 1.26; color: var(--on-navy2); }
main { column-count: 2; column-gap: .27in; padding: .15in .45in .07in; }
h2 { font: 600 10.2pt/1.15 'Avenir Next', 'Helvetica Neue', sans-serif; color: var(--ink); margin: .105in 0 .035in; padding-bottom: .02in; border-bottom: 1px solid var(--hair); break-after: avoid; }
main section:first-child h2 { margin-top: 0; }
figure, .sidebar, li { break-inside: avoid; }
p { margin: 0 0 .055in; }
strong { font-weight: 600; }
ul { margin: 0 0 .05in; padding-left: .135in; } li { margin-bottom: .032in; } li::marker { color: var(--ink2); }
figure { margin: .02in 0 .03in; } figure svg { display: block; width: 100%; height: auto; }
figcaption { font-size: 7.7pt; line-height: 1.28; color: var(--ink2); margin-top: .03in; } figcaption strong { color: var(--ink); }
.sidebar { background: var(--tint); border-radius: 5px; padding: .08in .1in .05in; margin-top: .09in; }
.sidebar h2 { border: 0; margin: 0 0 .02in; font-size: 8.9pt; padding: 0; }
.sidebar p { font-size: 7.75pt; line-height: 1.29; margin: 0 0 .03in; }
footer { margin: 0 .45in; padding: .07in 0 .14in; border-top: 1px solid var(--hair); font-size: 7.45pt; line-height: 1.29; color: var(--ink2); }
footer strong { color: var(--ink); } footer p { margin: 0 0 .03in; }
code { font: 6.7pt 'SF Mono', Menlo, monospace; }
"""


def build_html(md):
    title, intro, table, sections = parse(md)
    chips = {"Model looked strong": "good", "Model lost or tied": "mixed", "Inconclusive": "mixed", "PARTIAL": "mixed", "NO-GO": "stop"}
    tiles = "".join(f'<div class="tile"><b>{inline(t)}</b><span class="chip {chips[r]}">{inline(r)}</span><span class="note">{inline(n)}</span></div>'
                    for t, r, n in table)
    def section(name, cls=""):
        parts = "".join(blocks(sections[name]))
        return f'<section class="{cls}"><h2>{inline(name)}</h2>{parts}</section>'
    evidence = blocks(sections["What the evidence showed"])
    figure_index = evidence.index("@@FIGURE@@")
    caption = evidence[figure_index + 1].replace("<p>", "").replace("</p>", "")
    figure = f"<figure>{figure_svg()}<figcaption>{caption}</figcaption></figure>"
    evidence_html = "".join(evidence[:figure_index])
    limits = blocks(sections["Limits and AI disclosure"])
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>catmech case study</title><style>{CSS}</style></head>
<body><div class="page">
<header><h1>{inline(title)}</h1><p class="intro">{inline(intro)}</p>
<div class="strip">{tiles}</div></header>
<main>
{section("The problem")}{section("What was built")}<section><h2>What the evidence showed</h2>{evidence_html}</section>
<section><h2>Typical error, simple method against new model</h2>{figure}</section>{section("The decision")}{section("Technical sidebar", "sidebar")}
</main>
<footer><p><strong>Limits and AI disclosure.</strong> {limits[0].replace("<p>", "").replace("</p>", "")}</p>{limits[1]}</footer>
</div></body></html>
"""


def main():
    md = (HERE / "CASE_STUDY.md").read_text()
    (HERE / "case-study-figure.svg").write_text(figure_svg() + "\n")
    html_path, pdf_path = HERE / "catmech-case-study.html", HERE / "catmech-case-study.pdf"
    html_path.write_text(build_html(md))
    with tempfile.TemporaryDirectory() as profile:
        subprocess.run([CHROME, "--headless=new", "--disable-gpu", "--no-pdf-header-footer", f"--user-data-dir={profile}",
                        "--no-first-run", "--disable-extensions", f"--print-to-pdf={pdf_path}", html_path.as_uri()],
                       check=True, capture_output=True, timeout=120)
    info = subprocess.run(["pdfinfo", str(pdf_path)], check=True, capture_output=True, text=True).stdout
    pages = int(re.search(r"Pages:\s+(\d+)", info).group(1))
    size = re.search(r"Page size:\s+(.+)", info).group(1)
    print(f"{pdf_path.name}: {pages} page(s), {size}")
    if pages != 1:
        sys.exit("The case study must be exactly one page.")


if __name__ == "__main__":
    main()
