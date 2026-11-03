"""Draws docs/figures/how-it-works.svg: five equal cards, captions in short lines that fit the cards.
Usage: python3 -I tools/how-it-works-figure.py docs/figures/how-it-works.svg
Caption lines stay short so they fit with room to spare in wider fonts too. After a change, regenerate the
guide's inline copy (docs/learn/index.html, between the HOW-FIGURE comments); the verifier compares them."""
import sys
from pathlib import Path

W, H, CARD, GAP, M = 1000, 330, 180, 22, 8
xs = [M + i * (CARD + GAP) for i in range(5)]
out = []
add = out.append
add(f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" role="img" aria-labelledby="t d" font-family="-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif">
  <title id="t">How catmech works</title>
  <desc id="d">Five panels from left to right. One: a table where every row has a recorded value and only audited rows have a checked value. Two: on audited rows, a regression of the checked value on the recorded value and the other variables. Three: each unaudited row's checked value is drawn 20 times from that model. Four: the analysis runs on each of the 20 completed copies. Five: Rubin's rules combine them into one estimate whose interval includes the spread between copies.</desc>
  <style>
    .bg {{ fill: #ffffff; }} .card {{ fill: #f4f6f9; stroke: #d7dee7; }}
    .ink {{ fill: #121a24; }} .muted {{ fill: #4d5a69; }} .cell {{ fill: #c3ccd8; }} .empty {{ fill: none; stroke: #9aa6b5; stroke-dasharray: 3 2; }}
    .check {{ fill: #0b7268; }} .line {{ stroke: #0b7268; }} .axis {{ stroke: #9aa6b5; }} .arrow {{ stroke: #7d8a9b; fill: none; }} .head {{ fill: #7d8a9b; }}
    .copy {{ fill: #ffffff; stroke: #c3ccd8; }} .bar {{ stroke: #4d5a69; }} .dot {{ fill: #4d5a69; }} .final {{ stroke: #0b7268; }} .finaldot {{ fill: #0b7268; }}
    @media (prefers-color-scheme: dark) {{
      .bg {{ fill: #0d1117; }} .card {{ fill: #161e2a; stroke: #2e3a4d; }}
      .ink {{ fill: #e8eef6; }} .muted {{ fill: #98a5b8; }} .cell {{ fill: #3a475b; }} .empty {{ stroke: #62718a; }}
      .check {{ fill: #3ad7c8; }} .line {{ stroke: #3ad7c8; }} .axis {{ stroke: #62718a; }} .arrow {{ stroke: #8593a8; }} .head {{ fill: #8593a8; }}
      .copy {{ fill: #111721; stroke: #3a475b; }} .bar {{ stroke: #98a5b8; }} .dot {{ fill: #98a5b8; }} .final {{ stroke: #3ad7c8; }} .finaldot {{ fill: #3ad7c8; }}
    }}
    text {{ font-size: 13px; }} .h {{ font-size: 14px; font-weight: 600; }} .s {{ font-size: 12px; }}
  </style>
  <rect class="bg" x="0" y="0" width="{W}" height="{H}"/>''')

def card(i, title, caption):
    x = xs[i]
    add(f'\n  <rect class="card" x="{x}" y="8" width="{CARD}" height="{H - 16}" rx="8"/>')
    add(f'  <text class="ink h" x="{x + 14}" y="32">{title}</text>')
    for k, line in enumerate(caption):
        add(f'  <text class="muted s" x="{x + 14}" y="{248 + 16 * k}">{line}</text>')
    if i < 4:
        a = x + CARD + 3
        add(f'  <path class="arrow" d="M{a} 140 H {a + 9}" stroke-width="2"/><path class="head" d="M{a + 9} 135 L{a + 16} 140 L{a + 9} 145 Z"/>')

# 1. Records
x = xs[0]
card(0, "1. Your records", ["Every row has a", "recorded value.", "Audited rows also", "have a checked one."])
add(f'  <text class="muted s" x="{x + 24}" y="58">recorded</text>')
add(f'  <text class="muted s" x="{x + 100}" y="58">checked</text>')
for k, audited in enumerate([0, 1, 0, 0, 1, 0, 0, 1]):
    y = 68 + 20 * k
    add(f'  <rect class="cell" x="{x + 24}" y="{y}" width="60" height="14" rx="2"/><rect class="{"check" if audited else "empty"}" x="{x + 100}" y="{y}" width="60" height="14" rx="2"/>')

# 2. Learn
x = xs[1]
card(1, "2. Learn", ["On audited rows, fit", "checked on recorded", "and other variables."])
add(f'  <line class="axis" x1="{x + 30}" y1="208" x2="{x + 164}" y2="208"/><line class="axis" x1="{x + 30}" y1="208" x2="{x + 30}" y2="56"/>')
add(f'  <text class="muted s" x="{x + 72}" y="226">recorded</text>')
add(f'  <text class="muted s" x="0" y="0" transform="translate({x + 22} 166) rotate(-90)">checked</text>')
add(f'  <line class="line" x1="{x + 38}" y1="196" x2="{x + 156}" y2="70" stroke-width="2"/>')
for px, py in [(46, 182), (60, 176), (74, 150), (90, 144), (106, 116), (120, 112), (136, 88), (98, 132)]:
    add(f'  <circle class="check" cx="{x + px}" cy="{py}" r="4"/>')

# 3. Fill in
x = xs[2]
card(2, "3. Fill in 20 ways", ["Draw each missing", "checked value from", "the model, with its", "uncertainty."])
for dx, dy in [(46, 50), (38, 58), (30, 66)]:
    add(f'  <rect class="copy" x="{x + dx}" y="{dy}" width="100" height="132" rx="4"/>')
for k, solid in enumerate([0, 1, 0, 0, 1, 0]):
    y = 78 + 18 * k
    add(f'  <rect class="cell" x="{x + 42}" y="{y}" width="36" height="12" rx="2"/><rect class="check" x="{x + 86}" y="{y}" width="36" height="12" rx="2"{"" if solid else " opacity=\".45\""}/>')
add(f'  <text class="ink h" x="{x + 136}" y="218">×20</text>')

# 4. Analyze each
x = xs[3]
card(3, "4. Analyze each", ["Run your analysis", "on every completed", "copy."])
for k, (lo, hi, mid) in enumerate([(48, 98, 72), (56, 108, 82), (40, 92, 66), (52, 102, 77), (46, 96, 71), (60, 112, 86), (44, 94, 69), (50, 102, 76)]):
    y = 62 + 18 * k
    add(f'  <line class="bar" x1="{x + lo + 10}" y1="{y}" x2="{x + hi + 10}" y2="{y}" stroke-width="1.5"/><circle class="dot" cx="{x + mid + 10}" cy="{y}" r="3"/>')

# 5. Combine
x = xs[4]
card(4, "5. Combine", ["One estimate. Its", "interval adds how", "much the copies", "disagree (Rubin's", "rules)."])
add(f'  <line class="final" x1="{x + 46}" y1="132" x2="{x + 134}" y2="132" stroke-width="3"/>')
add(f'  <line class="final" x1="{x + 46}" y1="124" x2="{x + 46}" y2="140" stroke-width="2"/><line class="final" x1="{x + 134}" y1="124" x2="{x + 134}" y2="140" stroke-width="2"/>')
add(f'  <circle class="finaldot" cx="{x + 90}" cy="132" r="6"/>')
add("</svg>\n")
Path(sys.argv[1]).write_text("\n".join(out), encoding="utf-8")
print("wrote", sys.argv[1])
