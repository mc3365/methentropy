"""Draw docs/figures/methentropy_workflow.svg (render to PNG with rsvg-convert).

Panel colours follow the original manual schematic of the v1 pipeline:
blue = raw data, green = simulation / references, peach = results,
pink = visualisation.
"""

from pathlib import Path

import numpy as np

W, H = 1600, 1010
FONT = "Inter, Helvetica, Arial, sans-serif"
INK, INK2, LINE = "#111111", "#444444", "#333333"
BLUE, GREEN, PEACH, PINK, GREY = "#dbeaf6", "#dcf0d8", "#fbe7d9", "#f2def0", "#f4f3f0"
out = []


def rect(x, y, w, h, fill="white", stroke=LINE, sw=1.6, rx=10, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    out.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" '
               f'stroke="{stroke}" stroke-width="{sw}"{d}/>')


def text(x, y, lines, size=17, weight=400, anchor="middle", color=INK, lh=1.3):
    if isinstance(lines, str):
        lines = [lines]
    for i, ln in enumerate(lines):
        ln = ln.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        out.append(f'<text x="{x}" y="{y + i * size * lh:.1f}" font-family="{FONT}" '
                   f'font-size="{size}" font-weight="{weight}" fill="{color}" '
                   f'text-anchor="{anchor}">{ln}</text>')


def arrow(x1, y1, x2, y2, color=LINE):
    out.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" '
               f'stroke-width="2" marker-end="url(#ah)"/>')


def path_arrow(d):
    out.append(f'<path d="{d}" fill="none" stroke="{LINE}" stroke-width="2" '
               f'marker-end="url(#ah)"/>')


def grid(x, y, M, r=5.5, gap=15, label=None):
    """Reads (rows) x CpGs (columns); filled = methylated."""
    for i, row in enumerate(M):
        out.append(f'<line x1="{x - 4}" y1="{y + i * gap}" x2="{x + 3 * gap + 4}" '
                   f'y2="{y + i * gap}" stroke="#9a9993" stroke-width="1"/>')
        for j, v in enumerate(row):
            fill = INK if v else "white"
            out.append(f'<circle cx="{x + j * gap}" cy="{y + i * gap}" r="{r}" fill="{fill}" '
                       f'stroke="{INK}" stroke-width="1.3"/>')
    if label:
        text(x + 1.5 * gap, y + len(M) * gap + 8, label, size=13, color=INK2)


def panel(x, y, w, h, fill, num, title):
    rect(x, y, w, h, fill=fill, stroke="none", rx=14)
    text(x + 18, y + 32, f"{num}  {title}", size=19, weight=600, anchor="start")


out.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
           f'viewBox="0 0 {W} {H}">')
out.append('<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" '
           f'markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" '
           f'fill="{LINE}"/></marker></defs>')
rect(0, 0, W, H, fill="white", stroke="none", rx=0)
text(30, 44, "methentropy — is read-level methylation more ordered than chance?", size=24,
     weight=600, anchor="start")
text(30, 72, "Every reference is computed at the window's own read depth, so finite-sample "
     "entropy bias is shared by observation and null.", size=15, color=INK2, anchor="start")

# ---------------------------------------------------------------- 1 input
panel(30, 100, 350, 470, BLUE, "1", "Input")
for i, (lines, yy) in enumerate([(["WGBS BAM", "(sorted, deduplicated)"], 160),
                                 (["RLM, entropy mode", "(Hetzel et al.)"], 285),
                                 (["Per 4-CpG window", "16 epiallele counts,", "reads N, methylation",
                                   "of each CpG"], 410)]):
    h = 34 + 22 * len(lines)
    rect(65, yy, 280, h)
    text(205, yy + 30, lines, size=16)
for y1, y2 in [(238, 283), (363, 408)]:
    arrow(205, y1, 205, y2)

# ------------------------------------------------------------ 2 references
panel(420, 100, 1150, 470, GREEN, "2", "References at the window's own depth N")
bx = [445, 820, 1195]
titles = ["Minimum (exact)", "Bimodal / allelic", "Independence null"]
bodies = [["Split the N reads into classes", "n0 … n4 (j methylated CpGs),",
           "Σ j·nj = k methylated calls;", "lowest attainable entropy"],
          ["Each read fully methylated", "or fully unmethylated at the",
           "window's level m: h(m) / 4,", "two-population expectation"],
          ["Keep each CpG's methylated", "count, shuffle reads per CpG;",
           "1,000 draws, cached per", "(N, counts): mean, SD, p"]]
for x, t, b in zip(bx, titles, bodies):
    rect(x, 150, 350, 395)
    text(x + 175, 182, t, size=18, weight=600)
    text(x + 175, 214, b, size=15, color=INK2)
# mini illustrations: 8 reads x 4 CpGs
rng = np.random.default_rng(4)
same = [[1, 1, 0, 0]] * 4 + [[0, 0, 0, 0]] * 4
grid(bx[0] + 152, 345, same, label="identical patterns")
bim = [[1, 1, 1, 1]] * 4 + [[0, 0, 0, 0]] * 4
grid(bx[1] + 152, 345, bim, label="all-or-none reads")
obs = [[1, 1, 1, 1], [1, 1, 1, 0], [1, 1, 1, 1], [0, 0, 0, 0],
       [0, 0, 0, 1], [1, 1, 1, 1], [0, 0, 0, 0], [0, 0, 0, 0]]
cols = np.array(obs).T
shuf = np.array([rng.permutation(c) for c in cols]).T.tolist()
grid(bx[2] + 70, 345, obs, label="observed")
grid(bx[2] + 235, 345, shuf, label="shuffled per CpG")
arrow(bx[2] + 140, 397, bx[2] + 218, 397)
text(bx[2] + 175, 522, "below null = CpGs co-methylated on molecules", size=13, color=INK2)
text(bx[0] + 175, 522, "0 at 25 / 50 / 75 % methylation", size=13, color=INK2)
text(bx[1] + 175, 522, "0.25 at 50 % methylation", size=13, color=INK2)
path_arrow("M 345 470 L 395 470 L 395 345 L 418 345")

# -------------------------------------------------------------- 3 results
panel(420, 600, 1150, 165, PEACH, "3", "Per-window statistics  (methentropy run)")
stats_items = [("Observed entropy H", "from the epiallele counts"),
               ("Δ = H − null mean", "< 0: more ordered than chance"),
               ("Scaled = (H − min) / (null − min)", "0 = minimum, 1 = chance"),
               ("z,  p_low / p_high,  BH q", "Monte Carlo, genome-wide FDR")]
for i, (a, b) in enumerate(stats_items):
    x = 445 + i * 280
    rect(x, 645, 265, 100, fill="white")
    text(x + 132, 685, a, size=15, weight=600)
    text(x + 132, 712, b, size=13.5, color=INK2)
arrow(995, 547, 995, 598)

# calibration side box
panel(30, 600, 350, 165, GREY, "", "Calibration  (methentropy calibrate)")
text(205, 670, ["Linkage-free synthetic reads with", "each window's real N and CpG levels:",
                "expect ≤ 5 % at p < 0.05"], size=15, color=INK2)
arrow(205, 533, 205, 598)

# -------------------------------------------------------------- 4 outputs
panel(30, 795, 1540, 195, PINK, "4", "Region-level summaries and figures")
outs = [("metaplot", ["Regions scaled to 47 bins", "+ 10 flanking bins per side;",
                      "median and IQR across regions"]),
        ("compare", ["Per region, per replicate;", "Wilcoxon over regions, Welch",
                     "over replicates, BH per region"]),
        ("qc", ["Entropy vs methylation, Δ", "distribution, depth dependence,",
                "p-value histogram"]),
        ("aggregate", ["Mean Δ, coverage-weighted Δ", "and fraction of significant",
                       "windows per BED region"])]
for i, (name, body) in enumerate(outs):
    x = 60 + i * 380
    rect(x, 840, 350, 130)
    text(x + 175, 872, f"methentropy {name}", size=17, weight=600)
    text(x + 175, 900, body, size=14, color=INK2)
arrow(995, 767, 995, 838)

out.append("</svg>")
p = Path(__file__).with_name("methentropy_workflow.svg")
p.write_text("\n".join(out))
print(p)
