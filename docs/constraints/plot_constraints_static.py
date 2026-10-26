"""
plot_constraints_static.py
==========================

Render the structured constraint YAML as a compact static figure (PNG *and*
vector PDF) suitable for embedding in README.md, the docs site, and papers.
The PDF keeps every string searchable / copy-able.

Layout (top → bottom)
---------------------
Rows are grouped into two blocks:

* **Include bounds — allowed / favoured regions** (green): ``m > x`` /
  ``m < x`` / range bounds marking where wave DM is still allowed, including
  the 1σ point estimates of the literature (they are *favoured* values, so
  they sit with this block);
* **Exclude bounds — ruled-out regions** (orange, hatched).

Within a tone every *dataset family* (rotation curves, dwarf galaxies, Lyα
forest, 21 cm, CMB, lensing, black holes, …) gets its own shade, and a dashed
horizontal guide runs through every row so a label can be traced across the
figure.

Inside each block the rows follow the order of the source table
``ma_constraints.tex`` (the ``table_index`` field of the YAML), so the figure
can be read side by side with the paper table.  Every *paper* is one row: a
work quoting several intervals (e.g. Bar et al. 2019 with ``m ≲ 4×10⁻²² eV``
and ``2×10⁻²⁰ ≲ m ≲ 8×10⁻¹⁹ eV``) draws several bars on its shared row instead
of repeating the label three times.

Row label: ``FirstAuthor [et al.] Year``; two *different* papers that
truncate to the same label get ``a`` / ``b`` suffixes (Bar et al. 2019a/b,
Maleki et al. 2020a/b).

Annotation: ``data — method`` using the abbreviations of
``constraint_labels.py`` (``rotation curve`` → ``RC``, ``soliton–halo
relation`` → ``SH relation``, …), with the method half in italics.  The width
of every annotation is measured with ``TextPath`` *before* drawing and the
text is placed in the first slot of its own row that is wide enough:

1. **inside a band that has room for it** (widest bar first) — the ``m > x`` /
   ``m ≤ x`` bands on either side of the figure and the wide interior ranges
   (Bar et al. 2022) alike;
2. **hugging the right end of the bars**;
3. in the free **gap between two bars** of the same row;
4. **hugging the left end of the bars** (right-aligned, so the label stays
   glued to the bar instead of floating at the axis edge).

Nothing can spill past the axes (the first version clipped ``Bar et al.
2022`` and ``Bar et al. 2019`` at the left edge); text over a hatched fill
gets a translucent plate.

Bars whose bound is open (``m > x`` runs to the right edge, ``m < x`` to the
left edge) or that are cut by the axis range carry a small arrowhead at the
open end.  A bound whose quoted value is itself a range (Marsh & Niemeyer
2019: ``m ≳ 0.6–1 ×10⁻¹⁹ eV``) also carries a white hatched *uncertainty
band* on the bar.

Outputs
-------
* ``output/ma_constraints.png`` — 220 DPI raster (README / docs preview).
* ``output/ma_constraints.pdf`` — searchable vector PDF (for papers).
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patches as mpatches  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.textpath import TextPath  # noqa: E402
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter  # noqa: E402
import yaml  # noqa: E402

from constraint_labels import assign_row_labels, data_label, method_label  # noqa: E402

HERE = Path(__file__).resolve().parent
YAML_IN = HERE / "output" / "ma_constraints.yaml"
OUT_PNG = HERE / "output" / "ma_constraints.png"
OUT_PDF = HERE / "output" / "ma_constraints.pdf"

# Axis range and layout constants (y units = one row slot).
X_MIN, X_MAX = -26.0, -18.0
NDEC = X_MAX - X_MIN
FIG_W = 19.0          # inches
IN_PER_SLOT = 0.40    # vertical inches per row
# Font sizes.  The figure geometry above is fixed, so these are the single
# knob for "how large is the text relative to the plot".  A row whose
# annotation no longer fits its slot is shrunk by ``_place_annotation``.
ANNOT_FS = 17.5       # data — method annotation
ROW_LABEL_FS = 18.0   # left-hand row labels
TICK_FS = 21.0        # x-axis decade ticks
AXLABEL_FS = 25.0
TITLE_FS = 30.0
BLOCK_FS = 21.0
LEGEND_FS = 19.0
ROW_H, HEADER_H = 1.0, 0.85
# Visible gap between the two block bands: the bands are tinted, so what the
# reader sees is the untinted strip in between.  ``0.3`` slot units of layout
# spacing leaves ~1 mm of white between the include band and the exclude band.
BLOCK_GAP = 0.3
BAR_H = 0.80
ARROW_W = 0.085       # decades
PAD = 0.05            # decades of breathing room around bars
# Minimum text-free margin from an arrowhead (data-coord decades).  Text laid
# out inside a band next to an arrow would otherwise overlap the triangle.
PAD_ARROW = ARROW_W + 0.025

# Two tones carry the meaning: green = allowed / favoured, orange (hatched) =
# ruled out.  Within each tone every *dataset family* gets its own shade, so
# rows that share a dataset read as a group (the shades are decorative — the
# legend explains the tones, not the shades).
TONE = {
    "include": dict(fill="#a9d18e", edge="#3d6b1f", hatch=None),
    "exclude": dict(fill="#f6b26b", edge="#8a4b12", hatch="///"),
}
FAMILY_TINTS: list[tuple[tuple[str, ...], str, str]] = [
    # (substrings of the data label, include shade, exclude shade)
    (("m87", "bh spin"),                 "#a3cfa8", "#f2a95f"),   # black holes
    (("lens", "cozmic", "satellites", "subhalo"),
                                         "#c3d79b", "#f7c58a"),   # (sub)halo / lensing
    (("lyα", "lyman"),                   "#bfd98d", "#f5c07c"),   # Lyα forest
    (("21cm",),                          "#96cd9d", "#f0b46b"),   # 21 cm
    (("cmb",),                           "#cdd88f", "#f4c78e"),   # CMB / LSS
    (("pulsar",),                        "#b6d09a", "#efb87f"),   # pulsar timing
    (("dsph", "ufd", "segue", "eri ii", "fornax", "carina", "leo ii",
      "stellar kin", "gc timing", "stellar dynamics"),
                                         "#a9d6a0", "#f1bb72"),   # dwarf galaxies
    (("mw", "stellar streams"),          "#9dce8d", "#f4ac66"),   # Milky Way
    (("sparc", "m31", "ngc", "spiral", "rc"),
                                         "#a9d18e", "#f6b26b"),   # rotation curves
]
BAND_STYLE = dict(fill="#ffffff", edge="#3d6b1f", hatch="\\\\\\\\")
TEXT_DARK = "#1b1b1b"

BLOCKS = [
    ("include", "Include bounds — allowed / favoured regions", "#3d6b1f", "#f2f8ec"),
    ("exclude", "Exclude bounds — ruled-out regions", "#8a4b12", "#fdf3ea"),
]


# --------------------------------------------------------------------------- #
# Data helpers                                                                #
# --------------------------------------------------------------------------- #

def _load() -> dict:
    if not YAML_IN.exists():
        raise FileNotFoundError(f"Run ma_constraints_data.py first; missing {YAML_IN}")
    with open(YAML_IN, encoding="utf-8") as f:
        return yaml.safe_load(f)


def _tone_of(iv: dict) -> str:
    """Which of the two blocks/tone an interval belongs to."""
    return "exclude" if iv["type"].startswith("exclude") else "include"


def _bar_style(iv: dict, data_text: str) -> dict:
    """Bar colours: the tone of the interval type, shaded by dataset family."""
    tone = _tone_of(iv)
    style = dict(TONE[tone])
    label = (data_text or "").lower()
    for keys, include_fill, exclude_fill in FAMILY_TINTS:
        if any(k in label for k in keys):
            style["fill"] = include_fill if tone == "include" else exclude_fill
            break
    return style


def _is_open_lo(v) -> bool:
    return v is None or (isinstance(v, float) and math.isinf(v) and v < 0)


def _is_open_hi(v) -> bool:
    return v is None or (isinstance(v, float) and math.isinf(v) and v > 0)


def _span(iv: dict) -> tuple[float, float]:
    """Interval clipped to the drawn range, in log10(m_a / eV)."""
    lo, hi = iv["log10_min"], iv["log10_max"]
    left = X_MIN if _is_open_lo(lo) else max(float(lo), X_MIN)
    right = X_MAX if _is_open_hi(hi) else min(float(hi), X_MAX)
    return left, right


def _in_view(iv: dict) -> bool:
    left, right = _span(iv)
    return right > left


def _extends_left(iv: dict) -> bool:
    return _is_open_lo(iv["log10_min"]) or float(iv["log10_min"]) < X_MIN


def _extends_right(iv: dict) -> bool:
    return _is_open_hi(iv["log10_max"]) or float(iv["log10_max"]) > X_MAX


# --------------------------------------------------------------------------- #
# Text measurement + annotation placement                                      #
# --------------------------------------------------------------------------- #

_WIDTH_CACHE: dict[tuple[str, float], float] = {}


def _text_width_in(text: str, fontsize: float) -> float:
    """Rendered width of ``text`` in inches (mathtext aware, 3% head-room)."""
    key = (text, fontsize)
    if key not in _WIDTH_CACHE:
        try:
            ext = TextPath((0, 0), text, size=fontsize, usetex=False).get_extents()
            width = ext.width / 72.0
        except Exception:  # pragma: no cover - font metrics fallback
            width = 0.52 * fontsize * len(text) / 72.0
        _WIDTH_CACHE[key] = width * 1.03
    return _WIDTH_CACHE[key]


def _annotation(data_txt: str, method_txt: str) -> str:
    """``data — method`` as one string; the method half is mathtext-italic."""
    if data_txt and method_txt:
        esc = method_txt.replace(" ", r"\ ")
        return f"{data_txt}   $\\mathit{{{esc}}}$"
    return data_txt or method_txt


def _annotation_slots(bars: list[tuple[float, float]],
                      arrows: list[tuple[bool, bool]] | None = None) -> list[dict]:
    """Free slots of a row, in decreasing order of preference.

    ``inside``  — inside the bars (widest first).  Whenever the label fits in
                  the band it belongs *in* the band — the ``m > x`` /
                  ``m ≤ x`` bands on either side of the figure, the wide
                  interior ranges (Bar et al. 2022), all of them.  The slot
                  stops *before* an open-end arrowhead, so the label never
                  collides with the triangle;
    ``right``   — hugging the right end of the bars;
    ``middle``  — the free gap between two bars of the same row;
    ``left``    — hugging the *left* end of the bars, right-aligned, so the
                  label stays glued to the bar instead of floating at the
                  axis edge.

    ``arrows`` is one ``(has_left_arrow, has_right_arrow)`` flag per bar.
    Every slot carries the text anchor (``x``) and alignment (``ha``) to use;
    ``x0``/``x1`` only describe the free width the slot offers.
    """
    arrows = list(arrows) if arrows else [(False, False)] * len(bars)
    spans: list[tuple[float, float]] = []
    cur = X_MIN
    for left, right in sorted(bars):
        if left - cur > 0:
            spans.append((cur, left))
        cur = max(cur, right)
    if X_MAX - cur > 0:
        spans.append((cur, X_MAX))

    right_spans = [s for s in spans if abs(s[1] - X_MAX) < 1e-9]
    left_spans = [s for s in spans if abs(s[0] - X_MIN) < 1e-9]
    middle = [s for s in spans if s not in right_spans and s not in left_spans]
    middle.sort(key=lambda s: -(s[1] - s[0]))

    inner = sorted(
        ((l + (PAD_ARROW if al else PAD), r - (PAD_ARROW if ar else PAD))
         for (l, r), (al, ar) in zip(bars, arrows)
         if (r - l) > 2 * PAD),
        key=lambda s: -(s[1] - s[0]),
    )

    def _inside(l: float, r: float) -> dict:
        return dict(where="inside", x0=l, x1=r, x=l, ha="left", inside=True)

    slots = [_inside(l, r) for l, r in inner]
    slots += [dict(where="right", x0=a + PAD, x1=b - 0.02, x=a + PAD, ha="left", inside=False)
              for a, b in right_spans]
    slots += [dict(where="middle", x0=a + PAD, x1=b - PAD, x=a + PAD, ha="left", inside=False)
              for a, b in middle]
    slots += [dict(where="left", x0=a + PAD, x1=b - PAD, x=b - PAD, ha="right", inside=False)
              for a, b in left_spans]
    return slots


def _truncate(text: str, max_in: float, fontsize: float,
              width_fn=None) -> str:
    """Shorten ``text`` with an ellipsis until it fits into ``max_in`` inches."""
    width_fn = width_fn or _text_width_in
    if max_in <= 0:
        return ""
    lo, hi = 0, len(text)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if width_fn(text[:mid] + "…", fontsize) <= max_in:
            lo = mid
        else:
            hi = mid - 1
    return text[:lo] + "…"


def _place_annotation(text: str, bars: list[tuple[float, float]],
                      base_fs: float, in_per_dec: float,
                      width_fn=None,
                      arrows: list[tuple[bool, bool]] | None = None) -> dict | None:
    """Choose the first slot of the row that fits ``text``, never overflowing.

    ``width_fn(text, fontsize)`` returns the rendered width expressed in the
    same unit as ``in_per_dec`` (inches for the matplotlib output, pixels for
    the SVG output) — so the same placement rules serve every renderer.
    """
    width_fn = width_fn or _text_width_in
    slots = _annotation_slots(bars, arrows)
    if not slots:
        return None
    # Try the full size first, then step down in ~7 % increments so a crowded
    # row degrades gently instead of jumping straight to the last resort.
    ladder = [base_fs * f for f in (1.0, 0.93, 0.86, 0.79, 0.72, 0.65)]
    for fs in ladder:
        for slot in slots:
            if (slot["x1"] - slot["x0"]) * in_per_dec >= width_fn(text, fs):
                return dict(text=text, x=slot["x"], ha=slot["ha"], fs=round(fs, 1),
                            inside=slot["inside"], where=slot["where"])
    # Last resort (should not happen): widest slot, shrink, else truncate.
    slot = max(slots, key=lambda s: s["x1"] - s["x0"])
    avail = (slot["x1"] - slot["x0"]) * in_per_dec
    fs, txt = base_fs, text
    if width_fn(txt, fs) > avail:
        fs = max(8.0, round(base_fs * avail / width_fn(txt, fs), 1))
    if width_fn(txt, fs) > avail:
        txt = _truncate(txt, avail, fs, width_fn)
    return dict(text=txt, x=slot["x"], ha=slot["ha"], fs=fs,
                inside=slot["inside"], where=slot["where"] + "*")


# --------------------------------------------------------------------------- #
# Layout                                                                       #
# --------------------------------------------------------------------------- #

def _build_layout(data: dict) -> tuple[list[dict], list[dict], float]:
    """Two blocks (include / exclude), rows in ``ma_constraints.tex`` order."""
    per_paper: list[dict] = []
    for key, row in data.items():
        ivs = [iv for iv in row["constraints"] if _in_view(iv)]
        if not ivs:
            continue
        per_paper.append(dict(citekey=key, row=row, ivs=ivs,
                              block=_tone_of(ivs[0]),
                              table_index=row.get("table_index", 10 ** 6)))

    # ``a`` / ``b`` suffixes are assigned over the canonical (table) order.
    per_paper.sort(key=lambda e: (e["table_index"], e["citekey"]))
    labels = assign_row_labels([e["row"] for e in per_paper])

    blocks: list[dict] = []
    slots: list[dict] = []
    y = 0.0
    for block_key, header, colour, tint in BLOCKS:
        entries = [e for e in per_paper if e["block"] == block_key]
        if not entries:
            continue
        block = dict(key=block_key, header=header, colour=colour, tint=tint,
                     top=y, bottom=None)
        slots.append(dict(kind="header", y=y, block=block))
        y += HEADER_H
        for entry in entries:
            slots.append(dict(kind="row", y=y, entry=entry,
                              label=labels[entry["citekey"]]))
            y += ROW_H
        block["bottom"] = y - 0.10
        y += BLOCK_GAP
        blocks.append(block)
    return slots, blocks, y


# --------------------------------------------------------------------------- #
# Drawing                                                                      #
# --------------------------------------------------------------------------- #

def _draw_arrow(ax, x_tip: float, y: float, side: str, style: dict) -> None:
    """Triangle at an open end: tip on the axis edge, base inside the bar."""
    base = x_tip * 10 ** (-ARROW_W if side == "right" else ARROW_W)
    verts = [(base, y - BAR_H * 0.30), (base, y + BAR_H * 0.30), (x_tip, y)]
    ax.add_patch(mpatches.Polygon(verts, closed=True, facecolor=style["fill"],
                                  edgecolor=style["edge"], linewidth=0.6, zorder=3))


def main() -> None:
    data = _load()
    slots, blocks, total_y = _build_layout(data)
    rows = [s for s in slots if s["kind"] == "row"]

    fig_h = IN_PER_SLOT * total_y + 2.8
    fig, ax = plt.subplots(figsize=(FIG_W, fig_h))
    fig.patch.set_facecolor("white")

    ax.set_xscale("log")
    ax.set_xlim(10 ** X_MIN, 10 ** X_MAX)
    ax.set_ylim(total_y - 0.25, -0.55)
    ax.set_yticks([s["y"] + ROW_H / 2 for s in rows])
    ax.set_yticklabels([s["label"] for s in rows], fontsize=ROW_LABEL_FS)
    # The extra pad leaves room for the block ribbon drawn in the left margin.
    ax.tick_params(axis="y", which="both", length=0, pad=26)

    major_locs = [10 ** p for p in range(int(X_MIN), int(X_MAX) + 1)]
    minor_locs = [k * 10 ** p for p in range(int(X_MIN), int(X_MAX)) for k in range(2, 10)]
    ax.xaxis.set_major_locator(FixedLocator(major_locs))
    ax.xaxis.set_minor_locator(FixedLocator(minor_locs))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.xaxis.set_major_formatter(
        FuncFormatter(lambda v, _p=None: rf"$10^{{{math.log10(v):.0f}}}$" if v > 0 else "")
    )
    ax.tick_params(axis="x", which="major", labelsize=TICK_FS, length=5.5, width=0.9)
    ax.tick_params(axis="x", which="minor", labelsize=TICK_FS, length=3.0, width=0.7)
    ax.set_xlabel(r"$\log_{10}(m_a/\mathrm{eV})$", fontsize=AXLABEL_FS, labelpad=10)
    ax.set_title(r"Observational constraints on the wave-CDM particle mass $m_a$",
                 fontsize=TITLE_FS, pad=26)
    ax.grid(True, axis="x", alpha=0.22, linestyle="--", zorder=0.4)
    # Horizontal guides per row make it easy to trace a row across the figure.
    ax.grid(True, axis="y", alpha=0.20, linestyle="--", zorder=0.4)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#bbbbbb")

    # Block bands + ribbon markers in the left margin.
    for block in blocks:
        ax.axhspan(block["top"] - 0.30, block["bottom"], facecolor=block["tint"],
                   edgecolor="none", zorder=0.5)
        ax.add_patch(mpatches.Rectangle(
            (10 ** (X_MIN - 0.095), block["top"] - 0.22),
            10 ** (X_MIN - 0.020) - 10 ** (X_MIN - 0.095),
            block["bottom"] - block["top"] + 0.12,
            facecolor=block["colour"], edgecolor="none", alpha=0.85,
            clip_on=False, zorder=1))
    for slot in slots:
        if slot["kind"] == "header":
            block = slot["block"]
            ax.text(10 ** (X_MIN + 0.06), slot["y"] + HEADER_H * 0.45,
                    block["header"], fontsize=BLOCK_FS, fontweight="bold",
                    color=block["colour"], ha="left", va="center", zorder=5)

    # ``tight_layout`` fixes the axes box; measure text afterwards so the
    # annotation fitting uses the real inches-per-decade of the figure.
    fig.tight_layout()
    ax_pos = ax.get_position()
    in_per_dec = fig.get_size_inches()[0] * ax_pos.width / NDEC

    report: list[str] = []
    for slot in slots:
        if slot["kind"] != "row":
            continue
        entry, row, y = slot["entry"], slot["entry"]["row"], slot["y"] + ROW_H / 2
        bars: list[tuple[float, float]] = []
        arrows: list[tuple[bool, bool]] = []
        data_text = data_label(row)
        for iv in entry["ivs"]:
            left, right = _span(iv)
            if right <= left:
                continue
            bars.append((left, right))
            arrows.append((_extends_left(iv), _extends_right(iv)))
            style = _bar_style(iv, data_text)
            ax.add_patch(mpatches.Rectangle(
                (10 ** left, y - BAR_H / 2), 10 ** right - 10 ** left, BAR_H,
                facecolor=style["fill"], edgecolor=style["edge"],
                hatch=style["hatch"], linewidth=0.6, alpha=0.95, zorder=2))
            band_hi = iv.get("log10_band_max")
            if band_hi is not None:
                band_right = min(float(band_hi), X_MAX)
                if band_right > left:
                    ax.add_patch(mpatches.Rectangle(
                        (10 ** left, y - BAR_H / 2),
                        10 ** band_right - 10 ** left, BAR_H,
                        facecolor=BAND_STYLE["fill"], edgecolor=BAND_STYLE["edge"],
                        hatch=BAND_STYLE["hatch"], linewidth=0.7,
                        alpha=0.85, zorder=2.5))
            if _extends_right(iv):
                _draw_arrow(ax, 10 ** right, y, "right", style)
            if _extends_left(iv):
                _draw_arrow(ax, 10 ** left, y, "left", style)

        text = _annotation(data_text, method_label(row))
        placed = _place_annotation(text, bars, ANNOT_FS, in_per_dec, arrows=arrows)
        if placed is None:
            continue
        bbox = None
        if placed["inside"]:
            host = next((iv for iv in entry["ivs"]
                         if _span(iv)[0] <= placed["x"] <= _span(iv)[1]), None)
            if host is not None and _bar_style(host, data_text)["hatch"]:
                # A translucent plate keeps dark text legible on the hatched
                # fill.  Deliberately *not* a path-effect halo: path effects
                # would turn these labels into outlines in the PDF and they
                # would stop being searchable.
                bbox = dict(facecolor="white", alpha=0.75, edgecolor="none", pad=1.5)
        ax.text(10 ** placed["x"], y, placed["text"], fontsize=placed["fs"],
                ha=placed["ha"], va="center", color=TEXT_DARK, zorder=4,
                clip_on=True, bbox=bbox)
        report.append(f"{slot['label']:30s} slot={placed['where']:7s} fs={placed['fs']:4.1f}")

    # Legend below the axes.
    handles = [
        mpatches.Patch(facecolor=TONE["include"]["fill"], edgecolor=TONE["include"]["edge"],
                       label="Include"),
        mpatches.Patch(facecolor=TONE["exclude"]["fill"], edgecolor=TONE["exclude"]["edge"],
                       hatch=TONE["exclude"]["hatch"], label="Exclude"),
        Line2D([], [], color="#666666", marker=">", linestyle="none", markersize=8,
               label="Open end"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.048),
              ncol=3, frameon=False, fontsize=LEGEND_FS, handlelength=2.0,
              columnspacing=1.6)

    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=220, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(f"[i] {len(rows)} rows in {len(blocks)} blocks")
    print("[i] annotation placement:")
    for line in report:
        print("      " + line)
    print(f"[i] Wrote {OUT_PNG}")
    print(f"[i] Wrote {OUT_PDF}")


if __name__ == "__main__":
    main()