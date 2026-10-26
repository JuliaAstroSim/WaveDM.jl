"""
plot_constraints_svg.py
=======================

Render the constraint figure as a **self-contained SVG** in which every
constraint rectangle is wrapped in an ``<a href="https://doi.org/…">``
hyperlink: clicking a bar opens the corresponding paper.  This is the
no-JavaScript fallback for platforms that strip ``<iframe>`` (GitHub README,
static documentation), so it carries the same layout as the PNG/PDF figure:

* the same two blocks (include / exclude), each ordered like
  ``ma_constraints.tex``;
* one row per paper, several bars when a work quotes several intervals;
* the same abbreviated ``data — method`` annotations, placed by the *same*
  slot-fitting code as the matplotlib output (``plot_constraints_static``),
  with a text-width estimate so a label can never run off the plot;
* arrowheads on open bounds, uncertainty bands on range-quoted bounds.

Output
------
``docs/constraints/output/ma_constraints.svg``
"""

from __future__ import annotations

import math
from pathlib import Path

import plot_constraints_static as P          # layout + placement, shared with the PNG
from constraint_labels import data_label, method_label, strip_tex

HERE = Path(__file__).resolve().parent
OUT_SVG = HERE / "output" / "ma_constraints.svg"

PX_PER_DECADE = 152.0
ROW_PX = 34.0
# Font sizes, scaled ×1.6 to match the PNG/PDF output (the SVG geometry above
# is in pixels, so a larger font is exactly "bigger relative to the plot").
FS_TITLE = 34
FS_TICK = 24
FS_TICK_SUB = 18
FS_AXLABEL = 27
FS_AXLABEL_SUB = 19
FS_BLOCK = 24
FS_ROW_LABEL = 21.5
FS_ANNOT = 21.5
FS_LEGEND = 20
MARGIN_LEFT = 250.0
MARGIN_RIGHT = 70.0
MARGIN_TOP = 72.0
MARGIN_BOTTOM = 150.0
BAR_PX = 0.80 * ROW_PX
ARROW_PX = 0.085 * PX_PER_DECADE

# Rough per-character advance widths (fraction of the font size) for a
# sans-serif face.  Deliberately on the generous side: the SVG is rendered
# with whatever font the browser has, and an over-estimate only costs a little
# white space while an under-estimate would let text escape the plot.
_CHAR_W = {
    "iljtf.,:;'|!()[]": 0.32,
    "mwMW": 0.92,
    "ABCDEFGHKNOPQRSTUVXYZ": 0.70,
    "0123456789": 0.58,
    " ": 0.34,
}


def _char_width(ch: str) -> float:
    for chars, w in _CHAR_W.items():
        if ch in chars:
            return w
    return 0.56


def _svg_text_width_px(text: str, fontsize: float) -> float:
    """Estimated rendered width of ``text`` in SVG pixels."""
    plain = text.replace("$", "").replace("\\mathit{", "").replace("}", "")
    return sum(_char_width(c) for c in plain) * fontsize * 1.04


def _escap(s: str) -> str:
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def main() -> None:
    data = P._load()
    slots, blocks, total_y = P._build_layout(data)

    width = MARGIN_LEFT + P.NDEC * PX_PER_DECADE + MARGIN_RIGHT
    height = MARGIN_TOP + total_y * ROW_PX + MARGIN_BOTTOM
    axes_right = MARGIN_LEFT + P.NDEC * PX_PER_DECADE

    def x_px(log10: float) -> float:
        return MARGIN_LEFT + (log10 - P.X_MIN) * PX_PER_DECADE

    def y_px(slot_y: float) -> float:
        return MARGIN_TOP + slot_y * ROW_PX

    out: list[str] = []
    add = out.append
    add(f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'xmlns:xlink="http://www.w3.org/1999/xlink" width="{width:.0f}" '
        f'height="{height:.0f}" viewBox="0 0 {width:.0f} {height:.0f}" '
        f'font-family="Helvetica, Arial, DejaVu Sans, sans-serif">')
    add(f'<rect x="0" y="0" width="{width:.0f}" height="{height:.0f}" fill="white"/>')

    # Hatch patterns so the SVG uses the same texture language as the PNG
    # (diagonal hatch = ruled out, reverse hatch = uncertainty band).
    add('<defs>'
        '<pattern id="hatch-exclude" width="9" height="9" patternUnits="userSpaceOnUse" '
        'patternTransform="rotate(45)"><line x1="0" y1="0" x2="0" y2="9" '
        'stroke="#8a4b12" stroke-width="1.3" stroke-opacity="0.45"/></pattern>'
        '<pattern id="hatch-band" width="7" height="7" patternUnits="userSpaceOnUse" '
        'patternTransform="rotate(-45)"><line x1="0" y1="0" x2="0" y2="7" '
        'stroke="#3d6b1f" stroke-width="1.2" stroke-opacity="0.55"/></pattern>'
        '</defs>')

    # Title.
    add(f'<text x="{width / 2:.0f}" y="{FS_TITLE + 14}" text-anchor="middle" '
        f'font-size="{FS_TITLE}" '
        f'font-weight="bold">Observational constraints on the wave-CDM particle mass '
        f'<tspan font-style="italic">m</tspan><tspan baseline-shift="sub" '
        f'font-size="{FS_TITLE * 0.67:.0f}">a</tspan></text>')

    # Block bands, ribbons and headers.
    for block in blocks:
        top = y_px(block["top"]) - 0.30 * ROW_PX
        bottom = y_px(block["bottom"])
        add(f'<rect x="{MARGIN_LEFT:.0f}" y="{top:.1f}" width="{axes_right - MARGIN_LEFT:.0f}" '
            f'height="{bottom - top:.1f}" fill="{block["tint"]}"/>')
        add(f'<rect x="{MARGIN_LEFT - 22:.1f}" y="{y_px(block["top"]) - 0.22 * ROW_PX:.1f}" '
            f'width="12" height="{bottom - top:.1f}" fill="{block["colour"]}" opacity="0.85"/>')

    # Decade gridlines + ticks.
    for p in range(int(P.X_MIN), int(P.X_MAX) + 1):
        x = x_px(float(p))
        add(f'<line x1="{x:.1f}" y1="{MARGIN_TOP:.0f}" x2="{x:.1f}" '
            f'y2="{MARGIN_TOP + total_y * ROW_PX:.1f}" stroke="rgba(0,0,0,0.13)" '
            f'stroke-width="1" stroke-dasharray="4 4"/>')
        add(f'<text x="{x:.1f}" y="{MARGIN_TOP + total_y * ROW_PX + 34:.1f}" '
            f'text-anchor="middle" font-size="{FS_TICK}">10<tspan baseline-shift="super" '
            f'font-size="{FS_TICK_SUB}">{p}</tspan></text>')
    add(f'<text x="{(MARGIN_LEFT + axes_right) / 2:.0f}" '
        f'y="{MARGIN_TOP + total_y * ROW_PX + 76:.1f}" text-anchor="middle" '
        f'font-size="{FS_AXLABEL}">'
        f'log<tspan baseline-shift="sub" font-size="{FS_AXLABEL_SUB}">10</tspan>(m'
        f'<tspan baseline-shift="sub" font-size="{FS_AXLABEL_SUB}">a</tspan> / eV)</text>')

    # Rows.
    report: list[str] = []
    overhang = -1e9
    for slot in slots:
        if slot["kind"] == "header":
            block = slot["block"]
            add(f'<text x="{MARGIN_LEFT + 10:.0f}" '
                f'y="{y_px(slot["y"]) + 0.45 * 0.85 * ROW_PX:.1f}" font-size="{FS_BLOCK}" '
                f'font-weight="bold" fill="{block["colour"]}">'
                f'{_escap(block["header"])}</text>')
            continue

        entry, row = slot["entry"], slot["entry"]["row"]
        y = y_px(slot["y"]) + ROW_PX / 2
        # Horizontal dashed guide through the row (same as the PNG/PDF).
        add(f'<line x1="{MARGIN_LEFT:.0f}" y1="{y:.1f}" x2="{axes_right:.0f}" y2="{y:.1f}" '
            f'stroke="rgba(0,0,0,0.16)" stroke-width="1" stroke-dasharray="5 5"/>')
        add(f'<text x="{MARGIN_LEFT - 30:.0f}" y="{y + 6.5:.1f}" text-anchor="end" '
            f'font-size="{FS_ROW_LABEL}">{_escap(slot["label"])}</text>')

        data_text = data_label(row)
        bars: list[tuple[float, float]] = []
        arrows: list[tuple[bool, bool]] = []
        for iv in entry["ivs"]:
            left, right = P._span(iv)
            if right <= left:
                continue
            bars.append((left, right))
            arrows.append((P._extends_left(iv), P._extends_right(iv)))
            style = P._bar_style(iv, data_text)
            x0, x1 = x_px(left), x_px(right)
            bib = row.get("bib", {}) or {}
            url = bib.get("url") or ""
            tip = (f'{row["citekey"]} ({bib.get("year", "?")}) | '
                   f'{strip_tex(bib.get("title", ""))[:90]} | '
                   f'{data_label(row)} — {method_label(row)} | '
                   f'{strip_tex(row.get("constraint_raw", ""))}')
            rect = (f'<rect x="{x0:.1f}" y="{y - BAR_PX / 2:.1f}" width="{x1 - x0:.1f}" '
                    f'height="{BAR_PX:.1f}" fill="{style["fill"]}" stroke="{style["edge"]}" '
                    f'stroke-width="0.8"/>')
            if style["hatch"]:
                pattern = "hatch-exclude" if P._tone_of(iv) == "exclude" else "hatch-mass"
                rect += (f'<rect x="{x0:.1f}" y="{y - BAR_PX / 2:.1f}" width="{x1 - x0:.1f}" '
                         f'height="{BAR_PX:.1f}" fill="url(#{pattern})" stroke="none"/>')
            if url:
                add(f'<a xlink:href="{_escap(url)}" target="_blank"><title>{_escap(tip)}</title>'
                    f'{rect}</a>')
            else:
                add(f'{rect}<title>{_escap(tip)}</title>')

            # Uncertainty band (a bound whose quoted value is a range).
            band_hi = iv.get("log10_band_max")
            if band_hi is not None:
                band_right = min(float(band_hi), P.X_MAX)
                if band_right > left:
                    bx0, bx1 = x_px(left), x_px(band_right)
                    add(f'<rect x="{bx0:.1f}" y="{y - BAR_PX / 2:.1f}" width="{bx1 - bx0:.1f}" '
                        f'height="{BAR_PX:.1f}" fill="white" fill-opacity="0.85" '
                        f'stroke="{P.BAND_STYLE["edge"]}" stroke-width="0.9" '
                        f'stroke-dasharray="3 2"/>'
                        f'<rect x="{bx0:.1f}" y="{y - BAR_PX / 2:.1f}" width="{bx1 - bx0:.1f}" '
                        f'height="{BAR_PX:.1f}" fill="url(#hatch-band)" stroke="none"/>')
            # Open ends.
            if P._extends_right(iv):
                add(f'<polygon points="{x1 - ARROW_PX:.1f},{y - 0.30 * BAR_PX:.1f} '
                    f'{x1 - ARROW_PX:.1f},{y + 0.30 * BAR_PX:.1f} {x1:.1f},{y:.1f}" '
                    f'fill="{style["fill"]}" stroke="{style["edge"]}" stroke-width="0.8"/>')
            if P._extends_left(iv):
                add(f'<polygon points="{x0 + ARROW_PX:.1f},{y - 0.30 * BAR_PX:.1f} '
                    f'{x0 + ARROW_PX:.1f},{y + 0.30 * BAR_PX:.1f} {x0:.1f},{y:.1f}" '
                    f'fill="{style["fill"]}" stroke="{style["edge"]}" stroke-width="0.8"/>')

        # Annotation: identical slot-fitting logic, SVG width estimate.
        data_txt, method_txt = data_text, method_label(row)
        ann = f"{data_txt}   {method_txt}"
        placed = P._place_annotation(ann, bars, FS_ANNOT, PX_PER_DECADE,
                                     width_fn=_svg_text_width_px, arrows=arrows)
        if placed is None:
            continue
        anchor = x_px(placed["x"])
        anchor_attr = ' text-anchor="end"' if placed["ha"] == "right" else ""
        text_w = _svg_text_width_px(ann, placed["fs"])
        # Left edge of the rendered text depends on the alignment.
        if placed["ha"] == "right":
            left_edge = anchor - text_w
            right_edge = anchor
        else:
            left_edge = anchor
            right_edge = anchor + text_w
        overhang = max(overhang, MARGIN_LEFT - left_edge, right_edge - axes_right)
        halo = ""
        if placed["inside"]:
            host = next((iv for iv in entry["ivs"]
                         if P._span(iv)[0] <= placed["x"] <= P._span(iv)[1]), None)
            if host is not None and P._bar_style(host, data_text)["hatch"]:
                # Translucent plate behind the label (same treatment as the
                # PNG/PDF) so the text stays crisp over the hatched fill.
                plate_x = anchor - text_w if placed["ha"] == "right" else anchor
                add(f'<rect x="{plate_x - 3:.1f}" y="{y - 0.44 * ROW_PX:.1f}" '
                    f'width="{text_w + 6:.1f}" height="{0.88 * ROW_PX:.1f}" '
                    f'fill="white" fill-opacity="0.75"/>')
        add(f'<text x="{anchor:.1f}" y="{y + 4.8:.1f}" font-size="{placed["fs"]:.1f}" '
            f'fill="{P.TEXT_DARK}"{anchor_attr}>{_escap(data_txt)}'
            f'<tspan font-style="italic">   {_escap(method_txt)}</tspan></text>')
        report.append(f'{slot["label"]:30s} slot={placed["where"]:11s} ha={placed["ha"]:5s} '
                      f'fs={placed["fs"]:4.1f}')

    # Legend.
    legend = [
        (P.TONE["include"]["fill"], P.TONE["include"]["edge"], None, None, "Include"),
        (P.TONE["exclude"]["fill"], P.TONE["exclude"]["edge"], "hatch-exclude", None,
         "Exclude"),
    ]
    ly = MARGIN_TOP + total_y * ROW_PX + 108
    lx = MARGIN_LEFT
    for face, edge, pattern, dash, label in legend:
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        add(f'<rect x="{lx:.0f}" y="{ly:.0f}" width="34" height="17" fill="{face}" '
            f'stroke="{edge}" stroke-width="0.8"{dash_attr}/>')
        if pattern:
            add(f'<rect x="{lx:.0f}" y="{ly:.0f}" width="34" height="17" '
                f'fill="url(#{pattern})" stroke="none"/>')
        add(f'<text x="{lx + 42:.0f}" y="{ly + 15:.0f}" font-size="{FS_LEGEND}">'
            f'{_escap(label)}</text>')
        lx += 60 + 0.58 * FS_LEGEND * len(label)
    add(f'<polygon points="{lx + 12:.0f},{ly:.0f} {lx + 12:.0f},{ly + 17:.0f} '
        f'{lx + 30:.0f},{ly + 8.5:.0f}" fill="#666666"/>')
    add(f'<text x="{lx + 38:.0f}" y="{ly + 15:.0f}" font-size="{FS_LEGEND}">Open end</text>')

    add("</svg>")
    OUT_SVG.parent.mkdir(parents=True, exist_ok=True)
    OUT_SVG.write_text("\n".join(out), encoding="utf-8")
    print(f"[i] Wrote {OUT_SVG} ({OUT_SVG.stat().st_size / 1024:.0f} KiB)")
    print(f"[i] max annotation overhang past the axes: {overhang:.1f} px "
          f"({'OK' if overhang <= 0 else 'OVERFLOW'})")
    for line in report:
        print("      " + line)


if __name__ == "__main__":
    main()