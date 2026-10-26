"""
plot_constraints_html.py
========================

Produce the *interactive* version of the constraint figure with ``plotly``.
Hovering any bar reveals the citekey, year, title, authors, journal, DOI /
arXiv link, the full (un-abbreviated) ``data`` / ``method`` strings, the raw
constraint expression and the bounding interval in ``log10(m_a / eV)``.

Layout
------
Identical to the PNG / SVG figures, so the three outputs of the pipeline can
be read interchangeably: the same two blocks (include / exclude), the same row
order (``table_index``, i.e. the order of ``ma_constraints.tex``), one row per
paper with one bar per interval, the same colours and the same arrowheads /
uncertainty bands.

Implementation notes
--------------------
* The x axis is **linear in the exponent** ``log10(m_a/eV)`` with tick labels
  rendered as ``10^-26`` … — this keeps bar ``base`` / ``x`` values exact (a
  plotly ``log`` axis would misinterpret a negative ``base``).
* The y axis is **numeric** (one unit per row slot) so that the block bands,
  arrowheads and uncertainty bands can be drawn as data-coordinate shapes
  (category axes cannot host shapes).

Output
------
``docs/constraints/output/ma_constraints.html``
"""

from __future__ import annotations

from pathlib import Path

import plotly.graph_objects as go

import plot_constraints_static as P          # shared layout / styles
from constraint_labels import data_label, method_label, strip_tex

HERE = Path(__file__).resolve().parent
OUT_HTML = HERE / "output" / "ma_constraints.html"

TONES = [
    ("include", "Include"),
    ("exclude", "Exclude"),
]


def _safe(s) -> str:
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def _hover(row: dict, iv: dict) -> str:
    bib = row.get("bib", {}) or {}
    authors = bib.get("authors") or []
    a_short = ", ".join(authors[:4]) + (", et al." if len(authors) > 4 else "")
    url = bib.get("url") or ""
    lo, hi = iv["log10_min"], iv["log10_max"]
    lo_s = "—∞" if P._is_open_lo(lo) else f"10<sup>{lo:.2f}</sup>"
    hi_s = "+∞" if P._is_open_hi(hi) else f"10<sup>{hi:.2f}</sup>"
    band = iv.get("log10_band_max")
    band_s = (f"<br><span style='color:#555'>uncertainty band:</span> "
              f"10<sup>{lo:.2f}</sup>–10<sup>{band:.2f}</sup>" if band is not None else "")
    return (
        f"<b>{_safe(row['citekey'])} ({_safe(bib.get('year', '?'))})</b><br>"
        f"<i>{_safe(strip_tex(bib.get('title', '')))}</i><br>"
        f"{_safe(a_short)}<br>"
        f"<span style='color:#555'>journal:</span> {_safe(bib.get('journal', '?'))}<br>"
        f"<span style='color:#555'>data:</span> {_safe(strip_tex(row.get('data', '')))}<br>"
        f"<span style='color:#555'>method:</span> {_safe(strip_tex(row.get('method', '')))}<br>"
        f"<span style='color:#555'>constraint:</span> {_safe(strip_tex(row.get('constraint_raw', '')))}<br>"
        f"<span style='color:#555'>interval:</span> {lo_s} … {hi_s} eV  ({iv['type']}){band_s}<br>"
        + (f"<a href='{url}' target='_blank'>open paper</a>" if url else "")
    )


def main() -> None:
    data = P._load()
    slots, blocks, total_y = P._build_layout(data)

    fig = go.Figure()

    # Block bands (drawn as a light background per block).
    for block in blocks:
        fig.add_shape(type="rect", xref="x", yref="y",
                      x0=P.X_MIN, x1=P.X_MAX,
                      y0=-(block["top"] - 0.30), y1=-block["bottom"],
                      fillcolor=block["tint"], line=dict(width=0), layer="below")
        fig.add_annotation(x=P.X_MIN + 0.06, y=block_header_y(block, slots),
                           xref="x", yref="y", showarrow=False,
                           xanchor="left", font=dict(size=22, color=block["colour"]),
                           text=f"<b>{_safe(block['header'])}</b>")

    # One trace per tone; the dataset-family shade is set per bar.
    for key, legend_label in TONES:
        xs, bases, ys, hovers, labels, colors = [], [], [], [], [], []
        for slot in slots:
            if slot["kind"] != "row":
                continue
            entry, row = slot["entry"], slot["entry"]["row"]
            y = -(slot["y"] + P.ROW_H / 2)
            data_text = data_label(row)
            for iv in entry["ivs"]:
                if P._tone_of(iv) != key:
                    continue
                left, right = P._span(iv)
                if right <= left:
                    continue
                style = P._bar_style(iv, data_text)
                xs.append(right - left)
                bases.append(left)
                ys.append(y)
                colors.append(style["fill"])
                labels.append(slot["label"])
                hovers.append(_hover(row, iv))
                # Arrowheads at open ends (data-coordinate path shapes).
                if P._extends_right(iv):
                    fig.add_shape(type="path", xref="x", yref="y",
                                  path=(f"M {right - P.ARROW_W} {y - P.BAR_H * 0.30} "
                                        f"L {right - P.ARROW_W} {y + P.BAR_H * 0.30} "
                                        f"L {right} {y} Z"),
                                  fillcolor=style["fill"], line=dict(color=style["edge"], width=1))
                if P._extends_left(iv):
                    fig.add_shape(type="path", xref="x", yref="y",
                                  path=(f"M {left + P.ARROW_W} {y - P.BAR_H * 0.30} "
                                        f"L {left + P.ARROW_W} {y + P.BAR_H * 0.30} "
                                        f"L {left} {y} Z"),
                                  fillcolor=style["fill"], line=dict(color=style["edge"], width=1))
                # Uncertainty band of a range-quoted bound.
                band_hi = iv.get("log10_band_max")
                if band_hi is not None:
                    band_right = min(float(band_hi), P.X_MAX)
                    if band_right > left:
                        fig.add_shape(type="rect", xref="x", yref="y",
                                      x0=left, x1=band_right,
                                      y0=y - P.BAR_H / 2, y1=y + P.BAR_H / 2,
                                      fillcolor="rgba(255,255,255,0.85)",
                                      line=dict(color=P.BAND_STYLE["edge"], width=1,
                                                dash="dash"))
        if not xs:
            continue
        fig.add_trace(go.Bar(
            x=xs, y=ys, base=bases,
            orientation="h", width=P.BAR_H, name=legend_label,
            marker=dict(color=colors, line=dict(color=P.TONE[key]["edge"], width=1)),
            hovertext=hovers, hoverinfo="text", showlegend=True,
            text=labels, textposition="none",
        ))

    # Font sizes follow the PNG/PDF output (×1.6 of the original values), with
    # the row pitch and the margins scaled to match so nothing collides.
    row_px = 42
    fig.update_layout(
        title=dict(
            text=("Observational constraints on the wave-CDM particle mass m<sub>a</sub>"
                  "<br><span style='font-size:19px;color:#5a5a5a'>hover a bar for the "
                  "paper, data, method and constraint</span>"),
            x=0.5, xanchor="center", font=dict(size=29),
        ),
        height=max(1100, row_px * len([s for s in slots if s["kind"] == "row"]) + 320),
        barmode="overlay",
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=290, r=60, t=140, b=110),
        legend=dict(orientation="h", yanchor="bottom", y=-0.07, xanchor="center", x=0.5,
                    font=dict(size=19), title=None),
        xaxis=dict(
            title=dict(text="log<sub>10</sub>(m<sub>a</sub> / eV)", font=dict(size=22)),
            showgrid=True, gridcolor="rgba(0,0,0,0.15)", gridwidth=1,
            zeroline=False, range=[P.X_MIN, P.X_MAX],
            tickvals=[float(p) for p in range(int(P.X_MIN), int(P.X_MAX) + 1)],
            ticktext=[f"10^{p}" for p in range(int(P.X_MIN), int(P.X_MAX) + 1)],
            tickfont=dict(size=21),
        ),
        yaxis=dict(
            showgrid=True, gridcolor="rgba(0,0,0,0.16)", gridwidth=1, griddash="dash",
            zeroline=False, range=[-total_y, 0],
            tickvals=[-(s["y"] + P.ROW_H / 2) for s in slots if s["kind"] == "row"],
            ticktext=[s["label"] for s in slots if s["kind"] == "row"],
            tickfont=dict(size=19),
        ),
    )

    OUT_HTML.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(
        str(OUT_HTML), include_plotlyjs="cdn", full_html=True,
        config={"displaylogo": False,
                "modeBarButtonsToRemove": ["lasso2d", "select2d"],
                "toImageButtonOptions": {"format": "png", "filename": "ma_constraints"}},
    )
    print(f"[i] Wrote {OUT_HTML} ({OUT_HTML.stat().st_size / 1024:.1f} KiB)")


def block_header_y(block: dict, slots: list[dict]) -> float:
    """Numeric y of a block header row (negative, matching the bar rows)."""
    for slot in slots:
        if slot["kind"] == "header" and slot["block"] is block:
            return -(slot["y"] + P.HEADER_H * 0.45)
    return 0.0


if __name__ == "__main__":
    main()