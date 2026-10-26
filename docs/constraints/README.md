# WaveDM.jl — observational constraints on the ultralight dark-matter mass

This directory holds the data pipeline that turns the literature-constraint
table in [`ma_constraints.tex`](./ma_constraints.tex) into a single figure
that can be embedded both in the README (as a PNG) and in the GitHub Pages
documentation site (as an interactive HTML, plus a no-JS SVG fallback).

## Files

| file                                          | role                                                                                 |
| --------------------------------------------- | ------------------------------------------------------------------------------------ |
| [`ma_constraints.tex`](./ma_constraints.tex) | Source of truth: 42 rows × 4 columns (citekey / data / method / constraint).         |
| [`temp/references.bib`](./temp/references.bib) | Bibliography with title, authors, year, DOI / arXiv id for every citekey.         |
| `ma_constraints_data.py`                      | Parses both files into a single structured YAML/JSON.                                |
| `constraint_labels.py`                        | Shared text layer of the three renderers: TeX → Unicode cleaning, `data` / `method` abbreviations, `FirstAuthor et al. Year` labels with `a`/`b` suffixes. |
| `build_constraints_dataset.jl`                | Loads that YAML into a `ParameterSpace.DataSet{ConstraintRow}` (catalog schema).     |
| `plot_constraints_static.py`                  | Matplotlib PNG for README / docs.                                                    |
| `plot_constraints_julia.jl`                   | CairoMakie / `ParameterSpaceVis.plot_constraints_catalog` → PNG/PDF/SVG.               |
| `plot_constraints_html.py`                    | Plotly HTML for the docs site (hover + click → DOI).                                  |
| `plot_constraints_svg.py`                     | Pure-SVG output with `<a href="...">` wrappers (click → DOI, no JS).                 |
| `output/ma_constraints.yaml` / `.json`        | Canonical structured dataset.                                                         |
| `output/ma_constraints.pspace.yaml`           | Same data serialised through `ParameterSpace.save_dataset` (ConstraintRow schema).    |
| `output/ma_constraints_filtered.yaml`         | Filtered subset (2018–2025) demonstrating `ParameterSpace.DataSet` filtering.       |
| `output/ma_constraints.png`                   | Static figure (matplotlib, 220 DPI).                                                  |
| `output/ma_constraints.pdf`                   | Searchable vector PDF (matplotlib, text can be copy-pasted into papers).              |
| `output/ma_constraints_julia.png/pdf/svg`     | Same figure via Julia Makie (`ParameterSpaceVis`).                                    |
| `output/ma_constraints.svg`                   | Self-contained SVG with clickable rectangles (matplotlib + hand-written).             |
| `output/ma_constraints.html`                  | Self-contained interactive HTML (loads plotly.js from CDN).                          |

## Regeneration recipe

```powershell
# 1) Extract structured YAML / JSON from the LaTeX + bib files.
D:\Programs\miniconda3\python.exe docs/constraints/ma_constraints_data.py

# 2) Render the three figures.
D:\Programs\miniconda3\python.exe docs/constraints/plot_constraints_static.py
D:\Programs\miniconda3\python.exe docs/constraints/plot_constraints_html.py
D:\Programs\miniconda3\python.exe docs/constraints/plot_constraints_svg.py

# 3) (Optional) Julia pipeline via ParameterSpace.jl + ParameterSpaceVis.jl.
#    Requires in your default v1.13 env:
#        Pkg.develop(path="E:/JuliaAstroSim/ParameterSpaceVis.jl")
#        Pkg.add("YAML")
#    Then either invoke directly:
$env:JULIA_DEPOT_PATH = "E:\.julia"
julia --project=@. docs/constraints/build_constraints_dataset.jl
julia --project=@. docs/constraints/plot_constraints_julia.jl
#    Or use the bundled helper:
powershell -ExecutionPolicy Bypass -File docs/constraints/run_julia_pipeline.ps1

# 4) Copy artefacts into the Documenter assets folder so `julia docs/make.jl`
#    pulls them into the docs site build.
Copy-Item docs/constraints/output/ma_constraints.png  docs/src/assets/ma_constraints.png  -Force
Copy-Item docs/constraints/output/ma_constraints.pdf  docs/src/assets/ma_constraints.pdf  -Force
Copy-Item docs/constraints/output/ma_constraints.svg  docs/src/assets/ma_constraints.svg  -Force
Copy-Item docs/constraints/output/ma_constraints.html docs/src/assets/ma_constraints.html -Force

# 5) Build docs (only needed if you change *.md).
julia --project=docs docs/make.jl
```

The assets list in `docs/make.jl` already lists
`assets/ma_constraints.png`, `…svg`, and `…html`, so Documenter copies them
into `docs/build/assets/` automatically once they exist.

## Figure layout (shared by PNG / PDF / SVG / HTML)

* **Two blocks, two tones.** Rows are split into *Include bounds — allowed /
  favoured regions* (green) and *Exclude bounds — ruled-out regions* (orange,
  hatched).  Within a tone every *dataset family* (rotation curves, dwarf
  galaxies, Lyα forest, 21 cm, CMB, lensing, black holes, pulsar timing, …)
  gets its own shade, so rows sharing a dataset read as a group; the legend
  explains the two tones, not the shades.  A dashed horizontal guide runs
  through every row so a label can be traced across the figure.
* **Row order = table order.** ``table_index`` (written by
  `ma_constraints_data.py`) keeps the figure rows in the same sequence as
  `ma_constraints.tex`, so the figure can be read side by side with the paper
  table.
* **One row per paper.** A work quoting several intervals draws several bars on
  its single row (Bar et al. 2019, Marsh et al. 2019, Wardana et al. 2026, …)
  instead of repeating its label.  Two *different* papers that share an
  `Author et al. Year` label get `a` / `b` suffixes.
* **Annotations.** `data — method` is written next to each bar in the
  abbreviated form defined in `constraint_labels.py` (`rotation curve` → `RC`,
  `soliton–halo relation` → `SH relation`, `subhalo mass function` → `SHMF`,
  `Lyα flux power spectrum` → `Lyα P(k)`, …).  The width of every annotation is
  measured before drawing, and the text goes into the first slot of its row
  that is wide enough: **inside a band that has room for it** (widest bar
  first — the `m > x` / `m ≤ x` bands on either edge of the figure and wide
  interior ranges like Bar et al. 2022 alike) → **hugging the right end** of
  the bars → the free **gap between two bars** → **hugging the left end** of
  the bars (right-aligned, so the label stays glued to its bar instead of
  floating at the axis edge).  Nothing can spill past the axes, and text over
  a hatched fill gets a translucent plate.  The full, un-abbreviated strings
  stay in the YAML and in the HTML/SVG hover tooltips.
* **Open ends.** A bound that continues beyond the plotted range
  (1e-26 … 1e-18 eV) carries an arrowhead at that end; a bound whose *quoted
  value* is itself a range (Marsh & Niemeyer 2019: `m ≳ 0.6–1 ×10⁻¹⁹ eV`)
  carries a hatched uncertainty band, stored as `log10_band_max` on the
  interval.

## Constraint expression syntax

`ma_constraints_data.py` parses each cell of the ``limit`` column into one or
more intervals on $\log_{10}(m_a/\mathrm{eV})$. The recognised patterns are:

* `$a \sim b \times 10^{-E}$`  → `[log10(a) − E, log10(b) − E]` (range)
* `$10^{-X} \sim 10^{-Y}$`     → `[−X, −Y]`
* `$m = a^{+u}_{-l} \times 10^{-E}$`  → point-estimate range with 1σ errors,
  recorded as an ordinary (favoured) include range
* `$a < m < b$` / `$a \le m \le b$` (with `\le`, `\leq`, `\lesssim`, `~`, etc.)
* `$m \gtrsim a \to b \times 10^{-E}$` → lower bound anchored at `a` plus an
  uncertainty band up to `b` (the quoted limit is itself a range, e.g. Marsh &
  Niemeyer 2019)
* `$m > b$` / `$m \ge b$` / `$m \gtrsim b$` → lower-bound
* `$m < a$` / `$m \le a$` / `$m \lesssim a$` → upper-bound
* `$a - b$` style ranges for BHSR (e.g. `7×10⁻¹⁴ – 2×10⁻¹¹`)
* `constrain m ≤ b: ρ < 6 GeV/cm³` (PTA-style with extra annotation)
* multi-interval cells separated by `\\&`, `&`, `;`, or `,`

Note that `&` inside a *data* / *method* cell must be escaped as `\&` in the
LaTeX table (`Segue I \& II UFD`); the row splitter only cuts on unescaped
`&`, and `\&` stays a separator when it appears in the *limit* cell
(`7×10⁻²⁰–1×10⁻¹⁶ eV \& …`).

If parsing fails for a row, the script logs a warning and emits an empty
`constraints: []` entry — the offending row is still listed in the YAML so
the table layout in the README figure does not change.

## Why ParameterSpace.jl?

The pipeline uses `ParameterSpace.ConstraintRow` / `DataSet` from the
[JuliaAstroSim/ParameterSpace.jl](https://github.com/JuliaAstroSim/ParameterSpace.jl)
package as a *data-modelling* layer — one citekey = one `ConstraintRow`
carrying 0..N `ConstraintInterval`s plus the free-text `data` / `method`
strings.  The output `ma_constraints.pspace.yaml` is round-trippable through
`load_dataset` / `save_dataset`, and the same `ConstraintRow` objects can be
fed directly into `ParameterSpaceVis.plot_constraints_catalog` to render the
catalog figure (see `plot_constraints_julia.jl`).  Both extensions (the
extended `DataSource` with `authors`/`doi`/`eprint`/`bibcode`/`url` fields, and
the new `ConstraintRow` / `ConstraintInterval` types) are **backward-compatible
additions** — older `ConstraintData`-style YAMLs still load and round-trip.