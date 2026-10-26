"""
ma_constraints_data.py
======================

Parse the LaTeX table ``ma_constraints.tex`` together with the bibliography
``references.bib`` and emit ``ma_constraints.yaml`` / ``ma_constraints.json``
that describe every observation-derived constraint on the ultralight (wave /
fuzzy) dark-matter particle mass.

Usage
-----
    D:\\Programs\\miniconda3\\python.exe docs/constraints/ma_constraints_data.py

YAML schema::

    citekey:
      citekey: str
      data:    str
      method:  str
      constraint_raw: str
      category: str
      table_index: int          # 0-based row position in ma_constraints.tex
      constraints: list of intervals on log10(m_a / eV)
        - type: 'include_range' | 'exclude_range' | 'include_above' | 'include_below' | 'exclude_above' | 'exclude_below'
          log10_min: float
          log10_max: float
          log10_band_max: float   # optional: stronger end of a *quoted range*
                                  # that anchors a lower bound (e.g. "m ≳ 0.6–1 ×10⁻¹⁹ eV"); the
                                  # figures draw it as an uncertainty band on top of the bound.
          label: str
      bib:
        title:   str
        authors: list[str]
        year:    int
        journal: str
        doi:     str | ''
        eprint:  str | ''
        url:     str

The constraint parser is hand-written and runs against the *raw* LaTeX snippet
so it preserves the ``\\times`` / ``\\sim`` / ``\\le`` / ``\\ge`` / ``10^{...}``
syntax.  It is robust to:

* ``$a \\sim b \\times 10^{-E}$``  (range with leading coefficient)
* ``$10^{-X} \\sim 10^{-Y}$``     (range with pure power-of-ten)
* ``$m \\sim 2.5^{+3.6}_{-2.0} \\times 10^{-E}$``  (point estimate with errors)
* ``exclude $a < m < b$``         (excluding a finite range)
* ``exclude $a \\le m \\le b$``
* ``$m > b$`` / ``$m \\ge b$`` / ``$m \\gtrsim b$``
* ``$m < a$`` / ``$m \\le a$`` / ``$m \\lesssim a$``
* ``$10^{-X} \\le m \\le 10^{-Y}$``
* ``$a-b$`` style BH superradiance ranges (e.g. ``7\\times10^{-14}-2\\times10^{-11}``)
* multi-interval cells separated by ``\\&`` or ``;``

All logs are base-10 with respect to eV, i.e. the column ``log10_min`` is
``\\log_{10}(m_a / eV)`` for the lower edge of the interval.
"""

from __future__ import annotations

import json
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import yaml

# --------------------------------------------------------------------------- #
# Paths                                                                       #
# --------------------------------------------------------------------------- #

HERE = Path(__file__).resolve().parent
TEX_FILE = HERE / "ma_constraints.tex"
BIB_FILE = HERE / "temp" / "references.bib"
OUT_DIR = HERE / "output"
OUT_DIR.mkdir(parents=True, exist_ok=True)
YAML_OUT = OUT_DIR / "ma_constraints.yaml"
JSON_OUT = OUT_DIR / "ma_constraints.json"

# --------------------------------------------------------------------------- #
# Bibliography parsing                                                        #
# --------------------------------------------------------------------------- #


@dataclass
class BibEntry:
    citekey: str
    title: str = ""
    authors: list[str] = field(default_factory=list)
    year: int = 0
    journal: str = ""
    doi: str = ""
    eprint: str = ""

    @property
    def url(self) -> str:
        if self.doi:
            return f"https://doi.org/{self.doi}"
        if self.eprint:
            return f"https://arxiv.org/abs/{self.eprint}"
        return ""

    @property
    def short_authors(self) -> str:
        if not self.authors:
            return "Unknown"
        if len(self.authors) <= 3:
            return ", ".join(self.authors)
        return f"{self.authors[0]} et al."


def _strip_braces(s: str) -> str:
    return s.replace("{", "").replace("}", "").replace("\\textasciitilde", "~")


def _parse_authors(author_field: str) -> list[str]:
    author_field = _strip_braces(author_field).strip()
    parts = re.split(r"\s+and\s+", author_field, flags=re.IGNORECASE)
    out: list[str] = []
    for part in parts:
        part = part.strip()
        if not part:
            continue
        out.append(re.sub(r"\s+", " ", part))
    return out


def parse_bib(path: Path) -> dict[str, BibEntry]:
    text = path.read_text(encoding="utf-8")
    entries: dict[str, BibEntry] = {}
    pattern = re.compile(r"@\w+\{([^,]+),(.*?)\n\}\s*\n", re.DOTALL)
    # Match a bib field with either ``field = {value},`` or ``field = value,`` (no braces).
    field_pattern = re.compile(
        r"\s*(\w+)\s*=\s*(?:\{(.*?)\}|(\d+))\s*,?\s*\n",
        re.DOTALL,
    )
    for m in pattern.finditer(text):
        key, body = m.group(1), m.group(2)
        entry = BibEntry(citekey=key.strip())
        for line_match in field_pattern.finditer(body):
            field_name = line_match.group(1).lower()
            raw = line_match.group(2) if line_match.group(2) is not None else line_match.group(3)
            value = _strip_braces(raw).strip() if line_match.group(2) is not None else raw.strip()
            if field_name == "title":
                entry.title = value
            elif field_name == "author":
                entry.authors = _parse_authors(value)
            elif field_name == "year":
                try:
                    entry.year = int(value)
                except ValueError:
                    m2 = re.search(r"(\d{4})", value)
                    entry.year = int(m2.group(1)) if m2 else 0
            elif field_name in ("journal", "shorttitle"):
                entry.journal = value
            elif field_name == "doi":
                entry.doi = value.replace("\\", "")
            elif field_name == "eprint":
                entry.eprint = value
        entries[entry.citekey] = entry
    return entries


# --------------------------------------------------------------------------- #
# Numeric helpers                                                             #
# --------------------------------------------------------------------------- #


def _log10(x: float) -> float:
    return math.log10(x)


# Regex pieces
_COEF = r"(?P<coef>[0-9]+(?:\.[0-9]+)?)"  # like "2" or "2.5" (no leading dot)
_COEF_PLAIN = r"(?:[0-9]+(?:\.[0-9]+)?)"
_POWER = r"10\^\{-(?P<exp>[0-9]+(?:\.[0-9]+)?)\}"  # like "10^{-21}" or "10^{-21.5}"
# Same shape as _POWER, but no named group — safe to nest in larger regex.
_POWER_PLAIN = r"10\^\{-(?:[0-9]+(?:\.[0-9]+)?)\}"
# A "value expression" can be either a pure power-of-ten, or coef x 10^{-E}.
# Use NON-capturing wrappers so we can re-use the same pieces multiple times in one
# regex without triggering "redefinition of group name" errors.
_VAL_COEF_PART = rf"(?:{_COEF_PLAIN}\s*(?:\\times|×|x)\s*)?"
_VAL_FULL_NAMED = rf"{_VAL_COEF_PART}{_POWER}"  # for top-level use only
_VAL_FULL = rf"{_VAL_COEF_PART}{_POWER_PLAIN}"   # safe to embed in larger regex
_VAL = rf"(?:{_VAL_FULL}|{_POWER_PLAIN})"


def _coef_to_log10(coef: str | None, exp: str | None) -> float | None:
    """Convert (coefficient, optional exponent text) to log10."""
    if exp is None:
        # Pure number or pure 10^{-X}.
        if coef is None:
            return None
        try:
            return _log10(float(coef))
        except ValueError:
            return None
    # ``10^{-X}`` form: log10 = -X; with leading coef: log10(coef) - X.
    if coef is None:
        return -float(exp)
    return _log10(float(coef)) - float(exp)


def _value_match_to_log10(value_text: str) -> float | None:
    """Best-effort: turn a single value-text like '2.5 x 10^{-21}' into log10."""
    s = value_text.replace(" ", "").replace("\\,", "")
    m = re.match(rf"{_COEF}\s*(?:\\times|×|x)\s*{_POWER}\s*$", s)
    if m:
        return _coef_to_log10(m.group("coef"), m.group("exp"))
    m = re.match(rf"{_POWER}\s*$", s)
    if m:
        return -float(m.group("exp"))
    m = re.match(rf"{_COEF}\s*$", s)
    if m:
        try:
            return _log10(float(m.group("coef")))
        except ValueError:
            return None
    return None


def _extract_first_value(text: str) -> tuple[float | None, float | None]:
    """
    Return (lo, hi) of the first numeric expression in ``text`` after LaTeX
    cleanup.  Both values are log10 of m_a in eV.

    Recognises::

        10^{-X}
        a \\times 10^{-X}
        a x 10^{-X}
        a (just a plain coefficient — gives log10(a))

    Returns (None, None) on failure.
    """
    m = re.search(rf"{_COEF}\s*(?:\\times|×|x)\s*{_POWER}", text)
    if m:
        return _coef_to_log10(m.group("coef"), m.group("exp")), None
    m = re.search(_POWER, text)
    if m:
        return _coef_to_log10(None, m.group("exp")), None
    m = re.search(_COEF, text)
    if m:
        try:
            return _log10(float(m.group("coef"))), None
        except ValueError:
            return None
    return None, None


def _safe_log10(s: str) -> float | None:
    s = s.strip()
    # Match 10^{-X} first (preferred); fall back to plain float.
    m = re.match(_POWER + r"\s*$", s.replace(" ", ""))
    if m:
        return -float(m.group("exp"))
    m = re.match(rf"{_COEF}\s*(?:\\times|×|x)\s*{_POWER}\s*$", s.replace(" ", ""))
    if m:
        return _coef_to_log10(m.group("coef"), m.group("exp"))
    m = re.match(rf"({_COEF})\s*$", s.replace(" ", ""))
    if m:
        try:
            return _log10(float(m.group(1)))
        except ValueError:
            return None
    return None


# --------------------------------------------------------------------------- #
# Constraint expression parsing                                               #
# --------------------------------------------------------------------------- #


def _split_segments(raw: str) -> list[str]:
    """Split the constraint cell on ``\\&``, ``;``, ``,`` and the word 'and'."""
    text = re.sub(r"\s*\\&\s*", ";", raw)
    text = re.sub(r"\s*&\s*", ";", text)
    text = re.sub(r"\s+and\s+", ";", text, flags=re.IGNORECASE)
    text = re.sub(r"\s*,\s*", ";", text)
    return [p.strip() for p in re.split(r";", text) if p.strip()]


def _normalise_segment(s: str) -> str:
    """Normalise whitespace and unify common TeX shorthands for matching."""
    s = re.sub(r"\s+", " ", s)
    s = s.replace("\\,", "").replace("\\!", "").replace("\u2009", "").replace("\u202f", "")
    # Strip surrounding $ ... $ delimiters from inline math segments.
    s = s.replace("$", "")
    # Strip \mathrm{eV} / \mathrm{cm^{-2}} / \mathrm{...} noise from units.
    s = re.sub(r"\\mathrm\{[^}]*\}", "", s)
    # Collapse common LaTeX relation tokens so our regex character classes match.
    # Order matters — replace the longer variants first.
    s = s.replace("\\gtrsim", ">=")
    s = s.replace("\\lesssim", "<=")
    s = s.replace("\\geq", ">=")
    s = s.replace("\\leq", "<=")
    s = s.replace("\\ge", ">=")
    s = s.replace("\\le", "<=")
    # Em-dash separator ``--`` is used as a range marker for BHSR constraints.
    s = s.replace("--", "-")
    # ``m_a`` -> ``m``.
    s = re.sub(r"\bm\s*_\s*[a-zA-Z]\b", "m", s)
    return s.strip()


def parse_constraint(raw: str) -> list[dict]:
    """Convert the TeX snippet from the 'limit' column to log10(m_a/eV) intervals."""
    segments = _split_segments(raw)
    is_exclude = bool(re.search(r"\bexclude\b|\bexcludes\b", raw, re.IGNORECASE))
    kind_pref = "exclude" if is_exclude else "include"

    out: list[dict] = []
    for seg in segments:
        norm = _normalise_segment(seg)

        # Pattern A: $m = a^{+u}_{-l} x 10^{-E}$  (point estimate with asymmetric errors)
        # A *measured* mass with a 1σ confidence interval.  It is recorded as
        # an ordinary (favoured) range interval so it is grouped and drawn with
        # the include bounds — the 1σ character lives in the label.
        m = re.search(
            rf"m\s*=\s*{_COEF}_\{{-(?P<err_lo>[0-9]+(?:\.[0-9]+)?)\}}\^\{{\+?(?P<err_hi>[0-9]+(?:\.[0-9]+)?)\}}\s*(?:\\times|×|x)\s*{_POWER}",
            norm,
        )
        if m:
            cen = float(m.group("coef")) * (10 ** -float(m.group("exp")))
            lo = cen - 10 ** -float(m.group("exp")) * float(m.group("err_lo"))
            hi = cen + 10 ** -float(m.group("exp")) * float(m.group("err_hi"))
            cen_log = _log10(cen)
            out.append(
                dict(
                    type="include_range",
                    log10_min=_log10(max(lo, 1e-30)),
                    log10_max=_log10(max(hi, 1e-30)),
                    label=f"m ≈ 10^{cen_log:.2f} eV (1σ: {m.group('coef')} × 10^-{m.group('exp')})",
                )
            )
            continue

        # Pattern B: exclude/require a < m < b  OR a ≤ m ≤ b  (or ≲ / ≳)
        # Two side expressions joined via 'm'.
        m = re.search(
            rf"({_VAL}|10\^\{{-[0-9]+\}})\s*(?:<=|<|≲|≥|>|≳|~)\s*m\s*(?:<=|<|≲|≥|>|≳|~)\s*({_VAL}|10\^\{{-[0-9]+\}})",
            norm,
        )
        if m:
            a = _value_match_to_log10(m.group(1).strip())
            b = _value_match_to_log10(m.group(2).strip())
            if a is not None and b is not None:
                lo, hi = sorted([a, b])
                # Decide direction: "exclude a < m < b" or "include a < m < b"
                if is_exclude:
                    t = f"exclude_range"
                else:
                    t = f"include_range"
                out.append(
                    dict(
                        type=t,
                        log10_min=lo,
                        log10_max=hi,
                        label=f"10^{lo:.2f} ≲ m/eV ≲ 10^{hi:.2f}",
                    )
                )
                continue

        # Pattern E (reordered): 'a - b' style ranges (e.g. $7\times10^{-14}-2\times10^{-11} \mathrm{eV}$)
        m = re.search(
            rf"({_VAL_FULL})\s*-\s*({_VAL_FULL})",
            norm,
        )
        if m:
            a = _value_match_to_log10(m.group(1).strip())
            b = _value_match_to_log10(m.group(2).strip())
            if a is not None and b is not None:
                lo, hi = sorted([a, b])
                t = "exclude_range" if is_exclude else "include_range"
                out.append(
                    dict(
                        type=t,
                        log10_min=lo,
                        log10_max=hi,
                        label=f"10^{lo:.2f} ≲ m/eV ≲ 10^{hi:.2f}",
                    )
                )
                continue

        # Pattern F0: ``m \gtrsim a \to b \times 10^{-E}`` — a *lower bound*
        # whose quoted value is itself a range (e.g. Marsh & Niemeyer 2019:
        # ``m \gtrsim 0.6 \to 1\times10^{-19}\mathrm{eV}``).  It is NOT an
        # allowed band between a and b: the allowed region is m ≳ a, and
        # [a, b] only marks the spread of the quoted limit.  We therefore emit
        # an ``include_above`` interval anchored at the *weakest* value a and
        # record the stronger end b as an uncertainty band.  Must be tested
        # before Pattern F, which also accepts ``\to`` as a range dash.
        m = re.search(
            rf"m\s*(?:>=|>|≳|≥)\s*({_COEF_PLAIN})\s*(?:\\to|->|→|\\sim|~)\s*"
            rf"({_COEF_PLAIN})\s*(?:\\times|×|x)\s*{_POWER}",
            norm,
        )
        if m:
            lo = _log10(float(m.group(1))) - float(m.group("exp"))
            band = _log10(float(m.group(2))) - float(m.group("exp"))
            if band > lo:
                out.append(
                    dict(
                        type="include_above",
                        log10_min=lo,
                        log10_max=math.inf,
                        log10_band_max=band,
                        label=f"include m/eV ≥ 10^{lo:.2f} (quoted as 10^{lo:.2f}–10^{band:.2f})",
                    )
                )
                continue

        # Pattern F (reordered): a \\sim b \\times 10^{-e}  (range with shared exponent)
        m = re.search(
            rf"({_COEF_PLAIN})\s*(?:\\sim|\\to|~)\s*({_COEF_PLAIN})\s*(?:\\times|×|x)\s*{_POWER}",
            norm,
        )
        if m:
            a_val = _log10(float(m.group(1))) - float(m.group("exp"))
            b_val = _log10(float(m.group(2))) - float(m.group("exp"))
            lo, hi = sorted([a_val, b_val])
            t = "exclude_range" if is_exclude else "include_range"
            out.append(
                dict(
                    type=t,
                    log10_min=lo,
                    log10_max=hi,
                    label=f"10^{lo:.2f} ≲ m/eV ≲ 10^{hi:.2f}",
                )
            )
            continue

        # Pattern C: $m > b$ / $m \ge b$ / $m \gtrsim b$
        m = re.search(
            rf"m\s*(?:>=|>|~|≳|≥)\s*({_VAL}|{_COEF_PLAIN})",
            norm,
        )
        if m:
            b = _value_match_to_log10(m.group(1).strip())
            if b is not None:
                t = "exclude_above" if is_exclude else "include_above"
                out.append(
                    dict(
                        type=t,
                        log10_min=b,
                        log10_max=math.inf,
                        label=f"{'exclude' if is_exclude else 'include'} m/eV ≥ 10^{b:.2f}",
                    )
                )
                continue

        # Pattern D: $m < a$ / $m \le a$ / $m \lesssim a$
        m = re.search(
            rf"m\s*(?:<=|~|≲|<|≤)\s*({_VAL}|{_COEF_PLAIN})",
            norm,
        )
        if m:
            a = _value_match_to_log10(m.group(1).strip())
            if a is not None:
                t = "exclude_below" if is_exclude else "include_below"
                out.append(
                    dict(
                        type=t,
                        log10_min=-math.inf,
                        log10_max=a,
                        label=f"{'exclude' if is_exclude else 'include'} m/eV ≤ 10^{a:.2f}",
                    )
                )
                continue

        # Pattern G: 10^{-X} \\sim 10^{-Y}   (range with pure powers of ten)
        m = re.search(
            rf"(10\^\{{-[0-9]+(?:\.[0-9]+)?\}})\s*\\sim\s*(10\^\{{-[0-9]+(?:\.[0-9]+)?\}})",
            norm,
        )
        if m:
            exp_a = _value_match_to_log10(m.group(1))
            exp_b = _value_match_to_log10(m.group(2))
            if exp_a is not None and exp_b is not None:
                lo, hi = sorted([exp_a, exp_b])
                t = "exclude_range" if is_exclude else "include_range"
                out.append(
                    dict(
                        type=t,
                        log10_min=lo,
                        log10_max=hi,
                        label=f"10^{lo:.2f} ≲ m/eV ≲ 10^{hi:.2f}",
                    )
                )
                continue

        # Pattern H: 10^{-X} \\le m \\le 10^{-Y}
        m = re.search(
            rf"(10\^\{{-[0-9]+(?:\.[0-9]+)?\}})\s*(?:<=|<|≥|>)\s*m\s*(?:<=|<|≥|>)\s*(10\^\{{-[0-9]+(?:\.[0-9]+)?\}})",
            norm,
        )
        if m:
            exp_a = _value_match_to_log10(m.group(1))
            exp_b = _value_match_to_log10(m.group(2))
            if exp_a is not None and exp_b is not None:
                lo, hi = sorted([exp_a, exp_b])
                t = "exclude_range" if is_exclude else "include_range"
                out.append(
                    dict(
                        type=t,
                        log10_min=lo,
                        log10_max=hi,
                        label=f"10^{lo:.2f} ≲ m/eV ≲ 10^{hi:.2f}",
                    )
                )
                continue

        # Pattern I: a \\lesssim b \\times 10^{-E}  (upper bound only — treated as below)
        m = re.search(
            rf"m\s*\\lesssim\s*({_VAL_FULL})",
            norm,
        )
        if m:
            a = _value_match_to_log10(m.group(1).strip())
            if a is not None:
                t = "exclude_below" if is_exclude else "include_below"
                out.append(
                    dict(
                        type=t,
                        log10_min=-math.inf,
                        log10_max=a,
                        label=f"{'exclude' if is_exclude else 'include'} m/eV ≤ 10^{a:.2f}",
                    )
                )
                continue

        # Pattern J: 'constrain m ≤ b: ρ<6 GeV/cm³' — extract the mass bound only.
        m = re.search(rf"constrain\s+m\s*(?:<=|<|≥|>)\s*({_VAL_FULL})", norm)
        if m:
            a = _value_match_to_log10(m.group(1).strip())
            if a is not None:
                t = "include_below" if not is_exclude else "exclude_below"
                out.append(
                    dict(
                        type=t,
                        log10_min=-math.inf,
                        log10_max=a,
                        label=f"m/eV ≤ 10^{a:.2f}",
                    )
                )
                continue

        # Pattern K: 'constrain a ≤ m ≤ b: Ω_a/Ω_d ≤ 0.05'
        m = re.search(
            rf"constrain\s+({_VAL}|10\^\{{-[0-9]+\}})\s*(?:<=|<|≥|>)\s*m\s*(?:<=|<|≥|>)\s*({_VAL}|10\^\{{-[0-9]+\}})",
            norm,
        )
        if m:
            a = _value_match_to_log10(m.group(1).strip())
            b = _value_match_to_log10(m.group(2).strip())
            if a is not None and b is not None:
                lo, hi = sorted([a, b])
                out.append(
                    dict(
                        type="include_range",
                        log10_min=lo,
                        log10_max=hi,
                        label=f"10^{lo:.2f} ≤ m/eV ≤ 10^{hi:.2f}",
                    )
                )
                continue

    # De-duplicate.
    seen = set()
    uniq: list[dict] = []
    for c in out:
        key = (c["type"], round(c["log10_min"], 4), round(c["log10_max"], 4))
        if key in seen:
            continue
        seen.add(key)
        uniq.append(c)

    return uniq


# --------------------------------------------------------------------------- #
# TeX table parsing                                                           #
# --------------------------------------------------------------------------- #


@dataclass
class TableRow:
    citekey: str
    data: str
    method: str
    raw_constraint: str
    category: str = ""
    # Position of the row in the .tex table (0-based).  The figures use it to
    # keep their vertical order identical to the table the reader sees.
    table_index: int = 0


_CATEGORY_KEYWORDS: list[tuple[str, list[str]]] = [
    ("galaxy",  ["rotation curve", "spiral galaxy", "m31", "m33", "mw rotation", "sparc", "stellar streams", "ngc", "leo ii", "draco", "crater", "fornax", "sextans", "eridanus", "carina", "segue", "ufds", "dsph", "stellar kinematics", "dwarf"]),
    ("bh",      ["bh superradiance", "bh spin", "superradiance", "m87*", "sgr a*", "sgr", "eht", "event horizon", "black hole"]),
    ("lensing", ["lensing", "weak lensing", "cosmic-shear", "subhalo", "shmf", "strong lens", "hmf", "concentration"]),
    ("pta",     ["pulsar timing", "ppak", "parkes", "pulsar", "ppak"]),
    ("cmb",     ["cmb", "bbn", "lss", "wigglez", "kde"]),
    ("lya",     ["lyman", "ly$\\alpha$", "lyalpha", "lya"]),
    ("21cm",    ["21cm", "edges", "cosmic dawn", "21-cm"]),
    ("heating", ["heating", "granule", "disk heating", "dynamical heating"]),
]


def _classify(data: str, method: str, raw: str) -> str:
    blob = " ".join([data.lower(), method.lower(), raw.lower()])
    for cat, kws in _CATEGORY_KEYWORDS:
        for kw in kws:
            if kw in blob:
                return cat
    return "other"


def _split_tex_row(body: str) -> list[str]:
    """Split a TeX table row into cells on *unescaped* ``&`` only.

    ``\\&`` is an escaped ampersand — literal text inside a cell (e.g. the
    ``Segue I \\& II UFD`` / ``MW dSphs \\& UFDs`` data cells).  Splitting on
    every ``&`` silently swapped the ``data`` / ``method`` columns of those
    rows and leaked the method text into the constraint cell.

    The escaped ampersands are restored for ``data`` / ``method`` (they are
    part of the printed text), but deliberately *kept* escaped in the
    constraint cell, where ``\\&`` separates two intervals of one work
    (Stott 2018: ``7×10⁻¹⁴–2×10⁻¹¹ eV \\& 7×10⁻²⁰–1×10⁻¹⁶ eV``).
    """
    cells = [c.strip() for c in re.split(r"(?<!\\)&", body)]
    # The row body starts with the ``&`` that ends the \cite cell, hence the
    # leading empty cell.  Trailing empties come from trailing whitespace.
    cells = [c for c in cells if c]
    return [c.replace("\\&", "&") if i < len(cells) - 1 else c
            for i, c in enumerate(cells)]


def parse_table(path: Path) -> list[TableRow]:
    text = path.read_text(encoding="utf-8")
    rows: list[TableRow] = []
    pattern = re.compile(
        r"\\cite\{(?P<key>[^}]+)\}\s*(?P<body>.*?)\\\\",
        re.DOTALL,
    )
    for m in pattern.finditer(text):
        cells = _split_tex_row(m.group("body"))
        if len(cells) != 3:
            print(
                f"[!] Skipping malformed row for {m.group('key').strip()} "
                f"({len(cells)} cells instead of 3)",
                flush=True,
            )
            continue
        data, method, limit = cells
        rows.append(
            TableRow(
                citekey=m.group("key").strip(),
                data=re.sub(r"\s+", " ", data).strip(),
                method=re.sub(r"\s+", " ", method).strip(),
                raw_constraint=limit.strip(),
            )
        )
    for idx, row in enumerate(rows):
        row.table_index = idx
        row.category = _classify(row.data, row.method, row.raw_constraint)
    return rows


# --------------------------------------------------------------------------- #
# Driver                                                                      #
# --------------------------------------------------------------------------- #


def main() -> None:
    if not TEX_FILE.exists():
        sys.exit(f"Cannot find {TEX_FILE}")
    if not BIB_FILE.exists():
        sys.exit(f"Cannot find {BIB_FILE}")

    print(f"[i] Parsing {TEX_FILE} …", flush=True)
    rows = parse_table(TEX_FILE)
    print(f"    Found {len(rows)} constraint rows.", flush=True)

    print(f"[i] Parsing {BIB_FILE} …", flush=True)
    bib = parse_bib(BIB_FILE)
    print(f"    Found {len(bib)} bibliography entries.", flush=True)

    out: dict[str, dict] = {}
    for row in rows:
        intervals = parse_constraint(row.raw_constraint)
        b = bib.get(row.citekey)
        out[row.citekey] = dict(
            citekey=row.citekey,
            data=row.data,
            method=row.method,
            constraint_raw=row.raw_constraint,
            category=row.category,
            table_index=row.table_index,
            constraints=intervals,
            bib=(
                dict(
                    title=b.title if b else "",
                    authors=b.authors if b else [],
                    year=b.year if b else 0,
                    journal=b.journal if b else "",
                    doi=b.doi if b else "",
                    eprint=b.eprint if b else "",
                    url=b.url if b else "",
                )
                if b is not None
                else dict(title="", authors=[], year=0, journal="", doi="", eprint="", url="")
            ),
        )

    with open(YAML_OUT, "w", encoding="utf-8") as f:
        yaml.safe_dump(out, f, allow_unicode=True, sort_keys=True, width=200)
    with open(JSON_OUT, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)

    print(f"[i] Wrote {YAML_OUT}", flush=True)
    print(f"[i] Wrote {JSON_OUT}", flush=True)
    bad = [k for k, v in out.items() if not v["constraints"]]
    if bad:
        print(f"[!] WARNING: {len(bad)} entries produced no parsed intervals:", flush=True)
        for k in bad:
            print(f"    - {k}: {out[k]['constraint_raw']}", flush=True)
    else:
        print("[OK] All entries produced at least one parsed interval.", flush=True)


if __name__ == "__main__":
    main()