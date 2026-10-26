"""
constraint_labels.py
====================

Shared text helpers for the three ``ma_constraints`` renderers (static
PNG/PDF, clickable SVG, interactive HTML).

The figure is a visual companion to ``ma_constraints.tex``, so the row labels,
the TeX cleaning and the abbreviation table must be *identical* across all
output formats — keeping them in one module prevents the three renderers from
drifting apart.

Contents
--------
* :func:`strip_tex` — reduce inline-math / BibTeX noise to plain Unicode
  (handles accents written both as ``\\'{e}`` and as ``\\v si\\v c``, i.e. with
  the space that BibTeX allows between an accent macro and its letter).
* :func:`abbreviate` — the curated ``data`` / ``method`` short forms used in
  the figures (e.g. ``rotation curve`` → ``RC``), so that the
  ``data — method`` annotation of every row fits next to its bar without
  running off the axes.
* :func:`assign_row_labels` — ``FirstAuthor [et al.] Year`` labels, with
  ``a`` / ``b`` suffixes when two *different* papers by the same first author
  share a year (Bar et al. 2019a/2019b, Maleki et al. 2020a/2020b).
"""

from __future__ import annotations

import re

# --------------------------------------------------------------------------- #
# TeX → Unicode                                                               #
# --------------------------------------------------------------------------- #

# (accent macro, base letter) -> accented character.  Covers the full BibTeX
# accent set actually used by the bibliography (acute, grave, umlaut, tilde,
# circumflex, cedilla, double-acute, macron, dot, caron, breve).
ACCENT_MAP: dict[tuple[str, str], str] = {
    ("'", "a"): "á", ("'", "A"): "Á", ("'", "c"): "ć", ("'", "C"): "Ć",
    ("'", "e"): "é", ("'", "E"): "É", ("'", "i"): "í", ("'", "I"): "Í",
    ("'", "n"): "ń", ("'", "N"): "Ń", ("'", "o"): "ó", ("'", "O"): "Ó",
    ("'", "s"): "ś", ("'", "S"): "Ś", ("'", "u"): "ú", ("'", "U"): "Ú",
    ("'", "y"): "ý", ("'", "Y"): "Ý", ("'", "z"): "ź", ("'", "Z"): "Ź",
    ("`", "a"): "à", ("`", "A"): "À", ("`", "e"): "è", ("`", "E"): "È",
    ("`", "i"): "ì", ("`", "I"): "Ì", ("`", "o"): "ò", ("`", "O"): "Ò",
    ("`", "u"): "ù", ("`", "U"): "Ù",
    ('"', "a"): "ä", ('"', "A"): "Ä", ('"', "e"): "ë", ('"', "E"): "Ë",
    ('"', "i"): "ï", ('"', "I"): "Ï", ('"', "o"): "ö", ('"', "O"): "Ö",
    ('"', "u"): "ü", ('"', "U"): "Ü",
    ("~", "a"): "ã", ("~", "A"): "Ã", ("~", "n"): "ñ", ("~", "N"): "Ñ",
    ("~", "o"): "õ", ("~", "O"): "Õ",
    ("^", "a"): "â", ("^", "A"): "Â", ("^", "e"): "ê", ("^", "E"): "Ê",
    ("^", "i"): "î", ("^", "I"): "Î", ("^", "o"): "ô", ("^", "O"): "Ô",
    ("^", "u"): "û", ("^", "U"): "Û",
    ("c", "c"): "ç", ("c", "C"): "Ç", ("c", "s"): "ş", ("c", "S"): "Ş",
    ("c", "t"): "ţ", ("c", "T"): "Ţ",
    ("H", "o"): "ő", ("H", "O"): "Ő", ("H", "u"): "ű", ("H", "U"): "Ű",
    ("=", "a"): "ā", ("=", "A"): "Ā", ("=", "e"): "ē", ("=", "E"): "Ē",
    ("=", "i"): "ī", ("=", "I"): "Ī", ("=", "o"): "ō", ("=", "O"): "Ō",
    ("=", "u"): "ū", ("=", "U"): "Ū",
    (".", "z"): "ż", (".", "Z"): "Ż",
    ("v", "c"): "č", ("v", "C"): "Č", ("v", "n"): "ň", ("v", "N"): "Ň",
    ("v", "s"): "š", ("v", "S"): "Š", ("v", "z"): "ž", ("v", "Z"): "Ž",
    ("u", "o"): "ŏ", ("u", "O"): "Ŏ", ("u", "u"): "ŭ", ("u", "U"): "Ŭ",
}

# Stand-alone letter macros (not accents) that show up in author lists.
LETTER_MACROS: dict[str, str] = {
    "\\ss": "ß", "\\aa": "å", "\\AA": "Å", "\\ae": "æ", "\\AE": "Æ",
    "\\oe": "œ", "\\OE": "Œ", "\\o": "ø", "\\O": "Ø", "\\l": "ł", "\\L": "Ł",
    "\\i": "i", "\\j": "j",
}

# Relation / symbol macros that appear in the ``data`` / ``method`` cells.
SYMBOL_MAP: dict[str, str] = {
    "\\gtrsim": "≥", "\\lesssim": "≤", "\\geq": "≥", "\\leq": "≤",
    "\\ge": "≥", "\\le": "≤", "\\sim": "~", "\\to": "→", "\\approx": "≈",
    "\\times": "×", "\\cdot": "·", "\\ldots": "…", "\\pm": "±",
    "\\alpha": "α", "\\beta": "β", "\\gamma": "γ", "\\delta": "δ",
    "\\epsilon": "ε", "\\sigma": "σ", "\\mu": "μ", "\\nu": "ν",
    "\\Omega": "Ω", "\\rho": "ρ", "\\Lambda": "Λ",
}


def strip_tex(s: str | None) -> str:
    """Reduce TeX / BibTeX noise to plain Unicode text.

    Handles ``$...$`` delimiters, ``\\mathrm{...}`` unit wrappers, accent
    macros (both ``\\'{e}`` and ``\\v si\\v c`` forms), letter macros
    (``\\l`` → ł), super/subscripts (``M87^*`` → ``M87*``), relation macros
    (``\\leq`` → ≤) and the TeX non-breaking space ``~`` used inside names
    such as ``NGC~6240`` (only outside math, so ``$z\\sim15$`` keeps its
    "approximately" tilde).
    """
    if not s:
        return ""
    # 1. ``~`` is a non-breaking space in text but may be a "≈" glyph inside
    #    ``$...$`` (from \sim).  Replace it only in the text segments.
    parts = s.split("$")
    for i in range(0, len(parts), 2):
        parts[i] = parts[i].replace("~", " ")
    s = "".join(parts)
    # 2. ``\mathrm{eV}`` / ``\text{...}`` → contents.
    s = re.sub(r"\\(?:mathrm|text|ensuremath)\{([^}]*)\}", r"\1", s)
    # 3. Accent macros.  BibTeX allows a space between the macro and its
    #    letter (``Ir\v si\v c``), which the naive ``\v s`` pattern misses.
    for (acc, base), repl in ACCENT_MAP.items():
        s = re.sub(r"\\" + re.escape(acc) + r"\s*\{" + re.escape(base) + r"\}",
                   repl, s)
        s = re.sub(r"\\" + re.escape(acc) + r"\s*" + re.escape(base), repl, s)
    # 4. Letter macros.
    for macro, repl in LETTER_MACROS.items():
        s = re.sub(re.escape(macro) + r"(?![a-zA-Z])", repl, s)
    # 5. ``\^{...}`` / ``\_{...}`` → contents; bare ``\^x`` / ``\_x`` → x.
    s = re.sub(r"\\[\^_]\{([^}]*)\}", r"\1", s)
    s = re.sub(r"\\[\^_](\S)", r"\1", s)
    # 6. Any ``\X{...}`` accent we did not map: keep the base letter.
    s = re.sub(r"\\[`'\"~=\^\.vucHkdHrtb]\{([^}]*)\}", r"\1", s)
    # 7. Bare super/subscript markers: ``M87^*`` → ``M87*``, ``x_a`` → ``xa``.
    #    Only the marker is dropped — the ``*`` of ``M87*`` / ``Sgr A*`` is a
    #    meaningful part of the source name.
    s = re.sub(r"[\^_]", "", s)
    # 8. Relation / symbol macros.
    for macro, repl in SYMBOL_MAP.items():
        s = s.replace(macro, repl)
    # 9. Braces, leftover backslashes, en-dash.
    s = re.sub(r"\{([^{}]*)\}", r"\1", s)
    s = s.replace("\\", "")
    s = s.replace("--", "-")
    s = re.sub(r"\s+", " ", s).strip()
    # 10. Trailing punctuation that comes from the source table cells.
    s = re.sub(r"[,;]+$", "", s).strip()
    s = re.sub(r",\s*\.$", "", s).strip()
    return s


# --------------------------------------------------------------------------- #
# Abbreviations used in the figures                                            #
# --------------------------------------------------------------------------- #
# Applied in order to the *stripped* label text.  Only standard, unambiguous
# short forms are used (the ones a reader of the m_a literature already knows);
# anything that would need a glossary is left alone.
ABBREV_RULES: list[tuple[str, str]] = [
    (r"\bLyman-\s*α\b", "Lyα"),
    (r"\bflux power spec\.", "P(k)"),
    (r"\bpower spectrum\b", "P(k)"),
    (r"\bsoliton-halo relation\b", "SH relation"),
    (r"\bsubhalo mass function\b", "SHMF"),
    (r"\bsubhalo function\b", "SHMF"),
    (r"\brotation curves\b", "RCs"),
    (r"\brotation curve\b", "RC"),
    (r"\bEDGES High-Band", "EDGES HB"),
    (r"\bEridanus II\b", "Eri II"),
    (r"\bstellar kinematics\b", "stellar kin."),
    (r"\bdisk-star kinematics\b", "disk-star kin."),
    (r"\bgranule fluctuations\b", "granule fluct."),
    (r"\bglobular-cluster timing\b", "GC timing"),
    (r"\bconcentration-mass relation\b", "c-M relation"),
    (r"\bSPARC LSB galaxies\b", "SPARC LSBs"),
    (r"\bstatistical search\b", "stat. search"),
    (r"\bstar-cluster heating\b", "SC heating"),
    (r"\banalytical\b", "anal."),
]


def abbreviate(text: str) -> str:
    """Apply the curated abbreviation table to an already-stripped label."""
    out = text
    for pattern, repl in ABBREV_RULES:
        # ``lambda`` replacement: the short forms contain ``$...$`` math and
        # backslashes that must survive verbatim (``re`` would otherwise
        # interpret e.g. ``\s`` in the template as an escape).
        out = re.sub(pattern, lambda _m, r=repl: r, out, flags=re.IGNORECASE)
    out = re.sub(r"\s+", " ", out).strip()
    out = re.sub(r"\s+([,;)])", r"\1", out)
    return out


def data_label(row: dict) -> str:
    """Plain, abbreviated ``data`` label of a table row."""
    return abbreviate(strip_tex(row.get("data") or ""))


def method_label(row: dict) -> str:
    """Plain, abbreviated ``method`` label of a table row."""
    return abbreviate(strip_tex(row.get("method") or ""))


# --------------------------------------------------------------------------- #
# Row labels                                                                  #
# --------------------------------------------------------------------------- #

def short_author(name: str) -> str:
    """Family name of the first author, TeX-cleaned (``Ir\\v si\\v c`` → Iršić)."""
    return strip_tex(name)


def row_label(row: dict) -> str:
    """``FirstAuthor et al. Year`` (or ``FirstAuthor Year`` for single authors)."""
    bib = row.get("bib", {}) or {}
    year = bib.get("year") or "?"
    authors = bib.get("authors") or []
    if not authors:
        return f"Unknown {year}"
    first = short_author(authors[0].split(",")[0]).strip()
    if len(authors) > 1:
        return f"{first} et al. {year}"
    return f"{first} {year}"


def assign_row_labels(rows: list[dict]) -> dict[str, str]:
    """Map citekey → row label, disambiguating different papers by a/b.

    ``rows`` must be in display order; the first citekey that shares an
    ``Author et al. Year`` label gets the ``a`` suffix, the next ``b``, etc.
    Two *separate* papers that truncate to the same label (e.g. Bar et al.
    2019: "Looking for ULDM near SMBHs" and "ULDM in disk galaxies") stay
    distinguishable, while the multiple intervals of a *single* paper share
    one row and therefore one label.
    """
    groups: dict[str, list[str]] = {}
    for row in rows:
        groups.setdefault(row_label(row), []).append(row["citekey"])
    labels: dict[str, str] = {}
    for label, keys in groups.items():
        if len(keys) == 1:
            labels[keys[0]] = label
        else:
            for i, key in enumerate(keys):
                labels[key] = f"{label}{chr(ord('a') + i)}"
    return labels