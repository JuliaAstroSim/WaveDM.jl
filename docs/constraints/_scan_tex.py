"""Scan ma_constraints.json for raw TeX characters in data/method fields."""
import json
import re

d = json.load(open(r'e:\JuliaAstroSim\WaveDM.jl\docs\constraints\output\ma_constraints.json', encoding='utf-8'))
for k, v in d.items():
    for f in ['data', 'method']:
        s = v.get(f, '') or ''
        if '$' in s or re.search(r'\\[a-zA-Z]+', s):
            print(f'{k}.{f}: {s!r}')