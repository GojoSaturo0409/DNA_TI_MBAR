#!/usr/bin/env python3
"""Round every 'value +- uncertainty' in LaTeX files to the editor's rule: uncertainty to one significant figure
(two if its leading digit is 1), value to the same decimal place. Handles 'a $\\pm$ b', 'a\\pm b', '$-$a $\\pm$ b'.
Usage: round_sigfig.py file.tex [...]   (edits in place, prints changed pairs)"""
import re, sys, math
from decimal import Decimal, ROUND_HALF_UP
PAT = re.compile(r"(?P<sign>\$-\$|-)?(?P<val>\d+(?:\.\d+)?)(?P<sep>\s*(?:\$\\pm\$|\\pm)\s*)(?P<err>\d+(?:\.\d+)?)")
def rnd(x, nd):
    q = Decimal(1).scaleb(-nd) if nd > 0 else Decimal(1)
    d = Decimal(x).quantize(q, rounding=ROUND_HALF_UP) if nd >= 0 else (Decimal(x) / Decimal(10) ** (-nd)).quantize(Decimal(1), rounding=ROUND_HALF_UP) * Decimal(10) ** (-nd)
    return d
def fix(m):
    v, e = m.group("val"), m.group("err"); ef = float(e)
    if ef == 0: return m.group(0)
    p = -int(math.floor(math.log10(ef)))           # decimal place of leading digit
    lead = int(Decimal(e).scaleb(p))                # leading digit
    nd = p + 1 if lead == 1 else p                  # keep two figures if leading digit is 1
    er = rnd(e, nd)
    if er != 0 and er.adjusted() > Decimal(e).adjusted():   # rounding carried (e.g. 0.96 -> 1.0): recompute
        lead2 = int(er.scaleb(-er.adjusted())); nd = -er.adjusted() + (1 if lead2 == 1 else 0); er = rnd(e, nd)
    vr = rnd(v, nd)
    fmt = lambda d: (f"{d:.{nd}f}" if nd > 0 else f"{int(d)}")
    return f"{m.group('sign') or ''}{fmt(vr)}{m.group('sep')}{fmt(er)}"
for f in sys.argv[1:]:
    s = open(f).read(); out = PAT.sub(fix, s)
    if out != s:
        for a, b in zip(PAT.finditer(s), PAT.finditer(out)):
            if a.group(0) != b.group(0): print(f"{f}: {a.group(0)!r} -> {b.group(0)!r}")
        open(f, "w").write(out)
