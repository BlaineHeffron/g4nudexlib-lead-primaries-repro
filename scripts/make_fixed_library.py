"""Make a G4NUDEXLIB tree with the three lead blocks of PrimaryCaptureGammas.dat put on a per-capture basis.

Usage: python3 scripts/make_fixed_library.py G4NUDEXLIB1.0_DIR reference/lead_reference.tsv OUT_DIR

OUT_DIR gets symlinks to every entry of the shipped tree except PrimaryCaptureGammas.dat, which is rewritten:
- ZA 82205 (204Pb + n): the block from the shipped PrimaryCaptureGammas01.dat (ENSDF "204PB(N,G) E=THERMAL",
  NR = 0.1), which is already per capture.
- ZA 82207 (206Pb + n) and ZA 82208 (207Pb + n): the same primary energies, with intensities replaced by the EGAF
  per-capture values in the reference table (sigma_gamma / sigma_E, energy-normalised).
The header line keeps Sn, sets the line count and sum, and writes "1 0" for the two trailing numbers, which NuDEX
does not read.
"""
import os
import sys

src, ref_path, out = sys.argv[1:4]
os.makedirs(out)
for name in os.listdir(src):
    if name != "PrimaryCaptureGammas.dat":
        os.symlink(os.path.abspath(os.path.join(src, name)), os.path.join(out, name))

egaf = {}
for l in open(ref_path).read().splitlines()[1:]:
    prod, e, s, y = l.split("\t")[:4]
    if s == "EGAF":
        egaf.setdefault({"207PB": 82207, "208PB": 82208}.get(prod), []).append((float(e), float(y)))


def blocks(path):
    lines, i, out = open(path).read().splitlines(), 0, {}
    while i < len(lines):
        f = lines[i].split()
        if f and f[0] == "ZA":
            n = int(lines[i + 1].split()[1])
            out[int(f[1])] = (i, i + 2 + n, lines[i:i + 2 + n])
            i += 2 + n
        else:
            i += 1
    return lines, out


lines, idx = blocks(os.path.join(src, "PrimaryCaptureGammas.dat"))
_, idx01 = blocks(os.path.join(src, "PrimaryCaptureGammas01.dat"))
new = {82205: idx01[82205][2]}
for za in (82207, 82208):
    head, sn = idx[za][2][0], idx[za][2][1].split()[0]
    prim = []
    for row in idx[za][2][2:]:
        e = float(row.split()[0])
        match = [y for ee, y in egaf[za] if abs(ee - e) <= 1.0]
        assert len(match) == 1, (za, e, match)
        prim.append((e, 100 * match[0]))
    new[za] = [head, f"{sn}  {len(prim)}  {sum(y for _, y in prim):.6g}  1  0"] + \
              [f"      {e:g}   {y:.6g}" for e, y in prim]
for za in sorted(new, key=lambda z: -idx[z][0]):  # replace from the bottom so earlier indices stay valid
    a, b, _ = idx[za]
    lines[a:b] = new[za]
open(os.path.join(out, "PrimaryCaptureGammas.dat"), "w").write("\n".join(lines) + "\n")
for za in sorted(new):
    print("\n".join(new[za]))
