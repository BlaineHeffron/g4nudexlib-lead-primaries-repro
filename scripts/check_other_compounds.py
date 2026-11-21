"""Compare NuDEX's strongest primaries with ENSDF and EGAF for non-lead compounds flagged by a survey.

Usage: python3 scripts/check_other_compounds.py PRIMARY_CAPTURE_GAMMAS_DAT ENSDF_ZIP NGLIST_E_DAT

For each compound: NuDEX per capture (library % / 100); ENSDF per capture (RI * NR * BR / 100, gammas placed from
the capture state of the named thermal dataset); EGAF per capture as sigma_gamma / sigma0, with sigma0 read from
the block header, sigma0 = 100 / (4th number * theta), theta the natural abundance of the target in Geant4
(G4NistElementBuilder.cc). Columns 'rel' divide each source by its own value for the reference line (the strongest
NuDEX line present in both ENSDF and EGAF), so they compare shapes independently of any normalisation.
"""
import re
import sys
import zipfile

# product ZA, target abundance theta (Geant4 G4NistElementBuilder), ENSDF dataset ID
CASES = [
    (30068, 0.041, "67ZN(N,G) E=THERMAL"),
    (32073, 0.2754, "72GE(N,G) E=THERMAL"),
    (44100, 0.1276, "99RU(N,G) E=TH"),
    (44103, 0.3155, "102RU(N,G) E=TH"),
    (50116, 0.0034, "115SN(N,G) E=TH"),
    (52123, 0.0255, "122TE(N,G) E=TH"),
    (52125, 0.0474, "124TE(N,G) E=THERMAL"),
    (55134, 1.0, "133CS(N,G) E=THERMAL"),
    (78197, 0.2524, "196PT(N,G) E=THERMAL"),
]


def num(s):
    s = s.strip()
    return float(s) if re.fullmatch(r"[-+]?\d*\.?\d+(E[-+]?\d+)?", s) else None


def main():
    lib = open(sys.argv[1]).read().splitlines()
    zf = zipfile.ZipFile(sys.argv[2])
    egaf = [l.split() for l in open(sys.argv[3]).read().splitlines()[1:] if len(l.split()) >= 6]
    for za, theta, dsid in CASES:
        z, a = za // 1000, za % 1000
        i = next(k for k, l in enumerate(lib) if l.split()[:2] == ["ZA", str(za)])
        head = lib[i + 1].split()
        prim = [tuple(map(float, lib[i + 2 + k].split()[:2])) for k in range(int(head[1]))]
        sigma0 = 100 / (float(head[3]) * theta)
        nucid = lib[i].split()[2]
        text = zf.read(f"ensdf.{a:03d}").decode("latin-1").splitlines()
        k = next(j for j, l in enumerate(text) if l[:5].strip() == nucid and l[5:8] == "   " and l[9:39].strip() == dsid)
        block = []
        for l in text[k:]:
            if not l.strip():
                break
            block.append(l)
        nr = br = 1.0
        for l in block:
            if l[5:8] == "  N":
                nr, br = num(l[9:19]) or 1.0, num(l[31:39]) or 1.0
        level, levels, gam = None, [], []
        for l in block:
            if l[5:8] == "  L":
                level = num(l[9:19])
                levels.append(level)
            elif l[5:8] == "  G":
                e, ri = num(l[9:19]), num(l[21:29])
                if e is not None and ri is not None:
                    gam.append((e, ri * nr * br / 100, level))
        cap = max(v for v in levels if v is not None)
        ens = lambda e: sum(y for ee, y, lv in gam if lv == cap and abs(ee - e) <= 1.0)
        eg = lambda e: sum(float(r[4]) for r in egaf if r[1] == str(z) and int(re.match(r"\d+", r[0]).group()) == a
                           and abs(float(r[2]) - e) <= 1.0) / sigma0
        top = sorted(prim, key=lambda p: -p[1])[:6]
        ref = next((e for e, _ in top if ens(e) and eg(e)), None)
        r0 = {s: (f(ref) if ref else 0) for s, f in (("e", ens), ("g", eg))}
        y_ref = next(y for e, y in prim if e == ref) / 100 if ref else 0
        print(f"\n{za} {lib[i].split()[2]} {' '.join(lib[i].split()[3:])}  sigma0(header)={sigma0:.4g} b; "
              f"ENSDF '{dsid}' NR={nr}; reference line {ref} keV")
        print(f"  {'E_keV':>9} {'NuDEX':>8} {'ENSDF':>8} {'EGAF':>8} | {'rel NuDEX':>9} {'rel ENSDF':>9} {'rel EGAF':>9}")
        for e, y in top:
            y /= 100
            en, g = ens(e), eg(e)
            rel = lambda v, v0: f"{v / v0:9.3g}" if v0 and v else f"{'-':>9}"
            print(f"  {e:9.2f} {y:8.4g} {en:8.4g} {g:8.4g} | {rel(y, y_ref)} {rel(en, r0['e'])} {rel(g, r0['g'])}")


if __name__ == "__main__":
    main()
