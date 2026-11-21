"""Extract per-capture primary-gamma intensities for the three lead compounds from ENSDF and EGAF.

Usage: python3 scripts/extract_reference.py ENSDF_ZIP NGLIST_E_DAT > reference/lead_reference.tsv

ENSDF_ZIP     NNDC ENSDF archive (zip of ensdf.001 ... ensdf.300).
NGLIST_E_DAT  EGAF isotopic list (columns: product, Z, E keV, dE, sigma_gamma b, d sigma_gamma b).

ENSDF: intensity per capture = RI * NR * BR / 100 for gammas placed from the capture state (the highest level).
Where the dataset has no NR record (208Pb, "207PB(N,G)"), it is normalised by energy: sum(E I) = Sn.
EGAF: per capture = sigma_gamma / sigma_E with sigma_E = sum(E sigma_gamma) / Sn over every listed line of
the product (energy normalisation, independent of any sigma0). The 206Pb + n compound has no thermal ENSDF
dataset (only E=RES), so only EGAF is given for it.
"""
import re
import sys
import zipfile

# product nucid, ENSDF dataset ID (exact); Sn is the ENSDF capture level, or for 207Pb (no thermal ENSDF
# dataset) the strongest EGAF primary plus its recoil, E + E^2 / (2 M c^2), the ground-state transition.
CASES = [("205PB", "204PB(N,G) E=THERMAL"), ("207PB", None), ("208PB", "207PB(N,G)")]
AMU_KEV = 931494.1


def num(s):
    s = s.strip()
    return float(s) if re.fullmatch(r"[-+]?\d*\.?\d+(E[-+]?\d+)?", s) else None


def ensdf(zf, nucid, dsid):
    a = int(nucid[:3])
    text = zf.read(f"ensdf.{a:03d}").decode("latin-1").splitlines()
    starts = [i for i, l in enumerate(text)
              if l[:5].strip() == nucid and l[5:8] == "   " and l[9:39].strip() == dsid]
    assert len(starts) == 1, (nucid, dsid, len(starts))
    block = []
    for l in text[starts[0]:]:
        if not l.strip():
            break
        block.append(l)
    nr = br = None
    for l in block:
        if l[5:8] == "  N":
            nr, br = num(l[9:19]), num(l[31:39])
    level, levels, gam = None, [], []
    for l in block:
        if l[5:8] == "  L":
            level = num(l[9:19])
            levels.append(level)
        elif l[5:8] == "  G":
            e, ri = num(l[9:19]), num(l[21:29])
            if e is not None and ri is not None:
                gam.append((e, ri, level))
    sn = max(v for v in levels if v is not None)
    if nr:
        scale, how = nr * (br or 1.0) / 100, f"NR={nr}"
    else:
        scale, how = sn / sum(e * ri for e, ri, _ in gam), "energy"
    return sn, how, block[0][9:39].strip(), {e: ri * scale for e, ri, lev in gam if lev == sn}


def egaf(path, nucid):
    rows = []
    for l in open(path).read().splitlines()[1:]:
        f = l.split()
        if len(f) >= 6 and f[0] == nucid:
            rows.append((float(f[2]), float(f[4])))
    return rows


def main():
    zf = zipfile.ZipFile(sys.argv[1])
    print("product\tline_keV\tsource\tper_capture\tsigma_gamma_b\tnote")
    for nucid, dsid in CASES:
        eg = egaf(sys.argv[2], nucid)
        if dsid:
            sn, how, did, prim = ensdf(zf, nucid, dsid)
            for e, y in sorted(prim.items(), key=lambda p: -p[1]):
                print(f"{nucid}\t{e}\tENSDF\t{y:.5g}\t\t{did}; {how}; Sn={sn}")
        else:
            e_max = max(e for e, _ in eg)
            sn = e_max + e_max ** 2 / (2 * int(nucid[:3]) * AMU_KEV)
        sig_e = sum(e * s for e, s in eg) / sn
        for e, s in sorted(eg, key=lambda p: -p[1]):
            print(f"{nucid}\t{e}\tEGAF\t{s / sig_e:.5g}\t{s:.4g}\tsigma_E={sig_e:.5g} b; Sn={sn:.2f}")


if __name__ == "__main__":
    main()
