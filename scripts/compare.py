"""Print the lead primary-gamma intensities side by side: shipped G4NUDEXLIB, ENSDF, EGAF, NuDEX emitted.

Usage: python3 scripts/compare.py PRIMARY_CAPTURE_GAMMAS_DAT reference/lead_reference.tsv [RUN_TSV ...]

PRIMARY_CAPTURE_GAMMAS_DAT  G4NUDEXLIB1.0/PrimaryCaptureGammas.dat (the file NuDEX reads)
RUN_TSV                     output of nudex_lines (one per target); the emitted yield of a line is the sum of
                            the 1 keV bins within 3 keV of E (NuDEX emits a primary at the energy of the level it
                            matched within 10 keV, so an emitted line can sit a keV or two from the tabulated E;
                            bins below the run's MIN_YIELD are not written and show as '-')
"""
import re
import sys

PRODUCTS = {82205: ("205PB", "204Pb + n", 204), 82207: ("207PB", "206Pb + n", 206), 82208: ("208PB", "207Pb + n", 207)}


def library(path):
    out, lines, i = {}, open(path).read().splitlines(), 0
    while i < len(lines):
        f = lines[i].split()
        if f and f[0] == "ZA" and int(f[1]) in PRODUCTS:
            head = lines[i + 1].split()
            n = int(head[1])
            prim = [tuple(map(float, lines[i + 2 + k].split()[:2])) for k in range(n)]
            out[int(f[1])] = (lines[i].strip(), lines[i + 1].strip(), prim)
            i += 2 + n
        else:
            i += 1
    return out


def reference(path):
    ref = {}
    for l in open(path).read().splitlines()[1:]:
        prod, e, src, y = l.split("\t")[:4]
        ref.setdefault((prod, src), []).append((float(e), float(y)))
    return ref


def run(path):
    meta, bins = None, {}
    for l in open(path):
        if l.startswith("# target"):
            meta = l.strip()
        elif re.match(r"^\d+\t", l):
            b, c, y = l.split("\t")[:3]
            bins[int(b)] = float(y)
    a = int(re.search(r"Pb-(\d+)", meta).group(1))
    caps = int(re.search(r"captures (\d+)", meta).group(1))
    return a, caps, meta, bins


def near(rows, e, tol=1.0):
    return sum(y for ee, y in rows if abs(ee - e) <= tol)


def main():
    lib, ref = library(sys.argv[1]), reference(sys.argv[2])
    runs = {}
    for p in sys.argv[3:]:
        a, caps, meta, bins = run(p)
        runs[a] = (caps, meta, bins)
    for za, (prod, label, a) in PRODUCTS.items():
        head, nums, prim = lib[za]
        h = nums.split()
        print(f"\n{label} -> {prod[:3]}Pb   library block: '{head}'")
        print(f"  header: Sn={h[0]} keV, {h[1]} lines, sum={h[2]} %, 4th={h[3]}, 5th={h[4]}"
              f"  (sum of primaries {sum(y for _, y in prim) / 100:.4f} per capture)")
        r = runs.get(a)
        if r:
            print(f"  run: {r[1][2:]}")
        print(f"  {'E_keV':>9} {'library':>9} {'ENSDF':>9} {'EGAF':>9} {'emitted':>9} {'+-':>7}")
        for e, y in sorted(prim, key=lambda p: -p[1])[:6]:
            en, eg = near(ref.get((prod, "ENSDF"), []), e), near(ref.get((prod, "EGAF"), []), e)
            em = err = None
            if r:
                caps, _, bins = r
                em = sum(bins.get(b, 0.0) for b in range(int(e) - 3, int(e) + 4))
                err = (em / caps) ** 0.5 if em else 0.0
            f = lambda v: f"{v:9.4f}" if v else f"{'-':>9}"
            print(f"  {e:9.2f} {y / 100:9.4f} {f(en)} {f(eg)} {f(em)} {'' if err is None else f'{err:7.4f}'}")


if __name__ == "__main__":
    main()
