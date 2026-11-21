# Reproduction: G4NUDEXLIB-1.0 lead primary-gamma intensities

Requirements:
- Geant4 11.3.2 or 11.4.x, built with its HP data (G4NDL) and sourced via `geant4.sh`;
- the G4NUDEXLIB-1.0 data set;
- CMake 3.16 or later, a C++17 compiler, and Python 3.

Tested with Geant4 11.3.2 and 11.4.2, gcc 13.3 and CMake 3.28 on Linux.

The reference table `reference/lead_reference.tsv` is already in the repository. Steps 3 and 7 regenerate it and run the survey of other compounds from the public ENSDF archive and the EGAF list; both steps are optional.

```sh
# 1. Clone and build (point Geant4_DIR at your install if CMake cannot find it)
git clone https://github.com/BlaineHeffron/g4nudexlib-lead-primaries-repro.git
cd g4nudexlib-lead-primaries-repro
source /path/to/geant4-install/bin/geant4.sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# 2. Get the NuDEX data set, the one Geant4 11.3 and 11.4 declare
mkdir -p data
curl -L -o data/G4NUDEXLIB.1.0.tar.gz https://cern.ch/geant4-data/datasets/G4NUDEXLIB.1.0.tar.gz
sha256sum data/G4NUDEXLIB.1.0.tar.gz   # cac7d65e9c5af8edba2b2667d5822e16aaf99065c95f805e76de4cc86395f415
tar -xzf data/G4NUDEXLIB.1.0.tar.gz -C data
# Geant4 11.3.x up to 11.4.1 builds NuDEX file names in a 100-character buffer
# (G4NuDEXLevelDensity.cc:68), so keep the data path short, e.g. through a symlink:
ln -sfn "$PWD/data" /tmp/g4nx
export G4NUDEXLIBDATA=/tmp/g4nx/G4NUDEXLIB1.0

# 3. (optional) Regenerate the reference table from ENSDF and EGAF
#    ensdf.zip: NNDC ENSDF archive (https://www.nndc.bnl.gov/ensdfarchivals/), tested with the 2026-03-31
#    archive, sha256 572c0b94ce2de8e5d35a6384425aeff1d18ecf1d770f6627288c52a80b8f37ce
#    nglist_e.dat: IAEA EGAF isotopic list (https://www-nds.iaea.org/pgaa/), tested with sha256
#    b3147d245ad28333134264d97438499e7529a1f30d611176aed54dc087a4fc4d
python3 scripts/extract_reference.py /path/to/ensdf.zip /path/to/nglist_e.dat > reference/lead_reference.tsv
git diff --stat reference/   # no change with the files above

# 4. Run NuDEX cascades with the shipped library: pure Pb-204, Pb-206 and Pb-207, 0.0253 eV neutrons
#    (each about 1-2 min; Pb-206 is slower because neutrons diffuse longer before capture)
mkdir -p runs
./build/nudex_lines 204 100000 204 0.0005 > runs/shipped-Pb204.tsv
./build/nudex_lines 207 100000 207 0.0005 > runs/shipped-Pb207.tsv
./build/nudex_lines 206  20000 206 0.0005 > runs/shipped-Pb206.tsv
grep -h "Primary thermal gammas" runs/shipped-Pb*.tsv

# 5. Print library, ENSDF, EGAF and emitted yields side by side
python3 scripts/compare.py $G4NUDEXLIBDATA/PrimaryCaptureGammas.dat reference/lead_reference.tsv runs/shipped-Pb*.tsv

# 6. Same runs with the three blocks rewritten per capture (the suggested fix)
python3 scripts/make_fixed_library.py $G4NUDEXLIBDATA reference/lead_reference.tsv data/fixed/G4NUDEXLIB1.0
ln -sfn "$PWD/data/fixed" /tmp/g4nxf
G4NUDEXLIBDATA=/tmp/g4nxf/G4NUDEXLIB1.0 ./build/nudex_lines 204 100000 204 0.0005 > runs/fixed-Pb204.tsv
G4NUDEXLIBDATA=/tmp/g4nxf/G4NUDEXLIB1.0 ./build/nudex_lines 207 100000 207 0.0005 > runs/fixed-Pb207.tsv
G4NUDEXLIBDATA=/tmp/g4nxf/G4NUDEXLIB1.0 ./build/nudex_lines 206  20000 206 0.0005 > runs/fixed-Pb206.tsv
python3 scripts/compare.py /tmp/g4nxf/G4NUDEXLIB1.0/PrimaryCaptureGammas.dat reference/lead_reference.tsv runs/fixed-Pb*.tsv

# 7. (optional) Other compounds: NuDEX vs ENSDF and EGAF, absolute and relative
python3 scripts/check_other_compounds.py $G4NUDEXLIBDATA/PrimaryCaptureGammas.dat /path/to/ensdf.zip /path/to/nglist_e.dat
```

## Expected output

Step 4 prints NuDEX's own count of the tabulated primaries it accepted:

```
 NuDEX: Primary thermal gammas for ZA=82205 found in the database: 63.5557 %
 NuDEX: Primary thermal gammas for ZA=82207 found in the database: 94.3288 %
 NuDEX: Primary thermal gammas for ZA=82208 found in the database: 87.3573 %
```

Step 5 output with Geant4 11.3.2 (`expected/compare-g4-11.3.2-shipped.txt`; the 11.4.2 run is in
`expected/compare-g4-11.4.2-shipped.txt`). All values are per capture:

```
204Pb + n -> 205Pb   library block: 'ZA   82205  205PB  204PB(N,G) E=THERMAL: {~EGAF}'
  header: Sn=6731.8 keV, 15 lines, sum=63.5558 %, 4th=10757.6, 5th=757.576  (sum of primaries 0.6356 per capture)
      E_keV   library     ENSDF      EGAF   emitted      +-
    6729.38    0.3442    0.4850    0.4790    0.3445  0.0019
    4613.46    0.0732    0.1100    0.1087    0.0717  0.0008
    4919.53    0.0656    0.0940    0.0928    0.0659  0.0008

206Pb + n -> 207Pb   library block: 'ZA   82207  207PB  206PB(N,G) E=THERMAL: {~EGAF}'
  header: Sn=6737.74 keV, 2 lines, sum=94.3287 %, 4th=13558.8, 5th=55555.6  (sum of primaries 0.9433 per capture)
      E_keV   library     ENSDF      EGAF   emitted      +-
    6737.62    0.9369         -    0.9882    0.9352  0.0068

207Pb + n -> 208Pb   library block: 'ZA   82208  208PB  207PB(N,G) E=THERMAL: {~EGAF}'
  header: Sn=7367.92 keV, 8 lines, sum=87.3573 %, 4th=635.534, 5th=2808.99  (sum of primaries 0.8736 per capture)
      E_keV   library     ENSDF      EGAF   emitted      +-
    7367.78    0.8707    0.9969    0.9969    0.8692  0.0029
```

The emitted yields follow the library and not the references. In `runs/shipped-Pb207.tsv`, the bins at 6512 keV (0.028 per capture) and 3917 keV (0.014) are statistical-model primaries that fill the missing 12.6 %. EGAF has no 208Pb line there.

Step 6 output (`expected/compare-g4-11.3.2-fixed.txt`):

```
    6729.36    0.4850    0.4850    0.4790    0.4697  0.0022    (ENSDF primaries sum to 1.034; NuDEX normalizes)
    6737.62    0.9882         -    0.9882    0.9929  0.0070
    7367.78    0.9969    0.9969    0.9969    0.9972  0.0032
```

Step 7 output: `expected/check-other-compounds.txt`.

Numbers vary within the quoted statistical errors with seed and Geant4 version. A capture count slightly above the event count, such as 100001 captures from 100000 events, comes from photoneutrons: 7.37 MeV gammas exceed S_n(207Pb).
