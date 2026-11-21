# g4nudexlib-lead-primaries-repro

This repository reproduces a data problem in G4NUDEXLIB-1.0 (the data set behind Geant4's NuDEX capture model, in Geant4 11.3 through 11.5-beta). In `PrimaryCaptureGammas.dat`, the `{~EGAF}` blocks for Pb-204 + n, Pb-206 + n and Pb-207 + n do not give intensities per capture, and NuDEX emits them as tabulated. For example, the 7368 keV primary from Pb-207 + n is emitted at 0.869 per capture, while ENSDF and EGAF both give 0.997.

- `src/nudex_lines.cc`: Geant4 application. It builds a pure-isotope Pb cube, fires 0.0253 eV neutrons with QGSP_BERT_HP and NuDEX enabled, and writes the gamma yield per capture in 1 keV bins.
- `scripts/extract_reference.py`: extracts per-capture ENSDF and EGAF values into `reference/lead_reference.tsv`.
- `scripts/compare.py`: prints the library, ENSDF, EGAF and emitted values side by side.
- `scripts/make_fixed_library.py`: writes a data tree with the three blocks put on a per-capture basis.
- `scripts/check_other_compounds.py`: checks other flagged `{~EGAF}` blocks against ENSDF and EGAF.
- `expected/`: outputs from Geant4 11.3.2 and 11.4.2.

See REPRO.md for the steps.
