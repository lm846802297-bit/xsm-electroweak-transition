# Research log

## 2026-11-23: A2/A3 basis and matching milestone

### Completed

- Implemented exact singlet translation of the general renormalizable xSM potential, including the Higgs quadratic coefficient cH.
- Mapped the Liu--Xie b1=0 basis to the zero-temperature <S>=0 basis.
- Implemented LO scalar matching.
- Implemented NLO scalar matching Eqs. (B.5)--(B.10) of JHEP 09 (2026) 260.
- Cross-checked the NLO hard-mode tadpole structure against Schicho et al., JHEP 06 (2021) 130, Eq. (4.42).
- Added exact symbolic potential-translation tests, stationarity tests, LO covariance tests, symbolic NLO covariance tests, a 60-digit thermal-constant check and a deterministic 100-point numerical covariance ensemble.
- Inspected the independent DRalgo xSM example and recorded an exact external source snapshot for a later numerical cross-check.
- Flagged that the original path-classification ratio is not invariant under singlet field translations.

### Validation

- python -m pytest -q: 11 tests passed.
- Deterministic ensemble: 100 points, seed 20261123.
- Maximum relative covariance residual: 1.0374193435774441e-14.
- Maximum absolute covariance residual in the archived ensemble: 2.3283064365386963e-10.
- Acceptance threshold: 2e-11 relative residual.
- Native outputs: results/validation/covariance-ensemble.csv, covariance-summary.json, covariance-run.log, pytest.txt.

## 2026-11-23: A4 preliminary temporal-sector audit

- Derived the finite-shift transformation of the temporal operators in Eq. (B.1).
- Compared the required x3 -> x3 + 2 y3 delta3 with Eqs. (B.18)--(B.21) of arXiv:2409.17554. The published x3 relation gives a shift twice as large.
- Derived the relative portal normalization from the hard-mode Higgs determinant.
- Inspected the DRalgo xSM temporal-coupling test and tensor conventions at commit e00c28e3d01c83f2f2f205ba508b21709eef2d41; source-level values support a factor-one-half candidate for the x3 action coefficient.
- Added symbolic diagnostic tests. These tests document the mismatch and candidate normalization; they do not claim that the published formula is definitively wrong.
- Mathematica/DRalgo runtime execution and a complete independent gauge-covariant 1PI derivation remain open.

### Not yet completed

- DRalgo point-by-point comparison (source inspected; Mathematica runtime unavailable).
- Full temporal-gauge-field integration and its field-shift covariance audit.
- Zero-temperature one-loop physical-parameter conversion.
- Two-loop 3D effective potential, bounce, nucleation/percolation and gravitational-wave calculations.
- Full five-dimensional scan and muon-collider re-projection.
