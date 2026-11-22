# Research log

## 2026-11-23: A2/A3 basis and matching milestone

### Completed

- Implemented exact singlet translation of the general renormalizable xSM potential, including the Higgs quadratic coefficient `cH`.
- Mapped the Liu--Xie `b1=0` basis to the zero-temperature `<S>=0` basis.
- Implemented LO scalar matching.
- Implemented NLO scalar matching Eqs. (B.5)--(B.10) of JHEP 09 (2026) 260.
- Cross-checked the NLO hard-mode tadpole structure against Schicho et al., JHEP 06 (2021) 130, Eq. (4.42).
- Added exact symbolic potential-translation tests, stationarity tests, LO covariance tests, symbolic NLO covariance tests, a 60-digit thermal-constant check and a deterministic 100-point numerical covariance ensemble.
- Inspected the independent DRalgo xSM example and recorded an exact external source snapshot for a later numerical cross-check.
- Flagged that the original path-classification ratio is not invariant under singlet field translations.

### Validation

- `python -m pytest -q`: 9 tests passed.
- Deterministic ensemble: 100 points, seed `20261123`.
- Maximum relative covariance residual: `1.0374193435774441e-14`.
- Maximum absolute residual: `2.3283064365386963e-10`.
- Acceptance threshold: `2e-11` relative residual.
- Native outputs are archived under `results/validation/`.

### Not yet completed

- DRalgo point-by-point comparison (source inspected; Mathematica runtime unavailable).
- Full temporal-gauge-field integration and its field-shift covariance audit.
- Zero-temperature one-loop physical-parameter conversion.
- Two-loop 3D effective potential, bounce, nucleation/percolation and gravitational-wave calculations.
- Full five-dimensional scan and muon-collider re-projection.
