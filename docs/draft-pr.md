# Draft PR: A2/A3 exact xSM basis map and NLO 3D EFT matching audit

## Summary

- Derive the exact Liu--Xie `b1=0` to zero-singlet-VEV parameter map, including the Higgs quadratic coefficient and the opposite `mu^2` sign convention.
- Implement LO scalar matching and NLO xSM scalar matching equations B.5--B.10 of JHEP 09 (2026) 260.
- Verify the zero-temperature stationarity conditions and exact potential identity under singlet translations.
- Verify that matching and field translation commute at LO and NLO, both symbolically and on a deterministic 100-point numerical ensemble.
- Record a candidate second-stage temporal-sector covariance issue as an open question, not a confirmed error.

## Validation

- `python -m pytest -q`: 9 passed.
- `python scripts/run_covariance_scan.py`: PASS, 100 points.
- Maximum relative covariance residual: `1.0374193435774441e-14`.
- Thermal constant checked against its zeta-function definition at 60 decimal digits.
- DRalgo xSM source inspected; execution remains pending because no Mathematica runtime is available in the current environment.

## Scope boundary

This PR does not claim a full xSM parameter scan, two-loop 3D effective potential, bounce, nucleation/percolation, gravitational-wave spectrum, or muon-collider re-projection.
