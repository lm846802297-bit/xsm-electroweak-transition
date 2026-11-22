# Validation artifacts

covariance-ensemble.csv is the native per-point output for the deterministic 100-point NLO field-shift covariance ensemble.

covariance-summary.json records the seed, acceptance threshold, maximum absolute and relative residuals, worst component, and SHA-256 hashes of the input configuration, source script, source modules and CSV output.

covariance-run.log is the stdout record of the ensemble run.

pytest.txt records the full automated test result. The suite includes exact symbolic potential translation, zero-temperature stationarity, LO matching, numerical and symbolic NLO covariance, a 60-digit check of the thermal constant, and symbolic diagnostics for the temporal-sector factor-of-two discrepancy.

The temporal-sector tests document a candidate matching discrepancy; they do not certify a corrected published coefficient. No Mathematica/DRalgo runtime execution, bounce, nucleation, percolation or gravitational-wave calculation is represented by these files.
