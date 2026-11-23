# xSM Electroweak Transition: EFT Matching Audit

Public research repository for an independent audit of the real-singlet extension of the Standard Model (xSM), focusing on:

1. exact field-shift maps between the Liu--Xie b1 = 0 basis and the <S> = 0 basis used by Ramsey-Musolf, Tenkanen and Tran;
2. leading- and next-to-leading-order (NLO) high-temperature dimensional-reduction matching;
3. covariance of the 3D EFT matching under the singlet field redefinition;
4. an audit of field-shift covariance in the second-stage temporal-gauge-field matching;
5. subsequent phase-transition, nucleation and gravitational-wave calculations.

## Current status

- Analytic tree-level basis map: implemented.
- LO 3D scalar matching: implemented.
- NLO scalar matching equations B.5--B.10 of JHEP 09 (2026) 260: implemented.
- Validation: 11 tests passed, including exact symbolic NLO covariance, a deterministic 100-point parameter ensemble, a 60-digit thermal constant check and symbolic diagnostics for the temporal-sector factor-of-two discrepancy.
- Maximum relative NLO covariance residual: 1.04e-14 (acceptance threshold 2e-11).
- DRalgo xSM example and temporal-coupling test inspected; execution awaits a Mathematica runtime.
- Second-stage temporal-field matching: an algebraic factor-of-two discrepancy in the displayed x3 matching relation has been isolated. A complete 1PI derivation and executed DRalgo check remain pending.
- Full xSM scan, bounce, percolation and GW pipeline: not started.

## Conventions

The code stores the coefficient cH multiplying H^dagger H directly:

V contains cH * (H^dagger H) + lambda * (H^dagger H)^2.

Thus:

- Liu--Xie convention: cH = -mu2_LX.
- Ramsey-Musolf--Tenkanen--Tran convention: cH = mu2_RMT.

This avoids silently mixing the opposite signs used for the Higgs quadratic term in the two papers.

The basis map uses S_LX = S_RMT + delta, with delta = v_s when mapping the Liu--Xie basis to the zero-temperature <S> = 0 basis.

## Quick start

    python -m pip install -e '.[test]'
    pytest

## Research policy

Every numerical result must record source commit, input SHA-256, software versions, hardware, numerical tolerances, native output, convergence evidence and an independent physics check. Unfinished work stays on a topic branch and enters main only after PR review.

## Primary references

- W. Liu and K.-P. Xie, Probing electroweak phase transition with multi-TeV muon colliders and gravitational waves, JHEP 04 (2021) 015, arXiv:2101.10469.
- P. M. Schicho, T. V. I. Tenkanen and J. Österman, Robust approach to thermal resummation: Standard Model meets a singlet, JHEP 06 (2021) 130, arXiv:2102.11145.
- L. Niemi, M. J. Ramsey-Musolf and G. Xia, Nonperturbative study of the electroweak phase transition in the real scalar singlet extended standard model, Phys. Rev. D 110, 115016 (2024), arXiv:2405.01191.
- M. J. Ramsey-Musolf, T. V. I. Tenkanen and V. Q. Tran, Refining gravitational wave and collider physics dialogue via singlet scalar extension, JHEP 09 (2026) 260, arXiv:2409.17554.

See docs/theory/ for derivations and conventions.
