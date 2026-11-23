# Research Task Brief: neutral weak top pole block

Task ID: 20261123-top-neutral-weak. Pre-execution definition, 2026-11-23 UTC.

1. **Question.** Evaluate only the Z and neutral Goldstone contributions to the one-loop top pole self-energy in NST (A61) and their isolated insertion into (A30). Exclude all W terms, other top self-energies, physical Higgs states, 1PR tadpoles, full matching or thermal observables.
2. **Primary material.** NST, *Phys. Rev. D* 103 115035, arXiv:2103.07467v3, A20–A23, A26, A30, A36–A39, A42, A61, and APS erratum 109 039902. The pinned, already audited inputs are in the previous charged weak audit at the baseline below.
3. **Baseline.** xSM branch research/top-charged-weak-audit exact 8f368bf32a9fc7ebbd43f64adfb9b576baf363e1, clean detached checkout; its main is 5d354669a2c46da1d53b4560476f3981493ed7be. Draft #1 remains b1ba61e6a27224696822db20367da2e364e9e33b; no changes there. Track new source hashes and local dirty paths.
4. **Conventions.** Zero-temperature 4D electroweak theory with real singlet, tree singlet VEV set to zero, NST positive inverse-propagator self-energy, MS-bar, Feynman +i0, D=4−2ε, xi=1. Internal top and Z/Goldstone masses are the real pole inputs; Goldstone mass equals Z mass; couplings are tree-level Fermi-constant values. Masses and scales in GeV; dimensionless self-energies.
5. **Assumptions.** NST's reference input at mu=MZ, two diagnostic scales MW and Mt. This subset has no singlet parameters or basis change. Gauge-fixed intermediate class only.
6. **Known.** The charged weak result at 8f368bf is independent; confirm identical input file SHA. Massless vectors elsewhere are not used as loop data.
7. **Derivation.** Check the exact printed neutral Goldstone parentheses, derive chiral vector and pseudoscalar numerators, retain D-times-pole term, compare direct shifted integral with analytic scalar PV and vector reduction.
8. **Acceptance.** Complex paths agree within 1e-55 (65 digits) at the three scales. Drift at 45 vs 65 digits below 1e-40. Physical pole has zero absorptive part; below, exact, and above synthetic equal-mass threshold distinguish causal branch; evanescent sign and coefficients checked. Native byte replay, actual lock and task fingerprints, central adoption and RCRPS 1.0.1 checks pass. No CI requested.
9. **Outputs.** Scope-limited derivation, code, pinned input, native signed JSON, tests/checksums, run manifest/logs and linked machine-readable task record. Stage only.
10. **Remaining.** Full running Yukawa still requires massive physical Higgs and tadpole terms and remaining W/EW projections, with a consistent gauge-complete convention.

Central source db503f0eb579bdbcac702e44237a6610c34690e7. Canonical lock research/standards-lock.yaml SHA-256 6585db254d00d803c17307b03cffb2986f51897c6ac98f8a4f8e9fbbf6e412f2. Applied profiles: task-handoff, research-compute, research-quality-gate, calculation-result-archive, rcrps-v1.0.1. Applied local extension CONTRIBUTING.md for reproducibility and independent physics validation.
