# Research Task Brief: mixed physical-Higgs Z self-energy bubbles

1. Question: verify precisely the two physical-h_i mixed Z and neutral-Goldstone 1PI bubbles in NST (A59) at the 350 GeV fixture, and their signed contribution to (A29) holding all W inputs fixed.
2. Sources: NST Phys. Rev. D 103 115035, arXiv:2103.07467v3, §II, (A26)–(A29), (A36)–(A44), (A59), and APS erratum. JHEP09(2026)260 §4.1, via exact pinned existing input.
3. Source: clean detached b7ac8daf45fafb8f95ac720390fee3fcf7660840 at preflight, with main 5d354669a2c46da1d53b4560476f3981493ed7be and Draft #1 b1ba61e6a27224696822db20367da2e364e9e33b untouched. Track source, signed native data and dirty paths.
4. Conventions: zero temperature, Feynman gauge, NST plus-Π transverse inverse propagator and causal +i0, MS-bar at three pinned scales, tree-level W and Z pole inputs, zero singlet VEV basis and journal sin(theta)=0.1 for squared doublet projections.
5. Scope: Z mixed bubbles only, not seagulls, 1PR, W/fermion parts, complete gauge coupling or EFT matching. A29 is an additive Z-only contribution / relative to the tree hypercharge value, not a prediction of a full running coupling.
6. Known: prior W mixed and W seagull, W/Z fermions and common 1PR audited separately. Input hashes bind those fixtures without copying any observable or computing those diagrams.
7. Derive both neutral couplings directly from the real and imaginary components of the doublet kinetic term and verify graph symmetry factors. Compare independent twice-subtracted spectral derivative from a physical h_i+Z two-body cut versus the D-dimensional shifted analytic loop; check absolute signed values versus the analytic PV reduction and preserve its finite D times UV term.
8. Acceptance: at both zero and Z pole, real direct and reduced values agree within 1e−48 GeV² at 65 working digits; both masses satisfy independent spectral derivative within 1e−48 (GeV^-2); physical below-threshold and signed sampled above each physical threshold; early D=4 control and scale slope; byte replay and central/RCRPS pass.
9. Outputs: derivation, replayable code, pinned inputs, signed native output, tests, checksum and task/run provenance. Stage only; no GitHub write or CI.
10. Remaining: scalar seagull and other bosonic Z sectors, full pole-to-running assembly, DRalgo, thermal EFT, transition and GW outside scope.

Task 20261123-z-mixed-physical-scalar. Central db503f0eb579bdbcac702e44237a6610c34690e7 and lock SHA-256 6585db254d00d803c17307b03cffb2986f51897c6ac98f8a4f8e9fbbf6e412f2; profiles task-handoff, research-compute, research-quality-gate, calculation-result-archive, rcrps-v1.0.1; project extension CONTRIBUTING.md.
