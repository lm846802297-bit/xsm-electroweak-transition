# xSM Electroweak Transition — Repository Instructions

Shared research governance, compute execution, task handoff, scientific-change tracking, quality gates, result archival, structured records, document review, and Git integration are defined by the central research infrastructure repository and pinned in [research/standards-lock.yaml](research/standards-lock.yaml). Start from the [Central Research Standards Index](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/governance/CENTRAL_STANDARDS_INDEX.md).

This repository remains authoritative for xSM physics, scalar-field basis conventions, dimensional-reduction matching, thermal effective-field-theory calculations, perturbative order counting, and project-specific scientific acceptance criteria.

## 1. Preflight

Before substantive work:

1. Read this file, README.md, CONTRIBUTING.md, the PR template, and research/standards-lock.yaml.
2. Read the pinned Central Research Standards Index.
3. For a new derivation, matching calculation, numerical run, parameter scan, or multi-agent handoff, prepare a concise Research Task Brief using the [central Task Handoff Standard](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/RESEARCH_TASK_HANDOFF_STANDARD.md) and [template](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/templates/research-task-brief.md). Use the existing project research layout; do not restructure solely to host a brief.
4. Select the appropriate [central Research Quality Gate](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/governance/RESEARCH_QUALITY_GATE.md) profile. Define checks and acceptance criteria before starting the substantive calculation.
5. Inspect open PRs, exact source/input revisions, existing calculations, raw outputs, and provenance before repeating work.

## 2. xSM and EFT-specific requirements

The following project-specific scientific requirements remain local and authoritative:

- State the xSM field content, scalar potential, vacuum convention, gauge/Yukawa conventions, and parameter basis used in each derivation.
- For every A2/A3 or other basis translation, provide explicit forward and inverse maps, parameter definitions, dimensions, and a symbolic round-trip check.
- Distinguish four-dimensional zero-temperature parameters from three-dimensional EFT parameters after dimensional reduction. State the normalization and temperature factors in the 3D fields and couplings.
- Keep LO and NLO matching contributions separate. Identify loop order, counterterm scheme, renormalization scale, thermal resummation prescription, and the treatment of zero modes.
- Track scalar tadpoles, thermal masses, wave-function/temporal-spatial kinetic terms, and higher-dimensional operators consistently with the selected EFT power counting.
- Do not assume that temporal and spatial kinetic coefficients coincide at finite temperature. Verify the relevant tensor structures and field redefinitions explicitly.
- Keep matching coefficients, phase-transition thermodynamics, bounce/nucleation/percolation, gravitational-wave spectra, and collider observables as distinct stages with explicit dependencies.
- Compare perturbative EFT results with nonperturbative xSM benchmarks only under matched conventions and parameter points.
- A basis identity or symbolic test does not by itself validate a thermal matching coefficient or a physical phase-transition prediction.

## 3. Computation, provenance, and validation

Follow the pinned central [Research Compute Policy](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/RESEARCH_COMPUTE_POLICY.md). For every substantive calculation, preserve:

1. Exact source commit SHA and clean/dirty working-tree state.
2. SHA-256 hashes of source, input data, model parameters, and configuration files.
3. Actual OS/runtime, CPU/GPU, memory, compiler/interpreter, and dependency versions.
4. Exact commands, physical parameters, renormalization scales, perturbative orders, and solver settings.
5. Native symbolic/numerical outputs and the scripts that generate figures.
6. Numerical precision, scale dependence, convergence and uncertainty analysis.
7. Independent analytical or physical cross-checks.
8. A run manifest linking source, inputs, environment, outputs, and validation evidence.

New substantial numerical runs use central RCRPS v1.0.1, including the complete 14-field execution.source_control object. Historical records remain unchanged and are not retrofitted with missing provenance.

## 4. Research lifecycle and Git

- Keep active calculations and matching audits on topic branches. Do not develop directly on main.
- main contains completed and reviewed research records; a Draft PR is not an accepted result.
- Preserve the active A2/A3 and NLO matching work on its existing branch and PR. Governance migration must not retarget, rebase, or rewrite that scientific branch.
- Follow the pinned central [Git merge policy](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/governance/GIT_MERGE_POLICY.md). Substantive research changes default to Create a merge commit.
- After merge, verify the integrated main state and archive. Report software tests, symbolic identities, numerical convergence, and physical validation separately.

## 5. Reporting

Follow the central [Research Task Handoff Standard](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/RESEARCH_TASK_HANDOFF_STANDARD.md) and [Scientific Diff reporting standard](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/scientific-diff-workflow.md#11-research-reporting-and-communication-standard).

Report the conclusion, evidence, quantitative results, and remaining work. Keep completed matching results separate from pending DRalgo comparison, temporal-sector covariance audit, bounce/nucleation/percolation, gravitational-wave spectra, and collider re-projection.
