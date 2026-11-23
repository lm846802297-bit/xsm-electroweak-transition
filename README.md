# xSM Electroweak Transition Research

A reproducible research project studying the real-singlet extension of the Standard Model (xSM), its finite-temperature electroweak phase transition, dimensional-reduction EFT matching, gravitational-wave predictions and multi-TeV muon-collider probes.

## Research scope

- Audit the thermal tadpole and field-basis conventions in Liu and Xie (JHEP 04 (2021) 015, arXiv:2101.10469).
- Reproduce and validate high-temperature dimensional-reduction matching.
- Compare perturbative EFT predictions with nonperturbative xSM results.
- Reassess nucleation, percolation, gravitational-wave and collider observables.

## Shared research standards

This repository adopts the central research standards pinned in [research/standards-lock.yaml](research/standards-lock.yaml). See the [Central Research Standards Index](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/governance/CENTRAL_STANDARDS_INDEX.md). The xSM potential, A2/A3 basis maps, dimensional-reduction matching, thermal EFT conventions and physics-specific validation remain local.

## Repository policy

`main` is reserved for completed, reviewed research results. Active derivations and calculations live on topic branches and enter `main` only through reviewed pull requests. Every numerical result must preserve exact source and input hashes, runtime environment, hardware and compute parameters, native outputs, convergence evidence and an independent physical check.

The active development branch is `research/eft-matching-audit`.

## Primary references

See [`metadata/references.bib`](metadata/references.bib).
