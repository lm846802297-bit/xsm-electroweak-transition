# Contributing and review policy

1. `main` contains only completed, reviewed research artifacts.
2. Active calculations belong on topic branches. The current topic branch is `research/eft-matching-audit`.
3. Every physics change must include tests and a short derivation or source-equation reference.
4. Numerical results must include source commit, input SHA-256, environment manifest, compute parameters, native output, convergence evidence and an independent physical cross-check.
5. CI is an independent regression check. The primary calculation and validation must first run in the active local research environment.
6. Do not merge a Draft PR or promote a provisional numerical result to `main`.
