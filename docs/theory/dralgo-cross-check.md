# Independent DRalgo cross-check plan

- Repository: https://github.com/DR-algo/DRalgo
- Snapshot commit: `e00c28e3d01c83f2f2f205ba508b21709eef2d41`
- xSM example: `examples/15_smRealSinglet.m`
- Example blob SHA: `49a6150ac7d5a13807668e27b09ab9267b788cad`
- Soft-to-ultrasoft module blob SHA: `3776d90243ad4d95102dbf1781e3a8982c6d86ce`
- License: GPL-3.0; this project links to DRalgo and does not copy its source.

The DRalgo conventions are `mu_m -> a1`, `lambda_m -> a2`, `c3 -> b3`, `lambda_sigma -> b4`, `mu1 -> b1`, `mu_sigma -> b2`, and `lambda_1H -> lambda`.

Validation procedure:
1. Use identical renormalized 4D input, scale, temperature, gauge couplings and top Yukawa.
2. Run `PerformDRhard[]` and extract NLO tadpole, scalar masses and couplings.
3. Compare against this implementation after matching conventions.
4. Run `PerformDRsoft[...]` and compare second-stage EFT separately.
5. Preserve Mathematica/GroupMath versions, scripts, native outputs and hashes.

Execution status: source inspected; current Linux environment has no Mathematica, WolframKernel or wolframscript. No DRalgo numerical result is claimed.
