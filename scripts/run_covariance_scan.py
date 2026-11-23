#!/usr/bin/env python3
"""Generate a deterministic raw-output audit of NLO field-shift covariance."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import random
from pathlib import Path

from xsm_ewpt.basis import XSMParameters, shift_singlet_basis
from xsm_ewpt.matching import MatchingScales, match_scalar_sector_nlo, shift_3d_basis

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "metadata" / "covariance-ensemble.json"
OUT = ROOT / "results" / "validation" / "covariance-ensemble.csv"
SUMMARY = ROOT / "results" / "validation" / "covariance-summary.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    cfg = json.loads(CONFIG.read_text())
    rng = random.Random(cfg["seed"])
    rows = []
    max_abs = 0.0
    max_rel = 0.0
    worst = None
    for i in range(cfg["n_points"]):
        p = XSMParameters(
            cH=rng.uniform(*cfg["ranges"]["cH"]),
            lam=rng.uniform(*cfg["ranges"]["lam"]),
            a1=rng.uniform(*cfg["ranges"]["a1"]),
            a2=rng.uniform(*cfg["ranges"]["a2"]),
            b1=rng.uniform(*cfg["ranges"]["b1"]),
            b2=rng.uniform(*cfg["ranges"]["b2"]),
            b3=rng.uniform(*cfg["ranges"]["b3"]),
            b4=rng.uniform(*cfg["ranges"]["b4"]),
        )
        T = rng.uniform(*cfg["ranges"]["temperature"])
        parent = rng.uniform(*cfg["ranges"]["parent_scale"])
        eft = rng.uniform(*cfg["ranges"]["eft_scale"])
        delta = rng.uniform(*cfg["ranges"]["delta"])
        s = MatchingScales(T, parent, eft, cfg["couplings"]["g2"], cfg["couplings"]["gp2"], cfg["couplings"]["yt2"])

        path_a = shift_3d_basis(match_scalar_sector_nlo(p, s), delta / math.sqrt(T))
        path_b = match_scalar_sector_nlo(shift_singlet_basis(p, delta), s)
        row = {
            "point": i,
            "cH": p.cH,
            "lambda": p.lam,
            "a1": p.a1,
            "a2": p.a2,
            "b1": p.b1,
            "b2": p.b2,
            "b3": p.b3,
            "b4": p.b4,
            "temperature": T,
            "parent_scale": parent,
            "eft_scale": eft,
            "delta": delta,
            "g2": cfg["couplings"]["g2"],
            "gp2": cfg["couplings"]["gp2"],
            "yt2": cfg["couplings"]["yt2"],
        }
        for field in path_a.__dataclass_fields__:
            x, y = getattr(path_a, field), getattr(path_b, field)
            abs_err = abs(x - y)
            rel_err = abs_err / max(1.0, abs(x), abs(y))
            row[f"{field}_path_a"] = x
            row[f"{field}_path_b"] = y
            row[f"{field}_abs_error"] = abs_err
            row[f"{field}_rel_error"] = rel_err
            if abs_err > max_abs:
                max_abs = abs_err
                worst = {"point": i, "field": field, "path_a": x, "path_b": y}
            max_rel = max(max_rel, rel_err)
        rows.append(row)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "n_points": len(rows),
        "seed": cfg["seed"],
        "max_absolute_residual": max_abs,
        "max_relative_residual": max_rel,
        "worst_component": worst,
        "acceptance_relative_tolerance": cfg["acceptance_relative_tolerance"],
        "pass": max_rel <= cfg["acceptance_relative_tolerance"],
        "input_config_sha256": sha256(CONFIG),
        "source_script_sha256": sha256(Path(__file__)),
        "source_module_sha256": {
            "basis.py": sha256(ROOT / "src" / "xsm_ewpt" / "basis.py"),
            "matching.py": sha256(ROOT / "src" / "xsm_ewpt" / "matching.py"),
        },
        "raw_csv_sha256": sha256(OUT),
    }
    SUMMARY.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if not summary["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

# The archive-covariance GitHub Actions workflow watches this source path.
