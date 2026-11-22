#!/usr/bin/env python3
"""Refresh provenance and write a SHA-256 manifest for repository artifacts."""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "metadata" / "SHA256SUMS"
PROVENANCE = ROOT / "metadata" / "provenance.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# Refresh provenance from the exact files present in this checkout.
prov = json.loads(PROVENANCE.read_text())
prov.setdefault("execution", {})["environment_sha256"] = sha256(ROOT / "metadata" / "environment.json")
prov.setdefault("execution", {})["test_result"] = (ROOT / "results" / "validation" / "pytest.txt").read_text().strip()
prov.setdefault("covariance_validation", {})["source_script_sha256"] = sha256(ROOT / "scripts" / "run_covariance_scan.py")
prov["covariance_validation"]["raw_csv_sha256"] = sha256(ROOT / "results" / "validation" / "covariance-ensemble.csv")
prov["covariance_validation"]["source_module_sha256"] = {
    "basis.py": sha256(ROOT / "src" / "xsm_ewpt" / "basis.py"),
    "matching.py": sha256(ROOT / "src" / "xsm_ewpt" / "matching.py"),
}
prov["source_sha256"] = {
    rel: sha256(ROOT / rel)
    for rel in [
        "src/xsm_ewpt/basis.py",
        "src/xsm_ewpt/matching.py",
        "scripts/run_covariance_scan.py",
        "scripts/write_sha256_manifest.py",
        "tests/test_basis.py",
        "tests/test_matching.py",
        "tests/test_matching_symbolic.py",
        "tests/test_temporal_covariance.py",
    ]
}
prov["archive_workflow"] = {
    "path": ".github/workflows/archive-covariance.yml",
    "sha256": sha256(ROOT / ".github" / "workflows" / "archive-covariance.yml"),
    "manifest_script": "scripts/write_sha256_manifest.py",
    "manifest_script_sha256": sha256(ROOT / "scripts" / "write_sha256_manifest.py"),
    "purpose": "independent GitHub Actions rerun; archives raw CSV, summary, logs, pytest output and complete SHA-256 manifest to the research branch",
}
try:
    prov["validated_checkout_commit"] = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
except subprocess.CalledProcessError:
    prov["validated_checkout_commit"] = "unknown"
PROVENANCE.write_text(json.dumps(prov, indent=2) + "\n")

# Include tracked and non-ignored untracked files, excluding the manifest itself.
raw = subprocess.check_output(
    ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"],
    cwd=ROOT,
)
paths = sorted(
    {
        p.decode()
        for p in raw.split(b"\0")
        if p and p.decode() != "metadata/SHA256SUMS"
    }
)
lines = []
for rel in paths:
    path = ROOT / rel
    if path.is_file():
        lines.append(f"{sha256(path)}  {rel}")
MANIFEST.write_text("\n".join(lines) + "\n")
print(f"Wrote {len(lines)} SHA-256 entries to {MANIFEST.relative_to(ROOT)}")

# Keep provenance and manifest synchronized whenever the archive workflow changes.
