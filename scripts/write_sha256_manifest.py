#!/usr/bin/env python3
"""Write a SHA-256 manifest for tracked and non-ignored repository files."""
from __future__ import annotations
import hashlib
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "metadata" / "SHA256SUMS"
raw = subprocess.check_output(
    ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"],
    cwd=ROOT,
)
paths = sorted({p.decode() for p in raw.split(b"\0") if p and p.decode() != "metadata/SHA256SUMS"})
lines = []
for rel in paths:
    path = ROOT / rel
    if path.is_file():
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        lines.append(f"{digest}  {rel}")
MANIFEST.write_text("\n".join(lines) + "\n")
print(f"Wrote {len(lines)} SHA-256 entries to {MANIFEST.relative_to(ROOT)}")
