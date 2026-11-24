# xSM repository storage audit — 2026-11-24

**Scope:** read-only current-tree audit and adoption of the central three-tier scientific data storage policy.  
**Repository:** `lm846802297-bit/xsm-electroweak-transition`  
**Audited branch:** `main`

## 1. Repository baseline

- GitHub repository metadata size: **360 KiB**.
- Complete recursive main tree: **10 blobs**, approximately **0.017 MiB** summed file bytes.
- Recursive tree response was complete (`truncated=false`).
- No file exceeds 1 MiB.
- No numerical output or data bundle is present in the current tree.

The repository currently contains its project instructions, contribution policy, bibliography, Python project configuration, and central-standards validation workflow. GitHub repository metadata size and current-tree byte sum are different measures; this audit does not claim a historical Git object breakdown.

## 2. Storage-tier decision

- **Ordinary Git:** appropriate for all current contents, including derivations, bibliography, configuration, validation code, and future small/medium native result tables.
- **Git LFS:** not enabled; no current dataset justifies it.
- **External object storage:** not configured; no bulk or long-term data object currently requires it.

For future numerical results, retain native files, immutable run directories, input hashes, environment and hardware records, compute parameters, convergence evidence, and independent physical validation. Consider LFS or external storage only when actual data size/access/retention requirements justify it.

## 3. Integrity and migration boundary

- No existing file was modified, deleted, or moved by the audit.
- No Git history was rewritten.
- No calculation was run or rerun.
- No LFS or external storage was configured.

## 4. Standards adoption

The downstream lock pins central commit `14155487afc075d301fa7f6cf452c7bff137e869`, which includes central PR #42's Git / Git LFS / external object storage policy. The xSM-specific EFT matching, bibliography, derivation, and independent-validation requirements remain local.

## 5. Conclusion

This is a compact theory repository with no current numerical data bundle. Ordinary Git is appropriate. **No migration and no LFS activation are required.** Apply the central three-tier policy when substantive numerical outputs are added.
