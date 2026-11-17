# dartRverse 1.3.0: state of the packages (18 Nov 2026)

Read this first after a break. Details of the FBM work are in
`fbm-step2-plan.md` (step 2, WP0-WP6) and `fbm-step3-status.md` (step 3:
decisions, downstream audit, readers).

## Versions and branches

All `dev_fbm` branches were merged into `dev` on 17 Nov 2026. Everything is
committed and pushed; all local copies are on `dev`, clean and in sync.

| Package | Version on `dev` | Merge / version commit | Notes |
|---|---|---|---|
| dartR.base | 1.3.0 | e219dff / 64f8a74 | includes Luis's PR #464 (PLINK via `system2()`) |
| dartR.popgen | 1.3.0 | 8fd5e3f / 87020aa | includes Luis's PR #136 (`gl.LDNe` long names) |
| dartR.sim | 1.3.0 | 2895d11 / aa49af3 | |
| dartR.spatial | 1.3.0 | 5d3af66 / b8399ba | |
| dartR.captive | 1.3.0 | 9b168d6 / 984bd04 | |
| dartR.sexlinked | 1.3.0 | ba8279b / 5b038db | |
| dartR.data | **1.2.7** | (no `dev_fbm`) | stays as on CRAN; a 1.3.0 bump was reverted (c77c859) |
| dartRverse | 1.3.0 | 3badc31, 338181a | NEWS entry |

- Downstream packages require `dartR.base (>= 1.3.0)`; none needs more than
  `dartR.data (>= 1.2.7)`.
- `Date` and `Revision` in the DESCRIPTIONs were **not** changed (still
  2026-11-12, "Generous Gecko"): set them for the release.
- The `dev_fbm` branches still exist on GitHub; they are fully merged and
  can be deleted.
- `beta` and `main` were not touched.

## Continuous integration (GitHub, 1.3.0 commits)

All runs finished green: R-CMD-check on macOS, Windows, Ubuntu release and
Ubuntu devel (plus oldrel-1 where configured) for all eight packages, and
the new **FBM-examples** workflow in dartR.base (every help-page example run
in memory and file-backed; fails only when an example fails file-backed but
not in memory). FBM-examples runs on pushes to `dev` and `dev_fbm` and can be
started by hand (workflow_dispatch).

## Local R CMD check (devtools::check, --no-manual, R 4.3.1, Windows)

| Package | Result | Comment |
|---|---|---|
| dartR.base | 1 error, 1 warning, 1 note | error: `gl.report.excess.het` example needed ggtern (now installed, example runs); re-run in progress |
| dartR.data | 0 / 0 / 2 notes | installed size (data package), time check |
| dartR.popgen | 0 / 1 warning / 2 notes | warning: adegenet built under R 4.3.3; notes: `assignPOP` not installed, time check |
| dartR.sexlinked | 0 / 1 warning / 1 note | same R-version warning, time check |
| dartRverse | 0 / 0 / 1 note | time check |
| dartR.sim, dartR.spatial, dartR.captive | running | started 17 Nov 22:30 (paused overnight) |

All issues found so far come from the local machine (R 4.3.1 while
dependencies were built under 4.3.3; optional packages missing; "unable to
verify current time"), none from the packages. A local check on an R
matching the dependencies (e.g. R 4.6.2, now installed) would be clean of
the version warnings.

## Open pull requests

| Repo | PR | State |
|---|---|---|
| dartR.captive | #154 `standalone-engine` (Luis): `gl.relatedness` / `gl.diagnostics.relatedness` run the standalone program `dartR_coancestry` instead of the package `dartR.coancestry`; new `coancestry.path`; runs can be interrupted (`processx`) | **mergeable** (conflict with 1.3.0 `dev` resolved on the branch, 717a2c4: both DESCRIPTION lines kept). **Blocker:** `coancestry_{mac,linux,windows}.zip` are not yet in `dartRverse/binaries` (checked `main` and `dev`); Luis asks not to release `gl.relatedness` before they are; not tested on Windows R |
| dartR.base | #409 `dev_chingching` (2 months), #227 (7 months) | old |
| dartR.popgen | #26 `dev_carlo_DNADot` -> `dev_jesus` (2 years) | stale |

Note: on `dev` (1.3.0) `gl.relatedness` still uses `dartR.coancestry`, which
users cannot install (no releases); #154 fixes that once the zips exist.

## What 1.3.0 contains (summary; see each package's NEWS)

- File-backed (FBM) objects are read from disk in blocks of loci across
  dartR.base and the downstream packages; results identical to in-memory
  objects (tested for both storage types and one-locus blocks).
- Readers: `gl.read.PLINK` (.bed) and `gl.read.dart` read straight into a
  file-backed matrix; `gl.read.vcf(parser = "cpp")` optional fast route via
  tidypopgen (Suggests; diploid SNPs only).
- Functions that need the whole matrix decode explicitly with a memory
  limit (`options(dartR_fbm_max_decode_gb)`, default 4 Gb).
- Exported helpers `utils.fbm.*` for the downstream packages.
- Bug fixes found on the way, e.g. `gl.read.PLINK` ploidy inference.

## Decisions and open items

- Done: step 2, step 3 items 1, 2 (PLINK, DArT, VCF), 5, 7 (CI), downstream
  groups A, B, C.
- Dropped: #3 (save/load without converting; genlight is the storage
  format), #4 (filter without copying; not worth a class change).
- Low priority: #6 (read counts / likelihoods on disk).
- Skipped: CSV reader, group D (writers for external programs), group E
  (exports).
- Pre-existing quirk, not changed: `gl.read.vcf` drops multi-allelic records
  by INFO `AC` when an `AC` key exists in any record, otherwise by ALT.
- Before the CRAN release: set Date/Revision; decide on #154 (needs the
  coancestry zips); local check on a current R.

## Environment notes (this machine)

- R 4.3.1 (main) and R 4.6.2 (installed 17 Nov).
- ggtern 4.0.0 installed in the R 4.3 library (17 Nov).
- tidypopgen has no CRAN binary for R 4.3; it was built from source into a
  scratch library for testing; R 4.6.2 installs it as a binary.
