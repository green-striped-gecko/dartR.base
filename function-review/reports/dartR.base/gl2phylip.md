# Review: gl2phylip (dartR.base)

- Family mode: io
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203
- Datasets: testset.gl (SNP), testset.gs (SilicoDArT)
- Baseline: tests/testthat/test-gl2phylip.R (snapshot captured pre-review; 13 assertions, all passing against current code)

## Verdict

**Standards: Needs work** — verbosity gating is clean (`verbose = 0` is fully silent, verified), but the function bypasses `gl.check.wd`, uses a non-conforming build tag, and its roxygen block predates the ratified order.
**Spec: Needs work** — the phylip file itself is well formed (10-character name field, count header, appended bootstrap matrices), but with `bstrap > 1` the documented return value is wrong, the bootstrap resampling path inherits a known genotype-corruption defect, and a third of the frequency table is silently discarded on the reference data.

## Findings

**F1 [HIGH, confidence: high] — `bstrap > 1` returns the last bootstrap replicate, not the observed matrix (DOC5)**
`R/gl2phylip.r:129,152` — the loop reassigns `d` each replicate and the function returns `d` after the loop. `@return` promises the "Matrix of Euclidean distances between populations".
Failure scenario: `result <- gl2phylip(x, bstrap = 100)` hands the caller the 100th resampled matrix; on testset.gl the returned values differ from the observed-data matrix by up to 1.18 on distances of ~7 (17%). Verified; snapshotted in the baseline test.
Proposed change: hold the observed matrix in its own variable and return that.

**F2 [HIGH, confidence: medium] — bootstrap resampling hits the adegenet SNPbin duplicate-index defect (DAT2; known package-wide bug)**
`R/gl2phylip.r:120-122` — `x[, sample(h, size = nLoc(x), replace = TRUE)]` subsets a genlight with repeated locus indices. The dartR `[` method delegates to `SNPbin`'s `[`, which was confirmed in isolation (2026-08-18 sweep, flagged package-wide, deliberately not yet fixed) to silently turn `NA` into `0` for every duplicate occurrence of a repeated index.
Failure scenario: every bootstrap replicate on data with missing genotypes (testset.gl has substantial missingness) computes frequencies on genotypes where resampled-duplicate loci lost their NAs — bootstrap distances are biased toward the reference allele. The call-site pattern is confirmed; the end-to-end magnitude through the distance matrix was not separately quantified here (hence medium confidence on impact, high on mechanism).
Proposed change: rebuild the resampled object from `as.matrix(x)[, idx]` (matrix subsetting handles duplicates correctly) rather than genlight column subsetting — or compute the frequency matrix once and resample its columns, which also removes the per-replicate conversion cost (STY2).

**F3 [MEDIUM, confidence: high] — frequency cells with any missing genotype are discarded (DAT2)**
`R/gl2phylip.r:87-88,125-126` — `mean(e) / 2` without `na.rm = TRUE` makes a population×locus frequency NA whenever a single individual is unscored. On testset.gl, 35.2% of the frequency table is NA; `stats::dist` then drops those coordinates pairwise (with rescaling), so each pairwise distance rests on a different, undisclosed locus subset.
Failure scenario: two populations sharing few completely scored loci get a distance computed from a small rescaled subset, with no warning at any verbosity.
Proposed change: use `mean(e, na.rm = TRUE) / 2`, and warn (VRB4) when cells remain NA (a population entirely unscored at a locus). **Consequence: numerical output changes for any dataset with missing genotypes.**

**F4 [MEDIUM, confidence: high] — SilicoDArT admitted but frequencies halved (DAT7, DAT1)**
`R/gl2phylip.r:55,87` — the datatype check keeps the permissive default while `mean(e) / 2` assumes diploid dosage; for ploidy-1 presence/absence data all frequencies (and hence distances) are halved. The roxygen explicitly claims support for "SNP or presence/absence (SilicoDArT) data". Verified: testset.gs runs and yields max frequency 0.5.
Failure scenario: SilicoDArT distances are uniformly scaled by 0.5 — harmless for tree topology, wrong for any absolute use or comparison with SNP-derived matrices.
Proposed change: divide by ploidy (2 for SNP, 1 for SilicoDArT), or restrict with `accept = "SNP"` and fix the docs — custodian's call on which contract is intended.

**F5 [MEDIUM, confidence: high] — outpath bypasses the package working-directory convention (FS7)**
`R/gl2phylip.r:39` — `outpath = tempdir()` is hardcoded rather than the family idiom `outpath = NULL` resolved via `gl.check.wd()` (used by the other four converters in this batch).
Failure scenario: a user who set a session directory with `gl.set.wd()` gets every other converter's output there, but `gl2phylip` writes to tempdir and the file evaporates with the session.
Proposed change: adopt `outpath = NULL` + `gl.check.wd(outpath, verbose = 0)`.

**F6 [LOW, confidence: high] — unprotected sink (no rule; robustness principle)**
`R/gl2phylip.r:109-141` — the bootstrap loop (including the F2-affected subsetting) runs inside an open `sink` with no `on.exit(sink(NULL))`. Any error mid-loop leaves the R session's output diverted.
Proposed change: `on.exit(if (sink.number() > 0) sink(NULL), add = TRUE)` after opening the sink.

**F7 [LOW, confidence: medium] — 10-character truncation can collide population names (DOC5)**
`R/gl2phylip.r:94-95` — names are truncated to phylip's 10-character field with no uniqueness check. No collision occurs on testset.gl (31 populations remain unique), but two populations sharing their first 10 characters become identical taxa to phylip.
Proposed change: warn when `substr(..., 1, 10)` produces duplicates.

**F8 [LOW, confidence: high] — build tag and roxygen drift (STY1, DOC1, DOC2, DOC7 (proposed rule))**
`R/gl2phylip.r:51` — `build = "Jody"` where the package convention is a `v.<year>.<n>` tag. Roxygen: title says "(SNP)" while the parameters claim SilicoDArT support (align with the F4 decision); `@family linker` vs `linkers` in the paup converters (one taxonomy term per concept); `@return` after `@export` (DOC1, old order — flag the file); `verbose` default clause not the ratified DOC2 wording; custodian-only author block (DOC7 — default fix: add `Author(s): Arthur Georges.`).

## Proposed changes

1. Return the observed-data distance matrix regardless of `bstrap` (F1). **Consequence: the returned object changes for `bstrap > 1` callers.**
2. Resample via matrix indexing (or resample the frequency matrix) to avoid the SNPbin duplicate-index corruption (F2). **Consequence: bootstrap replicate values change on data with missing genotypes.**
3. Add `na.rm = TRUE` to the frequency computation with a VRB4 warning for all-NA cells (F3). **Consequence: numerical output changes for datasets with missing genotypes.**
4. Resolve the SilicoDArT contract: divide by ploidy or restrict to SNP and amend the docs (F4).
5. Adopt the `outpath = NULL` + `gl.check.wd` idiom (F5).
6. Protect the sink with `on.exit` (F6).
7. Warn on 10-character name collisions (F7).
8. Fix the build tag; align title/family/verbose/author roxygen and re-document (F8, DOC4).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT (n/a — no plots), STY — run
- Spec: `bstrap = 1` and `bstrap = 3` executed on testset.gl; file inspected against phylip distance-matrix format (count header, 10-char names, appended replicate matrices); SilicoDArT run on testset.gs subset — run
- verbose = 0 silence and sink release: run — clean
- F2 end-to-end bias quantification through the distance matrix: SKIPPED — mechanism confirmed in the 2026-08-18 isolation test; per-replicate magnitude not re-measured here
- Downstream execution in phylip (neighbor/fitch): SKIPPED — no phylip binary on this machine; format checked against the published spec
- FBM path (DAT6): SKIPPED — no FBM fixture; converter densifies via `as.matrix` by design

## Approval (Phase B)

All findings at every severity approved by Arthur Georges, 2026-09-05, via
the formal approval boxes -- a blanket class approval covering this batch,
explicitly acknowledging that converted outputs change where they were
wrong. The F4 custodian's-call is resolved by the approved DAT7 policy:
a fatal `accept = "SNP"` gate wherever SilicoDArT is silently admitted,
with the docs amended to match.

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges | 2026-09-05, blanket class approval; return changes for bstrap > 1 callers |
| 2 | Approved | Arthur Georges | 2026-09-05; matrix-level resampling per the report, adegenet not touched |
| 3 | Approved | Arthur Georges | 2026-09-05; numerical output changes acknowledged |
| 4 | Approved | Arthur Georges | 2026-09-05; DAT7 policy: fatal accept = "SNP" gate, docs amended |
| 5 | Approved | Arthur Georges | 2026-09-05, blanket class approval |
| 6 | Approved | Arthur Georges | 2026-09-05, blanket class approval |
| 7 | Approved | Arthur Georges | 2026-09-05, blanket class approval |
| 8 | Approved | Arthur Georges | 2026-09-05, blanket class approval |

## Outcome (Phase C)

All eight changes applied on branch `review-gl2phylip` (base
`upstream/dev` = ddaed27), 2026-09-05. Implementation notes: the genlight
is densified once (`mat <- as.matrix(x)`) and bootstrap replicates
resample `mat`'s columns, avoiding the SNPbin duplicate-index defect and
the per-replicate conversion cost; the sink guard is depth-based
(`sink.number()` recorded before opening) so it cannot disturb a sink
already open in the caller such as `capture.output`.

Verification: baseline characterization test updated -- every diff maps
to an approved finding (F1 observed-matrix return, F3 recomputed
distances d[1,2] 6.7673 -> 7.7752 under na.rm = TRUE, F4 SilicoDArT
rejection) -- 13 assertions pass. End-to-end on testset.gl at
verbose = 3 with bstrap = 3: clean, including the new F3 warning (815
unscored population x locus cells reported). verbose = 0 captures zero
lines and leaves no sink open. Sibling grep across the 8 clones: no
package-code callers of `gl2phylip` -- all clear.
PR: https://github.com/green-striped-gecko/dartR.base/pull/PRNUM

```json
{
  "function": "gl2phylip",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "medium", "rule": "DAT2", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DAT2", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DAT7", "status": "applied", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "FS7", "status": "applied", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "none (robustness)", "status": "applied", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "medium", "rule": "DOC5", "status": "applied", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 8}
  ],
  "coverage_skipped": ["F2 end-to-end bias magnitude: mechanism confirmed separately, not re-measured", "phylip execution: no binary", "DAT6: no FBM fixture"],
  "status": "phase-c-complete",
  "pr": null
}
```
