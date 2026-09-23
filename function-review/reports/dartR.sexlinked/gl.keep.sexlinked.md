# Review: gl.keep.sexlinked (dartR.sexlinked)

- Family mode: modify
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: f293fba (version 1.2.2, loaded with `devtools::load_all()`), branch `review-gl.keep.sexlinked` (stacked on PR #33)
- Datasets: LBP, testset.gl, testset.gs, platypus.gl (dartR.data 1.2.5); platypus_10KSNPs.gl (local, 200 ind x 10,000 loci)
- Baseline: `tests/testthat/test-gl.keep.sexlinked.R` (31 expectations, all pass on f293fba)
- Related: `gl.report.sexlinked` review, PR dartR.sexlinked#33 (open)

## Verdict

**Standards: Needs work**: the function carries the same classification code as `gl.report.sexlinked()`, and with it the same defects. It also returns genlight objects whose locus metadata is out of sync when the input is a plain `genlight`.

**Spec: Needs work**: with a missing `id` column or no males, it returns four `NULL` subsets and no error. `gl.infer.sex()` then crashes on that output with "invalid 'times' argument".

What works well: empty categories are guarded (lines 720-742), so the zero-loci crash described in the local note `BUG-gl.keep.sexlinked-zero-loci.md` no longer occurs. Returned dartR objects keep ploidy and in-sync `loc.metrics`, and the input object is not modified.

## Findings

The results table is identical to `gl.report.sexlinked()` on `origin/dev` for LBP in both `xy` and `zw` (`all.equal` TRUE), so findings F1-F6, F10 and F11 of the report review apply here unchanged. They are restated briefly, with keep-specific consequences.

**K1 [HIGH, confidence: high]: `loc.metrics` not subset for plain genlight input (DAT2, DAT3)**
`R/gl.keep.sexlinked.r:720-742`: the subsets are made with `x[, a]` and rely on the dartR `[` method to subset `loc.metrics`.
Failure scenario: `as(LBP, "genlight")` returns `$x.linked` with 66 loci but 1,000 `loc.metrics` rows. Any later function that reads locus metrics by row mismatches loci and metadata.
Proposed change: after each subset, re-subset from the original object: `A@other$loc.metrics <- x@other$loc.metrics[a, , drop = FALSE]` (the DAT3 idiom).

**K2 [HIGH, confidence: high]: individuals matched through an unchecked `id` column (DAT5); report F1**
`R/gl.keep.sexlinked.r:193-211`.
Failure scenario: `testset.gl` with no `id` column matches no individuals and returns four `NULL` subsets, with no error.
Proposed change: as report change 1: read sex from the `ind.metrics` rows in order, and error if the row count differs from `nInd(x)`.

**K3 [HIGH, confidence: high]: no guard for a missing sex (FS5); report F2**
`R/gl.keep.sexlinked.r:186-194`.
Failure scenario: `testset.gl` with no males returns `results.table` plus four `NULL`s. Passing that to `gl.infer.sex()` fails with "invalid 'times' argument".
Proposed change: as report change 2: error when either sex has no individuals.

**K4 [HIGH, confidence: high]: crash when a sex has exactly one individual (DAT5); report F3**
`R/gl.keep.sexlinked.r:210-211`.
Failure scenario: `testset.gl` with one female fails with `'x' must be an array of at least two dimensions`.
Proposed change: `drop = FALSE`.

**K5 [MEDIUM, confidence: high]: parallel workers left running after an error (FS5); report F4**
`R/gl.keep.sexlinked.r:137-140, 766-768`.
Failure scenario: `system = NULL, ncores = 2` leaves 2 connections open (0 → 2).
Proposed change: start the cluster after validation and stop it with `on.exit()`.

**K6 [MEDIUM, confidence: high]: SilicoDArT accepted (DAT1); report F5**
`R/gl.keep.sexlinked.r:130-135`.
Failure scenario: `testset.gs` returns subsets of presence/absence loci, classified partly on a "heterozygosity" that has no meaning for these data.
Proposed change: accept SNP data only.

**K7 [MEDIUM, confidence: medium]: zeros replaced by 1 before Fisher's test (numerical correctness); report F6**
`R/gl.keep.sexlinked.r:270, 313, 490, 539`.
Failure scenario: once PR #33 merges, the report uses observed counts and keep does not. On platypus_10KSNPs.gl the report then finds 95 x-linked loci, while keep returns 85.
Proposed change: as report change 6: Fisher's test on the observed counts; the chi-square branch keeps the replacement.

**K8 [LOW, confidence: high]: no history on the returned objects (FS8)**
`R/gl.keep.sexlinked.r:720-742`: the returned genlight objects are modified (loci removed), but no `match.call()` is appended to `@other$history`.
Failure scenario: `r$x.linked@other$history` ends at the entry that created LBP, so the object's provenance does not show that non-sex-linked loci were dropped.
Proposed change: append the call to the history of each non-`NULL` subset.

**K9 [LOW, confidence: high]: documentation disagrees with behaviour (DOC5 (proposed rule), DOC1, DOC2, DOC7 (proposed rule))**
`R/gl.keep.sexlinked.r:1-88`:
- `@return` says "a list of 5 elements and 4 plots", but the plots are not returned.
- Empty categories are returned as `NULL`, which is undocumented.
- The `id` column requirement is stated (it changes with K2).
- The thresholds, the sex codes and the chi-square switch are undocumented.
- There is no `@family` and no `Author(s):` line.

Failure scenario: a user writes `nLoc(r$w.linked)` and gets an error on datasets without W-linked loci.
Proposed change: a docs rewrite matching the report's, plus a note that empty categories are `NULL`.

**K10 [LOW, confidence: high]: `plot.theme` ignored (PLT1); report F9**
Lines 391-668: none of the four plots uses `plot.theme`.
Proposed change: add `plot.theme` to the four plots.

**K11 [LOW, confidence: high]: column 11 referenced by position (STY1); report F10**
Lines 371, 462, 511.
Proposed change: refer to the column by name.

**K12 [LOW, confidence: high]: row-by-row loops (STY2); report F11**
Failure scenario: platypus_10KSNPs.gl takes 6.4 s serial (the report took 6.3 s before its fix and 2.9 s after).
Proposed change: vectorise, as in the report.

**K13 [INFO, confidence: high]: stale local notes and a downstream bug**
- The local, untracked `BUG-gl.keep.sexlinked-zero-loci.md` and the "Known defect" section of the local `CLAUDE.md` describe the zero-loci crash as open. It was fixed upstream, and the guard is at lines 720-742.
- `gl.infer.sex()` fails with "invalid 'times' argument" when every subset is `NULL`. K3 stops that input from being produced, but the defect belongs to the `gl.infer.sex` review.
- Report F8 (gate the duplicate-sex-column warning) applies here too. It is not proposed, because it was rejected in the report review.

## Proposed changes

Changes 1-7 and 9-11 port the fixes approved for `gl.report.sexlinked()`. The numbering of the report review is kept where it applies.

1. Take sex from the `ind.metrics` rows in order (K2). **Consequence: objects whose `id` is missing or differs from `indNames(x)` get real subsets instead of four `NULL`s.**
2. Error when either sex has no individuals (K3). **Consequence: calls that currently return four `NULL`s now error.**
3. `drop = FALSE` when subsetting by sex (K4).
4. Start the cluster after validation and stop it with `on.exit()` (K5).
5. Restrict to SNP data (K6). **Consequence: SilicoDArT input errors.**
6. Fisher's test on the observed counts (K7). **Consequence: p-values, `ratio`/`stat` and the loci kept in each category change, mostly for small sample sizes (10K platypus set: 85 → 95 x-linked loci), matching `gl.report.sexlinked()` after PR #33.**
7. Documentation rewrite, including `NULL` for empty categories (K9). Docs only.
9. Apply `plot.theme` to the four plots (K10).
10. Refer to the w/y-linked column by name (K11).
11. Vectorise the classification (K12).
12. Re-subset `loc.metrics` from the original object for each returned subset (K1). **Consequence: for plain genlight input, the returned objects carry only their own loci's metrics.**
13. Append the call to `@other$history` of each returned subset (K8).

How to implement 1-6 and 9-11 is a separate choice, asked for in the approval boxes:
- (a) Copy the fixed code from PR #33 into this file. This keeps the PR self-contained, but leaves three copies to maintain.
- (b) Move the classification into one internal function (`utils.sexlinked.classify()`) used by report, keep and filter. This waits for #33 to merge and also touches `gl.report.sexlinked()`. The conventions (STY5) call for the custodian to be consulted first.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API: run.
- Modify-family checks: genotype–metadata sync on dartR and plain genlight input (run), flags (all `FALSE` on LBP already; keeping loci does not invalidate locus metrics, so no finding), history (run: K8), filter keeps what its name says (run: returned loci match the table categories).
- Spec: behaviour vs roxygen on LBP, platypus.gl, testset.gl, testset.gs, platypus_10KSNPs.gl: run.
- Downstream: `gl.infer.sex()` on keep output for LBP xy and platypus.gl zw (with `NULL` subsets): run, works. On all-`NULL` output: fails (K13).
- Engine equality with `gl.report.sexlinked()` on `origin/dev`: run, identical.
- Serial vs parallel: run, identical.
- Chi-square branch (1,000 or more sexed individuals): SKIPPED, because no fixture that size is available.
- FBM path (DAT6): SKIPPED, because no FBM fixture is available.
- dartR Google Group: not searched. Known complaint: OneDArT issue 8822 (zero-loci crash), already fixed upstream.
- Baseline timing: written after the first read of the code, before any change.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence approved |
| 2 | approved | Luis | consequence approved |
| 3 | approved | Luis | |
| 4 | approved | Luis | |
| 5 | approved | Luis | consequence approved |
| 6 | approved | Luis | consequence approved |
| 7 | approved | Luis | |
| 9 | approved | Luis | |
| 10 | approved | Luis | |
| 11 | approved | Luis | |
| 12 | approved | Luis | consequence approved |
| 13 | approved | Luis | |
| Approach | (a) copy into keep | Luis | shared helper not pursued now |

## Outcome

All approved changes are applied in `R/gl.keep.sexlinked.r`, with `man/gl.keep.sexlinked.Rd` regenerated. The classification code is copied from `gl.report.sexlinked()` at b077e6d (PR #33). The keep-specific parts are the four BEFORE/AFTER plots and the subsetting, which is now a `keep.loci()` helper.

- Equivalence: with the zero replacement put back, the new output (history removed) is identical to f293fba on LBP, testset.gl, platypus.gl and platypus_10KSNPs.gl, in both `xy` and `zw`. So every numerical diff comes from change 6.
- Consistency: `results.table` is identical to the fixed `gl.report.sexlinked()` on all four datasets and both systems. On the 10K platypus set, keep returns 95 x-linked loci, the same as the report (before: 85).
- 1, 2, 3 and 5: checked in the updated baseline. A missing `id` gives a table equal to the unmodified run. No males errors with "Found 114 females and 0 males". One female runs. `testset.gs` errors.
- 12: plain genlight input returns `$x.linked` with 66 `loc.metrics` rows whose `AlleleID` values match the kept loci (before: 1,000 rows).
- 13: each returned object gains one history entry, the `gl.keep.sexlinked()` call.
- 9: the 4 saved plots carry `theme_dartR()`.
- 11: serial time on 10K loci went from 6.4 s to 2.7 s.
- Downstream: `gl.infer.sex()` on the LBP output returns a 376 x 11 table.
- Tests: 57 expectations across both test files, 0 failures. The serial/parallel comparison removes history first, because the recorded call includes `ncores`.
- R CMD check (`--no-manual`): 0 errors. 1 WARNING comes from the environment (dependencies built under R 4.4.3). 1 NOTE (timestamp).
- Branching: the branch sits on `dev_luis` (on top of PR #33), because both PRs add `NEWS.md`, `tests/testthat.R` and the testthat Suggests line. Merge #33 first.
- API3 caller grep: dartr_shiny and dartr2shiny call `gl.keep.sexlinked()` (see PR). The list names and element order are unchanged.
- PR: green-striped-gecko/dartR.sexlinked#34 (commit 2ead5c5).

```json
{
  "function": "gl.keep.sexlinked",
  "package": "dartR.sexlinked",
  "family": "modify",
  "skill_version": "2.0.0",
  "commit": "f293fba",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "K1", "severity": "HIGH", "confidence": "high", "rule": "DAT3", "status": "approved", "change": 12},
    {"id": "K2", "severity": "HIGH", "confidence": "high", "rule": "DAT5", "status": "approved", "change": 1},
    {"id": "K3", "severity": "HIGH", "confidence": "high", "rule": "FS5", "status": "approved", "change": 2},
    {"id": "K4", "severity": "HIGH", "confidence": "high", "rule": "DAT5", "status": "approved", "change": 3},
    {"id": "K5", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 4},
    {"id": "K6", "severity": "MEDIUM", "confidence": "high", "rule": "DAT1", "status": "approved", "change": 5},
    {"id": "K7", "severity": "MEDIUM", "confidence": "medium", "rule": "none", "status": "approved", "change": 6},
    {"id": "K8", "severity": "LOW", "confidence": "high", "rule": "FS8", "status": "approved", "change": 13},
    {"id": "K9", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 7},
    {"id": "K10", "severity": "LOW", "confidence": "high", "rule": "PLT1", "status": "approved", "change": 9},
    {"id": "K11", "severity": "LOW", "confidence": "high", "rule": "STY1", "status": "approved", "change": 10},
    {"id": "K12", "severity": "LOW", "confidence": "high", "rule": "STY2", "status": "approved", "change": 11},
    {"id": "K13", "severity": "INFO", "confidence": "high", "rule": "none", "status": "noted", "change": null}
  ],
  "coverage_skipped": ["chi-square branch: no >=1000-individual fixture", "DAT6: no FBM fixture", "dartR Google Group not searched"],
  "status": "pr-open",
  "pr": 34
}
```
