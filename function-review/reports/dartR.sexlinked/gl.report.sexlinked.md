# Review: gl.report.sexlinked (dartR.sexlinked)

- Family mode: report
- Date: 2026-09-23
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: f293fba (version 1.2.2, loaded with `devtools::load_all()`)
- Datasets: LBP, testset.gl, testset.gs, platypus.gl (dartR.data 1.2.5); platypus_10KSNPs.gl (local, 200 ind x 10,000 loci) for timing
- Baseline: `tests/testthat/test-gl.report.sexlinked.R` (21 expectations, all pass on f293fba)

## Verdict

**Standards: Needs work**: the structure follows the dartR anatomy, but a parallel cluster leaks on error, one warning ignores `verbose`, and several documented parameters do nothing.

**Spec: Needs work**: when the function cannot match individuals to a sex, it still returns a table that classifies every locus as autosomal, with no error. With exactly one individual of a sex it crashes with an unhelpful error.

What works well: the classification logic gives identical results in serial and parallel runs, and the input object comes back unchanged (FS8).

## Findings

**F1 [HIGH, confidence: high]: individuals matched through an `id` column the function never checks (DAT5)**
`R/gl.report.sexlinked.r:216-231`: sexes are assigned by matching `ind.metrics$id` against `indNames(x)`. If `id` is missing, or its values differ from `indNames(x)`, no individual matches.
Failure scenario: `testset.gl` with the `id` column removed, or with `id` values renamed, gives `count.F.scored == 0` and `count.M.scored == 0` at every locus. The printed summary reads "0 sex-linked loci ... 255 autosomal loci", and the call does not error.
Proposed change: take sex from the `ind.metrics` rows in individual order, since they track individuals 1:1 (DAT2), and stop matching through `id`.

**F2 [HIGH, confidence: high]: no guard for a missing sex (FS5)**
`R/gl.report.sexlinked.r:209-217`: the only check is that at least one F or M exists overall.
Failure scenario: `testset.gl` with every male set to `Unknown` prints "Detected 114 females and 0 males", then classifies all 255 loci as autosomal. `scoringRate.M` is `NaN` at every locus.
Proposed change: `stop(error(...))` when either sex has zero individuals, and name the counts found.

**F3 [HIGH, confidence: high]: crash when a sex has exactly one individual (DAT5)**
`R/gl.report.sexlinked.r:230-231`: `gen[, cols]` without `drop = FALSE` turns a one-column selection into a vector.
Failure scenario: `testset.gl` with one female fails with `'x' must be an array of at least two dimensions`, raised from `rowSums()`.
Proposed change: add `drop = FALSE` to both subsets. Whether one individual per sex should be allowed at all depends on the minimum decided under change 2.

**F4 [MEDIUM, confidence: high]: parallel workers left running after an error (FS5, FS6)**
`R/gl.report.sexlinked.r:159-162, 705-707`: the cluster starts before `system` is validated and stops only at the end of the function.
Failure scenario: `gl.report.sexlinked(testset.gl, system = NULL, ncores = 2)` errors as designed but leaves 2 worker connections open (`showConnections()` 0 → 2). Any error after line 162 does the same.
Proposed change: start the cluster after parameter validation and register `on.exit(parallel::stopCluster(cl))`.

**F5 [MEDIUM, confidence: high]: SilicoDArT data accepted and processed as SNP (DAT1)**
`R/gl.report.sexlinked.r:152-157`: `accept` includes `"SilicoDArT"`, but phase 2 counts `1` as heterozygous. In presence/absence data, `1` means presence.
Failure scenario: `testset.gs` returns a 255 x 23 table with "heterozygosity" columns (mean `heterozygosity.F` 0.344) that have no genetic meaning, and gives no warning.
Proposed change: accept SNP data only and error clearly for SilicoDArT.

**F6 [MEDIUM, confidence: medium]: zeros replaced by 1 before every test (no rule fits: numerical correctness)**
`R/gl.report.sexlinked.r:275, 288, 330, 482, 493, 542`: every 0 in the 2x2 tables becomes 1 before `fisher.test()`/`chisq.test()`. The comment says this avoids an error, but `fisher.test()` accepts zeros.
Failure scenario: 3 females all scored and 3 males all missing gives p = 0.10 on the real counts and p = 0.49 after replacement. With small sample sizes, loci that are truly sex-linked lose significance. On LBP (373 sexed individuals) the replacement changes nothing: 10 loci with `p.adjusted <= 0.01` either way.
Proposed change: needs the custodian. Either keep the replacement for the chi-square branch only and run Fisher on the real counts, or keep it and document it as intended. It was not checked whether the published method (Robledo-Ruiz et al. 2023) specifies the replacement.

**F7 [LOW, confidence: high]: documentation disagrees with behaviour (DOC5 (proposed rule), DOC1, DOC2, DOC7 (proposed rule))**
`R/gl.report.sexlinked.r:1-110`:
- `@title` says "Filters loci", but the function reports and filters nothing.
- `@return` promises "a dataframe and 2 plots", but only the data frame is returned; the plots are printed or saved.
- `ratio` and `stat` are described as Fisher estimates, but they hold chi-square statistics when 1,000 or more sexed individuals are present.
- `MALE`/`FEMALE` values are accepted but not documented.
- The 0.1 call-rate threshold and the 0.01 FDR threshold are not documented.
- `@family` is missing, the `verbose` text is non-standard, and there is no `Author(s):` line.

Failure scenario: a user reading the help page expects to get the plots back, or expects `ratio` to be an odds ratio in a large study.
Proposed change: a docs-only rewrite of these items, then `devtools::document()`.

**F8 [LOW, confidence: high]: warning not gated by `verbose` (VRB3)**
`R/gl.report.sexlinked.r:194-200`: the multiple-"sex"-columns warning prints at `verbose = 0` and has no trailing newline, so the next output runs onto the same line.
Failure scenario: `platypus.gl` has both `Sex` and `sex` columns and prints the warning at `verbose = 0`.
Proposed change: wrap in `if (verbose >= 2)` and add `\n`.

**F9 [LOW, confidence: high]: `plot.theme` and `plot.colors` accepted but ignored (PLT1)**
`R/gl.report.sexlinked.r:116-117, 131-143, 409-652`: the colour-validation code runs, and then neither argument reaches the plots.
Failure scenario: `plot.theme = theme_bw()` has no effect, and no message says so.
Proposed change: add `plot.theme` to both plots. Leave `plot.colors` documented as not implemented, because the plots use five category colours.

**F10 [LOW, confidence: high]: column referenced by position (STY1)**
`R/gl.report.sexlinked.r:389, 465, 514`: `table[i, 11]` is the w/y-linked flag, which the baseline confirms is correct today. Inserting a column before position 11 would silently change which loci get excluded.
Proposed change: refer to the column by name.

**F11 [LOW, confidence: high]: row-by-row loops; parallel slower than serial (STY2)**
Lines 307-340, 359-394, 512-609: each locus is handled with `table[i, ...]` reads and writes in a loop.
Failure scenario: platypus_10KSNPs.gl (10,000 loci) takes 6.2 s serial and 9.9 s with `ncores = 2`.
Proposed change: vectorise the classification steps (`w/y.linked`, `sex.biased`, `x/z.linked`, `gametolog`). The results should be unchanged, and the baseline test checks that. Recommend deferring until F1-F6 are settled.

**F12 [INFO, confidence: high]: engine duplicated in `gl.keep.sexlinked()` and `gl.filter.sexlinked()`**
About 300 lines of the classification code are copied into both functions (77 diff lines between the keep and filter copies). F1-F4, F6 and F10 probably apply to them too, but that was not tested here. Fixes made here do not reach those copies.
Proposed change: none in this PR. Carry these findings into their reviews, or extract a shared internal helper later (a team decision).

## Proposed changes

1. Take sex from `ind.metrics` in row order instead of matching `id` (F1). **Consequence: numerical output changes for objects whose `id` column is missing or differs from `indNames(x)`: they get real counts instead of all-zero counts.**
2. Error when either sex has no individuals (F2). **Consequence: calls that currently return an all-autosomal table now error.**
3. Add `drop = FALSE` when subsetting genotypes by sex (F3).
4. Start the cluster after validation and stop it with `on.exit()` (F4).
5. Restrict to SNP data (F5). **Consequence: `gl.report.sexlinked()` on SilicoDArT data errors instead of returning a table.**
6. Zero replacement: run Fisher on the real counts, or document the replacement as intended (F6). **Consequence (if the replacement is removed): p-values, `ratio`/`stat` and locus classifications change, most for small sample sizes.**
7. Documentation rewrite (F7). Docs only.
8. Gate the multiple-"sex"-column warning by `verbose` and add a newline (F8).
9. Apply `plot.theme` to both plots (F9).
10. Refer to the w/y-linked column by name (F10).
11. Vectorise the classification loops (F11). Recommended to defer.

Adding the test file also needs `Suggests: testthat (>= 3.0.0)` and `Config/testthat/edition: 3` in `DESCRIPTION`, which goes in the PR with the test.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API: run.
- Report-family checks: input untouched (confirmed with `identical()`), no history append (confirmed), results independent of plotting (plots are built but do not feed the table): run.
- Spec: behaviour vs roxygen on LBP, testset.gl, testset.gs, platypus.gl: run.
- Edge cases (one individual of a sex, no males, no `id`, `id` differing from `indNames`, SilicoDArT, missing `system` with `ncores = 2`): run.
- Serial vs parallel equality: run, identical.
- Chi-square branch (1,000 or more sexed individuals): SKIPPED, because no fixture that size is available. Its behaviour was read from the code, not run.
- FBM path (DAT6): SKIPPED, because no FBM fixture is available. `as.matrix()` plus `t()` plus `as.data.frame()` densifies the full matrix, which was not measured on large data.
- Published-method check for F6: SKIPPED, because the paper was not consulted.
- Known complaints: GitHub issues on green-striped-gecko/dartR.sexlinked are empty or disabled. The dartR Google Group was not searched.
- Baseline timing: the baseline was written after the first read of the code, not before, but before any change.
- Environment note: Homebrew `Rscript` (R 4.6) has no dartR packages. All runs used `/usr/local/bin/Rscript` (R 4.4.2).

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence approved |
| 2 | approved | Luis | consequence approved |
| 3 | approved | Luis | |
| 4 | approved | Luis | |
| 5 | approved | Luis | consequence approved |
| 6 | approved | Luis | "Use real counts": Fisher on the observed counts; chi-square branch keeps the replacement. Consequence approved |
| 7 | approved | Luis | |
| 8 | rejected | Luis | not selected in the approval boxes |
| 9 | approved | Luis | |
| 10 | approved | Luis | |
| 11 | approved | Luis | approved despite the recommendation to defer |

## Outcome

All approved changes are applied in `R/gl.report.sexlinked.r`, with `man/gl.report.sexlinked.Rd` regenerated.

- 1: sex is read from the `ind.metrics` rows, with a new error when the row count differs from `nInd(x)`. Evidence: removing or renaming `id` in `testset.gl` gives a result equal to the unmodified run.
- 2: errors with "Found 114 females and 0 males; at least one of each is needed" (`testset.gl` with no males).
- 3: one female runs (`testset.gl`).
- 4: the cluster starts after validation and stops through `on.exit()`. Evidence: an error injected right after the cluster starts leaves 0 connections open (before: 2).
- 5: `testset.gs` errors with "found SilicoDArT expecting dartR or genlight or SNP".
- 6: Fisher's test uses the observed counts. Evidence: on the 10K platypus set, 85 → 95 x-linked loci (10 gained, 0 lost). y-linked, sex-biased and gametolog counts are unchanged. Category counts on LBP, testset.gl and platypus.gl are unchanged; `ratio`, `p.value` and `stat` values differ.
- 7: docs rewritten: title, return value, sex coding, row-order requirement, classification thresholds, test switch at 1,000 individuals, `@family matched report`, `Author(s):` line.
- 9: `plot.theme` is added to both plots. Evidence: the saved RDS carries `theme_dartR()` by default and `theme_bw()` when passed.
- 10 and 11: classification is vectorised and columns are addressed by name. Evidence: with the zero replacement put back into the new code, its output is identical to f293fba on LBP, testset.gl, platypus.gl and platypus_10KSNPs.gl under both `xy` and `zw`. So every numerical diff comes from change 6. Serial time on 10K loci: 6.3 s → 2.9 s. `ncores = 2` still takes 8.9 s, slower than serial, because each locus is a separate foreach task.
- Snapshot: the baseline edge-case expectations were updated only for changes 1, 2, 3 and 5. The LBP, testset.gl and platypus.gl category snapshots pass unchanged. New tests cover change 6 (the p-value equals `fisher.test()` on the observed counts) and change 4.
- R CMD check (`--no-manual`): 0 errors. 1 WARNING comes from the environment (ade4, ggplot2 and dplyr built under R 4.4.3). Notes: timestamp, and the top-level `function-review` directory (resolved when the report moves).
- API3 caller grep: dartr_shiny and dartr2shiny call `gl.report.sexlinked(x, system, ncores, plot.theme)` and read output column 11 by position. The column order is unchanged, and `plot.theme` now takes effect in the app. No sibling dartR.* package calls it.
- `DESCRIPTION`: `Suggests: testthat (>= 3.0.0)`, `Config/testthat/edition: 3`. Adds `NEWS.md`.
- Correction to step 1: the campaign manifest is `dartR.base/function-review/manifest.csv`, where this function was already listed as `pending`. A duplicate manifest created in dartR.sexlinked during the review was not committed.
- PR: green-striped-gecko/dartR.sexlinked#33 (commit b077e6d on `dev_luis`).

```json
{
  "function": "gl.report.sexlinked",
  "package": "dartR.sexlinked",
  "family": "report",
  "skill_version": "2.0.0",
  "commit": "f293fba",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "FS5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "DAT5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DAT1", "status": "approved", "change": 5},
    {"id": "F6", "severity": "MEDIUM", "confidence": "medium", "rule": "none", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "VRB3", "status": "rejected", "change": 8},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "PLT1", "status": "approved", "change": 9},
    {"id": "F10", "severity": "LOW", "confidence": "high", "rule": "STY1", "status": "approved", "change": 10},
    {"id": "F11", "severity": "LOW", "confidence": "high", "rule": "STY2", "status": "approved", "change": 11},
    {"id": "F12", "severity": "INFO", "confidence": "high", "rule": "none", "status": "noted", "change": null}
  ],
  "coverage_skipped": ["chi-square branch: no >=1000-individual fixture", "DAT6: no FBM fixture", "F6 published method not consulted", "dartR Google Group not searched"],
  "status": "pr-open",
  "pr": 33
}
```
