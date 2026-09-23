# Review: gl.filter.sexlinked (dartR.sexlinked)

- Family mode: modify
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: `R/gl.filter.sexlinked.r` as of f293fba (version 1.2.2); branch `review-gl.filter.sexlinked`, stacked on PR #34 (which is stacked on #33)
- Datasets: LBP, testset.gl, testset.gs, platypus.gl (dartR.data 1.2.5); platypus_10KSNPs.gl (local, 200 ind x 10,000 loci)
- Baseline: `tests/testthat/test-gl.filter.sexlinked.R` (20 expectations, all pass on the unchanged function)
- Related: `gl.report.sexlinked` (PR #33) and `gl.keep.sexlinked` (PR #34) reviews

## Verdict

**Standards: Needs work**: the function carries the same classification code as report and keep (a line diff against keep at f293fba shows only formatting, plot titles and the output block). It has the same defects, plus unsynced locus metadata for plain `genlight` input and no history entry.

**Spec: Needs work**: with no males or a missing `id` column, it returns the whole input unfiltered and raises no error. Every sex-linked locus then flows silently into downstream population-genetic analyses.

What works well: the returned loci are exactly the autosomal loci of the classification table. The all-sex-linked case returns `NULL` instead of crashing (b97c1cf). The input object is not modified.

## Findings

Findings F1-F6 and F9-F11 of the report review, and K1 and K8 of the keep review, apply here. They are restated with filter-specific consequences.

**L1 [HIGH, confidence: high]: nothing filtered when sexes cannot be matched (FS5, DAT5); report F1-F2**
`R/gl.filter.sexlinked.r:180-209`: with no males, or no matching `id`, every locus is classified as autosomal and returned.
Failure scenario: `testset.gl` with no `id` column, or with every male set to `Unknown`, returns all 255 loci. The message reads "Filtered out 0 sex-linked loci", and the call does not error.
Proposed change: as report changes 1 and 2: read sex by row, and error when either sex has no individuals.

**L2 [HIGH, confidence: high]: crash when a sex has exactly one individual (DAT5); report F3**
Failure scenario: `testset.gl` with one female fails with `'x' must be an array of at least two dimensions`.
Proposed change: `drop = FALSE`.

**L3 [HIGH, confidence: high]: `loc.metrics` not subset for plain genlight input (DAT2, DAT3); keep K1**
`R/gl.filter.sexlinked.r:624-628`.
Failure scenario: `as(LBP, "genlight")` returns 923 loci with 1,000 `loc.metrics` rows.
Proposed change: `gl.autosomal@other$loc.metrics <- x@other$loc.metrics[autosomal, , drop = FALSE]`.

**L4 [MEDIUM, confidence: high]: parallel workers left running after an error (FS5); report F4**
Failure scenario: `system = NULL, ncores = 2` leaves 2 connections open (0 → 2).
Proposed change: `on.exit()` after validation.

**L5 [MEDIUM, confidence: high]: SilicoDArT accepted (DAT1); report F5**
Failure scenario: `testset.gs` returns a filtered presence/absence object, with loci removed partly on a meaningless "heterozygosity".
Proposed change: SNP data only.

**L6 [MEDIUM, confidence: medium]: zeros replaced by 1 before Fisher's test (numerical correctness); report F6**
Failure scenario: after #33 and #34, report and keep use the observed counts and filter does not. On platypus_10KSNPs.gl, filter keeps 9,893 loci, while report and keep class 117 loci as sex-linked, so 10 of them would stay in the "autosomal" output.
Proposed change: Fisher's test on the observed counts.

**L7 [LOW, confidence: high]: no history entry (FS8); keep K8**
Failure scenario: `r@other$history` has the same length as the input's, although loci were removed.
Proposed change: append `match.call()` to the returned object.

**L8 [LOW, confidence: high]: documentation disagrees with behaviour (DOC5 (proposed rule), DOC1, DOC2, DOC7 (proposed rule))**
- `@return` says "A genlight object and 4 plots", but the plots are not returned, and the function returns `NULL` when every locus is sex-linked.
- The `id` requirement is stated, and the sex codes and thresholds are undocumented.
- There is no `@family` and no `Author(s):` line.

Failure scenario: a script does `nLoc(gl.filter.sexlinked(x, "xy"))` and errors on an all-sex-linked panel, which the help page does not warn about.
Proposed change: a docs rewrite matching keep's, stating the `NULL` case.

**L9 [LOW, confidence: high]: `plot.theme` ignored (PLT1); report F9**
Proposed change: apply it to the four plots.

**L10 [LOW, confidence: high]: column 11 referenced by position, row-by-row loops (STY1, STY2); report F10-F11**
Failure scenario: 10K loci take 6.3 s serial (keep after the fix: 2.7 s).
Proposed change: port the vectorised, by-name code.

**L11 [INFO, confidence: high]: return contract for the all-sex-linked case**
`NULL` was chosen in b97c1cf to match keep's empty-category contract. That is a recent deliberate decision, so it is not flagged as a defect. It is documented under L8.

## Proposed changes

The numbering follows the keep review.

1. Take sex from the `ind.metrics` rows in order (L1). **Consequence: objects whose `id` is missing or differs from `indNames(x)` are filtered, instead of returned whole.**
2. Error when either sex has no individuals (L1). **Consequence: calls that currently return the unfiltered object now error.**
3. `drop = FALSE` (L2).
4. Start the cluster after validation and stop it with `on.exit()` (L4).
5. Restrict to SNP data (L5). **Consequence: SilicoDArT input errors.**
6. Fisher's test on the observed counts (L6). **Consequence: the loci removed change, mostly for small sample sizes (10K platypus set: 9,893 → 9,883 loci kept), matching report and keep.**
7. Documentation rewrite, including `NULL` when every locus is sex-linked (L8). Docs only.
9. Apply `plot.theme` to the four plots (L9).
10. Refer to the w/y-linked column by name (L10).
11. Vectorise the classification (L10).
12. Re-subset `loc.metrics` from the input (L3). **Consequence: for plain genlight input, the returned object carries only its own loci's metrics.**
13. Append the call to `@other$history` (L7).

Implementation: copy the code from keep (#34), as decided for keep.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API: run.
- Modify-family checks: metadata sync on dartR and plain genlight input (run), flags (already `FALSE` on the inputs; keeping individuals intact does not invalidate locus metrics), history (run: L7), the filter removes what its name says (run: the returned loci equal the autosomal loci of the table).
- Edge cases: no males, no `id`, one female, SilicoDArT, all loci sex-linked, missing `system` with `ncores = 2`: run.
- Serial vs parallel: run, identical.
- Chi-square branch: SKIPPED, because no fixture with 1,000 or more individuals is available.
- FBM path (DAT6): SKIPPED, because no FBM fixture is available.
- Caller grep (API3): see Outcome.
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

## Outcome

All approved changes are applied in `R/gl.filter.sexlinked.r`, with `man/gl.filter.sexlinked.Rd` regenerated. The code is copied from `gl.keep.sexlinked()` at 2ead5c5 (PR #34). The filter-specific parts are the AFTER plots, which show autosomal loci only, and the single-object output.

- Equivalence: with the zero replacement put back, the new output (history removed) is identical to f293fba on LBP, testset.gl, platypus.gl and platypus_10KSNPs.gl, in both `xy` and `zw`.
- Consistency: on all four datasets and both systems, the returned loci are exactly the input minus the loci returned by `gl.keep.sexlinked()`. On the 10K platypus set, 9,893 → 9,883 loci are kept (change 6).
- 1, 2, 3 and 5: checked in the updated baseline. A missing `id` gives the same loci as the unmodified run. No males errors. One female runs. `testset.gs` errors.
- 12: plain genlight input gives 923 `loc.metrics` rows for 923 loci (before: 1,000).
- 13: one history entry, the `gl.filter.sexlinked()` call.
- 9: the 4 saved plots carry `theme_dartR()`.
- 11: serial time on 10K loci went from 6.3 s to 2.9 s.
- Tests: 78 expectations across three test files, 0 failures.
- R CMD check (`--no-manual`): 0 errors. 1 WARNING comes from the environment (dependencies built under R 4.4.3). 1 NOTE (timestamp).
- API3 caller grep: dartr_shiny and dartr2shiny call `gl.filter.sexlinked(x, system, ncores, plot.theme)` and assign the result to `MyData`. The signature is unchanged. Datasets with no males now raise an error in the app instead of silently replacing `MyData` with an unfiltered copy.
- Branching: stacked on #34 (itself stacked on #33). Merge in the order #33, #34, then this PR.
- PR: green-striped-gecko/dartR.sexlinked#35 (commit d0cfd22).

```json
{
  "function": "gl.filter.sexlinked",
  "package": "dartR.sexlinked",
  "family": "modify",
  "skill_version": "2.0.0",
  "commit": "f293fba",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "L1", "severity": "HIGH", "confidence": "high", "rule": "FS5", "status": "approved", "change": 1},
    {"id": "L2", "severity": "HIGH", "confidence": "high", "rule": "DAT5", "status": "approved", "change": 3},
    {"id": "L3", "severity": "HIGH", "confidence": "high", "rule": "DAT3", "status": "approved", "change": 12},
    {"id": "L4", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 4},
    {"id": "L5", "severity": "MEDIUM", "confidence": "high", "rule": "DAT1", "status": "approved", "change": 5},
    {"id": "L6", "severity": "MEDIUM", "confidence": "medium", "rule": "none", "status": "approved", "change": 6},
    {"id": "L7", "severity": "LOW", "confidence": "high", "rule": "FS8", "status": "approved", "change": 13},
    {"id": "L8", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 7},
    {"id": "L9", "severity": "LOW", "confidence": "high", "rule": "PLT1", "status": "approved", "change": 9},
    {"id": "L10", "severity": "LOW", "confidence": "high", "rule": "STY2", "status": "approved", "change": 11},
    {"id": "L11", "severity": "INFO", "confidence": "high", "rule": "none", "status": "noted", "change": null}
  ],
  "coverage_skipped": ["chi-square branch: no >=1000-individual fixture", "DAT6: no FBM fixture"],
  "status": "pr-open",
  "pr": 35
}
```
