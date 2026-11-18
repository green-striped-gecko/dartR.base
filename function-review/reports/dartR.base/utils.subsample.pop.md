# Review: utils.subsample.pop (dartR.base)

- Family mode: analysis (internal utility; with its helper `het_rep()`)
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: origin/dev after PR #435 (which added the `ploidy` argument and gametic Ho for polyploids)
- Datasets: testset.gl and subsets; FBM copy via `gl.gen2fbm`
- Baseline: `tests/testthat/test-utils.subsample.pop.R` (11 expectations, captured pre-review, all pass)
- Callers: `gl.report.heterozygosity` and `gl.report.polyploid_heterozygosity` (`subsample.pop = TRUE`); none in sibling packages

## Verdict

**Standards: Needs work** — no description of arguments or return value; one subsetting call lacks `drop = FALSE`.
**Spec: Needs work** — correct where it runs (a whole-population subsample reproduces the population's Ho; FBM input gives identical results), but any `n.limit` below 10, which users set through the report functions, stops with an error.

## Findings

**F1 [HIGH, confidence: high] — `n.limit` below 10 errors (spec)**
`R/utils.het.report.r` `het_rep()` — every population that passes `n.limit` is subsampled at every size in `subsamples = c(10, 5, 4, 3, 2)`, drawing without replacement. With `n.limit = 5`, a population of 7 is kept but cannot provide a subsample of 10: "cannot take a sample larger than the population when 'replace = FALSE'". Verified through `gl.report.heterozygosity(testset.gl, subsample.pop = TRUE, n.limit = 5)`; testset.gl has 7 populations of 5-9 individuals.
Failure scenario: a user lowers `n.limit` to include small populations, as the parameter invites, and the report stops.
Proposed change: for each population, use only the subsample sizes it can provide (at most its number of individuals); a population that can provide none is skipped like one below `n.limit`.
**Consequence: calls with `n.limit` below 10 now return results (small populations get fewer subsample sizes) instead of an error; results for `n.limit >= 10` are unchanged.**

**F2 [LOW, confidence: high] — single-locus input errors (DAT3)**
`het_rep()`, diploid branch — `mat[rows, ] == 1` drops to a vector when there is one locus, and `colMeans()` fails ("'x' must be an array of at least two dimensions"). The polyploid branch already uses `drop = FALSE`.
Proposed change: add `drop = FALSE`. Results are unchanged for more than one locus.

**F3 [LOW, confidence: high] — undocumented internal interface (STY1)**
Neither `utils.subsample.pop()` nor `het_rep()` states its arguments, return columns (`res.mean`, `res_SE`, `pop`, `subsample`) or that `res_SE` is the standard error across the 10 replicates.
Proposed change: add a comment header to each. No code change.

## Proposed changes

1. Subsample each population only at the sizes it can provide (F1). **Consequence: `n.limit < 10` returns results instead of an error; unchanged for `n.limit >= 10`.**
2. `drop = FALSE` in the diploid branch of `het_rep()` (F2).
3. Comment headers for `utils.subsample.pop()` and `het_rep()` (F3).

## Coverage

- Standards walk: DAT, STY, DEP — run. FS, VRB, DOC, PLT: not applicable (internal, no verbosity, no roxygen, no plot).
- Spec: whole-population subsample vs independent Ho — run, equal; SE 0 as expected
- FBM path (DAT6): run, `identical()` to dense for the same seed
- Diploid output unchanged by #435 — verified in that PR (identical to the pre-#435 helper for the same seed)
- Edge cases: `n.limit` 5 (error, F1), 1000 (empty table, correct), single locus (error, F2)
- Polyploid path: covered by the #435 test (tetraploid subsample of all individuals equals gametic Ho)

## Report notes (not findings)

- A subsample whose size equals the population size (for example 10 of 10) draws the same individuals in every replicate, so that row is the population's Ho with SE 0. This follows from the Schmidt et al. (2021) design; no change proposed.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | |
| 2 | approved | Luis | |
| 3 | approved | Luis | |

## Outcome

- Changes 1-3 applied on branch `review-utils.subsample.pop` (from `origin/dev`, which includes #435 and #436).
- Characterization test: 15 expectations pass. Diffs from baseline map to approved changes only: `n.limit = 5` returns all 27 populations of at least 5 individuals, each with the sizes it can supply (1); a single locus returns 100 rows (2).
- For `n.limit = 10`, output is `identical()` to the pre-review function for three seeds on testset.gl and for platypus.gl.
- Both reports run end to end with plots at `verbose = 3` with `n.limit` 5 (`gl.report.heterozygosity`) and 3 (`gl.report.polyploid_heterozygosity`).
- Related tests pass: `gl.report.heterozygosity` 37, `gl.report.polyploid_heterozygosity` 37, `gl.test.heterozygosity` 38.
- Full `R CMD check` not run.
- PR: pending.

## Machine block

```json
{
  "function": "utils.subsample.pop",
  "package": "dartR.base",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "origin/dev after #435",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "spec", "status": "approved", "change": 1},
    {"id": "F2", "severity": "LOW", "confidence": "high", "rule": "DAT3", "status": "approved", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "STY1", "status": "approved", "change": 3}
  ],
  "coverage_skipped": [],
  "status": "pr-open",
  "pr": null
}
```
