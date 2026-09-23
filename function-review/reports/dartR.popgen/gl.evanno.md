# Review: gl.evanno + utils.structure.evanno (dartR.popgen)
- Family mode: analysis (Evanno delta K from STRUCTURE runs)
- Scope: one review for the exported wrapper `gl.evanno` (3 lines) and the
  helper that holds all the logic, `utils.structure.evanno`; two manifest
  rows, one report, one PR (agreed with Luis, 2026-09-23)
- Date: 2026-09-23
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 95fde36 (origin/dev)
- Datasets: hand-built `structure.result` objects with controlled LnP(K);
  a real STRUCTURE 2.3.4 run (`~/programs/structure`, testset.gl three
  populations, K = 1:4, five replicates)
- Baseline: tests/testthat/test-gl.evanno.R (snapshot captured pre-review;
  runs without the STRUCTURE binary)

## Verdict

**Standards: Needs work** — the wrapper has no `verbose`, flags or checks;
the helper uses the deprecated `aes_string()`, draws with `gridExtra`
without a guard, and ignores `theme_dartR()` and the `plot.file` idiom.
**Spec: Needs work** — delta K is wrong or unusable in two common cases:
replicates that converge to the same LnP(K) give `Inf` (it happened on a
real STRUCTURE run of testset.gl), and a `k.range` with gaps is differenced
as if the K values were consecutive.

What works well: on consecutive K with varying replicates, `mean.ln.k`,
`sd.ln.k`, LnP'(K), |LnP''(K)| and delta K match a hand calculation
exactly, and K values of 10 and more keep numeric order.

## Findings

**F1 [MEDIUM, confidence: high] — identical replicates give delta K = `Inf` (DOC5)**
`R/utils.structure.evanno.r:33-35` — delta K divides by `sd.ln.k`, which is
0 when all replicates of a K report the same LnP(K).
Failure scenario: `gl.run.structure` on three populations of testset.gl
(K = 1:4, five replicates, 200 burn-in and 200 MCMC): every replicate of
K = 2, 3 and 4 has the same LnP(K), so `delta.k` is `Inf` at K = 2 and 3
and the delta K panel has no finite values. Nothing tells the user why.
Proposed change: set delta K to `NA` where `sd.ln.k` is 0 and say at
`verbose >= 1` that delta K is undefined for those K because the replicates
agree (and that LnP(K) should be read instead).

**F2 [MEDIUM, confidence: high] — K values with gaps are differenced as if consecutive (DOC5)**
`R/utils.structure.evanno.r:31-35` — `diff()` and the `i - 1`, `i + 1`
neighbours work on positions, not on K values.
Failure scenario: `k.range = c(1, 3, 5, 7)` reports LnP'(3) = LnP(3) −
LnP(1) and delta K at K = 3 from K = 1 and K = 5 (36 in the fixture), all
labelled as if they were the Evanno statistics at K = 3.
Proposed change: compute LnP'(K), |LnP''(K)| and delta K only where K − 1
(and K + 1) are present, `NA` elsewhere, with a `verbose >= 1` warning
naming the K values that lack neighbours.

**F3 [LOW, confidence: high] — delta K silently missing with one replicate (VRB4, proposed rule)**
`R/utils.structure.evanno.r:27-30, 83-86` — with one replicate, `sd` is
`NA`, so delta K is `NA` and the delta K plot is dropped without a message.
Failure scenario: `num.k.rep = 1` (the `gl.run.structure` default is 1)
returns a three-panel plot and a `delta.k` column of NAs; the user is not
told that delta K needs replicates.
Proposed change: warn at `verbose >= 1` that delta K needs at least two
replicates per K.

**F4 [LOW, confidence: high] — misleading error messages (VRB2)**
`R/utils.structure.evanno.r:15-21` — fewer than three K values stops with
"must have at least two values of k"; a wrong class says "not a result from
'structure.run'" (a strataG function name).
Failure scenario: `k.range = 2:3` gives the message "at least two values"
although two were supplied.
Proposed change: "needs at least three values of K" and "is not a
structure.result object returned by gl.run.structure", checked in
`gl.evanno` before work starts.

**F5 [LOW, confidence: high] — wrapper lacks the house structure (FS2, FS3, FS9, PLT1, PLT2)**
`R/gl.evanno.r:41-45` — no `verbose`, flags, `plot.theme`, `plot.dir` or
`plot.file`; plots use the default ggplot theme.
Failure scenario: the Evanno plot cannot be saved with the dartR
`plot.file` idiom, and its look differs from every other dartR plot
(including the Evanno panel `gl.run.structure` draws with `theme_dartR()`).
Proposed change: add `plot.theme = theme_dartR()`, `plot.dir`,
`plot.file` and `verbose` after `plot.out`; standard flags;
`utils.plot.save()` for the combined plot.

**F6 [LOW, confidence: high] — `gridExtra` used without a guard; plot not returned (DEP1, PLT1)**
`R/utils.structure.evanno.r:95-104` — the four panels are aligned with
`gridExtra::grid.arrange` (Suggests, unguarded), and the combined figure is
not returned.
Failure scenario: without `gridExtra` installed, `gl.evanno(sr)` stops with
"there is no package called 'gridExtra'" after the statistics are computed.
Proposed change: build the combined figure with patchwork (Imports), as
`gl.run.structure` already does, and return it as `plots$combined`;
`gridExtra` is then no longer used here. The four per-panel plots keep
their names, which `gl.run.structure` relies on.

**F7 [LOW, confidence: high] — deprecated `aes_string()` (PLT1)**
`R/utils.structure.evanno.r:53-84` — ggplot2 deprecated `aes_string()` in
3.0.0.
Failure scenario: the first call in a session prints a lifecycle warning
asking the user to report the issue to the dartR group.
Proposed change: `aes(x = .data$k, y = .data$mean.ln.k)` and so on.

**F8 [LOW, confidence: high] — documentation describes another function (DOC1, DOC5, DOC7)**
`R/gl.evanno.r:1-39`, `R/utils.structure.evanno.r:1-11` —
- `@description` says the function "takes a genlight object and runs a
  STRUCTURE analysis"; it takes a `structure.result` and computes the
  Evanno statistics.
- `@details` calls it a wrapper around strataG's `evanno`; strataG is no
  longer on CRAN and the code is now in `utils.structure.evanno`.
- `@return` says "a list of all four plots"; it returns `df` and `plots`,
  and there are three plots when delta K cannot be computed.
- The `df` columns are not documented; the method (means of LnP(K) per K,
  as in strataG) and the need for consecutive K and replicates are not
  stated.
- `@seealso` ends with a dangling comma; a stale strataG `@importFrom`
  comment remains; `@family` missing; `@author` lacks Author(s)/Custodian
  (DOC7, proposed rule). The helper's roxygen has no `@name`/`@title` in
  house order.
Failure scenario: a user reading the help page expects a genlight input
and four plots.
Proposed change: rewrite both headers; regenerate Rd.

**F9 [INFO, confidence: medium] — delta K is computed from mean LnP(K) (principle: method transparency)**
`R/utils.structure.evanno.r:31-35` — the second difference is taken on the
mean LnP(K) of each K, as in strataG, rather than averaging per-replicate
|L''(K)| values.
Failure scenario: none observed; values can differ from programs that
average per replicate.
Proposed change: none in code; state the method in `@details` (part of
change 8).

**F10 [INFO, confidence: high] — helper exported (FS1)**
`NAMESPACE:36` — `utils.structure.evanno` is exported although it is a
`utils.*` helper; `gl.run.structure` is its only internal caller.
Failure scenario: none; users can call an undocumented helper.
Proposed change: none — unexporting is an API change (API3) to decide
across the `utils.*` exports, not here.

## Proposed changes

1. `delta.k` is `NA` where `sd.ln.k` is 0, with a `verbose >= 1` warning
   (F1). **Consequence: `delta.k` changes from `Inf` to `NA` whenever all
   replicates of a K agree.**
2. LnP'(K), |LnP''(K)| and delta K only from true neighbours K − 1 and
   K + 1, `NA` otherwise, with a warning (F2). **Consequence: for a
   `k.range` with gaps, `ln.pk`, `ln.ppk` and `delta.k` become `NA` where
   they were computed across the gap.**
3. Warning when delta K cannot be computed for lack of replicates (F3).
4. Clear error messages; class and K-count checks in `gl.evanno` (F4).
5. `plot.theme`, `plot.dir`, `plot.file`, `verbose` arguments; flags;
   plots in `theme_dartR()`; combined plot saved with `utils.plot.save()`
   (F5). **Consequence: the Evanno plot looks different (dartR theme).**
6. Combined figure built with patchwork and returned as
   `plots$combined`; `gridExtra` dropped here (F6).
7. `aes()` with `.data` instead of `aes_string()` (F7).
8. Documentation rewrite for both functions, Rd regenerated, NEWS entry
   (F8, F9).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run on both
  functions.
- Spec: statistics checked against a hand calculation; edge cases (gaps in
  K, one replicate, identical replicates, unequal replicates, K >= 10,
  fewer than three K) — run.
- Real STRUCTURE output: one run (K = 1:4, five replicates) — run; it
  produced the `Inf` case in F1.
- Plot appearance: SKIPPED as a snapshot (no vdiffr); plots checked by
  class and panel names.
- DAT1–DAT6: not applicable (no genlight).
- dartR Google Group / GitHub issues search: not run (no search access).
- Downstream callers: `gl.run.structure` calls
  `utils.structure.evanno(sr, plot = FALSE)` and uses `plots$mean.ln.k`,
  `plots$ln.pk`, `plots$ln.ppk`, `plots$delta.k` — names kept by every
  proposed change; dartr2shiny (`template_report.csv`, `variables_matrix.csv`)
  calls `gl.evanno` — new arguments go after `plot.out`. No callers in other
  `dartR.*` packages.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | |
| 2 | approved | Luis | |
| 3 | approved | Luis | |
| 4 | approved | Luis | |
| 5 | approved | Luis | |
| 6 | approved | Luis | |
| 7 | approved | Luis | |
| 8 | approved | Luis | |

## Outcome

- Changes 1–8 applied in commit 0462be5 on `review-gl.evanno`, PR #96 to `dev`.
- Characterization test: 38 expectations pass; every diff from baseline is
  tagged `[approved n]` (1–7). Old vs new on consecutive K (1:4, 8:11, one
  replicate, random 1:8 with 4 replicates): identical `df` and panel data.
- Real STRUCTURE run: `delta.k` Inf at K = 2, 3 → NA with a warning.
- `test-gl.run.structure.R` with the STRUCTURE binary: 51 expectations pass.
- `R CMD check`: unrelated failure in `test-gl.ld.haplotype.R`; no new NOTE.
- Addendum (requested by Luis after PR #96): `gl.run.structure` called the
  helper without `verbose`, so the new warnings never reached its users.
  Fixed in commit 601ef54, PR #97 (base `review-gl.evanno`, to be retargeted
  to `dev` after #96 merges). Test with the STRUCTURE binary: warning at
  `verbose = 1`, silent at `verbose = 0`; 53 expectations pass.

## Machine block

```json
{
  "function": "gl.evanno",
  "also_covers": ["utils.structure.evanno"],
  "package": "dartR.popgen",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "95fde36",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "VRB4", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "FS2", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DEP1", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "PLT1", "status": "approved", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 8},
    {"id": "F9", "severity": "INFO", "confidence": "medium", "rule": "principle: method transparency", "status": "approved", "change": 8},
    {"id": "F10", "severity": "INFO", "confidence": "high", "rule": "FS1", "status": "no-change", "change": null}
  ],
  "coverage_skipped": [
    "plot snapshot: no vdiffr",
    "Google Group / issues search: no access"
  ],
  "status": "pr-open",
  "pr": 96
}
```
