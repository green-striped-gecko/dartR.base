# Review: gl.report.polyploid_heterozygosity (dartR.base)

- Family mode: report
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 2bd61c5 (dev_luis, level with origin/dev)
- Datasets: testset.gl; simulated autotetraploid dosage data (40 individuals in 2 populations, 200 loci, allele frequencies U(0.05, 0.95), binomial dosages 0-4, 200 missing genotypes, `set.seed(42)`); the same simulation at ploidy 2; FBM copy of testset.gl
- Baseline: `tests/testthat/test-gl.report.polyploid_heterozygosity.R` (25 expectations, captured pre-review, all pass)
- Custodian: Ching Ching Lau (author, 2024-08). The diploid sibling review was held for author consultation (STY5); the same may apply here.

## Verdict

**Standards: Needs work** — house structure is in place, but the crash fixes, verbose gating and return shape applied to `gl.report.heterozygosity` in its review were never carried to this copy.
**Spec: Rework** — on polyploid data, every reported statistic is wrong: He is computed from diploid genotype counts, Ho counts heterozygous individuals instead of the documented gametic heterozygosity, and the individual method counts only dosage 1. Bootstrap intervals are wrong at every ploidy. The function works correctly only where it reproduces the diploid function (Ho and He on diploid data).

## Findings

**F1 [BLOCKER, confidence: high] — He wrong for polyploid data (DAT1)**
`R/gl.report.polyploid_heterozygosity.r:555-563` — allele frequencies come from `utils.recalc.freqhomref/freqhomsnp/freqhets`, which count only codes 0, 1 and 2 (`colMeans(t == 0)` etc.). For a tetraploid, dosages 3 and 4 are dropped and `p + q` is not 1. On the simulated tetraploid, He is 0.723 and 0.717; the independent value, `1 - p^2 - q^2` with `q = sum(dosage) / (4 * n)`, is 0.356 and 0.358.
Failure scenario: any polyploid dataset — the function's reason to exist — reports He, uHe and FIS roughly doubled or worse.
Proposed change: see change 1.

**F2 [BLOCKER, confidence: high] — Ho is not the documented gametic heterozygosity (DOC5, spec)**
`R/gl.report.polyploid_heterozygosity.r:440-456` — `gamete_het()` computes the gametic heterozygosity of each genotype, `d(k - d) / choose(k, 2)`, and then the code thresholds it with `> 0`, so Ho becomes the share of individuals that are heterozygous at all. The documentation (Details, and the Moody et al. 1993 reference) defines Ho as the proportion of heterozygous gametes. On the simulated tetraploid, reported Ho is 0.645; gametic Ho is 0.357, which matches He (0.356) as expected under Hardy-Weinberg. FIS combines the two inconsistent quantities. For diploids both definitions give the same Ho.
Failure scenario: a user compares Ho with He or reads FIS on polyploid data; Ho looks nearly twice He, suggesting a heterozygote excess that does not exist.
Proposed change: see change 1.

**F3 [HIGH, confidence: high] — bootstrap resamples individuals, not loci (spec)**
`R/gl.report.polyploid_heterozygosity.r:638-649` — `boot::boot()` receives the individuals × loci matrix, so it resamples individuals; `pop.het()` then transposes and treats loci as individuals. `gl.report.heterozygosity` transposes before calling `boot()`. The statistic is also the diploid `pop.het_fun` (`colMeans(df == 1)`, `q = colMeans / 2`). On the simulated diploid data, He = 0.349 with a 95% interval of 0.498–0.500, and FIS = 0.028 with an interval of 0.28–0.34; on the tetraploid, the He interval is −0.103 to −0.030.
Failure scenario: every `nboots > 0` run, at any ploidy, reports intervals that do not describe the estimate.
Proposed change: transpose before `boot()` as the sibling does, and use the ploidy-aware statistic of change 1.

**F4 [HIGH, confidence: high] — `method = "ind"` counts only dosage 1 (DAT1)**
`R/gl.report.polyploid_heterozygosity.r:998-1006` — Ho per individual is `m == 1`, homozygous reference `m == 0`, homozygous alternative `m == 2`. For a tetraploid, dosages 2 and 3 are also heterozygous and dosage 4 is the alternative homozygote. On the simulation, mean individual Ho is 0.194 (share heterozygous: 0.647), and the three proportions sum to 0.58 instead of 1. The documentation defines individual Ho as the proportion of heterozygous gametes.
Failure scenario: individual heterozygosity for polyploids is underestimated about threefold; `gl.filter.heterozygosity`-style outlier calls on it are wrong.
Proposed change: see change 2.

**F5 [HIGH, confidence: high] — crashes already fixed in the diploid sibling (PLT3, spec)**
All verified on testset.gl:
- `method = "ind"` with `verbose = 0` or `plot.display = FALSE`: "object 'p1' not found" (outliers read from the built plot, line 1042).
- `plot.file` with `plot.display = FALSE`: "object 'p3' not found" (line 1087).
- `method = "ind", subsample.pop = TRUE`: "object 'res_sub' not found" (line 1078).
- `subsample.pop = TRUE` with plotting: "Supplied 150 items to be assigned to 100 items of column 'color'" — populations below `n.limit` are skipped, but colours are recycled with `rep(colors_pops, each = 5)` (line 818).
- `nboots > 0` where a replicate statistic is constant: "replacement has length zero" (lines 669-680).
Failure scenario: the default verbose-0 call on the individual method, used in pipelines, errors.
Proposed change: port the fixes from `gl.report.heterozygosity` (its F1-F3 and addendum, plus its `boot.ci` NA guard): outliers from `boxplot.stats()`, colours keyed by population name, `subsample.pop` ignored with a warning under `method = "ind"`, `plot.file` guarded, NA limits when `boot.ci` cannot form an interval.

**F6 [MEDIUM, confidence: high] — diploid uHe and FIS diverge from `gl.report.heterozygosity` (spec)**
`R/gl.report.polyploid_heterozygosity.r:571, 584` — the sample-size correction uses the population's mean sample size and is hard-coded `2n / (2n - 1)`; the sibling uses each locus's own sample size. On testset.gl, FIS differs by up to 0.144 (EmmacMDBSanf: −0.334 here, −0.190 in the sibling); uHe by up to 0.002.
Failure scenario: the same diploid data give two different FIS values depending on which function is called.
Proposed change: covered by change 1.

**F7 [MEDIUM, confidence: high] — repeated densification (DAT6)**
`gamete_het()` is called three times per population, each calling `as.matrix()` twice; `n_loc` and the monomorphic-locus loop densify again. For large datasets this multiplies memory and run time.
Proposed change: covered by change 1 (one matrix per population).

**F8 [LOW, confidence: high] — ungated warnings, empty error message (VRB3)**
Lines 354-395 — the method, negative `n.invariant` and `gl.filter.secondaries` warnings print at `verbose = 0`; `error.bar = "CI"` with `nboots = 0` prints a message and then calls `stop()` with an empty message.
Proposed change: gate warnings at `verbose >= 1`; single `stop(error(...))`, as in the sibling.

**F9 [LOW, confidence: high] — `subsample.pop` returns an unnamed list (DOC5)**
Line 1104 — `list(res_sub, df)`, while `@return` promises a data frame.
Proposed change: `list(subsample = res_sub, results = df)` and document both shapes, as in the sibling.
**Consequence: code that indexes the result by position keeps working; code expecting no names sees names.**

**F10 [LOW, confidence: high] — documentation (DOC1, DOC3, DOC5)**
- No runnable examples: the block is commented out (`# @examples`) and refers to `package = 'dartR'` and diploid `platypus.gl`.
- Title says "from SNP data" with no mention of dosage/polyploid input; `@param x` wording is ungrammatical.
- Nei 1978 and Schmidt et al. 2021 are cited in Details but missing from `@references`; "Schimdt"/"Schimddt" misspelt.
- Details equations describe diploid formulas (`2 * n_Ind`); after change 1 they need the ploidy `k`.
- `@family unmatched report`.
Proposed change: rewrite the header; add a runnable example that builds a small tetraploid genlight.

## Proposed changes

1. Ploidy-aware population statistics, shared by the point estimates and the bootstrap: per locus, `q = sum(dosage) / sum(ploidy)`, He `= 1 - p^2 - q^2`, gametic Ho `= mean(d(k - d) / choose(k, 2))`, uHe `= He * kn / (kn - 1)` with that locus's own sample size `n`, FIS `= 1 - Ho / uHe`; one genotype matrix per population. Implemented by adding a `ploidy` argument to the internal `pop.het_fun` (default 2 gives exactly its current diploid result, so `gl.report.heterozygosity` is unchanged). (F1, F2, F6, F7)
   **Consequence: for polyploid data Ho, He, uHe, FIS and their SD/SE change (simulated tetraploid: Ho 0.645 to about 0.357, He 0.723 to about 0.356); for diploid data Ho and He are unchanged and uHe/FIS change to match `gl.report.heterozygosity`.**
2. Ploidy-aware individual method: Ho per individual is the mean gametic heterozygosity; `f.hom.ref` is dosage 0 and `f.hom.alt` is dosage `k`. (F4)
   **Consequence: individual Ho, `f.hom.ref` and `f.hom.alt` change for polyploid data; diploid output unchanged.**
3. Bootstrap loci, not individuals: transpose before `boot()`, use the change-1 statistic. (F3)
   **Consequence: all bootstrap intervals change, at every ploidy.**
4. Port the sibling's crash fixes (outliers from data, colour keying, `subsample.pop` under `method = "ind"`, `plot.file` guard, `boot.ci` NA guard). (F5)
5. Gate warnings at `verbose >= 1`; single `stop(error(...))`. (F8)
6. Named list return for `subsample.pop = TRUE`. (F9) **Consequence: the returned list gains names.**
7. Rewrite the roxygen header with a runnable tetraploid example. (F10) Docs only.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: behaviour vs roxygen on testset.gl and simulated tetraploid — run; independent re-computation of He, gametic Ho and uHe for the tetraploid
- Diploid cross-check against `gl.report.heterozygosity` — run (Ho, He identical; uHe, FIS diverge, F6)
- Report-family checks: input object returned untouched — run (identical); no history append — run
- FBM path (DAT6): run on `gl.gen2fbm(testset.gl)`; Ho identical to dense
- Real polyploid data: SKIPPED — no packaged dosage-mode dataset; the simulation stands in
- Mixed-ploidy objects: not run; `utils.check.datatype` warns "ploidy is not uniformly 2; treating as SNP data"; change 1 handles per-individual ploidy by construction
- Plain `genlight` not passed through `gl.compliance.check`: fails inside `gl.allele.freq` ("'names' attribute [8] must be the same length as the vector [7]") — a `gl.allele.freq` issue, noted, not a finding here
- Known complaints: none found on GitHub; Google Group not searched (no access from this session)
- Callers: none in sibling `dartR.*` packages or dartr2shiny

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | |
| 2 | approved | Luis | |
| 3 | approved | Luis | |
| 4 | approved | Luis | |
| 5 | approved | Luis | |
| 6 | approved | Luis | |
| 7 | approved | Luis | `@family` kept as `unmatched report`: `gl.filter.heterozygosity` is diploid-only |
| A1 | approved | Luis | |
| A2 | approved | Luis | keep corrected counts; sibling to be fixed separately |

## Addendum (discovered during apply)

**A1 [HIGH, confidence: high] — `subsample.pop` Ho counts only dosage 1 (DAT1)**
`R/utils.het.report.r` `het_rep()` — the subsampled Ho is `colMeans(mat == 1)`, the same defect as F4, so polyploid subsampling results are wrong. The helper is shared with `utils.subsample.pop` (pending review) and `gl.report.heterozygosity`.
Proposed change: pass per-individual ploidy to `het_rep()` and use gametic heterozygosity when ploidy is not 2; the diploid path stays unchanged.
**Consequence: polyploid `subsample.pop` Ho values change; diploid unchanged.**

**A2 [MEDIUM, confidence: high] — all-NA loci counted as polymorphic (spec)**
Original lines 485-505 — loci with no data in a population are counted in `polyLoc` and subtracted from `monoLoc`. testset.gl population EmmacBrisWive: reported 21 polymorphic / 224 monomorphic / 10 all-NA; the independent count is 11 / 234 / 10. The diploid sibling `gl.report.heterozygosity` has the same defect (not fixed here). The new one-matrix count (change 1) gives 11 / 234 / 10.
Proposed change: keep the corrected count here and raise the sibling separately.
**Consequence: `polyLoc` and `monoLoc` change wherever a population has all-NA loci; the polyploid function then differs from the diploid function in these two columns until the sibling is fixed.**

## Outcome

- Changes 1-7 and addenda A1, A2 applied on branch `review-gl.report.polyploid_heterozygosity` (from `origin/dev`). The function body is rebuilt on the reviewed `gl.report.heterozygosity` code; `pop.het_fun`, `pop.het`, `het_rep` and `utils.subsample.pop` in `R/utils.het.report.r` gain a `ploidy` argument whose default (2) runs the unchanged diploid code.
- Characterization test: 37 expectations pass. Every diff from baseline maps to an approved change: tetraploid Ho/He/uHe/FIS equal an independent computation (1); tetraploid individual Ho/f.hom.ref/f.hom.alt equal an independent computation (2); bootstrap intervals contain their estimates at ploidy 2 and 4 (3); five former crashes run (4); silence at `verbose = 0` and error message (5); named list (6); polyLoc/monoLoc corrected (A2); tetraploid subsample Ho equals gametic Ho (A1).
- Diploid: on testset.gl every column except `polyLoc`/`monoLoc` equals `gl.report.heterozygosity()` (population and individual methods).
- Shared-helper safety: `gl.report.heterozygosity` (32), `gl.test.heterozygosity` (38), `gl.filter.heterozygosity` (19), `gl.alf` (68) and subsample tests pass unchanged; diploid `utils.subsample.pop` output identical to the pre-review helper for the same seed.
- `verbose = 3` end to end on the example tetraploid with `nboots = 200`, `error.bar = "CI"`, `subsample.pop = TRUE`: runs; returns the named list.
- `devtools::document()` run; NEWS entry added. Callers: none in sibling packages or dartr2shiny.
- Full `R CMD check` not run.
- PR: pending.

## Machine block

```json
{
  "function": "gl.report.polyploid_heterozygosity",
  "package": "dartR.base",
  "family": "report",
  "skill_version": "2.0.0",
  "commit": "2bd61c5",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "BLOCKER", "confidence": "high", "rule": "DAT1", "status": "approved", "change": 1},
    {"id": "F2", "severity": "BLOCKER", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "spec", "status": "approved", "change": 3},
    {"id": "F4", "severity": "HIGH", "confidence": "high", "rule": "DAT1", "status": "approved", "change": 2},
    {"id": "F5", "severity": "HIGH", "confidence": "high", "rule": "PLT3", "status": "approved", "change": 4},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "spec", "status": "approved", "change": 1},
    {"id": "F7", "severity": "MEDIUM", "confidence": "high", "rule": "DAT6", "status": "approved", "change": 1},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "VRB3", "status": "approved", "change": 5},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 6},
    {"id": "F10", "severity": "LOW", "confidence": "high", "rule": "DOC3", "status": "approved", "change": 7},
    {"id": "A1", "severity": "HIGH", "confidence": "high", "rule": "DAT1", "status": "approved", "change": "A1"},
    {"id": "A2", "severity": "MEDIUM", "confidence": "high", "rule": "spec", "status": "approved", "change": "A2"}
  ],
  "coverage_skipped": ["real polyploid dataset: none packaged", "mixed ploidy: not run", "Google Group: no access"],
  "status": "pr-open",
  "pr": null
}
```
