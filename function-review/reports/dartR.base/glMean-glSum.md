# Review: glMean and glSum (dartR.base)

- Family mode: analysis (utility overrides of the adegenet functions of the same name)
- Reviewed as a pair: both live in `R/utils.dartR.class.def.r` (lines 1071-1170) and `glMean` calls `glSum`
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 2bd61c5 (`R/utils.dartR.class.def.r` identical to origin/dev)
- Datasets: testset.gl and subsets; FBM copies via `gl.gen2fbm`; simulated 1000 × 20000 SNP matrix (5% missing) for timing
- Baseline: `tests/testthat/test-glMean-glSum.R` (34 expectations, captured pre-review, all pass)
- Author: Bernd Gruber (FBM support, 2025-10 to 2026-02); `@author` tag absent

## Verdict

**Standards: Needs work** — the roxygen headers are incomplete (no author/custodian, no examples, one wrong parameter description) and the FBM path makes two row-by-row passes over a column-stored file.
**Spec: Ready** — on dense, FBM-backed, subset and all-NA inputs, both functions return the same values, types and names as adegenet and as an independent computation.

## Findings

**F1 [MEDIUM, confidence: high] — FBM path reads the file one individual at a time, twice (DAT6)**
`R/utils.dartR.class.def.r:1094-1109, 1160-1167` — `glSum` loops over individuals reading `x@fbm[e, ]`, a row of a file stored column by column; `glMean` then calls the dartR `glNA` method, which scans the missing-data positions again. On 1000 individuals × 20000 loci: `glSum` 0.28 s + `glNA` 0.22 s. `bigstatsr::big_counts()` reads each column once and returns the count of each genotype code per locus, from which both the allele sum and the non-missing count follow; it took 0.037 s and gave identical sums and missing counts.
Failure scenario: FBM exists for datasets too large for memory; on such data (for example 10000 × 1 000 000), every `glMean` call — `gl.impute` calls it per population — costs minutes rather than seconds.
Proposed change: compute both functions' FBM results from one `big_counts()` pass (sum = count of 1 + 2 × count of 2; non-missing = total − count of NA), keeping the current return types (integer sums for `alleleAsUnit = TRUE`) and locus names. The dense path is unchanged.

**F2 [LOW, confidence: high] — documentation (DOC1, DOC5, DOC7)**
Lines 1071-1080, 1143-1150:
- `glSum` `@param alleleAsUnit` says "the mean is calculated"; it is the sum, counted in alleles (`TRUE`) or in individuals, each genotype divided by its ploidy (`FALSE`).
- `@description` for both reads as an internal note ("we need one for fbm projects") rather than telling the user what the function returns and when it differs from adegenet.
- `@return` omits that loci with no data give `NaN` (`glMean`) and that `glSum` returns integers for `alleleAsUnit = TRUE`.
- No `@examples`, no `@author` with custodian.
Proposed change: rewrite both headers; add an example on `testset.gl` (dense and FBM-backed).

## Proposed changes

1. Compute FBM-backed results from a single `bigstatsr::big_counts()` pass in both functions; dense path unchanged; same values, types and names (F1).
2. Rewrite both roxygen headers with examples and author/custodian (F2). Docs only. Custodian to be confirmed (author: Bernd Gruber).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run. FS/VRB structure (verbosity, flags, history) does not apply: these are low-level accessors called inside other functions, without `verbose`.
- Spec: values vs adegenet (dense) and vs an independent allele-frequency computation — run, identical
- FBM path (DAT6): full object, `[1:20, 1:50]`, reversed individual order, 5 individuals with 15 all-NA loci — run, identical to dense (values, types, names, NaN count)
- Plain `genlight` input (DAT5) — run, identical to adegenet
- Input object not modified by the `class(x) <- "genlight"` relabel (it acts on a local copy) — run
- Callers: `gl.impute` (3 calls); sibling packages import `dartR.base` and so resolve to these versions; dartR.popgenomics imports adegenet then dartR.base, so dartR.base's versions win there too
- Known complaints: none on GitHub; Google Group not searched (no access)

## Report notes (not findings here)

- `adegenet::glPca()` on an FBM-backed object fails with "subscript out of bounds": adegenet's internals read the empty `@gen` slot and call adegenet's own `glMean`, not this one. This is a general FBM limitation with adegenet functions, outside these two functions.
- The dartR `glNA` method (between the two functions) has the same row-wise pattern; change 1 removes `glMean`'s call to it but does not change `glNA` itself.
- `class(x) <- "genlight"` relabels the object rather than converting it (`as(x, "genlight")`); it works because adegenet only reads `@gen` and `@ploidy`, so no change is proposed.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | |
| 2 | approved | Luis | |

## Outcome

- Changes 1 and 2 applied on branch `review-glMean-glSum` (from `origin/dev`).
- Characterization test: 34 expectations pass with no change to any assertion (no behavioural diff).
- Old vs new on FBM-backed testset.gl, `[1:5, ]` and `[20:1, 1:50]`, both `alleleAsUnit` values: `identical()` TRUE for values, types and names.
- Timing, 1000 × 20000 FBM, mean of 3 runs: `glMean` 0.563 s to 0.037 s; `glSum` 0.311 s to 0.036 s.
- Mixed-ploidy FBM objects keep the previous row-wise code (`gl.gen2fbm` accepts SNP data only, so this path is not reached today).
- Related tests pass: `gl.impute` 75, `gl.gen2fbm` 37, `gl.alf` 68, `gl.allele.freq` 34, `utils.dartR.class.def` 17, `utils.impute` 9.
- `devtools::document()` run; examples run. Custodian set to Bernd Gruber (author) — assumption, to be confirmed.
- Full `R CMD check` not run.
- PR: pending.

## Machine block

```json
{
  "function": "glMean, glSum",
  "package": "dartR.base",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "2bd61c5",
  "verdict_standards": "needs_work",
  "verdict_spec": "ready",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "DAT6", "status": "approved", "change": 1},
    {"id": "F2", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "approved", "change": 2}
  ],
  "coverage_skipped": ["Google Group: no access"],
  "status": "pr-open",
  "pr": null
}
```
