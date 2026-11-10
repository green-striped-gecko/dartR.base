# Review: gl.randomize.snps (dartR.base)

- Family mode: analysis (behaves as `modify`: returns a recoded genlight; `modify` checks applied)
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 2bd61c5 (dev_luis, level with origin/dev)
- Datasets: testset.gl (after `gl.filter.monomorphs`, 250 ind × 111 loci), testset.gs, FBM copy of testset.gl via `gl.gen2fbm`
- Baseline: `tests/testthat/test-gl.randomize.snps.R` (19 expectations, captured pre-review, all pass)

## Verdict

**Standards: Needs work** — structure follows the house order and flags and history are handled; a dead duplicate source file, unused plot arguments and a full densification remain.
**Spec: Rework** — the returned object is wrong in three cases: SilicoDArT data gains invalid `2` codes, FBM-backed input comes back unrandomized, and the allele labels no longer describe the genotypes.

## Findings

**F1 [HIGH, confidence: high] — SilicoDArT data corrupted (DAT1)**
`R/gl.randomize.snps.r:74` — `utils.check.datatype()` accepts both datatypes, and the recoding turns every `0` into `2` in half the loci regardless. On `testset.gs`, 17,389 of 52,801 scored cells become `2` in a ploidy-1 object.
Failure scenario: a user runs the function on presence/absence data; downstream functions receive values outside 0/1 and either error or report nonsense.
Proposed change: restrict to SNP data with `utils.check.datatype(x, accept = "SNP", verbose = verbose)`, which stops with a clear message for SilicoDArT.

**F2 [HIGH, confidence: high] — FBM-backed input silently not randomized (DAT6)**
`R/gl.randomize.snps.r:92-94` — the recoded genotypes are written to `x@gen`, but for FBM-backed objects genotype reads come from `x@fbm`, which is untouched. On the FBM copy of testset.gl, `as.matrix()`, `glMean()` and `gl.fbm2gen()` on the result all equal the original. The object also now carries a full dense `@gen` (250 `SNPbin`) that disagrees with its FBM.
Failure scenario: with `options(dartR_fbm = TRUE)`, the function completes without warning, adds a history entry and resets the flags, but nothing is recoded.
Proposed change: for FBM-backed input, copy the FBM first (`bigstatsr::big_copy`, as `gl.sample` does), write the recoded columns into the copy, and leave `@gen` alone. The copy is needed because FBM storage is shared by reference: writing into `x@fbm` directly would also recode the caller's original object.

**F3 [HIGH, confidence: high] — allele labels not swapped (DAT2)**
`R/gl.randomize.snps.r:84-94` — the function swaps which homozygote is coded `0` and which `2`, but leaves `loc.all` unchanged. At locus `100049698-16-G/A`, 36 individuals were coded `0` (homozygous G) before and 1 after, while `loc.all` still reads `G/A`. The genotypes now claim the individuals are homozygous for the other allele.
Failure scenario: a user randomizes, then exports with `gl2vcf`/`gl2plink` or joins with another dataset by allele; half the loci carry the wrong allele for every homozygote. Allele frequencies reported per allele are inverted for those loci.
Proposed change: reverse `loc.all` for the recoded loci (`"G/A"` becomes `"A/G"`), so the object remains a correct description of the same individuals, with only the counted allele switched. The DArT report columns in `loc.metrics` (`SNP`, `AlleleID`, `TrimmedSequence`) keep the original reference; their flags are already reset.
**Consequence: `loc.all` changes for the recoded half of the loci.**

**F4 [MEDIUM, confidence: high] — duplicate definition in `R/gl.random.snp.r` (FS1)**
`R/gl.random.snp.r` holds an older copy of the whole function, with its own roxygen block under `@name gl.randomize.snps`. Files load alphabetically, so `gl.randomize.snps.r` overwrites it at build time; roxygen merges both blocks into `man/gl.randomize.snps.Rd`, which now runs two example sets (lines 80-85).
Failure scenario: a developer edits `gl.random.snp.r` (the file that matches a search for "random snp") and the change has no effect.
Proposed change: delete `R/gl.random.snp.r` and re-document.

**F5 [MEDIUM, confidence: high] — two full dense copies of the genotype matrix (DAT6)**
`R/gl.randomize.snps.r:80-90` — `as.matrix(x)` is called twice and both copies are kept alive, plus `hold <- x`. For large objects the function needs about three times the genotype matrix in memory.
Failure scenario: a large dataset exhausts memory in a function whose job is a column relabel.
Proposed change: build one matrix and recode only the sampled columns with `m[, idx] <- 2 - m[, idx]` (maps 0↔2, keeps 1 and `NA`). Keep `hold` only when `plot.display = TRUE`. Output is identical for the same seed.

**F6 [LOW, confidence: high] — `plot.theme` and `plot.colors` ignored (PLT1)**
`R/gl.randomize.snps.r:45-65, 124, 133` — both arguments are processed but never passed to `gl.smearplot()`. The warning says "More than 2 colors specified" while the code keeps 4, and the documented default (`c("#2171B5","#6BAED6")`) differs from the code default (4 colours).
Failure scenario: a user sets `plot.colors` or `plot.theme` and the plot does not change.
Proposed change: pass `plot.theme` and `plot.colors` to both `gl.smearplot()` calls; fix the warning text and the `@param` defaults.

**F7 [LOW, confidence: high] — documentation disagrees with behaviour (DOC5, proposed rule)**
Roxygen, lines 2-33 — `@details` says plots are saved to the session's temporary directory (they are saved only when `plot.file` is set); `@description` says 0s are recoded to 2s (the swap is both ways); `@return` omits the reset metric flags and, after F3, the swapped `loc.all`; `plot.dir` mentions RDS files only.
Failure scenario: a user looks for saved plots in `tempdir()` and finds none.
Proposed change: rewrite the header to match behaviour. Docs only.

## Proposed changes

1. Restrict input to SNP data; SilicoDArT stops with a message (F1). **Consequence: calls on SilicoDArT objects that ran before now error.**
2. Write the recoding into a copy of the FBM for FBM-backed input (F2). **Consequence: FBM-backed input is now actually recoded.**
3. Reverse `loc.all` for recoded loci (F3). **Consequence: `loc.all` changes for the recoded half of the loci.**
4. Delete the duplicate `R/gl.random.snp.r` and re-document (F4).
5. Single matrix, recode sampled columns only; keep the pre-plot copy only when plotting (F5).
6. Pass `plot.theme`/`plot.colors` to the smear plots; fix warning text and defaults (F6).
7. Rewrite the roxygen header to match behaviour (F7). Docs only.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: behaviour vs roxygen on testset.gl, testset.gs — run
- FBM path (DAT6): run on `gl.gen2fbm(testset.gl)`
- Plain `genlight` input (DAT5): run; returns a `genlight`, no error
- Plot/result independence (PLT3): run; same seed gives identical genotypes with and without plotting
- Edge case, 1 locus: run; `floor(1/2) = 0` loci recoded, returns unchanged object silently (not raised as a finding)
- Known complaints: GitHub issues for dartR.base and dartR searched for "randomize" — none found. Google Group: not searched (no access from this session)
- Callers: none in sibling `dartR.*` packages; exposed in dartr2shiny (`shiny_fun/Fun_gl.randomize.snps.R`, config `Data_Refinement`)

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

## Outcome

- Changes 1-7 applied on branch `review-gl.randomize.snps` (from `origin/dev`).
- Characterization test: 26 expectations pass. Diffs from baseline map to approved changes only: SilicoDArT now errors (1), FBM result equals the dense result and the input FBM is untouched (2), `loc.all` reversed at the 55 recoded loci only (3), plot arguments reach `gl.smearplot` (6).
- SNP genotypes: the pre-review function (from `git show HEAD:R/gl.randomize.snps.r`) and the new one give identical `as.matrix()` output for `set.seed(7)`; `loc.metrics` and `ind.metrics` identical; `loc.all` differs at 55 of 111 loci (change 3).
- `verbose = 3` end to end on FBM-backed testset.gl: `validObject()` TRUE, FBM kept, 55 loci recoded. Examples run with `dartR_fbm` FALSE and TRUE.
- `devtools::document()` run; `man/gl.randomize.snps.Rd` now carries one example set.
- NEWS entry added. Callers: none in sibling `dartR.*` packages; dartr2shiny exposes the function (SilicoDArT input there will now error).
- PR: pending.

## Machine block

```json
{
  "function": "gl.randomize.snps",
  "package": "dartR.base",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "2bd61c5",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT1", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DAT6", "status": "approved", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "DAT2", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS1", "status": "approved", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DAT6", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "PLT1", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 7}
  ],
  "coverage_skipped": ["Google Group search: no access from session"],
  "status": "pr-open",
  "pr": null
}
```
