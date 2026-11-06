# Review: gl.outflank (dartR.popgen)

Covers `gl.outflank` and the helpers it runs: `utils.outflank`,
`utils.outflank.MakeDiploidFSTMat`, `utils.outflank.plotter` and the
internal functions in `utils.outflank.diploids.r`, `utils.outflank.fst.r`
and `utils.outflank.likelihood.r`.

- Family mode: analysis
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 4b177c9 (origin/dev, reviewed state)
- Datasets: bandicoot.gl (dartR.data, the roxygen example: 96 x 1,000, 5
  populations); Balding-Nichols simulation (5 populations x 30 diploids,
  2,000 neutral loci at F = 0.05 and 20 selected loci at F = 0.4, 5%
  missing calls); testset.gs. Reference: the authors' `OutFLANK` package
  (`MakeDiploidFSTMat` + `OutFLANK`) run on the same genotypes.
- Baseline: tests/testthat/test-gl.outflank.R (9 tests on a smaller
  simulation, snapshot captured pre-review; defects marked
  `BASELINE (F<n>)`)

**Standards: Needs work** — the wrapper has no verbosity, no start/end
flags and no input checks, and its dependency guard returns `-1` instead
of stopping.
**Spec: Needs work** — the Fst calculation and the likelihood fit are
faithful copies of OutFLANK, but the wrapper feeds every SNP in twice.
That doubles the outlier count it reports, shifts the inferred df and
q-values away from OutFLANK, and makes locus names with a dot fail.

What works well: per-locus Fst, He and mean Fst agree with the `OutFLANK`
package to machine precision (max difference 2.2e-16), and when the same
helpers get one column per SNP they reproduce `OutFLANK` exactly (q-value
difference 0).

## Findings

**F1 [HIGH, confidence: high] — every SNP is analysed twice, once per
allele (principle: agreement with the reference implementation; DOC5)**
`R/gl.outflank.r:61-87` — the genlight is converted with `gl2gi()`, and
`as.matrix()` of the genind has one column per allele, two per SNP. Both
columns go to `utils.outflank.MakeDiploidFSTMat()` as separate loci. The
trimming, the ML fit of df and the q-values then run on 2L entries, and
the duplicates are removed only afterwards by stripping everything after
the first dot of the name.
Failure scenario (simulation, 20 selected loci):

| | gl.outflank | OutFLANK package |
|---|---|---|
| `numberHighFstOutliers` | 26 | 13 |
| loci flagged | 13 | 13 |
| `dfInferred` | 4.808 | 4.796 |
| max q-value difference | 0.0037 | — |

The count is doubled because it is taken before the duplicates are
removed. Here the flagged set is the same, but a locus near `qthreshold`
can move either way. `results$meanAlleleFreq` refers to whichever allele
genind lists first, which changes from locus to locus (correlation with the
alternate-allele frequency 0.35). The Fst step also takes 2.5 times as
long (10,100 loci: 7.5 s against 2.9 s).
Proposed change: pass the 0/1/2 genotype matrix of the genlight
(`as.matrix(x)`, NA as 9) to the helpers, as OutFLANK expects; convert
genind input with `gi2gl()`; drop the name stripping and de-duplication
(change 1).

**F2 [HIGH, confidence: high] — loci missing in every individual are
removed silently, so `index` no longer lines up with the input (DAT2)**
`R/gl.outflank.r:62` — `gl.filter.allna()` drops all-NA loci before the
analysis. `index` and `results` then have fewer rows than `nLoc(x)`, and
every locus after a dropped one sits one position earlier.
Failure scenario: one all-NA locus at position 5 of 2,020. `index` has
2,019 values; `sel-1-A/G` is at position 2,000 in `results` and 2,001 in
the genlight. `x[, !res$index]` then selects the wrong loci without any
warning.
Proposed change: keep every input locus in `results` and `index`; loci
without data get `NA` (as monomorphic loci already do) (change 2).

**F3 [MEDIUM, confidence: high] — `index` is TRUE for loci that are not
outliers, while `@return` calls it "an index of outliers" (DOC5)**
`R/gl.outflank.r:18, 93` — `index.outflank <- !(outf$results$OutlierFlag)`.
Loci that could not be tested (monomorphic, He below `Hmin`) are `NA`.
Failure scenario: a user following the help page runs
`x[, res$index]` expecting the outliers and gets every neutral locus
(2,007 of 2,020 in the simulation). Monomorphic loci are `NA`, and `NA`
in a logical index causes an error when subsetting.
Proposed change: keep the values (flipping them would silently invert
existing scripts) and document them exactly: `TRUE` = not an outlier,
`FALSE` = outlier, `NA` = not tested; show `x[, which(!res$index)]` in
the example (change 3, docs only).

**F4 [MEDIUM, confidence: high] — locus names that contain a dot make the
function fail (DAT5)**
`R/gl.outflank.r:63, 85` — `gl2gi()` builds genind columns named
`<locus>.<allele>`, and adegenet refuses names with more than one dot.
Failure scenario: VCF-derived names such as `scaf.1-A/G`: "more than one
'.' in column names; please name column as [LOCUS].[ALLELE]". Resolved by
change 1, which no longer goes through genind or strips names.

**F5 [MEDIUM, confidence: high] — no input checks; failures surface as
unrelated errors or as return values (FS4, FS5, DEP1)**
`R/gl.outflank.r:49-58, 69-73`; `R/utils.outflank.r:129-145, 198-225` and
`R/utils.outflank.diploids.r:88-100`.
- A missing `qvalue` prints a message and returns `-1` (DEP1 requires
  `stop(error())`).
- One population: "missing value where TRUE/FALSE needed".
- No population: `MakeDiploidFSTMat` prints "your population names do not
  match your SNP matrix" and carries on.
- SilicoDArT input fails inside `gl2gi()`.
- When `utils.outflank` cannot trim or fit, it prints and returns `NULL`,
  `0` or the string "FAIL". `gl.outflank` then fails with
  "$ operator is invalid for atomic vectors".
Failure scenario: `gl.outflank()` on a single-population object gives an
error that names no argument and no cause.
Proposed change: datatype check (SNP only) and a check for at least two
populations before any work. Stop with `error()` when `qvalue` is
missing. In `utils.outflank` and `MakeDiploidFSTMat`, replace
print-and-return with `stop(error(...))`, keeping the original advice text
(change 4).

**F6 [MEDIUM, confidence: high] — no `verbose` argument; output cannot be
silenced (FS2, FS3, FS9, VRB1, VRB2)**
`R/gl.outflank.r:42-48`; `R/utils.outflank.diploids.r:102, 108-110`;
`R/utils.outflank.r:164, 233` — there is no `verbose` argument, no
`utils.flag.start()` and no end message. "Calculating FSTs, may take a
few minutes..." and the `gl2gi` messages print on every call, and
`print(paste(i, "done of", nloci))` every 10,000 loci.
Failure scenario: a script looping over data sets cannot run quietly.
Proposed change: add `verbose = NULL` with the standard blocks; progress
at 2 and a summary at 3 (loci tested, outliers, df, mean Fst). The helpers
print only when the wrapper's verbosity allows (change 5).
**Consequence: the signature gains a `verbose` argument (additive; existing
calls keep working).**

**F7 [LOW, confidence: high] — `...` is documented as "additional
parameters" but is ignored (API2, proposed rule)**
`R/gl.outflank.r:17, 48` — nothing reads `...`.
Failure scenario: `gl.outflank(x, NumberOfSamples = 2)` runs with the
number of populations, and nothing tells the user that the argument was
ignored.
Proposed change: remove `...` so R rejects unknown arguments (change 6).
**Consequence: calls that pass extra arguments, which currently have no
effect, stop with "unused argument".**

**F8 [LOW, confidence: high] — the plot ignores the user's `Hmin` (PLT1)**
`R/gl.outflank.r:90` — `utils.outflank.plotter(outf)` uses its default
`Hmin = 0.1`, so with `Hmin = 0.2` the histogram shows loci the analysis
excluded. The plot is base graphics, not ggplot (PLT1). A ggplot rewrite
is outside the scope of this review.
Proposed change: pass `Hmin` to the plotter (change 7).

**F9 [LOW, confidence: high] — roxygen gaps (DOC1, DOC2, DOC5, DOC7
(proposed rule))**
`R/gl.outflank.r:1-40` — no `@name`, `@title`, `@family` or
`@author`/custodian. The title says "per population" but the test is one
Fst scan across all populations. `@details` says OutFLANK must be
installed from GitHub; the package bundles the code and needs only
`qvalue`. `@return` does not describe `index` (F3) or the `outflank`
elements.
Proposed change: rewrite the header in house order (change 3).

**F10 [INFO] — behaviour inherited from OutFLANK**
`R/utils.outflank.diploids.r:16-18` — a polymorphic locus with identical
frequencies in every population (`s2 == 0`) returns a single 0, which is
recycled into every column, so its He is recorded as 0 and the locus is
excluded as low-He. The `OutFLANK` package does the same; no change
proposed, so results stay comparable with it.

## Proposed changes

1. Analyse the 0/1/2 genotype matrix directly (one column per SNP, NA as
   9); convert genind input with `gi2gl()`; drop the name stripping and
   de-duplication; count populations after dropping empty levels (F1, F4).
   **Consequence: numerical output changes. `numberHighFstOutliers` and
   `numberLowFstOutliers` halve to the true counts; `dfInferred` and
   q-values change to match the OutFLANK package exactly; outlier calls
   can change for loci near `qthreshold`; `meanAlleleFreq` becomes the
   reference-allele frequency for every locus.**
2. Keep all-NA loci in `results` and `index` as `NA` so both have one row
   per input locus, in input order (F2).
   **Consequence: `index` and `results` get one row per input locus when
   all-NA loci are present (previously shorter and shifted).**
3. Roxygen rewrite: document `index` as TRUE = not an outlier, FALSE =
   outlier, NA = not tested; list the `outflank` elements; fix the title and
   the installation note; add author/custodian (F3, F9). Docs only.
4. Input checks and hard errors: SNP datatype, at least two populations,
   `stop(error())` for missing `qvalue`; the helpers stop with their
   original advice text instead of returning `NULL`/`0`/"FAIL" (F5).
5. Add `verbose` with the standard structure; the helpers' output follows
   it (F6).
   **Consequence: the signature gains `verbose = NULL`.**
6. Remove `...` (F7).
   **Consequence: calls passing extra arguments now error instead of
   ignoring them.**
7. Pass `Hmin` to the plot (F8).

Callers: dartr2shiny generator copy (`input_generator/dartR.popgen/gl.outflank.r`,
takes the roxygen header); the old `dartR` package (unmaintained) has its
own copy. No `dartR.*` sibling calls these functions.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run. FS8 not
  applicable (returns a list, does not modify `x`).
- Spec: behaviour vs roxygen on bandicoot.gl and the simulation — run.
- Numerical check against the `OutFLANK` package on the same genotypes
  — run (Fst exact; df, q-values and outlier counts differ as in F1).
- Power sanity check: 12 of 20 simulated selected loci flagged, 1 neutral
  false positive, in both implementations.
- Plot: drawn to a null device; result unchanged by `plot` — run; visual
  check of the histogram not done.
- dartR Google Group: dartrverse Gmail searched for "outflank", no
  messages. GitHub issues not searched.
- FBM path (DAT6): SKIPPED — no FBM fixture; the function densifies the
  genotype matrix (the genind route already did).
- `utils.outflank.plotter` arguments `withOutliers = FALSE` (index
  misalignment after NA removal) and `Zoom`: not reachable from
  `gl.outflank`; not reviewed further.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence approved: numerical output changes |
| 2 | approved | Luis | consequence approved: index/results one row per input locus |
| 3 | approved | Luis |  |
| 4 | approved | Luis |  |
| 5 | approved | Luis | consequence approved: signature gains verbose |
| 6 | approved | Luis | consequence approved: extra arguments now error |
| 7 | approved | Luis |  |

## Outcome

- Changes 1-7 applied (branch `review-outflank`): `R/gl.outflank.r`,
  `R/utils.outflank.r`, `R/utils.outflank.diploids.r`; three Rd files
  regenerated; NEWS entry added. `@family selection` is new (one member),
  chosen because a selection scan does not belong in "population
  structure". Custodian recorded as Bernd Gruber (author of the wrapper
  and custodian of `utils.outflank`), since the file named none.
- Snapshot diffs against the pre-review baseline: 8, all mapped. Outlier
  count no longer doubled, df and q-values now equal OutFLANK (change 1,
  three expectations); names with a dot no longer fail (change 1); index
  one value per input locus and positions unshifted (change 2, two
  expectations); one-population error now names the cause (change 4);
  progress text silent by default (change 5). The index-values test is
  unchanged (change 3 is docs only).
- Tests rewritten for the approved behaviour: 12 tests, 29 expectations,
  all pass, including exact equality with `OutFLANK::OutFLANK` for df,
  q-values, flags and outlier counts.
- Simulation (2,020 loci): df 4.7962, 13 high-Fst outliers, q-values
  identical to OutFLANK (before: df 4.8082, 26 reported). Run time for
  10,100 loci: 1.2 s (before: 7.5 s). bandicoot.gl at `verbose = 3`: 942
  of 1,000 loci tested, no outliers, df 4.258.
- Correction to F3 made during Phase C: loci below `Hmin` are not tested
  but OutFLANK gives them `OutlierFlag = FALSE`, so their index is TRUE,
  not NA as first drafted in the help text. The help page and the verbose
  summary now say so; behaviour unchanged and identical to OutFLANK.
- The `MakeDiploidFSTMat` value check now accepts any subset of 0/1/2/9
  (it required all three genotype classes to be present) and stops
  instead of printing; part of change 4.
- PR: dartR.popgen#111 (commit a4ef513, branch `review-outflank`).

## Machine block

```json
{
  "function": "gl.outflank",
  "package": "dartR.popgen",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "4b177c9",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DAT2", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DAT5", "status": "approved", "change": 1},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 4},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "FS2", "status": "approved", "change": 5},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "API2", "status": "approved", "change": 6},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "PLT1", "status": "approved", "change": 7},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 3},
    {"id": "F10", "severity": "INFO", "confidence": "high", "rule": "principle: parity with OutFLANK", "status": "no_change", "change": null}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "GitHub issues not searched", "plot visual check"],
  "status": "done",
  "pr": 111
}
```
