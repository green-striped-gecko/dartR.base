# Review: gl.gen2fbm + gl.fbm2gen (dartR.base)

The pair converts genotypes between the in-memory SNPbin list (`@gen`) and
a file-backed matrix (`@fbm`, bigstatsr FBM.code256). They are reviewed
together; each finding names the file it applies to.

## Provenance

- Family mode: io (both)
- Date: 2026-09-24
- Reviewer: Claude (Opus 5.5, claude-opus-5-5, Claude Code),
  dartr-function-review v2.0.0
- Package commit: f9f1be8 (origin/dev), worktree branch `review-fbm`
- Datasets: testset.gl, testset.gs, a plain adegenet genlight built from
  testset.gl, and a simulated 2000 x 20000 dartR object with 1e5 missing
  genotypes
- Baseline: `tests/testthat/test-gl.gen2fbm.R` (20 expectations), new;
  passes at the reviewed state
- Release state: both are on `dev` only (FBM support is not on
  `main`/CRAN).
- Callers: `gl.read.dart`, `gl.read.csv`, `gl.read.PLINK` and
  `gl.read.fasta` call `gl.gen2fbm` when `fbm = TRUE`; `gl.save` and
  `gl.fst.pop` call `gl.fbm2gen`; about 115 help-page examples call
  `gl.gen2fbm` when `options(dartR_fbm = TRUE)`.

## Verdicts

**Standards: Needs work** — the error and message idioms and the roxygen
structure do not follow the conventions, and there is dead code.
**Spec: Needs work** — the round trip is exact (genotypes, NA, ploidy,
names, pop, position and all of `@other` identical). However,
`gl.fbm2gen` loads the whole matrix into memory while its documentation
promises block-wise streaming, and `gl.gen2fbm` rejects the plain genlight
it says it accepts.

## Findings

**F1 [MEDIUM, confidence: high] — gl.fbm2gen loads the full matrix despite promising blocks (DAT6 proposed rule, DOC5 proposed rule)**
`R/gl.fbm2gen.r:1-13,56` — the documentation says the conversion is
"column-chunked via bigstatsr::big_apply" with a `chunk` argument, but
the code runs `new("genlight", gen = x@fbm[])`. That decodes the entire
FBM into a dense double matrix. `chunk` is never used.
Failure scenario: 2000 individuals x 20000 loci: peak R heap 822 Mb
(baseline 260 Mb) to build a 13 Mb SNPbin list, about twice the size of
the 305 Mb dense matrix. The FBM format exists for data that do not fit
in memory. `gl.save()` calls `gl.fbm2gen()` before saving, so saving a
large FBM object can fail for lack of memory.
Proposed change: convert in blocks of `chunk` individuals, decoding
`x@fbm[rows, ]` for each block and concatenating the SNPbin lists, so the
full matrix is never decoded at once. The output is identical.
**Consequence: `chunk` becomes the number of individuals per block (it
was documented as loci per block but was ignored).**

**F2 [MEDIUM, confidence: high] — plain genlight input fails (DAT5, DOC5 proposed rule)**
`R/gl.gen2fbm.r:10-11,67-100`, `R/gl.fbm2gen.r:36` — the documentation
says a `genlight` is accepted and coerced to dartR, but there is no
coercion: `x@fbm <- G` fails on a genlight. `gl.fbm2gen` uses
`stopifnot(inherits(x, "dartR"))`.
Failure scenario: `gl.gen2fbm(g)` stops with "'fbm' is not a slot in
class 'genlight'". `gl.fbm2gen(g)` stops with 'inherits(x, "dartR") is
not TRUE' instead of returning the input unchanged, as it does for a
gen-backed dartR.
Proposed change: coerce with `methods::as(x, "dartR")`. The
`class(x) <- "dartR"` idiom used elsewhere produces an object without the
`fbm` slot that fails `validObject()` (see Notes). In `gl.fbm2gen`, return
a genlight without FBM unchanged.

**F3 [LOW, confidence: high] — error and message idioms (VRB2, VRB3, FS9)**
`R/gl.gen2fbm.r:49-51,62-65,71-76,104-106`: `stop("...")` without
`error()`; `message()` rather than `cat(report())` at `verbose > 2`. The
"already FBM-backed" path returns before "Completed:", and the
both-slots-populated case stops with "Object already FBM-backed; cannot
convert", which names the wrong problem. There is commented-out dead code
(lines 53-60) and a bare `x` (line 103). `R/gl.fbm2gen.r:47-49` prints its
"no FBM" note with `message()`.
Failure scenario: errors are not coloured like the rest of dartR; an
object with both slots set (invalid) gets a misleading message.
Proposed change: `stop(error(...))` with accurate text, `cat(report())`
for notes, "Completed:" on every path, and remove the dead code.

**F4 [LOW, confidence: high] — roxygen structure and examples (DOC1, DOC3)**
Neither file has `@name`, `@title` or `@family`, the `gl.gen2fbm` verbose
`@param` has no default, and both examples are `\dontrun{}` on undefined
objects (`gl`, `d_fbm`), so they never run.
Failure scenario: the help pages sit under no family, and nobody sees an
example break.
Proposed change: add the tags (`@family data manipulation`), the standard
verbose text, and runnable examples on testset.gl.

Notes (outside these files, not proposed here):
- `class(x) <- "dartR"` on a plain genlight gives an object with no `fbm`
  slot that fails `validObject()`; `methods::as(x, "dartR")` gives a valid
  one. The idiom appears in several functions (for example
  `gl.report.contamination`, `gl.subsample.ind`). It works today only
  because those functions do not touch `@fbm`. It is worth a
  package-wide fix.
- `as.matrix()` returns integer for gen-backed and double for FBM-backed
  objects (values identical). This is an `as.matrix` method issue, not one
  of these functions.

## Proposed changes

1. `gl.fbm2gen` converts in blocks of `chunk` individuals (F1).
   **Consequence: `chunk` means individuals per block; output identical;
   peak memory reduced.**
2. Plain genlight accepted by `gl.gen2fbm` (coerced to dartR) and returned
   unchanged by `gl.fbm2gen` (F2).
3. Error and message idioms, accurate messages, dead code removed (F3).
4. Roxygen tags, verbose text, runnable examples (F4). Docs only.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run on both
- Round trip on testset.gl: genotypes, NA count (8608), ploidy, names,
  pop, position, chromosome, n.loc and `@other` — run, all identical
- `chunk` variation in `gl.gen2fbm` (7 loci per block) — run, identical
- Memory and time of `gl.fbm2gen` at 2000 x 20000 — run (F1)
- SilicoDArT: `gl.gen2fbm` stops ("Only SNP data supported at this
  time"), as documented — run
- DEP1: bigstatsr/bigsnpr are Imports — no guard needed
- `gl.save` round trip with FBM input: not run separately; covered by
  change 1's identical output
- Google Group / GitHub issues: not searched (not available: no browser
  session); FBM support is unreleased

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence (chunk = individuals per block) approved explicitly |
| 2 | approved | Luis | |
| 3 | approved | Luis | |
| 4 | approved | Luis | |

## Outcome

- Change 1 (F1): row-blocked decoding. Output identical for chunk 1, 7,
  249, 250 and 1000 on testset.gl (test), and the SNPbin list is
  identical at 2000 x 20000. Peak R heap above baseline (separate process
  per run): 2000 x 20000, whole-matrix load +429 Mb vs chunk 256 +191 Mb;
  4000 x 20000, +884 Mb vs +334 Mb. Time is unchanged (2.4-4.4 s). The
  peak is reduced but not constant in n: the adegenet genlight constructor
  uses about 8x the block's size, and garbage accumulates between
  collections. A `gc()` after every block lowers it further (+153 / +276
  Mb) at about 50 % more run time, so it was not adopted. The default
  `chunk` is 256 individuals (it was 2048 and never used); at 256 the peak
  sits near the floor for every block size tried (16-256).
- Change 2 (F2): `methods::as(x, "dartR")`; the baseline expectations
  "not a slot" / "dartR" were flipped (approved diff). A plain genlight
  now converts and round-trips.
- Change 3 (F3): `stop(error())`; the both-slots case gives "Invalid
  object: both @fbm and @gen hold genotypes" (test); "Completed:" on every
  path; dead code removed.
- Change 4 (F4): `@name`/`@title`/`@family data manipulation`/`@seealso`,
  standard verbose text, runnable examples (both run).
- Callers' tests pass: gl.save 7, gl.fst.pop 141, gl.read.csv 37,
  gl.read.dart 40, gl.read.fasta 31, gl.read.PLINK 26, utils.read.dart 21,
  utils.read.fasta 32 expectations. The new file has 37 expectations.
- `devtools::document()` also rewrote unrelated Rd files (cross-link
  drift, PR #421); reverted.
- PR: pending.

```json
{
  "function": "gl.gen2fbm",
  "companion": "gl.fbm2gen",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "2.0.0",
  "commit": "f9f1be8",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "DAT6", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DAT5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 4}
  ],
  "coverage_skipped": ["gl.save FBM round trip not run separately", "forum/issues: no browser session"],
  "status": "awaiting-approval",
  "pr": null
}
```
