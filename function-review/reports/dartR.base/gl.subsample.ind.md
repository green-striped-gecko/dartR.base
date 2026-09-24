# Review: gl.subsample.ind + gl.subsample.loc + gl.subsample.loci (dartR.base)

Reviewed together as the subsampling trio. Each finding names the file it
applies to. `utils.subsample.pop` (in `utils.het.report.r`) is left for the
heterozygosity-report review.

## Provenance

- Family mode: modify (all three)
- Date: 2026-09-24
- Reviewer: Claude (Opus 5.5, claude-opus-5-5, Claude Code),
  dartr-function-review v2.0.0
- Package commit: f9f1be8 (origin/dev), worktree branch `review-subsample`
- Datasets: testset.gl (250 individuals, 255 loci, 30 populations of 1-11),
  platypus.gl (populations of 23, 17 and 41), testset.gs, a plain adegenet
  genlight built from testset.gl, and FBM conversions
- Baseline: `tests/testthat/test-gl.subsample.ind.R` (10 expectations),
  `test-gl.subsample.loc.R` (7) and `test-gl.subsample.loci.R` (4), new
  files; all pass at the reviewed state
- Release state: all three are on `main`/CRAN.
- Callers: `gl.subsample.loc` in dartR.popgen `gl.nhybrids` (with
  `replace = FALSE`) and in two dartR.captive examples. There are no
  callers of `gl.subsample.ind` or `gl.subsample.loci` in the siblings or
  in dartr2shiny's config.

## Verdicts

**Standards: Needs work** — history, flag resets and metadata sync follow
the conventions (recent fixes 3de3266 and 7887bf7); the remaining issues
are argument checks and stale metrics.
**Spec: Rework** — `gl.subsample.ind` returns the wrong number of
individuals when it samples across all individuals or upsamples by
population, and it crashes at its documented default.
`gl.subsample.loc` behaves as documented.

## Findings

**F1 [HIGH, confidence: high] — upsampling by population returns too many individuals (spec)**
`R/gl.subsample.ind.r:129` — `for (i in 1:trunc(n/nInd(x))-1)` parses as
`(1:k) - 1`, which is k iterations rather than k - 1. With the initial
draw, each population above its size receives one extra full copy of
itself.
Failure scenario: platypus.gl, `n = 70, by.pop = TRUE`: populations of
23/17/41 return 93/87/70 individuals; `n = 100` returns 123/117/141.
Every population should have n.
Proposed change: `for (i in seq_len(trunc(n / nInd(x)) - 1))`.

**F2 [HIGH, confidence: high] — sampling across all individuals is capped by the number of loci (spec)**
`R/gl.subsample.ind.r:99-109` — the `by.pop = FALSE` branch compares `n`
with `nLoc(x)` instead of `nInd(x)`, and sets `n <- nLoc(x)` on both
branches of the `replace` test.
Failure scenario: testset.gl (250 individuals, 255 loci):
`n = 300, replace = TRUE` and `n = 600, replace = TRUE` both return 255
individuals. `n = 252, replace = FALSE` stops with "Cannot upsample"
instead of being capped at 250.
Proposed change: compare with `nInd(x)`. With replacement, upsample to
`n`; without replacement, cap at `nInd(x)` with the existing warning.

**F3 [MEDIUM, confidence: high] — the documented default n = NULL crashes (spec, FS5)**
`R/gl.subsample.ind.r:71-86` — `if (n < 1)` runs before the NULL default
is resolved.
Failure scenario: `gl.subsample.ind(x)` stops with "argument is of length
zero". The documentation promises the smallest population size (by
population) or `nInd(x)`.
Proposed change: resolve the default first, then check that `n` is a
single number of at least 1, with a `stop(error(...))` naming it.

**F4 [MEDIUM, confidence: high] — `method = "pic"` ranks by stale information content (DAT4)**
`R/gl.subsample.loci.r:97-102,124-128` — `AvgPIC` (SNP) or `PIC`
(SilicoDArT) is read as stored, even when its flag says the value no
longer describes the individuals present.
Failure scenario: testset.gl with 20 of 30 populations dropped (flag
`AvgPIC = FALSE`): only 38 of the 50 loci chosen are in the true top 50
after recalculation.
Proposed change: when the flag is not TRUE, recalculate with
`utils.recalc.avgpic()` before ranking.
**Consequence: "pic" selections change for objects whose PIC values are
out of date.**

**F5 [LOW, confidence: high] — method matching and plain genlights (DAT5)**
`R/gl.subsample.loci.r:48-83` — `method = "PIC"` is rejected by the check
at line 74 although line 124 accepts it, so it silently becomes random.
The monomorph checks read `loc.metrics.flags$monomorphs`, which a plain
adegenet genlight lacks.
Failure scenario: `method = "PIC"` samples at random (warning only at
`verbose >= 1`). A plain genlight stops with "argument is of length zero".
Proposed change: match `method` case-insensitively and stop on an
unknown value; treat a missing flag as "not checked" (`isTRUE()`).

**F6 [LOW, confidence: high] — missing `n` in gl.subsample.loc gives an opaque error (FS5, DOC5 proposed rule)**
`R/gl.subsample.loc.r:10,31,51` — `n` has no default, but the
documentation says `[default NULL]`. A missing `n` fails inside the range
check with a bare call trace. Line 21 has a stray `#' #'` in the details.
Failure scenario: `gl.subsample.loc(x)` stops with an error naming no
argument.
Proposed change: document `n` as `[required]`, and stop with "n must be
supplied" when it is missing. Fix the stray tag.

**F7 [LOW, confidence: high] — error idiom and documentation (VRB2, DOC1)**
`R/gl.subsample.ind.r:72-73,142-143` use `cat(error(...)); stop()`, which
prints the message and then stops with an empty error.
`R/gl.subsample.loci.r:1-27` has no `@name`, `@title` or `@family`, and
its example labelled "Tag P/A data" runs on SNP data (testset.gl).
Failure scenario: `tryCatch()` receives an empty error message; the help
page is filed under no family.
Proposed change: `stop(error(...))`; add `@name`/`@title`/`@family data
manipulation`; run the P/A example on testset.gs. Mostly docs.

**F8 [INFO, confidence: high] — two functions for random locus subsampling (design)**
`gl.subsample.loci(method = "random")` does what
`gl.subsample.loc(replace = FALSE)` does; the only unique capability of
`gl.subsample.loci` is `method = "pic"` (and `mono.rm`).
Failure scenario: users choose between near-identical names with
different defaults (`replace = TRUE` in one, no replacement in the other)
and different argument checks (one caps `n`, the other stops).
Proposed change: add `method = c("random", "pic")` and `mono.rm` to
`gl.subsample.loc`, and make `gl.subsample.loci` a deprecated wrapper
that calls it, as was done for `gl.filter.excess.het`.
**Consequence: `gl.subsample.loci` warns that it is deprecated;
`gl.subsample.loc` gains two arguments with defaults that keep current
calls unchanged.**

## Proposed changes

1. Fix the upsampling loop count (F1).
   **Consequence: `by.pop = TRUE` with `n` above a population's size
   returns n individuals per population instead of more.**
2. Compare with `nInd`, upsample to n with replacement, cap at `nInd`
   without (F2).
   **Consequence: `by.pop = FALSE` returns n individuals where it
   returned nLoc(x).**
3. Resolve `n = NULL` before checking it; validate `n` (F3).
4. Recalculate stale PIC before ranking (F4).
   **Consequence: "pic" selections change for objects with out-of-date
   PIC values.**
5. Case-insensitive `method`, stop on an unknown value, `isTRUE()` flag
   checks (F5).
6. Clear error for missing `n` in `gl.subsample.loc`; docs (F6).
7. `stop(error(...))` idiom and `gl.subsample.loci` roxygen (F7).
8. Merge `gl.subsample.loci` into `gl.subsample.loc` as a deprecated
   wrapper (F8).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run on all three
- Spec: returned sizes against the requested `n` for every branch
  (global/by population × with/without replacement × below/above size) —
  run (F1, F2, F3)
- DAT2 metadata sync: loc.metrics aligned with loci after `loc` and
  `loci` sampling; ind.metrics rows = nInd after `ind` sampling — run,
  passes
- DAT1 ploidy: SilicoDArT through `gl.subsample.ind` keeps ploidy 1 — run
- FBM path (DAT6): run — `gl.subsample.loc`, `gl.subsample.loci` and
  `gl.subsample.ind` complete on FBM-backed testset.gl (the ind count
  shows F1)
- Plain genlight (DAT5): `loc` works; `loci` crashes (F5)
- Randomness: sizes and alignment are checked, not the sampling
  distribution (`sample()` is used directly)
- Google Group / GitHub issues: not searched (not available: no browser
  session in this run)

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence (per-population sizes change) approved explicitly |
| 2 | approved | Luis | consequence (global sizes change) approved explicitly |
| 3 | approved | Luis | |
| 4 | approved | Luis | consequence ("pic" selections change) approved explicitly |
| 5 | approved | Luis | |
| 6 | approved | Luis | |
| 7 | approved | Luis | |
| 8 | approved | Luis | consequence (gl.subsample.loci deprecated; two new arguments in gl.subsample.loc) approved explicitly |

## Outcome

- Change 1 (F1): `seq_len(k - 1)`. platypus.gl `by.pop = TRUE`: `n = 70`
  gives 70/70/70 (was 93/87/70), `n = 100` gives 100/100/100 (was
  123/117/141); FBM testset.gl `n = 10`: 300 (was 326). Baseline flipped
  (approved diff).
- Change 2 (F2): `nInd`. testset.gl `by.pop = FALSE`: `n = 300` and `600`
  with replacement give 300 and 600 (was 255); `n = 252` without
  replacement gives 250 (was an error). Baseline flipped (approved diff).
- Change 3 (F3): the default is resolved before the check;
  `gl.subsample.ind(testset.gl)` returns one individual per population
  (30; the smallest population has 1). Invalid `n` stops naming it.
- Change 4 (F4): stale PIC recalculated via `utils.recalc.avgpic()`; after
  dropping 20 populations, "pic" picks 50 of the true top 50 (was 38).
- Change 5 (F5): `tolower(method)`; an unknown method stops; plain genlight
  works (via `gl.subsample.loc`, which does not read the monomorph flag).
- Change 6 (F6): missing `n` stops with "n, the number of loci to
  subsample, must be supplied"; docs mark `n` [required].
- Change 7 (F7): `stop(error(...))` in `gl.subsample.ind`;
  `gl.subsample.loci` roxygen gains `@name`/`@title`/`@family`, and its P/A
  example uses testset.gs.
- Change 8 (F8): `gl.subsample.loc(method, mono.rm)` added after
  `error.check`, before `verbose`. `gl.subsample.loci` is a deprecated
  wrapper that keeps its stop on out-of-range `n` and records its own call
  in the history (tested: the entry replays to the same loci). Callers
  (gl.nhybrids, gl2snapper, gl.tree.fitch, gl2paup.parsimony, gl2bpp
  tests) all use named arguments; unaffected.
- `gl.subsample.loc` with `mono.rm = TRUE` records one history entry (the
  internal `gl.filter.monomorphs` entry is not kept).
- Tests: ind 14, loc 16, loci 12 expectations, all pass.
- `devtools::document()` also rewrote unrelated Rd files (cross-link
  drift, PR #421); reverted.
- PR: #426.

```json
{
  "function": "gl.subsample.ind",
  "companion": ["gl.subsample.loc", "gl.subsample.loci"],
  "package": "dartR.base",
  "family": "modify",
  "skill_version": "2.0.0",
  "commit": "f9f1be8",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "spec", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "spec", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DAT4", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "DAT5", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "FS5", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "approved", "change": 7},
    {"id": "F8", "severity": "INFO", "confidence": "high", "rule": "design", "status": "approved", "change": 8}
  ],
  "coverage_skipped": ["sampling distribution not tested", "forum/issues: no browser session"],
  "status": "pr-open",
  "pr": 426
}
```
