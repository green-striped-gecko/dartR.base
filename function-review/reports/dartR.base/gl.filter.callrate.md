# Review: gl.filter.callrate (dartR.base)

- Family mode: modify
- Date: 2026-08-28
- Reviewer: Claude (claude-fable-5), dartr-function-review v2.0.0
- Package commit: b04c5bb
- Datasets: testset.gl, testset.gs
- Baseline: tests/testthat/test-gl.filter.callrate.R (snapshot captured pre-review)

## Verdict

**Standards: Rework** — `method='pop'` silently corrupts the entire
`loc.metrics` data frame (every row all-NA) via name-based row indexing
against non-name rownames; plus the recurring monomorphs-flag NULL crash,
a stale monomorphs flag after individual filtering, and boundary/verbosity
defects.

**Spec: Needs work** — the filter itself filters correctly on all three
methods (verified against independent computations), but the
deleted-individuals listing contradicts the retention rule at the exact
threshold, and the docs mis-describe the input class, the return value,
and the threshold type.

## Findings

**F1 [HIGH, confidence: high] — `method='pop'` corrupts all locus metrics (DAT2/DAT3)**
`R/gl.filter.callrate.r:514-515` — `x2@other$loc.metrics <-
x@other$loc.metrics[locall,]` indexes data-frame ROWS by locus NAME, but
`rownames(loc.metrics)` are not locus names (confirmed:
`identical(rownames, locNames)` is FALSE on testset.gl). Character
indexing against non-matching rownames yields all-NA rows. Confirmed
empirically: after `method='pop', threshold=0.8`, all 127 retained loci
have every metric NA (AlleleID included). The dartR `[` operator on the
previous line had already subset `loc.metrics` correctly — the manual
reassignment then overwrites it with garbage. A correctly-computed
positional `index` sits unused on line 512.
Failure scenario: any pop-filtered object silently loses every locus
metric; all downstream metric-dependent analyses on it are garbage or
crash.
Proposed change: subset positionally
(`keep <- locNames(x) %in% locall`), matching DAT3 and the loc branch.

**F2 [HIGH, confidence: high] — monomorphs-flag NULL crash (DAT4/FS5)**
`R/gl.filter.callrate.r:136` — `x@other$loc.metrics.flags$monomorphs ==
FALSE` throws "argument is of length zero" when the flag is unset — the
same recurring defect class fixed across the `utils.recalc.*` family and
`gl.subsample.loci`. Confirmed empirically.
Proposed change: `!isTRUE(...)` idiom.

**F3 [MEDIUM, confidence: high] — at-threshold individuals retained but listed as deleted (Spec)**
`R/gl.filter.callrate.r:304 vs 310-311 (and 372 vs 396-397)` — retention
uses `>= threshold`, but the deleted-report and `x3` use `<= threshold`.
Confirmed: with a threshold set to an actual call-rate value, the
individual exactly at the threshold is retained AND printed in the
"Individuals deleted" list.
Proposed change: deleted-set becomes the exact complement
(`< threshold`), both branches. Listing only; retention unchanged.

**F4 [MEDIUM, confidence: high] — monomorphs flag left stale after individual filtering (DAT4)**
`R/gl.filter.callrate.r:333-346 (and 427-440)` — the manual flag-reset
lists (AvgPIC ... allna) omit `monomorphs`. Confirmed: input with
`monomorphs = TRUE`, `method='ind'`, `mono.rm=FALSE` — individuals
removed, monomorphic loci demonstrably present in the result, flag still
TRUE. The drop/keep siblings set this flag FALSE unconditionally after
individual removal.
Proposed change: add `monomorphs <- FALSE` to both ind-branch flag-reset
blocks.

**F5 [LOW, confidence: high] — console hygiene (VRB5/VRB2)**
Ungated `cat("\n")` at line 324 (confirmed: 1 line leaks at
`verbose = 0` when individuals are deleted) and ungated
`cat(report("  Removing monomorphic loci\n"))` at line 419 in the
recursive branch; several warning/report strings are multi-line literals
that print with embedded source indentation (lines 140-141, 163-164,
360-362, 492-493).
Proposed change: gate the two ungated calls; reflow the broken strings
to single lines.

**F6 [LOW, confidence: high] — documentation corrections (DOC1/DOC5/DOC7)**
`@param x` says "the genind object containing the SilocoDArT data"
(wrong class, misspelt SilicoDArT); `@return` says "genlight or genind
object, plus a summary" (no genind, no summary is returned); the
threshold warning says "must be an integer between 0 and 1" (it is a
proportion); `@details` still says results "are saved to the session's
temporary directory" (superseded by plot.file/plot.dir); `@author` lacks
the Author(s) part (DOC7); tags out of ratified DOC1 order.
Proposed change: correct all of the above (docs/messages only).

**F7 [notes, no fix proposed]**
(a) The loc branch's manual flag-reset invalidates per-locus metrics that
locus removal does not actually invalidate — conservative, harmless.
(b) Duplicate `plot.display` suppression (lines 109 and 155-157).
(c) `verbose` param wording variant (DOC2 canon note, as all siblings).
(d) The commented-out CallRate-flag check (line 147) means CallRate is
always recalculated — appropriate for a filter keyed on it.

## Cleared during verification (checked, not a defect)

- All three methods filter correctly (independent recomputation:
  loc-threshold set, ind-threshold set, and every retained locus passing
  the threshold in every population for method='pop').
- Genotype/metadata sync correct in loc and ind branches; ploidy
  preserved on both datatypes; history appended on the returned object.
- Too-stringent ind threshold fails fast with a clear error.
- The `for (i in 1:10)` recursion cap breaks on stability — no runaway.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: all three methods verified against independent computations on
  testset.gl / testset.gs — run
- FBM path (DAT6): SKIPPED — no FBM-backed fixture readily available
- Caller grep (for F1/F3/F4 behaviour changes): run in Phase C, results
  recorded in Outcome

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur | pop-branch metadata consequence approved |
| 2 | approved | Arthur | crash -> run consequence approved; on verification failure (crash re-fired in the callee) Arthur directed: "fix utils.recalc.callrate" — the identical one-line guard included in this PR as a declared dependency fix |
| 3 | approved | Arthur | listing-boundary consequence approved |
| 4 | approved | Arthur | flag-invalidation consequence approved |
| 5 | approved | Arthur | |
| 6 | approved | Arthur | |

F7 notes only, not offered.

## Outcome

- **F1**: applied — pop branch now subsets `loc.metrics` positionally.
  Verified: 0 NA rows (was 127/127); metrics rows are exactly the source
  rows for the retained loci; per-population threshold check passes
  independently.
- **F2**: applied in `gl.filter.callrate` AND, per Arthur's direction
  after the verification stop, the identical `!isTRUE()` guard applied
  to the callee `utils.recalc.callrate` (its pre-campaign fix existed
  only on an unpushed local branch). Verified: NULL-flag object now runs
  clean end to end.
- **F3**: applied — deleted-set and listing use `< threshold` in both
  branches. Verified: at-threshold individual retained and no longer
  listed.
- **F4**: applied — `monomorphs <- FALSE` added to both ind-branch
  flag-reset blocks. Verified with a monomorphs-flag-TRUE input.
- **F5**: applied — `verbose = 0` fully silent (was 1 leaked line);
  recursive-branch message gated; broken strings reflowed.
- **F6**: applied — all documentation/message corrections, tags
  reordered per ratified DOC1. `man/gl.filter.callrate.Rd` regenerated.
- Caller grep: only two genuine call sites anywhere (the function's own
  pop-branch self-call with method='loc', and `gl.mahal.assign` at
  method='loc'/verbose=0 reading only the genotype matrix) — neither
  affected. `utils.recalc.callrate`'s widening is unreachable by any
  in-repo caller (all five call sites guarantee a non-NULL flag). No
  dartr2shiny. NEWS.md entry added.
- Addendum note (scope rule): `utils.recalc.avgpic` upstream still
  carries the unguarded `== FALSE` idiom (its fix likewise stranded on
  the unpushed local branch) — one-line note for its own future review,
  not touched here.
- Characterization test: 53/53 pass. `verbose = 3` end-to-end clean on
  all three methods plus the recursive SilicoDArT path.
- PR: recorded below after opening.

```json
{
  "function": "gl.filter.callrate",
  "package": "dartR.base",
  "family": "modify",
  "skill_version": "2.0.0",
  "commit": "b04c5bb",
  "verdict_standards": "rework",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT2", "status": "approved"},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DAT4", "status": "approved"},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved"},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DAT4", "status": "approved"},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "VRB5", "status": "approved"},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved"},
    {"id": "F7", "severity": "LOW", "confidence": "medium", "rule": "DOC2", "status": "note_only"}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture available"],
  "status": "pr-open",
  "pr": null
}
```
