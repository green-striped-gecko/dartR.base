# Review: gl.edit.recode.ind (dartR.base)
- Family mode: modify (interactively recodes individual labels; deletes flagged individuals; optional recalc/mono.rm)
- Date: 2026-09-15
- Reviewer: Claude (Claude Opus 4.8), dartr-function-review v1.0.0
- Package commit: f02bd34 (upstream/dev; working copy verified identical by `git diff upstream/dev -- R/gl.edit.recode.ind.r`)
- Datasets: testset.gl (SNP), testset.gs (SilicoDArT) — dartR.data 1.2.5
- Baseline: tests/testthat/test-gl.edit.recode.ind.R (new file, snapshot captured pre-review; 10 assertions, all passing)

## Verdicts

**Standards: Needs work** — the FS backbone is present, but the object's
metric-validity flags depend on the verbosity level, the `outpath` parameter
is computed and then ignored, and a missing metrics flag crashes the run.

**Spec: Needs work** — the interactive recode/delete behaviour is sound, but
the documented defaults contradict the actual defaults, `outpath` does not
work as documented, and the roxygen carries copy-paste errors from the
population-recode sibling.

## Blast radius

`gl.edit.recode.ind` has **no internal callers** anywhere in the
dartRverse — it is a user-facing interactive function (it calls `edit()` to
open a spreadsheet editor).

## Independent verification (spec axis)

The `edit()` step no-ops under a non-interactive session, so the surrounding
logic (which runs regardless of what the user edits) was exercised
directly. The recode-application loop and the delete path were verified by
inspection; the flag-handling, outpath and default-argument findings were
confirmed empirically (tests below).

## Findings

**F1 [HIGH, confidence: high] — metric-validity flags depend on the verbosity level (DAT4; verbosity-dependent-results class)**
`R/gl.edit.recode.ind.r:197-202` — when `recalc = FALSE` (the default), the
`x <- utils.reset.flags(x, verbose = 0)` call that marks the locus metrics
as stale (individuals may have been deleted) sits **inside** the
`if (verbose >= 2)` block. So the flags are reset only at `verbose >= 2`.
Failure scenario: `gl.edit.recode.ind(x, recalc = FALSE)` at `verbose = 0`
(or 1) returns an object whose `loc.metrics.flags` still read TRUE — the
metrics claim to be current even though individuals were recoded/deleted —
whereas the same call at `verbose = 2` returns an object with the flags
correctly reset to FALSE. Confirmed: `CallRate` flag TRUE at verbose 0,
FALSE at verbose 2 (test 2). Downstream functions that trust the flags
(e.g. deciding whether to recalculate) then behave differently for the same
data depending only on how loudly this function was run.
Proposed change: move `utils.reset.flags()` out of the verbose gate so it
always runs when `recalc = FALSE`; keep only the message gated.

**F2 [MEDIUM, confidence: high] — a missing monomorphs flag crashes the run (DAT5)**
`R/gl.edit.recode.ind.r:183` —
`if (x@other$loc.metrics.flags$monomorphs == FALSE)` evaluates
`if (logical(0))` when the flag is absent (NULL), raising "argument is of
length zero".
Failure scenario: a genlight not built by dartR, or one whose flags were
not fully populated, crashes here regardless of the recode. Confirmed
(test 3). This is the same class fixed in gl.filter.ld.
Proposed change: `if (!isTRUE(x@other$loc.metrics.flags$monomorphs))`.

**F3 [MEDIUM, confidence: high] — documented defaults contradict the signature (DOC5 (proposed rule))**
`R/gl.edit.recode.ind.r:18-19` document `recalc` and `mono.rm` as
`[default TRUE]`, but the signature (lines 74-75) sets both to `FALSE`.
Failure scenario: a user reading the manual expects monomorphic loci
removed and metrics recalculated by default after recoding; neither
happens. Confirmed via `formals()` (test 4).
Proposed change: correct the documentation to `[default FALSE]` (the
signature default is preserved — changing the defaults themselves would be
an API change affecting existing callers, and is noted below as a separate
consideration).

**F4 [MEDIUM, confidence: high] — outpath is computed but the recode file ignores it (FS7)**
`R/gl.edit.recode.ind.r:93,132-138` — `outfilespec <- file.path(outpath,
out.recode.file)` is built but never used; `write.table(..., file =
out.recode.file)` writes to the bare filename, so the file lands in
`getwd()`, not the `outpath` the user set (and the docs at line 45-46 tell
them to set).
Failure scenario: `gl.edit.recode.ind(x, out.recode.file = "r.csv",
outpath = "some/dir")` writes `r.csv` to the working directory, not
`some/dir`. Confirmed (test 1).
Proposed change: `write.table(..., file = outfilespec)`.

**F5 [LOW, confidence: high] — roxygen copy-paste errors and stale text (DOC1, DOC2, DOC5 (proposed rule))**
`R/gl.edit.recode.ind.r:16` — `@param outpath "Directory to save the plot
RDS files"` (this function writes a recode CSV, not plot RDS files);
`:20-22` the `verbose` text is the pre-DOC2 form ("progress but not
results"); the in-body messages say "populations" for individuals
(`:152` "Assigning new populations to x", `:161` "Remove populations
flagged for deletion") — copied from `gl.edit.recode.pop`; the `@details`
(`:38-43`) states that mono.rm/recalc are "not available for Tag P/A data
(SilicoDArT)", but the code applies them to SilicoDArT without gating (runs
without error — the restriction is documented but not implemented).
Proposed change: correct the `outpath` description, the individual-vs-
population message wording, and the verbose text; either gate
mono.rm/recalc for SilicoDArT or drop the "not available" claim.

## Proposed changes

1. Move `utils.reset.flags()` out of the `if (verbose >= 2)` block so the
   metric flags are reset whenever `recalc = FALSE`, independent of
   verbosity (F1).
   **Consequence: for `recalc = FALSE` at verbose 0/1 the returned object's
   loc.metrics.flags are now reset to FALSE (as they already are at verbose
   2) — i.e. the object correctly reports its metrics as stale.**
2. Write the recode file to `outfilespec` so `outpath` is honoured (F4).
3. Null-safe monomorphs-flag check (F2).
4. Correct the documented `recalc`/`mono.rm` defaults to FALSE (F3).
5. Roxygen and message fixes: outpath description, individual-vs-population
   wording, verbose text, and the SilicoDArT restriction claim (F5).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP (n/a), PLT (n/a), STY — run
- Spec: outpath handling, verbose-dependent flag state, missing-flag crash,
  default arguments, history append — run; the interactive recode/delete
  application loop verified by inspection
- Interactive `edit()` path: SKIPPED for unit testing — it no-ops
  non-interactively; the recode/delete logic downstream of it was inspected,
  not driven with scripted edits
- Caller survey (API3): no internal callers; user-facing interactive
  function
- GitHub issues / Google Group: SKIPPED — not queried this session

Out-of-scope note (no action): the `recalc`/`mono.rm` defaults themselves
(FALSE) leave stale metrics after deletion unless the user opts in; whether
the defaults should be TRUE is an API decision for the sibling
recode/edit family as a whole, not this function alone.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur | consequence (flags reset at verbose 0/1 too) approved explicitly |
| 2 | approved | Arthur | |
| 3 | approved | Arthur | |
| 4 | approved | Arthur | docs corrected to FALSE (defaults unchanged) |
| 5 | approved | Arthur | |

## Outcome

(pending Phase C)

```json
{
  "function": "gl.edit.recode.ind",
  "package": "dartR.base",
  "family": "modify",
  "skill_version": "1.0.0",
  "commit": "f02bd34",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT4", "status": "proposed", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DAT5", "status": "proposed", "change": 3},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "proposed", "change": 4},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS7", "status": "proposed", "change": 2},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "proposed", "change": 5}
  ],
  "coverage_skipped": ["interactive edit() path not unit-tested", "GitHub issues not queried", "Google Group not queried"],
  "status": "awaiting-approval",
  "pr": null
}
```
