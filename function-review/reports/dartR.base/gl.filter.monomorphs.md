# Review: gl.filter.monomorphs (dartR.base)

- Family mode: modify
- Date: 2026-08-29
- Reviewer: Claude (claude-fable-5), dartr-function-review v2.0.0
- Package commit: 2739215
- Datasets: testset.gl, testset.gs, plus a crafted 3x4 genlight with
  all-heterozygous loci
- Baseline: tests/testthat/test-gl.filter.monomorphs.R (snapshot captured pre-review)

## Verdict

**Standards: Needs work** — the filter itself is correct and clean, but
an unassigned call dumps the object summary (missing `invisible()`), the
history records two entries per call (the internal `gl.drop.loc`
delegation leaks, carrying the full removed-locus list), and there is a
works-by-coincidence parenthesis slip plus routine doc items.

**Spec: Ready** — removal set verified exactly against independent
computations on both datatypes (including the subtle case: uniform-
heterozygous loci correctly retained as polymorphic); metadata sync,
ploidy, and the monomorphs flag all correct; idempotent on a clean
object.

## Findings

**F1 [MEDIUM, confidence: high] — missing `invisible()` on return (FS10)**
`R/gl.filter.monomorphs.r:141` — `return(x)`. Confirmed: an unassigned
call prints the 22-line object summary. Every reviewed sibling filter
uses the invisible convention.
Proposed change: `return(invisible(x))`.

**F2 [MEDIUM, confidence: high] — internal delegation leaks a second history entry (FS8)**
`R/gl.filter.monomorphs.r:110, 133-134` — the `gl.drop.loc()` call
appends its own `match.call()` (including the full `loc.list` of removed
locus names — 144 names on testset.gl) to `@other$history`, and the
function then appends its own entry: two history records for one user
action, the first being an internal implementation detail that bloats
the history and misrepresents what the user called. Confirmed: history
grows by 2 per call.
Proposed change: restore the pre-call history before appending the
function's own entry, so one user action records one entry.

**F3 [LOW, confidence: high] — parenthesis slip that works by coincidence (STY3)**
`R/gl.filter.monomorphs.r:105` — `if (length(loc.list > 0))` computes
the length of a comparison vector rather than comparing the length; it
happens to behave identically (`length(v > 0) == length(v)`, and
`if (n)` treats nonzero as TRUE) but reads as a defect and would break
under any future edit that relies on it being a logical.
Proposed change: `if (length(loc.list) > 0)`. No behaviour change.

**F4 [LOW, confidence: high] — documentation (DOC7/DOC1)**
`@author` lacks the Author(s) part; detail paragraphs sit inside
`@description` with no `@details` tag; `@return` last after `@export`
(ratified DOC1 order); `@import utils patchwork` and
`@importFrom plyr count` import things this function never uses; the
"Removing monomorphic loci..." report string is a multi-line literal
that prints with embedded source indentation.
Proposed change: correct all (docs/messages only).

**F5 [notes, no fix proposed]**
(a) The verbose>=3 summary prints "Monomorphic loci" and "No. of loci
deleted" as the same number — redundant but accurate.
(b) `verbose` doc clause is the "[default 2, unless...]" variant (DOC2
canon note, as all siblings).

## Cleared during verification (checked, not a defect)

- Removal set exact against independent computation on both datatypes;
  all-heterozygous loci correctly retained (crafted fixture).
- Metadata rows track loci 1:1 after filtering; ploidy preserved;
  individuals/populations untouched.
- `monomorphs` flag correctly set TRUE after the run (this is the one
  function entitled to do so).
- Idempotent: a second run removes nothing and keeps the flag TRUE.
- All progress/summary messages properly verbosity-gated; the only
  verbose=0 output is the unassigned-call dump (F1).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: removal set vs independent computation on testset.gl /
  testset.gs; all-het fixture — run
- FBM path (DAT6): SKIPPED — no FBM-backed fixture readily available
- Caller grep (for F1/F2 behaviour changes): run in Phase C, results
  recorded in Outcome

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur | invisible-return consequence approved |
| 2 | approved | Arthur | one-history-entry consequence approved |
| 3 | approved | Arthur | |
| 4 | approved | Arthur | |

F5 notes only, not offered.

## Outcome

- **F1**: applied -- `return(invisible(x))`. Verified: unassigned call
  prints 0 lines (was 22).
- **F2**: applied -- pre-call history restored before the function's own
  append. Verified: exactly one new entry, naming gl.filter.monomorphs,
  no gl.drop.loc entry.
- **F3**: applied -- `length(loc.list) > 0`. Removal set unchanged
  (verified against independent computation).
- **F4**: applied -- Author(s) line, details split out, ratified tag
  order, unused imports (`utils`, `patchwork`, `plyr::count`) dropped,
  broken string reflowed. `man/gl.filter.monomorphs.Rd` regenerated.
- Caller grep: ~30 internal callers plus 4 in siblings, ALL assigning
  the return (invisible change unobservable) and none parsing history
  for the removed internal entry (the three real history parsers all
  search for gl.filter.secondaries). No dartr2shiny. NEWS.md entry
  added.
- Characterization test: 17/17 pass, incl. exact removal-set checks on
  both datatypes and the all-heterozygous fixture. `verbose = 3`
  end-to-end clean on both datatypes.
- PR: [green-striped-gecko/dartR.base#251](https://github.com/green-striped-gecko/dartR.base/pull/251)

```json
{
  "function": "gl.filter.monomorphs",
  "package": "dartR.base",
  "family": "modify",
  "skill_version": "2.0.0",
  "commit": "2739215",
  "verdict_standards": "needs_work",
  "verdict_spec": "ready",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "FS10", "status": "approved"},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "FS8", "status": "approved"},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "STY3", "status": "approved"},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "approved"},
    {"id": "F5", "severity": "LOW", "confidence": "medium", "rule": "DOC2", "status": "note_only"}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture available"],
  "status": "pr-open",
  "pr": 251
}
```
