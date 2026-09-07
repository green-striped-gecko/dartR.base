# Review: gl.report.monomorphs (dartR.base)

- Family mode: report
- Date: 2026-08-29
- Reviewer: Claude (claude-fable-5), dartr-function-review v2.0.0
- Package commit: 2739215
- Datasets: testset.gl, testset.gs
- Baseline: tests/testthat/test-gl.report.monomorphs.R (snapshot captured pre-review)

## Verdict

**Standards: Needs work** — the results block ignores `verbose` (8 lines
at `verbose = 0`), plus the routine doc items.

**Spec: Rework** — the function's central contract is broken: documented
as returning "An unaltered genlight object", it actually returns the
FILTERED object (255 loci in, 111 out on testset.gl — 144 monomorphic
loci silently deleted, plus a `gl.drop.loc` history entry the user never
made). The counts themselves are correct (verified independently on both
datatypes).

## Findings

**F1 [HIGH, confidence: high] — a report function that silently filters the returned object (Spec / report-mode contract / FS8)**
`R/gl.report.monomorphs.r:99, 118` — the function drops the monomorphic
loci from `x` (to derive the polymorphic count) and then returns that
filtered `x`; the untouched copy saved in `hold` (line 50) is used only
for the "No. of loci" line. Confirmed: 255 loci in, 111 returned;
history grew by 1 (the internal `gl.drop.loc` call). The docs — and the
report-family convention stated in sibling docs ("To avoid issues from
inadvertent use of this function in an assignment statement, the
function returns the genlight object unaltered") — promise no
modification.
Failure scenario: `gl <- gl.report.monomorphs(gl)` — the natural,
docs-sanctioned idiom — silently deletes 144 loci from the user's
object. Filtering is `gl.filter.monomorphs`'s job.
Proposed change: return `invisible(hold)` — the unaltered input. The
internal filtered copy remains for the counts. This also removes the
spurious history append (FS8: report functions do not append history).

**F2 [MEDIUM, confidence: high] — results block prints unconditionally (VRB5/VRB3)**
`R/gl.report.monomorphs.r:103-108` — confirmed: 8 lines at
`verbose = 0`. Same defect class fixed in gl.report.bases (#244) and
gl.report.callrate (#247), both gated at `verbose >= 1`.
Proposed change: gate the block at `verbose >= 1`.

**F3 [LOW, confidence: high] — documentation (DOC7/DOC1)**
`@author` lacks the Author(s) part; tags out of ratified DOC1 order
(`@details` after `@param`, `@return` last after `@export`).
Proposed change: add `Author(s): Arthur Georges.`; reorder (docs only).

**F4 [notes, no fix proposed]**
(a) The printed "Monomorphic loci" count includes the all-NA loci (which
are also reported on their own line) — consistent with
`gl.filter.monomorphs`'s removal set, but the label could be read as
excluding them.
(b) `verbose` doc clause is the "[default NULL, unless...]" variant
(DOC2 canon note, as all siblings).
(c) Per-locus loop over `as.matrix(x)` — fine at report scale; DAT6
consideration only for FBM-backed objects.

## Cleared during verification (checked, not a defect)

- Monomorphic/all-NA counts match independent computations exactly on
  both datatypes (uniform dosage 1 correctly NOT counted as monomorphic).
- The SNP/SilicoDArT dispatch uses the correct dosage sets (0/2 vs 0/1).
- No plot bundle — the graphics half of VRB5 does not apply.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: counts vs independent computation on testset.gl / testset.gs —
  run
- FBM path (DAT6): SKIPPED — no FBM-backed fixture readily available
- Caller grep (for F1/F2 behaviour changes): run in Phase C, results
  recorded in Outcome

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur | unaltered-return consequence approved (hidden-filtering reliance must move to gl.filter.monomorphs) |
| 2 | approved | Arthur | verbose=0 silence consequence approved |
| 3 | approved | Arthur | |

F4 notes only, not offered.

## Outcome

- **F1**: applied — `invisible(hold)`. Verified: `expect_identical()`
  passes on both datatypes; no history growth; counts still correct
  (the internal filtered copy remains for the report only).
- **F2**: applied — results block gated at `verbose >= 1`. Verified:
  0 lines at `verbose = 0` (was 8); the two baseline count tests moved
  from verbose=0 to verbose=1 (expected diff, maps to F2).
- **F3**: applied — Author(s) line added, tags reordered per ratified
  DOC1. `man/gl.report.monomorphs.Rd` regenerated.
- Caller grep: zero call sites anywhere that assign the return or parse
  output — only unassigned roxygen examples and one `eval=FALSE`
  tutorial chunk whose surrounding prose already describes the function
  as a pure report paired with gl.filter.monomorphs for removal. No
  dartr2shiny. NEWS.md entry added.
- Characterization test: 7/7 pass. `verbose = 3` end-to-end clean on
  both datatypes.
- PR: [green-striped-gecko/dartR.base#249](https://github.com/green-striped-gecko/dartR.base/pull/249)

```json
{
  "function": "gl.report.monomorphs",
  "package": "dartR.base",
  "family": "report",
  "skill_version": "2.0.0",
  "commit": "2739215",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "FS8", "status": "approved"},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "approved"},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "approved"},
    {"id": "F4", "severity": "LOW", "confidence": "medium", "rule": "DOC2", "status": "note_only"}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture available"],
  "status": "pr-open",
  "pr": 249
}
```
