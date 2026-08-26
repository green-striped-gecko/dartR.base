# Review: gl.drop.pop (dartR.base)

- Family mode: modify
- Date: 2026-08-26
- Reviewer: Claude (claude-sonnet-5), dartr-function-review v1.4.0
- Package commit: ba7145a
- Datasets: testset.gl, testset.gs
- Baseline: tests/testthat/test-gl.drop.pop.R (snapshot captured pre-review)

## Setup note (not a finding)

This checkout had none of the campaign infrastructure the skill assumes yet:
no `dev_<name>` branch (work has been happening directly on `dev`), no
`tests/testthat/` scaffold, and no `function-review/manifest.csv`. Created
the `testthat` scaffold, the manifest (this function's row only), and this
`function-review/` tree fresh as part of this review rather than blocking on
their absence. Stayed on `dev` rather than inventing a personal branch
unasked.

## Verdict

**Standards: Needs work** — five specific, bounded issues: a missing
`invisible()` wrap, one roxygen tag out of house order, an uncoloured
console-output block, a documented `verbose` text deviation shared with much
of the codebase, and a locus-metric flag that understates confidence after
`mono.rm=TRUE`. Nothing structural.

**Spec: Ready** — every behavioural claim in the docs (population removal,
`mono.rm`, `recalc`, `as.pop`, dual SNP/SilicoDArT support) held up against
the characterization test; no discrepancy found between documented and
actual behaviour.

## Findings

**F1 [MEDIUM, confidence: high] — missing `invisible()` on return (FS10)**
`R/gl.drop.pop.r:213` — `return(x)` is not wrapped in `invisible()`, unlike
the house convention for `modify`-class functions. Confirmed: calling
`gl.drop.pop(x, pop.list = "EmsubRopeMata", verbose = 0)` without assignment
prints the full 22-line "DARTR OBJECT" summary to console.
Failure scenario: a user calling the function for its side effect inside a
loop or pipeline gets unexpected console spam on every iteration.
Proposed change: wrap the final `return(x)` in `invisible()`.

**F2 [MEDIUM, confidence: high] — `@return` out of house order (DOC1)**
`R/gl.drop.pop.r:46` — `@return` is the last roxygen tag, after `@export`,
rather than immediately following the `@param` block. Actual order: `@param`
→ `@author` → `@examples` → `@seealso` → `@export` → `@return`; expected:
`@param` → `@return` → `@author` → `@examples` → `@export`.
Failure scenario: none functional (roxygen2 renders regardless of tag
order), but it's a mechanical, checkable deviation from house style and
inconsistent with sibling functions in the same family.
Proposed change: move `@return` to directly follow the `@param` block.

**F3 [LOW, confidence: high] — uncoloured console output in the verbose>=3 block (VRB2)**
`R/gl.drop.pop.r:128, 168-192` — roughly 20 lines (the deletion notice and
the "Summary of recoded dataset" block) use plain `cat()`/`cat(paste(...))`
rather than the `report()`/`warn()` colour helpers every other message in
this file uses. Confirmed: `verbose = 0` correctly suppresses the block (0
output lines); `verbose = 3` prints 19 lines, none colour-coded.
Failure scenario: purely cosmetic — gating works correctly, only the colour
convention is inconsistent with the rest of the function and the package.
Proposed change: wrap each `cat()` call in `report()`.

**F4 [LOW, confidence: medium] — `verbose` param text deviates from DOC2's canonical wording (proposed rule)**
`R/gl.drop.pop.r:24-26` — reads "2, progress but not results" / "[default 2
or as specified using gl.set.verbosity]" versus DOC2's canonical "2, brief
progress messages" / "[default 2, unless specified using gl.set.verbosity]".
Note for the skill maintainer, not just this function: this precise
deviation (and a third variant, "2, progress log") is extremely common
across dartR.base functions reviewed informally outside this campaign so
far — worth auditing whether DOC2's canonical text is actually representative
of current code before treating deviations from it as per-function findings,
since [confirmed] status is supposed to mean "matches current code."
Proposed change: align wording with DOC2, but only as part of a wider pass —
not worth a one-off fix here if the canonical text itself needs revisiting.

**F5 [LOW, confidence: medium] — `monomorphs` flag left FALSE even when `mono.rm=TRUE` just confirmed none remain (DAT4-adjacent)**
`R/gl.drop.pop.r:137-164` — with `recalc = FALSE` (the default), the final
`utils.reset.flags()` call unconditionally resets the `monomorphs` flag to
FALSE, even on the `mono.rm = TRUE` path where `gl.filter.monomorphs()` was
just run. Confirmed empirically: after `mono.rm = TRUE, recalc = FALSE`, no
monomorphic loci remain in the output, but `loc.metrics.flags$monomorphs`
still reads FALSE.
Failure scenario: none in the data itself — `utils.reset.flags()` is doing
exactly what it says (blanket-invalidating all locus metrics after
individual/population deletion, since most other metrics genuinely did
become stale). This just means a user reading the flag alone would
under-trust a metric that was, in this one case, just correctly verified.
Proposed change: none proposed — flagging for awareness; may be acceptable
given `utils.reset.flags()`'s intentionally blanket scope.

## Cleared during verification (checked, not a defect)

- **The recurring `loc.metrics.flags$monomorphs == FALSE` crash** seen
  repeatedly elsewhere in this codebase (`NULL == FALSE` → "argument is of
  length zero") does **not** apply here: line 137 unconditionally assigns a
  definite `FALSE` literal immediately before the line 147 check, so the
  flag can never be `NULL` at that point. Verified empirically by forcing
  the incoming object's flag to `NULL` before calling — no crash.
- Ploidy, individual counts, and history-append all confirmed correct for
  both SNP and SilicoDArT via the characterization test.
- The `for (case in pop.list)` loop mutating `pop.list` mid-iteration
  (lines 112-119) is safe in R — the iteration sequence is fixed at loop
  start and unaffected by the later reassignment.
- No `accept =` restriction on `utils.check.datatype()` is correct and
  intentional: the function's own docs state it supports both SNP and
  SilicoDArT, and both pass in the characterization test.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: behaviour vs roxygen on testset.gl / testset.gs — run
- FBM path (DAT6): SKIPPED — no FBM-backed fixture readily available for
  this function
- Sibling `dartR.*` / dartr2shiny caller grep (API1-3): not run — no
  API-affecting change proposed, so the escalation gate does not apply

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur | |
| 2 | approved | Arthur | |
| 3 | approved | Arthur | |

F4 and F5 carried no proposed fix (flagged for awareness only) and were not
offered as approvable changes.

## Outcome

- **F1 (invisible)**: applied at `R/gl.drop.pop.r:214`. Verified: an
  unassigned call now prints 0 lines (was 22); assigned calls' return value
  unchanged (`class`, `nInd` confirmed identical).
- **F2 (@return order)**: applied — moved to directly follow `@param`.
  `devtools::document()` re-run; `man/gl.drop.pop.Rd` unchanged (roxygen2
  renders `.Rd` sections in a fixed order regardless of source comment
  order, so this was source-hygiene only, as expected).
- **F3 (colour the summary block)**: applied to all ~20 `cat()` calls at
  lines 129, 171-192. Verified: `verbose = 3` end-to-end run on both the
  `as.pop`-null and `as.pop`-set branches completes cleanly with the
  expected coloured output.
- Characterization test (`tests/testthat/test-gl.drop.pop.R`): all 13 cases
  still pass unchanged — no behavioural diff, matching expectations for
  print/doc-only changes.
- PR: not yet opened (pending pre-push confirmation).

```json
{
  "function": "gl.drop.pop",
  "package": "dartR.base",
  "family": "modify",
  "skill_version": "1.4.0",
  "commit": "ba7145a",
  "verdict_standards": "needs_work",
  "verdict_spec": "ready",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "FS10", "status": "pending"},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DOC1", "status": "pending"},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "pending"},
    {"id": "F4", "severity": "LOW", "confidence": "medium", "rule": "DOC2", "status": "pending"},
    {"id": "F5", "severity": "LOW", "confidence": "medium", "rule": "DAT4", "status": "pending"}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture available"],
  "status": "awaiting-approval",
  "pr": null
}
```
