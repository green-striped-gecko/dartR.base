# Review: gl.filter.allna (dartR.base)

- Family mode: modify
- Date: 2026-08-29
- Reviewer: Claude (claude-fable-5), dartr-function-review v2.0.0
- Package commit: 2739215
- Datasets: testset.gl, testset.gs, plus a crafted 3x4 genlight with one
  all-NA individual (testset objects have all-NA loci but no all-NA
  individuals)
- Baseline: tests/testthat/test-gl.filter.allna.R (snapshot captured pre-review)

## Verdict

**Standards: Needs work** — flag invalidation and history are gated on
locus-count change only, so removing an all-NA individual silently
modifies the object with stale CallRate flags and no history record;
the datatype check is commented out; the by.pop path leaks an internal
history entry; plus the family's routine items.

**Spec: Needs work** — the filtering itself is exact on all three paths
(verified independently, including the crafted all-NA individual), but
the individuals listing carries the same literal-"NULL" corruption as
its report sibling, and the individual-removal subsetting works only by
the coincidence that no individual is ever named "NULL".

## Findings

**F1 [HIGH, confidence: high] — individuals-only removal leaves stale flags and no history (DAT4/FS8)**
`R/gl.filter.allna.r:204, 219` — both the flag-reset block and the
history append are gated on `nLoc(x) != nLoc(x2)`. When only an all-NA
INDIVIDUAL is removed (no all-NA loci), neither fires. Confirmed on the
crafted fixture: individual removed, `CallRate` flag still TRUE (stale —
removing an individual changes every locus's CallRate denominator),
history grew by 0 (a modified object with no record of the call).
Proposed change: gate both on "anything removed"
(`nLoc changed OR nInd changed`).

**F2 [MEDIUM, confidence: high] — missing `invisible()` on return (FS10)**
`R/gl.filter.allna.r:229` — confirmed: unassigned call dumps the 22-line
object summary.
Proposed change: `return(invisible(x2))`.

**F3 [MEDIUM, confidence: high] — datatype check commented out (FS4/DAT5)**
`R/gl.filter.allna.r:65` — same defect as its report sibling (PR #250).
Confirmed: data.frame input dies with "unable to find an inherited
method..." instead of the standard message.
Proposed change: restore the check as written in the comment.

**F4 [MEDIUM, confidence: high] — by.pop path leaks an internal history entry (FS8)**
`R/gl.filter.allna.r:194, 219-221` — the `gl.drop.loc()` delegation
appends its own entry (with the full removed-locus list) before the
function appends its own: confirmed +2 history entries per by.pop call.
Same defect class fixed in gl.filter.monomorphs (PR #251).
Proposed change: restore the pre-call history before the function's own
append — one user action, one record.

**F5 [LOW, confidence: high] — individuals listing prints literal "NULL"s; subsetting correct only by coincidence (STY3)**
`R/gl.filter.allna.r:135, 145, 154, 162` — `ind.list <- vector("list",
nI)` (the same list-vs-array slip as gl.report.allna, PR #250).
Confirmed listing: "NULL, IND_ALLNA, NULL". The removal itself
(`!x2$ind.names %in% ind.list`) currently works only because `%in%`
coerces the list to c("NULL", ...) and no individual is literally named
"NULL".
Proposed change: `array(NA, nI)` and print the count, mirroring the
sibling fix.

**F6 [LOW, confidence: high] — documentation (DOC1/DOC5)**
Two `@family` tags ("unmatched filter" + "filter functions") — wrong
twice over, since the matched report `gl.report.allna` exists; garbled
`@title` ("populations with all NA across loci" — individuals); detail
paragraphs inside `@description`; two broken multi-line report strings
(lines 73-74, 80-81); `@return` after `@examples` block order issues per
ratified DOC1; `@import utils patchwork` unused.
Proposed change: correct all (docs/messages only); single
`@family matched filter`.

**F7 [notes, no fix proposed]**
(a) The flag-reset deliberately leaves `maf` alone (commented out) —
defensible, since removing all-NA loci/individuals cannot change any
remaining locus's allele frequencies; CallRate is the one metric that
does change.
(b) `verbose` doc clause is the "[default 2, unless...]" variant (DOC2
canon note, as all siblings).

## Cleared during verification (checked, not a defect)

- Removal sets exact on all three paths against independent
  computations: the 3 all-NA loci of testset.gl; the by.pop union (36
  loci); the crafted all-NA individual.
- Metadata rows track loci 1:1; ploidy preserved on both datatypes;
  individual metadata handled by the dartR `[` operator.
- When loci ARE removed, flags reset and history appended correctly
  (the +1 path works; only the individuals-only path is broken).
- All progress/summary messages properly verbosity-gated; the only
  verbose=0 output is the unassigned-call dump (F2).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: removal sets vs independent computation on testset.gl /
  testset.gs / crafted fixture, all three paths — run
- FBM path (DAT6): SKIPPED — no FBM-backed fixture readily available
- Caller grep (for F1-F4 behaviour changes): run in Phase C, results
  recorded in Outcome

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur | flags/history-on-any-removal consequence approved |
| 2 | approved | Arthur | invisible-return consequence approved |
| 3 | approved | Arthur | fail-fast consequence approved |
| 4 | approved | Arthur | single-history-entry consequence approved |
| 5 | approved | Arthur | |
| 6 | approved | Arthur | |

F7 notes only, not offered.

## Outcome

- **F1**: applied -- both gates now fire on any removal (loci OR
  individuals). Verified on the crafted fixture: CallRate flag FALSE,
  history +1 after an individuals-only removal.
- **F2**: applied -- `return(invisible(x2))`. Verified: 0 lines on an
  unassigned call (was 22).
- **F3**: applied -- datatype check restored. Verified: data.frame
  input fails with the standard message. Benign side effect noted: the
  check's own all-NA probe means a depth-2-terminating recursion; extra
  cost and an advisory warning only for verbose >= 2 callers.
- **F4**: applied -- pre-call history restored after the by.pop
  delegation. Verified: +1 entry naming gl.filter.allna only.
- **F5**: applied -- `array(NA, nI)` + count in the listing. Verified:
  one line, names IND_ALLNA, zero "NULL" strings; removal no longer
  depends on the no-individual-named-NULL coincidence.
- **F6**: applied -- single `@family matched filter` (cross-reference
  updates ripple into 12 sibling .Rd files, pure family-index links),
  fixed title, details split, strings reflowed, ratified tag order,
  unused imports dropped.
- Caller grep: 14 live call sites (10 dartR.base incl. the
  compliance/datatype probes, 4 siblings), ALL assigning the return,
  none reading flags/history from the result; the one by.pop caller
  (dartR.sim gl.report.nall) never inspects history. No dartr2shiny.
  NEWS.md entry added.
- Characterization test: 24/24 pass. `verbose = 3` end-to-end clean on
  both datatypes and both paths.
- PR: recorded below after opening.

```json
{
  "function": "gl.filter.allna",
  "package": "dartR.base",
  "family": "modify",
  "skill_version": "2.0.0",
  "commit": "2739215",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT4", "status": "approved"},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "FS10", "status": "approved"},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "FS4", "status": "approved"},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS8", "status": "approved"},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "STY3", "status": "approved"},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved"},
    {"id": "F7", "severity": "LOW", "confidence": "medium", "rule": "DOC2", "status": "note_only"}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture available"],
  "status": "pr-open",
  "pr": null
}
```
