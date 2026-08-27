# Review: gl.report.callrate (dartR.base)

- Family mode: report
- Date: 2026-08-27
- Reviewer: Claude (claude-fable-5), dartr-function-review v2.0.0
- Package commit: b04c5bb
- Datasets: testset.gl, testset.gs
- Baseline: tests/testthat/test-gl.report.callrate.R (snapshot captured pre-review)
- Prior art: Luis's v1.x trial review (later reverted, per the campaign
  README) found the verbose leak, the phantom `by.pop` examples, and the
  unvalidated `method`; all three are still present and re-confirmed here.

## Verdict

**Standards: Needs work** — the results tables ignore `verbose` entirely
(33 lines at `verbose = 0`), `method` is never validated, and a handful of
cosmetic defects (stray `)`, "3r quartile", doc/signature default
mismatch).

**Spec: Needs work** — the docs promise "returns the genlight object
unaltered" but the function returns the object with CallRate recalculated
in place (confirmed with a deliberately stale CallRate), and two
documented examples call a `by.pop` argument that does not exist
(silently swallowed by `...`).

## Findings

**F1 [HIGH, confidence: high] — results tables print unconditionally, ignoring `verbose` (VRB5/VRB3)**
`R/gl.report.callrate.r:156-195, 247-287` — both the `method='loc'` and
`method='ind'` statistics/table blocks are ungated. Confirmed: 33 lines at
`verbose = 0` (plot side is correctly suppressed). Same defect class fixed
in `gl.report.bases` (PR #244), where the ratified gate level was
`verbose >= 1`.
Proposed change: gate both printing blocks at `verbose >= 1`; default
`verbose = 2` output unchanged.

**F2 [HIGH, confidence: high] — documented examples call a nonexistent `by.pop` argument (DOC5)**
`R/gl.report.callrate.r:77-78` — two `@examples` lines pass
`by.pop=TRUE`; the function has no such parameter, and `...` (destined
for `ggsave`) swallows it silently. Confirmed: no error, no per-population
behaviour. A user copying the example believes a per-population report
ran.
Proposed change: remove the two example lines (docs only). The
per-population summary the examples imply already exists inside
`method='ind'`.

**F3 [MEDIUM, confidence: high] — `method` never validated (FS5)**
`R/gl.report.callrate.r:152, 245` — any value other than "loc"/"ind"
(e.g. `method = "pop"`) matches neither branch: confirmed 0 lines of
output, no plot, no error — the function silently does nothing and
returns. House style for recoverable input is warn-and-coerce.
Proposed change: validate up front; unknown values warn and coerce to
"loc".

**F4 [MEDIUM, confidence: high] — returned object is altered, docs say unaltered (Spec / report-mode contract)**
`R/gl.report.callrate.r:145, 341` — `x <- utils.recalc.callrate(x,
verbose = 0)` overwrites `x` before the report, and `invisible(x)` returns
that modified object. Confirmed: input with a deliberately stale
`CallRate` column comes back changed. The recalculation itself is right
(the report must not use stale metrics); returning the mutated object is
what breaks the documented contract and the report-mode rule.
Proposed change: keep the recalc for the report's internal computation
but hold the original input and return THAT unaltered.

**F5 [LOW, confidence: high] — stray `)` printed after the individual listing (VRB2)**
`R/gl.report.callrate.r:287` — `cat("\n)")`; confirmed in output as a
`)` glued to the next line.
Proposed change: `cat("\n")`.

**F6 [LOW, confidence: high] — "3r quartile" typo, both branches (DOC5)**
`R/gl.report.callrate.r:164, 256`.
Proposed change: "3rd quartile".

**F7 [LOW, confidence: high] — `ind.to.list = 0` lists one individual (FS5)**
`R/gl.report.callrate.r:148, 285` — only the upper bound is clamped;
`ind.means[1:0, ]` yields row 1. Confirmed.
Proposed change: clamp to at least 1 and skip the listing entirely when
the caller asks for 0.

**F8 [LOW, confidence: high] — `@param bins` says default 25, signature says 50 (DOC1)**
`R/gl.report.callrate.r:25 vs 103`.
Proposed change: doc corrected to 50 (docs only).

**F9 [LOW, confidence: high] — roxygen tag order vs ratified DOC1 (DOC1)**
`@details` follows `@param`, `@return` sits last after `@export`.
Proposed change: reorder per the ratified house order (docs only).

**F10 [notes, no fix proposed]**
(a) `verbose` param default clause is the "[default 2, unless...]"
variant (canon since v1.6.2 is the NULL -> global -> 2 cascade) — same
widespread-deviation note as all reviewed functions.
(b) `aggregate(ind.means ~ pop(x), x, mean)` passes the genlight as the
`data` argument — works via formula-environment lookup, fragile-looking
but verified correct in the baseline.
(c) No `nLoc == 0` guard (NaN "Missing Rate Overall" on a degenerate
object) — low-likelihood boundary, consistent with siblings.
(d) `@param plot.file` lacks the standard "[default NULL]" clause.

## Cleared during verification (checked, not a defect)

- No history append (correct for report mode); dimensions, ploidy, and
  metadata row counts preserved; both datatypes pass.
- Reported median matches an independent `summary()` of the recalculated
  CallRate exactly.
- Plot correctly suppressed at `verbose = 0` (graphics half of VRB5) and
  by `plot.display = FALSE`; `utils.plot.save()` decoupled from display.
- `@author` already carries both Author(s) and Custodian (DOC7
  compliant).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: behaviour vs roxygen on testset.gl / testset.gs; quantile table
  cross-checked against an independent computation — run
- FBM path (DAT6): SKIPPED — no FBM-backed fixture readily available
- Caller grep (for F1/F3/F4 behaviour changes): run in Phase C, results
  recorded in Outcome

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur | gate at verbose >= 1 (console-change consequence approved) |
| 2 | approved | Arthur | |
| 3 | approved | Arthur | silent no-op -> warn + 'loc' consequence approved |
| 4 | approved | Arthur | truly-unaltered-return consequence approved |
| 5 | approved | Arthur | |
| 6 | approved | Arthur | |
| 7 | approved | Arthur | |
| 8 | approved | Arthur | |
| 9 | approved | Arthur | |

F10 notes only, not offered.

## Outcome

- **F1**: applied — both printing blocks gated at `verbose >= 1`.
  Verified: 0 lines at `verbose = 0` for both methods (was 33); default
  output unchanged.
- **F2**: applied — the two `by.pop` example lines removed.
- **F3**: applied — unknown `method` warns and coerces to "loc".
  Verified with `method = "pop"`.
- **F4**: applied — input held before the internal
  `utils.recalc.callrate()` and returned untouched. Verified:
  `expect_identical()` passes even with a deliberately stale CallRate,
  and the printed report still uses the fresh values.
- **F5-F9**: applied (stray `)`, "3rd quartile", `ind.to.list = 0` skips
  the listing, bins doc -> 50, tags reordered per ratified DOC1).
- Caller grep: zero production callers of `gl.report.callrate()` in
  dartR.base or any sibling; only `eval=FALSE` tutorial/README chunks,
  and the dartRstartup tutorial already instructs users not to assign
  the return — F4 aligns behaviour with that guidance. No dartr2shiny.
  NEWS.md entry added.
- Characterization test: 15/15 pass. `verbose = 3` end-to-end clean on
  both datatypes and both methods. `man/gl.report.callrate.Rd`
  regenerated.
- PR: [green-striped-gecko/dartR.base#247](https://github.com/green-striped-gecko/dartR.base/pull/247)

```json
{
  "function": "gl.report.callrate",
  "package": "dartR.base",
  "family": "report",
  "skill_version": "2.0.0",
  "commit": "b04c5bb",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "VRB5", "status": "approved"},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved"},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved"},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS8", "status": "approved"},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "approved"},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved"},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "FS5", "status": "approved"},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved"},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved"},
    {"id": "F10", "severity": "LOW", "confidence": "medium", "rule": "DOC2", "status": "note_only"}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture available"],
  "status": "pr-open",
  "pr": 247
}
```
