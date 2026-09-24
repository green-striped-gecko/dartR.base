# Review: gl.print.history (dartR.base)

- Family mode: report (prints the `@other$history` of an object; read-only)
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 2bd61c5 (`dev_luis`, synced with `origin/dev` ac6b0bd)
- Custodian: Bernd Gruber (STY5: changes go through this report)
- Datasets: `testset.gl` filtered twice with `gl.filter.callrate`
  (3 history entries)
- Baseline: `tests/testthat/test-gl.print.history.R` (snapshot captured
  pre-review, 13 expectations pass)
- Callers: none in dartR.base, sibling packages or the dartr2shiny
  config.

## Verdicts

**Standards: Needs work**: an unguarded Suggests package, unused imports
and dead code, and the header has no `@description`.

**Spec: Needs work**: printing the full history works. Selecting entries
renumbers them, an out-of-range entry prints "NULL", three kinds of
invalid input fail with R-internal errors, and the function returns
nothing.

## Findings

**F1 [MEDIUM, confidence: high] — `knitr` used without a guard (DEP1)**
`R/gl.print.history.r:71`: `knitr::kable()` formats the table, but
`knitr` is in Suggests, not Imports, and there is no
`requireNamespace()` guard.
Failure scenario: on an installation without `knitr` (a minimal server
install, for example), every call at `verbose >= 1` stops with
"there is no package called 'knitr'". Not reproduced here, because
`knitr` is installed; this is standard R behaviour for an uninstalled
namespace.
Proposed change: print the table with base R (`cat()`/`format()`) and
drop the `knitr` call. Change 3 needs this anyway.

**F2 [MEDIUM, confidence: high] — selected entries renumbered; missing
entries print "NULL" (spec axis)**
`R/gl.print.history.r:50, 64`: `nr` is `1:nh` on the subset, not the
entry's position in the history. An index beyond the history returns
`NULL`, which prints as a row.
Failure scenario (verified): `gl.print.history(gl3, history = 3)` labels
the third call as entry 1. `history = c(1, 7)` on a 3-entry history
prints entry 1 and a row "2 | NULL" without a message.
Proposed change: number rows by their original position, and stop with
an error naming the valid range when an index is out of range.

**F3 [LOW, confidence: high] — wrapped calls break the table (spec
axis)**
`R/gl.print.history.r:67-69`: calls longer than 80 characters are
wrapped with a newline inside a markdown table cell, so the row spills
onto a new line outside the table.
Failure scenario (verified): the `gl.read.dart(...)` entry prints as
"| 1 |gl.read.dart(... ind.metafile = metadata," followed by a bare
line "probar = TRUE) |".
Proposed change: print each entry as "number, then the call", with
continuation lines indented under the call.

**F4 [LOW, confidence: high] — invalid inputs fail with R-internal
errors (FS5)**
`R/gl.print.history.r:45-53`: `hist2` is only assigned on two paths.
Failure scenario (verified): `gl.print.history()` and
`gl.print.history(data.frame(a = 1))` stop with "object 'hist2' not
found"; `gl.print.history(gl3, history = gl3@other$history)`, which
matches the documented "link to a history slot", stops with "invalid
subscript type 'list'".
Proposed change: check the inputs first. Accept a genlight (optionally
with numeric `history` indices), or a history list on its own, and stop
with a message saying so otherwise. When `x` is given and `history` is
a list, use the list.

**F5 [LOW, confidence: high] — returns nothing; `verbose = 0` does
nothing (DOC5 (proposed rule), FS10)**
`R/gl.print.history.r:70, 85`: the table prints only at `verbose >= 1`,
and the function returns `NULL` (the value of the final `if`). At
`verbose = 0` the call has no effect at all. `@return` says only that
it prints.
Failure scenario (verified): `h <- gl.print.history(gl3, verbose = 0)`
prints nothing and `h` is `NULL`, so a script cannot inspect the history
table.
Proposed change: return the table (entry number and call) as a data
frame, invisibly, at every verbosity; keep printing at `verbose >= 1`.
**Consequence: the return value changes from `NULL` to an invisible
data frame.**

**F6 [LOW, confidence: high] — unused imports and dead code (STY1,
DEP2)**
`R/gl.print.history.r:30, 62, 74-77`: `@importFrom gridExtra grid.table
ttheme_default` imports functions used only in commented-out code. No
other dartR.base function uses `gridExtra` (the one reference in
`gl.report.hwe` calls `ggtern::grid.arrange`). So `gridExtra` is an
Imports dependency that is never called.
Failure scenario: every dartR.base installation installs `gridExtra`
for nothing; removing the import without removing the dependency makes
`R CMD check` report "Namespace in Imports field not imported from".
Proposed change: remove the import and the commented code, and move
`gridExtra` out of `Imports` in `DESCRIPTION`. dartR.popgen and
dartR.captive declare `gridExtra` themselves.

**F7 [LOW, confidence: high] — header incomplete (DOC1)**
`R/gl.print.history.r:1-33`:
- no `@description` (the Rd has no description section);
- `@param history` says a full history "recreates the identical object
  x", a sentence copied from `gl.play.history`; this function only
  prints;
- "forth" for "fourth";
- `@return` describes printing, not a value;
- `@examples` sits before `@return`; the example is in `\donttest{}` and
  reads a CSV file where `testset.gl` would do and run fast;
- `warning(warn(...))` wraps a crayon string in an R warning; the house
  style is `cat(warn(...))` or `stop(error(...))`, which F4's checks
  replace anyway.
Proposed change: rewrite the header; the example uses `testset.gl` and
runs outside `\donttest{}`.

## Proposed changes

1. Print the table with base R, without `knitr` (F1). With change 3.
2. Number entries by their position in the history; stop with a message
   for out-of-range indices (F2).
   **Consequence: `history = c(1, 7)` on a 3-entry history stops with an
   error instead of printing a NULL row.**
3. Print each entry as number then call, with wrapped lines indented
   (F3). With change 1.
4. Check inputs first; accept a history list alongside `x` (F4).
5. Return the table as an invisible data frame at every verbosity (F5).
   **Consequence: the return value changes from `NULL` to a data frame.**
6. Remove the unused `gridExtra` imports and the dead code; move
   `gridExtra` out of `Imports` in `DESCRIPTION` (F6).
   **Consequence: dartR.base no longer installs `gridExtra`; any user
   code that relied on dartR.base to pull it in must load it itself.**
7. Rewrite the header and use a fast `testset.gl` example (F7).

## Coverage

- Standards walk: FS, DOC, VRB, DEP, STY — run. DAT not applicable (no
  genotypes touched). PLT not applicable (no plot). DAT6: history lives
  in `@other`, the same for FBM-backed objects — not run separately.
- Spec: behaviour against roxygen on 9 input cases — run.
- Report-mode check: the input object is not modified and no history is
  appended — confirmed by reading the code (the function never assigns
  to `x`).
- Missing `knitr`: not reproduced (installed); reasoning from R
  behaviour.
- Callers: none found in dartR.base, siblings or dartr2shiny config.
- Google Group / GitHub issues: not searched (not available: no browser
  session).

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | |
| 2 | approved | Luis | |
| 3 | approved | Luis | |
| 4 | approved | Luis | |
| 5 | approved | Luis | |
| 6 | approved | Luis | |
| 7 | approved | Luis | |

## Outcome

- Changes 1 and 3 applied: base-R output, one entry per call ("1
  gl.read.dart(...)"), wrapped lines indented by two spaces; no
  `knitr` call.
- Change 2 applied: `history = 3` prints as entry 3; `history = c(1, 7)`
  stops with "history must be entry numbers from 1 to 3".
- Change 4 applied: missing or non-genlight `x` stops with "provide a
  genlight object as x, or a history list as history"; a list alongside
  `x` is used.
- Change 5 applied: returns an invisible data frame (`nr`, `history`) at
  every verbosity; empty history warns at `verbose >= 2` and returns 0
  rows.
- Change 6 applied: `importFrom(gridExtra, ...)` gone from `NAMESPACE`;
  `gridExtra` removed from `Imports`. `R CMD check` (no tests, examples,
  vignettes or manual) raises no note or warning about it.
- Change 7 applied: header rewritten; example uses `testset.gl`, outside
  `\donttest{}`, and runs.
- Snapshot: 17 expectations pass; every diff from the baseline is marked
  `[approved diff change N]`. An FBM-backed object prints the same.
- Full suite: same 10 failing files and 21 failures as the run before
  this change (pre-existing; see the gl.fdsim report).
- NEWS entry added (escalation gate: return value and dependency change).
- PR: pending.

```json
{
  "function": "gl.print.history",
  "package": "dartR.base",
  "family": "report",
  "skill_version": "2.0.0",
  "commit": "2bd61c5",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "DEP1", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "spec", "status": "approved", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "spec", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "FS5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "FS10", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DEP2", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 7}
  ],
  "coverage_skipped": ["DAT/PLT: not applicable", "missing-knitr not reproduced", "Google Group: no browser session"],
  "status": "pr-open",
  "pr": null
}
```
