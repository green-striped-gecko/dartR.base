# Review: gl.check.verbosity (dartR.base)

- Family mode: analysis (environment helper; entry-point check called by
  159 files in `R/`)
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: a1440d7 (`dev_luis`, merged with `origin/dev` aa8b0f9)
- Custodian: Bernd Gruber (STY5: changes go through this report)
- Datasets: none required (no genlight input); `testset.gl` used for one
  downstream check through `gl.report.callrate`
- Baseline: `tests/testthat/test-gl.check.verbosity.R` (snapshot captured
  pre-review, 34 expectations pass)

## Verdicts

**Standards: Needs work**: the logic is short and correct for valid
input, but the header is missing `@details`, a complete `@param` and the
house tag order, and the warning text prints with a stray line break.

**Spec: Needs work**: three kinds of invalid `verbose` value crash with
an R-internal error instead of the documented "set to 2" coercion, and
the global option is returned without any check.

The valid path works: an explicit value in 0 to 5 overrides the global
option, `NULL` falls back to `options()$dartR_verbose` and then to 2.

## Findings

**F1 [MEDIUM, confidence: high] — invalid values crash instead of
coercing (spec axis, FS5)**
`R/gl.check.verbosity.r:31`: the check
`is.numeric(x) & x >= 0 & x <= 5` assumes `x` is a single non-missing
number. For `NA_real_` the condition is `NA`, for a length-0 vector it is
`logical(0)`, and for `c(1, 3)` it has length 2. `if` then stops with an
R-internal error. Because every `gl.*` function calls this helper first,
the user sees the error from inside the function they called, with no
mention of `verbose`.
Failure scenario (verified): `gl.report.callrate(testset.gl, verbose =
NA_real_)` stops with "missing value where TRUE/FALSE needed".
`gl.check.verbosity(c(1, 3))` stops with "the condition has length > 1";
`gl.check.verbosity(numeric(0))` stops with "argument is of length zero".
A computed value such as `verbose = mean(v)` over a vector with an `NA`
reaches this path.
Proposed change: guard with `!is.numeric(x) || length(x) != 1 ||
is.na(x) || x < 0 || x > 5`, the same test `gl.set.verbosity` already
uses, so these values warn and return 2 like every other invalid value.

**F2 [LOW, confidence: high] — the global option is not validated
(spec axis, DOC5 (proposed rule))**
`R/gl.check.verbosity.r:28`: `options()$dartR_verbose` is returned
as is. `gl.set.verbosity` validates what it writes, so the gap is only
reachable when the option is set directly with `options()`.
Failure scenario (verified): after `options(dartR_verbose = "loud")`,
`gl.check.verbosity()` returns `"loud"`; downstream `verbose >= 2`
becomes a string comparison (`"loud" >= "2"` is `TRUE`), so the function
runs at an arbitrary level. `options(dartR_verbose = 9)` returns 9.
Proposed change: pass the option value through the same validation as
an explicit argument (warn, return 2).

**F3 [LOW, confidence: high] — warning text contains a line break and
20 spaces (VRB2, STY1)**
`R/gl.check.verbosity.r:36-37`: the message string spans two source
lines, so the printed warning reads
"...in the range \n                    0 to 5, set to 2".
Failure scenario (verified): `gl.check.verbosity(6)` prints the warning
split across two lines, the second indented by 20 spaces.
Proposed change: build the message with `paste0()` or a single-line
string; name the value received, for example
"Warning: Parameter verbose must be a number from 0 to 5 (received 6); set to 2".

**F4 [LOW, confidence: high] — header incomplete and out of order
(DOC1, DOC5 (proposed rule))**
`R/gl.check.verbosity.r:1-18`:
- no `@details`;
- `@param x` does not say what values are valid, or that invalid values
  become 2 with a warning;
- `@title` says "Checks the current global verbosity", but the function
  returns the effective verbosity, and an explicit argument takes
  precedence over the global setting;
- `@return` says "in variable verbose", which describes the caller, not
  the return value;
- tag order is `@examples`, `@author`, `@export`, `@return` instead of
  the house order, and a bare (non-`#'`) blank line sits between
  `@description` and `@param` at line 9. roxygen2 renders the Rd
  correctly today (checked), but the block depends on that tolerance.
Failure scenario: a developer reading `?gl.check.verbosity` cannot tell
that 2.5 is accepted, that `TRUE` becomes 2, or which of the argument
and the option wins.
Proposed change: rewrite the header in house order with the missing
tags. Docs-only; `devtools::document()` in the same change (DOC4).

**F5 [INFO, confidence: high] — non-integer and reserved levels accepted
(VRB1)**
`gl.check.verbosity(2.5)` returns 2.5 and `gl.check.verbosity(4)`
returns 4. Callers gate with `>=`, so 2.5 behaves as level 2 and 4
behaves as level 3. Nothing breaks. No change proposed; F4 documents
the accepted range.

## Proposed changes

1. Guard against `NA`, length-0 and multi-value `verbose` so they warn
   and return 2 instead of stopping (F1).
   **Consequence: calls that currently stop with an R error now run at
   verbosity 2, after a warning.**
2. Validate the global option with the same rule, warning and returning
   2 on an invalid value (F2).
   **Consequence: a session with an invalid `dartR_verbose` option now
   runs at verbosity 2 with a warning at every function call, instead of
   running at an arbitrary level silently.**
3. Print the warning on one line and include the value received (F3).
4. Rewrite the roxygen header in house order: accurate `@title`,
   `@details` on precedence and coercion, complete `@param` and
   `@return` (F4). Docs-only.

## Coverage

- Standards walk: FS, DOC, VRB, DEP, STY — run. FS2–FS4 and FS8–FS9 do
  not apply: this is the helper those rules call, and a start/end banner
  inside it would print on every function call.
- DAT, PLT, DAT6 (FBM): not applicable, no genlight input or plot.
- Spec: behaviour against roxygen, 15 input cases probed with
  `devtools::load_all()` on R 4.4.2 — run.
- Downstream: one caller (`gl.report.callrate`) run with an invalid
  option and an `NA_real_` argument — run.
- Callers of a changed signature: not needed, the signature does not
  change.
- Google Group / GitHub issues: not searched (not available: no browser
  session in this review).

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | warn and use 2 (not stop) |
| 2 | approved | Luis | |
| 3 | approved | Luis | |
| 4 | approved | Luis | |

## Outcome

- Change 1 applied: `NA_real_`, `c(1, 3)` and `numeric(0)` warn and
  return 2; `gl.report.callrate(testset.gl, verbose = NA_real_)` runs.
  Test: `test-gl.check.verbosity.R`, blocks "NA_real_ and wrong-length
  values" and "a caller given verbose = NA_real_".
- Change 2 applied: option `"loud"` and `9` warn and return 2; an explicit
  valid argument still wins silently over an invalid option. Test block
  "an invalid global option".
- Change 3 applied: warning is one line, for example
  "Warning: verbose must be a single number from 0 to 5 (received 6); set
  to 2". The warning names its source (`verbose` or `option
  dartR_verbose`).
- Change 4 applied: header rewritten; `devtools::document()` run;
  `man/gl.check.verbosity.Rd` regenerated.
- Snapshot: 47 expectations pass. Every diff from the baseline is marked
  `[approved diff change N]` in the test file; unchanged: `NULL` fallback,
  argument precedence, silent pass-through of 0 to 5 and 2.5, and the
  coercion of 6, -1, `TRUE`, `"3"`, `NA` to 2.
- Full suite (`devtools::test()`): 9 other files fail
  (`gl.filter.hwe`, `gl.filter.ld`, `gl.fixed.diff`, `gl.read.vcf`,
  `gl.report.allelerich`, `gl.report.basics`, `gl.report.hwe`,
  `gl.report.ld.map`, `utils.heatmap`). The same files fail with the
  original `gl.check.verbosity`, and run file by file both versions give
  13 failures in them, so none come from this change.
- End to end: `gl.filter.callrate(testset.gl, threshold = 0.9, verbose = 3)`
  runs with the new helper.
- NEWS entry added (escalation gate: user-visible behaviour change).
- PR: pending.

```json
{
  "function": "gl.check.verbosity",
  "package": "dartR.base",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "a1440d7",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 4},
    {"id": "F5", "severity": "INFO", "confidence": "high", "rule": "VRB1", "status": "no-change", "change": null}
  ],
  "coverage_skipped": ["DAT/PLT/DAT6: no genlight input or plot", "Google Group: no browser session"],
  "status": "pr-open",
  "pr": null
}
```
