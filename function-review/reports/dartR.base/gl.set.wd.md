# Review: gl.set.wd (dartR.base)
- Family mode: analysis (utility; sets a global option from a path)
- Date: 2026-09-14
- Reviewer: Claude (Claude Opus 4.8), dartr-function-review v1.0.0
- Package commit: f02bd34 (upstream/dev; working copy verified identical by `git diff upstream/dev -- R/gl.set.wd.r`)
- Datasets: none required
- Baseline: tests/testthat/test-gl.set.wd.R (new file, snapshot captured pre-review; 14 assertions, all passing)

## Verdicts

**Standards: Needs work** — the FS backbone is present, but the path guard
crashes on non-character / multi-element input and the success message is
printed regardless of whether anything was set.

**Spec: Rework** — on an invalid path the function does not set the global
working directory, yet it returns that path and prints "Global working
directory set to <path>". The core contract of a setter — that it either
sets the value or tells you it did not — is broken. This is not a bounded
fix to one line; the invalid-path behaviour has to be redefined.

## Blast radius

`gl.set.wd` has **no internal callers** anywhere in the dartRverse — it is
a user-facing setter invoked directly from scripts (and taught in the eBook,
Ch 1 p18 and the save workflow p48). It writes `options(dartR_wd = )`, which
`gl.check.wd` reads; the round-trip works for a valid path (test 7). So no
sibling function is affected by any change here; the contract that matters
is the one presented to the user.

## Independent verification (spec axis)

For a valid directory the option is set and the path returned (test 1); the
default `gl.set.wd()` sets the option to a fresh `tempdir()` (test 2); and
`gl.check.wd()` then reads back exactly what was set (test 7). The failure
is entirely on the invalid-path branch, verified below.

## Findings

**F1 [HIGH, confidence: high] — an invalid path is silently ignored while success is reported (DOC5 (proposed rule); recurring silent-failure class)**
`R/gl.set.wd.r:36-42` — the guard sets `options(dartR_wd = wd)` only when
`wd` is an existing directory, but the success message
(`"  Global working directory set to <wd>"`, line 42) and the
`return(wd)` (line 48) run unconditionally.
Failure scenario: `gl.set.wd("C:/typo_dir")` prints "Global working
directory set to C:/typo_dir", returns that path, and leaves the real
`dartR_wd` option unchanged (still its prior value or `tempdir()`). The
user believes every subsequent plot/output goes to their chosen directory;
in fact it goes elsewhere, and there is no warning. Confirmed for a
non-existent path (test 3) and for an existing file that is not a directory
(test 4): the option stays at its prior value while the message claims it
was set.
Proposed change: make the invalid-path case honest. The custodian chooses
the contract (Phase B question) — the recommended option is to
`stop(error(...))` naming the missing directory, since a setter invoked
deliberately with a specific path should fail loudly rather than silently
redirect the user's explicit choice.

**F2 [MEDIUM, confidence: high] — a non-character or NULL wd throws an opaque error (STY3)**
`R/gl.set.wd.r:36-38` — `!is.null(wd) & is.character(wd) & dir.exists(wd)`
uses `&`, which evaluates every operand, so `dir.exists(wd)` runs on a NULL
or numeric `wd` and raises "invalid filename argument".
Failure scenario: `gl.set.wd(NULL)` or a programmatic numeric path errors
uninformatively instead of being rejected clearly. Confirmed for `NULL`
and `5` (test 5).
Proposed change: short-circuit with `&&` and validate `wd` is a single
character string; a non-character / NULL `wd` then takes the invalid-path
branch (F1) with a clear message.

**F3 [MEDIUM, confidence: high] — a length>1 wd raises a condition-length error (STY3)**
`R/gl.set.wd.r:36-38` — for a multi-element `wd`, `dir.exists()` returns a
logical vector and the `if()` condition has length > 1 (a hard error in
R >= 4.2). Confirmed (test 6). Folds into the F2 guard rewrite (single
character-string validation).

**F4 [LOW, confidence: high] — documentation gaps (DOC1, DOC2)**
`R/gl.set.wd.r:22` — `@return "path the the working directory"` (typo:
"path the the"); `:10-12` the `verbose` text is the pre-DOC2 form; `:9`
`@param wd` does not state the default (`tempdir()`); "explicitely" is
misspelled twice (lines 7, 9); there is no `@details`.
Proposed change: adopt the DOC2 verbose text, fix "path to the", state the
default, correct the typos.

## Proposed changes

1. Redefine the invalid-path behaviour so the function no longer reports
   success when nothing was set, and validate `wd` as a single character
   string so non-character / NULL / multi-element input is handled by the
   same branch rather than crashing (F1, F2, F3). The specific contract
   (error vs warn-and-fallback-to-tempdir vs warn-and-leave-unchanged) is
   the Phase B question.
   **Consequence: a call with a bad path changes behaviour — currently it
   returns the path and prints a false success message (string paths) or
   throws an opaque low-level error (non-character / vector); after the fix
   it produces one clear, truthful outcome.**
2. Documentation: DOC2 verbose text, `@return` typo, state the default,
   fix "explicitely" (F4).

## Coverage

- Standards walk: FS, DOC, VRB, DAT (n/a), DEP (n/a), PLT (n/a), STY — run
- Spec: valid dir, default (tempdir), non-existent path, file-not-dir,
  non-character, NULL, length>1, and the gl.set.wd -> gl.check.wd
  round-trip — run
- Caller survey (API3): no internal callers in any sub-package; user-facing
  setter only
- GitHub issues / Google Group: SKIPPED — not queried this session

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur | contract: **error loudly** on an invalid path (stop() naming the missing directory; global wd left unchanged). Consequence approved. |
| 2 | approved | Arthur | |

## Outcome

(pending Phase C)

```json
{
  "function": "gl.set.wd",
  "package": "dartR.base",
  "family": "environment",
  "skill_version": "1.0.0",
  "commit": "f02bd34",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "proposed", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "STY3", "status": "proposed", "change": 1},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "STY3", "status": "proposed", "change": 1},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "proposed", "change": 2}
  ],
  "coverage_skipped": ["GitHub issues not queried", "Google Group not queried"],
  "status": "awaiting-approval",
  "pr": null
}
```
