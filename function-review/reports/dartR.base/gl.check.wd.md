# Review: gl.check.wd (dartR.base)
- Family mode: analysis (utility; accepts a path, not a genlight)
- Date: 2026-09-14
- Reviewer: Claude (Claude Opus 4.8), dartr-function-review v1.0.0
- Package commit: f5e7b72 (upstream/dev; working copy verified identical by `git diff upstream/dev -- R/gl.check.wd.r`)
- Datasets: none required (path-resolution utility)
- Baseline: tests/testthat/test-gl.check.wd.R (new file, snapshot captured pre-review; 16 assertions, all passing)

## Verdicts

**Standards: Needs work** — the resolution logic is correct and the FS
backbone (verbosity, flag start/end, explicit return) is present, but the
fallback warning prints at `verbose = 0` and the path check crashes on
inputs it documents a graceful fallback for.

**Spec: Needs work** — the three documented resolution sources (explicit
`wd`, the `dartR_wd` option, `tempdir()`) all behave as described, and the
returned path is always correct for valid input; but the `@param wd`
default is documented as `tempdir()` while the signature default is `NULL`,
and the "fall back to tempdir on a bad path" contract is only honoured for
a character path that exists-or-not — a non-character or multi-element `wd`
errors instead.

## Blast radius

`gl.check.wd` is called at **87 sites across 86 files in every
dartRverse sub-package** (base, popgen, sim, spatial, captive, sexlinked).
Every internal call passes `verbose = 0` and uses only the return value
(`plot.dir <- gl.check.wd(plot.dir, verbose = 0)` or the `outpath`
equivalent). The signature is therefore effectively frozen: no
finding here changes it, and none changes the returned path for the input
every caller actually passes (a character path or NULL).

## Independent verification (spec axis)

All three sources resolved correctly (test 1): `wd = NULL` with no option →
`tempdir()`; `wd = NULL` with `options(dartR_wd=)` set → the option; an
explicit existing directory → returned unchanged. A non-existent path and
an existing file (not a directory) both fall back to `tempdir()` (tests 2,
3). The function does not create a missing directory (test 6) — it only
tests existence and falls back.

## Findings

**F1 [MEDIUM, confidence: high] — the fallback warning prints at verbose = 0 (VRB5, VRB3)**
`R/gl.check.wd.r:52-58` — when `wd` is a non-existent path the function
prints `cat(warn("Warning: The path to the working directory does not
exist! Set to tempdir().\n"))`, ungated by `verbose`.
Failure scenario: every one of the 87 callers passes `verbose = 0` and
relies on `gl.check.wd` being silent. A user who calls, say,
`gl.report.callrate(x, plot.dir = "typo/path", verbose = 0)` gets a stray
"Warning: The path to the working directory does not exist!" line in the
middle of an otherwise-silent run, attributed to no visible function.
Confirmed: one line captured at `verbose = 0` (test 2).
Proposed change: gate the warning at `verbose >= 1`. This keeps `verbose =
0` fully silent (VRB5) while a user running `gl.check.wd` directly, or any
caller at `verbose >= 1`, still learns their path was rejected (VRB4 —
the fallback changes where files are written).
Visible-output note: for the 87 internal callers (all `verbose = 0`) this
removes a leak and changes nothing else; the returned path is unchanged.

**F2 [MEDIUM, confidence: high] — a non-character wd throws instead of falling back (STY3, DAT5)**
`R/gl.check.wd.r:47` — `if (is.character(wd) & dir.exists(wd))` uses `&`,
which evaluates both operands, so `dir.exists(wd)` runs even when
`is.character(wd)` is already `FALSE`. `dir.exists()` requires a character
argument.
Failure scenario: a caller (or a mis-set `options(dartR_wd = 5)`) passes a
numeric or `NA` `wd` and gets an opaque "invalid filename argument" error
instead of the documented tempdir fallback. Confirmed for `5` and `NA`
(test 4).
Proposed change: short-circuit with `&&` and validate up front, so a
non-character `wd` takes the warn-and-fallback branch.

**F3 [MEDIUM, confidence: high] — a length>1 wd raises a condition-length error (STY3)**
`R/gl.check.wd.r:47` — for a multi-element `wd`, `dir.exists()` returns a
logical vector and the `if()` condition has length > 1, a hard error in
R >= 4.2 ("the condition has length > 1"). Confirmed (test 5).
Proposed change: validate that `wd` is a single character string before the
existence test (folds together with F2 as one input-validation guard).

**F4 [LOW, confidence: high] — documentation is stale and internally inconsistent (DOC1, DOC2, DOC5 (proposed rule))**
`R/gl.check.wd.r:10` — `@param wd ... [default: tempdir()]` but the
signature default is `NULL` (line 25); the actual cascade is NULL → the
`dartR_wd` option → `tempdir()`. `:8` refers to "gl.setwd" (the function is
`gl.set.wd`) and misspells "acccordance"; `:11-13` the `verbose` text is
the pre-DOC2 form ("[default 2, unless specified using gl.set.verbosity]").
Failure scenario: a reader of the manual believes the default is
`tempdir()` and misses that a globally-set `dartR_wd` is the effective
default.
Proposed change: `[default NULL; resolves to the dartR_wd option if set,
otherwise tempdir()]`; fix the two typos and the `gl.set.wd` reference;
adopt the DOC2 standard verbose text.

**F5 [INFO, confidence: high] — dead no-op assignment (STY1)**
`R/gl.check.wd.r:47-49` — the valid-path branch is `wd <- wd` (a no-op).
Harmless; remove for clarity when the surrounding guard is rewritten.

## Proposed changes

1. Gate the non-existent-path warning at `verbose >= 1` (F1).
2. Validate `wd` is a single character string, short-circuiting with `&&`,
   so a non-character or multi-element `wd` takes the warn-and-fallback
   branch instead of erroring (F2, F3).
3. Documentation: correct the `wd` default, the `gl.set.wd` reference and
   the two typos, and adopt the standard verbose text (F4).
4. Remove the `wd <- wd` no-op (F5).

None of these changes the signature or the returned path for a valid
character/NULL `wd`, so no caller among the 87 is affected in its result.

## Coverage

- Standards walk: FS, DOC, VRB, DAT (n/a — no genotypes), DEP (n/a),
  PLT (n/a), STY — run
- Spec: three resolution sources, non-existent path, file-not-dir,
  missing-dir (no create), non-character, length>1, NA — run
- FS4 datatype check: correctly ABSENT (this is a path utility, not a
  genlight function)
- Caller survey (API3): 87 call sites across 86 files in all six
  sub-packages — every one passes `verbose = 0` and consumes only the
  return value; signature frozen
- GitHub issues / Google Group: SKIPPED — not queried this session

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur | |
| 2 | approved | Arthur | |
| 3 | approved | Arthur | |
| 4 | approved | Arthur | |

Cross-package caller survey (API3): 87 call sites across 86 files in all six
sub-packages, every one `gl.check.wd(<path>, verbose = 0)` consuming only the
return value. No signature change and no change to the returned path for a
character/NULL wd, so no caller is affected. All clear.

## Outcome

Changes 1-4 applied on branch review-gl.check.wd (commit df042cf), PR
green-striped-gecko/dartR.base#395.

- Characterization suite green (16 assertions); every diff from the
  pre-review baseline maps to an approved change: silence at verbose 0 with
  a bad path (F1), tempdir fallback for non-character / NA / length>1 wd
  (F2, F3).
- Returned working directory unchanged for the three resolution sources and
  for any valid character / NULL wd.
- End-to-end: gl.report.callrate(x, plot.dir = "no/such/path", verbose = 0)
  now emits 0 lines (was 1) — F1 leak closed at a real call site;
  gl.check.wd(verbose = 3) runs clean.
- Signature unchanged; the 87-caller survey is all-clear.

```json
{
  "function": "gl.check.wd",
  "package": "dartR.base",
  "family": "environment",
  "skill_version": "1.0.0",
  "commit": "f5e7b72",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "STY3", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "STY3", "status": "approved", "change": 2},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 3},
    {"id": "F5", "severity": "INFO", "confidence": "high", "rule": "STY1", "status": "approved", "change": 4}
  ],
  "coverage_skipped": ["GitHub issues not queried", "Google Group not queried"],
  "status": "pr-open",
  "pr": 395
}
```
