# Review: gl.report.excess.het + gl.filter.excess.het (dartR.base)

Matched pair reviewed together. Both are deprecated wrappers: since
0187142 they call `gl.report.hwe()` / `gl.filter.hwe()` with the settings
of the Robledo-Ruiz et al. (2023) workflow. Each finding names the file it
applies to.

## Provenance

- Family mode: report (`gl.report.excess.het`), modify
  (`gl.filter.excess.het`)
- Date: 2026-09-24
- Reviewer: Claude (Opus 5.5, claude-opus-5-5, Claude Code),
  dartr-function-review v2.0.0
- Package commit: f9f1be8 (origin/dev), worktree branch `review-excess-het`
- Datasets: LBP, testset.gl, testset.gs (dartR.data)
- Baseline: `tests/testthat/test-gl.report.excess.het.R` (7 expectations)
  and `tests/testthat/test-gl.filter.excess.het.R` (6), new files; all pass
  at the reviewed state. Existing stub tests in `test-gl.report.hwe.R` and
  `test-gl.filter.hwe.R` also cover the wrappers.
- Release state: the wrapper form is on `dev` only. `main`/CRAN still ship
  the original stand-alone implementations.
- Callers: none in the sibling packages. dartr2shiny wires
  `gl.filter.excess.het` (config/functions.csv, input_generator).

## Verdicts

**Standards: Needs work** — the wrappers are minimal and correct in what
they call; the filter's history entry cannot be replayed and loses the
user's call, and the report silently discards four plot arguments.
**Spec: Needs work** — both wrappers reproduce the published result (6 loci
on LBP, matching the original); the deprecation message that tells users
how to migrate gives a call that removes 282 loci instead of 6.

## Findings

**F1 [HIGH, confidence: high] — migration advice gives a different result (API1 proposed rule, DOC5 proposed rule)**
`R/gl.filter.excess.het.r:6-9,74-82`, `R/gl.report.excess.het.r:6-9,81-89`
— the `.Deprecated()` message and `@description` give
`gl.filter.hwe(direction = 'excess', test.type = 'ChiSquare',
mult.comp.adj = TRUE, mult.comp.adj.method = 'fdr')`, which omits
`min.hobs = 0.5` and `cc.val`. `min.hobs = 0.5` is the screen that confines
the FDR adjustment to high-heterozygosity loci. Without `cc.val`, the
Yates correction defaults to on (0.5), whatever `Yates` was.
Failure scenario: on LBP the wrapper removes 6 loci. The call printed in
the deprecation warning removes 282, as does the equivalent
`gl.report.hwe()` call. Adding `cc.val = 0, min.hobs = 0.5` gives 6 again.
A user who migrates as told loses 276 extra loci without being told.
Proposed change: give the full equivalent call in both messages and both
`@description`s, including `min.hobs = 0.5` and a `cc.val` (`cc_val` for
the report) equal to the value actually used for this call's `Yates`.

**F2 [MEDIUM, confidence: high] — history entry cannot be replayed (FS8)**
`R/gl.filter.excess.het.r:84-106` — the history keeps the entry that
`gl.filter.hwe()` appends:
`gl.filter.hwe(x = x, ..., cc.val = ifelse(Yates, 0.5, 0), ..., verbose = verbose)`.
`Yates` and `verbose` exist only inside the wrapper, and the call the user
made (`gl.filter.excess.het(LBP, Yates = TRUE)`) is not recorded.
Failure scenario: re-running the history outside the wrapper fails with
"object 'Yates' not found"; the history also misrepresents which function
the user ran.
Proposed change: replace the last history entry with the wrapper's own
`match.call()`. The number of history entries does not change.

**F3 [LOW, confidence: high] — plot arguments accepted and discarded silently (API2 proposed rule)**
`R/gl.report.excess.het.r:68-71` — `plot.theme`, `plot.colors`, `plot.file`
and `plot.dir` are kept for backward compatibility and ignored. This is
documented, but nothing is printed at run time.
Failure scenario: `gl.report.excess.het(LBP, plot.file = "p", plot.dir = d)`
writes no file and gives no message; the user finds out only when the file
is missing.
Proposed change: at `verbose >= 1`, warn once naming the non-NULL ignored
arguments.

**F4 [LOW, confidence: high] — author block and ASCII (DOC7, DOC6 proposed rules)**
`R/gl.filter.excess.het.r:40`, both files lines 6-9 and 44 — the filter's
`@author` reads "Author(s): ... (Custodian: Ching Ching Lau)" rather than
the labelled `Custodian:` form; both files use "Jesús" and an em dash in
roxygen. `gl.report.excess.het` is tagged `@family unmatched report`
although `gl.filter.excess.het` is its matched filter.
Failure scenario: an inconsistent author block, and a PDF manual that
depends on non-ASCII handling across platforms.
Proposed change: labelled `Author(s): ... Custodian: ...` form, ASCII
("Jesus", "--"), and `@family matched report`. Docs only.

Notes (outside these files, not proposed here):
- At `verbose = 0`, one LBP run raises about 940 R warnings ("Expected
  counts below 5: chi-square approximation may be incorrect") from
  `HardyWeinberg::HWChisq` inside `gl.filter.hwe()`/`gl.report.hwe()`.
  `gl.report.hwe()` also prints "Starting gl.colors" at `verbose = 0`
  because its default `plot_colors = gl.colors("2c")` runs at the global
  verbosity. Both belong to the hwe pair (reviewed, done, PR #271) and
  should be raised with its custodian.
- dartr2shiny still generates a module for `gl.filter.excess.het`, so the
  deprecation warning reaches GUI users. Before the wrappers are removed,
  that module must be retargeted to `gl.filter.hwe` (API3).

## Proposed changes

1. Deprecation messages and `@description` give the full equivalent call,
   including `min.hobs = 0.5` and the `cc.val`/`cc_val` matching `Yates`
   (F1). Wrapper output unchanged.
2. The filter records its own call in the history in place of the inner
   `gl.filter.hwe()` call (F2).
   **Consequence: the last history entry changes from `gl.filter.hwe(...)`
   to `gl.filter.excess.het(...)`; the entry count is unchanged.**
3. The report warns when ignored plot arguments are supplied (F3).
4. Author block, ASCII and `@family` fixes (F4). Docs only.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run on both files
- Spec: wrapper results against the original's documented LBP result
  (6 loci, Yates TRUE and FALSE) — run, matches; migration advice against
  the wrapper — run, does not match (F1)
- Report leaves input untouched (FS8) — run, passes
- SilicoDArT input — run; stops with the hwe pair's "Cannot calculate HWE
  from fragment presence/absence data"
- Empty and non-empty results — run (testset.gl: 3 loci; LBP: 6)
- DAT2 loc.metrics sync after filtering — run, passes
- FBM path (DAT6): SKIPPED — the wrappers only delegate to the hwe pair,
  whose review covers it
- Independent numerical check: SKIPPED — the numbers come from
  `gl.report.hwe()`/`gl.filter.hwe()`, reviewed separately; the wrappers
  match the original's published result
- Google Group / GitHub issues: not searched (not available: no browser
  session in this run)

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | |
| 2 | approved | Luis | consequence (last history entry changes) approved explicitly |
| 3 | approved | Luis | |
| 4 | approved | Luis | |

## Outcome

- Change 1 (F1): new tests parse the call from each deprecation warning,
  run it, and compare: same 6 LBP loci as the wrapper for `Yates = TRUE`
  and `FALSE`, in both functions.
- Change 2 (F2): the baseline expectation "history records
  gl.filter.hwe" was flipped (approved diff). A new test evaluates the
  recorded entry and gets the same loci; the entry count is unchanged.
- Change 3 (F3): new test; `plot.file`/`plot.dir` are named in a warning,
  with no warning when none are supplied.
- Change 4 (F4): docs; both R files are ASCII only; `man/` regenerated.
- Baseline: removals and report loci unchanged (LBP 6, testset.gl 3/252).
  Report tests 11 expectations, filter tests 10, all pass. Both functions
  run end to end at `verbose = 3`.
- Pre-existing, unrelated: `test-gl.report.hwe.R:83` (245 vs 196) and
  `test-gl.filter.hwe.R:29` (249 vs 251) fail on clean origin/dev
  f9f1be8 as well (HardyWeinberg 1.7.9, dartR.data 1.2.5). Not caused by
  this change; to raise with the hwe pair's custodian.
- `devtools::document()` also rewrote unrelated Rd files (cross-link
  drift, PR #421); reverted.
- PR: pending.

```json
{
  "function": "gl.report.excess.het",
  "companion": "gl.filter.excess.het",
  "package": "dartR.base",
  "family": "report+modify",
  "skill_version": "2.0.0",
  "commit": "f9f1be8",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "API1", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "FS8", "status": "approved", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "API2", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "approved", "change": 4}
  ],
  "coverage_skipped": ["DAT6: delegated to hwe pair", "numerical check: delegated to hwe pair", "forum/issues: no browser session"],
  "status": "awaiting-approval",
  "pr": null
}
```
