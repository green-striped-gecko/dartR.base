# Review: gl.report.secondaries (dartR.base)

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.report.secondaries
- Datasets: platypus.gl (working dataset — 9 secondaries, TrimmedSequence present),
  testset.gl (no-secondaries case), testset.gs (SilicoDArT rejection)
- Family mode: report
- Checks skipped: plot rendering verified only for object construction, not
  visual output (no display in headless run); Google Group not searched
  (not available: no browser session in this run — GitHub issues in
  green-striped-gecko/dartR.base showed no open complaint naming this function).

## Verdicts

- **Standards: FAIL** — results and a warning print unconditionally at
  verbose = 0 (VRB2), ~330 lines of raw lambda iterations print at default
  verbose = 2, the return is visible, and the header breaks the ratified tag
  order and author format (DOC1/DOC2/DOC7).
- **Spec: FAIL** — the documented no-secondaries path ("Warning: No loci with
  secondaries…") is unreachable: an unused leftover subset line crashes the
  function on any dataset without secondaries. On data *with* secondaries, all
  eight returned values were verified correct against independent
  recomputation.

## Findings

### F1 — Crashes on any dataset with no secondaries (HIGH, confidence: certain)

Rule: spec axis (documented branch unreachable); STY3 (dead code).
Location: `R/gl.report.secondaries.r:216`.

`x.secondaries <- x[, duplicated(b)]` is never used anywhere in the function.
When the dataset has no secondaries, `duplicated(b)` is all-FALSE, the subset
produces zero loci, and dartR's genlight subset method throws
`"Subsetting resulted in zero loci."` — before the function's own graceful
branch at line 363 ("No loci with secondaries, no plot produced" + a
data.frame of zeros/NAs) can run. Verified empirically: `testset.gl`
(255 loci, all unique tags) crashes; so would any post-`gl.filter.secondaries`
object, which is a natural thing to feed a report function.

Proposed change: delete line 216. The graceful branch then works as
documented (verified on testset.gl with the line removed in a scratch copy).
Behaviour change: crash → documented result; escalation-gate class because
user-visible behaviour changes (error becomes a normal return).

### F2 — Results print unconditionally at verbose = 0 (MEDIUM, confidence: certain)

Rule: VRB2 (verbose = 0 fully silent); ratified gate for report results is
`verbose >= 1`.
Location: `R/gl.report.secondaries.r:378-417`.

The entire results block (`Total number of SNP loci scored…` through `Total
Number of invariant sites…`) is bare `cat()` with no verbosity gate. At
verbose = 0 the function emits 10–11 lines. Verified: probe at verbose = 0 on
platypus.gl emitted 11 lines.

Proposed change: wrap the results block in `if (verbose >= 1) { … }` (the
ratified level for report-function results).

### F3 — TrimmedSequence-absent warning prints at verbose = 0 (MEDIUM, confidence: certain)

Rule: VRB2. Location: `R/gl.report.secondaries.r:204-211`.

`cat(warn("The column 'TrimmedSequence' was not found…"))` is ungated.
Verified: with TrimmedSequence removed, 2 warning lines print at verbose = 0.

Proposed change: gate with `if (verbose >= 1)`.

### F4 — ~330 lines of raw lambda iterations at default verbosity (MEDIUM, confidence: certain)

Rule: VRB3 (verbose 2 = progress log, not debug dump); STY3.
Location: `R/gl.report.secondaries.r:282-284`.

Inside the fixed-point loop, `if (verbose >= 2) print(k)` prints the raw
lambda estimate every iteration. On platypus.gl the loop runs ~330 iterations,
so a default-verbosity call emits ~330 bare `[1] 0.0864…` lines before the
"Converged on Lambda" message. This is a debug leftover, not a progress log.

Proposed change: gate the per-iteration print at `verbose >= 5` (full report);
the existing "Converged on Lambda of …" message at verbose >= 2 already
reports the outcome.

### F5 — Return is visible (MEDIUM, confidence: certain)

Rule: VRB5/FS (report family returns invisibly).
Location: `R/gl.report.secondaries.r:434`.

`return(data.frame(…))` auto-prints the 8-row parameter table at the console
on a bare call, duplicating the cat() results block. Verified with
`withVisible()`: visible = TRUE. All documented usage assigns the result
(`n.inv <- gl.report.secondaries(test)`), which is unaffected.

Proposed change: `return(invisible(…))`. Escalation-gate class (console
behaviour changes for bare calls; assigned-use callers unaffected).

### F6 — Header: tag order, verbose wording, author format (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC2, DOC7. Location: `R/gl.report.secondaries.r:1-111`.

- Tag order is `@description, @param, @details, @author, @examples, @seealso,
  @references, @importFrom/@import, @export, @return` — ratified order is
  `@description, @details, @param, @return, @author, @examples, @export` (with
  seealso/references/imports before export).
- `@param verbose` uses the old wording — ratified canon: "2, progress log …
  [default NULL, adopting the global verbosity set by gl.set.verbosity(), or
  2 if no global is set]".
- `@author` has Custodian only — DOC7 requires both; default the missing
  author to match: "Author(s): Arthur Georges. Custodian: Arthur Georges…".
- `@param nsim` says "number of simulations to estimate the mean of the
  Poisson distribution" — no simulation occurs; it is the maximum number of
  iterations of the fixed-point convergence. Reword.
- The itemize list of returned parameters sits inside @details after the
  ggplot theme URLs; it belongs under @return. The last item is labelled
  `k` but the returned data.frame row is named `Lambda` — align.

Proposed change: rewrite the header to the ratified template (docs-only; no
behaviour change).

### F7 — FLAG SCRIPT START sits after wd/colour setup (LOW, confidence: certain)

Rule: FS4 (model-function section order: verbosity → flag start → datatype →
specific checks). Location: `R/gl.report.secondaries.r:127-144`.

SET WORKING DIRECTORY and SET COLOURS run before FLAG SCRIPT START, so the
"Starting gl.report.secondaries" banner appears after any colour warning.
Proposed change: move the FLAG SCRIPT START block up to directly follow SET
VERBOSITY (no behaviour change beyond message ordering).

## Report notes (other functions / not fixed here)

- `utils.check.datatype` on platypus.gl warns about all-NA loci at verbose 2 —
  expected, not a defect here.
- Dead commented-out code elsewhere in the function (lines 225, 255, 179,
  181, 317-318) left in place; only the crashing line 216 is touched by F1.
- `xlim(range = c(-1, 1))` at line 231 works only because ggplot2's `xlim(…)`
  accepts dots; harmless wart, not fixed.

## Coverage

Characterization baseline: `tests/testthat/test-gl.report.secondaries.R` —
17 assertions: return shape/names, all 8 values vs independent recomputation
(platypus.gl), input untouched + no history append, SilicoDArT rejected,
no-secondaries crash (baseline), verbose-0 output present (baseline), visible
return (baseline). The last three flip to post-fix assertions if F1/F2/F5 are
approved. All 17 pass on the pre-fix code.

## Approval

Decisions recorded 2026-08-29 via approval boxes:

- F1 (crash on no-secondaries datasets): **approved**, including the
  behaviour change (error becomes documented data.frame return).
- F2 (results gate at verbose >= 1): **approved**.
- F3 (TrimmedSequence warning gate): **approved**.
- F4 (lambda iteration dump to verbose >= 5): **approved**.
- F5 (invisible return): **rejected** — the visible return is retained by
  the member's decision; bare calls continue to auto-print the data.frame.
- F6 (header rewrite to ratified template): **approved**.
- F7 (flag-start moved to follow verbosity): **approved**.

## Outcome

Applied F1, F2, F3, F4, F6, F7. Verification: all 23 characterization
assertions pass post-fix; every diff from baseline maps to an approved
finding (no-secondaries crash → graceful return [F1]; verbose = 0 fully
silent [F2/F3]; verbose = 2 output 339 lines → 18 [F4]; message order [F7]).
End-to-end at verbose = 3 on platypus.gl reproduces all 8 values exactly;
testset.gl (no secondaries) returns the documented zeros/NA data.frame.
Caller grep across dartR.base + 7 sibling packages: doc/example references
only (gl.filter.secondaries, gl.report.heterozygosity,
gl.report.polyploid_heterozygosity, dartRstartup tutorial) — none affected
by the approved changes. dartr2shiny: not present in the workspace. NEWS
entry added. Note under F2: the previously line-wrapped result strings
("…removed on \n filtering:") were reflowed to single lines while gating —
cosmetic text change within the approved finding.

```json
{
  "function": "gl.report.secondaries",
  "package": "dartR.base",
  "family_mode": "report",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "HIGH", "rules": ["spec", "STY3"], "loc": "R/gl.report.secondaries.r:216", "status": "proposed"},
    {"id": "F2", "severity": "MEDIUM", "rules": ["VRB2"], "loc": "R/gl.report.secondaries.r:378-417", "status": "proposed"},
    {"id": "F3", "severity": "MEDIUM", "rules": ["VRB2"], "loc": "R/gl.report.secondaries.r:204-211", "status": "proposed"},
    {"id": "F4", "severity": "MEDIUM", "rules": ["VRB3", "STY3"], "loc": "R/gl.report.secondaries.r:282-284", "status": "proposed"},
    {"id": "F5", "severity": "MEDIUM", "rules": ["VRB5"], "loc": "R/gl.report.secondaries.r:434", "status": "proposed"},
    {"id": "F6", "severity": "LOW", "rules": ["DOC1", "DOC2", "DOC7"], "loc": "R/gl.report.secondaries.r:1-111", "status": "proposed"},
    {"id": "F7", "severity": "LOW", "rules": ["FS4"], "loc": "R/gl.report.secondaries.r:127-144", "status": "proposed"}
  ],
  "datasets": ["platypus.gl", "testset.gl", "testset.gs"],
  "baseline_test": "tests/testthat/test-gl.report.secondaries.R",
  "pr": 253
}
```
