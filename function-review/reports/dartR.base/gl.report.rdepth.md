# Review: gl.report.rdepth (dartR.base)

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.report.rdepth
- Datasets: testset.gl (SNP, rdepth metric), testset.gs (SilicoDArT,
  AvgReadDepth metric)
- Family mode: report
- Checks skipped: plot rendering verified only for object construction
  (headless run); FBM path not exercised (not available: no FBM fixture);
  Google Group not searched (not available: no browser session — GitHub
  issues showed no open complaint naming this function).

## Verdicts

- **Standards: FAIL** — the stats block and quantile table print
  unconditionally at verbose = 0 (VRB2), the `verbose == 0 →
  plot.display = FALSE` guard is missing, the "3r quartile" typo, a large
  commented-out dead block, and the header breaks DOC1/DOC2/DOC7.
- **Spec: FAIL** — `plot.file` with `plot.display = FALSE` crashes
  ("object 'p3' not found") because the plots are only built inside the
  display block, so saving is impossible without displaying; and the
  "Retained" counts in the quantile table are NA-inflated when the
  read-depth metric contains NAs. On clean data all printed values match
  independent recomputation; the returned object is the input, unaltered
  and invisible, with no history append (correct for the family).

## Findings

### F1 — plot.file + plot.display=FALSE crashes; saving requires displaying (MEDIUM, confidence: certain)

Rule: spec axis; PLT3 (results/outputs independent of plotting).
Location: `R/gl.report.rdepth.r:108-139,183-195`.

`p1`/`p2` are built only inside `if (plot.display)`, and `p3` only inside
the second `if (plot.display)` block, but the `plot.file` save block
references `p3` unconditionally. Verified:
`gl.report.rdepth(testset.gl, plot.display = FALSE, plot.file = "x")`
errors with "object 'p3' not found".

Proposed change: build the plots unconditionally (cheap — ggplot objects
are lazy); gate only `print(p3)` on plot.display. plot.file then saves the
RDS regardless of display, and the crash disappears.

### F2 — Stats block and quantile table print at verbose = 0 (MEDIUM, confidence: certain)

Rule: VRB2 (ratified gate for report results: verbose >= 1).
Location: `R/gl.report.rdepth.r:142-154,188`.

The summary statistics (`Reporting Read Depth by Locus`, min/quartiles/
mean/max, missing rate) and `print(df)` are ungated. Verified: 33 lines at
verbose = 0.

Proposed change: wrap both in `if (verbose >= 1)`.

### F3 — Missing `verbose == 0 → plot.display <- FALSE` guard (LOW, confidence: certain)

Rule: VRB2/PLT (ratified pattern present in the other reviewed report
functions). Location: `R/gl.report.rdepth.r:64-66`.

At verbose = 0 with the default plot.display = TRUE the plot still
renders. Proposed change: add the standard guard after SET VERBOSITY.

### F4 — Retained counts NA-inflated (LOW, confidence: certain)

Rule: spec axis (DAT NA-handling lens). Location:
`R/gl.report.rdepth.r:158-160`.

`length(rdepth[rdepth >= y])` counts NA comparisons as retained loci
(NA subscripts yield NA elements, which length() counts). Verified: with
10 NAs injected into rdepth, the 0% row reports 255 retained where 245 is
correct. testset.gl/testset.gs have no NAs in the metric, so current
printed output is unaffected — this bites only on datasets with NA
read-depth metrics.

Proposed change: `sum(rdepth >= y, na.rm = TRUE)`.

### F5 — "3r quartile" typo (LOW, confidence: certain)

Rule: STY. Location: `R/gl.report.rdepth.r:150`. Same typo already fixed
in gl.report.callrate (PR #247). Proposed change: "3rd quartile".

### F6 — Header: tag order, verbose wording, author, example naming (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC2, DOC7. Location: `R/gl.report.rdepth.r:1-57`.

- Tag order: `@return` after `@export`; `@details` after `@param`;
  ratified order is description, details, param, return, author, examples,
  export.
- `@param verbose` old wording — replace with ratified canon.
- `@author` Custodian only — add "Author(s): Arthur Georges." (DOC7).
- Examples assign to `df` (`df <- gl.report.rdepth(...)`) — misleading:
  the return is the unaltered genlight object, not the printed table.

Proposed change: rewrite header to the ratified template (docs-only).

### F7 — Structure: dead block, section order (LOW, confidence: certain)

Rules: STY3, FS4. Location: `R/gl.report.rdepth.r:68-80,197-225`.

The 29-line commented-out "SAVE INTERMEDIATES TO TEMPDIR" block is dead
code. SET WORKING DIRECTORY and SET COLOURS precede FLAG SCRIPT START.
Proposed change: delete the dead block; move FLAG SCRIPT START to directly
follow SET VERBOSITY.

## Report notes (other functions / not fixed here)

- Default plot.colors is obtained via
  `gl.select.colors(library="brewer", palette="Blues", select=c(7,5),
  verbose=0)` which returns exactly the documented `#2171B5`/`#6BAED6`
  (verified) — equivalent to the literal assignment the other reviewed
  functions use; left as-is.
- Missing-rate line densifies the genotype matrix (`as.matrix(x)`) — FBM
  densification lens; acceptable for a report; left as-is.

## Coverage

Characterization baseline: `tests/testthat/test-gl.report.rdepth.R` — 10
assertions: unaltered invisible return + no history, quantile-table rows
vs independent recomputation, SilicoDArT path, 33-line verbose-0 output
(baseline), p3 crash with plot.file + display FALSE (baseline), error when
metric absent. All 10 pass on the pre-fix code.

## Approval

Decisions recorded 2026-08-29 via approval boxes — all seven findings
**approved**, including both escalation-gate consequences:

- F1: plot.file + plot.display=FALSE crash becomes a working save —
  approved.
- F4: Retained/Filtered counts become NA-correct (values change only for
  datasets with NA read-depth metrics) — approved.
- F2, F3, F5, F6, F7 — approved.

## Outcome

Applied F1–F7. Verification: all 13 characterization assertions pass
post-fix; every diff from baseline maps to an approved finding (RDS saves
without displaying and the file lands on disk [F1]; verbose = 0 fully
silent, table at verbose >= 1 [F2/F3]; retained = 245 not 255 with 10 NA
metrics [F4]). End-to-end at verbose = 3 on testset.gl and testset.gs;
both return their input identical. Caller grep across dartR.base + 7
siblings: doc cross-references only (gl.filter.rdepth seealso/examples);
no code caller. dartr2shiny: not present in the workspace. NEWS entry
added.

```json
{
  "function": "gl.report.rdepth",
  "package": "dartR.base",
  "family_mode": "report",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "rules": ["spec", "PLT3"], "loc": "R/gl.report.rdepth.r:108-139,183-195", "status": "proposed"},
    {"id": "F2", "severity": "MEDIUM", "rules": ["VRB2"], "loc": "R/gl.report.rdepth.r:142-154,188", "status": "proposed"},
    {"id": "F3", "severity": "LOW", "rules": ["VRB2", "PLT"], "loc": "R/gl.report.rdepth.r:64-66", "status": "proposed"},
    {"id": "F4", "severity": "LOW", "rules": ["spec", "DAT"], "loc": "R/gl.report.rdepth.r:158-160", "status": "proposed"},
    {"id": "F5", "severity": "LOW", "rules": ["STY"], "loc": "R/gl.report.rdepth.r:150", "status": "proposed"},
    {"id": "F6", "severity": "LOW", "rules": ["DOC1", "DOC2", "DOC7"], "loc": "R/gl.report.rdepth.r:1-57", "status": "proposed"},
    {"id": "F7", "severity": "LOW", "rules": ["STY3", "FS4"], "loc": "R/gl.report.rdepth.r:68-80,197-225", "status": "proposed"}
  ],
  "datasets": ["testset.gl", "testset.gs"],
  "baseline_test": "tests/testthat/test-gl.report.rdepth.R",
  "pr": 255
}
```
