# Review: gl.report.reproducibility (dartR.base)

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.report.reproducibility
- Datasets: testset.gl (SNP, RepAvg), testset.gs (SilicoDArT,
  Reproducibility)
- Family mode: report
- Checks skipped: plot rendering verified only for object construction
  (headless run); FBM path not exercised (not available: no FBM fixture);
  Google Group not searched (not available: no browser session — GitHub
  issues showed no open complaint naming this function).

## Verdicts

- **Standards: FAIL** — stats block and quantile table print
  unconditionally at verbose = 0 (VRB2), missing `verbose == 0 →
  plot.display = FALSE` guard, "3r quartile" typo, 28-line dead commented
  block, header breaks DOC1/DOC2/DOC7.
- **Spec: FAIL** — `plot.file` with `plot.display = FALSE` crashes
  ("object 'p3' not found"), and the "Retained" counts are NA-inflated
  when the repeatability metric contains NAs. The same defect family as
  the rdepth/taglength reports (PRs #255/#257) — this file is a
  line-for-line sibling of gl.report.rdepth. On clean data all printed
  values match independent recomputation; the return is already
  unaltered, invisible, and history-free.

## Findings

### F1 — plot.file + plot.display=FALSE crashes (MEDIUM, confidence: certain)

Rule: spec axis; PLT3. Location: `R/gl.report.reproducibility.r:134-153,197-209`.

`p1`/`p2` are built only inside `if (plot.display)` and `p3` in the second
display block; the save block references `p3` unconditionally. Verified:
"object 'p3' not found". Proposed change: build the plots unconditionally;
gate only `print(p3)`.

### F2 — Stats block and quantile table print at verbose = 0 (MEDIUM, confidence: certain)

Rule: VRB2 (ratified gate: verbose >= 1). Location:
`R/gl.report.reproducibility.r:155-168,202`.

Verified: 33 lines at verbose = 0. Proposed change: wrap in
`if (verbose >= 1)`.

### F3 — Missing `verbose == 0 → plot.display <- FALSE` guard (LOW, confidence: certain)

Rule: VRB2/PLT. Location: `R/gl.report.reproducibility.r:73-75`.
Proposed change: add the standard guard.

### F4 — Retained counts NA-inflated (LOW, confidence: certain)

Rule: spec axis (DAT NA-handling lens). Location:
`R/gl.report.reproducibility.r:172-174`.

`length(repeatability[repeatability >= y])` counts NA metrics as retained
(verified: 255 vs correct 245 with 10 NAs). Proposed change:
`sum(repeatability >= y, na.rm = TRUE)`.

### F5 — "3r quartile" typo (LOW, confidence: certain)

Rule: STY. Location: `R/gl.report.reproducibility.r:164`. Fourth
occurrence of this typo in the campaign (PRs #247/#255/#257).
Proposed change: "3rd quartile".

### F6 — Header conformance and dead block (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC2, DOC7, STY3. Location:
`R/gl.report.reproducibility.r:1-64,211-239`.

- `@return` after `@export`; `@details` after `@param` (DOC1); DOC2
  verbose wording stale; `@author` Custodian only (DOC7).
- 28-line commented-out "SAVE INTERMEDIATES TO TEMPDIR" block is dead
  code (same block removed from gl.report.rdepth in PR #255).
- "repeatbility" typo in @details.

Proposed change: rewrite the header to the ratified template; delete the
dead block.

### F7 — Structure: section order (LOW, confidence: certain)

Rule: FS4. Location: `R/gl.report.reproducibility.r:76-88`. SET WORKING
DIRECTORY and SET COLOURS precede FLAG SCRIPT START. Proposed change:
move FLAG SCRIPT START to directly follow SET VERBOSITY.

## Report notes (other functions / not fixed here)

- gl.filter.reproducibility is the matched filter (next in the queue);
  expect the same NA/plot lenses to apply.
- Boxplot `ylim(c(min(repeatability), 1))` lacks na.rm — with NA metrics
  the lower bound becomes NA (ggplot treats it as unbounded); harmless,
  left as-is.

## Coverage

Characterization baseline:
`tests/testthat/test-gl.report.reproducibility.R` — 8 assertions:
unaltered invisible return + no history, quantile table vs independent
recomputation (verbose = 1), SilicoDArT path, metric-absent error,
33-line verbose-0 output (baseline), p3 crash (baseline). All 8 pass on
the pre-fix code.

## Approval

Decisions recorded 2026-08-29 via approval boxes — all seven findings
**approved**, including both escalation-gate consequences (F1: crash
becomes a working save; F4: retained counts NA-correct on affected data).

## Outcome

Applied F1–F7. Verification: all 11 characterization assertions pass
post-fix; every diff from baseline maps to an approved finding (RDS saved
without display [F1]; verbose = 0 fully silent [F2/F3]; retained = 245
not 255 with 10 NA RepAvg values [F4]; "3rd quartile" [F5]). The boxplot/
histogram limits also gain na.rm = TRUE so plots build on NA-bearing
metrics. End-to-end at verbose = 3 on testset.gl and testset.gs; both
return their input identical. Caller grep across dartR.base + 7 siblings:
docs, seealso, and tutorial references only — no code caller.
dartr2shiny: not present in the workspace. NEWS entry added.

```json
{
  "function": "gl.report.reproducibility",
  "package": "dartR.base",
  "family_mode": "report",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "rules": ["spec", "PLT3"], "loc": "R/gl.report.reproducibility.r:134-153,197-209", "status": "proposed"},
    {"id": "F2", "severity": "MEDIUM", "rules": ["VRB2"], "loc": "R/gl.report.reproducibility.r:155-168,202", "status": "proposed"},
    {"id": "F3", "severity": "LOW", "rules": ["VRB2", "PLT"], "loc": "R/gl.report.reproducibility.r:73-75", "status": "proposed"},
    {"id": "F4", "severity": "LOW", "rules": ["spec", "DAT"], "loc": "R/gl.report.reproducibility.r:172-174", "status": "proposed"},
    {"id": "F5", "severity": "LOW", "rules": ["STY"], "loc": "R/gl.report.reproducibility.r:164", "status": "proposed"},
    {"id": "F6", "severity": "LOW", "rules": ["DOC1", "DOC2", "DOC7", "STY3"], "loc": "R/gl.report.reproducibility.r:1-64,211-239", "status": "proposed"},
    {"id": "F7", "severity": "LOW", "rules": ["FS4"], "loc": "R/gl.report.reproducibility.r:76-88", "status": "proposed"}
  ],
  "datasets": ["testset.gl", "testset.gs"],
  "baseline_test": "tests/testthat/test-gl.report.reproducibility.R",
  "pr": null
}
```
