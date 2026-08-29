# Review: gl.report.taglength (dartR.base)

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.report.taglength
- Datasets: testset.gl (SNP), testset.gs (SilicoDArT — has TrimmedSequence)
- Family mode: report
- Checks skipped: plot rendering verified only for object construction
  (headless run); FBM path not exercised (not available: no FBM fixture);
  Google Group not searched (not available: no browser session — GitHub
  issues showed no open complaint naming this function).

## Verdicts

- **Standards: FAIL** — stats block and quantile table print
  unconditionally at verbose = 0 (VRB2), missing `verbose == 0 →
  plot.display = FALSE` guard, quartiles mislabelled "quantile" (plus the
  recurring "3r" typo), DOC2 wording stale.
- **Spec: FAIL** — `plot.file` with `plot.display = FALSE` crashes
  ("object 'p3' not found"), and the "Retained" counts are NA-inflated
  when tag lengths are NA — the same two defects just fixed in the rdepth
  pair (PRs #255/#256). On clean data all printed values match independent
  recomputation; the return is already unaltered, invisible, and
  history-free (correct for the family).

## Findings

### F1 — plot.file + plot.display=FALSE crashes (MEDIUM, confidence: certain)

Rule: spec axis; PLT3. Location: `R/gl.report.taglength.r:175-187`.

`p1`/`p2` are built unconditionally but `p3` only inside
`if (plot.display)`; the save block references `p3` unconditionally.
Verified: "object 'p3' not found".

Proposed change: build `p3` unconditionally; gate only `print(p3)`.

### F2 — Stats block and quantile table print at verbose = 0 (MEDIUM, confidence: certain)

Rule: VRB2 (ratified gate: verbose >= 1). Location:
`R/gl.report.taglength.r:132-145,180`.

Verified: 33 lines at verbose = 0. Proposed change: wrap the stats block
and `print(df)` in `if (verbose >= 1)`.

### F3 — Missing `verbose == 0 → plot.display <- FALSE` guard (LOW, confidence: certain)

Rule: VRB2/PLT (ratified pattern). Location:
`R/gl.report.taglength.r:74-76`. Proposed change: add the standard guard.

### F4 — Retained counts NA-inflated (LOW, confidence: certain)

Rule: spec axis (DAT NA-handling lens). Location:
`R/gl.report.taglength.r:150-152`.

`length(nchar.tags[nchar.tags >= y])` counts NA tag lengths as retained
(255 vs the correct 245 with 10 NAs injected — verified). NA tag lengths
arise when TrimmedSequence entries are NA. Proposed change:
`sum(nchar.tags >= y, na.rm = TRUE)`.

### F5 — Quartiles mislabelled, "3r" typo (LOW, confidence: certain)

Rule: STY. Location: `R/gl.report.taglength.r:138,141`.

"1st quantile" / "3r quantile" — these are quartiles (the rdepth sibling
says quartile), and "3r" repeats the typo fixed in PRs #247/#255.
Proposed change: "1st quartile" / "3rd quartile".

### F6 — Header conformance (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC2. Location: `R/gl.report.taglength.r:1-65`.

- `@return` after `@export` (DOC1); DOC2 verbose wording stale.
- The ggplot-themes reference block (lines 46-49) is commented with `# `
  instead of `#' `, so it is dead text that never reaches the rendered
  docs.
- The boxplot title hardcodes "SNP data - Tag Length" although SilicoDArT
  data (with TrimmedSequence, as in testset.gs) is accepted and works —
  title made datatype-aware in passing.
- @author already carries Author(s) + Custodian (compliant).

Proposed change: rewrite header to the ratified template; datatype-aware
plot title.

### F7 — Structure: section order (LOW, confidence: certain)

Rule: FS4. Location: `R/gl.report.taglength.r:74-89`. SET WORKING
DIRECTORY and SET COLOURS precede FLAG SCRIPT START. Proposed change:
move FLAG SCRIPT START to directly follow SET VERBOSITY.

## Report notes (other functions / not fixed here)

- gl.filter.taglength is the matched filter (next in the queue); expect
  the same NA/boundary/plot lenses to apply there.
- Return is already `invisible(x)` — no change needed (contrast with the
  siblings fixed in PRs #251-#254).

## Coverage

Characterization baseline: `tests/testthat/test-gl.report.taglength.R` —
9 assertions: unaltered invisible return + no history, quantile table vs
independent recomputation (at verbose = 1, valid pre- and post-fix),
SilicoDArT path, TrimmedSequence-absent error, 33-line verbose-0 output
(baseline), p3 crash (baseline). All 9 pass on the pre-fix code.

## Approval

Decisions recorded 2026-08-29 via approval boxes — all seven findings
**approved**, including both escalation-gate consequences (F1: crash
becomes a working save; F4: retained counts NA-correct on affected data).

## Outcome

Applied F1–F7. Verification: all 12 characterization assertions pass
post-fix; every diff from baseline maps to an approved finding (RDS saved
without display [F1]; verbose = 0 fully silent [F2/F3]; retained = 245
not 255 with 10 NA TrimmedSequence entries [F4]; quartile labels [F5]).
End-to-end at verbose = 3 on testset.gl and testset.gs; both return their
input identical; the SilicoDArT plot title now says "Fragment P/A data".
Caller grep across dartR.base + 7 siblings: doc cross-references and a
dartRstartup tutorial bare call at default verbosity only — no code
caller. dartr2shiny: not present in the workspace. NEWS entry added.

```json
{
  "function": "gl.report.taglength",
  "package": "dartR.base",
  "family_mode": "report",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "rules": ["spec", "PLT3"], "loc": "R/gl.report.taglength.r:175-187", "status": "proposed"},
    {"id": "F2", "severity": "MEDIUM", "rules": ["VRB2"], "loc": "R/gl.report.taglength.r:132-145,180", "status": "proposed"},
    {"id": "F3", "severity": "LOW", "rules": ["VRB2", "PLT"], "loc": "R/gl.report.taglength.r:74-76", "status": "proposed"},
    {"id": "F4", "severity": "LOW", "rules": ["spec", "DAT"], "loc": "R/gl.report.taglength.r:150-152", "status": "proposed"},
    {"id": "F5", "severity": "LOW", "rules": ["STY"], "loc": "R/gl.report.taglength.r:138,141", "status": "proposed"},
    {"id": "F6", "severity": "LOW", "rules": ["DOC1", "DOC2"], "loc": "R/gl.report.taglength.r:1-65", "status": "proposed"},
    {"id": "F7", "severity": "LOW", "rules": ["FS4"], "loc": "R/gl.report.taglength.r:74-89", "status": "proposed"}
  ],
  "datasets": ["testset.gl", "testset.gs"],
  "baseline_test": "tests/testthat/test-gl.report.taglength.R",
  "pr": null
}
```
