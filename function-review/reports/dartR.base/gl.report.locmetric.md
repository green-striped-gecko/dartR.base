# Review: gl.report.locmetric (dartR.base)

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.report.locmetric
- Datasets: testset.gl (SNP, metric = SnpPosition), testset.gs
  (SilicoDArT, metric = AvgReadDepth)
- Family mode: report
- Checks skipped: plot rendering verified only for object construction
  (headless run); FBM path not exercised (not available: no FBM fixture);
  Google Group not searched (not available: no browser session — GitHub
  issues showed no open complaint naming this function).

## Verdicts

- **Standards: FAIL** — stats block and quantile table print
  unconditionally at verbose = 0 (VRB2); quartiles mislabelled "quantile"
  with the recurring "3r" typo; header breaks DOC1/DOC2/DOC7. Otherwise
  notably compliant: the verbose-0 plot guard, unconditional plot build,
  working plot.file-without-display, and invisible return are all already
  present.
- **Spec: FAIL** — the stats lines are garbled: `summary()` is applied to
  the one-column data.frame instead of the vector, so every line prints a
  doubled label ("Minimum      :  Min.   : 5.0"); and the "Retained"
  counts are NA-inflated when the chosen metric contains NAs. Quantile
  table values on clean data match independent recomputation; the return
  is unaltered, invisible, history-free; both datatypes and all error
  paths work.

## Findings

### F1 — Stats block and quantile table print at verbose = 0 (MEDIUM, confidence: certain)

Rule: VRB2 (ratified gate: verbose >= 1). Location:
`R/gl.report.locmetric.r:163-173,207`.

Verified: 32 lines at verbose = 0. Proposed change: wrap in
`if (verbose >= 1)`.

### F2 — Stats lines carry doubled labels (LOW, confidence: certain)

Rule: spec axis (output quality). Location:
`R/gl.report.locmetric.r:164`.

`stats <- summary(metric_df)` summarises the data.frame, yielding
formatted strings ("Min.   : 5.0   "); each cat line then prepends its
own label, producing "    Minimum      :  Min.   : 5.0   " (verified).
Proposed change: `summary(metric_df$field)` — numeric values print
cleanly, matching the other report functions.

### F3 — Retained counts NA-inflated (LOW, confidence: certain)

Rule: spec axis (DAT NA lens). Location:
`R/gl.report.locmetric.r:178-180`.

`length(field[field >= y])` counts NA metric values as retained
(verified: 255 vs correct 245 with 10 NAs). User-supplied custom metrics
make NAs particularly plausible here. Proposed change:
`sum(field >= y, na.rm = TRUE)`.

### F4 — Quartiles mislabelled, "3r" typo (LOW, confidence: certain)

Rule: STY. Location: `R/gl.report.locmetric.r:169,172`.

"1st quantile" / "3r quantile" — quartiles; fifth occurrence of the "3r"
typo in the campaign. Proposed change: "1st quartile" / "3rd quartile".

### F5 — Header conformance, section order, ploidy-based title (LOW, docs-only + tidy, confidence: certain)

Rules: DOC1, DOC2, DOC7, FS4. Location:
`R/gl.report.locmetric.r:1-93,103-124,144-148`.

- Tag order: `@details` after `@param`, `@return` after `@export` (DOC1);
  DOC2 verbose wording stale.
- `@author` "Luis Mijangos (Post to …)" — no Author(s)/Custodian labels;
  per DOC7 default the missing role to match: "Author(s): Luis Mijangos.
  Custodian: Luis Mijangos".
- SET WORKING DIRECTORY / SET COLOURS precede FLAG SCRIPT START (FS4).
- The plot title keys on `all(x@ploidy == 2)` although `datatype` is
  already computed — switched to datatype in passing.

Proposed change: rewrite header to the ratified template; reorder
sections; datatype-keyed title.

## Report notes (other functions / not fixed here)

- This file already embodies several campaign fixes (plot guard,
  unconditional plot build, invisible return) — a useful contrast with
  its older siblings.
- gl.filter.locmetric is the matched filter (next in the queue).

## Coverage

Characterization baseline: `tests/testthat/test-gl.report.locmetric.R` —
10 assertions: unaltered invisible return + no history, quantile table vs
independent recomputation, SilicoDArT path, both error paths, plot.file
save without display (already working), 32-line verbose-0 output
(baseline), doubled-label stats formatting (baseline). All 10 pass on the
pre-fix code.

## Approval

Decisions recorded 2026-08-30 via approval boxes — all five findings
**approved**, including the escalation-gate consequence (F3: retained
counts NA-correct on affected metrics).

## Outcome

Applied F1–F5. Verification: all 14 characterization assertions pass
post-fix; every diff from baseline maps to an approved finding (verbose =
0 fully silent [F1]; stats lines single-labelled numeric [F2]; retained =
245 not 255 with 10 NA metric values [F3]; quartile labels [F4]).
End-to-end at verbose = 3 on testset.gl (SnpPosition) and testset.gs
(AvgReadDepth); both return their input identical. Caller grep across
dartR.base + 7 siblings: no references at all — not even docs (noted:
gl.filter.locmetric's header does not link its report; to be addressed in
that function's own review). dartr2shiny: not present in the workspace.
NEWS entry added.

```json
{
  "function": "gl.report.locmetric",
  "package": "dartR.base",
  "family_mode": "report",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "rules": ["VRB2"], "loc": "R/gl.report.locmetric.r:163-173,207", "status": "proposed"},
    {"id": "F2", "severity": "LOW", "rules": ["spec"], "loc": "R/gl.report.locmetric.r:164", "status": "proposed"},
    {"id": "F3", "severity": "LOW", "rules": ["spec", "DAT"], "loc": "R/gl.report.locmetric.r:178-180", "status": "proposed"},
    {"id": "F4", "severity": "LOW", "rules": ["STY"], "loc": "R/gl.report.locmetric.r:169,172", "status": "proposed"},
    {"id": "F5", "severity": "LOW", "rules": ["DOC1", "DOC2", "DOC7", "FS4"], "loc": "R/gl.report.locmetric.r:1-93", "status": "proposed"}
  ],
  "datasets": ["testset.gl", "testset.gs"],
  "baseline_test": "tests/testthat/test-gl.report.locmetric.R",
  "pr": 263
}
```
