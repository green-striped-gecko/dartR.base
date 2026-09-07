# Review: gl.report.overshoot (dartR.base)

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.report.overshoot
- Datasets: testset.gl (SNP — carries 21 genuine overshoot loci), crafted
  variant (SnpPosition forced out of range), testset.gs (SilicoDArT
  rejection)
- Family mode: report
- Checks skipped: FBM path not exercised (not available: no FBM fixture);
  Google Group not searched (not available: no browser session — GitHub
  issues showed no open complaint naming this function).

## Verdicts

- **Standards: FAIL** — the results (count + locus listing, or the
  no-overshoot message) print unconditionally at verbose = 0 (VRB2);
  header breaks DOC1/DOC2/DOC7 and the @details falsely claim ggplots and
  tabulations are saved to the tempdir (this function produces no plots
  and saves nothing — stale copied text).
- **Spec: PASS (logic)** — the overshoot test `(SnpPosition + 1) >
  nchar(TrimmedSequence)` matches independent recomputation exactly
  (21/21 on testset.gl; 24 with 3 crafted extras, correctly named); the
  return is unaltered, invisible, history-free; SilicoDArT and
  missing-metric inputs error properly. The FAIL elements are
  presentation-layer only, so the overall spec verdict on behaviour vs
  docs is dragged down solely by the false tempdir claim.

## Findings

### F1 — Results print unconditionally at verbose = 0 (MEDIUM, confidence: certain)

Rule: VRB2 (ratified gate: verbose >= 1). Location:
`R/gl.report.overshoot.r:84-92`.

Both branches (count + listing, and "There were no loci…") are bare
`cat()`. Verified: 2 lines at verbose = 0. Proposed change: gate the
results block at `verbose >= 1`.

### F2 — Wasteful subset and malformed listing (LOW, confidence: certain)

Rule: STY; FBM densification lens. Location:
`R/gl.report.overshoot.r:82-89`.

`xx <- x[, os]` constructs a full genlight subset only to call `nLoc(xx)`
and `locNames(xx)` — `length(os)` and `locNames(x)[os]` suffice and avoid
densifying FBM-backed data. The listing `cat(paste0(locNames(xx), sep =
","))` passes `sep` into `paste0`'s dots (paste0 has no sep argument), so
each name gets a comma glued on and cat joins them with spaces, leaving a
stray trailing comma (visible in output). Proposed change: drop the
subset; list with `paste(locNames(x)[os], collapse = ", ")`.

### F3 — Header: false tempdir claim, tag order, wording, author (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC2, DOC7, DOC5 (docs must not promise what the function
does not do). Location: `R/gl.report.overshoot.r:1-34`.

- @details: "Resultant ggplot(s) and the tabulation(s) are saved to the
  session's temporary directory" — no plots exist and nothing is saved;
  remove.
- Tag order: `@details` after `@param`, `@return` after `@export` (DOC1).
- `@param verbose` wording near-canon but missing the "adopting the
  global verbosity … or 2 if no global is set" clause (DOC2).
- `@author` Custodian only (DOC7).

Proposed change: rewrite the header to the ratified template.

## Report notes (other functions / not fixed here)

- testset.gl genuinely contains 21 overshoot loci — useful as a live
  fixture for the matched filter's review (gl.filter.overshoot, next).
- No plot parameters exist on this report — appropriate for its output;
  no plot.file/plot.display findings apply.

## Coverage

Characterization baseline: `tests/testthat/test-gl.report.overshoot.R` —
11 assertions: count vs independent recomputation (21), crafted detection
(24, names listed), unaltered invisible return + no history, SilicoDArT
rejected, both missing-metric errors, 2-line verbose-0 output (baseline).
All 11 pass on the pre-fix code.

## Approval

Decisions recorded 2026-08-30 via approval boxes — all three findings
**approved** (F1 verbosity gate, F2 subset + listing, F3 header rewrite).

## Outcome

Applied F1–F3. Verification: all 13 characterization assertions pass
post-fix; every diff from baseline maps to an approved finding (verbose =
0 fully silent [F1]; listing comma-separated with no trailing comma and
no genlight subset [F2]). End-to-end at verbose = 3 reproduces the count
of 21 and returns the input identical. Caller grep across dartR.base + 7
siblings: tutorial bare calls at default verbosity only — no code caller.
dartr2shiny: not present in the workspace. NEWS entry added.

```json
{
  "function": "gl.report.overshoot",
  "package": "dartR.base",
  "family_mode": "report",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "rules": ["VRB2"], "loc": "R/gl.report.overshoot.r:84-92", "status": "proposed"},
    {"id": "F2", "severity": "LOW", "rules": ["STY"], "loc": "R/gl.report.overshoot.r:82-89", "status": "proposed"},
    {"id": "F3", "severity": "LOW", "rules": ["DOC1", "DOC2", "DOC5", "DOC7"], "loc": "R/gl.report.overshoot.r:1-34", "status": "proposed"}
  ],
  "datasets": ["testset.gl", "testset.gs"],
  "baseline_test": "tests/testthat/test-gl.report.overshoot.R",
  "pr": 261
}
```
