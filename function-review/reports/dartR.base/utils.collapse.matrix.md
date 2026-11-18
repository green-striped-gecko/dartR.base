# Review: utils.collapse.matrix (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.collapse.matrix. Reviewed in the
  infrastructure wave (seven files, one approval round, per-function
  PRs), under the standing member directive that utility functions
  are not for end users.
- Datasets: testset.gl subsets; constructed matrices; ggplot fixtures;
  composed-command inspection (no PLINK binary required).
- Family mode: analysis/infrastructure utility.
- Checks skipped: Google Group not searched (not available: no
  browser session); dartr2shiny not present in the local workspace.

## Verdicts

- **Standards: FAIL** — two cat(error()) + bare stop() splits (empty
  fatal messages); an indented @details tag; @return says "distances
  between individuals" for a function whose whole purpose is
  populations; the dartR-coercion warning is gated at verbose > 2
  (nonstandard).
- **Spec: FAIL (diagonal only)** — off-diagonal population means
  verified exact against manual recomputation, dist-in/dist-out
  preserved, silent at verbose 0; but the within-population diagonal
  includes the zero self-distances, deflating it (1.6195 including vs
  1.7995 excluding on the fixture). The diagonal is only visible in
  the matrix output (as.dist drops it), and gl.dist.pop consumes the
  dist form, so exposure is limited to direct matrix use.

## Findings

### I9 — Diagonal excludes self-distances (MEDIUM) [escalation: matrix diagonal changes]

Within-population means computed over distinct pairs only.

### I10 — Stops, gates, docs (LOW)

stop(error(...)) x2; coercion warning gated at >= 2; @return
corrected; @details indentation fixed; @keywords internal (STAYS
exported — called by gl.dist.pop and captive/popgen functions).

## Coverage

test-utils.collapse.matrix.R — 4 assertions: off-diagonal anchor,
diagonal self-inclusion (baseline), verbose-0 silence + dist class
(anchor). All pass pre-fix.

## Approval

All findings approved via the approval boxes (2026-09-01).

## Outcome

All findings applied (I9 diagonal over distinct pairs, I10 stops/gates/docs). Suite: 4/4; off-diagonal anchors unchanged. PR #325.

```json
{"function": "utils.collapse.matrix", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [{"id": "I9", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/utils.collapse.matrix.r diagonal", "status": "applied"},
  {"id": "I10", "severity": "LOW", "rules": ["STY", "DOC"], "loc": "R/utils.collapse.matrix.r", "status": "applied"}],
 "datasets": ["testset.gl", "constructed"],
 "baseline_test": "tests/testthat/test-utils.collapse.matrix.R",
 "pr": 325}
```
