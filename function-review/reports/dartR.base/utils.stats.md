# Review: utils.stats (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.stats. Reviewed in the
  infrastructure wave (seven files, one approval round, per-function
  PRs), under the standing member directive that utility functions
  are not for end users.
- Datasets: testset.gl subsets; constructed matrices; ggplot fixtures;
  composed-command inspection (no PLINK binary required).
- Family mode: analysis/infrastructure utility.
- Checks skipped: Google Group not searched (not available: no
  browser session); dartr2shiny not present in the local workspace.

## Verdicts

- **Standards: FAIL** — a four-line file defining top-level std.error
  with no documentation of any kind.
- **Spec: PASS** — the computation matches the standard-error
  definition (verified); three in-package callers (utils.het.report,
  gl.report.heterozygosity, gl.report.polyploid_heterozygosity).

## Findings

### I5 — Documentation (LOW)

A brief @noRd header (purpose, formula, NA handling). Unexported
status unchanged.

## Coverage

test-utils.stats.R — 2 assertions. No flips.

## Approval

All findings approved via the approval boxes (2026-09-01).

## Outcome

Documentation added (I5). Suite: 2/2. PR recorded below.

```json
{"function": "utils.stats", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "PASS",
 "findings": [{"id": "I5", "severity": "LOW", "rules": ["DOC"], "loc": "R/utils.stats.r", "status": "applied"}],
 "datasets": ["testset.gl", "constructed"],
 "baseline_test": "tests/testthat/test-utils.stats.R",
 "pr": null}
```
