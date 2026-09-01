# Review: utils.flag.start (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.flag.start. Reviewed in the
  infrastructure wave (seven files, one approval round, per-function
  PRs), under the standing member directive that utility functions
  are not for end users.
- Datasets: testset.gl subsets; constructed matrices; ggplot fixtures;
  composed-command inspection (no PLINK binary required).
- Family mode: analysis/infrastructure utility.
- Checks skipped: Google Group not searched (not available: no
  browser session); dartr2shiny not present in the local workspace.

## Verdicts

- **Standards: PASS (with notes)** — no @description section; the
  verbose == 5 branch collapses to the else branch whenever build is
  NULL (redundant nesting).
- **Spec: PASS** — verbosity contract verified (silent at 0, start
  line at 1+, build info at 5); NULL func stops informatively;
  returns the function name invisibly. The family's most-called
  utility behaved correctly throughout the campaign.

## Findings

### I1 — Tidy and docs (LOW)

@description added; redundant nesting simplified; @keywords internal
(STAYS exported — called by every gl.* function across the family).

## Coverage

test-utils.flag.start.R — 7 assertions. No flips.

## Approval

All findings approved via the approval boxes (2026-09-01).

## Outcome

Docs/tidy applied (I1). Suite: 7/7; verbosity contract anchors unchanged. PR recorded below.

```json
{"function": "utils.flag.start", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "PASS", "verdict_spec": "PASS",
 "findings": [{"id": "I1", "severity": "LOW", "rules": ["DOC", "STY"], "loc": "R/utils.flag.start.r", "status": "applied"}],
 "datasets": ["testset.gl", "constructed"],
 "baseline_test": "tests/testthat/test-utils.flag.start.R",
 "pr": null}
```
