# Review: utils.transpose (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.transpose. Reviewed in the
  infrastructure wave (seven files, one approval round, per-function
  PRs), under the standing member directive that utility functions
  are not for end users.
- Datasets: testset.gl subsets; constructed matrices; ggplot fixtures;
  composed-command inspection (no PLINK binary required).
- Family mode: analysis/infrastructure utility.
- Checks skipped: Google Group not searched (not available: no
  browser session); dartr2shiny not present in the local workspace.

## Verdicts

- **Standards: FAIL** — every line is narrated by an explanatory
  comment (tutorial-style noise); the file's Rd has been a source of
  formatting churn all campaign.
- **Spec: PASS** — verified: dimensions, individual/locus names and
  loc/ind metrics all swap correctly, and a double transpose
  reproduces the original genotypes exactly. Depends on matrix2gen
  from utils.impute.R (noted). Mixed-ploidy input would fail at the
  ploidy assignment (edge, noted).

## Findings

### I6 — Tidy and docs (LOW)

Narration comments reduced to the load-bearing ones; @keywords
internal (stays unexported).

## Coverage

test-utils.transpose.R — 5 assertions including the round-trip
anchor. No flips.

## Approval

All findings approved via the approval boxes (2026-09-01).

## Outcome

Tidy applied (I6). Suite: 5/5; round-trip anchor unchanged. PR recorded below.

```json
{"function": "utils.transpose", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "PASS",
 "findings": [{"id": "I6", "severity": "LOW", "rules": ["STY", "DOC"], "loc": "R/utils.transpose.R", "status": "applied"}],
 "datasets": ["testset.gl", "constructed"],
 "baseline_test": "tests/testthat/test-utils.transpose.R",
 "pr": null}
```
