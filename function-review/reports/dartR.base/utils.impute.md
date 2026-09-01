# Review: utils.impute (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.impute. Reviewed in the
  infrastructure wave (seven files, one approval round, per-function
  PRs), under the standing member directive that utility functions
  are not for end users.
- Datasets: testset.gl subsets; constructed matrices; ggplot fixtures;
  composed-command inspection (no PLINK binary required).
- Family mode: analysis/infrastructure utility.
- Checks skipped: Google Group not searched (not available: no
  browser session); dartr2shiny not present in the local workspace.

## Verdicts

- **Standards: FAIL** — the doc block is unfilled boilerplate
  ("[Custodian to provide a title]" etc.) attached to @name
  utils.impute, a function that does not exist — the file actually
  defines three helpers (matrix2gen, s_alleles, sample_genotype) and
  its ghost Rd has churned all campaign.
- **Spec: FAIL** — matrix2gen's parallel branch assigns to i@gen where
  no object i exists: parallel = TRUE crashes with "object 'i' not
  found" (verified) — and gl.impute exposes parallel to end users, so
  this is a live user-facing crash. The serial branch and both
  sampling helpers verify (SNPbin construction; q=0/q=1 degenerate
  probabilities; NA propagation; Hardy-Weinberg probabilities).

## Findings

### I7 — matrix2gen parallel branch (MEDIUM)

The stray i@gen assignment removed; the mclapply result is returned
directly (mirroring the serial branch). Windows note: mclapply falls
back to serial with mc.cores = 1.

### I8 — Real documentation (LOW)

The placeholder block replaced by @noRd headers on each helper; the
ghost utils.impute Rd removed (ends the churn).

## Coverage

test-utils.impute.R — 9 assertions: serial anchor, parallel crash
(baseline), sampling contracts. All pass pre-fix.

## Approval

All findings approved via the approval boxes (2026-09-01).

## Outcome

All findings applied (I7 parallel branch returns the mclapply result, I8 real @noRd docs, ghost Rd removed). Suite: 9/9; gl.impute caller smoke clean. PR recorded below.

```json
{"function": "utils.impute", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [{"id": "I7", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/utils.impute.R matrix2gen", "status": "applied"},
  {"id": "I8", "severity": "LOW", "rules": ["DOC"], "loc": "R/utils.impute.R placeholders", "status": "applied"}],
 "datasets": ["testset.gl", "constructed"],
 "baseline_test": "tests/testthat/test-utils.impute.R",
 "pr": null}
```
