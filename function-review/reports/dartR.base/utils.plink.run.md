# Review: utils.plink.run (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.plink.run. Reviewed in the
  infrastructure wave (seven files, one approval round, per-function
  PRs), under the standing member directive that utility functions
  are not for end users.
- Datasets: testset.gl subsets; constructed matrices; ggplot fixtures;
  composed-command inspection (no PLINK binary required).
- Family mode: analysis/infrastructure utility.
- Checks skipped: Google Group not searched (not available: no
  browser session); dartr2shiny not present in the local workspace.

## Verdicts

- **Standards: FAIL** — dead commented blocks; the on.exit/setwd
  workaround documented only by an inline aside.
- **Spec: FAIL** — the composed command is malformed twice (verified:
  "path/nonexistent_exe_999 --file hapmap1--out hapmap1"): (1) the
  default plink.path = "path" — documented as "plink is on the PATH,
  no path needed" — is pasted literally, producing a bogus "path/"
  prefix; (2) there is no space between syntax and "--out", gluing
  the last flag ("hapmap1--out"). Any call that does not both supply
  an explicit plink.path and end syntax with a trailing space runs a
  broken command. One caller in the family (gl.read.PLINK).

## Findings

### I11 — Command composition (MEDIUM) [escalation: composed command changes]

plink.path == "path" now means bare plink.cmd (PATH lookup as
documented); explicit paths joined with file.path; a space guaranteed
before --out. Callers that already worked around the gluing (trailing
space in syntax) produce a harmless double space.

### I12 — Tidy (LOW)

Dead commented blocks removed; @keywords internal (STAYS exported).

## Coverage

test-utils.plink.run.R — 2 assertions on the composed command
(baseline both defects). All pass pre-fix.

## Approval

All findings approved via the approval boxes (2026-09-01).

## Outcome

All findings applied (I11 command composition, I12 tidy). Suite: 2/2; the composed command is now well formed under the default path. PR recorded below.

```json
{"function": "utils.plink.run", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [{"id": "I11", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/utils.plink.run.r command", "status": "applied"},
  {"id": "I12", "severity": "LOW", "rules": ["STY"], "loc": "R/utils.plink.run.r", "status": "applied"}],
 "datasets": ["testset.gl", "constructed"],
 "baseline_test": "tests/testthat/test-utils.plink.run.R",
 "pr": null}
```
