# Review: utils.plot.save (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.plot.save. Reviewed in the
  infrastructure wave (seven files, one approval round, per-function
  PRs), under the standing member directive that utility functions
  are not for end users.
- Datasets: testset.gl subsets; constructed matrices; ggplot fixtures;
  composed-command inspection (no PLINK binary required).
- Family mode: analysis/infrastructure utility.
- Checks skipped: Google Group not searched (not available: no
  browser session); dartr2shiny not present in the local workspace.

## Verdicts

- **Standards: FAIL** — the "No plot saved" note is ungated (1 line at
  verbose 0, verified); a duplicated getwd() line; @param ... claims
  the dots are passed to ggplot2::ggsave but they are never used
  (inert-parameter doc class); NULL returned visibly.
- **Spec: FAIL** — (1) verbose is never normalized: the documented
  default NULL crashes with "argument is of length zero" the moment
  any gated branch is evaluated (verified; latent because every
  in-family caller passes a numeric). (2) The nonexistent-directory
  fallback assigns tempfile() — a path that does not exist — as the
  directory, so the rescue path itself crashes saveRDS with "cannot
  open the connection" (verified). The happy path verifies: RDS
  written, reloadable.

## Findings

### I2 — Normalize verbose (MEDIUM)

verbose <- gl.check.verbosity(verbose) at entry.

### I3 — Fallback to tempdir() (MEDIUM)

The rescue path uses tempdir() (which exists) and reports the actual
location at verbose >= 2.

### I4 — Gates, docs, tidy (LOW)

"No plot saved" gated at verbose >= 2; duplicate line removed; the
ggsave claim dropped from @param ... (dots retained for signature
compatibility); invisible(NULL); @keywords internal (STAYS exported —
called across the family).

## Coverage

test-utils.plot.save.R — 6 assertions: happy path (anchor),
default-verbose crash (baseline), broken fallback (baseline), ungated
note (baseline). All pass pre-fix.

## Approval

All findings approved via the approval boxes (2026-09-01).

## Outcome

All findings applied (I2 verbose normalized, I3 tempdir fallback, I4 gates/docs). Suite: 6/6; the documented default no longer crashes and the rescue path saves to a real directory. PR recorded below.

```json
{"function": "utils.plot.save", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [{"id": "I2", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/utils.plot.save.r verbose", "status": "applied"},
  {"id": "I3", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/utils.plot.save.r fallback", "status": "applied"},
  {"id": "I4", "severity": "LOW", "rules": ["VRB", "DOC"], "loc": "R/utils.plot.save.r", "status": "applied"}],
 "datasets": ["testset.gl", "constructed"],
 "baseline_test": "tests/testthat/test-utils.plot.save.R",
 "pr": null}
```
