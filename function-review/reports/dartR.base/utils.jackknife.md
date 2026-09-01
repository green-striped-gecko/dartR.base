# Review: utils.jackknife (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.jackknife. Reviewed in the
  distance wave (five files, one approval round, per-function PRs),
  under the standing member directive that utility functions are not
  for end users.
- Datasets: constructed 10x30 SNP fixture with missing data (plus its
  allele-relabelled mirror), testset.gs, raw-vector sequence fixtures,
  brute-force references.
- Family mode: analysis (distance/resampling kernel).
- Checks skipped: Google Group not searched (not available: no
  browser session); dartr2shiny not present in the local workspace.

## Verdicts

- **Standards: FAIL** — 6 lines leak at verbose 0 from the
  gl.set.verbosity save/restore performed inside every replicate
  (verified); `recal =` relies on partial matching for `recalc`;
  "to conducts" typo; a stale developer TODO comment.
- **Spec: FAIL** — `length(unit == 1)` takes the length of the
  comparison result rather than of unit (the luck-idiom class), so a
  unit vector of length > 1 crashes with "the condition has length
  > 1" instead of reaching the function's own informative stop
  (verified). The resampling itself verifies: one result per unit,
  correct subsetting (pop jackknife returns the per-population
  complements).

## Findings

### D10 — length(unit == 1) luck idiom (MEDIUM)

length(unit) == 1; the informative stop becomes reachable.

### D11 — Verbose-0 leak from the verbosity dance (MEDIUM)

The gl.set.verbosity save/restore calls are silenced (their output
captured); FUN's own output continues to be suppressed as before.

### D12 — Tidy and docs (LOW)

recal -> recalc; unique(pop(x)) coerced to character; typos; a note
that FUN must be available on worker processes for n.cores > 1;
@keywords internal (stays exported, called by gl.diagnostics.hwe).

## Coverage

test-utils.jackknife.R — 6 assertions: pop-jackknife anchor,
verbose-0 leak (baseline), length crash (baseline), invalid-unit
stop. All pass pre-fix.

## Approval

All findings approved via the approval boxes (2026-09-01); utils.hamming
disposition: REMOVE.

## Outcome

All findings applied (D10 length check, D11 silent verbosity dance, D12 tidy/docs). Suite: 6/6; the pop-jackknife anchor unchanged and now silent at verbose 0. PR recorded below.

```json
{"function": "utils.jackknife", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [{"id": "D10", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/utils.jackknife.R unit check", "status": "applied"},
  {"id": "D11", "severity": "MEDIUM", "rules": ["VRB"], "loc": "R/utils.jackknife.R verbosity dance", "status": "applied"},
  {"id": "D12", "severity": "LOW", "rules": ["DOC", "STY"], "loc": "R/utils.jackknife.R", "status": "applied"}],
 "datasets": ["constructed", "testset.gs"],
 "baseline_test": "tests/testthat/test-utils.jackknife.R",
 "pr": null}
```
