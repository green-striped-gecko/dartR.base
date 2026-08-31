# Review: gl.reassign.pop (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Package commit: 2739215 (dev, synced with upstream/dev);
  Branch: review-gl.reassign.pop; Datasets: testset.gl, testset.gs
- Family mode: modify (population assignment; genotypes untouched)
- Checks skipped: Google Group not searched (not available: no browser
  session).

## Verdicts

- **Standards: FAIL** — minor header drift only.
- **Spec: FAIL** — `as.pop` is never validated: a metric name that does
  not exist in ind.metrics assigns NULL to pop(x), silently destroying
  every population assignment (verified: nPop = 0, zero output at
  verbose 0). An object lacking ind.metrics entirely fails the same
  way. The valid path verified correct on both SNP and SilicoDArT data
  (sex → Female/Male/Unknown; history appended; input untouched).

## Findings

### F1 — No validation of as.pop; silent destruction of pop slot (MEDIUM, confidence: certain)

Rule: spec/DAT. Proposed: fatal error when ind.metrics is absent or
as.pop is not one of its columns, naming the available columns.
Behaviour change: previously a silent no-op that emptied pop(x).

### F2 — NA handling and header (LOW, confidence: certain)

Rule: VRB/DOC. When the chosen metric contains NAs the resulting pop
factor carries NA assignments with no notice — a gated (verbose >= 2)
warning reporting the NA count is added. Header: minor spacing;
datatype line in @param x says "SNP genotypes" though SilicoDArT is
accepted.

## Coverage

`tests/testthat/test-gl.reassign.pop.R` — 8 assertions (sex
reassignment on gl and gs, history, input untouched, silent-destruction
baseline). All pass pre-fix.

## Approval

Both findings approved via the approval boxes (2026-08-31); F1 approved as a behaviour change (fatal instead of silent destruction).

## Outcome

Both applied. Suite 7/7 (the two-assertion baseline block became one expect_error); the flipped assertion maps to F1 (bogus as.pop now fatal, naming available metrics). Verbose-3 end-to-end clean on SNP and SilicoDArT.

```json
{"function": "gl.reassign.pop", "package": "dartR.base", "family_mode": "modify",
 "commit": "2739215", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [
  {"id": "F1", "severity": "MEDIUM", "rules": ["spec", "DAT"], "loc": "R/gl.reassign.pop.r as.pop", "status": "applied"},
  {"id": "F2", "severity": "LOW", "rules": ["VRB", "DOC"], "loc": "R/gl.reassign.pop.r", "status": "applied"}],
 "datasets": ["testset.gl", "testset.gs"],
 "baseline_test": "tests/testthat/test-gl.reassign.pop.R", "pr": 282}
```
