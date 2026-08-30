# Review: gl.reassign.ind (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Package commit: 2739215 (dev, synced with upstream/dev);
  Branch: review-gl.reassign.ind; Datasets: testset.gl
- Family mode: modify (population assignment; genotypes untouched)
- Checks skipped: Google Group not searched (not available: no browser
  session).

## Verdicts

- **Standards: FAIL (minor)** — the five fatal exits use plain
  `stop("Fatal error: ...")` without the house `error()` styling
  (verified: no crayon codes in the condition); the empty-selection
  path uses R's `warning()` instead of a gated `cat(warn())`; header
  otherwise conforms (this is a recent v.2025.1 function).
- **Spec: PASS** — character, numeric, and logical selections verified
  correct; invalid names/indices/lengths fatal; reassignment to an
  existing population works; others retain assignments; single history
  call appended.

## Findings

### F1 — Error styling (LOW, confidence: certain)

Rule: STY. Five `stop("Fatal error: ...")` → `stop(error(...))`.

### F2 — Empty-selection path and dedup (LOW, confidence: certain)

Rule: VRB/STY. `warning()` → `cat(warn())` gated at verbose >= 2 (the
unchanged object is still returned); numeric indices deduplicated with
`unique()` so a repeated index is not double-processed (currently
harmless but untidy in the reported count).

## Report notes (other functions / not fixed here)

- Functional overlap with gl.define.pop — governance question for the
  coordinator (recorded in gl.define.pop's report too).

## Coverage

`tests/testthat/test-gl.reassign.ind.R` — 12 assertions (three
selection modes, fatal guards, existing-pop reassignment, history).
All pass pre-fix.

## Approval

Both findings approved via the approval boxes (2026-08-31).

## Outcome

Both applied. Suite 12/12 (no flips needed; behaviour unchanged for valid input). Verbose-3 end-to-end clean; empty selection silent at verbose 0.

```json
{"function": "gl.reassign.ind", "package": "dartR.base", "family_mode": "modify",
 "commit": "2739215", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "PASS",
 "findings": [
  {"id": "F1", "severity": "LOW", "rules": ["STY"], "loc": "R/gl.reassign.ind.r stops", "status": "applied"},
  {"id": "F2", "severity": "LOW", "rules": ["VRB", "STY"], "loc": "R/gl.reassign.ind.r selection", "status": "applied"}],
 "datasets": ["testset.gl"],
 "baseline_test": "tests/testthat/test-gl.reassign.ind.R", "pr": null}
```
