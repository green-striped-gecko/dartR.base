# Review: gl.define.pop (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Package commit: 2739215 (dev, synced with upstream/dev);
  Branch: review-gl.define.pop; Datasets: testset.gl
- Family mode: modify (population assignment; genotypes untouched)
- Checks skipped: Google Group not searched (not available: no browser
  session).

## Verdicts

- **Standards: FAIL** — the not-present-individual warning prints at
  verbose 0 (verified); a full gl.filter.monomorphs run (0.4 s on the
  small testset) is executed solely to warn that monomorphs exist —
  irrelevant to defining a population (the preamble class removed from
  gl.filter.locmetric, PR #264); header drift.
- **Spec: PASS** — assignment verified correct: listed individuals move
  to the new population, others untouched, empty old levels drop,
  single history call appended, input object unmodified.

## Findings

### F1 — Not-present warning leaks at verbose 0 (MEDIUM, confidence: certain)

Rule: VRB. The per-individual "not present -- ignored" warning is
ungated (verified: 1 line at verbose 0). Proposed: gate at
verbose >= 2 (the gl.recode.ind pattern, PR #267).

### F2 — Irrelevant monomorphs preamble (MEDIUM, confidence: certain)

Rule: STY/perf (irrelevant-preamble class, PR #264 precedent). A full
monomorph filter runs on every call solely to emit a warning that has
nothing to do with population definition. Proposed: remove the scan.

### F3 — Header and style conformance (LOW, confidence: certain)

Rules: DOC1, STY. `@return` after `@export`; dead
`is.na(length(pop(x)))` condition in the pop guard; per-case loop
filter rewritten `ind.list[ind.list %in% indNames(x)]`; the
"Assigned..." progress message unstyle (no report()) and printed
before the assignment happens; array→factor idiom.

## Report notes (other functions / not fixed here)

- gl.define.pop and gl.reassign.ind (v.2025.1) implement the same
  operation; governance may wish to deprecate one. Flagged for the
  coordinator — no action here.

## Coverage

`tests/testthat/test-gl.define.pop.R` — 10 assertions (assignment,
input untouched, history call, verbose-0 warning baseline, all-absent
fatal). All pass pre-fix.

## Approval

All three findings approved via the approval boxes (2026-08-31).

## Outcome

All applied. Suite 10/10; the flipped assertion maps to F1. Verbose-3 end-to-end clean.

```json
{"function": "gl.define.pop", "package": "dartR.base", "family_mode": "modify",
 "commit": "2739215", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "PASS",
 "findings": [
  {"id": "F1", "severity": "MEDIUM", "rules": ["VRB"], "loc": "R/gl.define.pop.r warning", "status": "applied"},
  {"id": "F2", "severity": "MEDIUM", "rules": ["STY"], "loc": "R/gl.define.pop.r preamble", "status": "applied"},
  {"id": "F3", "severity": "LOW", "rules": ["DOC1", "STY"], "loc": "R/gl.define.pop.r", "status": "applied"}],
 "datasets": ["testset.gl"],
 "baseline_test": "tests/testthat/test-gl.define.pop.R", "pr": null}
```
