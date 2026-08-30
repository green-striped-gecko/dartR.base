# Review: gl.merge.pop (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Package commit: 2739215 (dev, synced with upstream/dev);
  Branch: review-gl.merge.pop; Datasets: testset.gl
- Family mode: modify (population assignment; genotypes untouched)
- Checks skipped: Google Group not searched (not available: no browser
  session).

## Verdicts

- **Standards: FAIL** — the empty-old validation sits INSIDE the
  `verbose >= 1` gate: `old = character(0)` is fatal at verbose >= 1
  but a silent no-op at verbose 0 (verified — the
  verbosity-dependent-validation class); a redundant
  `is(x, "genlight")` check follows the datatype check; header drift.
- **Spec: FAIL** — old populations that do not exist are silently
  ignored (verified: no error, no message, object unchanged even at
  verbose 0) — a user who mistypes a population name believes the
  merge happened. gl.rename.pop received a fatal error for the same
  class (PR #268, FORBID chosen). Correct merges verified: multi-pop
  merge, single-pop rename, counts and history all correct.

## Findings

### F1 — Verbosity-dependent validation (MEDIUM, confidence: certain)

Rule: VRB (ratified class). The `length(old)` check lives in the
verbose >= 1 announcement block. Proposed: validate old/new upfront,
ungated; announcement stays gated.

### F2 — Nonexistent old populations silently ignored (MEDIUM, confidence: certain)

Rule: spec axis. Proposed: fatal error naming the missing populations
(mirrors gl.rename.pop, PR #268). Behaviour change: previously a
silent no-op.

### F3 — Tidy and header (LOW, confidence: certain)

Rules: DOC1, STY. Redundant genlight check removed; duplicate
is.null(old) validation consolidated; description opening paragraph is
a copy-paste about csv metadata files (belongs to gl.reassign.pop);
`@return` after `@export`.

## Coverage

`tests/testthat/test-gl.merge.pop.R` — 14 assertions (merge, rename,
counts, history, silent-nonexistent baseline, verbosity-dependent
empty-old baseline, NULL guards). All pass pre-fix.

## Approval

(to be recorded)

## Outcome

(to be recorded)

```json
{"function": "gl.merge.pop", "package": "dartR.base", "family_mode": "modify",
 "commit": "2739215", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [
  {"id": "F1", "severity": "MEDIUM", "rules": ["VRB"], "loc": "R/gl.merge.pop.r validation", "status": "proposed"},
  {"id": "F2", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl.merge.pop.r old check", "status": "proposed"},
  {"id": "F3", "severity": "LOW", "rules": ["DOC1", "STY"], "loc": "R/gl.merge.pop.r", "status": "proposed"}],
 "datasets": ["testset.gl"],
 "baseline_test": "tests/testthat/test-gl.merge.pop.R", "pr": null}
```
