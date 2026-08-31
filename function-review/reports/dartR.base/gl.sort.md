# Review: gl.sort (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Package commit: 2739215 (dev, synced with upstream/dev);
  Branch: review-gl.sort; Datasets: bandicoot.gl, testset.gl
- Family mode: modify (reordering; genotype values untouched)
- Custodian: Bernd Gruber
- Checks skipped: Google Group not searched (not available: no browser
  session).

## Verdicts

- **Standards: FAIL** — the history entry is appended as
  `c(match.call())`, which coerces the call to a LIST and corrupts the
  history chain (verified: entry class "list", is.call FALSE — replay
  via gl.print.history breaks on it); there is no FLAG SCRIPT END
  block, so "Completed:" never prints at any verbosity (verified); the
  no-chromosome warning under order.by.chr.pos is ungated (prints at
  verbose 0); the dartR-conversion notice is gated `verbose > 2`
  where the house convention is `>= 2`.
- **Spec: PASS** — default pop sort, explicit order.by pop sort
  (levels honoured), ind sort by name and by external vector all
  verified correct; ind.metrics and latlon travel on the same
  permutation; a missing latlon slot is tolerated (NULL subsetting
  yields NULL); length mismatches fatal. Note: bandicoot.gl's
  ind.metrics$id column does not match indNames even before sorting —
  a fixture quirk, not a gl.sort defect.

## Findings

### F1 — History entry appended as a list (MEDIUM, confidence: certain)

Rule: FS8. `c(match.call())` → `match.call()`.

### F2 — Missing script-end block; warning gates (MEDIUM, confidence: certain)

Rule: FS/VRB. Add the standard FLAG SCRIPT END (`Completed:` at
verbose >= 1); gate the chr-pos warning at verbose >= 2; align the
dartR-conversion notice to verbose >= 2.

### F3 — Duplicate check and header (LOW, confidence: certain)

Rules: DOC2, STY. The sort.by='ind' length re-check duplicates the
upfront length validation and carries a misleading message ("does not
contain all levels" for a length mismatch) — removed; verbose param
canon wording; sort.by doc quote style; @return prose tidied.

## Coverage

`tests/testthat/test-gl.sort.R` — 13 assertions (three sort modes,
order.by honoured, permutation integrity of ind.metrics, fatal guards,
list-history baseline, no-Completed baseline). All pass pre-fix.

## Approval

All three findings approved via the approval boxes (2026-08-31).

## Outcome

All applied. Suite 13/13; flipped assertions map to F1 (history entry is a call) and F2 (Completed: prints at verbose >= 1). Verbose-3 end-to-end clean.

```json
{"function": "gl.sort", "package": "dartR.base", "family_mode": "modify",
 "commit": "2739215", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "PASS",
 "findings": [
  {"id": "F1", "severity": "MEDIUM", "rules": ["FS8"], "loc": "R/gl.sort.r history", "status": "applied"},
  {"id": "F2", "severity": "MEDIUM", "rules": ["FS", "VRB"], "loc": "R/gl.sort.r end block", "status": "applied"},
  {"id": "F3", "severity": "LOW", "rules": ["DOC2", "STY"], "loc": "R/gl.sort.r", "status": "applied"}],
 "datasets": ["bandicoot.gl", "testset.gl"],
 "baseline_test": "tests/testthat/test-gl.sort.R", "pr": 283}
```
