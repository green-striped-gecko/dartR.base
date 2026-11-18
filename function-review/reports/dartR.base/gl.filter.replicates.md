# Review: gl.filter.replicates (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Package commit: 2739215 (dev, synced with upstream/dev);
  Branch: review-gl.filter.replicates; Datasets: platypus.gl subset
  with two exact duplicate individuals injected
- Family mode: modify (individual removal via gl.drop.ind)
- Custodian: Luis Mijangos
- Checks skipped: Google Group not searched (not available: no browser
  session).

## Verdicts

- **Standards: FAIL** — no datatype/genlight check; no validation of
  replicates.report at all; the history of the returned object records
  gl.drop.ind's call rather than gl.filter.replicates' (verified —
  the history-delegation class, PR #251 pattern); gl.drop.ind runs at
  the caller's verbosity producing a second "Completed:" banner; the
  return is visible where filters return invisibly; typo "functionv".
- **Spec: FAIL** — three verified defects, two shared with the report.
  (1) Tied-missingness pairs lose BOTH members: the doubled table plus
  the orientation-dependent drop rule removed 5 individuals from the
  test object where 3 (one per replicate pair) is correct — for exact
  duplicates this discards both copies of the sample. (2) Feeding the
  report's no-pairs string return crashes ("$ operator is invalid for
  atomic vectors"). (3) Re-thresholding to an empty set crashes via
  gl.drop.ind ("no individuals to drop") instead of returning the
  object unchanged.

## Findings

### F1 — Both tied replicates dropped; empty re-threshold crash (HIGH, confidence: certain)

Rule: spec axis. Location: drop-rule block.

Proposed change: deduplicate to unordered pairs before applying the
drop rule, with the same deterministic tie-break as the report (drop
ind2); when re-thresholding leaves no pairs, return the object
unchanged with a gated (verbose >= 1) message instead of crashing.

### F2 — replicates.report unvalidated (MEDIUM, confidence: certain)

Rule: spec/API. Location: input handling.

Proposed change: fatal error with a clear message unless
replicates.report is a list carrying a table.rep data frame; an empty
table.rep (the report's fixed no-pairs return) is handled as
nothing-to-drop.

### F3 — History delegation and verbosity (MEDIUM, confidence: certain)

Rule: FS8/VRB. Location: gl.drop.ind call.

Proposed change: call gl.drop.ind at verbose 0, itemise the dropped
individuals at verbose >= 2, snapshot the input's history and append a
single gl.filter.replicates entry (the PR #251/#266 pattern); datatype
check added.

### F4 — Return visibility and header (LOW, confidence: certain)

Rule: VRB/DOC. Invisible return (house standard for filters); typo
"functionv"; docs note that thresholds can only be tightened relative
to those used in gl.report.replicates (the table only contains pairs
that passed the original thresholds).

## Report notes (other functions / not fixed here)

- Legacy underscore parameter naming shared with the report pair —
  left to the custodian (recorded in the report review too).

## Coverage

`tests/testthat/test-gl.filter.replicates.R` — 8 assertions: 5-drop
tied baseline (correct member kept for the unequal pair), history
entry baseline, string-report crash baseline, empty-re-threshold crash
baseline. All 8 pass on the pre-fix code.

## Approval

All four findings approved via the approval boxes (2026-08-31).

## Outcome

All four findings applied. Pairs are canonicalised (ind1/ind2 in
alphabetical order) and deduplicated before the drop rule, so the fix
works with report tables generated both before (doubled) and after
(deduplicated) the companion gl.report.replicates fix (PR #287) —
verified against the OLD-format table. Characterization suite: 11/11
pass; flips map to F1 (29 individuals retained, one member per tied
pair, empty re-threshold returns unchanged), F2 (clear validation
error), F3 (history entry is gl.filter.replicates, single Completed
banner). End-to-end at verbose 3 clean with the dropped individuals
itemised. NEWS entry added.

```json
{"function": "gl.filter.replicates", "package": "dartR.base", "family_mode": "modify",
 "commit": "2739215", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [
  {"id": "F1", "severity": "HIGH", "rules": ["spec"], "loc": "R/gl.filter.replicates.r drop rule", "status": "applied"},
  {"id": "F2", "severity": "MEDIUM", "rules": ["spec", "API"], "loc": "R/gl.filter.replicates.r input", "status": "applied"},
  {"id": "F3", "severity": "MEDIUM", "rules": ["FS8", "VRB"], "loc": "R/gl.filter.replicates.r delegation", "status": "applied"},
  {"id": "F4", "severity": "LOW", "rules": ["VRB", "DOC"], "loc": "R/gl.filter.replicates.r", "status": "applied"}],
 "datasets": ["platypus.gl (constructed duplicates)"],
 "baseline_test": "tests/testthat/test-gl.filter.replicates.R", "pr": 288}
```
