# Review: gl.filter.hwe (dartR.base) — post excess.het migration

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev) + migration
  commits on branch migrate-excess.het-to-hwe
- Datasets: testset.gl, LBP (stub-equivalence fixture)
- Family mode: modify
- Context: reviewed together with gl.report.hwe immediately after the
  approved migration (direction + min.hobs parameters; deprecation
  stubs).
- Checks skipped: as for gl.report.hwe.

## Verdicts

- **Standards: FAIL** — same defect family as the report sibling: the
  population-skipping blocks sit inside `verbose >= 2` gates, the
  HardyWeinberg check returns -1, the alpha warning is ungated and
  mis-worded, `@family filter functions` instead of "matched filter",
  DOC1/DOC2/DOC7 drift.
- **Spec: FAIL** — the tested population pool depends on verbosity
  (structurally identical to the report sibling's verified defect; on
  testset.gl with default settings the removed-locus set happens to
  coincide, but the FDR pool provably differs and other datasets will
  diverge). Otherwise sound: filtering removes exactly the loci the
  report flags (verified for direction = 'excess', 61/61), metadata
  stays in sync, exactly one history entry, invisible return, and the
  gl.filter.excess.het stub removes exactly the original's 6 loci on
  LBP.

## Findings

### F1 — Tested population set depends on verbosity (HIGH, confidence: certain)

Rule: VRB. Location: `R/gl.filter.hwe.r` (population-skipping blocks).
Identical structure to gl.report.hwe R1 (verified there: 30 vs 23
populations tested). The filter's removed-locus set therefore depends on
verbosity whenever a small/monomorphic population contributes a
significant test or shifts the adjustment pool. Proposed change: same
fix — removals unconditional, warnings gated.

### F2 — HardyWeinberg check returns -1 (MEDIUM, confidence: certain)

Rule: spec axis (return contract). `cat(error(...)); return(-1)` from a
function documented to return a genlight object. Proposed change:
`stop(error(...))`.

### F3 — Alpha warning ungated and mis-worded (LOW, confidence: certain)

Rule: VRB2. Same as report R4. Proposed change: gate at verbose >= 1;
"a value between 0 and 1".

### F4 — Header conformance (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC2, DOC7. `@family filter functions` corrected to
"matched filter"; tag order; DOC2 wording; Author(s)/Custodian labels;
boundary rows with Prob exactly equal to alpha receive Sig = NA (noted).

## Report notes (other functions / not fixed here)

- Migration verified: direction filter removes exactly the report's
  flagged loci; min.hobs reproduces the published excess.het pool; the
  stub matches the original's removals on LBP (6/6) and appends exactly
  one history entry.
- The internal per-population gl.filter.monomorphs calls operate on
  discarded copies — no history leak to the returned object (verified:
  +1 entry only).

## Coverage

Characterization baseline: `tests/testthat/test-gl.filter.hwe.R` — 11
assertions: excess-direction filter/report agreement (61/61), default
'both' behaviour preserved (nLoc 249), sync, single history entry,
invisible return, input untouched, verbose-0 silence, stub equivalence
(6 loci, history +1). All 11 pass post-migration, pre-review-fix.

## Approval

Decisions recorded 2026-08-30 via approval boxes — all findings
**approved**, including the escalation-gate consequences (verbosity-
independent tested set; fatal errors for missing packages; ggtern only
when plotting). A follow-up decision during the caller grep: the
one-line dependency fix in gl.diagnostics.hwe (pass plot.out = FALSE)
was approved and applied.

## Outcome

Applied all approved findings on the migration branch. Verification: all
24 characterization assertions across the pair pass post-fix; every diff
from baseline maps to an approved finding (tested set identical at every
verbosity — 245 rows / 23 populations on testset.gl; missing packages
now stop; alpha warning gated; families corrected). The migration
invariants held throughout: direction partitions exactly, Het.exp
matches independent recomputation, filter removes exactly the report's
flagged loci, and the deprecation stubs reproduce the original
excess.het output exactly on LBP (6 loci, 0 extras, one history entry).
Caller grep: gl.diagnostics.hwe was the only code caller — it relied on
the old return(-1) contract and has been given the approved one-line fix
(plot.out = FALSE), which also spares it from rendering ternary plots
mid-diagnostics; its now-dead -1 guard is noted for its own review.
dartr2shiny: not present in the workspace. NEWS entry added.

```json
{
  "function": "gl.filter.hwe",
  "package": "dartR.base",
  "family_mode": "modify",
  "commit": "2739215+migration",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "HIGH", "rules": ["VRB"], "loc": "R/gl.filter.hwe.r pop-skipping blocks", "status": "proposed"},
    {"id": "F2", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl.filter.hwe.r package check", "status": "proposed"},
    {"id": "F3", "severity": "LOW", "rules": ["VRB2"], "loc": "R/gl.filter.hwe.r alpha check", "status": "proposed"},
    {"id": "F4", "severity": "LOW", "rules": ["DOC1", "DOC2", "DOC7"], "loc": "R/gl.filter.hwe.r header", "status": "proposed"}
  ],
  "datasets": ["testset.gl", "LBP"],
  "baseline_test": "tests/testthat/test-gl.filter.hwe.R",
  "pr": null
}
```
