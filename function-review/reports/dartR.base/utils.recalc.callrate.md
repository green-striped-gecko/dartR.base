# Review: utils.recalc.callrate (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.recalc.callrate. Reviewed as part
  of the recalc-battery run (seven functions, one approval round,
  per-function PRs), with the member directive that utility functions
  are not for end users.
- Datasets: testset.gl (10 populations dropped, metrics recomputed),
  testset.gs.
- Family mode: analysis (metric-maintenance utility).
- Custodian: Luis Mijangos.
- Checks skipped: Google Group not searched (not available: no browser
  session); dartr2shiny not present in the local workspace (noted).

## Verdicts

- **Standards: PASS (with notes)** — the NULL-safe isTRUE monomorphs
  check is the battery's good example; docs carry the internal-use
  warning; the seealso block references nonexistent functions
  (utils.recalc.metrics, gl.recalc.maf, gl.recalc.rdepth as gl.*).
- **Spec: PASS** — CallRate verified exact against an independent
  recomputation on a 20-population subset; datatype-agnostic (works
  for SilicoDArT, verified); silent at verbose 0.

## Findings

### B1 — Visibility (LOW) [battery-wide, member directive]

Already unexported (the @export tag is commented out). Proposed:
@keywords internal so the man page leaves the user-facing index.

### B4 — Docs (LOW)

seealso corrected to existing functions; title covers individuals as
well as populations; creation-message gate standardised to verbose >= 3.

## Coverage

test-utils.recalc.callrate.R — 6 assertions (recomputation anchor,
flags-slot-absent safety, silico support). All pass pre-fix.

## Approval

All findings approved via the approval boxes (2026-09-01), including
the battery-wide visibility shape (B1: @keywords internal on all
seven; utils.reset.flags stays exported for dartR.sim) and the
silico-fatal escalation (B3).

## Outcome

All findings applied. Battery suite: 40/40 assertions pass; the
metric recomputation anchors are unchanged (arithmetic verified
exact pre- and post-fix); gl.recalc.metrics runs end-to-end on SNP
and SilicoDArT. PR #302.

```json
{"function": "utils.recalc.callrate", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "PASS", "verdict_spec": "PASS",
 "findings": [{"id": "B1", "severity": "LOW", "rules": ["DOC"], "loc": "R/utils.recalc.callrate.r visibility", "status": "applied"},
  {"id": "B4", "severity": "LOW", "rules": ["DOC"], "loc": "R/utils.recalc.callrate.r docs", "status": "applied"}],
 "datasets": ["testset.gl", "testset.gs"],
 "baseline_test": "tests/testthat/test-utils.recalc.callrate.R",
 "pr": 302}
```
