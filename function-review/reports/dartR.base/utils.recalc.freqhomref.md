# Review: utils.recalc.freqhomref (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.recalc.freqhomref. Reviewed as part
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

- **Standards: FAIL** — the monomorphs check `flag == FALSE` is
  `if (logical(0))` when the flags slot is absent (verified crash:
  "argument is of length zero"); the seealso block references
  nonexistent functions; the @title contains a stray "#'" that renders into the docs.
- **Spec: FAIL** — SilicoDArT input runs through the diploid 0/1/2
  arithmetic silently (absences counted as homozygous reference); the metric itself is verified exact
  against an independent recomputation on a 20-population subset, and
  verbose 0 is silent.

## Findings

### B2 — NULL-unsafe monomorphs check (MEDIUM)

Proposed: `!isTRUE(...)` (the idiom utils.recalc.callrate already
uses).

### B3 — SilicoDArT accepted (MEDIUM) [escalation: silico now fatal]

Proposed: utils.check.datatype(x, accept = "SNP"). All in-package
callers verified SNP-side (gl.recalc.metrics dispatches silico to
avgpic/callrate only).

### B1 — Visibility (LOW) [battery-wide, member directive]

Already unexported; @keywords internal added.

### B4 — Docs (LOW)

seealso corrected; the stray "#'" removed from the title; creation-message gate standardised.

## Coverage

test-utils.recalc.freqhomref.R — 5 assertions (recomputation anchor, flags-crash
baseline, silico baseline). All pass pre-fix.

## Approval

All findings approved via the approval boxes (2026-09-01), including
the battery-wide visibility shape (B1: @keywords internal on all
seven; utils.reset.flags stays exported for dartR.sim) and the
silico-fatal escalation (B3).

## Outcome

All findings applied. Battery suite: 40/40 assertions pass; the
metric recomputation anchors are unchanged (arithmetic verified
exact pre- and post-fix); gl.recalc.metrics runs end-to-end on SNP
and SilicoDArT. PR #304.

```json
{"function": "utils.recalc.freqhomref", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [{"id": "B2", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/utils.recalc.freqhomref.r monomorphs check", "status": "applied"},
  {"id": "B3", "severity": "MEDIUM", "rules": ["spec", "DAT"], "loc": "R/utils.recalc.freqhomref.r datatype", "status": "applied"},
  {"id": "B1", "severity": "LOW", "rules": ["DOC"], "loc": "R/utils.recalc.freqhomref.r visibility", "status": "applied"},
  {"id": "B4", "severity": "LOW", "rules": ["DOC"], "loc": "R/utils.recalc.freqhomref.r docs", "status": "applied"}],
 "datasets": ["testset.gl", "testset.gs"],
 "baseline_test": "tests/testthat/test-utils.recalc.freqhomref.R",
 "pr": 304}
```
