# Review: utils.recalc.avgpic (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.recalc.avgpic. Reviewed as part
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

- **Standards: FAIL** — the NULL-unsafe monomorphs check (crash class
  verified); `c` shadows base::c; seealso references nonexistent
  functions.
- **Spec: PASS** — both datatype branches verified exact against
  independent recomputations (SNP: OneRatioRef/Snp, PICRef/Snp,
  AvgPIC; SilicoDArT: OneRatio, PIC); proper datatype dispatch (the
  battery's only member that handles silico correctly); silent at
  verbose 0. All-NA loci yield NaN metrics (noted; consistent with
  division by a zero count).

## Findings

### B2 — NULL-unsafe monomorphs check (MEDIUM)

`!isTRUE(...)` as in callrate.

### B1 — Visibility (LOW) [battery-wide, member directive]

Already unexported; @keywords internal added.

### B4 — Docs and tidy (LOW)

seealso corrected; `c` renamed to ctot.

## Coverage

test-utils.recalc.avgpic.R — 6 assertions (SNP + silico anchors,
flags-crash baseline). All pass pre-fix.

## Approval

All findings approved via the approval boxes (2026-09-01), including
the battery-wide visibility shape (B1: @keywords internal on all
seven; utils.reset.flags stays exported for dartR.sim) and the
silico-fatal escalation (B3).

## Outcome

All findings applied. Battery suite: 40/40 assertions pass; the
metric recomputation anchors are unchanged (arithmetic verified
exact pre- and post-fix); gl.recalc.metrics runs end-to-end on SNP
and SilicoDArT. PR #307.

```json
{"function": "utils.recalc.avgpic", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "PASS",
 "findings": [{"id": "B2", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/utils.recalc.avgpic.r monomorphs check", "status": "applied"},
  {"id": "B1", "severity": "LOW", "rules": ["DOC"], "loc": "R/utils.recalc.avgpic.r visibility", "status": "applied"},
  {"id": "B4", "severity": "LOW", "rules": ["DOC", "STY"], "loc": "R/utils.recalc.avgpic.r tidy", "status": "applied"}],
 "datasets": ["testset.gl", "testset.gs"],
 "baseline_test": "tests/testthat/test-utils.recalc.avgpic.R",
 "pr": 307}
```
