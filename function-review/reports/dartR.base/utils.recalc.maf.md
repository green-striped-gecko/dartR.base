# Review: utils.recalc.maf (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.recalc.maf. Reviewed as part
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
  verified); a dead commented gl.allele.freq call; seealso references
  nonexistent functions.
- **Spec: FAIL** — SilicoDArT accepted despite the doc's own "only
  applies to SNP genotype data"; maf itself verified exact
  (pmin(alf, 1-alf)) against an independent recomputation; silent at
  verbose 0. Note: the function recalculates FreqHets/FreqHomRef/
  FreqHomSnp as a side effect although the maf computation does not
  use them — documented behaviour, retained.

## Findings

### B2 — NULL-unsafe monomorphs check (MEDIUM)

`!isTRUE(...)` as in callrate.

### B3 — SilicoDArT accepted (MEDIUM) [escalation: silico now fatal]

accept = "SNP" (the doc already promises this). Callers verified.

### B1 — Visibility (LOW) [battery-wide, member directive]

Already unexported; @keywords internal added.

### B4 — Docs and tidy (LOW)

seealso corrected; dead commented line removed; creation-message gate
standardised.

## Coverage

test-utils.recalc.maf.R — 5 assertions. All pass pre-fix.

## Approval

All findings approved via the approval boxes (2026-09-01), including
the battery-wide visibility shape (B1: @keywords internal on all
seven; utils.reset.flags stays exported for dartR.sim) and the
silico-fatal escalation (B3).

## Outcome

All findings applied. Battery suite: 40/40 assertions pass; the
metric recomputation anchors are unchanged (arithmetic verified
exact pre- and post-fix); gl.recalc.metrics runs end-to-end on SNP
and SilicoDArT. PR recorded below.

```json
{"function": "utils.recalc.maf", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [{"id": "B2", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/utils.recalc.maf.r monomorphs check", "status": "applied"},
  {"id": "B3", "severity": "MEDIUM", "rules": ["spec", "DAT"], "loc": "R/utils.recalc.maf.r datatype", "status": "applied"},
  {"id": "B1", "severity": "LOW", "rules": ["DOC"], "loc": "R/utils.recalc.maf.r visibility", "status": "applied"},
  {"id": "B4", "severity": "LOW", "rules": ["DOC", "STY"], "loc": "R/utils.recalc.maf.r docs", "status": "applied"}],
 "datasets": ["testset.gl", "testset.gs"],
 "baseline_test": "tests/testthat/test-utils.recalc.maf.R",
 "pr": null}
```
