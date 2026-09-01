# Review: utils.reset.flags (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.reset.flags. Reviewed as part
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

- **Standards: FAIL** — the out-of-range `value` warning prints at
  verbose 0 (verified); the @title carries a stray "#'"; the SNP
  branch's verbose-2 message lists OneRatio/PIC twice with
  contradictory targets; a duplicated monomorphs flag assignment; a
  commented-out dead block.
- **Spec: FAIL** — (1) the `value` parameter ("set the default
  verbosity") is validated and then IGNORED: the verbosity slot is
  hardcoded to 2 (verified: value = 4 yields @other$verbose == 2) —
  the inert-parameter class. (2) The SNP branch invents a
  loc.metrics$monomorphs COLUMN for what is a flag, not a locus metric
  (the silico branch has this very block commented out); testset
  objects already carry the bogus column. Flag semantics themselves
  verified: SNP/silico branches reset their own metrics to `set` and
  disable the other datatype's flags.

## Findings

### R1 — value parameter inert (MEDIUM)

Honour it where the slot is created: x@other$verbose <- value. No
current caller passes value, so behaviour is unchanged for existing
code paths.

### R2 — Ungated coercion warning (LOW)

Gate at verbose >= 1.

### R3 — Bogus loc.metrics$monomorphs column; duplicate assignment (LOW)

Stop creating the column (flag-only, matching the silico branch);
remove the duplicated flag set and the dead commented block; fix the
contradictory verbose-2 message.

### B1 — Visibility (LOW) [battery-wide, member directive]

KEPT EXPORTED — dartR.sim calls it (utils.sims.r line 449) and it is
the one battery member in NAMESPACE. @keywords internal added so it
leaves the user-facing index while remaining callable cross-package.

### B4 — Docs (LOW)

Stray "#'" in the title; value/verbose param docs disentangled;
seealso corrected.

## Coverage

test-utils.reset.flags.R — 9 assertions (flag semantics both
datatypes, value-inert baseline, warning-leak baseline, bogus-column
baseline). All pass pre-fix.

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
{"function": "utils.reset.flags", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [{"id": "R1", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/utils.reset.flags.r value", "status": "applied"},
  {"id": "R2", "severity": "LOW", "rules": ["VRB"], "loc": "R/utils.reset.flags.r warning", "status": "applied"},
  {"id": "R3", "severity": "LOW", "rules": ["spec", "STY"], "loc": "R/utils.reset.flags.r monomorphs column", "status": "applied"},
  {"id": "B1", "severity": "LOW", "rules": ["DOC"], "loc": "R/utils.reset.flags.r visibility", "status": "applied"},
  {"id": "B4", "severity": "LOW", "rules": ["DOC"], "loc": "R/utils.reset.flags.r docs", "status": "applied"}],
 "datasets": ["testset.gl", "testset.gs"],
 "baseline_test": "tests/testthat/test-utils.reset.flags.R",
 "pr": null}
```
