# Review: gl.report.hwe (dartR.base) — post excess.het migration

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev) + migration
  commits on branch migrate-excess.het-to-hwe
- Datasets: testset.gl (SNP, 30 pops incl. 7 with <= 5 individuals), LBP
  (stub-equivalence fixture)
- Family mode: report
- Context: reviewed together with gl.filter.hwe immediately after the
  approved migration of gl.report.excess.het/gl.filter.excess.het into
  this pair (direction + min.hobs parameters; deprecation stubs).
- Checks skipped: ternary-plot rendering verified only for object
  construction (headless); Exact-test p-values not re-derived (delegated
  to HardyWeinberg::HWExactStats; the ChiSquare path was verified against
  HWChisq directly); Google Group not searched (not available: no browser
  session).

## Verdicts

- **Standards: FAIL** — no `@family` tag at all (the function sits in no
  family index), alpha warning ungated at verbose = 0 with wrong wording
  ("must be an integer"), legacy save2tmp/tempfile output instead of
  plot.file/plot.dir, DOC1/DOC2/DOC7 drift.
- **Spec: FAIL** — the tested set depends on verbosity: the skipping of
  monomorphic and small-sample populations sits INSIDE `verbose >= 2`
  gates, so verbose 0 tests all 30 populations (311 rows at
  sig_only = FALSE) while verbose 2 tests 23 (245 rows); any
  multiple-comparison adjustment pool differs accordingly. Missing
  packages return -1 from a function documented to return a dataframe,
  and ggtern is demanded even when plot.out = FALSE (verified: blocked
  non-plotting use on this machine until installed). The core statistics
  are sound: ChiSquare p-values match HWChisq directly; the new
  direction partition and Het.exp column were verified exact; the
  migration stub reproduces the original excess.het locus set on LBP
  (6/6).

## Findings

### R1 — Tested population set depends on verbosity (HIGH, confidence: certain)

Rule: VRB (object/results must not depend on verbosity — the
gl.recode.pop class, here affecting RESULTS). Location:
`R/gl.report.hwe.r` (population-skipping blocks).

`poplist <- poplist[-pops_to_remove]` (monomorphic pops) and the
small-sample equivalent are inside `if (verbose >= 2)`. Verified:
sig_only = FALSE returns 311 rows / 30 pops at verbose 0 vs 245 rows /
23 pops at verbose 2. With multi_comp = TRUE the adjustment pool differs
between verbosity levels.

Proposed change: perform both removals unconditionally; keep only the
warnings gated. Consequence: verbose 0/1 results change to match
verbose >= 2 (small and monomorphic populations always skipped, as the
messages claim).

### R2 — Missing packages return -1; ggtern demanded without plotting (MEDIUM, confidence: certain)

Rule: spec axis (return contract); DEP. Location: package checks.

`cat(error(...)); return(-1)` for both HardyWeinberg and ggtern — a
value of -1 from a function documented to return a dataframe, invisible
to scripts. And ggtern is required even when plot.out = FALSE (verified
empirically). Proposed change: `stop(error(...))` for HardyWeinberg;
move the ggtern check inside `if (plot.out)`.

### R3 — npop column absent when sig_only = FALSE (LOW, docs, confidence: certain)

`@return` documents the npop column unconditionally; it exists only when
sig_only = TRUE (verified). Proposed change: docs state the condition.

### R4 — Alpha warning ungated and mis-worded (LOW, confidence: certain)

Rule: VRB2. The out-of-range alpha warning prints at verbose = 0 and
says alpha "must be an integer between 0 and 1". Proposed change: gate
at verbose >= 1; "a value between 0 and 1".

### R5 — Header conformance (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC2, DOC7. No `@family` tag (added: "matched report");
tag order; DOC2 wording; Author(s)/Custodian labels; boundary rows with
Prob exactly equal to alpha receive Sig = NA (noted in docs). The legacy
save2tmp/plot.out API is retained as-is in this pass (renaming to
plot.file/plot.dir/plot.display is an API change beyond this migration's
scope — recorded as a deferred note).

## Report notes (other functions / not fixed here)

- Migration verified: direction excess/deficit partition the departure
  set exactly; Het.exp matches independent recomputation; min.hobs
  screens before the adjustment so the FDR pool equals the published
  workflow's; the gl.report.excess.het stub reproduces the original's 6
  flagged loci on LBP exactly.
- plot.out/save2tmp modernisation deferred (API change).
- gl.colors banner ("Starting gl.colors...") leaks from the default
  plot_colors argument evaluation even at verbose = 0 — a gl.colors
  defect, noted for that function's review.

## Coverage

Characterization baseline: `tests/testthat/test-gl.report.hwe.R` — 13
assertions: direction partition, Het.exp recomputation, min.hobs screen,
invisible return + input untouched, ChiSquare vs HWChisq, verbosity
dependence 311/245 (baseline), stub shape + 6-locus LBP equivalence. All
13 pass post-migration, pre-review-fix.

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
  "function": "gl.report.hwe",
  "package": "dartR.base",
  "family_mode": "report",
  "commit": "2739215+migration",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "R1", "severity": "HIGH", "rules": ["VRB"], "loc": "R/gl.report.hwe.r pop-skipping blocks", "status": "proposed"},
    {"id": "R2", "severity": "MEDIUM", "rules": ["spec", "DEP"], "loc": "R/gl.report.hwe.r package checks", "status": "proposed"},
    {"id": "R3", "severity": "LOW", "rules": ["DOC5"], "loc": "R/gl.report.hwe.r @return", "status": "proposed"},
    {"id": "R4", "severity": "LOW", "rules": ["VRB2"], "loc": "R/gl.report.hwe.r alpha check", "status": "proposed"},
    {"id": "R5", "severity": "LOW", "rules": ["DOC1", "DOC2", "DOC7"], "loc": "R/gl.report.hwe.r header", "status": "proposed"}
  ],
  "datasets": ["testset.gl", "LBP"],
  "baseline_test": "tests/testthat/test-gl.report.hwe.R",
  "pr": null
}
```
