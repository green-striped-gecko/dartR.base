# Review: utils.is.fixed (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.is.fixed. Reviewed in the
  numerical-kernel wave (six functions, one approval round,
  per-function PRs), under the standing member directive that utility
  functions are not for end users.
- Datasets: testset.gl subsets; constructed 3-population objects
  (complete, 25% missing, per-population all-NA loci, monomorphic);
  platypus.gl; hierfstat and HardyWeinberg as external references.
- Family mode: analysis (numerical kernel).
- Checks skipped: Google Group not searched (not available: no browser
  session); dartr2shiny not present in the local workspace.

## Verdicts

- **Standards: PASS (with notes)** — trailing commented test calls;
  "to tests" and "tollerance" typos.
- **Spec: PASS (docs wrong)** — the truth table and tolerance
  boundaries verify exactly as documented (95,5 at tloc = 0.05 is
  fixed; 94,6 is not; NA propagates), but @return and @details claim
  TRUE/FALSE where the function returns numeric 1/0/NA — and callers
  rely on the numeric contract (gl.nhybrids tests `> 0`).

## Findings

### K4 — Docs corrected to the numeric contract (LOW)

@return: 1 (fixed difference), 0 (alleles shared), NA (missing
frequency); typos fixed; commented test block tidied; @keywords
internal added (STAYS exported — dartR.popgen::gl.nhybrids calls it).

## Coverage

test-utils.is.fixed.R — 10 assertions (truth table, tolerance
boundaries, NA propagation, numeric type). No flips (docs-only).

## Approval

All findings approved via the approval boxes (2026-09-01).

## Outcome

All findings applied. Kernel suite: 30/30 assertions pass; flips map
to approved findings only. utils.basic.stats verified to match
hierfstat::basic.stats exactly on complete data, 25% missingness,
monomorphic loci, per-population absent loci, and the real testset
subset; gl.report.fstat runs on the fixed kernel (its own one-line
verbose-0 leak is noted for its review). PR recorded below.

```json
{"function": "utils.is.fixed", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "PASS", "verdict_spec": "PASS",
 "findings": [{"id": "K4", "severity": "LOW", "rules": ["DOC"], "loc": "R/utils.is.fixed.r header", "status": "applied"}],
 "datasets": ["testset.gl", "constructed", "platypus.gl"],
 "baseline_test": "tests/testthat/test-utils.is.fixed.R",
 "pr": null}
```
