# Review: utils.hwe (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.hwe. Reviewed in the
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

- **Standards: FAIL** — three top-level functions (GenerateSamples,
  CritSam, CritSam_Chi) with no documentation at all; `ncomp`
  assigned and never used in both CritSam variants.
- **Spec: PASS** — GenerateSamples(n) enumerates exactly the
  choose(n+2, 2) genotype compositions summing to n (verified);
  CritSam/CritSam_Chi return critical-sample frequency matrices with
  rows summing to 1 (verified); the HardyWeinberg dependency is
  guarded by the callers (gl.report.hwe / gl.diagnostics.hwe).

## Findings

### K6 — Minimal documentation and tidy (LOW)

Each helper gains a brief @noRd roxygen header (purpose, params,
return); unused `ncomp` removed; unexported status unchanged.

## Coverage

test-utils.hwe.R — 6 assertions (composition enumeration, critical
matrices). No flips.

## Approval

All findings approved via the approval boxes (2026-09-01).

## Outcome

All findings applied. Kernel suite: 30/30 assertions pass; flips map
to approved findings only. utils.basic.stats verified to match
hierfstat::basic.stats exactly on complete data, 25% missingness,
monomorphic loci, per-population absent loci, and the real testset
subset; gl.report.fstat runs on the fixed kernel (its own one-line
verbose-0 leak is noted for its review). PR #312.

```json
{"function": "utils.hwe", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "PASS",
 "findings": [{"id": "K6", "severity": "LOW", "rules": ["DOC", "STY"], "loc": "R/utils.hwe.r", "status": "applied"}],
 "datasets": ["testset.gl", "constructed", "platypus.gl"],
 "baseline_test": "tests/testthat/test-utils.hwe.R",
 "pr": 312}
```
