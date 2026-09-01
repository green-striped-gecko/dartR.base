# Review: utils.allelic.richness (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.allelic.richness. Reviewed in the
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

- **Standards: FAIL** — @return documents "calling function name"
  (wrong); param docs are placeholders ("indices indices"); the
  `allelicRichness <- function(){}` R CMD check hack is undocumented.
- **Spec: PASS** — the hypergeometric rarefaction arithmetic was
  validated against an independent recomputation in the
  gl.report.allelerich review (PR #286); kernel re-anchored here
  (fixture value 1.2; monomorphic input returns exactly 1).
  Rcpp::cppFunction compiles once per session (cached), noted as a
  performance consideration for bootstrap-heavy callers.

## Findings

### K7 — Docs (LOW)

@return corrected (mean rarefied allelic richness across sites);
params documented (df orientation, boot::boot usage); the check hack
commented; @keywords internal added (retains its NAMESPACE export to
avoid hand-editing NAMESPACE; full un-export deferred - no known
external callers, dartr2shiny unverified).

## Coverage

test-utils.allelic.richness.R — 2 assertions (fixture anchor,
monomorphic edge). No flips.

## Approval

All findings approved via the approval boxes (2026-09-01).

## Outcome

All findings applied. Kernel suite: 30/30 assertions pass; flips map
to approved findings only. utils.basic.stats verified to match
hierfstat::basic.stats exactly on complete data, 25% missingness,
monomorphic loci, per-population absent loci, and the real testset
subset; gl.report.fstat runs on the fixed kernel (its own one-line
verbose-0 leak is noted for its review). PR #313.

```json
{"function": "utils.allelic.richness", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "PASS",
 "findings": [{"id": "K7", "severity": "LOW", "rules": ["DOC"], "loc": "R/utils.allelic.richness.r header", "status": "applied"}],
 "datasets": ["testset.gl", "constructed", "platypus.gl"],
 "baseline_test": "tests/testthat/test-utils.allelic.richness.R",
 "pr": 313}
```
