# Review: utils.n.var.invariant (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.n.var.invariant. Reviewed in the
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

- **Standards: FAIL** — the secondaries-in-history warning prints at
  verbose 0 (verified); the `@details` tag is indented (fragile
  roxygen); a stray developer aside in the code; verbose doc says
  "[default NULL]".
- **Spec: PASS** — verified on platypus.gl: n.variant equals the
  CloneID multiplicity for every locus and n.invariant equals
  lenTrimSeq - n.variant; TrimmedSequence requirement enforced with
  an informative stop; silent at verbose 0 on a clean history. The
  CloneID fallback assumes 3 pipe-separated AlleleID fields (holds
  for DArT data; noted). History append retained (the function
  returns a modified object).

## Findings

### K8 — Warning gate and docs (LOW)

Secondaries warning gated at verbose >= 1; @details indentation
fixed; verbose doc standardised; developer aside removed; @keywords
internal (stays unexported).

## Coverage

test-utils.n.var.invariant.R — 6 assertions (count verification,
history-is-call, warning-leak baseline). All pass pre-fix.

## Approval

All findings approved via the approval boxes (2026-09-01).

## Outcome

All findings applied. Kernel suite: 30/30 assertions pass; flips map
to approved findings only. utils.basic.stats verified to match
hierfstat::basic.stats exactly on complete data, 25% missingness,
monomorphic loci, per-population absent loci, and the real testset
subset; gl.report.fstat runs on the fixed kernel (its own one-line
verbose-0 leak is noted for its review). PR #314.

```json
{"function": "utils.n.var.invariant", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "PASS",
 "findings": [{"id": "K8", "severity": "LOW", "rules": ["VRB", "DOC"], "loc": "R/utils.n.var.invariant.r", "status": "applied"}],
 "datasets": ["testset.gl", "constructed", "platypus.gl"],
 "baseline_test": "tests/testthat/test-utils.n.var.invariant.R",
 "pr": 314}
```
