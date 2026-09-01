# Review: utils.het.pop (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.het.pop. Reviewed in the
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

- **Standards: FAIL** — near-empty doc block (no description/details);
  the file defines a second top-level helper (ind.count) with a
  generic name; `t` shadows base::t.
- **Spec: PASS (with a methods note)** — the computation verifies
  exactly against its own definition: unbiased expected
  heterozygosity 2n/(2n-1)*(1-p^2-q^2) averaged over loci, with n the
  MEAN number of genotyped individuals across loci (not per-locus n).
  The mean-n shortcut differs from a per-locus correction (0.014200
  vs 0.014841 on the fixture) — recorded as a methods note for the
  custodian, NOT changed here: the heterozygosity outputs were
  settled with collaborators in the gl.report.heterozygosity review
  (PR #291).

## Findings

### K5 — Docs, names, tidy (LOW)

Return vector gains population names; doc block filled out
(description, details incl. the mean-n definition); @keywords
internal (stays unexported); `t` renamed.

## Report notes

- Mean-n vs per-locus-n correction: methods decision for the
  custodian (would change published numbers).

## Coverage

test-utils.het.pop.R — 2 assertions (computation anchor, unnamed
baseline). All pass pre-fix.

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
{"function": "utils.het.pop", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "PASS",
 "findings": [{"id": "K5", "severity": "LOW", "rules": ["DOC", "STY"], "loc": "R/utils.het.pop.r", "status": "applied"}],
 "datasets": ["testset.gl", "constructed", "platypus.gl"],
 "baseline_test": "tests/testthat/test-utils.het.pop.R",
 "pr": null}
```
