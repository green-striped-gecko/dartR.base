# Review: utils.dist.binary (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.dist.binary. Reviewed in the
  distance wave (five files, one approval round, per-function PRs),
  under the standing member directive that utility functions are not
  for end users.
- Datasets: constructed 10x30 SNP fixture with missing data (plus its
  allele-relabelled mirror), testset.gs, raw-vector sequence fixtures,
  brute-force references.
- Family mode: analysis (distance/resampling kernel).
- Checks skipped: Google Group not searched (not available: no
  browser session); dartr2shiny not present in the local workspace.

## Verdicts

- **Standards: FAIL** — the scale warning prints at verbose 0
  (verified); "N11" leftovers in the type doc and two report
  messages; an indented @details tag; large commented-out reference
  loops.
- **Spec: FAIL (one gap)** — Jaccard and Sorensen verified exact
  against manual N-count recomputations; Euclidean/simple standard;
  but bray-curtis is documented AND implemented in the C++ yet
  missing from the validation list, so it silently falls back to
  simple matching (verified). Bray-Curtis and Sorensen are the same
  formula here (the docs already note the synonymy).

## Findings

### D5 — bray-curtis accepted (MEDIUM) [escalation: previously silently simple]

Added to the validation list; requests that used it now get the
documented Sorensen/Dice dissimilarity instead of simple matching.

### D6 — Scale warning gated (LOW)

verbose >= 1.

### D7 — Docs and tidy (LOW)

"N11" leftovers removed; @details indentation; dead commented loops
removed; @keywords internal (stays unexported).

## Coverage

test-utils.dist.binary.R — 4 assertions: Jaccard/Sorensen anchors,
bray-curtis rejection (baseline), scale-warning leak (baseline). All
pass pre-fix.

## Approval

All findings approved via the approval boxes (2026-09-01); utils.hamming
disposition: REMOVE.

## Outcome

All findings applied (D5 bray-curtis accepted, D6 warning gated, D7 docs/tidy). Suite: 4/4; Jaccard/Sorensen anchors unchanged. gl.dist.ind caller smoke clean on SilicoDArT. PR #316.

```json
{"function": "utils.dist.binary", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [{"id": "D5", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/utils.dist.binary.r method list", "status": "applied"},
  {"id": "D6", "severity": "LOW", "rules": ["VRB"], "loc": "R/utils.dist.binary.r scale warning", "status": "applied"},
  {"id": "D7", "severity": "LOW", "rules": ["DOC", "STY"], "loc": "R/utils.dist.binary.r docs", "status": "applied"}],
 "datasets": ["constructed", "testset.gs"],
 "baseline_test": "tests/testthat/test-utils.dist.binary.R",
 "pr": 316}
```
