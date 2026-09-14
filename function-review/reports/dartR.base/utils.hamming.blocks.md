# Review: utils.hamming.blocks (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.hamming.blocks. Reviewed in the
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

- **Standards: PASS** — modern, documented, @keywords internal
  already, session-cached compilation, caller guards documented.
- **Spec: PASS** — the dedup block-hashing detector matches a
  brute-force worst-to-best scan exactly on a 40-sequence fixture,
  and the pairwise mismatch counts are exact (both verified). No
  changes proposed; the only note is cosmetic (file
  utils.hamming.blocks.r hosts function utils.hamming.engine).

## Findings

None.

## Coverage

test-utils.hamming.blocks.R — 5 assertions (dedup vs brute force,
cap flag, pairwise exactness). No flips.

## Approval

All findings approved via the approval boxes (2026-09-01); utils.hamming
disposition: REMOVE.

## Outcome

No changes (PASS/PASS). Characterization tests added: 5/5 (dedup vs brute force, pairwise exactness). PR #318.

```json
{"function": "utils.hamming.blocks", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "PASS", "verdict_spec": "PASS",
 "findings": [],
 "datasets": ["constructed", "testset.gs"],
 "baseline_test": "tests/testthat/test-utils.hamming.blocks.R",
 "pr": 318}
```
