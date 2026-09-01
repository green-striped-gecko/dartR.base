# Review: utils.dist.ind.snp (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.dist.ind.snp. Reviewed in the
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

- **Standards: FAIL** — the file carries a botched merge: an ungated
  progress cat (its gating if is commented out, the body live -
  1 line leaks at verbose 0), the matrix converted twice, two copies
  of the tail (Returning/Completed print TWICE at verbose >= 2, and
  the dist object is converted dist -> matrix -> dist); large
  commented-out reference loops.
- **Spec: FAIL** — the Simple and Absolute methods are
  reference-allele ASYMMETRIC (verified: relabelling which allele is
  reference changes the distances, max discrepancy 0.167 on the
  fixture; Euclidean is invariant as a control) and contradict their
  own documentation: both-homozygous-REFERENCE pairs are scored as
  sharing NO alleles while both-homozygous-ALTERNATE pairs share
  both. Euclidean and Manhattan/Czekanowski verified exact against
  manual recomputations. The Euclidean scale doc claims the result
  falls in [0,1]; with genotypes 0..2 the scaled maximum is 2.

## Findings

### D1 — Simple/Absolute allele-sharing scoring (HIGH) [escalation: numbers change]

Proposed (doc-conform): per locus shared alleles
S = 2 - |g1 - g2| (0/1/2 as documented); Simple distance
= 1 - mean(S)/2; Absolute rv = (S > 0), distance = proportion of
opposite-homozygote loci. Both become invariant under allele
relabelling. Note: doc-conform Simple is then arithmetically
identical to Czekanowski (recorded in the docs).

### D2 — Merge-artifact cleanup (MEDIUM)

One clean path: single matrix extraction, gated progress message,
single tail (convert once, report once). Fixes the verbose-0 leak
and the doubled messages.

### D3 — Docs (LOW)

Euclidean scale range corrected; dead commented loops removed;
@keywords internal (stays unexported).

## Coverage

test-utils.dist.ind.snp.R — 6 assertions: Euclidean/Manhattan
anchors, relabelling asymmetry (baseline), leak + doubled messages
(baseline). All pass pre-fix.

## Approval

All findings approved via the approval boxes (2026-09-01); utils.hamming
disposition: REMOVE.

## Outcome

All findings applied (D1 allele-sharing scoring, D2 single clean path, D3 docs). Suite: 6/6; Euclidean/Manhattan anchors unchanged; Simple/Absolute now relabel-invariant. gl.dist.ind caller smoke clean on SNP data. One apply-time note: Rcpp's cppFunction signature scanner mis-parses parenthesised text in leading C++ comments (it generated a wrapper named from comment text); commentary moved inside the function body. PR recorded below.

```json
{"function": "utils.dist.ind.snp", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [{"id": "D1", "severity": "HIGH", "rules": ["spec"], "loc": "R/utils.dist.ind.snp.r simple/absolute", "status": "applied"},
  {"id": "D2", "severity": "MEDIUM", "rules": ["VRB", "STY"], "loc": "R/utils.dist.ind.snp.r merge artifacts", "status": "applied"},
  {"id": "D3", "severity": "LOW", "rules": ["DOC"], "loc": "R/utils.dist.ind.snp.r docs", "status": "applied"}],
 "datasets": ["constructed", "testset.gs"],
 "baseline_test": "tests/testthat/test-utils.dist.ind.snp.R",
 "pr": null}
```
