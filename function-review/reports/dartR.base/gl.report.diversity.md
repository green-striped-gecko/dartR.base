# Review: gl.report.diversity (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-gl.report.diversity.
- Datasets: testset.gl (30 pops), testset.gs, pair subset via
  gl.keep.pop, constructed 2-population object with a per-population
  all-NA locus.
- Family mode: report (analysis lens on the entropy/Hill-number
  computations, Sherwin et al. 2017/2021).
- Author: Bernd Gruber; Contributors: William B. Sherwin, Alexander
  Sentinella.
- Checks skipped: Google Group not searched (not available: no browser
  session).

## Verdicts

- **Standards: FAIL** — the kable tables print unconditionally (356
  lines at verbose 0 on testset); the gl.colors banner leaks at
  verbose 0 when plot.display = TRUE; plot.file is documented twice
  with contradictory text; a large commented-out SAVE INTERMEDIATES
  block; the "number of missing loci" comment describes a non-missing
  count; `((...) > 0) + 1 - 1` logical coercion; the example passes
  table = FALSE where the doc offers 'N'.
- **Spec: FAIL** — (1) The q = 1 (Shannon) computations mishandle
  per-population all-missing loci: they enter the alpha mean as ZERO
  diversity (verified: reported 0.5844 equals the zeros-version, the
  correct exclusion gives 0.6151 on the constructed object), and the
  full-length internal vector misaligns with the beta indexing by
  logical recycling whenever a population carries an all-NA locus —
  which is every one of testset's 30 populations (up to 36 such loci).
  q = 0 handles the case via na.rm and q = 2 by exclusion, so q = 1 is
  also internally inconsistent. (2) one_H_beta and two_H_beta pool ALL
  populations rather than the pair, and two_H_beta's correction factor
  uses the dataset-wide npops: the same pair's beta reads 0.0396 in
  the full dataset vs 0.0103 pair-only (q = 1) and 0.0211 vs 0.0108
  (q = 2) — a pairwise measure that changes when unrelated populations
  are added (the same defect class as the Chao estimator fixed in
  PR #299). zero_H_beta is already pair-based (verified identical).
  (3) plot.file with plot.display = FALSE crashes ("object 'p3' not
  found"). (4) SilicoDArT input runs silently through diploid 0/1/2
  arithmetic and produces meaningless indices (the title says "for
  SNPs"; the datatype is checked and ignored). (5) A mistyped table
  argument silently suppresses all tables even at verbose 2. The q = 0
  and q = 2 alpha computations are verified exact against independent
  recomputations.

## Findings

### F1 — Shannon (q=1) missing-locus handling (HIGH, confidence: certain) [escalation: numbers change]

Location: one_H_alpha_es (dummys keeps zeros for all-missing loci);
one_H_beta_es / two_H_beta_es indexing.

Proposed change: drop per-population all-missing loci from the Shannon
per-locus vector (as q = 2 already does). This corrects the alpha
deflation and simultaneously aligns the beta indexing (the logical
subscripts then match the vector length, as they already do for
q = 2). one_H_alpha, one_H_alpha_sd, one_D_alpha and one_H/D_beta
values change for any dataset with per-population missing loci — on
testset, all populations.

### F2 — Pairwise betas pooled over all populations (HIGH, confidence: certain) [escalation: numbers change]

Location: one_H_beta (one_H_alpha_all from the full x); two_H_beta
(two_H_alpha_all from the full x; npops/(npops-1) with dataset npops).

Proposed change: compute the pooled term from the two populations of
the pair, and use the pairwise correction factor 2 (= npops/(npops-1)
with npops = 2). Pairwise beta values then depend only on the pair
(zero_H_beta already behaves this way). Values change for every
dataset with more than two populations.

### F3 — Tables ungated and table argument unvalidated (MEDIUM, confidence: certain) [escalation: user-visible behavior]

Proposed change: print the tables only at verbose >= 1 (the table
argument remains the format selector); validate table against
c("D", "H", "DH", "HD", "N") with an informative stop (currently a
typo — or the old example's table = FALSE — silently suppresses all
tables); the example switched to table = "N".

### F4 — plot.file crashes without plot.display (MEDIUM, confidence: certain)

Proposed change: exists-guard on p3 with a gated note (the established
idiom).

### F5 — SilicoDArT produces silent garbage (MEDIUM, confidence: certain) [escalation: silico now fatal]

Proposed change: utils.check.datatype(x, accept = "SNP") — the
entropy formulas assume diploid 0/1/2 genotypes; the @param x text
drops the presence/absence claim (the @title already says SNPs).

### F6 — gl.colors banner leaks at verbose 0 (LOW, confidence: certain)

Proposed change: the caller-side hand-patch used across the campaign —
default plot.colors.pop = gl.colors("dis", verbose = 0).

### F7 — Tidy and docs (LOW, confidence: certain)

Duplicate contradictory @param plot.file removed; verbose doc wording
standardised; "20 tables" corrected (10 printed tables; the return
list has 20 components); the "number of missing loci" comment
corrected; logical coercion tidied; the commented-out SAVE
INTERMEDIATES block removed; example table = "N".

## Report notes (other functions / not fixed here)

- gl.filter.allna is called internally with verbose = verbose — its
  own gating is correct, so no leak; noted only.
- The author TODO ("adjust calculation of betas for population sizes")
  is retained above the function.

## Coverage

`tests/testthat/test-gl.report.diversity.R` — 12 assertions: q0/q2
alpha vs independent recomputation (anchor), verbose-0 leak
(baseline), Shannon zeros-bias (baseline), Shannon beta misalignment
(baseline), beta pair-dependence with zero_H control (baseline), p3
crash (baseline), silico runs (baseline), table typo silent
(baseline), gl.colors banner (baseline). All 12 pass pre-fix.

## Approval

All findings approved via the approval boxes (2026-09-01).

## Outcome

All findings applied. Suite: 12/12 pass; flips map to F1 (Shannon
alpha equals the exclusion recomputation; beta equals the pair-pooled
recomputation), F2 (betas identical between full-dataset and pair-only
runs; zero_H_beta control unchanged), F3 (silent at verbose 0;
informative stop on an invalid table), F4 (plot.file without
plot.display completes), F5 (SilicoDArT rejected), F6 (no banner at
verbose 0 when plotting). End-to-end verbose 3 on a 3-population
subset clean with sane Hill numbers. PR #301.

```json
{"function": "gl.report.diversity", "package": "dartR.base", "family_mode": "report",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [
  {"id": "F1", "severity": "HIGH", "rules": ["spec"], "loc": "R/gl.report.diversity.r one_H", "status": "applied"},
  {"id": "F2", "severity": "HIGH", "rules": ["spec"], "loc": "R/gl.report.diversity.r betas", "status": "applied"},
  {"id": "F3", "severity": "MEDIUM", "rules": ["VRB", "spec"], "loc": "R/gl.report.diversity.r tables", "status": "applied"},
  {"id": "F4", "severity": "MEDIUM", "rules": ["spec", "PLT"], "loc": "R/gl.report.diversity.r plot.file", "status": "applied"},
  {"id": "F5", "severity": "MEDIUM", "rules": ["spec", "DAT"], "loc": "R/gl.report.diversity.r datatype", "status": "applied"},
  {"id": "F6", "severity": "LOW", "rules": ["VRB"], "loc": "R/gl.report.diversity.r palette default", "status": "applied"},
  {"id": "F7", "severity": "LOW", "rules": ["DOC", "STY"], "loc": "R/gl.report.diversity.r", "status": "applied"}],
 "datasets": ["testset.gl", "testset.gs", "constructed"],
 "baseline_test": "tests/testthat/test-gl.report.diversity.R",
 "pr": 301}
```
