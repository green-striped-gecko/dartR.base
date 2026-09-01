# Review: utils.basic.stats (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.basic.stats. Reviewed in the
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

- **Standards: FAIL** — `cat(error()) + stop()` split leaves the fatal
  error with an EMPTY message (verified); the `@family utilities` tag
  is indented so roxygen treatment is fragile; exported with @examples
  despite the internal-use warning.
- **Spec: FAIL** — the doc claims results "match exactly"
  hierfstat::basic.stats. Verified exactly TRUE for complete data,
  25% random missingness, and monomorphic loci — but FALSE whenever a
  locus is entirely missing from any population: the harmonic mean
  sample size includes the zero count (1/mean(1/np) with np = 0 gives
  mn = 0), so Hs/Ht/Fst and the rest go NaN for that locus and drop
  silently from the overall averages, where hierfstat computes them
  from the populations that have data. On a 3-population testset
  subset: overall Fst 0.2944 (dartR) vs 0.2483 (hierfstat), Fis
  -0.358 vs -0.112. Every testset population carries all-NA loci, so
  this is live on real data.

## Findings

### K1 — Per-pop absent loci break the hierfstat equivalence (HIGH) [escalation: numbers change]

Proposed: harmonic mean over the populations that have data
(`1/mean(1/y[y > 0])`). Verified plan: post-fix output must match
hierfstat exactly on the per-pop-all-NA construction and the testset
subset.

### K2 — Empty fatal error (MEDIUM)

`stop(error(...))` instead of cat + bare stop().

### K3 — Docs/visibility (LOW)

@family indentation fixed; @keywords internal added (STAYS exported —
dartR.popgen::gl.ld.haplotype calls it); single-population edge
(Dstp/Gst_max division by zero) documented as not meaningful.

## Report notes

- The commented "option 2 / Nei 1987" blocks retained (author's
  working notes).

## Coverage

test-utils.basic.stats.R — 6 assertions: exact hierfstat equivalence
on complete data (anchor), per-pop-absent divergence (baseline),
empty error message (baseline). All pass pre-fix.

## Approval

All findings approved via the approval boxes (2026-09-01).

## Addendum (discovered during apply)

Loci carried by a SINGLE population made Dstp = Dst/0 = -Inf, which
poisoned the overall means even after the harmonic-mean fix. Under
K1's approved consequence (exact hierfstat equivalence) Dstp is set
NaN for such loci so they drop from the overall averages, matching
hierfstat. Verified: all nine overall statistics now match hierfstat
exactly on the real testset subset (Fst 0.2483 = 0.2483).

## Outcome

All findings applied. Kernel suite: 30/30 assertions pass; flips map
to approved findings only. utils.basic.stats verified to match
hierfstat::basic.stats exactly on complete data, 25% missingness,
monomorphic loci, per-population absent loci, and the real testset
subset; gl.report.fstat runs on the fixed kernel (its own one-line
verbose-0 leak is noted for its review). PR #309.

```json
{"function": "utils.basic.stats", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [{"id": "K1", "severity": "HIGH", "rules": ["spec"], "loc": "R/utils.basic.stats.r harmonic mean", "status": "applied"},
  {"id": "K2", "severity": "MEDIUM", "rules": ["STY"], "loc": "R/utils.basic.stats.r stop", "status": "applied"},
  {"id": "K3", "severity": "LOW", "rules": ["DOC"], "loc": "R/utils.basic.stats.r docs", "status": "applied"}],
 "datasets": ["testset.gl", "constructed", "platypus.gl"],
 "baseline_test": "tests/testthat/test-utils.basic.stats.R",
 "pr": 309}
```
