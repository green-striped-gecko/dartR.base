# Review: gl.allele.freq (dartR.base)

## Provenance

- Model: Claude Fable 5 (claude-fable-5, Claude Code) via dartr-dev agent;
  Skill: dartr-function-review v2.0.0; Base: upstream/dev at ddaed27
  (`git diff upstream/dev -- R/gl.allele.freq.r` empty — the loaded code
  is the reviewed code); working branch integration-local at ed99203.
- Reviewed in the population-distance chain wave (four functions
  reviewed coherently: gl.allele.freq, gl.dist.pop, gl.dist.ind,
  utils.collapse.matrix v2), because gl.allele.freq feeds every SNP
  method in gl.dist.pop.
- Datasets: testset.gl, testset.gs, constructed plain-genlight fixture.
- Family mode: analysis (frequency table generator; unmatched report).
- Baseline: tests/testthat/test-gl.allele.freq.R (10 tests, all pass at
  the reviewed state).
- Checks skipped: Google Group not searched (not available: no browser
  session); FBM path (DAT6) not exercised (no FBM fixture in this wave —
  the `@examples` gen2fbm line was not run).

## Verdicts

**Standards: Needs work** — preamble, gating, and silence at `verbose = 0`
all conform (verified empirically: zero lines); the defects are an
unguarded flags access, an unvalidated `by` argument, and documentation
drift.
**Spec: Needs work** — `by = 'popxloc'` is cell-exact against hand
computation (every population x locus cell on a 4-population subset,
plus `sum`/`nobs`/`nmissing`), but the `by = 'loc'` branch breaks two of
the function's own contracts: it ignores `percent = TRUE`, and for
SilicoDArT it returns half the presence frequency.

Empirical NA policy (this is the chain's reference statement): each
population x locus cell is `mean(genotypes, na.rm = TRUE)/2` — missing
genotypes are dropped within the cell; a cell with no scored genotypes
is NaN (545 of 7650 cells on testset.gl). `by = 'pop'` then averages
cells across loci with NaN cells dropped (`na.rm = TRUE` in the
aggregate). Frequencies are rounded to 2 dp on the percentage scale
before any output or aggregation (max error 1.1e-05 on the proportion
scale — benign, but it propagates into every gl.dist.pop method).

Cross-note (gl.tree.nj, PR #370): the PR's claim that `na.rm = TRUE` is
"the policy gl.allele.freq/gl.dist.pop use" is correct at the frequency
level — verified here cell-by-cell.

## Findings

**F1 [HIGH, confidence: high] — `by = 'loc'` ignores `percent = TRUE` (DOC5 (proposed rule))**
`R/gl.allele.freq.r:194` — the branch overwrites `frequency` with
`colMeans(as.matrix(x), na.rm = TRUE)/2` after the percent scaling has
been applied, so the output is always a proportion.
Failure scenario: `gl.allele.freq(testset.gl, percent = TRUE, by = 'loc')`
returns frequencies with maximum 1.0 (verified); a user comparing
`by = 'loc'` with `by = 'popxloc'` output at `percent = TRUE` sees values
differing by a factor of 100.
Proposed change: scale the overwritten column by 100 when
`percent = TRUE` (or compute it before the percent branch).

**F2 [HIGH, confidence: high] — SilicoDArT `by = 'loc'` and `simple = TRUE` return half the presence frequency (DOC5 (proposed rule), DAT1 lens)**
`R/gl.allele.freq.r:93-96,194` — the SilicoDArT 1-to-2 recode is applied
to the per-population split list only; the `by = 'loc'` overwrite uses
the raw 0/1 matrix and still divides by 2 (the SNP ploidy divisor).
Failure scenario: on testset.gs, `by = 'popxloc'` returns the presence
percentage but `by = 'loc'` (and therefore `simple = TRUE`, which forces
`by = 'loc'`) returns exactly half the presence proportion (verified:
equal to `colMeans/2`, not `colMeans`); `alf1`/`alf2` from
`simple = TRUE` are wrong for Tag P/A data.
Proposed change: skip the /2 divisor (or apply the recode) when
`datatype == "SilicoDArT"` in the `by = 'loc'` overwrite.
**Consequence: numerical output changes for SilicoDArT `by = 'loc'` and
`simple = TRUE` callers.**

**F3 [MEDIUM, confidence: high] — unguarded flags access breaks plain genlight input (DAT5)**
`R/gl.allele.freq.r:70` — `x@other$loc.metrics.flags$monomorphs` is
indexed with no existence check.
Failure scenario: a genlight not built by dartR (no
`loc.metrics.flags`) fails with "argument is of length zero" (verified)
instead of a compliance message.
Proposed change: guard with `is.null()` (or route through
`gl.compliance.check`), per the DAT5 idiom.

**F4 [MEDIUM, confidence: high] — `by` is never validated (FS5)**
`R/gl.allele.freq.r:172-213` — the if/else chain treats any string that
is not `'pop'` or `'loc'` as `'popxloc'`.
Failure scenario: `by = 'population'` (or any typo) silently returns the
full population x locus table instead of erroring — the same silent-
fallback class as gl.tree.nj's "ugpma" precedent.
Proposed change: validate `by` against `c('pop','loc','popxloc')` and
stop (or warn and default) on anything else.

**F5 [LOW, confidence: high] — documentation drift (DOC5 (proposed rule), DOC1)**
`R/gl.allele.freq.r:38-39,84-106` — `@return` says "A matrix"; the
function returns a data.frame (verified). The `verbose >= 2` progress
messages say "Tag allele frequencies" for SNP data ("Tag" is the P/A
vocabulary). `simple = TRUE` silently overrides user-supplied `percent`
and `by` (verified identical output) without a note in `@param simple`.
`by = 'loc'` averages `nobs`/`nmissing`/`n` across populations (mean
per-population nobs, verified 8.3 vs total 249), while the doc says
"averaged across individuals".
Proposed change: docs-only — correct `@return`, the message text, the
`simple` and `by` parameter descriptions.

**F6 [INFO, confidence: high] — 2 dp percentage rounding is baked into every output**
`R/gl.allele.freq.r:129` — frequencies are rounded to 2 dp on the
percentage scale at source; `percent = FALSE` is exactly the rounded
percentage /100 (verified). Max distortion 1.1e-05 as a proportion.
Not proposed for change (behaviour is long-standing and benign); noted
because the chain review verifies gl.dist.pop against these rounded
values.

## Proposed changes

1. Scale the `by = 'loc'` frequency column when `percent = TRUE` (F1).
   **Consequence: numerical output changes for `percent = TRUE, by='loc'`
   callers (values x100).**
2. Correct the SilicoDArT divisor in the `by = 'loc'` overwrite (F2).
   **Consequence: numerical output changes for SilicoDArT `by='loc'` and
   `simple=TRUE` callers (values x2).**
3. Guard the monomorphs-flag access for non-dartR genlights (F3).
4. Validate `by` (F4).
5. Docs-only corrections: `@return`, progress-message wording, `simple`
   override note, `by='loc'` column semantics (F5).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT (n/a — no plot), STY — run.
- Spec, cell-level: every popxloc cell vs hand computation from
  `as.matrix` on a 4-pop subset (SNP) and 3-pop subset (SilicoDArT);
  `sum`/`nobs`/`nmissing` columns; all-NA cells; percent/proportion
  relation; `by='pop'` aggregation; `by='loc'` overwrite; `simple=TRUE`;
  monomorphic loci retained (verified present in output) — run.
- verbose=0 silence: capture.output empty — run.
- FBM path (DAT6): SKIPPED — no FBM fixture in this wave.
- Google Group search: SKIPPED — not available, no browser session.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur Georges (2026-09-06) | consequence acknowledged: values x100 for `percent = TRUE, by = 'loc'` |
| 2 | approved | Arthur Georges (2026-09-06) | consequence acknowledged: values x2 for SilicoDArT `by = 'loc'`/`simple = TRUE` |
| 3 | approved | Arthur Georges (2026-09-06) | |
| 4 | approved | Arthur Georges (2026-09-06) | |
| 5 | approved | Arthur Georges (2026-09-06) | |

## Outcome

All 5 approved changes applied on branch `review-gl.allele.freq`
(base `ddaed27`, upstream/dev) and submitted as
[PR #374](https://github.com/green-striped-gecko/dartR.base/pull/374),
covering findings F1-F5; F6 is INFO/no
action. Verification:

- `by = 'popxloc'` (percent TRUE/FALSE, SNP and SilicoDArT) and
  `by = 'pop'` outputs are byte-identical to the pre-change state
  (`identical()` on saved snapshots), and `gl.dist.pop()` euclidean/nei
  on `testset.gl` are byte-identical through this branch — the chain
  input is untouched.
- F1 cell-verified: `percent = TRUE, by = 'loc'` equals
  `round(colMeans(as.matrix(x), na.rm = TRUE)/2*100, 4)` exactly
  (max |diff| 0 on testset.gl).
- F2 cell-verified: SilicoDArT `by = 'loc'` equals
  `round(colMeans(as.matrix(x), na.rm = TRUE), 4)` exactly (max |diff|
  0 on testset.gs) — the presence frequency, no longer half of it.
- Baseline characterization test: 34 assertions pass; the four flipped
  expectations are tagged `[approved F1]`-`[approved F4]` (plus two
  `[approved F5]` comment updates); no unexplained diffs.
- `verbose = 0` empirically silent (capture.output empty); `verbose = 3`
  end-to-end run clean on both datatypes.
- Caller grep across the 8 dartRverse clones: every consumer that feeds
  distances uses `by = 'popxloc'` (gl.dist.pop, gl.fdsim, gl.fixed.diff,
  gl2bayescan, gl2snapper, gl2treemix, dartR.popgen gl.nhybrids) —
  unaffected. The two `by = 'loc'` consumers (gl.sim.genotypes,
  gl.report.polyploid_heterozygosity via `simple = TRUE`) are
  SNP-context (the latter `accept = "SNP"`), for which `percent = FALSE`
  behaviour is unchanged. No breaking caller.

Integration probe (all three chain branches loaded together): see the
probe record in gl.dist.pop.md — `gl.allele.freq` popxloc feeds
`gl.dist.pop` identically in the combined state.

```json
{
  "function": "gl.allele.freq",
  "package": "dartR.base",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "ddaed27",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DAT5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "applied", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 5},
    {"id": "F6", "severity": "INFO", "confidence": "high", "rule": "STY1", "status": "noted", "change": null}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "Google Group: no browser session"],
  "baseline_test": "tests/testthat/test-gl.allele.freq.R",
  "status": "pr-open",
  "pr": 374
}
```
