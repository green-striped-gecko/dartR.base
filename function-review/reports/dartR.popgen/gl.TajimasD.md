# Review: gl.TajimasD + utils.get.allele.freq (dartR.popgen)
- Family mode: analysis (helper: utility)
- Scope: one review, one report, one PR for the function and the helper it
  calls; two manifest rows
- Custodian: Ching Ching Lau (gl.TajimasD, author Renee Catullo); Arthur
  Georges (utils.get.allele.freq). STY5: this report is the discussion
  record; changes approved by Luis
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 88b9c22 (origin/dev)
- Datasets: possums.gl populations A and B at full call rate (30 + 30
  individuals, 200 loci); platypus.gl SEVERN_ABOVE and SEVERN_BELOW (6.8%
  missing calls); testset.gs
- External software: ms and sample_stats (`~/programs`); pegas 
  `tajima.test` as the independent reference for D
- Baseline: tests/testthat/test-gl.TajimasD.R (8 tests, 16 expectations; 2
  need MS_DIR)

## Verdict

**Standards: Needs work** — the helper prints at every verbosity and fails
on a genlight without dartR flags, SilicoDArT data is accepted, `cleanup`
never runs, and ms problems surface as ms's usage text or unrelated R
errors.
**Spec: Rework** — D itself is right (matches pegas), but all three
p-values assume the SNPs lie on one non-recombining locus. DArT SNPs are
largely unlinked, so the p-values are far too large and real departures
from neutrality are reported as non-significant; the ms simulation also
uses the number of individuals as the number of sequences.

What works well: D, pi and S match `pegas::tajima.test` on complete data
(D to 4-5 decimals); per-locus sample sizes are used for pi, so missing
calls do not misclassify loci as in `gl.sfs` before #106.

## Findings

**F1 [HIGH, confidence: high] — p-values assume one non-recombining locus (DOC5)**
`R/gl.TajimasD.r:143-157,216-279` — the variance in D's denominator, the
beta and normal p-values (Tajima 1989) and the ms null (`ms nsam rep -s S`,
one locus with S linked sites) all assume complete linkage. For unlinked
sites the variance of pi - S/a1 has no S(S-1) covariance term, so D's null
spread is much narrower.
Failure scenario (probe): possums.gl population A (60 sequences, S = 176):
the null standard deviation of D is 0.87 when simulated as the code does
and 0.20 for 176 unlinked sites. A D of 0.5 has a normal p-value of 0.62
but is well outside the unlinked null (|D| > 0.5 in under 2% of
simulations at sd 0.20). Users testing selection or demography with DArT
data get non-significant results for real departures. In addition, DArT
loci are ascertained as polymorphic, which shifts D upwards (D = 2.6 and
3.1 in the possums populations); the documentation does not say so.
Proposed change: simulate the null as S unlinked sites (`ms 2N rep*S -s
1`, one segregating site per locus, grouped into `rep` replicates of S
sites, D computed in R), report its p-value as `sim_pval`, and document
that `Pval.normal`/`Pval.beta` assume a single non-recombining locus and
that ascertainment inflates D. `sample_stats` is no longer needed.
**Consequence: `sim_pval` changes (it now tests against an unlinked
null); `sample_stats` is no longer required.**

**F2 [MEDIUM, confidence: high] — ms is given the number of individuals as the number of sequences (DOC5)**
`R/gl.TajimasD.r:227,243` — `nsam` is `N`, the mean number of genotyped
individuals, while D is computed on 2N sequences; with missing data `N` is
non-integer (21.5) and ms truncates it.
Failure scenario (probe): population A, 30 individuals: ms simulates 30
sequences instead of 60 (null sd 0.83 instead of 0.87).
Proposed change: `nsam = round(2N)` (part of change 1).

**F3 [MEDIUM, confidence: high] — `cleanup` never runs (DOC5)**
`R/gl.TajimasD.r:361-367` — the cleanup code comes after `return()`; the
working directory is switched to the temporary folder for the whole
function even without ms.
Failure scenario (baseline test 7): each ms run leaves its folder, with
copies of ms and sample_stats, in `tempdir()`.
Proposed change: remove the folder before returning when `cleanup = TRUE`;
change directory only around the ms calls.

**F4 [MEDIUM, confidence: high] — SilicoDArT data is accepted (FS4)**
`R/gl.TajimasD.r:74` — `datatype` is never used; the helper recodes
presence/absence as homozygotes.
Failure scenario (baseline test 3): `testset.gs` returns a D for each of 29
populations.
Proposed change: stop on SilicoDArT input.
**Consequence: SilicoDArT input errors.**

**F5 [MEDIUM, confidence: high] — utils.get.allele.freq: fails without dartR flags; prints regardless of verbosity (DAT5, VRB1)**
`R/utils.get.allele.freq.r:24,38,53,67` — `x@other$loc.metrics.flags$monomorphs`
is read unguarded; `verbose = 2` is a fixed default (no
`gl.check.verbosity`) and `gl.TajimasD` calls the helper without passing
its own `verbose`; the messages name `gl.percent.freq`.
Failure scenario (baseline tests 2, 4, 6): a genlight without
`loc.metrics.flags` stops with "argument is of length zero";
`gl.TajimasD(..., verbose = 0)` prints five helper lines.
Proposed change: guard the flag check; `verbose = NULL` with
`gl.check.verbosity`; `gl.TajimasD` passes its verbosity; correct the
messages.

**F6 [LOW, confidence: high] — D computed from rounded, averaged inputs (DOC5)**
`R/gl.TajimasD.r:94,130-143` — allele frequencies come from the helper
rounded to 4 decimals; the constants a1, a2, e1, e2 use the mean sample
size (truncated in `1:(n - 1)` when non-integer), and Watterson's term is
S/a1(mean n) for all sites.
Failure scenario: complete data: D differs from pegas in the 5th decimal
(2.592523 vs 2.592547). platypus SEVERN_BELOW (missing calls): S/a1(mean
n) = 113.14 against 112.10 from per-site sample sizes.
Proposed change: compute p from the allele counts and number of called
sequences per locus, and Watterson's term as the sum over segregating
sites of 1/a1(n_i); keep the constants for the variance at the mean
sample size and document the approximation.
**Consequence: D changes in the 5th decimal on complete data and by up to
about 1% with missing calls.**

**F7 [MEDIUM, confidence: high] — ms problems surface as unrelated errors (FS5, VRB2)**
`R/gl.TajimasD.r:173-213` — a wrong `ms.path` fails in `normalizePath`;
`rep = NULL` passes an empty count to ms, which prints its usage text and
the function then fails with "length of 'dimnames' [2] not equal to array
extent"; missing binaries print with `cat()` then `stop()` with no
message.
Failure scenario (baseline test 8): `ms.path` set, `rep` left NULL.
Proposed change: require `rep` when `ms.path` is set; clear
`stop(error())` messages for the path and binaries (pointing to
`gl.download.binary("ms")`); check the ms exit status.

**F8 [LOW, confidence: high] — verbosity default and documentation (FS2, DOC1, DOC5, DOC6 (proposed rule), DOC7 (proposed rule))**
`R/gl.TajimasD.r:1-63` — `verbose = 2` ignores `gl.set.verbosity`; no
`@family`; `@author` has no Author(s) label; `@return` and `@export`
twice; a stray `#' #'`; non-ASCII dashes; `theta_per_site` is pi divided by
the number of loci (not base pairs); `N` is the mean number of genotyped
individuals; `plot.dir` default documented as the working directory.
`plot_theme` keeps its name (dartR Shiny passes it).
Failure scenario: a user reads `theta_per_site` as per base pair.
Proposed change: `verbose = NULL`; fix the items above.

## Proposed changes

1. ms null of S unlinked sites with 2N sequences; document the linkage
   assumption of `Pval.normal`/`Pval.beta` and ascertainment (F1, F2).
   **Consequence: `sim_pval` changes; `sample_stats` no longer
   required.**
2. `cleanup` runs; directory changed only around ms (F3).
3. Reject SilicoDArT input (F4).
   **Consequence: SilicoDArT input errors.**
4. utils.get.allele.freq: guard the flag, honour verbosity, correct
   messages; gl.TajimasD passes `verbose` (F5).
5. D from exact counts with per-site sample sizes in Watterson's term
   (F6).
   **Consequence: D changes in the 5th decimal on complete data and by up
   to about 1% with missing calls.**
6. ms input checks and clear errors (F7).
7. `verbose = NULL` and documentation (F8).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run
- Spec: D, pi, S, p-values against `pegas::tajima.test` (complete data);
  null spread of D simulated linked (as the code) and unlinked with ms;
  missing data on platypus.gl; SilicoDArT; no flags; one population;
  `plot.file`; wrong `ms.path`; `ms.path` without `rep` — run
- Callers (API3): dartR Shiny (`rep`, `ms.path`, `simulation.out`,
  `plot_theme`); no sibling dartR.* package calls either function — run
- Windows (`shell()` path): SKIPPED — not run
- Recombination between some DArT loci (partial linkage): not modelled;
  the unlinked null is the closer approximation for DArT data, documented
- DAT6 (FBM): SKIPPED — no FBM fixture; the helper densifies per
  population
- Google Group / GitHub issues: not searched in this session

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | addendum: null sampled in R, not ms (see Outcome) |
| 2 | approved | Luis |  |
| 3 | approved | Luis |  |
| 4 | approved | Luis |  |
| 5 | approved | Luis |  |
| 6 | approved | Luis |  |
| 7 | approved | Luis |  |

## Outcome

Branch `review-tajimasd` (from origin/dev 88b9c22). Evidence:
`tests/testthat/test-gl.TajimasD.R`, 8 tests, 22 expectations, all passing
(no external software). Every changed assertion is tagged with its approved
change.

- Addendum to change 1 (approved by Luis in Phase C): simulating S
  unlinked sites with `ms -s 1` prints nsam lines per site (about 320
  million lines for rep = 1000, S = 5000) and ms's fixed-site conditioning
  departs from the neutral 1/k distribution (chi-square p = 0.001 on
  20,000 loci). The null is sampled in R instead: derived counts with
  P(k) proportional to 1/k from round(2N) sequences, S sites per replicate.
  ms and sample_stats are no longer needed; `ms.path` is accepted and
  ignored with a warning; `cleanup` has nothing to clean.
- 1 applied: null sd 0.20 for 60 sequences and 176 sites (0.195 from ms
  with unlinked loci); on simulated neutral unlinked data sim_pval < 0.05
  in 5% of 200 datasets; seeded runs reproducible, and the user's random
  number stream is restored. Documentation of the linkage assumption of
  Pval.normal/Pval.beta and of ascertainment added.
- 2 applied: no temporary folder is created any more (moot after the
  addendum).
- 3 applied: SilicoDArT errors.
- 4 applied: helper works without flags, `verbose = NULL`, messages
  corrected; `gl.TajimasD` calls it with `verbose = 0`; `verbose = 0`
  prints nothing. Helper output identical to before.
- 5 applied: D equals pegas to 5.8e-15 on complete data (was 2.5e-5).
  **Correction to the approved consequence:** with missing calls D changed
  by 1.1% (SEVERN_ABOVE) and 5.3% (SEVERN_BELOW), not "up to about 1%";
  Watterson's term moves under 1%, but D divides the difference pi - theta_W.
  Pval.beta differs from pegas by up to 9e-6 (pegas's own beta bounds).
- 6 applied: `rep` required with `ms.path`; `rep < 1` rejected.
- 7 applied: `verbose = NULL`, `plot.dir` resolved with `gl.check.wd`,
  documentation rewritten, `plot_theme` name kept for dartR Shiny.
  `devtools::document()` also rewrote stale Rd files on dev; reverted.
  The "demographic history" family now links gl.TajimasD from gl.sfs,
  gl.run.epos and gl.run.stairway2.
- `devtools::check()`: 0 errors; install WARNING (packages built under
  R 4.4.3) and 3 NOTEs present before this change.
- NEWS entry added. Caller: dartR Shiny (passes `ms.path`, `rep`,
  `simulation.out`, `plot_theme`; still works, `ms.path` ignored).
- PR: #108

## Machine block

```json
{
  "function": "gl.TajimasD",
  "functions_in_scope": ["gl.TajimasD", "utils.get.allele.freq"],
  "package": "dartR.popgen",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "88b9c22",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS4", "status": "approved", "change": 3},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DAT5", "status": "approved", "change": 4},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 5},
    {"id": "F7", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 6},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "FS2", "status": "approved", "change": 7}
  ],
  "coverage_skipped": ["Windows shell() path not run", "partial linkage not modelled", "DAT6: no FBM fixture", "Google Group / issues not searched"],
  "status": "pr-open",
  "pr": 108
}
```
