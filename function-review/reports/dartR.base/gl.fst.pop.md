# Review: gl.fst.pop (dartR.base)

- Family mode: analysis
- Date: 2026-09-08
- Reviewer: Claude Opus 5 (`claude-opus-5`), dartr-function-review v2.0.0
- Package commit: `ddaed27` (upstream/dev; `git diff upstream/dev -- R/gl.fst.pop.r` is empty, so
  the reviewed state is upstream/dev unmodified. Working tree is `integration-local` at `ed99203`
  with unrelated uncommitted test edits.)
- Dependency versions: StAMPP 1.6.3, hierfstat 0.5.11, dartR.data 1.2.5, R 4.4.2
- Datasets: platypus.gl, possums.gl, testset.gl, testset.gs, plus four constructed fixtures
  (exchangeable split, all-monomorphic, all-NA loci, single-individual populations)
- Baseline: `tests/testthat/test-gl.fst.pop.R` (new file, 115 assertions, all passing against the
  reviewed state)

## Verdict

**Standards: Needs work** — the roxygen block is missing `@title`, `@description`, `@name`,
`@family` and `@details`, so `?gl.fst.pop` renders a four-sentence title; there is no FS5
validation block and no `accept = "SNP"` on the datatype check.

**Spec: Rework** — the point estimate is exactly right (it reproduces
`hierfstat::pairwise.WCfst` to 1.4e-17), but the inference layer around it is not usable as
documented: bootstrap confidence intervals and p-values are not reproducible under `set.seed()`,
and the p-value is a one-tailed bootstrap fraction, not the "probability of being different from
zero" the description promises.

What works: the delegation feeds StAMPP the right representation, the pairwise matrix dimnames
track their content under non-alphabetical population levels, the input object comes back
untouched, `verbose = 0` is fully silent, and the FBM path returns the dense answer exactly.

## Independent verification

`gl.fst.pop` wraps `StAMPP::stamppFst`, which implements the Weir & Cockerham (1984) theta.

| Check | Result |
|---|---|
| `gl.fst.pop` vs `StAMPP::stamppFst` called directly (platypus.gl) | max abs diff 0, dimnames identical |
| `gl.fst.pop` vs `hierfstat::pairwise.WCfst` (platypus.gl, `gl.filter.callrate(threshold = 1)` + `gl.filter.monomorphs`) | max abs diff 1.4e-17 over 3 pairs, r = 1.000000 |
| `gl.fst.pop` vs `hierfstat::pairwise.WCfst` (possums.gl, monomorphs filtered) | max abs diff 0.000000 over 45 pairs, r = 1.000000 |
| `gl.fst.pop` vs `hierfstat::pairwise.WCfst` (platypus.gl unfiltered) | max abs diff 2.96e-4, r = 0.99989 |

The 2.96e-4 residual on unfiltered platypus.gl is a fixture artefact, not a disagreement between
the estimators: `gl2gi` drops "markers with no scored alleles" on the hierfstat side while StAMPP
keeps them. Once monomorphic and all-missing loci are removed the two implementations agree to
machine precision on every pair of both datasets.

**Right engine, right input.** `stamppFst` requires `class(geno) == "genlight"` to take its
genlight branch; on a `dartR` object that test is FALSE and the object would fall through to the
data-frame path and error. Line 55's `class(x) <- "genlight"` is the reason it does not, and the
inline comment says so. The conversion StAMPP then applies — `as.matrix(x) * (1/ploidy)`,
`NA -> NaN` — receives 0/1/2 dosages at ploidy 2 and produces the correct allele frequencies.
Verified end to end by the exact hierfstat agreement above.

## The three inherited leads

| Lead (from the gl.report.fstat review, PR #384) | Verdict | Evidence |
|---|---|---|
| (i) `@return` says `dist`, a matrix is returned | **Confirmed** | `class(out)` is `c("matrix", "array")`; `length(as.numeric(out))` is `nPop^2`, not `nPop*(nPop-1)/2`. See F3. |
| (ii) SilicoDArT accepted | **Confirmed, and worse than "accepted"** | No error, no warning, and the returned number is driven by the ploidy slot. See F2. |
| (iii) bootstraps not seed-reproducible | **Confirmed, and it reaches the p-values** | Six calls under the same `set.seed(99)` returned six different CI limits and p-values of 0.02, 0.01, 0.01, 0.03, 0.02, 0.02. See F1. |

Lead (iii) was recorded in the sibling review as "`Fsts` reproduce under `set.seed`, `Bootstraps`
do not". That is understated. The p-values do not reproduce either; the earlier observation was
made on a saturated fixture where every p-value was 0.

## Findings

**F1 [HIGH, confidence: high] — bootstrap inference is not reproducible under `set.seed()`
(no rule fits; see the skill-maintainer note)**

`R/gl.fst.pop.r:57-61` — when `nboots > 1`, `StAMPP::stamppFst` runs the locus bootstrap inside
`foreach(...) %dopar%` on a PSOCK cluster it creates itself with `makeCluster(nclusters)` and
registers with `registerDoParallel(cl)`. The workers' RNG streams are never seeded from the master
session, and the cluster handle is internal to StAMPP, so `set.seed()` in the user's session has
no effect on the bootstrap.

Evidence: 25 loci, two exchangeable populations of 15 possums, `nboots = 100`, `set.seed(99)`
called immediately before each of six identical calls:

| Run | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| p-value | 0.02 | 0.01 | 0.01 | 0.03 | 0.02 | 0.02 |
| Lower CI | 0.00026 | 0.00305 | 0.00350 | -0.00019 | 0.00164 | 0.00088 |

The point estimates are unaffected — `Fsts` is computed before the bootstrap branch and is
identical across runs. There is no `seed` parameter in the signature.

Failure scenario: a p-value quoted in a manuscript cannot be reproduced from the same script and
the same data. On this fixture the 95 per cent CI covers zero in one run out of six while the
p-value stays below 0.05 in all six, so the reported conclusion changes with the run.

Proposed change: two options, the second subsuming the first.
(a) Document the limitation in `@details` and in the `nboots` parameter text — cheap, no
behaviour change, and honest about what the numbers are.
(b) Move the locus bootstrap into `gl.fst.pop`. StAMPP's bootstrap resamples column indices from
the precomputed per-population `p`, `oh` and `ninds` matrices; reproducing it locally is about 30
lines, removes the parallel round-trip for small `nboots`, and makes a `seed` argument possible.
This is also the only route to fixing F5 and the CI indexing in F6.

**F2 [HIGH, confidence: high] — SilicoDArT is accepted and returns a ploidy-driven artefact
(DAT7)**

`R/gl.fst.pop.r:48` — `utils.check.datatype(x, verbose = verbose)` is called without `accept`,
which under DAT7 admits SilicoDArT. The algorithm is dosage-based and has no presence/absence
meaning: `stamppFst` computes `geno * (1/ploidy)`, so 0/1 scores at ploidy 1 become allele
frequencies of exactly 0 or 1; the observed-heterozygosity term is identically zero (verified: sum
of the term over four populations of `testset.gs` is 0); and per-population sample size is
accumulated as `ploidy/2` per individual, i.e. exactly half the real n.

Evidence: on four populations of `testset.gs` the function returns 0.3032 to 0.6394 with no error
and no warning at any verbosity. The same 0/1 matrix with `ploidy` forced to 2 returns 0.0636 to
0.1391 — a five- to six-fold change driven entirely by a slot that carries no allelic information
for presence/absence data.

Failure scenario: a user runs `gl.fst.pop` on a SilicoDArT object, gets values in the range Fst
normally occupies, and reports them as Fst.

Proposed change: `datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)`.

**F3 [MEDIUM, confidence: high] — `@return` states class `dist`; a matrix is returned
(DOC5, proposed rule)**

`R/gl.fst.pop.r:18-22` — "A matrix of distances between populations (class dist)". The returned
object at `nboots = 1` is `c("matrix", "array")`: `nPop x nPop`, lower triangle populated, upper
triangle and diagonal NA. `length()` is `nPop^2`; a `dist` of the same data has length
`nPop*(nPop-1)/2`. Any code written to the documented class — `hclust()`, `ape::nj()`,
`labels()` — fails or misreads it. The `@return` hint `as.matrix(as.dist(fsts))` does work and
returns the symmetric matrix as promised.

The single live caller is unaffected: `dartR.popgen::gl.check.panel:60,62` uses
`as.numeric()` then `complete.cases()`, which yields the same 45 paired values under either class
(verified on possums.gl — 100 flattened values with 45 non-NA versus 45 dist values, same value
set).

Proposed change: correct `@return` to describe a lower-triangular base matrix. Do not change the
return class; that would be an API1 change for no benefit to the one caller.

**F4 [MEDIUM, confidence: high] — the estimator is never named (DOC5, proposed rule)**

`R/gl.fst.pop.r:1-6` — the documentation says only "based on the implementation in the StAMPP
package (?stamppFst)". It never states that the statistic is Weir & Cockerham's (1984) theta, and
there is no `@references`. The sibling `gl.report.fstat` returns Nei (1987) Gst-family statistics
under the same word "Fst"; the two diverge by up to 0.30 on `testset.gl` (per the PR #384 review).

Failure scenario: a user compares `gl.report.fstat$Fst` against `gl.fst.pop` on unevenly sampled
data, sees numbers that differ by a factor of two, and has nothing in either help page that
explains why.

Proposed change: name the estimator in `@details`, add `@references` for Weir & Cockerham (1984),
and add `@seealso gl.report.fstat` with one sentence on which estimator each returns.

**F5 [MEDIUM, confidence: high] — the p-value is a one-tailed bootstrap fraction, not the
documented quantity (DOC5, proposed rule)**

`R/gl.fst.pop.r:3-6` — "run bootstrap to estimate probability of Fst values to be different from
zero". The implemented quantity is
`pval = (number of bootstrap replicates with Fst <= 0) / nboots` — a one-tailed bootstrap achieved
significance level for H0: theta <= 0. Verified exactly: on an exchangeable split of possums
population A (`nboots = 200`), `Pvalues[2,1]` equals `mean(bootstraps <= 0)` to full precision.

Three consequences follow from the definition, none of them documented:
- Undifferentiated populations give p = 1, not a large-but-not-1 value. On the exchangeable
  fixture theta is -0.0104 and all 200 replicates are at or below zero.
- p can be exactly 0 (it is 0 for every pair of platypus.gl at `nboots = 20`). A bootstrap ASL
  should not report 0; `(k + 1)/(nboots + 1)` is the standard floor.
- No multiple-testing correction is applied or mentioned across the `nPop*(nPop-1)/2` pairs — 435
  simultaneous tests on `testset.gl`.

Failure scenario: p = 0 is quoted as "p < 0.001" when the resolution of a 20-replicate bootstrap
is 0.05, and p = 1 for an undifferentiated pair reads as positive evidence of no differentiation.

Proposed change: state the definition, the one-tailed direction, the absence of any correction and
the `1/nboots` resolution in `@details`; apply the `(k + 1)/(nboots + 1)` floor if F1(b) is taken.

**F6 [MEDIUM, confidence: high] — no parameter validation; out-of-range arguments produce
R-internals errors and silently degenerate confidence intervals (FS5)**

`R/gl.fst.pop.r:47-56` — there is no FUNCTION SPECIFIC ERROR CHECKING section between the datatype
check and the job. Observed:

| Input | Result |
|---|---|
| one population | `Error: subscript out of bounds` |
| `percent = 150` | `Error: subscript out of bounds` |
| `percent = "abc"` | `Error: non-numeric argument to binary operator` |
| `nboots = -5` | `Error: invalid 'length' argument` |
| `nboots = 0` | returns the matrix silently; `@return` says the matrix arrives "if nboots = 1" |
| `nboots = 2` | lower and upper CI limits are the same number |

The CI case is the one that returns a wrong answer rather than an unhelpful message. StAMPP indexes
the percentile as `ceiling(0.025 * nboots)`, which is 1 for every `nboots` up to 40 — so for any
`nboots <= 40` the reported "lower bound CI limit" is the sample minimum of the bootstrap
distribution, not its 2.5th percentile. At `nboots = 2` both bounds resolve to the same order
statistic.

Failure scenario: a user runs `nboots = 20` for speed and reads the resulting interval as a 95 per
cent CI; it is the observed range of 20 replicates.

Proposed change: add the FS5 block — require `nPop(x) >= 2` with a message naming the object,
`nboots` a whole number >= 1, `0 < percent < 100`, `nclusters >= 1` — and warn at
`verbose >= 1` (VRB4) when `nboots < 40`, because the requested percentile cannot be reached.

**F7 [MEDIUM, confidence: high] — undifferentiable pairs return NaN with no warning (VRB4)**

`R/gl.fst.pop.r:57-61` — populations of n = 1 make `n.bar - 1` zero and every locus term
non-finite; all-monomorphic data makes every term NaN. Both cases return a NaN cell and say
nothing at any verbosity.

Evidence: `testset.gl` ships two populations of n = 1 (`EmmacNormLeic`, `EmmacNormSalt`); their
pair is the one NaN in the 435-cell lower triangle, and `verbose = 2` output on that call is three
lines, none mentioning it. An all-monomorphic 30 x 20 fixture returns NaN for all three pairs
silently.

Failure scenario: a NaN cell propagates into a downstream mean or heatmap and is read as missing
data rather than as a population that cannot support the estimator.

Proposed change: after the call, count non-finite cells in the lower triangle and report them at
`verbose >= 1`, naming the population pairs and the likely cause (n = 1, or no polymorphic loci
shared).

**F8 [LOW, confidence: high] — the roxygen block has no `@title`, `@description`, `@name`,
`@family` or `@details` (DOC1)**

`R/gl.fst.pop.r:1-6` — the block opens with a bare title line, a blank line, then four more bare
lines. roxygen2 folds all of it into `\title`, so `man/gl.fst.pop.Rd` carries a four-sentence
title and repeats the identical text in `\description`. `?gl.fst.pop` shows that title, and so
does the package index entry.

Proposed change: split into `@title` (one line, no full stop), `@description` (one paragraph) and
`@details`; add `@name gl.fst.pop` and an `@family` matching the sibling F-statistic functions.

**F9 [LOW, confidence: high] — `@author` lacks the Author(s)/Custodian structure
(DOC7, proposed rule)**

`R/gl.fst.pop.r:25-26` — "Bernd Gruber (bugs? Post to \url{...})". Neither label is present.

Proposed change: `Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to
\url{https://groups.google.com/d/forum/dartr}`.

**F10 [LOW, confidence: high] — the `verbose` parameter text is not the canonical wording
(DOC2)**

`R/gl.fst.pop.r:15-17` — reads "[default 2, unless specified using gl.set.verbosity]" where the
canon is "[default NULL, adopting the global verbosity set by gl.set.verbosity(), or 2 if no
global is set]", and carries a stray space in "progress log ;".

Proposed change: replace with the DOC2 text verbatim.

**F11 [LOW, confidence: high] — verbosity levels 3 and 5 promise output that never appears
(VRB1, DOC5 proposed rule)**

The documented scale promises "3, progress and results summary; 5, full report". Verified by
`capture.output` on possums.gl: `verbose = 2` and `verbose = 3` produce the same three lines;
`verbose = 5` adds only the `[dartR.base vers. 1.2.3 Build = Jody ]` banner from
`utils.flag.start`. No results summary exists at any level.

Proposed change: print a summary at `verbose >= 3` — number of populations, number of pairs, range
of theta, count of non-finite pairs (this is also F7's warning) — or drop the promise from the
parameter text.

**F12 [INFO, confidence: high] — `utils.flag.start(build = "Jody")` (FS3)**

`R/gl.fst.pop.r:43-45` — the `build` argument is still accepted by `utils.flag.start` but is
outdated under FS3, and it is what surfaces the "Build = Jody" banner at `verbose = 5`.

Proposed change: drop `build = "Jody"`.

## Proposed changes

1. Add `accept = "SNP"` to the `utils.check.datatype` call (F2).
   **Consequence: SilicoDArT objects that today return a number will error.**
2. Correct `@return` to describe a lower-triangular base matrix, and note that `nboots > 1`
   (not `nboots != 1`) selects the list form (F3, part of F6). Docs only; return class unchanged.
3. Name the estimator: add `@details` stating Weir & Cockerham (1984) theta, add `@references`, add
   `@seealso gl.report.fstat` distinguishing the two (F4).
4. Document the p-value precisely in `@details` — one-tailed fraction of bootstrap replicates at or
   below zero, resolution `1/nboots`, no multiple-testing correction (F5).
5. Document the bootstrap's non-reproducibility in `@details` and in the `nboots` parameter text
   (F1a). Docs only.
6. Reimplement the locus bootstrap inside `gl.fst.pop`, add a `seed` parameter, fix the percentile
   indexing and apply the `(k+1)/(nboots+1)` floor to the p-value (F1b, F5, part of F6).
   **Consequence: p-values and confidence intervals change numerically for every call with
   `nboots > 1`.** Supersedes change 5 if taken.
7. Add the FS5 validation block: `nPop >= 2`, whole-number `nboots >= 1`, `0 < percent < 100`,
   `nclusters >= 1`, plus a `verbose >= 1` warning when `nboots < 40` (F6).
   **Consequence: single-population input and out-of-range arguments now error with a dartR
   message instead of an R-internals message.**
8. Warn at `verbose >= 1` when any pair returns a non-finite theta, naming the pairs (F7).
9. Add a results summary at `verbose >= 3` (F11).
10. Restructure the roxygen header: `@title`, `@description`, `@name`, `@family`, `@details`,
    canonical `@param verbose` text, `Author(s):`/`Custodian:` in `@author` (F8, F9, F10). Run
    `devtools::document()` in the same change (DOC4).
11. Drop `build = "Jody"` from `utils.flag.start` (F12).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: behaviour vs roxygen on platypus.gl, possums.gl, testset.gl, testset.gs — run
- Independent numerical verification: `StAMPP::stamppFst` direct and `hierfstat::pairwise.WCfst`
  on platypus.gl and possums.gl — run
- Label integrity (non-alphabetical population levels, possums.gl relevelled to reverse order,
  45 pairs recomputed pairwise) — run, 0 mismatches
- Bootstrap determinism, six same-seed runs — run
- Inference semantics on an exchangeable-populations fixture — run
- NA policy: entirely-NA loci, all-NA pop x locus cells, all-monomorphic data, single locus,
  single-individual populations — run
- Edge cases: 1, 2 and 30 populations; identical genotypes; `nboots` 0/1/2/20/50/100/200 — run
- Contract: input untouched, no history appended, `verbose = 0` silence, `nclusters = 2` — run
- FBM path (DAT6): run — `gl.gen2fbm` fixture on three possums populations matches the dense
  answer exactly. The `.fbm_or_null` branch at line 53 is exercised.
- Caller impact (API3): grep across nine live clones — one caller,
  `dartR.popgen/R/gl.check.panel.r:60,62`, both `nboots = 1`. Its consumption pattern is pinned in
  the baseline test.
- PLT rules: SKIPPED — the function does not plot
- DEP1 guard: SKIPPED as not applicable — StAMPP is in `Imports`, not `Suggests`, so no
  `requireNamespace` guard is required
- dartR Google Group / GitHub issue search: SKIPPED — no network access in this session

## Notes (outside the nominated function; recorded, not fixed)

- `dartR.popgen::gl.check.panel:60,62` pairs `as.numeric(gl.fst.pop(xorig))` with
  `as.numeric(gl.fst.pop(x))` and relies on both matrices having identical dimnames and identical
  NA patterns. That holds only because the function sorts both objects by population first
  (`x[order(pop(x)),]`) and checks the populations match. If a panel ever drops a population
  entirely, `nPop` differs, the flattened vectors are different lengths, and `complete.cases`
  would pair unrelated cells. Not a `gl.fst.pop` defect; recorded for the popgen review.
- `dartR.spatial::gl.ibd` calls `StAMPP::stamppFst` directly rather than through `gl.fst.pop`, so
  it carries the same non-reproducible bootstrap and would not pick up any fix made here.
- Recorded follow-up (custodian, 2026-09-08): both `gl.fst.pop` and `gl.report.fstat` move to
  dartR.popgen after the `utils.basic.stats` home question is settled. "Fix now, move later" — the
  move is not proposed as a finding in this review.

## Note for the skill maintainer

Two rule gaps, both hit again here.

1. **Method correctness and method naming.** F4 and F5 are statistical-method defects — the
   estimator is not named, and the inferential quantity is described as something it is not. The
   only rule that reaches them is DOC5, which is `[proposed]` and therefore cannot carry more than
   a MEDIUM. This is the third review in the distance/ordination cluster to land in the same place
   (`utils.dist.ind.snp`, `gl.report.fstat`, now `gl.fst.pop`). A `[confirmed]` rule in a new
   NUM/method-correctness group would carry these directly.
2. **Reproducibility of stochastic output.** F1 has no rule at all. A function whose documented
   output is a p-value must produce the same p-value from the same seed, and nothing in FS, DOC,
   VRB, DAT, DEP, PLT, STY, API or TST says so. Suggested: a REP group covering seed handling,
   parallel RNG streams, and the requirement that any stochastic result be reproducible from a
   documented seed argument.

## Approval

Approved 2026-09-08 by Arthur Georges via the formal approval boxes, with the consequences
acknowledged as stated against each change.

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges, 2026-09-08 | F2 HIGH. `accept = "SNP"`. Consequence acknowledged: SilicoDArT objects that return a number today will error. The one live caller, `gl.check.panel`, is SNP-context |
| 2 | Approved | Arthur Georges, 2026-09-08 | F3 MEDIUM. `@return` corrected to the matrix actually returned; the return class is unchanged and the correction is verified safe for `gl.check.panel` |
| 3 | Approved | Arthur Georges, 2026-09-08 | F4 MEDIUM. The estimator is named as Weir & Cockerham (1984) theta |
| 4 | Approved | Arthur Georges, 2026-09-08 | F5 MEDIUM. P-value semantics documented as the one-tailed `mean(replicates <= 0)`, with the absence of any multiple-testing correction across up to 435 pairs stated plainly. The definition itself is unchanged |
| 5 | Superseded | — | F1a, documentation of the non-reproducibility. Change 6 removes the defect instead |
| 6 | Approved | Arthur Georges, 2026-09-08 | F1 HIGH. Own bootstrap over loci under the caller's seed. Consequence acknowledged: p-values and confidence intervals change numerically; they were never reproducible anyway. Requirement set with the approval: two runs under the same `set.seed()` must give identical `Bootstraps` and `Pvalues`, different seeds must differ, and the point estimate must stay byte-identical |
| 7 | Approved | Arthur Georges, 2026-09-08 | F6 MEDIUM. Parameter validation, plus a guard warning for the `nboots <= 40` degenerate-CI regime. Consequence acknowledged: single-population input and out-of-range arguments now error with a dartR message |
| 8 | Approved | Arthur Georges, 2026-09-08 | F7 MEDIUM. Non-finite pairs named at `verbose >= 1` |
| 9 | Approved | Arthur Georges, 2026-09-08 | F11 LOW. Results summary at `verbose >= 3` |
| 10 | Approved | Arthur Georges, 2026-09-08 | F8, F9, F10 LOW. Roxygen structure so `?gl.fst.pop` renders, DOC7 author line, DOC2 verbose wording |
| 11 | Approved | Arthur Georges, 2026-09-08 | F12. Build tag dropped, named in the approval line alongside the LOWs |

**Relocation deferred.** The custodian chose "fix now, move later". Moving `gl.fst.pop` and
`gl.report.fstat` to dartR.popgen, and settling where `utils.basic.stats` lives, is a separate
follow-up job. Nothing in this change touches the location of either function.

## Outcome

Applied 2026-09-08 on branch `review-gl.fst.pop`, cut from `upstream/dev` at `ddaed27`, by
Claude Fable 5 (claude-fable-5, Claude Code). All 10 approved changes are in; change 5 is
superseded by change 6.

**Applied**

1. (F2) `utils.check.datatype(x, accept = "SNP", verbose = verbose)`.
2. (F3) `@return` describes the lower-triangular base matrix of class `c("matrix", "array")`,
   states that it is not a `dist`, notes that dimnames follow order of appearance rather than
   factor level, and documents the three-element list returned when `nboots > 1`.
3. (F4) `@details` names Weir and Cockerham's (1984) theta, `@references` carries the citation,
   and `@seealso` distinguishes it from the Nei (1987) Gst-family statistics of
   `gl.report.fstat`.
4. (F5) `@details` states the p-value definition (`mean(replicates <= 0)`, one-tailed), its
   `1/nboots` resolution, that p = 1 is a failure to reject rather than evidence of no
   differentiation, and that no multiple-testing correction is applied across the
   `nPop*(nPop-1)/2` pairs — 435 for a 30-population object.
6. (F1) The locus bootstrap is this function's own. `StAMPP::stamppFst` is still called, with
   `nboots = 1`, so the point estimate is byte-identical by construction. Per-population,
   per-locus allele frequency, observed heterozygosity and sample size are then computed from
   the DECODED genotype matrix, exactly as StAMPP computes them, and each replicate resamples
   the COLUMNS of those three matrices with replacement, nLoc at a time, under the calling
   session's RNG. Theta is recomputed per replicate by a local implementation of the same
   estimator. Nothing is subsetted out of the genlight with repeated indices, so adegenet's
   `SNPbin` `[` method — which drops NA when an index repeats, the SNPbin trap the sibling
   `gl.report.fstat` fix had to work around with `new("genlight", ...)` — is never reached at
   all.
7. (F6) A FUNCTION SPECIFIC ERROR CHECKING block validates `nPop(x) >= 2`, `nboots` as a whole
   number of 1 or more, `0 < percent < 100` and `nclusters >= 1`, all before any work is done,
   and warns at `verbose >= 1` when `ceiling(alpha * nboots) <= 1` — every `nboots` up to 40 at
   the default `percent = 95` — naming the number of replicates needed to reach the requested
   percentile.
8. (F7) Non-finite pairs are counted and named at `verbose >= 1`, capped at ten names with a
   count of the remainder, with the two likely causes stated.
9. (F11) A results summary prints at `verbose >= 3`: populations, pairs, theta min/max/mean,
   non-finite pair count, and for `nboots > 1` the replicate count and the no-correction
   reminder.
10. (F8, F9, F10) `@name`, `@title`, `@family distance`, `@description`, `@details` split out of
    the bare title block; `Author(s): Bernd Gruber. Custodian: Bernd Gruber` per DOC7; DOC2
    verbose wording verbatim. `devtools::document()` run.
11. (F12) `build = "Jody"` dropped from `utils.flag.start`.

**Three deliberate departures from the proposed change 6**, all following the approval text
rather than the finding text:

- **No `seed` parameter.** The approval requires reproducibility under the caller's
  `set.seed()`, which the session RNG delivers. Adding an argument is an API change that was
  not among the acknowledged consequences, and the sibling `gl.report.fstat` (PR #384) took the
  same decision. The `set.seed()` guidance is documented in `@details` instead.
- **No `(k + 1)/(nboots + 1)` floor on the p-value.** The approval names the documented quantity
  as `mean(replicates <= 0)`, so the definition is documented, not changed.
- **The percentile indexing is unchanged.** The approval asks for a guard/warning for the
  `nboots <= 40` regime, not a re-index. Changing the index would have moved the confidence
  limits for a second, unapproved reason.

**Verification** (R 4.4.2, dartR.data 1.2.5, StAMPP 1.6.3, hierfstat 0.5.11, `pdf(NULL)`):

- Point estimates byte-identical. The pinned values all hold: platypus.gl lower triangle
  0.0823972148 / 0.0745767976 / 0.0600481604, possums.gl lower-triangle sum 13.4991617506,
  `[2,1]` 0.2757253000, `[10,1]` 0.2916881324, testset.gl lower-triangle sum 176.3633320934.
  Maximum absolute difference against `StAMPP::stamppFst` called directly is 0 on both
  platypus.gl and possums.gl, with identical dimnames. Against
  `hierfstat::pairwise.WCfst` on the filtered fixtures: 1.39e-17 on platypus.gl (3 pairs) and
  5.55e-17 on possums.gl (45 pairs). `Fsts` at `nboots = 10` is `identical()` to the
  `nboots = 1` matrix.
- F1 reproducibility. On three possums populations at `nboots = 50`: two runs under
  `set.seed(123)` give `identical()` `Bootstraps` and `identical()` `Pvalues`; `set.seed(124)`
  gives different replicates and different lower CI limits. On a marginal fixture (an
  exchangeable split of possums population A, 40 loci, `nboots = 100`) the p-value is 0.96
  under seed 201 on both runs and 0.99 under seed 202 — the p-value itself now moves with the
  seed and is stable within it.
- F1 resampling unit. Intercepting `sample()` lexically and recording every draw across 50
  replicates on a 30-individual, 200-locus fixture: 10000 indices drawn, maximum index 200,
  200 distinct values — the draw spans the full locus set, not the first nInd loci.
- F1 estimator agreement. With `sample()` replaced by the identity, every replicate equals the
  reported point estimate to a maximum absolute difference of 0, on a two-population fixture,
  a three-population fixture and a fixture carrying NA cells.
- F1 NA preservation. On the same fixture, `as.matrix(g[, c(1,1,1,2,2,3:30)])` carries 7 NA
  where the decoded matrix indexed identically carries 17 — the SNPbin duplicate-index defect,
  reproduced. The function's replicate path resamples matrix columns and never touches that
  method; all replicates on the NA fixture are finite.
- F1 known answer. An exchangeable split of possums population A gives theta -0.01037 and
  p = 1.000 at `nboots = 200` — no differentiation, no signal — and the p-value equals
  `mean(replicates <= 0)` to full precision. The differentiated pair A vs B gives theta 0.2757,
  p = 0.000 and a 95 per cent interval of [0.2298, 0.3117] that brackets the point estimate.
- F2. `gl.fst.pop(testset.gs)` stops with "found SilicoDArT expecting SNP".
- F6. `nboots = 0`, `nboots = -5`, `nboots = 2.5`, `percent = 150`, `percent = "abc"`,
  `nclusters = 0` and single-population input all stop with a dartR `Fatal Error` message.
  `nboots = 20` and `nboots = 40` warn at `verbose = 1`; `nboots = 41` does not.
- F7. `testset.gl` reports "1 of 435 population pairs returned a non-finite Fst ...
  EmmacNormLeic vs EmmacNormSalt" at `verbose = 1`, and zero lines at `verbose = 0`. An
  all-monomorphic fixture reports 3 of 3.
- Caller contract (API3). `dartR.popgen::gl.check.panel:60,62` runs its
  `as.numeric()`/`complete.cases()` pattern unchanged: identical dimnames, identical NA
  pattern, 45 paired values, and the value set equal to `as.numeric(as.dist(fo))`.
- Contract. Input object byte-identical after the call, no history entry, `verbose = 0` zero
  lines at both `nboots = 1` and `nboots = 20`, `nclusters = 1` and `2` agree, and the FBM
  fixture matches the dense answer.
- Console lines on possums.gl: 2 at `verbose 1`, 3 at `verbose 2`, 7 at `verbose 3`, 7 at
  `verbose 5`. Levels 2 and 3 were identical at 3 lines before, and only the "Build = Jody"
  banner separated 5.
- Characterization tests: 141 assertions pass, no warnings, no skips. Every changed assertion
  carries an `# [approved Fn]` comment. Seven blocks flipped, each mapping to an approved
  finding: `nboots = 0` and the two error-path blocks to F6, the SilicoDArT block to F2, the
  reproducibility block to F1, the NaN block to F7, and the verbosity block to F11 and F12.
  No unexplained diff.
- Caller grep across all eight live clones plus dartR.data: `gl.check.panel` only.

The `@family distance` retag adds a `gl.fst.pop` cross-reference to `man/gl.dist.ind.Rd`,
`man/gl.dist.pop.Rd`, `man/gl.fdsim.Rd` and `man/utils.dist.ind.snp.Rd`. Only
`man/gl.fst.pop.Rd` is committed, to keep the change to one function. `devtools::document()`
also raised unrelated pre-existing drift in NAMESPACE and 34 other `man/*.Rd` files; that was
discarded.

```json
{
  "function": "gl.fst.pop",
  "package": "dartR.base",
  "family": "analysis",
  "skill_version": "2.0.0",
  "model": "claude-opus-5",
  "commit": "ddaed27",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "none (proposed REP group)", "status": "proposed", "change": 5},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DAT7", "status": "proposed", "change": 1},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "proposed", "change": 2},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "proposed", "change": 3},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "proposed", "change": 4},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "proposed", "change": 7},
    {"id": "F7", "severity": "MEDIUM", "confidence": "high", "rule": "VRB4", "status": "proposed", "change": 8},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "proposed", "change": 10},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "proposed", "change": 10},
    {"id": "F10", "severity": "LOW", "confidence": "high", "rule": "DOC2", "status": "proposed", "change": 10},
    {"id": "F11", "severity": "LOW", "confidence": "high", "rule": "VRB1", "status": "proposed", "change": 9},
    {"id": "F12", "severity": "INFO", "confidence": "high", "rule": "FS3", "status": "proposed", "change": 11}
  ],
  "verification": {
    "engine": "StAMPP::stamppFst, Weir & Cockerham (1984) theta",
    "vs_stampp_direct": {"dataset": "platypus.gl", "max_abs_diff": 0},
    "vs_hierfstat_wcfst": [
      {"dataset": "platypus.gl (cr=1, monomorphs filtered)", "pairs": 3, "r": 1.0, "max_abs_diff": 1.4e-17},
      {"dataset": "possums.gl (monomorphs filtered)", "pairs": 45, "r": 1.0, "max_abs_diff": 0.0},
      {"dataset": "platypus.gl (unfiltered)", "pairs": 3, "r": 0.99989, "max_abs_diff": 2.96e-4}
    ],
    "label_integrity": {"fixture": "possums.gl, levels reversed", "pairs": 45, "mismatches": 0}
  },
  "leads": {
    "return_class_dist_vs_matrix": "confirmed",
    "silicodart_accepted": "confirmed",
    "bootstrap_not_seed_reproducible": "confirmed, extends to p-values"
  },
  "callers": {"gl.fst.pop": ["dartR.popgen/R/gl.check.panel.r:60", "dartR.popgen/R/gl.check.panel.r:62"]},
  "caller_impact": "none for changes 2-5 and 9-11; change 1 (accept=SNP) does not affect gl.check.panel, which passes SNP data; changes 6 and 7 alter numerical output and error paths only",
  "coverage_skipped": [
    "PLT rules: function does not plot",
    "DEP1: StAMPP is in Imports, not Suggests",
    "Google Group / GitHub issue search: no network access in this session"
  ],
  "baseline_test": "tests/testthat/test-gl.fst.pop.R",
  "approval": {
    "by": "Arthur Georges",
    "date": "2026-09-08",
    "approved_changes": [1, 2, 3, 4, 6, 7, 8, 9, 10, 11],
    "superseded_changes": [5],
    "consequences_acknowledged": [
      "SilicoDArT input now errors (change 1)",
      "p-values and confidence intervals change numerically; they were never reproducible (change 6)",
      "single-population input and out-of-range arguments now error with a dartR message (change 7)"
    ],
    "deferred": "relocation of gl.fst.pop and gl.report.fstat to dartR.popgen"
  },
  "applied": {
    "date": "2026-09-08",
    "branch": "review-gl.fst.pop",
    "base_commit": "ddaed27",
    "model": "claude-fable-5",
    "departures": [
      "no seed parameter added; reproducibility comes from the session RNG, as in PR #384",
      "no (k+1)/(nboots+1) floor on the p-value; the approval documents mean(replicates <= 0)",
      "percentile indexing unchanged; the approval asks for a warning, not a re-index"
    ],
    "verification": {
      "point_estimates": "byte-identical; max abs diff 0 vs StAMPP::stamppFst direct, 1.39e-17 vs hierfstat on platypus.gl, 5.55e-17 on possums.gl",
      "seed_reproducible": true,
      "replicate_span": "10000 draws over 50 replicates, 200 distinct of 200 loci",
      "identity_replicate_equals_point_estimate": true,
      "exchangeable_split": {"theta": -0.01037, "p": 1.0},
      "tests": {"assertions": 141, "failed": 0, "warnings": 0, "skipped": 0}
    }
  },
  "status": "pr-open",
  "pr": null
}
```
