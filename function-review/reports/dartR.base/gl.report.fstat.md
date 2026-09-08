# Review: gl.report.fstat (dartR.base)

## Provenance

- Model: Claude Opus 5 (claude-opus-5, Claude Code) via the
  dartr-function-review skill; Skill: dartr-function-review v2.0.0;
  Base: upstream/dev at ddaed27 (`git diff upstream/dev --
  R/gl.report.fstat.r` empty, and the same for `R/gl.fst.pop.r` — the
  loaded code is the reviewed code); working branch integration-local at
  ed99203.
- Family mode: report (read-only contract) with analysis-grade numerical
  verification, because the F-statistics are the substance of the output.
- Datasets: platypus.gl (81 x 1000 x 3), possums.gl (300 x 200 x 10),
  testset.gl (250 x 255 x 30), testset.gs (218 x 255 x 29), plus
  constructed fixtures: a "stacked" 40 x 300 two-population set whose
  first 40 loci are undifferentiated and whose remaining 260 are fixed
  differences; a 200-locus set with 40 loci entirely absent from one
  population; a 300-locus set with 30 per cent missing at random plus 30
  loci absent from one population; a monomorphic-only set; an
  identical-frequency set; a non-alphabetical population-level relabel of
  platypus.gl. dartR.data 1.2.5, R 4.4.2, hierfstat 0.5.11, StAMPP 1.6.3,
  boot 1.3.31, gplots 3.1.3.1, adegenet 2.1.10.
- Baseline: `tests/testthat/test-gl.report.fstat.R` (26 blocks, 79
  assertions, all pass at the reviewed state; 24 `boot.ci` warnings are
  themselves a symptom of F1).
- Checks skipped: the rendered heatmap could not be inspected —
  `dendextend` is not installed on this machine, so `gl.plot.heatmap`
  returns `-1` after printing its guard message and no graphic is
  produced. The plot-side half of VRB5 is therefore argued from the code
  (there is no `verbose == 0` gate on `plot.display`), not from a
  rendered window. Parallel execution with `ncpus > 1` was run once
  (completes) but not benchmarked. The dartR Google Group was not
  searched (no browser session).

## Verdicts

**Standards: Needs work** — the house preamble is present and in order
(verbosity, working directory, start flag, datatype check), the report
contract is honoured exactly (input object returns byte-identical, no
history entry, results identical with and without plotting), and the
pairwise matrices are labelled correctly even when population levels are
not alphabetical. The gaps are a missing `verbose == 0` plot gate, no
validation of any parameter, a `@return` that describes two of the four
return shapes, and an inlined 100-line duplicate of `utils.basic.stats`.

**Spec: Rework** — the point estimates are exactly right: all four
statistics reproduce `hierfstat::basic.stats` to the function's own
4-decimal rounding on every pair tested. The confidence intervals are
not. `boot::boot` is given a data frame of individuals-by-loci and draws
its resample indices over the rows (individuals), while the statistic
applies those indices to the columns (loci). Only the first `nInd` loci
of a pair can ever enter a bootstrap replicate — 40 of 1000 on
platypus.gl — and the replicates are computed by a stale copy of
`utils.basic.stats` rather than the current one. On a fixture where the
leading loci are unrepresentative the reported 95 per cent interval
excludes the point estimate entirely. Everything the function documents
under "Confidence Intervals" needs rebuilding.

What works well: the pairwise matrix assembly. `mat_pops[[i]]` is filled
from `t(combn(npops, 2))` into `lower.tri()` and mirrored with
`t(m)[rev(lower.tri(m))]`, which is `upper.tri(m)` in column-major order;
the two orderings agree, so cells and dimnames line up. This is not a
repeat of the `gl.dist.pop` F1 positional-dimnames defect — verified
against a non-alphabetical relabel and against direct per-pair
computation on 10 populations.

## What the function actually computes

| Property | Value |
|---|---|
| Estimator family | Nei (1987) Gst family, as re-implemented in `utils.basic.stats` (a genlight-native port of `hierfstat::basic.stats`) |
| `Fst` | Nei's Gst, biased: `overall(Dst) / overall(Ht)` |
| `Fstp` | Nei's sample-size-corrected Gst: `overall(Dstp) / overall(Htp)` |
| `Dest` | Jost's D: `overall(Dstp) / (1 - overall(Hs))` |
| `Gst_H` | Hedrick's G'st: `Fstp / Gst_max` |
| Not computed | Weir & Cockerham theta, Hudson's Fst, Fis (computed by `utils.basic.stats`, dropped here) |
| Level | Pairwise only. Every pair is rebuilt with `rbind.dartR` and passed to `utils.basic.stats`; the "overall" figure of each pair is a locus-average, not a global multi-population Fst |
| Per-locus output | none — `utils.basic.stats$perloc` is discarded |
| Confidence intervals | `boot::boot` + `boot::boot.ci`, `nboots` replicates, `CI.type` one of norm/basic/perc/bca, `conf` level. No seed argument; reproducible only by setting the global RNG before the call |

Independent verification, `platypus.gl`, all three pairs, against
`hierfstat::basic.stats` on the same pair converted through
`gl2gi()` -> `hierfstat::genind2hierfstat()`:

| Pair | Fst (dartR / hierfstat) | Fstp | Dest | Ht |
|---|---|---|---|---|
| SEVERN_ABOVE vs SEVERN_BELOW | 0.0307 / 0.0307 | 0.0599 / 0.0599 | 0.0101 / 0.0101 | 0.1410 / 0.1410 |
| SEVERN_ABOVE vs TENTERFIELD | 0.0392 / 0.0392 | 0.0754 / 0.0754 | 0.0134 / 0.0134 | 0.1471 / 0.1471 |
| SEVERN_BELOW vs TENTERFIELD | 0.0434 / 0.0434 | 0.0836 / 0.0836 | 0.0150 / 0.0150 | 0.1478 / 0.1478 |

Maximum absolute difference 0.0000 on all four statistics. Tolerance:
the function rounds to 4 decimal places inside `utils.basic.stats`, so
agreement is exact at the reported precision; no looser tolerance was
needed. Cross-checked cell-for-cell against direct
`utils.basic.stats()` calls on four pairs of `possums.gl` (10
populations) — every matrix cell equals its independent per-pair value.

Estimator identity is documented correctly. The roxygen names Nei (1987)
for Fst and Fstp, Jost (2008) for Dest and Hedrick (2005) for Gst_H, and
that is what the code computes. `hierfstat::wc()` on the same pair gives
theta = 0.0598 against `Fst` = 0.0307 and `Fstp` = 0.0599: the
documented statistic is not W&C, and the docs do not claim it is.

NA policy: inherited from `utils.basic.stats`, which is per-cell
`na.rm = TRUE` (consistent with the `gl.allele.freq` policy, not the
`gl.tree.nj` defect). Loci absent from one population of a pair are
carried by the surviving population, the harmonic mean of sample sizes
skips the zero count, and `Dstp` is set to `NaN` where fewer than two
populations carry the locus so the locus drops out of the average.
Monomorphic loci are not dropped; on monomorphic-only data every cell
returns `NaN` rather than erroring. Populations with identical
frequencies give small negative Fst (-0.011 to -0.015 on the fixture),
which is the expected behaviour of a ratio-of-averages Gst.

## Findings

**F1 [BLOCKER, confidence: high] — the bootstrap resamples loci using
individual indices (spec axis; no catalogue rule covers "documented
numerical method not implemented" — see the note to the skill maintainer
at the end of this section)**

`R/gl.report.fstat.r:455` — the statistic passed to `boot` is
`df <- x[, indices]`, i.e. `indices` selects columns (loci). The data
handed to `boot::boot` at `:502` and `:550` is
`as.data.frame(as.matrix(tpop))`, which is individuals-by-loci, and
`boot` draws its resample from `1:nrow(data)` — the individuals.

Two consequences, both verified:

1. Only the first `nInd` loci of the pair can ever be drawn. Spying on
   the indices `boot` generates for the SEVERN_ABOVE/SEVERN_BELOW pair of
   platypus.gl (40 individuals, 1000 loci): every replicate is 40 indices
   long, and the maximum index over 20 replicates is 40. Loci 41 to 1000
   — 96 per cent of the data — never enter any interval. Each replicate
   is therefore a 40-locus dataset, and the intervals are far wider than
   a locus bootstrap of the full data would give. `boot.ci` emits
   "extreme order statistics used as endpoints" repeatedly for the same
   reason.
2. If the leading loci are unrepresentative, the interval does not cover
   the estimate. On the stacked fixture (first 40 loci undifferentiated,
   remaining 260 fixed differences) the function reports
   `Value = 0.8622` with a 95 per cent percentile interval of
   `[-0.0175, -0.0175]` — degenerate, and entirely below the reported
   value. The bootstrap mean on platypus.gl is 0.0260 against a point
   estimate of 0.0307, tracking the 40-locus subset value of 0.0244.

Failure scenario: any user who reports a confidence interval from this
function is reporting the sampling variability of the first `nInd` loci
in the file, in file order. Reordering the loci changes the interval;
adding loci beyond the `nInd`-th does not narrow it.

Proposed change: bootstrap over an explicit locus index. Hand `boot` a
one-column data frame of locus positions and have the statistic subset
the genotype matrix by `d$loc[i]`, so `R` replicates each draw
`nLoc` loci with replacement.

**F2 [HIGH, confidence: high] — the bootstrap kernel is a stale copy of
`utils.basic.stats` (STY3, single source of truth)**

`R/gl.report.fstat.r:360-462` — `pop.diff_fun` duplicates the body of
`utils.basic.stats` inline so it can run on a matrix. The copy predates
the review of that helper (commit 62d5525) and is missing both
corrections it received:

- `utils.basic.stats.r:66` guards the harmonic mean of per-population
  sample sizes with `y[y > 0]`; the copy at `:390` does not, so a locus
  absent from one population of the pair drives `mn` to zero and poisons
  every cross-population statistic for that locus.
- `utils.basic.stats.r:131` sets `Dstp[n.pop < 2] <- NaN`; the copy at
  `:421` divides by `n.pop - 1 = 0` instead.

Consequences, verified: on clean platypus.gl data the two kernels already
disagree (`Fst` 0.0307 point estimate against 0.0309 in the bootstrap
kernel; `Gst_H` 0.0754 against 0.0753) for two of the three pairs, so
the reported `Value` and the distribution the interval is built from are
different estimators. On a 200-locus fixture with 40 loci absent from one
population the point estimate is `Fst` 0.0219 and the bootstrap kernel
gives 0.0256, `Gst_H` 0.0734 against 0.0651 — and asking for a bootstrap
aborts the whole call with `missing value where TRUE/FALSE needed` from
`boot.ci`, because the kernel returns `NaN`.

Failure scenario: a user with loci that failed in one population gets no
result at all; a user without them gets an interval around a slightly
different statistic than the one printed beside it.

Proposed change: extract the estimator from `utils.basic.stats` into a
matrix-level kernel that both the helper and this function call, so
there is one implementation. Failing that, call `utils.basic.stats` on
the resampled subset.

**F3 [HIGH, confidence: high] — SilicoDArT data is admitted to
diploid-dosage arithmetic (DAT7)**

`R/gl.report.fstat.r:342` — `utils.check.datatype(x, verbose = verbose)`
without `accept = "SNP"`. The default admits both datatypes. The
arithmetic is SNP-dosage-specific throughout: heterozygosity is
`sgl_mat == 1` (`:395`, and `utils.basic.stats.r:73`) and allele
frequency is `colMeans(...) / 2` (`:402`). For SilicoDArT a score of 1
means the tag is present, so every present-tag individual is counted as a
heterozygote and the "allele frequency" is half the presence rate.

Failure scenario: `gl.report.fstat(testset.gs)` returns a full set of
pairwise matrices with `Fst[2,1] = 0.06` and no warning. The numbers are
meaningless but indistinguishable from a valid SNP result.

Proposed change: `utils.check.datatype(x, accept = "SNP", verbose =
verbose)`. Note the same gap exists in `gl.fst.pop` (sibling note).

**F4 [MEDIUM, confidence: high] — `verbose = 0` is not silent and does
not suppress the plot (VRB5)**

`R/gl.report.fstat.r:675` — `if (plot.display & npops > 2)` has no
verbosity gate, and the preamble never applies the house idiom
`if (verbose == 0) plot.display <- FALSE`. With the default
`plot.display = TRUE`, `capture.output(gl.report.fstat(platypus.gl,
verbose = 0))` returns one line. On this machine that line is
`gl.plot.heatmap`'s dependency guard, because `dendextend` is absent;
where `dendextend` is installed the same code path opens a heatmap
window at `verbose = 0` instead (argued from the code, not observed —
see Coverage).

Failure scenario: a script looping over datasets at `verbose = 0` still
paints a heatmap per iteration.

Proposed change: add `if (verbose == 0) plot.display <- FALSE` to the
preamble.

**F5 [MEDIUM, confidence: high] — no parameter or precondition validation
(FS5)**

`R/gl.report.fstat.r:328-358` — the error-checking block is absent
entirely. Observed failures:

| Call | Result |
|---|---|
| one population | `Error: n < m` (from `combn(1, 2)`) |
| all populations singletons | `Fatal Error: no populations listed to keep!` (from `gl.keep.pop`) |
| `CI.type = "rubbish"` | `Error: subscript out of bounds` after the full bootstrap has run |
| `nboots = 1` | `Error: replacement has length zero` |
| `nboots = 0.5` | `Error: dims [product 0] do not match the length of object [20]` |
| `nboots = -5` | no error; returns a malformed object carrying the string "All values of t are equal to ..." |

Failure scenario: a user who passes `nboots = 1` — the default of the
sibling `gl.fst.pop` — gets an R internals message with no indication of
which argument is wrong; a user who passes a negative count gets an
object that looks like a result.

Proposed change: validate after the datatype check — `nPop(x) >= 2`
after singleton removal, `nboots` a non-negative whole number and not 1,
`CI.type %in% c("norm","basic","perc","bca")` and length 1,
`conf` in (0, 1), `plot.stat %in% c("Fst","Fstp","Dest","Gst_H")` — each
with `stop(error(...))`.

**F6 [MEDIUM, confidence: high] — `@return` describes two of the four
return shapes (DOC1, DOC5 proposed rule)**

`R/gl.report.fstat.r:306-310` states "Two lists, the first list contains
matrices ... the second list contains tables". The four branches at
`:720-739` return:

| `nboots` | populations | returned |
|---|---|---|
| `> 0` | `> 2` | `list(Stat_matrices, Confidence_Intervals)` — matches |
| `> 0` | `<= 2` | `list(Stat_tables, Confidence_Intervals)` — first element is a data frame, not matrices |
| `0` | `> 2` | `list(Stat_matrices, <unnamed data frame>)` — second element has no name |
| `0` | `<= 2` | a bare `data.frame`, not a list at all |

Failure scenario: code written against the documented contract —
`res$Stat_matrices` or `res[[2]]$something` — fails on a two-population
dataset, which is the commonest case for a pairwise comparison.

Proposed change: document all four shapes, and name the second element
`Stat_tables` in the `nboots = 0`, `> 2` populations branch.

**F7 [LOW, confidence: high] — every heatmap failure is reported as a
plot-pane size problem (FS5)**

`R/gl.report.fstat.r:676-690` — the `tryCatch` around `p3()` maps all
errors to "Your plot was not shown in full because your 'Plots' pane is
too small". `plot.stat` is never validated and is dereferenced lazily
inside `create_heatmap` at `:662`, so `plot.stat = "Nonsense"` produces
`mat_pops[["Nonsense"]]` -> `NULL` -> an error -> the plot-pane message.

Failure scenario: a typo in `plot.stat` is reported as a window-size
problem, and the user resizes the pane instead of fixing the argument.

Proposed change: validate `plot.stat` in the F5 block, and either narrow
the `tryCatch` or include `conditionMessage(e)` in the warning.

**F8 [LOW, confidence: high] — the default `CI.type = "bca"` fails at the
replicate count the help recommends (DOC5, proposed rule)**

`R/gl.report.fstat.r:318` — `CI.type = "bca"` is the default. The help
(`:199`) says "Consider increasing the number of bootstrap replicates to
at least 200". At `nboots = 30` the call aborts with `estimated
adjustment 'a' is NA`; at `nboots = 200` it completes but warns "extreme
order statistics used as endpoints" repeatedly (itself a symptom of F1).

Failure scenario: the first thing a user tries — a small `nboots` to
check the call works — is a hard error rather than a wide interval.

Proposed change: either default to `"perc"`, or emit an informative
`stop(error(...))` when `CI.type == "bca"` and `nboots < 200`.

**F9 [LOW, confidence: high] — the results summary prints at
`verbose >= 2` (VRB1)**

`R/gl.report.fstat.r:629` — the whole result object is printed under
`if (verbose >= 2)`. VRB1 reserves level 2 for a progress log and level 3
for the results summary. Measured on platypus.gl: verbose 1 gives 2
lines, verbose 2 gives 48, verbose 3 gives the identical 48, verbose 5
gives 49. Levels 2 and 3 are indistinguishable.

Failure scenario: a user asking for a progress log on a 30-population
dataset gets 435 pairwise rows and four 30 x 30 matrices dumped to the
console.

Proposed change: gate the print block on `verbose >= 3`.

**F10 [LOW, confidence: high] — removal of one-individual populations is
announced only at `verbose >= 2` (VRB4, proposed rule)**

`R/gl.report.fstat.r:345-354` — populations of one individual are
dropped, together with their individuals, before anything is computed.
The warning is inside `if (verbose >= 2)`. Verified: at verbose 0 and
verbose 1 nothing is said; at verbose 2 the message appears. On
`testset.gs` this silently reduces 29 populations to 25.

Failure scenario: a user running at verbose 1 sees a 25 x 25 matrix from
a 29-population object and no explanation.

Proposed change: print the warning at `verbose >= 1`, and name the
populations removed.

**F11 [LOW, confidence: high] — the pairwise table is named differently
in each branch (STY3)**

`R/gl.report.fstat.r:734` wraps the table as
`data.frame(Stat_tables = pairpop_res)`, which prefixes every column:
`Stat_tables.A_vs_B`. At `:738` the same wrap on a single-column frame
produces the bare pair name, `SEVERN_ABOVE_vs_SEVERN_BELOW`. The verbose
print at `:647` and `:652` reproduces both.

Failure scenario: code that selects `res[[2]][["A_vs_B"]]` works for two
populations and returns `NULL` for three or more.

Proposed change: drop the `Stat_tables =` wrap and set the list element
name instead.

**F12 [LOW, confidence: high] — roxygen defects (DOC7, DOC2, DOC1,
DOC5)**

- `:259` — `@author Custodian: Luis Mijangos` with no `Author(s):` line
  (DOC7). Proposed: add `Author(s): Luis Mijangos.`
- `:44-46` — the `verbose` text ends "[default NULL, unless specified
  using gl.set.verbosity]" rather than the DOC2 canonical clause.
- `:133-135` — Jost's D is illustrated with `\figure{Dstequation.jpg}`,
  the same figure already used for `Dst` at `:105` and for `Dstp` at
  `:118`. Three different quantities share one image; `man/figures/` has
  no `Destequation.jpg` or `Dstpequation.jpg`.
- `:5` — `@family matched reports`. The house tag is `matched report`
  (15 files use it, 3 use the plural). There is also no `gl.filter.fstat`
  for this to be matched to.
- `:306` — `@return` sits after `@export`, the historical order ruled
  outdated under DOC1.

Failure scenario: the pkgdown reference index splits the report family
across two headings, and the manual shows the wrong equation for Jost's D.

Proposed change: fix each as listed and run `devtools::document()`.

**F13 [LOW, confidence: high] — `utils.flag.start` is called with the
outdated `build` argument (FS3)**

`R/gl.report.fstat.r:337-339` — `build = "v.2023.2"`. FS3 records
`build=` as outdated. Proposed change: drop the argument.

**F14 [INFO, confidence: medium] — the genotype matrix is densified per
pair (DAT6, proposed rule)**

`R/gl.report.fstat.r:499` and `:546` build
`as.data.frame(as.matrix(tpop))` for every pair, and `utils.basic.stats`
calls `lapply(seppop(x), as.matrix)`. FBM-backed objects work — verified
with `gl.gen2fbm(platypus.gl)`, which returns results identical to the
dense path — but they are materialised in full. Unlike `gl.fst.pop`
(`gl.fst.pop.r:53`) there is no `.fbm_or_null` branch. No proposed change
beyond recording it; densifying one pair at a time is the bounded case.

*Note to the skill maintainer*: F1 has no rule ID. The catalogue has no
rule for "the numerical method the documentation describes is not the
method implemented" — DOC5 covers documentation-behaviour agreement but
is `[proposed]` and so cannot carry a BLOCKER, and no DAT/STY rule
reaches a statistical-method defect. A confirmed rule in a new
"NUM"/method-correctness group would have caught this directly.

## Overlap analysis

### Is `gl.report.fstat` the same statistic as `gl.fst.pop`?

**No — two different estimators, which coincide numerically on balanced
data.** `gl.report.fstat` computes Nei's (1987) Gst family;
`gl.fst.pop` is a five-line wrapper over `StAMPP::stamppFst`, which
computes Weir & Cockerham's (1984) theta.

Cell-for-cell, upper triangle, on three datasets:

| Dataset | pops | `Fst` (Nei, biased) vs W&C | `Fstp` (Nei, unbiased) vs W&C |
|---|---|---|---|
| platypus.gl | 3 | r = 0.9998, max abs 0.0390, mean -0.0346 | r = 1.0000, max abs 0.0012, mean +0.0006 |
| possums.gl | 10 | r = 0.9970, max abs 0.1652, mean -0.1195 | r = 1.0000, max abs 0.0000, mean +0.0000 |
| testset.gl | 30 | r = 0.9828, max abs 0.3041, mean -0.1496 | r = 0.9979, max abs 0.0724, mean -0.0039 |

Side by side on platypus.gl:

| Pair | Nei Fst | Nei Fstp | W&C (StAMPP) |
|---|---|---|---|
| SEVERN_ABOVE vs SEVERN_BELOW | 0.0307 | 0.0599 | 0.0600 |
| SEVERN_ABOVE vs TENTERFIELD | 0.0392 | 0.0754 | 0.0746 |
| SEVERN_BELOW vs TENTERFIELD | 0.0434 | 0.0836 | 0.0824 |

On possums.gl, which has ten equal-sized populations, `Fstp` and
`stamppFst` agree to four decimal places on all 45 pairs. The two
estimators separate where sample sizes are uneven — testset.gl, with 30
populations of 2 to 30 individuals, reaches a 0.0724 discrepancy. The
biased `Fst` is systematically below theta everywhere, by 0.03 to 0.15
on average.

This is not a disagreement to be resolved: `Fst` and `Fstp` are the two
Nei estimators, `stamppFst` is the W&C estimator, and each function
documents what it computes. The practical point for the custodian is
that **`gl.report.fstat$Fstp` and `gl.fst.pop` answer the same question
and, on well-sampled data, return the same number.**

On a missing-data fixture (30 per cent missing at random plus 30 loci
absent from one population) both return finite values on every pair and
stay within 0.0015 of each other:

| Pair | Nei Fst | Nei Fstp | W&C |
|---|---|---|---|
| SEVERN_ABOVE vs SEVERN_BELOW | 0.0229 | 0.0489 | 0.0484 |
| SEVERN_ABOVE vs TENTERFIELD | 0.0345 | 0.0742 | 0.0731 |
| SEVERN_BELOW vs TENTERFIELD | 0.0419 | 0.0804 | 0.0790 |

### How they differ in everything except the number

| | `gl.report.fstat` | `gl.fst.pop` |
|---|---|---|
| Engine | `utils.basic.stats` (dartR-native port of `hierfstat::basic.stats`) | `StAMPP::stamppFst` |
| Statistics | Fst, Fstp, Dest (Jost), Gst_H (Hedrick) | theta only |
| Output shape | list of four `nPop x nPop` matrices plus a 4 x nPairs table; two other shapes for 2 populations or `nboots = 0` (F6) | lower-triangular matrix at `nboots = 1`; `list(Fsts, Pvalues, Bootstraps)` above that. The `@return` says class `dist`; the observed class is `matrix` (sibling note) |
| Inference | bootstrap over loci with `boot`, four CI types — broken (F1, F2) | bootstrap over loci inside StAMPP, percentile CI plus a p-value per pair. Working: `Fsts` reproduce under `set.seed`, `Bootstraps` do not (sibling note) |
| NA policy | `na.rm = TRUE` per cell, `n.pop < 2` loci excluded from the average | delegated to StAMPP |
| Datatype gate | none (F3) — `testset.gs` returns 0.06 | none — `testset.gs` returns 0.6484 (sibling note) |
| Singleton populations | dropped silently (F10) | not handled |
| Speed | platypus.gl 0.14 s, possums.gl 0.85 s, testset.gl 5.64 s | 1.30 s, 0.82 s, 2.71 s |
| Plot | heatmap of one chosen statistic | none |
| Callers | none | 1 |

`gl.report.fstat` is faster on few populations and many loci, and slower
as populations multiply — it rebuilds a two-population genlight with
`rbind.dartR` and re-runs the whole `utils.basic.stats` pipeline for each
of the `nPop*(nPop-1)/2` pairs.

### Other F-statistic paths in the family

| Path | What it computes | Relation |
|---|---|---|
| `dartR.base::utils.basic.stats` | the Nei engine itself; Ho, Hs, Fis, Fst, Fstp, Dest, Gst_H per locus and overall | the substance of `gl.report.fstat`; also called by `gl.diagnostics.hwe` (Fst/Fis scatter and jackknife standard errors) and `dartR.popgen::gl.ld.haplotype` |
| `dartR.popgen::utils.outflank.fst` / `utils.outflank.diploids` | W&C Fst and FstNoCorr per locus, for OutFLANK outlier detection | a third, independent W&C implementation, per-locus not pairwise |
| `dartR.spatial::gl.ibd` | calls `StAMPP::stamppFst` directly for its `distance = "Fst"` option | bypasses `gl.fst.pop` entirely — the same computation, wired twice |
| `dartR.sim::gl.diagnostics.sim` | `hierfstat::pairwise.neifst` | a fourth path, and the only one that calls hierfstat rather than porting it |
| `dartR.base::gl.dist.pop` | methods `euclidean`, `reynolds`, `nei`, `chord`, `fixed-diff`, `simple` | no Fst method. `reynolds` is a coancestry distance related to Fst but is not Fst; `fixed-diff` counts fixed differences; `nei` is Nei's standard genetic distance D |
| `dartR.captive::gl.report.kin.sets` | `Fst[s,t] = 1 - GDt/GDb` from kinship | a different quantity (PMx's kinship-based analogue). Note only — no overlap |

Five separate implementations of pairwise or per-locus Fst live in the
family: one Nei port, two W&C wrappers over StAMPP, one W&C
re-implementation inside OutFLANK, and one hierfstat call. That is a
consolidation question for the team, beyond this review.

### Caller inventory (eight live clones; dated backup copies excluded)

Clones searched: `dartR.base`, `dartR.captive`, `dartR.data`,
`dartR.popgen`, `dartR.sexlinked`, `dartR.sim`, `dartR.spatial`,
`dartRstartup`, `dartRverse`.

- **`gl.report.fstat`: no callers.** The only occurrences outside its own
  definition are its own `@examples` (`:264`) and `@details` (`:232`) and
  the generated `man/gl.report.fstat.Rd`. No sibling package, vignette or
  test calls it.
- **`gl.fst.pop`: one caller.** `dartR.popgen/R/gl.check.panel.r:60` and
  `:62`, both `nboots = 1`, inside the `parameter == "Fst"` branch.
- **`utils.basic.stats`: three callers.**
  `dartR.base/R/gl.diagnostics.hwe.r:203` and `:307` (the latter through
  `utils.jackknife`), `dartR.popgen/R/gl.ld.haplotype.r:436`, and
  `gl.report.fstat` itself (as the inlined copy, F2).

### Recommendation

**Keep both, fix `gl.report.fstat` in place, then move both to
dartR.popgen.**

*Why not merge.* They deliver different things. Merging into `gl.fst.pop`
loses the Jost and Hedrick statistics, which have no other home in the
family; merging into `gl.report.fstat` loses the per-pair p-value, which
`gl.check.panel` depends on. A single function taking an `estimator`
argument would have to reconcile two output shapes and two bootstrap
schemes for no gain.

*Why not deprecate `gl.report.fstat`.* Zero callers and a broken
bootstrap make it the obvious candidate, but the Nei/Jost/Hedrick panel
is the only implementation of Dest and Gst_H in the dartRverse, and the
point estimates are exactly right. The defect is confined to the
inference layer.

*Why not deprecate `gl.fst.pop`.* It is the only working inferential Fst
in the family today and it has a live caller.

*On relocation.* Both belong in dartR.popgen. Neither has anything to do
with dartR.base's remit of reading, filtering and reshaping genlight
objects; both are population-genetic estimators, which is exactly the
line the custodian has been drawing (`gl.amova` to popgen, `gl.Ho`/
`gl.He` to sim, `gl.propShared` to spatial). The cost:

- `gl.report.fstat` — 0 callers, so nothing breaks. The move carries the
  function, `man/gl.report.fstat.Rd`, the nine `man/figures/*.jpg`
  equation images, and the `boot` and `gplots` dependencies. It also
  needs `gl.plot.heatmap` and `gl.colors`, both exported from
  dartR.base, so popgen would import them — it already depends on
  dartR.base.
- `gl.fst.pop` — 1 caller, and that caller (`gl.check.panel`) is already
  in dartR.popgen, so the move removes a cross-package call rather than
  adding one. `StAMPP` moves from dartR.base's dependencies to
  dartR.popgen's. `dartR.spatial::gl.ibd` calls `StAMPP::stamppFst`
  directly and is unaffected.
- The one real decision: `utils.basic.stats` is used by
  `gl.diagnostics.hwe`, which stays in dartR.base. Either that helper
  stays in base and popgen imports it, or `gl.diagnostics.hwe` moves too.
  Settle the helper's home before moving `gl.report.fstat`, because F2's
  fix (a shared kernel) should land on one side of the boundary, not
  both.

*Sequence.* Fix F1, F2 and F3 in dartR.base first, so the relocation
commit is a pure move with no behaviour change in it.

## Proposed changes

1. Bootstrap over loci: hand `boot::boot` an explicit locus index and
   subset the genotype matrix by it, so every replicate draws `nLoc`
   loci with replacement instead of the first `nInd` (F1).
   **Consequence: every confidence interval this function has ever
   produced changes, and typically narrows.**
2. Replace the inlined `pop.diff_fun` with a call to the current
   `utils.basic.stats` estimator — extracting a matrix-level kernel that
   both call (F2).
   **Consequence: bootstrap values shift slightly on all data, and calls
   that currently abort on loci absent from one population start
   returning intervals.**
3. Restrict the datatype to SNP: `utils.check.datatype(x, accept =
   "SNP", verbose = verbose)` (F3).
   **Consequence: `gl.report.fstat` on SilicoDArT data now errors where
   it previously returned numbers.**
4. Add `if (verbose == 0) plot.display <- FALSE` to the preamble (F4).
5. Add the missing error-checking block: `nPop(x) >= 2` after singleton
   removal, `nboots` a non-negative whole number other than 1, `CI.type`
   a single valid choice, `conf` in (0, 1), `plot.stat` one of the four
   statistic names (F5, F7).
   **Consequence: `nboots = 1` and negative `nboots` now error instead of
   failing obscurely or returning a malformed object.**
6. Document all four return shapes in `@return`, and name the second
   element `Stat_tables` in the `nboots = 0` branch (F6, F11).
7. Include the underlying error message in the heatmap `tryCatch`
   warning instead of always reporting a pane-size problem (F7).
8. Default `CI.type` to `"perc"`, or refuse `"bca"` below 200 replicates
   with an informative message (F8).
   **Consequence: the default interval type changes, so default output
   changes for `nboots > 0` callers.**
9. Move the results-summary print from `verbose >= 2` to `verbose >= 3`
   (F9).
10. Announce singleton-population removal at `verbose >= 1` and name the
    populations dropped (F10).
11. Roxygen: add the `Author(s):` line, adopt the DOC2 verbose wording,
    correct the Jost's D figure, change `@family matched reports` to
    `matched report`, move `@return` before `@export`, drop `build=`
    from `utils.flag.start`; run `devtools::document()` (F12, F13).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run.
- Report contract (input untouched, no history, plot-decoupled) — run,
  all three hold.
- Spec: behaviour against roxygen on platypus.gl, possums.gl,
  testset.gl, testset.gs — run.
- Independent numerical verification against `hierfstat::basic.stats`
  (all three platypus.gl pairs, all four statistics) and against direct
  `utils.basic.stats` per-pair calls (four possums.gl pairs) — run,
  exact agreement at 4 decimal places.
- Estimator identity against `hierfstat::wc` and `StAMPP::stamppFst` —
  run.
- Label integrity on a non-alphabetical population-level relabel and on
  10 populations — run, no defect.
- NA policy, monomorphic loci, identical frequencies, single population,
  singleton populations, two populations, monomorphic-only data — run.
- Bootstrap: index provenance, seed determinism, `nboots`/`conf`/
  `CI.type` honoured, CI coverage on a stacked fixture — run.
- FBM path (DAT6) — run with `gl.gen2fbm(platypus.gl)`; results identical
  to the dense path.
- Every documented parameter exercised: `nboots`, `conf`, `CI.type`,
  `ncpus`, `plot.stat`, `plot.display`, `palette.divergent` (default
  only), `font.size` (default only), `plot.dir`, `plot.file`, `verbose`,
  `...` (not exercised) — `palette.divergent`, `font.size` and `...`
  reach `gl.plot.heatmap` only, which could not run here.
- Rendered heatmap: SKIPPED — `dendextend` is not installed, so
  `gl.plot.heatmap` returns `-1` without drawing. The VRB5 plot-side
  claim in F4 is from code reading; the text-side claim is measured.
- `ncpus > 1`: run once (completes on Windows with `parallel = "snow"`);
  not benchmarked, and the parallel workers' access to the statistic was
  not inspected.
- dartR Google Group / GitHub issues: SKIPPED — no browser session.

## Approval

Approved 2026-09-08 by Arthur Georges via the formal approval boxes, with
the consequences acknowledged as stated against each change.

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges, 2026-09-08 | F1 BLOCKER. Resample loci, not individuals. Consequence acknowledged: every confidence interval the function has produced was meaningless, so all CI output changes |
| 2 | Approved | Arthur Georges, 2026-09-08 | F2 HIGH. Delete the inlined copy, call `utils.basic.stats`. Consequence acknowledged: bootstrap values shift, and the absent-loci abort disappears |
| 3 | Approved | Arthur Georges, 2026-09-08 | F3 HIGH. `accept = "SNP"`. Consequence acknowledged: SilicoDArT now errors. No callers anywhere, so no migration cost |
| 4 | Approved | Arthur Georges, 2026-09-08 | F4 MEDIUM |
| 5 | Approved | Arthur Georges, 2026-09-08 | F5, F7 MEDIUM |
| 6 | Approved | Arthur Georges, 2026-09-08 | F6, F11 MEDIUM |
| 7 | Approved | Arthur Georges, 2026-09-08 | F7 LOW |
| 8 | Approved | Arthur Georges, 2026-09-08 | F8 LOW. Applied as the second of the two options the finding offered: `CI.type` keeps its `"bca"` default and fewer than 200 replicates is refused. Changing the default would have altered numerical output for every existing CI caller, which was not among the acknowledged consequences |
| 9 | Approved | Arthur Georges, 2026-09-08 | F9 LOW |
| 10 | Approved | Arthur Georges, 2026-09-08 | F10 LOW |
| 11 | Approved | Arthur Georges, 2026-09-08 | F12, F13 LOW |

F14 (INFO, densification per pair) is a note and was not actioned.

**Relocation deferred.** The custodian chose "fix now, move later". Moving
`gl.report.fstat` and `gl.fst.pop` to dartR.popgen, and settling where
`utils.basic.stats` lives, is a separate follow-up job. Nothing in this
change touches `gl.fst.pop` or the location of `utils.basic.stats`.

## Outcome

Applied 2026-09-08 on branch `review-gl.report.fstat`, cut from
`upstream/dev` at `ddaed27`, by Claude Fable 5 (claude-fable-5, Claude
Code). All 11 approved changes are in; F14 remains a note.

**Applied**

1. (F1) `boot::boot` is given `data.frame(loc = seq_len(nLoc(tpop)))` and
   the statistic subsets the genotype matrix by `loc.index$loc[indices]`,
   so each replicate draws `nLoc` loci with replacement. The replicate
   genlight is rebuilt from the decoded matrix with `new("genlight", ...)`
   rather than by subsetting the genlight, because `adegenet`'s `SNPbin`
   `[` method drops `NA` on repeated indices and a with-replacement draw
   always produces repeated indices.
2. (F2) `pop.diff_fun`, the inlined copy of the estimator, is deleted. The
   statistic calls `utils.basic.stats()`, the function that produces the
   point estimate.
3. (F3) `utils.check.datatype(x, accept = "SNP", verbose = verbose)`.
4. (F4) `if (verbose == 0) plot.display <- FALSE` in the preamble.
5. (F5, F7) A function-specific error-checking block validates `nboots`,
   `CI.type`, `conf`, `plot.stat`, the `bca`/200 combination and the
   two-population precondition, all before any work is done.
6. (F6, F11) `@return` documents all four shapes; the `nboots = 0`,
   `> 2` populations branch names its second element `Stat_tables`; the
   `data.frame(Stat_tables = ...)` wrap is gone from both branches, so
   columns are the bare pair names everywhere.
7. (F7) The heatmap `tryCatch` reports `conditionMessage(e)`.
8. (F8) `CI.type = "bca"` with fewer than 200 replicates stops with an
   informative message. The default is unchanged.
9. (F9) The results summary prints at `verbose >= 3`.
10. (F10) Singleton-population removal is announced from `verbose >= 1`
    and names the populations dropped.
11. (F12, F13) `Author(s):` line added; DOC2 `verbose` wording; the two
    quantities illustrated with the wrong figure (`Dstp` and `Dest`, both
    carrying `Dstequation.jpg`) now carry their own `\deqn{}` formulae,
    since `man/figures/` has no image for either; `@family matched
    report`; `@return` before `@export`; `build =` dropped from
    `utils.flag.start`. `devtools::document()` run.

**No seed argument was added** — the report did not propose one. The
`set.seed` guidance is documented in the Confidence Intervals section of
the roxygen instead.

**Verification** (R 4.4.2, dartR.data 1.2.5, hierfstat 0.5.11,
`pdf(NULL)`):

- Point estimates byte-identical before and after the change:
  `Fst` 0.0309 / 0.0392 / 0.0436, `Fstp` 0.0599 / 0.0754 / 0.0836,
  `Dest` 0.0101 / 0.0134 / 0.0150, `Gst_H` 0.0753 / 0.0958 / 0.1062 on
  `platypus.gl`, and exactly equal to `utils.basic.stats()` called
  directly on each of the three pairs.
- These are not the values pinned by the Phase A baseline
  (`Fst` 0.0307 / 0.0392 / 0.0434, `Gst_H` 0.0754 / 0.0958 / 0.1063),
  because the branch is cut from `upstream/dev`, which does not yet carry
  the reviewed `utils.basic.stats` (PR #309, open). Copying that helper
  into the tree restores the pinned values exactly and brings agreement
  with `hierfstat::basic.stats` to a maximum absolute difference of
  0.0000 over three pairs and three statistics. On `dev` as it stands,
  `Fstp` and `Dest` already agree with hierfstat exactly and `Fst` sits
  2e-4 away on two of the three pairs. The difference is in
  `utils.basic.stats`, not in `gl.report.fstat`.
- (F1) Bootstrap indices now span the full locus set: over 21 calls on the
  SEVERN_ABOVE/SEVERN_BELOW pair (40 individuals, 1000 loci) each
  replicate is 1000 indices long and all 1000 loci appear, against a
  maximum index of 40 before. The stacked fixture that produced
  `Value = 0.8622` with a 95 per cent percentile interval of
  `[-0.0175, -0.0175]` now returns `[0.8211, 0.9068]`; on all four
  statistics the interval brackets the value.
- (F2) A replicate built from the identity index reproduces the point
  estimate exactly (`identical()` on all four statistics). The 200-locus
  fixture with 40 loci absent from one population, which aborted with
  "missing value where TRUE/FALSE needed", now completes and returns
  `Fst` 0.0256 with `[0.0117, 0.0461]`.
- (F3) `gl.report.fstat(testset.gs)` stops with "Fatal Error:
  inappropriate object passed to function, found SilicoDArT expecting
  SNP".
- (F4) `capture.output(gl.report.fstat(platypus.gl, verbose = 0,
  plot.display = TRUE))` is zero lines. Before the change it was one line:
  `gl.plot.heatmap`'s `dendextend` guard. Its absence is direct evidence
  that the plot path is no longer entered at `verbose = 0`, which closes
  the plot-side check the review had to argue from the code.
- (F9) Console lines on `platypus.gl`: verbose 1 gives 2, verbose 2 gives
  5, verbose 3 gives 43. Levels 2 and 3 were identical at 48 before.
- (F5, F8) `nboots = 1`, `nboots = -5`, `nboots = 0.5`,
  `CI.type = "rubbish"`, `conf = 95`, `plot.stat = "Nonsense"`,
  `CI.type = "bca"` with `nboots = 30`, and a single-population object
  each stop with a message naming the argument.
- (F6) All four return shapes were produced and match the new `@return`.
- Report contract: input object byte-identical after the call, no history
  entry, results identical with plotting on and off.
- `ncpus = 2` completes on Windows with the shared estimator
  (`parallel = "snow"`).
- Baseline `tests/testthat/test-gl.report.fstat.R`: 99 assertions, all
  pass, no warnings. Every changed assertion carries an `[approved Fn]`
  comment. The baseline now pins the point estimates against
  `utils.basic.stats` rather than against literals, so it passes with or
  without PR #309 in the tree.

**One prediction did not hold.** The 24 `boot.ci` "extreme order
statistics used as endpoints" warnings the baseline raised are not an F1
symptom. They come from the two `nboots = 30` calls in the reproducibility
block: 30 replicates cannot supply the order statistics for a 95 per cent
percentile interval, whichever unit is resampled, and a correct locus
bootstrap raises them too. At `nboots = 200` there are none, which the
baseline now asserts.

**Left out of this change.** The `@family` retag from `matched reports` to
`matched report` moves `gl.report.fstat` between family cross-reference
blocks in 18 other `man/*.Rd` files. Those regenerated files are not in
this change, which carries only `man/gl.report.fstat.Rd`, to keep the PR
to one function; the committed `man/` on `dev` is already stale against
its roxygen for unrelated reasons.

Caller grep across all eight live clones plus `dartR.data`: no callers of
`gl.report.fstat` outside its own definition, its generated `.Rd`, the
`NAMESPACE` export and the family cross-links in other `.Rd` files. The
Phase A finding of zero callers is confirmed.

```json
{
  "function": "gl.report.fstat",
  "package": "dartR.base",
  "family": "report",
  "skill_version": "2.0.0",
  "commit": "ddaed27",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "BLOCKER", "confidence": "high", "rule": "none (spec axis: documented method not implemented)", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "STY3", "status": "applied", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "DAT7", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "applied", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "applied", "change": 5},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "FS5", "status": "applied", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 8},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "VRB1", "status": "applied", "change": 9},
    {"id": "F10", "severity": "LOW", "confidence": "high", "rule": "VRB4", "status": "applied", "change": 10},
    {"id": "F11", "severity": "LOW", "confidence": "high", "rule": "STY3", "status": "applied", "change": 6},
    {"id": "F12", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "applied", "change": 11},
    {"id": "F13", "severity": "LOW", "confidence": "high", "rule": "FS3", "status": "applied", "change": 11},
    {"id": "F14", "severity": "INFO", "confidence": "medium", "rule": "DAT6", "status": "note", "change": null}
  ],
  "approved_by": "Arthur Georges",
  "approved_date": "2026-09-08",
  "applied_branch": "review-gl.report.fstat",
  "applied_base": "ddaed27",
  "follow_up": "relocation of gl.report.fstat and gl.fst.pop to dartR.popgen, and the home of utils.basic.stats, deferred to a separate job",
  "overlap": {
    "sibling": "gl.fst.pop",
    "same_statistic": false,
    "detail": "gl.report.fstat computes Nei (1987) Gst (Fst biased, Fstp unbiased) plus Jost D and Hedrick G'st; gl.fst.pop wraps StAMPP::stamppFst, i.e. Weir & Cockerham theta. Fstp and theta agree to 1e-4 on possums.gl (45 pairs), 0.0012 on platypus.gl, 0.0724 on testset.gl.",
    "other_paths": ["utils.basic.stats", "dartR.popgen::utils.outflank.fst", "dartR.spatial::gl.ibd (StAMPP direct)", "dartR.sim::gl.diagnostics.sim (hierfstat::pairwise.neifst)", "dartR.captive::gl.report.kin.sets (different quantity)"],
    "callers": {"gl.report.fstat": [], "gl.fst.pop": ["dartR.popgen/R/gl.check.panel.r:60", "dartR.popgen/R/gl.check.panel.r:62"], "utils.basic.stats": ["dartR.base/R/gl.diagnostics.hwe.r:203", "dartR.base/R/gl.diagnostics.hwe.r:307", "dartR.popgen/R/gl.ld.haplotype.r:436"]},
    "recommendation": "keep both, fix gl.report.fstat, then move both to dartR.popgen"
  },
  "coverage_skipped": [
    "rendered heatmap: dendextend not installed",
    "ncpus > 1 not benchmarked",
    "dartR Google Group / GitHub issues: no browser session"
  ],
  "status": "pr-open",
  "pr": 384
}
```
