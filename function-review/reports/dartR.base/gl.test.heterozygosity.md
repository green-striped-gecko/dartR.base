# Review: gl.test.heterozygosity (dartR.base)

## Provenance

- Model: Claude Fable 5.1 (claude-fable-5-1)
- Skill: dartr-function-review v2.0.0 (Phase A, read-only review)
- Package commit: 98d3a98 (origin/dev, merged into dev_luis before the
  review; `R/gl.test.heterozygosity.r` and `R/utils.het.report.r`
  identical between the two)
- Date: 2026-09-15
- Family mode: analysis (population-genetic test with plotting; the input
  object must come back untouched, no history append)
- Datasets: platypus.gl (3 pops), testset.gl subsets (8 largest pops for
  the label/p-value comparison, 5 for paging, full set with two
  single-individual pops), testset.gs (SilicoDArT), gl.gen2fbm(platypus.gl),
  a plain `genlight` built with `new()` (no dartR metadata)
- Baseline: tests/testthat/test-gl.test.heterozygosity.R (new file; 32
  assertions, all passing at 98d3a98; defects pinned as-is and tagged
  with finding IDs; bootstraps seeded, nreps small, so numerical pins are
  seed-specific)
- Note: this function is a golden fixture for the PLT3 class (plot.out =
  FALSE once emptied the results, fixed in a018871). That defect is
  confirmed absent: results are identical with and without plotting.

## Verdicts

**Standards: Needs work** -- the datatype gate, working-directory check
and flag end conform, but the result table prints at every verbosity, two
input warnings are ungated, the "Starting" line is printed by hand after
the warnings instead of by `utils.flag.start`, and a `genlight` without
dartR flags crashes.

**Spec: Needs work** -- the observed differences equal an independent
per-locus uHe computation and results do not depend on plotting, but the
significance labels are computed from one-tailed quantiles, so
"sig @0.05" is a two-sided test at 0.10 and contradicts the p value and
the 95% interval reported in the same row; SilicoDArT data are accepted
although heterozygosity is undefined for presence/absence; fewer than two
populations crash.

What works well: the `pop.het` engine shared with
`gl.report.heterozygosity` gives the same uHe as a direct per-locus
Nei (1978) computation, and the FBM path returns byte-identical results.

## Findings

**F1 [HIGH, confidence: high] -- significance labels use one-tailed
quantiles at alpha, so "sig @0.05" is a two-sided test at 0.10 (DOC5;
statistical correctness)**
`R/gl.test.heterozygosity.r:159-167, 300-322` -- `upper1 = 1 - alpha1`
and `lower1 = alpha1` give the 5th and 95th percentiles of the bootstrap
difference; zero outside that 90% band is labelled "sig @0.05". The p
value (`:328-329`) is two-sided and the confidence interval (`:333-335`)
uses `(1 - conf) / 2` tails, so the three columns disagree.
Failure scenario: 8 largest testset.gl populations, `nreps = 1000`,
seed 7: 2 of 28 pairs are labelled "sig @0.05" with p = 0.074 and 0.092
and a 95% CI that includes zero; 2 pairs are "sig @0.01" with
0.01 < p <= 0.02 (reproduced; baseline test "F1"). On platypus.gl with
`nreps = 50` two pairs read "sig @0.01" next to p = 0.039. The plot's red
lines, which the `@details` text tells the reader to judge significance
from, sit at the same one-tailed quantiles.
Proposed change: take the label quantiles (and the plot lines) at
`alpha/2` and `1 - alpha/2`, so a label at `alpha` matches the two-sided
p value and the `conf = 1 - alpha` interval.

**F2 [MEDIUM, confidence: high] -- SilicoDArT objects are accepted
(FS4, DOC5)**
`R/gl.test.heterozygosity.r:131` -- `utils.check.datatype(x)` admits both
datatypes; `@param x` says "SNP genotypes". `pop.het` then treats 0/1
presence/absence scores as dosages (p = mean/2) and reports a
"heterozygosity" difference.
Failure scenario: `gl.test.heterozygosity(testset.gs)` returns a table
of significant "uHe" differences for presence/absence data (reproduced;
baseline test "FBM and SilicoDArT inputs").
Proposed change: `utils.check.datatype(x, accept = "SNP", verbose)`.

**F3 [MEDIUM, confidence: high] -- fewer than two populations crash with
"subscript out of bounds" (FS5)**
`R/gl.test.heterozygosity.r:170-182, 290` -- an object without
population assignments is given a single population "pop1", and then
`1:(nPop(x) - 1)` runs `1:0`, so `D[y, z]` is indexed at column 0 or 2.
Failure scenario: a one-population object or `pop(x) <- NULL` errors
after the bootstrap has run (reproduced; baseline test "F3/F4").
Proposed change: stop early with "at least two populations are needed"
when `nPop(x) < 2`, including the no-population case.

**F4 [MEDIUM, confidence: high] -- a genlight without dartR flags crashes
(DAT5)**
`R/gl.test.heterozygosity.r:186` -- `if (x@other$loc.metrics.flags$monomorphs
== FALSE)` is `if (logical(0))` when the flags are absent.
Failure scenario: `new("genlight", ...)` with populations errors with
"argument is of length zero" (reproduced; baseline test "F3/F4").
Proposed change: `if (isFALSE(x@other$loc.metrics.flags$monomorphs))`.

**F5 [LOW, confidence: high] -- table printed and warnings emitted at
`verbose = 0` (VRB1, VRB3)**
`R/gl.test.heterozygosity.r:137-157, 507` -- `print(df)` is
unconditional; the alpha and `boot.method` warnings are not gated.
Failure scenario: `verbose = 0` still prints the 8-line table (reproduced;
baseline test "platypus.gl"); scripts that assign the result get console
output they cannot silence.
Proposed change: print the table at `verbose >= 3`; gate the warnings at
`verbose >= 2`.

**F6 [LOW, confidence: high] -- legend labels swap when `alpha1 > alpha2`
(DOC5)**
`R/gl.test.heterozygosity.r:161-165, 413-451` -- the quantile levels are
swapped so `u1quantile` is always the looser line, but the legend labels
`paste("Sig.", alpha1)` keep the user's order.
Failure scenario: `alpha1 = 0.01, alpha2 = 0.05` draws the 0.05 lines
labelled "Sig. 0.01" and the 0.01 lines labelled "Sig. 0.05" (reproduced
by tracing the swap; table labels are correct).
Proposed change: swap `alpha1`/`alpha2` themselves when
`alpha1 > alpha2`, and derive the quantiles from the swapped values.

**F7 [LOW, confidence: medium] -- the result table is also saved as an
RDS named `table_<plot.file>` through `utils.plot.save`, undocumented
(DOC5, PLT2)**
`R/gl.test.heterozygosity.r:511-516` -- `@details` "Saving the plot"
mentions only the ggplot; the table lands beside it via the plot-saving
helper.
Failure scenario: `plot.file = "het"` writes `het_1_to_3.RDS` and
`table_het.RDS` (reproduced; baseline test "PLT3"); a user looking for
the table has no documentation pointing at the second file.
Proposed change: document the table file in `@details` and `@param
plot.file` (docs only).

**F8 [LOW, confidence: high] -- roxygen gaps (DOC1, DOC2, DOC7 proposed)**
`R/gl.test.heterozygosity.r:29-36, 92-93, 101` -- `plot.colors` is
documented as `gl.colors(2)` but defaults to
`gl.select.colors(ncolors = 2, verbose = 0)`; the `verbose` text reads
"progress log" and "[default NULL, ...]" instead of the DOC2 standard;
`@author` names a custodian but no `Author(s):`; `@family` sits after
`@examples` instead of after `@title`.
Proposed change: docs only, followed by `devtools::document()`.

**F9 [INFO, confidence: high] -- housekeeping (FS3, STY1)**
- `:196-204` the "Starting" line is printed by hand after the datatype
  check and the input warnings; the house idiom is `utils.flag.start`
  right after the verbosity check, before any other output.
- `:472-476` `match_call` is built and never used.
- `:453` `override.aes = list(size = 5)` on line layers (ggplot2 3.4
  deprecated `size` for lines; the warning is suppressed).
- `:300-307` label quantiles run without `na.rm` while the CI quantile
  uses `na.rm = TRUE`; a missing replicate would crash the former.
- `:232` `as.matrix(sgl[[i]])` densifies each population (DAT6, proposed
  rule); the engine needs a matrix, so this is noted, not proposed.
Proposed change: `utils.flag.start` at the top, drop `match_call`,
`linewidth` in `override.aes`, `na.rm = TRUE` on the label quantiles.

## Proposed changes

1. Compute the significance labels and the plot's red lines from the
   `alpha/2` and `1 - alpha/2` quantiles (F1). **Consequence: labels
   change for existing calls; pairs whose two-sided p lies between alpha
   and 2*alpha lose their "sig" label (2 of 28 pairs in the testset.gl
   check), and the red lines on the histograms move outward.**
2. Restrict input to SNP data with `accept = "SNP"` (F2). **Consequence:
   SilicoDArT objects now stop with the datatype error instead of
   returning a table.**
3. Stop with a clear message when the object has fewer than two
   populations (F3).
4. Tolerate a missing `loc.metrics.flags` slot (F4).
5. Print the table at `verbose >= 3` and gate the alpha/`boot.method`
   warnings at `verbose >= 2` (F5). Console output only.
6. Swap `alpha1`/`alpha2` when given in the wrong order so the legend
   labels match the lines (F6).
7. Documentation pass: table RDS file, `plot.colors` default, DOC2
   verbose text, `Author(s)`/`Custodian`, tag order; then
   `devtools::document()` (F7, F8).
8. Housekeeping: `utils.flag.start` before the checks, drop
   `match_call`, `linewidth` in `override.aes`, `na.rm` on the label
   quantiles (F9).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY -- run
- Spec: observed differences vs an independent per-locus uHe (Nei 1978
  eq. 2) on platypus.gl -- run, equal to 1e-10
- Label vs p value vs CI agreement at nreps = 1000 -- run (F1)
- PLT3 (results independent of plotting, same seed) -- run, identical
- `boot.method = "loc"`, paging with 10 pairs and `max_plots = 6`,
  `plot.file` output, `alpha1 > alpha2` -- run
- FBM path (DAT6): run on `gl.gen2fbm(platypus.gl)`; results identical
  to the dense object
- SilicoDArT, single population, no population, plain genlight -- run
- Two single-individual populations (full testset.gl): run, no crash
- Distributional correctness of the bootstrap itself: SKIPPED --
  delegated to package boot and to the gl.report.heterozygosity review
- Rendered plots: SKIPPED -- introspection stops at the returned table
  and the saved patchwork object
- Google Group / GitHub issues: searched "test.heterozygosity" in
  dartR.base and dartR issue trackers -- none found
- Callers: no sibling `dartR.*` package calls the function; dartr2shiny
  calls it with `nreps`, `alpha1`, `alpha2`, `max_plots`, `plot.theme`
  and shows the returned table (`shiny_fun/Fun_gl.test.heterozygosity.R`).
  No signature change is proposed; change 1 alters label values only.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | |
| 2 | approved | Luis | |
| 3 | approved | Luis | |
| 4 | approved | Luis | |
| 5 | approved | Luis | |
| 6 | approved | Luis | |
| 7 | approved | Luis | |
| 8 | approved | Luis | |

## Outcome

Applied on branch `review-gl.test.heterozygosity` (from origin/dev
98d3a98); PR number recorded in the machine block once opened.

- Change 1: quantiles at `alpha/2` and `1 - alpha/2` for the labels and
  the plot lines; labels read `alpha1`/`alpha2` directly. On the 8
  largest testset.gl populations at `nreps = 1000` (seed 7) no pair
  disagrees with its p value or its CI in either direction (was 2 of
  28); 10 pairs are "sig @0.01". A `@details` sentence records the
  `2 / (nreps + 1)` floor of the p value.
- Change 2: `utils.check.datatype(x, accept = "SNP")`; testset.gs stops
  with "found SilicoDArT expecting SNP".
- Change 3: `nPop(x) < 2` stops with "At least two populations are
  needed" before the bootstrap (one-population and no-population
  objects).
- Change 4: `isFALSE()` on the monomorphs flag; a `new("genlight")`
  object with populations returns the 3-row table.
- Change 5: `print(df)` at `verbose >= 3`; alpha/`boot.method` warnings
  at `verbose >= 2` (0 lines at `verbose = 0`, both warnings shown at 2).
- Change 6: `alpha1`/`alpha2` swapped when reversed; legend labels for
  `alpha1 = 0.01, alpha2 = 0.05` read "Sig. 0.05", "Sig. 0.01" against
  the matching lines (was reversed).
- Change 7: roxygen updated (`plot.file`/`plot.colors`/`alpha*`/verbose
  texts, table RDS in `@details`, `Author(s)`/`Custodian`, `@family`
  after `@title`); `devtools::document()` run, `devtools::check_man()`
  clean.
- Change 8: `utils.flag.start` right after the verbosity check ("Starting"
  is now the first line at `verbose >= 1`), hand-rolled start block and
  `match_call` removed, `linewidth` in `override.aes`, `na.rm = TRUE` on
  the label quantiles.
- Characterization test: 32 baseline assertions at 98d3a98; after the
  changes 38 assertions pass, the flipped ones tagged `[approved diff, change n]`.
  Observed differences still equal the independent uHe computation;
  results still identical with and without plotting; FBM identical to
  dense.
- NEWS.md entry added under 1.2.3 (development).
- Callers: no sibling package calls the function; dartr2shiny passes
  `nreps`, `alpha1`, `alpha2`, `max_plots`, `plot.theme` and displays the
  table, all unchanged in shape.

```json
{
  "function": "gl.test.heterozygosity",
  "package": "dartR.base",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "98d3a98",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "FS4", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DAT5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "VRB3", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "medium", "rule": "DOC5", "status": "approved", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 7},
    {"id": "F9", "severity": "INFO", "confidence": "high", "rule": "FS3", "status": "approved", "change": 8}
  ],
  "coverage_skipped": ["bootstrap distributional correctness: delegated to boot", "rendered plots: not reproducible in a test"],
  "status": "pr-open",
  "pr": null
}
```
