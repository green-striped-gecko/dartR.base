# dartR.base 1.2.3 (development)

* `gl.randomize.snps()`:
  - FLAGS: swapping the 0/2 homozygote coding of half the loci leaves the
    per-locus allele-frequency metrics (`OneRatioRef`/`OneRatioSnp`,
    `FreqHomRef`/`FreqHomSnp`, `PICRef`/`PICSnp`, `maf`, ...) out of step with
    the recoded genotypes. The locus-metric flags are now reset
    (`utils.reset.flags()`) so downstream recalculation knows to recompute
    them. (Applied in both `gl.randomize.snps.r` and the duplicate
    `gl.random.snp.r`, which define the same function.)

* `gl.plot.snp.density()` (function review):
  - PLOT CHANGE: bins containing no SNPs were absent from the plot data
    and drawn as background, although the help promised the lowest
    palette colour. Every bin from the start of a chromosome to its last
    SNP is now drawn, so SNP deserts are coloured and the fill scale
    starts at 0.
  - PLOT CHANGE: chromosomes were ordered by a plain string sort (`chr10`
    between `chr1` and `chr2`) while the help promised longest to
    shortest; the order is now alphabetical and numeric-aware, first name
    at the top, and documented as such.
  - `plot.display`, `plot.file` and `plot.dir` (with `utils.plot.save()`)
    replace `save2tmp`; a call passing `save2tmp` now stops with "unused
    argument".
  - an object whose `@chromosome`/`@position` slots are empty stops with a
    message naming the slots to fill instead of "arguments imply differing
    number of rows".
  - verbose >= 2 reports how many chromosomes each filter dropped; verbose
    >= 3 prints the per-chromosome table (SNPs, last SNP position, bins).
  - the `min.snps`/`min.length` messages say ">= 1", matching the check.
  - documentation: `@family graphics` added; "chromosome length" is
    described as the position of the last SNP.

* New function `gl.report.contamination()`: screens a genlight object for
  cross-contaminated samples from genotypes, population labels and, when
  present, the plate wells stored by `gl.read.dart()`. Tier 1 flags
  individuals whose heterozygosity is an outlier within their population
  (with the leave-one-out rare-allele burden as supporting evidence); tier
  2 reports each individual's strongest excess-kinship partner and whether
  it sits in an adjacent well. Validated on a DArT plate with seven
  contaminated samples confirmed by species-diagnostic loci and a lab
  note: 7/7 flagged, 0 false positives, source well named in 5/7.
* `gl.plot.heatmap()` (function review):
  - BEHAVIOUR CHANGE: a `matrix` was coerced with `as.dist()`, which kept
    the lower triangle only and set the diagonal to zero, so relatedness
    matrices (dartR.captive `gl.relatedness()`, `gl.run.EMIBD9()`) lost
    their self-values and asymmetric matrices lost one direction without
    a message. A matrix is now drawn as supplied; a matrix with one
    triangle entirely NA is mirrored from the other; a non-square matrix
    or one whose row and column names differ stops with a message.
    `dist` input is unchanged.
  - `verbose = 0` is silent: the default `palette.divergent` no longer
    prints the three `gl.colors()` lines (callers forwarding
    `verbose = 0`, such as `gl.report.fstat()`, inherit the fix).
  - the legend is built one row per population, so two populations
    sharing a colour no longer swap swatches.
  - an `fd` object with `x` supplied no longer stops with "arguments imply
    differing number of rows"; `x` is ignored for population-level input
    (`fd`, `gl.dist.pop()` distances) with a note at verbose >= 2.
  - individuals are matched to `x` by name: colours are drawn whenever
    every column of `D` is an individual of `x` (a subset of `x` works);
    otherwise a warning at verbose >= 1 replaces the silent drop or the
    "ColSideColors must be a character vector" error.
  - `par(mar)` is restored after the legend; a `palette_discrete` vector
    of the wrong length stops with a message; a missing dendextend stops
    with the install message instead of returning -1.
  - documentation: the function calls `utils.heatmap()` (not
    `gplots::heatmap.2`), returns that list invisibly, `legendy` default
    corrected, legend coordinates documented, unused gtools import removed.

* `gl.sim.crosses()`: the offspring object now carries the parents'
  `loc.all`, so a brood with no homozygous-alternate genotype (small broods,
  few loci, monomorphic parents) is no longer rejected by the
  content-vs-ploidy check in `utils.check.datatype()` as presence/absence
  data.

* Test suite: fixtures built by hand in the gl.sim.crosses, gl.fst.pop and
  gl2gi tests now carry SNP metadata (`loc.all`) so they pass the
  content-vs-ploidy check; six expectations that pinned since-fixed defects
  (gl.compliance.check F5/F9, the history entry leaked into gl.read.dart,
  gl.read.fasta and gl.read.vcf, gl.pcoa.plot F11, the gl.read.dart
  plate_location header overrun) now assert the fixed behaviour.

* `gl.test.heterozygosity()` (function review):
  - METHODS CHANGE: the significance labels ("sig @0.05", "sig @0.01")
    and the red lines on the histograms were taken from the alpha and
    1 - alpha quantiles of the bootstrap distribution, which is a
    two-sided test at 2 * alpha; they now use alpha/2 and 1 - alpha/2, so
    a label at alpha agrees with the two-sided p value and the (1 - alpha)
    confidence interval reported in the same row. Pairs whose p value lies
    between alpha and 2 * alpha lose their "sig" label (2 of 28 pairs on
    the 8 largest testset.gl populations at nreps = 1000).
  - the function is restricted to SNP data (accept = "SNP"); SilicoDArT
    objects now stop with the datatype error instead of returning
    "heterozygosity" differences for presence/absence scores.
  - fewer than two populations (including objects without population
    assignments) stop with a clear message instead of "subscript out of
    bounds" after the bootstrap has run; a genlight without
    `loc.metrics.flags` no longer crashes on the monomorphs check.
  - the result table prints only at verbose >= 3 and the alpha /
    boot.method warnings only at verbose >= 2 (the returned value is
    unchanged); when alpha1 > alpha2 the two levels are swapped so the
    legend labels match the lines they describe.
  - documentation: the table_<plot.file>.RDS file is documented, the
    plot.colors default corrected.

* `gl.map.interactive()` (function review):
  - `matrix` accepts `dist` objects (what `gl.dist.pop()` and
    `gl.dist.ind()` return) instead of crashing with "argument is of
    length zero".
  - BEHAVIOUR CHANGE: an individual-level matrix was reordered by row only
    (rows in alphabetical individual order, columns and coordinates in
    object order), so links joined the wrong individuals whenever the
    individuals were not already sorted. Matrices are now aligned to the
    object by row/column names when they match `popNames`/`indNames`, and
    taken in object order otherwise (with a warning at verbose >= 2 when
    names are present but do not match).
  - BEHAVIOUR CHANGE: symmetric links are drawn only for pairs with a
    non-missing value above 0 and never from a point to itself; the
    previous filter never fired, so every pair including self-links was
    drawn. Line width now follows the (standardised) value as documented,
    in addition to colour.
  - a single-population object placed its label with latitude and
    longitude swapped when `latlon` columns were ordered lat,lon.
  - an asymmetric matrix with `NA` cells crashed; missing pairs are now
    skipped in that direction and drawn in grey in the other.
  - `ind.circle.cols` with fewer colours than populations now errors
    instead of silently drawing the remaining individuals in leaflet's
    default colour; `latlon` stored as a matrix is accepted; missing
    `leaflet`/`leaflet.minicharts`/`terra`/`scales` stop with an error
    instead of returning -1; `gl.colors()` no longer prints at
    `verbose = 0`.

* `gl.report.heterozygosity()` (from Carlo Pacioni's PR #229, re-applied
  on the reviewed code):
  - the point estimates and the bootstrap replicates are now computed by
    the same helper (`pop.het_fun` in utils.het.report.r); previously two
    near-duplicate code paths existed and had drifted apart, so
    confidence intervals were bootstrapped around a slightly different
    estimator than the reported value.
  - METHODS CHANGE: unbiased expected heterozygosity (uHe), and FIS which
    derives from it, now apply Nei's sample-size correction per locus
    (2n/(2n-1) with n the individuals genotyped at that locus, as in
    GenAlEx and hierfstat) instead of one mean n per population. Values
    shift slightly wherever missingness varies across loci (on testset.gl
    up to 0.001 in uHe and 0.02 in FIS); the bootstrap already used the
    per-locus definition.
  - the standard errors of adjusted heterozygosity (Ho.adjSE, He.adjSE)
    were divided by the number of scored loci although the SD includes
    the invariant sites; they now use scored loci + n.invariant, so were
    previously too large by sqrt((n.Loc + n.invariant) / n.Loc).
  - a population with a single individual crashed both the plain and the
    bootstrap path ('x' must be an array of at least two dimensions):
    fixed, and where boot.ci cannot produce an interval the limits are NA.
  - the duplicate `std.error()` in utils.het.report.r was removed
    (utils.stats.r keeps the one definition).

* `gl.edit.recode.pop()`: `pop.recode` is now the read-only input recode
  table (loaded and applied when supplied, as documented and as the rest of
  the recode family use it) and `out.recode.file` is the output, written
  under `outpath`. Previously `pop.recode` was neither read nor applied;
  instead the function OVERWROTE it with a freshly generated identity table
  -- silently destroying a user's existing recode file -- while
  `out.recode.file` and `outpath` were ignored. Also: the locus-metric
  flags are now reset whenever `recalc = FALSE` regardless of verbosity (the
  `utils.reset.flags()` call sat inside an `if (verbose >= 2)` block); a
  missing `monomorphs` flag no longer crashes the run; the documented
  `recalc`/`mono.rm` defaults are corrected to FALSE; and roxygen
  descriptions are fixed. Callers who passed `pop.recode` to receive the
  output file must switch to `out.recode.file`.
* `gl.edit.recode.ind()`: the locus-metric flags are now reset whenever
  `recalc = FALSE`, regardless of verbosity. The `utils.reset.flags()` call
  sat inside an `if (verbose >= 2)` block, so the returned object's
  `loc.metrics.flags` were only reset (to signal that the metrics are stale
  after recoding/deleting individuals) at `verbose >= 2` -- the object's
  flag state depended on the reporting level. The `out.recode.file` is now
  written to `outpath` (it was written to the working directory, ignoring
  the computed `outfilespec`); a missing `monomorphs` flag no longer crashes
  the run; the documented `recalc`/`mono.rm` defaults are corrected to
  FALSE; and roxygen/message copy-paste errors from the population-recode
  sibling are fixed.
* `gl.add.indmetrics()`: a metadata file whose ids are a superset of (or
  only partly overlap) the genlight no longer crashes. The function subsets
  x to the matching individuals but had left the metadata frame at its full
  row count, so `ind.cov$pop_old <- x@pop` failed with "replacement has N
  rows, data has M" whenever the metadata carried any individual not present
  in x -- the very "a subset matches is fine" case the function's own
  warning describes. The metadata is now aligned to the matched individuals
  before use. Also: the duplicate-id check raises a proper
  `stop(error(...))`; the match-count message is corrected; and the roxygen
  gains a `@family` tag and a corrected datatype description (SNP or
  SilicoDArT genlight, not "genind"). Behaviour for exact-match and
  strict-subset metadata is unchanged.
* `gl.write.csv()`: now returns `invisible(NULL)` instead of a visible NULL
  (it previously printed a bare "NULL" at the console on every un-assigned
  call, even at `verbose = 0`); `outpath` is now resolved through
  `gl.check.wd()`, so a non-existent output directory falls back to
  tempdir() (with a warning at `verbose >= 1`) as the other io functions
  do, instead of failing with an opaque "cannot open the connection" error.
  Documentation now notes the SilicoDArT (0/1) coding and adopts the
  standard verbose text. The written file is otherwise unchanged.
* `gl.set.wd()`: an invalid working directory now raises a clear error and
  leaves the global working directory unchanged. Previously an invalid path
  was silently ignored -- the `dartR_wd` option was not set, yet the
  function returned the path and printed "Global working directory set to
  <path>", so a mistyped directory sent subsequent output elsewhere with no
  warning. A non-character, NULL or multi-element `wd` (which previously
  raised an opaque "invalid filename argument" / condition-length error) now
  takes the same clear error path. Behaviour for a valid directory is
  unchanged. Documentation corrected (default, verbose text, typos).
* `gl.check.wd()`: the "path does not exist" fallback warning now gates at
  `verbose >= 1` (it previously printed even at `verbose = 0`, leaking a
  stray line into the ~87 functions across the dartRverse that call
  `gl.check.wd(plot.dir, verbose = 0)`); a non-character, NA or
  multi-element `wd` now takes the documented tempdir fallback instead of
  raising an opaque error; documentation corrected (the `wd` default is
  NULL, resolving to the dartR_wd option then tempdir(); gl.set.wd
  reference and typos fixed). The returned working directory is unchanged
  for any valid character or NULL `wd`.

* `gl.dist.phylo()`: subst.model = "BH87" with the default
  pairwise.missing = TRUE crashed the R session (ape::dist.dna has no
  pairwise-deletion routine for BH87 and its C code faults on N/ambiguity
  codes, ape <= 5.8.1; this is what killed the macOS R CMD check job).
  The combination now falls back to global deletion with a warning.
  The same ape routine never writes the diagonal of the BH87 matrix it
  returns (uninitialised memory); gl.dist.phylo now sets it to 0, which
  also stops the garbage feeding within-population averages when
  by.pop = TRUE.
* `gl2fasta()`: the on.exit sink cleanup popped the caller's sink (e.g.
  capture.output or testthat) after the function had already closed its
  own; it now unwinds only the sinks it opened.
  R CMD check hygiene: `function-review/` added to .Rbuildignore, unused
  MASS import dropped from DESCRIPTION, gdsfmt declared in Suggests (used
  by test-gl2gds), utils.heatmap's hist() qualified, over-long example
  line in gl2paup.parsimony wrapped.

- gl.propShared has moved to dartR.spatial (green-striped-gecko/dartR.spatial#34), where its only callers (gl.ibd, gl.spatial.autoCorr) live.

- gl.Ho and gl.He have moved to dartR.sim (their only dartRverse caller is gl.diagnostics.sim); the reviewed versions and their test travel with them (dartR.sim PR #44).

* `gl.sim.genotypes()`: **behaviour change.** Allele frequencies are now
  estimated WITHIN each population and each population is simulated
  separately, so the population structure of the source object is carried
  through to the simulated object. Previously every individual was pooled
  into one panmictic gene pool without a check or a message: on
  `possums.gl` (10 populations) simulated mean Ho was 0.472, the pooled He,
  against a source Ho of 0.347. Per population, simulated Ho now tracks its
  own source population (for example A 0.328 -> 0.322, C 0.428 -> 0.404).
  `n.ind` is consequently the number of individuals simulated PER
  POPULATION, so a call on a 10-population object returns `10 * n.ind`
  individuals; this is documented in `@details`. A single-population source
  is unaffected.
* `gl.sim.genotypes()`: a locus with no calls -- and therefore no estimable
  allele frequency -- aborted the run inside `sample()` with "NA in
  probability vector", naming neither the locus nor the function. Such a
  locus is now returned as missing for the individuals of the population
  concerned, so the number of loci and the missingness structure of the
  source are preserved. The function ran on only one of the five packaged
  datasets (`possums.gl`); it now runs on all of them (`testset.gl` 3/255
  all-NA loci, `testset2.gl` 3/755, `platypus.gl` 6/1000).
* `gl.sim.genotypes()`: the announced `n.ind` cap is now applied. The
  function printed "Setting n.ind to <nLoc>" and then simulated the full
  requested `n.ind`; the branch contained no assignment.
  `gl.sim.genotypes(possums.gl, n.ind = 300)` on 200 loci returned 300
  individuals per the old code and now returns 200 per population.
* `gl.sim.genotypes()`: SilicoDArT input is refused
  (`accept = "SNP"`). The algorithm is diploid dosage arithmetic, so
  presence/absence data was silently returned as a ploidy-2 SNP object.
* `gl.sim.genotypes()`: the returned object records the call that created
  it as its history, in place of the two internal `gl.recalc.metrics()` and
  `gl.compliance.check()` entries it carried; `ind.metrics` now holds `id`
  and `pop`, and neither `loc.metrics` nor `loc.metrics.flags` carries the
  `array(NA, nLoc(x))` / `array(NA, 1)` column inherited from the
  compliance check. `n.ind` is validated (invalid values reached base R and
  failed there); the `n.ind` warning no longer prints at `verbose = 0`;
  `ploidy` is built per individual rather than per locus; a plain genlight
  is accepted; `verbose = 3` prints the results summary it promised; and
  the roxygen block gains `@details`, runnable `@examples` and an
  `Author(s):` line.
* `gl.sim.crosses()`: **simulated datasets produced by earlier versions
  carry fabricated correlation between loci and between parents, and
  should be regenerated.** Gametes were drawn with
  `ifelse(mmat == 1, sample(c(0, 2), mhet, replace = TRUE), mmat)`, which
  draws only `mhet` values -- the total number of heterozygous calls in
  the parent matrix -- and lets `ifelse()` recycle them across the whole
  matrix. Heterozygous calls at column-major positions congruent modulo
  `mhet` therefore always transmitted the same allele, within a parent
  and between unrelated parents. Per-locus segregation ratios were
  correct throughout, so the defect is invisible in any per-locus
  summary; the joint structure of the simulated cohort was not. Anything
  depending on that joint structure -- relatedness or kinship estimation,
  parentage assignment, LD, power analysis, Fst between simulated cohorts
  -- was invalid. Each heterozygous call now gets its own draw. Output
  changes for every call, under any seed.
  Also in this function: the documented `n` (offspring to retain) was
  read only by a warning test and is now applied, with `n = NULL`
  resolving to the documented lesser of 1000 and the brood total;
  `error.check = FALSE`, the path `@details` recommends for simulations,
  aborted with "object 'noff' not found" and now works; SilicoDArT input
  is refused instead of being crossed as though 1 meant heterozygote;
  `compliance.check = FALSE` now returns a structurally valid object
  (`ind.metrics` was assigned into a NULL and came back a bare list);
  non-positive or fractional `broodsize` and out-of-range `sexratio` now
  actually assign their documented fallbacks instead of only announcing
  them; unequal parent cohorts and mismatched locus panels raise
  informative errors; parental sequence-level locus metrics are carried
  onto the offspring; `ind.metrics` records `mother` and `father`;
  history is the call itself rather than the internal helper calls; the
  standard preamble and the history append no longer sit inside the
  optional `error.check` block; `verbose = 0` is silent; and the roxygen
  header is documented under `gl.sim.crosses` (was `gl.sim.cross`) with
  a runnable `@details` recipe and `@examples`.
* `gl.sample()`: any population with exactly one member was silently
  dropped from the sample and its slots filled with unrelated individuals
  drawn from elsewhere in the object. `sample(v, n)` treats a length-1 `v`
  as the range `1:v`, so a single-member population's index became a
  range to draw from. The returned individual count was always right,
  which is why the substitution went unnoticed. On `testset.gl` this hit
  every call, the default call included, because the two singleton
  populations set the default `nsample` to 1; over 20 consecutive runs
  those populations appeared zero times. Singleton populations now return
  their own member, repeated `nsample` times. **Sampled objects change
  for any dataset containing a population of one, seeded output
  included.** Second fix in the same call: the result was assembled with
  `do.call(rbind, ...)`, whose SNPbin fallback discarded the whole
  `@other` list, so `loc.metrics`, `ind.metrics`, `latlon`,
  `loc.metrics.flags` and `history` came back `NULL` and
  `gl.filter.callrate()` on the result failed with "incorrect number of
  dimensions". A single positive-index subset replaces the assembly, so
  individual metadata now tracks the drawn individuals (repeated rows
  included) and the FBM and SNPbin paths return the same object for the
  same call and seed. Also: all `loc.metrics.flags` are set FALSE because
  resampling individuals invalidates every locus metric computed across
  individuals (on one test object 248 of 255 stored `CallRate` values
  disagreed with the truth while the flag still read TRUE); the call is
  appended to `@other$history`; `ind.metrics$id` is kept in step with the
  renamed `indNames`; `nsample` is validated, so `nsample = 2.7` now
  errors instead of silently truncating to 2 and an over-large `nsample`
  under `replace = FALSE` names the smallest population; row order now
  follows draw order whether or not `nsample * nPop` exceeds `nInd`; and
  the roxygen block documents that `nsample` is per-population under
  `onepop = FALSE` and a whole-object total under `onepop = TRUE`.
* `gl.fst.pop()`: **the bootstrap confidence limits and p-values change
  numerically, and are now reproducible.** They were computed by
  `StAMPP::stamppFst()` inside `foreach %dopar%` on a PSOCK cluster StAMPP
  creates itself, whose worker RNG streams are never seeded from the calling
  session; six runs under the same `set.seed(99)` returned six different
  intervals and p-values of 0.02, 0.01, 0.01, 0.03, 0.02, 0.02, so no
  reported interval or p-value could be reproduced from the script that
  produced it. The locus bootstrap now runs in the calling session, so
  `set.seed()` immediately before the call fixes the result. Point
  estimates are unchanged: they still come from `StAMPP::stamppFst()` and
  match it to a maximum absolute difference of 0, and match
  `hierfstat::pairwise.WCfst()` to 1.4e-17 on filtered fixtures. Other
  changes: SilicoDArT input now errors instead of returning a
  ploidy-driven artefact in the range Fst normally occupies; `nboots = 0`,
  fractional or negative `nboots`, `percent` outside (0, 100), `nclusters`
  below 1 and single-population input now stop with dartR messages instead
  of R-internals messages; 40 or fewer bootstrap replicates now warn that
  the reported lower confidence limit is the smallest replicate rather than
  the 2.5th percentile; population pairs returning a non-finite Fst are
  named at `verbose >= 1`; the results summary promised at `verbose >= 3`
  now prints; and the help page documents the estimator as Weir and
  Cockerham's (1984) theta, states the p-value definition
  (`mean(replicates <= 0)`, one-tailed, no multiple-testing correction
  across up to 435 pairs) and describes the matrix actually returned
  instead of a `dist`.
* `gl.report.fstat()`: **every confidence interval this function has
  produced was computed on the wrong data and all CI output changes.**
  `boot::boot()` was handed an individuals-by-loci data frame, so it drew
  its resample indices over individuals, while the statistic applied those
  indices to the columns (loci). Only the first `nInd` loci of a pair could
  ever enter a replicate -- 40 of 1000 on `platypus.gl` -- and on a fixture
  whose leading loci are unrepresentative the reported 95 per cent interval
  was `[-0.0175, -0.0175]` around a value of 0.8622. The bootstrap now
  resamples loci, and the interval brackets the point estimate. Point
  estimates are unchanged. The bootstrap kernel was also a stale inlined
  copy of `utils.basic.stats()`, so the interval and the value beside it
  came from different estimators and data with loci absent from one
  population aborted the whole call; the copy is deleted and the shared
  helper is called. SilicoDArT data is now refused (`accept = "SNP"`)
  instead of being read as SNP dosages. Arguments are validated:
  `nboots = 1`, negative or fractional `nboots`, an unknown `CI.type` or
  `plot.stat`, `conf` outside (0, 1), fewer than two populations, and
  `CI.type = "bca"` below 200 replicates now stop with an informative
  message. `verbose = 0` is fully silent and draws no heatmap; the results
  summary moved from `verbose >= 2` to `verbose >= 3`; dropped
  single-individual populations are named from `verbose >= 1`. The
  `nboots = 0` return with more than two populations names its second
  element `Stat_tables` and no longer prefixes its columns with
  `Stat_tables.`.
* `gl.pcoa.plot()`: three documented parameters were accepted and then
  ignored, and one documented save location was wrong. All four now behave
  as the help page describes, so existing calls can produce different
  output. (1) `hadjust` and `vadjust` were validated and never used; they
  now set the horizontal and vertical justification of the point labels, so
  labels move -- including on default calls, since the defaults
  `hadjust = 1.5`, `vadjust = 1` are not the neutral 0.5. Pass
  `hadjust = 0.5, vadjust = 0.5` for the previous centred labels. (2)
  `pop.labels = 'ind'` passed validation but no branch built the plot, so
  the call died with "object 'plott' not found"; it now labels each point
  with its individual name. (3) `plot.file` with the default `plot.dir`
  wrote the RDS into the current working directory; the function now
  resolves `plot.dir` through `gl.check.wd()`, so the file lands in
  `tempdir()` (or the directory set by `gl.setwd()`) unless `plot.dir` is
  given. (4) The axis-range check reset an out-of-range `yaxis` to the
  constant 2 and `zaxis` to 3 without checking those against the ordination,
  so any ordination holding fewer axes -- `gl.pcoa(gl, nfactors = 1)`, a
  two-individual object, `zaxis = 5` on a two-factor ordination -- died with
  "subscript out of bounds". Axis choices are now bounds-checked, must
  differ from one another, and too-small ordinations stop with a message
  naming the shortfall.
* `gl.pcoa.plot()`: conformance and message fixes. A `plot.display`
  argument (default TRUE) was added and, per the house rule that
  `verbose = 0` is fully silent, `verbose = 0` no longer displays the plot;
  callers that relied on `verbose = 0` to silence messages while still
  showing the plot (including `gl.assign.pca()` in dartR.captive and
  dartR.popgen) will no longer see it. Twelve messages that printed at
  `verbose = 0` are now gated at `verbose >= 2`. A PCoA of a distance matrix
  corrected with `cailliez` or `lingoes` was labelled "PCA Axis" because the
  classification keyed on `$loadings` being absent; it is now identified
  from the shape of the loadings and labelled "PCoA Axis". `@return`
  documented NULL although a visible ggplot is returned. The `directlabels`,
  `plotly`, `gganimate` and `tibble` guards returned `-1` instead of
  stopping. An ordination whose entities do not match the genlight object
  (for example one built from `gl.dist.pop()`) failed with "arguments imply
  differing number of rows" and now stops with a message that names the
  mismatch. `pt.colors` and `pt.shapes` are honoured in the
  `interactive = TRUE` branch and documented as unavailable for shapes in
  3D. The `as.pop` error named `loc.metrics` where the lookup is in
  `ind.metrics`; the `legend` branch rendered its legend titled "pop"
  despite mapping `Population`; a `plot.theme` argument
  (default `theme_dartR()`) was added.
* `gl.report.shannon()`: **results change for `level = 'beta'` and
  `'gamma'`**. Both returned degenerate constants (beta identically 1,
  gamma identical to alpha, one row per individual) because a single
  individual's dosage vector was fed to machinery expecting an abundance
  matrix. They now compute the real population-level partition of Ma, Li &
  Zhang (2020): each population's individuals form an abundance matrix,
  gamma is the diversity of the pooled locus abundances, alpha the mean
  within-individual diversity, and beta = gamma/alpha the effective number
  of distinct individuals; the partition gamma = alpha x beta holds at
  every order and the output is one row per population (`pop`, q0 ...).
  `level = 'alpha'` results are unchanged (verified byte-identical). Also:
  `verbose = NULL` default, so the global `gl.set.verbosity()` setting is
  honoured; `verbose = 0` is now fully silent (reshape2 melt message
  removed, plot gated off); `level` and `order` are validated with
  informative errors instead of failing deep in the loop; a malformed,
  dead dependency guard removed; individuals with no non-missing,
  non-zero dosages return NA rows (was 0/1/Inf) with a gated warning
  naming them; roxygen header rewritten to describe the actual
  computation (per-individual Hill-number profiles q0, q1, q2 ...,
  natural-log convention, and the population partition semantics).
* DESIGN CHANGE - genome-only @position. The genlight @position and
  @chromosome slots are now reserved for GENOME coordinates; the
  position of the SNP within the sequence tag lives solely in
  gl@other$loc.metrics$SnpPosition (0-based). Concretely:
  - utils.dart2genlight()/gl.read.dart() no longer copy SnpPosition
    into @position: DArT-read objects carry @position = NULL (and
    @chromosome = NULL) until genome coordinates are assigned
    explicitly (e.g. x$position <- x$other$loc.metrics$ChromPos_...).
  - gl.compliance.check() no longer fills @position from SnpPosition;
    instead it CLEARS a provably stale copy (when @position is
    identical to SnpPosition), with a gated message. Genuine genome
    coordinates are never touched.
  - gl2bpp() now reads SnpPosition from loc.metrics (as gl2fasta
    already did) - output verified byte-identical; it previously read
    @position, which silently corrupted haplotypes once users assigned
    genome coordinates to the slot as documented elsewhere.
  - gl2vcf()/gl2hapmap() replace the fragile max(position) < 1000
    heuristic with an explicit NULL-slot test.
  - The show method no longer claims @position is "[within 69 base
    sequence]".
  Migration note: scripts that read tag positions from @position
  should use gl@other$loc.metrics$SnpPosition instead.
* `gl.pcoa()`: five fixes from the function-review campaign (findings
  verified empirically from an external code-read handed over by the
  custodian). (1) `nfactors` was silently ignored on the file-backed (FBM)
  big_SVD path -- scores/loadings came back with `nInd - 1` columns;
  they are now truncated to `nfactors` (clamped with a warning when fewer
  axes exist), `$eig` stays full length as on the glPca path. (2) The
  distance-matrix path crashed with "subscript out of bounds" when the
  matrix had fewer than `nfactors + 1` entities; `nfactors` is now clamped
  to the available axes with a warning. (3) The `verbose >= 3` summary
  printed "NA % of the total variance" lines when fewer than 2-3
  informative/positive axes existed (the dist-branch guard
  `length(eig.top >= 2)` tested the length of a logical, a precedence
  bug); axis-combination lines now print only when enough axes exist.
  (4) A dead assignment to `e` (immediately overwritten) was removed.
  (5) `@details` now documents the second algorithm: FBM-backed objects
  use `bigstatsr::big_SVD()` after neighbour imputation, in-memory
  objects use `adegenet::glPca()`. The FBM reconstruction rescaling and
  the Tracy-Widom criterion remain open custodian items.
* `gl.dist.phylo()`: numerics unchanged at defaults; the documentation
  now states what the numbers mean -- heterozygous sites carry no
  distance signal (heterozygotes are written as IUPAC ambiguity codes,
  which `ape::dist.dna()` treats as missing data), and a warning at
  `verbose >= 1` reports the fraction of heterozygous genotype calls,
  with a stronger warning when `pairwise.missing = FALSE` deletes every
  het-bearing site globally. The gamma corrections and variances the
  details promised are now requestable (new `gamma` and `variance`
  parameters forwarded to `ape::dist.dna()`; defaults reproduce the
  previous output exactly). Also: the `setwd()` round-trip is gone, so a
  mid-pipeline failure can no longer strand the session working
  directory in `tempdir()`; SilicoDArT input is rejected up front with
  the redirect to `gl.dist.pop`/`gl.dist.ind` (the previous check
  compared against the wrong case and never fired); a single-population
  object with `by.pop = TRUE` and an all-monomorphic object now fail
  with informative errors instead of "subscript out of bounds" and a
  bare subsetting message; the `by.pop = FALSE` label format
  (indName_pop), the BH87 asymmetric-matrix return and the
  TrimmedSequence/SnpPosition/loc.all prerequisites are documented.
* `utils.check.datatype()`: residual findings of the second-pass review
  applied. (1) The courtesy all-NA scan no longer densifies an FBM-backed
  object at `verbose >= 2`; it reads the FBM in column blocks instead of
  materialising the full genotype matrix. (2) Content-vs-ploidy
  consistency: an object of uniform ploidy 2 whose non-missing genotypes
  are all 0 or 1 AND that carries no SNP metadata (empty `loc.all` slot,
  no `SNP`/`SnpPosition` locus metrics) -- presence/absence content
  mislabelled as SNP -- is now rejected with an actionable fatal at any
  verbosity, instead of passing `accept = "SNP"` gates into dosage-based
  0/1/2 arithmetic. Clean SNP objects are unaffected: the scan exits at
  the first genotype of 2, and a subset that merely lacks the
  homozygous-alternate class is vouched for by its `loc.all` slot.
  (3) A data.frame is now classified `"data.frame"` rather than `"list"`.
  (4) The classification-affecting mixed-ploidy notice prints from
  `verbose >= 1` (was `>= 2`); the accept-gate fatal gained its missing
  trailing newline; the ploidy-slot-only classification contract, the
  case-sensitivity of `accept`, and the data.frame mapping are now
  documented; an `Author(s):` line was added to the header.
* `gl.allele.freq()`: two frequency corrections in the `by = 'loc'`
  breakdown -- `percent = TRUE` was silently ignored (the frequency
  column was always a proportion; it is now rescaled to the percentage
  scale), and for SilicoDArT presence/absence data the column was
  exactly half the presence frequency (the raw 0/1 matrix was divided
  by the SNP ploidy divisor; `simple = TRUE`, which routes through
  `by = 'loc'`, returned the same halved values in `alf1`/`alf2`).
  The `by = 'popxloc'` and `by = 'pop'` breakdowns are unchanged
  (verified byte-identical), so downstream consumers such as
  `gl.dist.pop()` are unaffected. Also: an unrecognised `by` now stops
  with an informative error instead of silently returning the full
  population x locus table; a genlight not built by dartR (no
  `loc.metrics.flags`) no longer crashes on the monomorphs-flag access;
  documentation corrected (`@return` is a data.frame; the
  `simple = TRUE` override of `percent`/`by` and the `by = 'loc'`
  column semantics are now stated).
* `gl.read.silicodart()`: four fixes. (1) A metafile without a pop
  column always crashed with "object 'out' not found" (a variable-name
  typo - the branch had never worked); it now mirrors the SNP path and
  defaults every individual to 'pop1' with a gated warning. (2)
  Duplicate non-numeric CloneIDs produced NA locus names (the
  uniquification loop assigned new strings into a factor column); the
  column is now converted to character first, so duplicates get the
  promised _1, _2 suffixes. (3) All messages are verbosity-gated, raw
  cat() calls routed through the colour helpers, and verbose is passed
  to the closing gl.compliance.check: `verbose = 0` is silent
  (previously 31 lines printed). (4) Individuals absent from the
  ind.metafile are still removed (the contract aligned with the SNP
  path), but the removal is now announced with a warning stating how
  many individuals were REMOVED, and the message no longer mislabels
  them "loci" ("Subsetting loci now!").
* `gl.compliance.check()`: two fixes from its function review, applied at
  a deliberately minimal scope. (1) The genotype-coding check was vacuous
  whenever the data contained a missing value (`max(mat)` without na.rm is
  NA, and `NA %in% c(0, 1, 2, NA)` is TRUE), so out-of-range genotypes
  passed as "confirmed"; the check now tests the value set exactly.
  Message gating and warn-not-stop semantics are unchanged: violations are
  reported at `verbose >= 1`, not repaired and not fatal. (2) Input in
  which every locus is monomorphic crashed the monomorph check
  ("Subsetting resulted in zero loci" via `gl.filter.monomorphs()`), and
  input in which every locus is all missing crashed the all-NA check the
  same way; both cases now complete with the relevant flag set FALSE and
  a gated message. With every locus all missing and `verbose >= 2`, the
  same root cause still errors earlier inside `utils.check.datatype()`
  (pre-existing, recorded in the review report). The review's remaining
  findings are recorded as deferred in
  `function-review/reports/dartR.base/gl.compliance.check.md`.
* `gl.tree.nj()`: THE DEFAULT TREE CHANGES on any dataset with missing
  data -- the default Euclidean distance is now computed from
  per-population allele frequencies with missing genotypes excluded
  (na.rm = TRUE), matching the policy of `gl.dist.pop()`; previously a
  single missing genotype nullified a population's frequency at that
  locus (23.8% of the frequency matrix on `testset.gl`, topology shift
  RF 32). Frequencies are also scaled by ploidy, so Tag P/A
  (SilicoDArT) branch lengths are no longer halved. Other repairs:
  `by.pop = FALSE` (individual-level tree) crashed on any ordinary
  dataset ("Vector length does no match number of populations") -- now
  functional; `method = "UPGMA"`/`"upgma"` is now accepted (previously
  only the misspelling "ugpma" selected UPGMA and the correct spelling
  silently returned an nj tree; "ugpma" is retained as a synonym), and
  the unknown-method fallback warning is gated (was printing at
  `verbose = 0`); plotting is decoupled from the result (new
  `plot.display` argument, forced FALSE at `verbose = 0`, and a plot
  failure no longer loses the computed tree); `type` is validated up
  front; outgroup rooting resolves the root so `ape::is.rooted()` is
  TRUE; a `dist.matrix` supplied as a plain matrix is coerced with
  `as.dist()` (the ugpma path previously failed opaquely); the `as.pop`
  error message now points at `ind.metrics` (was `loc.metrics`);
  documentation corrected (return value, parameter order, method
  spelling, plot types).
* `gl.read.dart()`: `verbose = 0` no longer leaks the full
  `gl.compliance.check()` progress log; the internal `utils.recalc.maf()`
  result is now assigned (previously discarded -- the final object was
  rescued by the compliance check); the lastmetric autodetection now runs
  after the standard preamble, checks that the file exists and fails with
  a clear message (naming the `lastmetric` parameter) when the file has no
  '*' header rows, instead of dying pre-banner with "argument of length
  0"; and DArT reports lacking CallRate/AvgCountRef/AvgCountSnp no longer
  crash the read-depth calculation ("replacement has length zero") --
  `rdepth` is set NA with a gated warning.
* `utils.vcfr2genlight.polyploid()` (the genotype-conversion engine of
  `gl.read.vcf()`): three fixes. (1) Genotype mode coded half-missing
  calls (`0/.`, `./1`) as heterozygous -- a call with no observed ALT
  allele returned 1, inflating heterozygosity on low-coverage vcfs;
  dosage mode handled the same calls inconsistently (1 and 0). Any call
  still carrying a "." after separator stripping is now NA in both
  modes. (2) Ploidy was whatever adegenet inferred from each
  individual's maximum dosage (a triploid individual with no 1/1/1 call
  was stamped diploid; genotype-mode objects came back haploid/diploid
  mixtures); ploidy is now derived from the allele count (arity) of the
  GT calls, which is exact. (3) Omitting `mode2` failed with an
  unrelated closure-coercion error (the signature default `mode2 = mode`
  resolved to `base::mode`); the default is now "genotype", matching the
  caller.
* `utils.dart2genlight()` (affects `gl.read.dart()` imports): five fixes.
  (1) The ind.metafile acts as a filter as well as a metadata source -
  individuals in the DArT file without a metafile row are removed; that
  behaviour is retained but is now documented in the help and announced
  with a warning stating how many individuals were REMOVED (previously
  the removal was silent beyond "Maybe this is fine if a subset
  matches"). (2) A DArT file lacking a SNP/Variant column now fails fast
  with a clear header error instead of "cannot coerce type 'closure'".
  (3) The two-row genotype translation table is exhaustive: 0/0 and
  partial-missing pairs map to NA deliberately (no more spurious "NAs
  introduced by coercion" warnings) and any unrecognised pattern in a
  corrupted file is counted and reported. (4) Duplicate ids in the
  ind.metafile abort with a proper error message (previously a bare
  stop() with an empty condition message). (5) The TrimmedSequence and
  id-mismatch warnings are verbosity-gated (the unmatched-id listing
  moved to `verbose >= 2`) and terminate with newlines; `verbose = 0`
  is silent.
* `gl.read.csv()`: five fixes from function review. (1) The documented
  `loc.metafile` feature crashed on every use ("argument is of length
  zero"); the metrics file's own AlleleID column is now validated against
  the locus names of the input data and then attached. (2) Out-of-range
  numeric genotype codes (e.g. a stray 5) were admitted silently and
  survived into the returned object; they now raise a fatal error naming
  the offending values. (3) Files with fewer than 5 loci or 5 individuals
  crashed on the unclamped 5x5 type-sniff window. (4) An ind.metafile
  with the id column anywhere but first was spuriously rejected; the id
  column is now compared by name, as documented. (5) Character-data
  allele pairs are now stored in `loc.all` (ref = most frequent allele)
  instead of the uniform "A/C" compliance placeholder, so downstream
  exports (gl2vcf, gl2fasta, gl2plink) see the real allele spellings.
* `gl.read.vcf()`: six fixes. (1) In dosage mode every individual was
  stamped diploid even though the matrix holds copy numbers up to the
  data's ploidy (a triploid vcf returned dosages of 3 under
  `ploidy == 2`); the ploidy is now set from the data's maximum copy
  number per the documented dosage semantics (uniform across
  individuals, as `gl.compliance.check()` requires a single ploidy
  level), and ploidy 2 is forced only in genotype mode. (2) INFO fields were split
  positionally with column names taken from the first record, so any vcf
  whose INFO keys vary in order or presence across records (merged or
  multi-caller vcfs) got silently swapped locus metrics; INFO is now
  parsed per record by key. (3) `verbose = 0` printed a full compliance
  transcript (~48 lines) because `gl.compliance.check()` and
  `gl.recalc.metrics()` were called without the user's verbosity; it is
  now passed through, and the metafile id-mismatch warnings are gated at
  `verbose >= 1`. (4) The presence of one multi-allelic record coerced
  every INFO column to numeric, destroying FILTER ("PASS" -> NA) and any
  character INFO field for all retained loci; only numeric-parseable
  columns are now coerced. (5) Individuals absent from the ind.metafile
  were silently dropped from the returned object; they are now retained
  with NA metadata and listed in a warning at `verbose >= 1`. (6) Haploid
  or polyploid calls read in genotype mode were silently recoded onto the
  diploid scale; a warning now states this at `verbose >= 1`, and the
  dosage orientation (counts ALT -- the opposite of `gl.read.PLINK()`,
  which counts allele.2) is documented.
* `gl2related()`: (1) The exported file no longer wraps individual names
  in quote characters (`write.table` ran with its default `quote = TRUE`),
  which COANCESTRY -- the file's documented consumer -- would read as part
  of the identifier; file output changes for every caller. (2) SilicoDArT
  input is now rejected (`accept = "SNP"`): presence scores were silently
  written as heterozygote pairs, making relatedness estimates meaningless.
  (3) The allele coding (1/3 alleles, 0 missing, two columns per locus,
  tab-separated, no header) is documented in `@details`, and the
  verbose/author/return roxygen brought to standard.
* `gl2hiphop()`: SilicoDArT input, which was silently recoded into 0/2
  pseudo-genotypes ('heterozygote' in hiphop coding) and would score
  meaningless parentage mismatches downstream, is now rejected
  (`accept = "SNP"`). The docs now state the dartR -> hiphop coding map
  (0 -> 0, 1 -> 2, 2 -> 1, NA preserved); the recode is numeric,
  dropping the character round trip through dplyr.
* `gl2demerelate()`: SilicoDArT (presence/absence) data was silently
  accepted and recoded as fake diploid genotypes (every present tag became
  a 1/2 "heterozygote"), feeding Demerelate meaningless relatedness input;
  the function now stops with an error unless the data are SNP. Roxygen
  header completed (the 1/2 allele coding and NA-as-missing convention are
  now documented) and a dead assignment removed.
* `gl2paup.svdquartets()`: (1) SilicoDArT output is written with
  `format datatype = standard` (0/1 symbols) instead of `dna`, which PAUP
  rejects for presence/absence characters. (2) Data is sorted by
  population once, ahead of the ploidy branch, so SilicoDArT taxpartition
  ranges correspond to the matrix rows; `gl.sort` runs quietly and the
  monomorph warning is gated, making `verbose = 0` fully silent. (3) The
  monomorph check tolerates genlight objects without dartR
  `loc.metrics.flags`. (4) Single-population input no longer crashes on
  the taxpartition construction. (5) `method = 0` (or any value other
  than 1/2) is coerced to 2 with a warning instead of silently accepted,
  and the documented "method=2 is assumed" for presence/absence data is
  now enforced (silico with `method = 1` previously crashed). (6) The
  PAUP log/tree filenames inside the paup block are derived from
  `outfile` instead of hardcoded `svd.*`. (7) Individual names are
  sanitised for parentheses like population names. (8) The return is
  `invisible(NULL)`.
* `gl2bpp()`: seven fixes. (1) On plain genlight input (any object built
  by adegenet/vcfR rather than dartR) the internal position sort did not
  carry `loc.metrics`, so locus names were paired with the wrong
  TrimmedSequence and the alignment was genetically wrong with no warning;
  the sort now re-subsets the metadata explicitly and dartR-class and plain
  genlight input produce identical output. (2) The `merge.secondaries`
  block produced a structurally invalid file (stale block-header lengths,
  orphan headers when a clone had three or more secondaries, a duplicated
  segment between the last two SNP positions, and labels read from the
  wrong lines); it is rewritten block-wise and verified on a synthetic
  secondary pair. (3) `method` is now validated -- an invalid value
  previously ran to "Completed" and silently wrote no alignment.
  (4) An object without `loc.metrics.flags` crashed with "argument is of
  length zero"; absent flags are now treated as monomorphs-not-confirmed
  (gated warning). (5) A vestigial dependency guard on seqinr (never used
  by the function; non-stopping `return(-1)` idiom) is removed. (6) Output
  lines no longer carry stray leading/trailing spaces. (7) The return is
  `invisible(NULL)` and the roxygen header is brought to the house
  template.
* `gl.mahal.assign()` REMOVED (family consolidation): one of three
  near-duplicate Mahalanobis assignment functions across the verse
  (with dartR.captive's `gl.assign.mahal()` and
  `gl.assign.mahalanobis()`). The reviewed and corrected
  `gl.assign.mahalanobis()` (envelope-dimension cap, MASS::ginv
  pseudo-inverse, rank-based chi-square df) becomes the single
  implementation, arriving with the captive-to-popgen assignment-suite
  migration. No callers existed in the family.
* `gl2gds()`: two silent data corruptions fixed. (1) With
  `snp.pos`/`snp.chr` supplied, the genotype/snp.id/snp.allele records
  were reordered with the INVERSE of the chromosome/position sort
  permutation while the coordinate columns took the forward-sorted
  values, so almost every record in the written gds carried another
  locus's position and chromosome (7 of 8 in the platypus
  demonstration); one permutation now reorders every per-locus field
  together, and each record's coordinates map 1:1 to the supplied
  loc.metrics fields (verified with SNPRelate across all 1000 platypus
  loci). (2) The stored genotype was the raw genlight dosage (count of
  the second allele) although SNPRelate defines it as the count of the
  first allele of snp.allele; the dosage is now written as 2 - dosage,
  so `snpgdsSNPRateFreq` and every ref/alt-aware export report correct
  allele frequencies (PCA/IBD were unaffected). Also: the gds structure
  dump no longer prints at `verbose < 3`, and SilicoDArT input is
  rejected (`accept = "SNP"`) instead of being written as a pseudo-SNP
  gds.
* `gl2phylip()`: (1) With `bstrap > 1` the function returned the LAST
  bootstrap replicate's distance matrix; it now returns the observed-data
  matrix as documented. (2) Bootstrap resampling subsetted the genlight
  with duplicated locus indices, hitting the adegenet SNPbin defect that
  silently converts missing genotypes to homozygous reference; loci are
  now resampled at the matrix level, so replicate values change on data
  with missing genotypes. (3) Allele frequencies are computed with
  `na.rm = TRUE` -- previously any missing genotype voided the whole
  population x locus cell (35% of cells on testset.gl), leaving each
  pairwise distance resting on an undisclosed locus subset; numerical
  output changes for datasets with missing genotypes, and a warning
  reports cells with no scored genotypes. (4) SilicoDArT input is now
  rejected (`accept = "SNP"`): the diploid-dosage math halved all
  presence/absence frequencies. (5) `outpath` follows the family idiom
  (`NULL` resolved via `gl.check.wd`). (6) The output sink is released
  on error. (7) A warning is given when 10-character truncation collides
  population names. (8) Build tag and roxygen conformance.
* `gl2paup.parsimony()`: (1) a bootstraps/ncpus count that does not divide
  evenly is now a fatal error instead of a printed "Fatal Error" that let
  the run continue and write a fractional bootstrap count into the PAUP
  block. (2) The three bash-mode generator scripts are written to
  `outpath` alongside the nexus file, not to the current working
  directory. (3) The generated PBS jobs cd into `base.dir.name` (a
  personal directory was hardcoded) and take their project code from the
  new `pbs.project` parameter [default 'xl04']; the storage directive
  lists that project, the project holding `base.dir.name`, and gdata/if89
  (PAUP's home). (4) The consensus script gains its missing shebang, all
  #PBS directives sit directly under it, and its job name no longer
  contains a literal `${i}`. (5) Single-population input no longer
  crashes on the taxpartition construction. (6) Advisory and coercion
  warnings are verbosity-gated and the return is invisible, so
  `verbose = 0` is fully silent. (7) The examples use the real
  `outfileprefix` parameter, and the test-mode subset respects a
  user-supplied `ncpus`.
* `gl2plink()`: (1) When `@chromosome`/`@position` are unset, the .map no
  longer fabricates coordinates (chromosome "1", positions 1..nLoc);
  unmapped loci are written with PLINK's own convention, chromosome 0 and
  position 0, so downstream position-aware operations (LD pruning,
  clumping) cannot silently run on fictitious coordinates. The fallback
  is documented, along with the hazard of within-tag offsets in
  `@position`. (2) `sex.code` is recoded to PLINK's 1/2/0 codes
  unconditionally -- a scalar (including the default "unknown") was
  previously written verbatim into the .ped sex column, so an all-female
  cohort supplied as `sex.code = "F"` silently lost its sex information.
  (3) SilicoDArT input is rejected at the datatype gate
  (`accept = "SNP"`) instead of crashing mid-write. (4) With
  `bed.files = TRUE`, the written allele list (now named
  `<outfile>_a2_alleles.txt`) is passed to PLINK via `--a2-allele`, so
  the .bed/.bim ref/alt orientation matches the genlight instead of
  being reassigned by minor-allele frequency. (5) The fallback warnings
  are gated at `verbose >= 1` (they printed at `verbose = 0`).
  (6) `@description` states the two-file default (.ped embeds the fam
  columns) rather than promising bed/bim/fam unconditionally.
  (7) Roxygen conformance and re-documentation.
* `gl2hapmap()`: the documented `pos` and `chrom` arguments were silently
  ignored whenever the corresponding slot (`@position`, `@chromosome`)
  was already populated -- a nominated loc.metrics field now takes
  precedence over the slot, so output positions change for callers who
  passed these arguments against an object with populated slots. A
  genlight without allele definitions (`loc.all` NULL) now fails fast
  with a clear message instead of an opaque replacement-length error.
  The zero-fill convention for positions is documented and messaged at
  `verbose >= 2`. The saved-file message printed at `verbose = 0` (now
  gated at `verbose >= 2`) and the visible NULL return auto-printed
  (now invisible). The SilicoDArT rejection raised a condition with an
  empty message (`cat(error()); stop()`); it now goes through
  `utils.check.datatype(accept = "SNP")`.
* `gl2faststructure()`: SilicoDArT (presence/absence) objects were
  silently admitted and written as meaningless pseudo-diploid genotypes;
  the function now accepts SNP data only and stops otherwise. The NULL
  return is now invisible (no stray "NULL" printed on unassigned calls).
* `gl2treemix()`: SilicoDArT objects were admitted with each ploidy-1
  individual counted as two allele copies, understating drift throughout
  a treemix run; the function now stops with a clear datatype error
  (`accept = "SNP"`). Also: the gz write restores the console and closes
  its connection if it fails midway; the `verbose > 2` record count now
  reports loci (and populations) instead of individuals; the NULL return
  is now invisible. SNP output is unchanged.
* `gl2gi()`: df2genind silently dropped loci with no scored alleles while
  the full `loc.metrics` table was copied wholesale, so converting data
  containing all-NA loci returned a genind whose locus metadata was out
  of register with the genotypes (and `gi2gl(gl2gi(x))` crashed on such
  data); all-NA loci are now removed up front (gated warning) and
  `loc.metrics` is subset to the surviving loci. SilicoDArT input, which
  previously converted without warning to a meaningless ploidy-2 genind,
  is now rejected (`accept = "SNP"`). Single-locus objects no longer
  crash ("X is not a matrix"). A warning is issued when placeholder
  alleles are fabricated for objects without allele definitions. The
  genotype recode is vectorised; the progress message is gated at
  `verbose >= 2`; docs corrected (`probar` default is FALSE).
* `gl2genalex()`: SilicoDArT (presence/absence) objects were silently
  exported as meaningless codominant GenAlEx genotypes; the function now
  accepts SNP data only. A missing poppr installation raises a proper
  error condition instead of returning -1 (scripts relying on the -1
  return now see an error). Also: poppr's conversion chatter is
  suppressed below `verbose = 2`, all-NA loci dropped before export are
  reported at `verbose >= 1`, the NULL return is invisible, and the
  `overwrite = FALSE` doc now states the actual raise-on-exists
  behaviour.
* `gl2structure()`: four fixes. (1) A re-run with
  `export.marker.names = FALSE` appended to the existing file instead of
  overwriting it, silently doubling every individual for STRUCTURE; the
  file is now truncated whenever no marker-name header precedes the data.
  (2) SilicoDArT objects were admitted and exported as fabricated diploid
  genotypes; the function now stops with a clear datatype error
  (`accept = "SNP"`). (3) Duplicated locus names collapsed genotype
  columns (header listed all loci, data rows carried fewer, misaligning
  every column after the collision); columns are now assigned
  positionally. (4) The NULL return is now invisible. Docs corrected: the
  object needs no lat/long location data, and the `ind.names` default is
  `indNames(x)`.
* `gl2bayescan()`: SilicoDArT (presence/absence) data was silently accepted
  and its tag counts doubled into a codominant (diploid gene-copy) BayeScan
  file; the function now stops with an error unless the data are SNP.
  Population x locus combinations with no genotyped individuals (written as
  samples of zero gene copies) are now counted and reported with a warning
  at `verbose >= 1`, with a recommendation to filter on call rate before
  export. The console sink is protected with `on.exit()` so an error during
  writing no longer leaves the session's console redirected, and the
  function returns NULL invisibly (an unassigned call no longer prints
  NULL).
* `gl2genepop()`: an all-monomorphic object wrote a silently malformed
  file (the header claimed all loci but every row carried exactly two
  "0000" fields); such objects now write one correctly coded genotype
  per locus. `pop.order` entries that were unlisted, misspelled or
  duplicated silently dropped populations from the exported file;
  pop.order is now validated against the object's population names and
  errors naming the offending entries. Also: the save-path message is
  gated at `verbose >= 2` (it printed at verbose 0, with a stray path
  separator), the per-locus column lookup is hoisted out of the
  individual loop (was O(nInd x nLoc^2)), the SilicoDArT rejection
  raises a proper error condition (was an empty stop()), a
  single-individual object fails fast with a clear message (was "X is
  not a matrix"), and all-NA locus drops are reported at
  `verbose >= 1`.
* `gl2snapper()`: the autapomorphy classifier (`rm.autapomorphies = TRUE`)
  now scores a population as polymorphic from allele counts (both alleles
  observed) instead of thresholds on rounded frequencies, so near-fixed
  polymorphic populations no longer cause informative loci to be dropped;
  the set of retained loci can change. Also: preprocessing chatter
  (name-mangling warnings, gl.allele.freq/gl.drop.loc progress) is gated
  so `verbose = 0` is silent; the nexus write restores the console if it
  fails midway; the NULL return is now invisible; the documented
  `rm.autapomorphies` default corrected to FALSE (the signature default,
  which is unchanged).
* `gl2gapit()`: now actually writes the hapmap table to
  `outpath/outfile` (tab-delimited) -- previously the parameters and the
  progress message promised output files but nothing was written. The
  roxygen header belonged to `gl2geno` (so `?gl2gapit` found nothing); a
  genuine header and man page now document the invisible hapmap
  data.frame return. Chromosome names are recoded to stable integer
  codes (alphabetical order of the distinct names, mapping reported at
  `verbose >= 2`) instead of unstable factor level indices; the
  `assembly` column is NA (was hardcoded "Oilpalm"); SilicoDArT input is
  rejected with a clear datatype error (previously an opaque crash); the
  empty-slot warnings are verbosity-gated; and the genotype mapping is
  vectorised.
* `gi2gl()`: a genind containing any locus with more than two alleles was
  silently corrupted -- the tab-column walk assumed at most two alleles
  per locus, so every locus after a multiallelic one read the wrong
  column; such input now raises a fatal error naming the offending loci.
  A single-locus genind crashed ("argument is of length zero"); fixed.
  Non-diploid genind objects are now rejected (their allele counts
  cannot be represented as 0/1/2 dosages). Allele spellings are now
  reconstructed from the genind allele names instead of a uniform
  placeholder, consistent with the returned dosages. The function's
  documentation page was merged into `gl2gi`'s (`@name gl2gi` in
  gi2gl.R); `gi2gl` now has its own manual page with examples.
* `gl2geno()`: the verbose >= 1 output-file message printed a single
  garbled name (`<outfile>.geno.lfmm.`); it now prints the two real
  output paths.
* `gl2bayesAss()`: SilicoDArT (presence/absence) data was silently accepted
  and written as fake diploid BA3 genotypes; the function now stops with an
  error unless the data are SNP. The `ploidy != 2` refusal now uses the
  standard fatal-error message format, and `ploidy` is documented as a
  diploid-confirmation guard rather than a settable option.
* `gl2vcf()`: coordinate and fidelity fixes. (1) Explicitly supplied
  `snp.pos`/`snp.chr` fields now take precedence over populated
  `@position`/`@chromosome` slots -- previously the arguments were
  silently ignored whenever the slots were already set (including stale
  tag-offset copies in legacy objects), emitting within-tag offsets as
  genome POS; the source used is announced at `verbose >= 1`, and a
  malformed (wrong-length) slot is discarded with a warning instead of
  silently. POS values change for callers passing `snp.pos` on objects
  with a populated position slot. (2) REF is now pinned from `loc.all`
  via PLINK's `--a2-allele`, so loci fixed for the alternate allele in
  the exported sample keep their recorded reference base instead of
  degrading to `N`. (3) SilicoDArT objects now stop with a clear datatype
  error (`accept = "SNP"`) instead of failing cryptically inside
  gl2plink. (4) A factor-typed `snp.pos` field previously exported its
  factor level codes as POS; values are coerced via character and
  non-integer position fields are a fatal error. (5) The PLINK binary is
  checked up front with a download URL in the error. (6) `verbose = 0`
  is silent: gl2plink runs at the passed verbosity and PLINK's captured
  log prints only at `verbose >= 2`. (7) The ped/map intermediates are
  written to tempdir() and, with PLINK's .log/.nosex by-products,
  removed on success. Docs corrected accordingly.

* `gl2fasta()` review (driven by a user report of single-locus output
  on a merged/co-analysed dataset - reproduced): the internal
  overshoot pre-filter is no longer silent - removals are warned at
  verbose >= 1, and fewer than 2 surviving loci is a fatal error
  explaining the SnpPosition-vs-TrimmedSequence inconsistency;
  factor-coded SnpPosition is now recovered correctly (it previously
  fed substr with factor level codes, silently corrupting method-1
  sequences) and genuinely non-numeric positions are fatal; sink() is
  guarded with on.exit (a mid-write crash could hijack the session
  console); the documented method=0 help listing works (a cat():cat()
  typo crashed it); NULL-safe flag check; @return corrected; FASTA
  lines written without trailing spaces.
* `utils.dartR.class.def` (the dartR S4 class layer) review: the show
  method's loc.metrics detail no longer vanishes when ind.metrics is
  absent (a copy-paste guard); subsetting with an unmatched locus name
  fails informatively instead of the cryptic "Cannot subset a SNPbin
  with mixed subscripts"; negative indices now work on FBM-backed
  objects; the dead .fbmsub_copy helper (broken free variable,
  superseded by big_copy) removed. Verified sound: XOR validity,
  subset/rbind/cbind round-trips exact, glSum/glMean match adegenet
  exactly, and the FBM layer matches the gen-backed reference exactly.
  Governance notes recorded: the duplicate-class cache message needs
  ecosystem-level resolution (retired dartR monolith defines the same
  class); adegenet-internal getFromNamespace fragility; glSum/glMean
  shadow-function dispatch limitation.
* `utils.heatmap()` review: documented as a deliberate fork of
  gplots::heatmap.2 (colored dendrogram leaf labels via dendextend,
  auto-sized margins, NULL side-color defaults); matrices without
  dimnames no longer crash on the auto-margin; clustering verified
  identical to gplots::heatmap.2. Note for the custodian: gplots
  remains in Imports but is no longer used anywhere in the package.
* `utils.plink.run()` review: the composed command is now well
  formed - plink.path="path" performs a bare PATH lookup as
  documented (previously a literal "path/" prefix), and a space is
  guaranteed before --out (previously glued onto the last syntax
  token); marked internal (stays exported).
* `utils.collapse.matrix()` review: the within-population diagonal
  now averages distinct pairs only (self-distances of zero deflated
  it; matrix output only - dist consumers such as gl.dist.pop are
  unaffected); empty fatal-error messages restored; off-diagonal
  means verified exact (unchanged); marked internal (stays
  exported). Addendum from the population-distance chain review: the
  name guard is two-directional (a D computed on a subset of the
  object's individuals, or an unnamed matrix, now fails with a clear
  message instead of a bare subscript/dimnames error).
* `utils.transpose()` review: verified exact (dimension/name/metric
  swaps; double transpose reproduces the original genotypes);
  narration comments tidied; marked internal.
* `utils.stats` (std.error): documented; computation verified.
* `utils.plot.save()` review: the documented default verbose = NULL
  no longer crashes (verbosity is normalized on entry); a nonexistent
  save directory now falls back to tempdir() instead of crashing on a
  tempfile() path; the "No plot saved" note respects verbose 0; the
  unused ggsave passthrough claim dropped from the docs; marked
  internal (stays exported - called family-wide).
* `utils.flag.start()` review: verbosity contract verified (no
  behaviour change); docs completed; marked internal (stays
  exported - called family-wide).

* `utils.jackknife()` review: a unit vector of length > 1 now reaches
  the informative stop instead of crashing ("the condition has length
  > 1" - the check took length(unit == 1), the length of the
  comparison); the gl.set.verbosity save/restore inside each
  replicate no longer leaks 6 lines at verbose 0; the recal= partial
  match corrected to recalc=; marked internal (stays exported).

* `utils.hamming()` REMOVED: an orphan since utils.hamming.engine
  superseded it (no callers anywhere in the family), carrying an
  off-by-one (comparison started at the last recognition-site base)
  and a wrong proportion denominator. The verified engine in
  utils.hamming.blocks.r is the live implementation.

* `utils.dist.binary()` review: method="bray-curtis" is now accepted
  (it was documented and implemented but missing from the validation
  list, silently falling back to simple matching); the scale warning
  respects verbose 0; Jaccard and Sorensen verified exact
  (unchanged); doc leftovers ("N11") removed.

* `utils.dist.ind.snp()` review: the Simple and Absolute distances
  are now invariant to which allele is scored as reference, computing
  shared alleles per locus as documented (previously both-homozygous-
  reference pairs scored as sharing NO alleles, so relabelling the
  reference allele changed the distances) - Simple and Absolute
  values change; a merge artifact is cleaned up (an ungated progress
  message printed at verbose 0, Completed/Returning printed twice,
  and the result was converted dist->matrix->dist); Euclidean and
  Manhattan verified exact and unchanged.

* `utils.n.var.invariant()` review: the secondaries-in-history
  warning now respects verbose 0; variant/invariant counts verified
  exact on platypus.gl; docs tidied; marked internal.

* `utils.allelic.richness()` review: @return corrected (it documented
  "calling function name"); parameters documented; rarefaction kernel
  re-anchored (validated in PR #286); marked internal.

* `utils.hwe` helpers (GenerateSamples, CritSam, CritSam_Chi):
  documented with @noRd headers; enumeration and critical-sample
  outputs verified; unused variable removed.

* `utils.het.pop()` review: the returned vector now carries
  population names; documentation added (including the mean-n
  definition of the unbiased correction); computation verified
  exact.

* `utils.is.fixed()` review: documentation corrected to the numeric
  1/0/NA contract callers rely on (docs claimed TRUE/FALSE); truth
  table and tolerance boundaries verified; marked internal (stays
  exported for dartR.popgen).

* `utils.basic.stats()` review: loci absent from one or more
  populations no longer NaN-poison the cross-population statistics
  (the harmonic mean sample size is now taken over the populations
  that carry the locus, and single-population loci drop from the
  overall averages as undefined) - output now matches
  `hierfstat::basic.stats` exactly in all verified cases, honouring
  the documented equivalence claim; overall Fst/Fis change on
  datasets with per-population absent loci; the single-individual-pop
  error regains its message; marked internal (stays exported for
  dartR.popgen).

* `utils.reset.flags()` review: the `value` argument is now honoured
  when creating @other$verbose (it was validated then ignored -
  hardcoded 2); the bogus loc.metrics$monomorphs COLUMN is no longer
  invented (monomorphs is a flag, not a locus metric); the
  out-of-range value warning is gated at verbose >= 1; marked internal
  but stays exported (dartR.sim calls it).

* `utils.recalc.avgpic()` review: NULL-safe monomorphs-flag check (no
  crash on flag-less objects); marked internal; SNP and SilicoDArT
  arithmetic both verified exact (unchanged).

* `utils.recalc.maf()` review: SilicoDArT input is now rejected (the
  doc already promised SNP-only); NULL-safe monomorphs-flag check;
  marked internal; maf arithmetic verified exact (unchanged).

* `utils.recalc.freqhomref()` review: SilicoDArT input is now rejected
  (was silently given diploid metrics); NULL-safe monomorphs-flag
  check; marked internal; arithmetic verified exact (unchanged).

* `utils.hamming.engine()` review: verified sound - the block-hashing
  dedup detector matches a brute-force scan exactly and pairwise
  mismatch counts are exact. No changes; characterization tests
  added.

* `utils.impute` helpers: matrix2gen(parallel = TRUE) no longer
  crashes ("object 'i' not found" - a live path via gl.impute's
  parallel argument); the placeholder documentation block replaced
  with real @noRd headers and the ghost utils.impute Rd removed.

* `utils.recalc.freqhomsnp()` review: SilicoDArT input is now rejected
  (was silently given diploid metrics); NULL-safe monomorphs-flag
  check; marked internal; arithmetic verified exact (unchanged).

* `utils.recalc.freqhets()` review: SilicoDArT input is now rejected
  (previously presences were silently counted as heterozygotes);
  NULL-safe monomorphs-flag check (no crash on flag-less objects);
  marked internal; FreqHets arithmetic verified exact (unchanged).

* `gl.filter.pa()` review: SilicoDArT frequencies are no longer halved
  as if diploid, so presence-fixed private alleles are found (on
  testset.gs the kept set grows from 45 to 66 loci); a bogus population
  name now fails informatively; the filtered object returns invisibly
  with a before/after summary at verbose >= 2.

* `utils.check.datatype()`: four fixes to the package's central datatype
  dispatcher. (1) The all-NA screen ran a full `gl.filter.allna()` pass on
  every function entry at `verbose >= 2` (~50 ms per call on the small
  test dataset); replaced with a direct single-pass check -- same
  warnings, strictly less work -- and the all-NA-individuals warning no
  longer carries the copy-pasted loci wording. (2) Mixed or non-diploid
  ploidy was silently classified as SNP with no notice (the
  misspecification error was unreachable); it is still classified as SNP
  (the polyploid paths depend on this) but now announces itself at
  `verbose >= 2`. (3) `accept = "genlight"` (or "dartR") alone rejected
  every genlight object because the returned datatype is never literally
  "genlight"; those entries now admit both genotype datatypes, unless a
  specific datatype is also listed, in which case the specific listing
  governs (so `c("genlight","SNP")` remains SNP-only). (4) The
  unknown-class warning printed at `verbose = 0` (gated at >= 1); dead
  code removed and the @return documentation aligned with the strings
  actually returned ("matrix", "glPca").

* `gl.filter.maf()`: six fixes; the global-path arithmetic was verified
  against an independent recomputation and is unchanged. (1) The by.pop
  path collapsed when NO loci qualified for removal -- `x[, -integer(0)]`
  subset to zero loci and crashed (e.g. whenever no population met
  ind.limit); the object is now returned unchanged with gated messages
  for both the no-qualifying-loci and no-qualifying-populations cases.
  (2) `plot.file` with `plot.display = FALSE` crashed ("object 'p3' not
  found"), and the by.pop display path with 0-1 qualifying populations
  referenced a plot built only in the global path; each path now prints
  the plots it builds and the save is guarded with a gated notice. (3)
  The per-population plot list was named positionally with the
  qualifying-population names, mislabeling plots whenever a
  non-qualifying population preceded a qualifying one; plots are now
  named by their source population. (4) The threshold-range warning
  printed at `verbose = 0` (gated), and the MAC interpretation applied at
  `threshold >= 1` while the documentation says "> 1" (aligned). (5) The
  by.pop tally used `maf <= threshold` against the documented "less than"
  (aligned to `<`); loci with undefined MAF removed by the global filter
  are itemised at `verbose >= 3` and the NA policy is documented; the
  return is invisible per the filter convention; dead assignments and an
  irrelevant monomorphs preamble removed. (6) SilicoDArT objects were
  accepted and run through the MAF machinery producing meaningless
  values; the function is now restricted to SNP data (fatal on
  presence/absence input), with the matched report to follow.

* `utils.recalc.callrate()` review: marked internal
  (@keywords internal, per the utility-function policy); seealso
  corrected; CallRate arithmetic verified exact (unchanged).

* `gl.report.pa()` review: the Chao1/Chao2 estimates of undetected
  private alleles are now computed from the pair being compared
  (previously they were based on minor allele frequencies recomputed
  across ALL populations in the dataset, pairs with zero private
  alleles silently fell back to the full locus set, and f1/f2 were
  taken by table position rather than by category) - Chao values
  change for most datasets; `plot.file` without `plot.display` no
  longer crashes; a mistyped `method` and a datatype-mismatched `x2`
  now fail informatively; `method="one2rest"` with the default palette
  no longer crashes and its orientation no longer depends on
  population names sorting before the internal sentinel; missing
  suggested packages are now fatal errors; `test.asym` documentation
  now describes the implemented permutation test and the test is
  skipped with a warning for SilicoDArT; verbose=0 is now silent.

* `gl.filter.overshoot()`: loci whose overshoot status cannot be
  assessed (missing TrimmedSequence or SnpPosition) previously passed
  through the filter silently; by direction of the module coordinator
  they are RETAINED (not removed), and their count is reported in the
  `verbose >= 3` summary so the user knows they were not assessed. The
  no-op message is gated at `verbose >= 1`, the return is invisible, and
  the verbose >= 3 locus listing no longer carries a stray trailing
  comma. Core removal logic verified exact against independent
  recomputation (testset.gl: 21 removed, 234 retained).
* `gl.report.diversity()` review: the Shannon (q=1) measures no longer
  count per-population all-missing loci as zero diversity, and the
  internal vector misalignment that corrupted one_H_beta via logical
  recycling is fixed; pairwise one_H_beta and two_H_beta are now pooled
  from the pair being compared with the pairwise correction factor
  (previously they pooled ALL populations, so a pair's beta changed
  when unrelated populations were present - zero_H_beta already
  behaved correctly) - q=1 values and the q=1/q=2 beta matrices change;
  tables now print at verbose >= 1 only and an invalid `table` value is
  a fatal error instead of silently suppressing output; `plot.file`
  without `plot.display` no longer crashes; SilicoDArT input is now
  rejected (the entropy formulas assume diploid SNP genotypes;
  previously it computed meaningless indices silently); no gl.colors
  banner at verbose 0.

* `gl.report.maf()`: five fixes; the MAF values themselves are unchanged
  (verified against a hand computation). (1) The verbosity contract was
  broken wholesale -- a hardcoded `verbose = 3` inside the per-population
  function overrode the caller's setting, and the overall statistics, the
  quantile table, the limit-coercion warnings and the singleton-drop
  notice were all ungated, printing 343 lines at `verbose = 0` on the
  test dataset. Per-population statistics now display at `verbose >= 3`,
  the overall statistics and quantile table at `verbose >= 1`, warnings
  at `verbose >= 1` and the singleton notice at `verbose >= 2`. (2) A
  second handwritten FLAG SCRIPT START block duplicated the "Starting"
  banner and referenced a `build` variable that resolved only by lexical
  accident to a stale package-level constant -- removed. (3) The
  singleton-drop notice reported population indices instead of names.
  (4) The overall statistics printed the unfiltered locus count while the
  statistics were computed on monomorph-filtered data (both counts now
  shown); the bad-as.pop error directed users to loc.metrics though the
  check is against ind.metrics. (5) Docs: @return corrected (the function
  returns an invisible dataframe of MAF by locus and population, not "an
  unaltered genlight object"); as.pop wording; "3r quantile" typo (x2);
  verbose canon. (Amendment, with gl.filter.maf F6:) restricted to SNP
  data -- MAF is undefined for presence/absence data, which was
  previously accepted and produced meaningless values.

* `gl2eigenstrat()`: the genotype file counted the wrong allele -- the raw
  dartR score (copies of the ALTERNATE allele, the sixth column of the .snp
  file) was written where EIGENSTRAT defines copies of the REFERENCE allele
  (the fifth column), so every genotype was allele-flipped relative to the
  declared ref/var and downstream statistics keyed to allele identity were
  inverted; the geno value is now `2 - score` (9 stays missing). Numerical
  output therefore changes for every exported dataset. Factor metadata
  fields nominated via `snp.chr`/`snp.pos` were coerced with `as.numeric()`
  directly, writing factor LEVEL CODES instead of the actual values; the
  coercion now goes through `as.character()`, chromosome labels 'X', 'Y',
  'MT'/'mtDNA' and 'XY' are mapped to the documented 23/24/90/91 encoding,
  and the documented (but never implemented) removal of loci with illegal
  chromosome values now happens, with a warning at `verbose >= 1`; if no
  locus is encodable the function stops (the roxygen example no longer
  nominates the un-encodable platypus chromosome field). SilicoDArT data
  (which produced a malformed 4-column .snp file) is now rejected.
  `sex.code`/`phen.value` lengths are validated instead of being silently
  recycled down the .ind file. Docs aligned with the signature (`pos.cM`
  default is 0; the numeric sentinel semantics of `snp.pos`/`snp.chr` are
  stated).
- gl.amova has moved to dartR.popgen (green-striped-gecko/dartR.popgen#88), carrying the PR #378 review fixes and its test; no dartRverse code calls it.

* `utils.read.fasta()` (the engine behind `gl.read.fasta()`): the
  genotype-classification core was redesigned. Previously anything that
  was not hom-ref or hom-alt fell through to heterozygous, so missing
  data (N, gaps, V/H/D/B) was coded 1 instead of NA, silently inflating
  heterozygosity; truly triallelic columns of homozygote classes escaped
  the more-than-2-alleles skip and their third allele came back as a fake
  het; lowercase (softmasked) bases registered as distinct alleles and
  fabricated loci; and the documented most-frequent-allele reference
  actually followed the modal GENOTYPE class, flipping polarity when the
  het was modal. Genotypes are now classified explicitly against the
  allele pool (NA for unrecognized codes), multiallelism is detected from
  the pool, sequences are upper-cased on read, and ref/alt follow summed
  allele counts. Note the output changes: NA at masked/missing sites,
  triallelic columns dropped, polarity corrected at het-modal loci.
  Also: unequal-length sequences now stop with an error naming the
  offending records (base-R recycling previously fabricated loci from
  ragged input); the self-referential signature defaults
  (`parallel = parallel`, `verbose = verbose`) are replaced with working
  ones; `parallel = TRUE` works (n.cores = NULL resolves to all cores;
  Windows falls back to serial with a gated warning); the
  multiallelic-skip and no-polymorphism messages are gated at
  `verbose >= 1` (they printed even at verbose 0); and `merge_gl_fasta()`
  refuses duplicate individual names per file, which `merge()` previously
  joined many-to-one, silently copying one record's genotypes onto
  several rows.
* `utils.read.dart()` (affects `gl.read.dart()` imports): five fixes.
  (1) A DArT report in which any locus ID had an unexpected row count
  (e.g. one allele row hand-deleted or duplicated) previously had EVERY
  locus removed by the duplicate-ID cleanup, and the read died with a
  misleading "must be either 1row or 2row" error; now only loci whose
  row count departs from the modal count are removed, and the removal
  is reported accurately. (2) A genuine 1-row report with no
  heterozygous calls was silently misread as 2-row, halving the locus
  count and pairing unrelated rows into fabricated genotypes; the
  genotype-range and locus-ID-count signals are now cross-checked and
  the read stops with a clear error when they disagree. (3) The
  service/plate extraction read past the header block when the file had
  fewer header rows than plate.row + 2, pasting the column-header and
  first data rows into ind.metrics (the canonical
  testset_SNPs_2Row.csv, with 3 header rows, got plate_location
  "UC_1-AA0109150"); out-of-range rows now yield NA with a gated
  warning. (4) The duplicate-individual-name uniquification was
  computed but never applied - the promised '_n' suffixes now actually
  reach the genotype column names (previously read.csv's '.1' suffixes
  leaked through, contradicting the printed warning). (5) All warnings
  are now verbosity-gated: `verbose = 0` is silent.
* `gl.read.fasta()`: the closing fbm block ignored the `fbm` argument --
  an empty `if (fbm) {}` guard meant `gl.gen2fbm()` always ran (clearing
  `@gen`), and at `verbose <= 2` a misattached else-branch wiped `@fbm`
  as well, so every default-verbosity call (including the documentation
  example) returned an object with no genotype data at all (nInd = 0),
  and `verbose >= 3` calls returned an FBM-backed object although fbm
  was declined. The conversion now runs only when `fbm = TRUE`. Note
  the output change: at default verbosity the function now returns the
  actual data. Also: input where no file yields SNPs (no polymorphism)
  now stops with a clear message instead of an opaque "argument of
  length 0" failure, and line-wrapped (multi-line) FASTA -- which the
  two-line reader silently mis-groups -- is rejected up front with an
  informative error, now documented.
* `utils.read.ped()` (vendored copy of `snpStats::read.pedfile`; used by
  `gl.report.ld.map()`): four fixes. (1) `lex.order = TRUE` swapped the
  map alleles but discarded the switched genotype matrix, so dosages at
  reordered loci counted the wrong allele relative to the returned map;
  the switched matrix is now assigned. (2) `show_warnings = FALSE` also
  disabled the multi-allelic NA reset, silently retaining genotypes at
  loci the function itself had classified as unreliable; the reset now
  runs regardless and only the warnings are gated. (3) Multi-allelic
  detection missed a novel allele whenever it was paired with a known
  one (e.g. an A/T carrier at an A/G locus was coded homozygous A/A,
  unflagged); detection is now per allele column. (4) The `split`
  argument was honoured when counting loci but ignored on data lines
  (a comma-separated file returned mangled ids and all-NA genotypes
  without error); it is now applied to every line. Real documentation
  replaces the placeholder header (the previous `@return` wrongly
  claimed a genlight object; the function returns a
  list(genotypes, fam, map)) and the helper is now exported with the
  internal-use warning, per the utils convention. F1/F3/F4 reproduce in
  the installed snpStats original and are candidates to offer upstream.
* `gl.dist.pop()`: distances landed on the wrong population labels
  whenever the population factor levels were not in alphabetical order
  -- the frequency matrix rows come back from `reshape2::dcast` in
  alphabetical order while the matrix labels were assigned positionally
  from `popNames(x)`. The rows are now re-anchored to `popNames(x)` by
  name, so euclidean, nei, reynolds and chord distances follow their
  labels for any level order (numerical output changes -- becomes
  correct -- for objects with non-alphabetical population levels;
  alphabetical-level objects, including all the packaged datasets, are
  unchanged, verified to the anchor values). Also: `plot.file` with
  `plot.display = FALSE` (including any `verbose = 0` call) no longer
  crashes with "object 'p3' not found" after computing the distances --
  the plot is built whenever it is displayed or saved; an unknown
  `method` now stops (previously a non-fatal "Fatal Error" printed at
  every verbosity and euclidean ran silently); `type = "matrix"` returns
  a full symmetric zero-diagonal matrix for every method (previously the
  SNP frequency methods returned only the lower triangle, upper triangle
  and diagonal NA); the missing-reshape2 guard stops instead of
  returning -1; fewer than two populations fails fast with a clear
  message; documentation corrected (the scaled SNP euclidean maximum is
  0.5, not 1; the SilicoDArT method list completed; the reynolds
  -log(1-D) linearised variant and the chord 2*sqrt(2)/pi scaling
  constant stated; the `as.pop` error message points at ind.metrics).
  Note: `method = "sorensen"` routes through `gl.dist.ind()`, whose
  companion fix (adding sorensen to its accepted methods) delivers true
  Sorensen distances here once merged.
* `gl.recalc.metrics()`: nine fixes from its function review. The values of
  every metric the function computes are unchanged. (1) An object with no
  `loc.metrics` data frame was corrupted rather than repaired: the helpers
  read the slot with `$`, which partial-matches to `loc.metrics.flags`, so
  with more than one locus the call died with "replacement has N rows, data
  has 1" and with exactly one locus it silently returned the flags table as
  the locus metrics, losing AlleleID, TrimmedSequence, rdepth and the read
  counts without a message. A conforming table is now built, with a gated
  message; a table whose row count does not match the number of loci is now
  a fatal error naming the condition. (2) History: `gl.recalc.metrics()` is
  an implementation step of 26 functions in this package and two in
  dartR.popgen, and appended an entry on each of their behalf, so nested
  calls multiplied entries (two successive `gl.compliance.check()` calls on
  `testset.gl` went 1 -> 3 -> 5). A call made from inside another dartRverse
  function now appends nothing, and `mono.rm = TRUE` no longer adds
  `gl.filter.monomorphs()`' entry as well as its own. A direct call still
  appends exactly one entry recording the call. **This changes the history
  contents of every function that calls `gl.recalc.metrics()` internally.**
  (3) The `monomorphs` flag was passed through untouched when
  `mono.rm = FALSE`, so a stale TRUE suppressed the monomorph warnings of
  every downstream report; the flag is now set from a check made on the
  metrics just recalculated. (4) `mono.rm = TRUE` on data in which every
  locus is monomorphic crashed with "Subsetting resulted in zero loci"; the
  case now completes with a gated warning and the flag set FALSE. (5)
  `mono.rm` is validated. (6) The six `utils.recalc.*` helpers are called
  with `verbose = 0` and this function reports on their behalf, so
  `verbose = 1` prints two lines rather than fourteen and a warning is not
  repeated once per helper. (7) `@family environment` was indented, so
  roxygen read it as part of the title; the function now appears in its
  family index. (8) Documentation corrected: `maf` listed as recalculated,
  the false "only RepAvg and TrimmedSeq are unaltered" claim replaced with
  the accurate split, `@param x` widened to SilicoDArT, the standard
  `verbose` wording adopted, an `Author(s):` line and `@details` added. (9)
  The outdated `build =` argument dropped from `utils.flag.start()`. The
  review's F1 (`rdepth` and silico `AvgReadDepth` are not recalculated and
  can contradict the refreshed metrics) is deferred by the custodian and
  remains open, with its evidence in
  `function-review/reports/dartR.base/gl.recalc.metrics.md`.
* `gl.report.ld.map()`: `plot.file` now works with `plot.display = FALSE`
  (previously crashed with "object 'p4' not found" after the full
  computation); fully silent at `verbose = 0` (warnings gated, plot
  suppressed); loci sharing a map position that are excluded from the
  analysis are now reported at `verbose >= 1`; SilicoDArT objects are
  rejected at the entry datatype check; documentation corrected
  (`ld.max.pairwise` semantics, `plot.display` text, `ind.limit`
  boundary, family "matched report") and the truncation to pairs with a
  positive statistic value is now documented. The returned data frame is
  unchanged.
* `gl.filter.ld()`: `ld.report` is now validated at entry (a clear error
  names `gl.report.ld.map` instead of an obscure downstream failure);
  objects without `loc.metrics.flags` no longer crash with "argument is of
  length zero"; the "No pair of loci" message gates at `verbose >= 1`; a
  single history entry is appended per call (previously two, one exposing
  internal variable names); the `pop.limit` default is computed explicitly
  as half of the populations represented in `ld.report` (same value as
  before, previously an accident of lazy evaluation) and documented as
  such; documentation corrected (threshold boundary "at or above", actual
  sequential pair-resolution rule described). Loci removed are unchanged.
* `gl.report.ld()`: crash-restart fixed — rerunning with the same
  `chunkname` at `verbose < 2` used to crash with "subscript out of
  bounds" (the cached-return sat inside a verbosity guard), and chunk
  discovery searched the working directory while chunks were saved to
  `outpath`, silently defeating the restart for the default
  `outpath = tempdir()`. Chunk files are now written only when
  `save = TRUE` (previously always, with collision-prone
  `LD_chunks__i.rdata` names). Dependency guards now `stop()` instead of
  printing and returning -1 (and the redundant guards for Imports
  packages are gone); the internal `gl2gi` conversion is silenced at
  `verbose = 0`; SilicoDArT is rejected (`accept = "SNP"`). Documentation
  now describes what the function computes — LD across all loci with all
  individuals POOLED, not per population (use `gl.report.ld.map` for
  within-population LD) — and drops the incorrect claim that genind
  input is accepted. The returned statistics are unchanged.
* `gl.filter.factorloadings()`: `retain = TRUE` results change — it
  previously returned exactly the same object as `retain = FALSE` (the
  high-loading loci were dropped instead of kept, i.e. the complement of
  what was documented); it now retains the loci with |loading| at or above
  the threshold, as documented. A pca that does not match the genlight
  (different locus count after monomorph removal) now raises a clear error
  instead of silently recycling loadings onto the wrong loci. A single
  history entry naming the function replaces the two internal entries;
  glPca and axis validation errors are labelled; the documented `...`
  (save parameters) is now forwarded; documentation corrected (`@return`,
  family, verbose text).

* `gl.tree.fitch()`: bootstrap support values change -- the previous
  values were extracted from the wrong edges of the wrong tree (an
  edge-order assumption that ape does not honour) and typically displayed
  full support regardless of the data; supports are now computed by clade
  matching (`ape::prop.clades`) against the replicate trees, returned on
  the tree as `node.label` together with the majority-rule consensus tree
  (attribute `consensus.tree`), and the replicate distances are computed
  under the same substitution-model settings as the main tree (new
  `subst.model`, `pairwise.missing`, `min.tag.len` parameters). Also: the
  function now works on Windows (the PHYLIP invocation performed no stdin
  redirection, so no run could complete); the working directory is
  restored after bootstrap runs; `out.path` actually receives the PHYLIP
  output files; a supplied `outgroup` roots the returned tree; the
  distance matrix is scaled to preserve branch-length precision through
  PHYLIP's 6-decimal tree file; guards added for fewer than 4 taxa,
  duplicate 10-character truncated labels and D/x population mismatch;
  `verbose = 0` is now fully silent (previously the whole PHYLIP menu
  dialogue printed); plotting is optional (`plot.display`) and decoupled
  from the returned result.
* `gl.set.verbosity()`: invalid values were a silent no-op that claimed
  success -- `gl.set.verbosity(7)` warned, left the global option
  untouched, and then printed "Global verbosity set to: 7"; a character
  value rode string comparisons through the same false echo; and
  `gl.set.verbosity(NULL)` crashed with "argument is of length zero".
  Invalid values (including NULL) now warn and coerce to the default 2,
  which is then genuinely set and honestly reported. The function also
  returns the value actually set (invisibly), honouring its @return
  contract (it returned NULL), validation runs before the start banner,
  and the header gains its @family tag.

* `gl.impute()`: **`method = "frequency"` changes numerical output.** It now
  does what its documentation always claimed -- a deterministic fill with
  the expected dosage 2q at the locus in the individual's population,
  rounded to the nearest valid genotype (exact ties, 2q = 0.5 or 1.5, go to
  the heterozygote). It was previously a random Bernoulli(q)-pair draw,
  distributionally identical to `method = "HW"`; every existing
  "frequency" result therefore changes, and repeated runs now give
  identical results without a seed. For presence-absence data the fill is
  the majority band state (ties to presence).
* `gl.impute()`: residual missing values on FBM-backed objects are now
  preserved as NA in the imputed object; the raw FBM write-back previously
  coerced them to genotype 0, fabricating homozygous-reference calls at
  all-NA loci. Dense and FBM-backed runs now return identical genotypes.
* `gl.impute()`: presence-absence (SilicoDArT) support is now real:
  "random" draws from 0/1, "HW" becomes a Bernoulli draw with the band
  frequency, and the residual fill draws Bernoulli from the global band
  frequency ("random" and "frequency" previously wrote genotype 2s into
  0/1 data, and "frequency"/"HW" then crashed with "negative
  probability"); "beagle" is blocked for presence-absence data with a
  clear error. Also: `method` is validated up front; all-missing-locus
  warnings are ploidy-aware and no longer fire for any locus above 50%
  missing; the random/beagle branches report every affected population
  rather than only the last; the beagle path checks the jar, java and
  PLINK up front with informative errors (instead of returning -1 or a
  raw system error) and no longer blanks singleton-scaffold chromosome
  names in the returned object.
* `gl.dist.ind()`: `method = "sorensen"` now computes the Sorensen (=
  Dice) distance -- previously "sorensen" was missing from the
  accepted-methods list, so it was silently coerced to simple matching
  with a warning that leaked at `verbose = 0` (numerical output changes
  for sorensen callers, including `gl.dist.pop(method = "sorensen")`,
  which routes through here). Also: the unknown-method fallback warnings
  are gated at `verbose >= 1`; `type` is normalised with `tolower()` and
  validated ("Matrix" now returns a matrix; an unrecognised type stops
  instead of silently returning a dist); NA distances (e.g. from an
  individual with no scored genotypes in common with another) are
  counted and warned at `verbose >= 1` instead of propagating silently;
  documentation corrected (`scale` default is FALSE and applies to
  euclidean only; the Sorensen/Bray-Curtis synonymy for binary data is
  stated; the `@author` line repaired).
* `gl.read.PLINK()`: the returned object now contains the genotypes at
  every verbosity. Previously `gl.gen2fbm()` always ran and the result
  was discarded at verbose <= 2 (the default), so default-settings calls
  returned an object with no genotype data; at verbose = 3 the object
  came back FBM-backed even with `fbm = FALSE`. The `fbm` argument now
  decides the backend, as documented. Individual metafile rows are now
  matched to the .fam individuals by `id` (previously bound in file
  order, silently misassigning metadata when the orders differed), and
  a `pop` column in the metafile is now applied to `pop(gl)`. A missing
  AlleleID column in the locus metafile now stops with the intended
  message instead of an unrelated subscript error. The .ped-to-.bed
  conversion now runs in a temporary directory (previously it wrote
  .bed/.bim/.fam/.log into the user's input directory) and PLINK's
  console chatter is suppressed below verbose 3. Documentation: the
  dosage orientation is now stated (counts allele.2, the PLINK 1.x
  major allele -- the opposite orientation to `gl.read.vcf()`, which
  counts ALT).
* R CMD check: silenced "no visible binding" NOTEs for ggplot aes
  variables in `gl.report.hamming()` (Threshold, Removed, current) and
  `gl.report.secondaries()` (count).
* CI repair: six review-campaign test files hardcoded expectations from a
  `testset.gl` that contained 3 all-NA loci; the CRAN dartR.data 1.2.2
  `testset.gl` (the one CI installs) has none, so every `dev` run since
  PR #252 failed. Expectations in test-gl.filter.allna,
  test-gl.fixed.diff, test-gl.report.allelerich, test-gl.report.basics,
  test-gl.report.hwe and test-gl.filter.hwe recomputed against the CRAN
  data.
* `gl.report.hwe()` / `gl.filter.hwe()`: crashed ("Subsetting resulted
  in zero loci") whenever a population was entirely monomorphic --
  `gl.filter.monomorphs()` cannot return a zero-locus object, so the
  functions' own skip-empty-populations logic was unreachable. Such
  populations are now caught and skipped as intended.
* `gl.report.allelerich()`: the bootstrap path crashed ("replacement has
  length zero") when `boot.ci()` returned no interval for a statistic
  constant across replicates; the CI cells are now left as NA.
* `gl.report.heterozygosity()`: three crash fixes and conformance, applied
  after consultation with the authors. (1) `subsample.pop = TRUE` crashed
  whenever ANY population was below `n.limit` (utils.subsample.pop stored
  NA placeholders that rbindlist rejects) although n.limit is documented
  as a skip threshold; small populations are now skipped and the
  subsample plot keys its colours by population name so skipped
  populations cannot desync the palette. (2) `method = 'ind'` with
  `plot.display = FALSE` and `verbose >= 3` crashed ("object 'outliers'
  not found") -- outliers are now computed from the data (Tukey boxplot
  statistics) independent of plotting. (3) `subsample.pop = TRUE` with
  `method = 'ind'` crashed ("object 'res_sub' not found"); it is now
  ignored with a gated warning. Also: the method-coercion,
  negative-n.invariant and secondaries-history warnings printed at
  `verbose = 0` (now gated); the nboots/CI check raises a proper error
  condition; the subsample return is a named list
  (`$subsample`, `$results`) and documented; `plot.file` with
  `plot.display = FALSE` no longer crashes on the unbuilt plot (gated
  warning, nothing saved); header conformance (matched report family,
  author credits, preliminaries order). Reported statistics are
  unchanged.
* `gl.filter.heterozygosity()`: five fixes. (1) Removing individuals
  left the locus-metric flags stale (CallRate, allele frequencies, PIC
  values and the like were invalidated by the removal but still flagged
  valid); the flags are now reset when individuals are removed, at
  every verbosity, matching gl.drop.ind. (2) An individual with all
  genotypes missing (undefined heterozygosity) crashed the filter; such
  individuals are now removed and itemised in the `verbose >= 3`
  summary. (3) The monomorphs warning printed at `verbose = 0` on every
  call whose monomorphs flag was FALSE, and the unguarded flag access
  crashed on flag-less objects; now isFALSE()-guarded and gated at
  `verbose >= 2`. (4) Reversed thresholds (t.lower > t.upper) fell
  through to a cryptic zero-individuals error and the t.lower range
  message named t.upper; thresholds are now swapped with a gated
  warning and the message corrected. (5) The return is invisible, and
  the indented `@family` tag no longer leaks into the rendered help
  title.
* `gl.load()`: the `compliance` parameter (documented [default FALSE]) was
  inert -- `gl.compliance.check()` ran unconditionally on every load,
  potentially modifying the loaded object and printing its full output
  even at `verbose = 0`. As directed by the module coordinator, the check
  now runs only when `compliance = TRUE` (with verbose passed through).
  The "Loaded object" and fbm-conversion messages are gated at
  `verbose >= 2` (they printed at verbose 0). A missing file and an RDS
  that does not contain a genlight object now give clear fatal errors
  (previously a raw connection error and a cryptic "no applicable method
  for `@`" failure). Docs corrected: `fbm` default is FALSE (was
  documented TRUE), `file` is the file to read (was "receive data"), an
  examples block added.
* `gl.filter.replicates()`: four fixes. (1) When the members of a
  replicate pair had tied missing-data rates -- the exact-duplicate case
  -- BOTH were removed (the doubled pair table evaluated the drop rule in
  both orientations); pairs are now canonicalised and deduplicated before
  the rule is applied, so one member per pair is removed (ties remove the
  alphabetically later individual), and this works with report tables
  generated both before and after the companion gl.report.replicates
  fix. (2) Re-thresholding to an empty set crashed via gl.drop.ind ("no
  individuals to drop"); the object is now returned unchanged with a
  gated message. (3) `replicates.report` was never validated -- the
  report's old no-pairs string return crashed with "$ operator is invalid
  for atomic vectors"; malformed input is now a clear fatal error, and an
  empty table is handled as nothing-to-drop. (4) The history of the
  returned object recorded gl.drop.ind's call instead of
  gl.filter.replicates' (gl.drop.ind now runs quietly, the dropped
  individuals are itemised at `verbose >= 2`, and a single
  gl.filter.replicates entry is appended); a datatype check was added and
  the object is returned invisibly per house standard.
* `gl.report.replicates()`: four fixes. (1) The pair table carried BOTH
  orderings of every pair, and the drop rule picked the opposite member in
  each ordering whenever missing-data rates tied -- which is precisely the
  exact-duplicate case -- so BOTH replicates landed in `ind.list.drop` and
  the histograms double-counted every pair. The table now holds one row
  per unordered pair (ind1 = the individual earlier in the object) and a
  tie deterministically drops ind2. (2) When no pairs passed the
  thresholds the function returned a bare character message instead of the
  documented 3-element list, crashing `gl.filter.replicates()` downstream;
  it now returns the documented structure with an empty table and a gated
  message. (3) `ind.list.rep` used `>= perc_geno` where every other output
  used `>` (aligned), and its NaN diagonal previously injected NA entries;
  a datatype check was added, Rcpp/RcppParallel are checked before
  compiling (RcppParallel added to Suggests), and the plot no longer
  renders at `verbose = 0`. (4) The result is returned invisibly (house
  standard) with the pair table printed at `verbose >= 3`; docs cleaned
  (@family, literal "##" headings, typos).
* `gl.report.hwe()` / `gl.filter.hwe()`: the functionality of the
  deprecated `gl.report.excess.het()` / `gl.filter.excess.het()` has been
  migrated into this pair. Both functions gain `direction` ('both',
  'excess', 'deficit') to restrict attention to heterozygote excess or
  deficit (an expected-heterozygote column, Het.exp, is now included in
  the report output), and `min.hobs` to restrict testing to loci with
  observed heterozygosity at or above a threshold, applied before the
  multiple-comparison adjustment. The published Robledo-Ruiz et al.
  (2023) excess-heterozygosity workflow is reproduced with
  direction='excess', min.hobs=0.5, ChiSquare test and fdr adjustment;
  the deprecated functions remain as thin wrappers that emit a
  deprecation warning and produce the same flagged/removed loci as
  before (verified exactly on the LBP dataset), and will be removed in a
  future release. Migrating also retires a defect in the old
  gl.filter.excess.het, which computed its per-population genotype
  counts on the wrong individuals (a recycled per-population index).

  Review fixes applied to the pair in the same change: (1) the skipping
  of monomorphic and small-sample populations happened only at
  `verbose >= 2`, so the tested set -- and any multiple-comparison
  pool -- depended on verbosity (verbose 0 tested all 30 testset.gl
  populations, verbose 2 only 23); populations are now skipped at every
  verbosity. (2) A missing HardyWeinberg package now raises a fatal
  error instead of returning -1, and ggtern is required only when
  ternary plots are requested (plot.out=TRUE). (3) The out-of-range
  alpha warning is gated at `verbose >= 1` and no longer calls alpha an
  "integer". (4) gl.report.hwe now carries a @family tag, and
  gl.filter.hwe's duplicated/indented family tags are repaired.
* `gl.reassign.pop()`: the `as.pop` metric name was never validated -- a
  name absent from `ind.metrics` assigned NULL to `pop(x)`, silently
  destroying every population assignment. A missing ind.metrics slot or an
  unknown metric name is now a fatal error naming the available metrics.
  A gated warning (verbose >= 2) reports how many individuals carry NA
  assignments when the chosen metric has missing values.
* `gl.reassign.ind()`: conformance fixes -- the five fatal exits now use
  the house `stop(error(...))` styling; the empty-selection notice is a
  gated `cat(warn())` (verbose >= 2) instead of an R `warning()`; repeated
  indices in a numeric `ind.list` are deduplicated.
* `gl.define.pop()`: the not-present-individual warning printed at
  `verbose = 0`; now gated at `verbose >= 2`. An irrelevant preamble that
  ran a full `gl.filter.monomorphs()` on every call solely to warn that
  monomorphic loci exist was removed. Style: the assignment message is now
  styled and printed after the assignment actually happens; dead
  `is.na(length())` condition removed; header tag order.
* `gl.merge.pop()`: two behaviour fixes. (1) Validation of `old` sat inside
  the `verbose >= 1` announcement block, so `old = character(0)` was fatal
  at `verbose >= 1` but a silent no-op at `verbose = 0`; validation now
  runs upfront at every verbosity. (2) Populations listed in `old` that do
  not exist in the dataset were silently ignored -- a mistyped population
  name left the object unchanged with no message; this is now a fatal
  error naming the missing populations (matching `gl.rename.pop`). Tidy:
  redundant genlight check and duplicate validation removed; the
  description opening (a copy-paste about csv metadata files) corrected;
  header tag order.
* `gl.save()`: the "Saved object" / "Load again" messages (and the
  FBM-conversion message) printed at `verbose = 0`; now gated at
  `verbose >= 2`, and the message no longer calls the RDS file an "RDA
  file". The @return contract ("the input object") is now honoured: the
  class-attribute stamping and any FBM-to-gen conversion apply only to
  the copy that is saved, and the input is returned unchanged. A
  nonexistent target directory gives a clear fatal error instead of a raw
  connection error; description wording corrected.
* `gl.report.allelerich()`: plumbing fixes; the rarefaction calculation
  itself was verified against an independent recomputation and is
  unchanged. (1) The plot rendered at `verbose = 0` (missing
  `plot.display` guard) and the lazy signature default
  `gl.colors("dis")` printed a 3-line banner at every verbosity -- both
  silenced (the default now passes `verbose = 0`). (2) An unrecognized
  `error.bar` value crashed downstream with "object 'max_val' not
  found"; unknown values now coerce to "SD" with a gated warning, and the
  silent override of the user's error-bar choice when `nboots > 0` is
  announced at `verbose >= 2`. (3) The package check called
  `requireNamespace()` on a vector, silently checking only dplyr, and
  returned -1 through a cat(); each package is now checked individually
  with `stop(error(...))`, and boot/Rcpp (needed for bootstrapping) are
  checked before the bootstrap path. (4) Dead code removed (an unused
  first plot built on every call, a commented parallel block, duplicated
  global declarations); `boot.method` validated; header canon.
* `gl.join()`: five fixes. (1) A join by shared loci LOST the individual
  metrics entirely -- plain `rbind()` returns NULL metadata and the
  function never rebuilt it; the combined object now carries the
  row-bound ind.metrics of both inputs (with the id column re-synced when
  duplicate names are made unique). (2) A join by shared individuals
  CRASHED ("replacement has 0 rows") on objects whose loc.metrics.flags
  data.frame lacks the OneRatio/PIC columns -- which is standard
  SNP-report data; only flags present in both objects are now combined.
  (3) The same flags block was triplicated (each path set the flags and a
  third copy ran again for both), warnings printed at `verbose = 0`
  (method deprecation, missing metrics/flags), and two
  `cat(error()) + stop()` splits returned no condition message; the
  duplicate block is gone, warnings are gated at `verbose >= 2`, and the
  exits use `stop(error(...))`. (4) SNP and SilicoDArT objects with
  matching names could be joined silently -- the datatypes were checked
  individually but never compared; now fatal. (4a, amendment) The
  historical legacy values `method='end2end'` and `method='sidebyside'`
  -- accepted by the pre-refactor implementation and described in the
  documentation ever since, and still used by real callers (the
  dartR.popgen gl.assign functions) -- were fatal because the legacy shim
  mapped only join.by.loc/join.by.ind; they are now mapped to their
  historical meanings, and any explicitly requested join is validated
  against the data so a mismatch fails with a clear message rather than
  a cryptic cbind/rbind error. (5) Messages used
  `substitute()` inside `cat()`, which printed garbage at `verbose = 2`
  and crashed at `verbose >= 3` whenever the arguments were expressions
  rather than names (e.g. `gl.join(x[1:7, ], y)`); arguments are now
  deparsed once for display. Docs: the description claimed the history was
  cleared (it is carried from the first object and appended); @details
  described method='sidebyside'/'end2end' values that never existed;
  duplicate/incorrect end-of-run summary removed; typos.
* `gl.sort()`: the history entry was appended as `c(match.call())`, which
  coerces the call to a list and corrupts the history chain -- now a
  proper call. A standard FLAG SCRIPT END block was added ("Completed:"
  never printed at any verbosity). The no-chromosome warning under
  `order.by.chr.pos` printed at `verbose = 0`; now gated at
  `verbose >= 2`, and the dartR-conversion notice gate aligned from
  `> 2` to `>= 2`. A redundant length re-check with a misleading message
  in the sort.by='ind' path was removed (the upfront validation already
  covers it); verbose param doc canon.
* `gl.select.shapes()`: four fixes. (1) The range validation was a
  parenthesis slip (`min(select < 0 | max(select > 25))`) that only fired
  when every element was negative -- a partially-negative `select` such as
  `c(-1, 5)` passed straight through to `pch`. It now correctly rejects any
  value outside 0-25. (2) The documented `x=` genlight argument was
  non-functional: `nPop(x)` was computed and discarded, so a `select` of
  the wrong length passed silently and a NULL `select` returned all 26
  shapes regardless of the number of populations. Now (matching
  `gl.select.colors`) a length mismatch with `nPop(x)` is a fatal error, a
  NULL `select` with `x` returns one shape per population, and more than
  26 populations without an explicit `select` is a fatal error (only 26
  distinct shapes exist). (3) New `plot.display` parameter (default TRUE);
  the palette chart was previously drawn unconditionally and can now be
  suppressed, and is suppressed automatically at `verbose = 0` (which
  previously also leaked the datatype banner). (4) Cosmetic: the genlight
  argument `x` was shadowed by the plot x-coordinates mid-function, a
  "Requires shapes" typo, and header conformance.
* `gl.colors()`: three fixes. (1) Both invalid-type exits used
  `cat(error(...))` followed by `stop(-1)`, so the error condition an
  upstream `tryCatch()` received carried the message "-1" while the real
  message printed to stdout where even `try(silent = TRUE)` could not
  suppress it; they now use `stop(error(...))`. (2) The return is now
  visible -- the documented example `gl.colors(2)` previously displayed
  nothing because the result was returned `invisible()`. (3) Documentation:
  the description listed a `"pal"` category that the function has never
  accepted (it was a fatal error) -- removed; the implemented but
  undocumented `"structure"` type (35 discrete colors, used by
  `dartR.popgen::gl.plot.snmf`) is now documented, along with header
  conformance fixes and an accurate `@return` (the four palette types
  return a function, not a vector). Note: `gl.colors()` evaluated as a
  default argument in a caller's signature still prints its banner at the
  caller's `verbose = 0` unless the default passes `verbose = 0`; changing
  the default verbosity was considered and deliberately not adopted.
* `gl.select.colors()`: five fixes. (1) An unrecognised `library` value
  silently returned base R's `colors()` FUNCTION as the "colour vector"
  (the internal variable was never assigned and lexical scoping found
  grDevices::colors); unknown libraries now coerce to the default
  scales/hue_pal with a warning at `verbose >= 1`. (2) Brewer requests
  are honoured: fewer than 3 colours are trimmed from the 3-colour pull
  (you get exactly what you asked for), and requests above the palette
  maximum return the maximum with a clear gated warning (previously 2
  requested delivered 3, and 12 requested from Blues silently delivered
  9). (3) Out-of-bounds `select` indices, which produced NA colours,
  are now a fatal error. (4) baseR palette='heat' passed validation but
  returned rainbow; it now dispatches to heat.colors. (5) The datatype
  banner no longer prints at `verbose = 0` when a genlight object is
  supplied, and missing packages raise a fatal error instead of
  returning -1. The internal default-colour contract used across the
  package (brewer, Blues, select=c(7,5)) is unchanged and
  regression-tested.
* `gl.report.basics()`: four fixes. (1) The function crashed on ALL
  SilicoDArT data ("length of 'dimnames' [2] not equal to array extent")
  because the composition table hard-coded four column names onto a
  three-class presence/absence table; the same crash occurred for SNP data
  in which any genotype class was entirely absent. The table is now
  tabulated over explicit per-datatype levels -- SilicoDArT is supported
  for the first time; SNP output values are unchanged. (2) The all-NA
  individuals listing printed the entire NA-padded array
  ("NA NA ind3 NA NA ..."); it now lists only the names. (3) Objects
  without an rdepth locus metric triggered a raw mean.default warning and
  printed NA for Average Read Depth; now "not available", and when
  present the mean is computed with na.rm and rounded. (4) Style and
  efficiency: per-locus for-loops replaced with vectorized counts, the
  two per-population gl.keep.pop subset loops merged into a single matrix
  pass, the datatype check passes verbose through, and the @return
  documentation states that NULL is returned invisibly.
* `gl.fixed.diff()`: four fixes. (1) The documented `mono.rm` parameter
  ([default TRUE], remove monomorphic loci before computation) was never
  referenced in the body, and the flag logic standing in its place was
  inverted -- `gl.filter.monomorphs()` ran only when the flag certified
  the data already monomorph-free, so monomorphic loci were effectively
  never removed. `$nloc` was inflated and `$pcfd` denominators with it
  (roughly fivefold understatement of percent fixed differences on a test
  subset with 208/255 monomorphic loci). `mono.rm` is now honoured:
  TRUE (default) removes monomorphs, FALSE retains them with the warning
  and reproduces the previous numbers. Raw `$fd` counts are unaffected
  (a monomorphic locus cannot carry a fixed difference), so
  `gl.fdsim()`/`gl.collapse()` amalgamation decisions are unchanged;
  `gl.dist.pop(method = "fixed-diff")` distances change with the
  corrected `$pcfd`. (2) The "false positives can only be simulated for
  tloc=0" warning printed at `verbose = 0`; now gated at `verbose >= 2`.
  A dead, unreachable `tloc.hold` block inside the pairwise loop (with
  its own ungated output) was removed. (3) The `verbose >= 4` return
  listing misnamed `$pval` as `$prob`, described `$sdfpos` with the
  copy-pasted expected-count text, and called `$gl` the input object;
  matrix dimnames and diagonals are now set once instead of on every
  pass of the pairwise loop. (4) Documentation: the details promised a
  per-comparison warning at sample sizes below 5 that never existed
  (actual behaviour: a global minimum-sample warning at n < 10,
  `verbose >= 3`); the progress bar's `verbose >= 2` gate documented;
  header conformance.
* `gl.Ho()` / `gl.He()`: the per-locus heterozygosity accessors now
  reject SilicoDArT (presence/absence) data with a fatal error --
  previously they silently returned meaningless values. Documentation
  brought to standard: gl.He's `@return` no longer claims to return
  observed heterozygosity, both functions state their semantics (pooled
  across all individuals; He is plain 2p(1-p) with no sample-size
  correction; NaN for all-NA loci; deliberately silent pure accessors)
  and cross-link gl.alf and gl.report.heterozygosity. Values verified
  exact against hand computation, gl.alf, and the population report.
* `gl.smearplot()`: six fixes. (1) The documented `plot.display` parameter
  was accepted and guarded but never used -- `print()` was unconditional,
  so the plot rendered even at `verbose = 0` or `plot.display = FALSE`;
  the print is now gated (saving via `plot.file` remains available either
  way). (2) The SilicoDArT legend has never shown its intended labels: a
  named assignment appended "Absence"/"Presence" to an unnamed c("0","1")
  vector instead of replacing, so the legend displayed "0" and "1"; it now
  reads Absence/Presence. (3) `het.only = TRUE` on SilicoDArT data warned
  "Set to FALSE" but rendered BOTH presence and absence in the het-only
  gray (#d3d3d3) because the palette was overridden before the datatype
  branch; the override now applies to SNP data only, and the warning is
  gated at `verbose >= 2` instead of printing at verbose 0. (4) The two
  package checks returned -1 after a cat() instead of raising an error;
  now `stop(error(...))`, and plotly (Suggests) is checked before
  `interactive = TRUE` use. (5) Dead "Missing data" label lines removed
  from both label blocks (NAs were stripped before the check, so the
  legend entry never appeared); documented that NA cells are colored
  without a legend entry. (6) Documentation and style: `group.pop` default
  documented as TRUE but is FALSE, stray duplicated doc line, grammar and
  quote fixes, `seq(1:n)` idiom, dead assignment in the dendrogram branch,
  a gated message when `den = TRUE` overrides `group.pop`, and the
  `loc.order` chromosome guard hardened against a zero-length chromosome
  slot (which would silently have dropped every locus).
* `gl.make.recode.pop()`: mirror of the gl.make.recode.ind fixes. The
  visible NULL return is now invisible with `@return` stating the actual
  contract (proforma written to file); the indented `@family` tag no
  longer leaks into the rendered help title; the `outpath` description
  no longer claims to save "plot RDS files"; and the @details wording
  now says populations where it said individuals. Proforma content and
  the round-trip through `gl.recode.pop()` were verified exact.
* `gl.make.recode.ind()`: the visible NULL return (a bare call printed
  "NULL") is now invisible, and the `@return` documentation -- which
  promised "A vector containing the new individual names" -- now states
  the actual contract: the proforma recode table is written to file and
  NULL is returned invisibly. The `outpath` parameter description no
  longer claims to save "plot RDS files" (a copy-paste from a plot
  function). Proforma content and the round-trip through
  `gl.recode.ind()` were verified exact on both datatypes.
* `gl.recode.ind()`: five fixes, mirroring `gl.recode.pop()`. (1) The
  `verbose = 3` deletions listing printed the literal word "Delete"
  instead of the deleted individuals' original identifiers (they had
  already been renamed before listing); it now lists the original names
  from the recode table, at `verbose >= 3` rather than only exactly 3.
  (2) verbose = 0 was not silent: the internal gl.drop.ind call always
  received both 'Delete' and 'delete' and warned about the absent one;
  only present spellings are now passed. (3) The object's locus-metric
  flags depended on verbosity (a pure renaming run invalidated them at
  `verbose >= 2` only); the misplaced reset is removed. (4) One call now
  appends exactly one history entry (the Delete path leaked the internal
  gl.drop.ind call). (5) The results summary gates at `verbose >= 3`
  (was >= 2), the return is invisible, the monomorphs-flag check is
  isFALSE()-guarded, and the @return no longer claims a genind can be
  returned.
* `gl.recode.pop()`: four fixes. (1) The `verbose = 3` deletions
  listing named the wrong individuals with the wrong count -- it indexed
  the individual names by a recycled per-population logical (e.g. 16
  arbitrary names listed when 20 individuals were deleted, none of them
  actual deletions); it now lists the true deletions, and appears at
  `verbose >= 3` rather than only exactly 3. (2) The object's
  locus-metric flags depended on verbosity: a pure renaming run (which
  leaves locus metrics valid) invalidated every flag at `verbose >= 2`
  but not below; the misplaced reset is removed (deletion runs are
  unaffected -- the internal gl.drop.pop already resets flags at every
  verbosity). (3) One call now appends exactly one history entry; the
  Delete path previously leaked the internal gl.drop.pop call as a
  second entry. (4) The return is invisible, and the monomorphs-flag
  check is isFALSE()-guarded against flag-less objects.
* `gl.report.factorloadings()`: five fixes. (1) The report line and
  top-N table printed at every verbosity level including `verbose = 0`;
  now gated at `verbose >= 1`. (2) `n.display` beyond the number of loci
  printed garbage NA rows (189 of them at n.display = 300 on a
  111-locus loadings table), and `n.display = 0` printed one row anyway
  (the `1:0` slip); the display is now clamped with `head()`. (3) An
  out-of-range `axis` produced a cryptic "subscript out of bounds"; now
  a clear fatal error. (4) The `@return` documentation claimed "The
  unchanged genlight object" -- the input is a glPca and the actual
  return is an invisible data.frame of the axis loadings; corrected.
  The `...` parameter, documented as passed to ggsave, was never
  forwarded; it now reaches the plot-save call. The plural
  `@family matched reports` no longer creates an orphan doc concept.
  (5) The glPca type check now uses inherits() with a proper error.
* `gl.filter.locmetric()`: major fixes. (1) `keep = "outside"` had
  never worked -- its condition (`metric <= lower AND >= upper`) is
  impossible whenever lower < upper, so every call crashed with
  "Subsetting resulted in zero loci". It now retains the exact
  complement of 'within' (metric < lower or > upper; 'within' keeps
  [lower, upper] inclusive). (2) An invalid `keep` value crashed with
  the cryptic "object 'x2' not found"; it now coerces to 'within' with
  a warning at `verbose >= 1`. A non-numeric metric produced factor
  warnings then a crash; it now stops with a clear fatal error (the
  check the report sibling already had). (3) Two irrelevant preambles
  removed: a full monomorphs scan run on every call solely to warn that
  monomorphic loci exist, and a block that silently stamped
  `pop = 'pop1'` onto objects without population assignments. (4) Loci
  with NA metric values (already removed, correctly) are now itemised
  in the `verbose >= 3` summary, and the return is invisible.
* `gl.report.locmetric()`: four fixes. (1) The summary statistics and
  quantile table printed at every verbosity level including
  `verbose = 0`; now gated at `verbose >= 1`. (2) The stats lines carried
  doubled labels ("Minimum      :  Min.   : 5.0") because `summary()`
  was applied to the one-column data.frame rather than the vector; they
  now print clean numeric values. (3) The "Retained" counts treated NA
  metric values as retained; they now exclude NAs (relevant for
  user-supplied custom metrics). (4) "1st quantile"/"3r quantile"
  corrected to quartile labels. This function already had the verbose-0
  plot guard, unconditional plot build, working plot.file-without-
  display, and invisible return.
* `gl.report.overshoot()`: the results (count and locus listing, or the
  no-overshoot message) printed at every verbosity level including
  `verbose = 0`; they are now gated at `verbose >= 1`. The locus listing
  no longer carries a stray trailing comma (a `paste0(..., sep = ",")`
  slip), and an unnecessary genlight subset used only for counting was
  removed. The help page no longer claims that plots and tabulations are
  saved to the tempdir (this function produces neither). The overshoot
  logic itself was verified correct against independent recomputation
  (testset.gl carries 21 genuine overshoot loci).
* `gl.filter.reproducibility()`: five fixes. (1) A dataset missing the
  RepAvg (SNP) or Reproducibility (SilicoDArT) metric was returned
  UNFILTERED with no error or warning; it now stops with the same fatal
  error as the report sibling. The function also demanded
  AlleleID/CloneID, metrics it never uses -- that check is removed, so
  datasets carrying the repeatability metric but no AlleleID now work.
  (2) One call appended two history entries (the internal `gl.drop.loc`
  delegation leaked its own); now exactly one, per the
  gl.filter.monomorphs precedent. (3) Loci with NA repeatability
  silently passed the filter; they are now removed and itemised in the
  `verbose >= 3` summary, matching gl.filter.rdepth/taglength. (4)
  `plot.file` with `plot.display = FALSE` crashed ("object 'p3' not
  found"); plots are now always built and the save works without
  displaying. (5) The out-of-range threshold warning is gated at
  `verbose >= 1`, the return is invisible, and the indented `@family`
  tag no longer leaks into the rendered help title.
* `gl.report.reproducibility()`: four fixes, same family as the
  rdepth/taglength reports (#255/#257). (1) `plot.file` with
  `plot.display = FALSE` crashed ("object 'p3' not found"); plots are now
  always built and the RDS save works without displaying. (2) The
  summary statistics and quantile table printed at every verbosity level
  including `verbose = 0`; now gated at `verbose >= 1`, and `verbose = 0`
  forces `plot.display = FALSE`. (3) The "Retained" counts treated loci
  with NA RepAvg/Reproducibility as retained; they now exclude NAs.
  (4) The "3r quartile" typo, and a 28-line dead commented block
  removed.
* `gl.filter.taglength()`: four fixes, mirroring `gl.filter.rdepth`
  (#256). (1) Loci with a missing (NA) TrimmedSequence silently corrupted
  the output -- genotypes were dropped while `loc.metrics` kept a
  garbage all-NA row for each, desyncing the object. NA-length loci are
  now removed cleanly and itemised in the `verbose >= 3` summary. (2)
  The progress message claimed loci *between* the thresholds are removed
  although those are the loci retained; it now reads "tag length < lower
  or > upper". (3) The threshold swap/range warnings printed at
  `verbose = 0` and the lower-range message named the wrong parameter
  ("'verbose'"); warnings now gate at `verbose >= 1` with corrected
  text. (4) The return is now invisible, and the indented `@family` tag
  no longer leaks into the rendered help title.
* `gl.report.taglength()`: four fixes, mirroring the rdepth pair
  (#255/#256). (1) `plot.file` with `plot.display = FALSE` crashed
  ("object 'p3' not found"); the composite plot is now always built and
  the RDS save works without displaying. (2) The summary statistics and
  quantile table printed at every verbosity level including
  `verbose = 0`; now gated at `verbose >= 1`, and `verbose = 0` forces
  `plot.display = FALSE`. (3) The "Retained" counts treated loci with NA
  tag length (NA TrimmedSequence) as retained; they now exclude NAs.
  (4) "1st quantile"/"3r quantile" corrected to "1st quartile"/"3rd
  quartile", and the plot title is now datatype-aware instead of always
  saying "SNP data".
* `gl.filter.rdepth()`: five fixes. (1) Loci with a missing (NA) read
  depth metric silently corrupted the output -- the genotype subset
  dropped them but the locus-metrics subset kept an all-NA row for each,
  leaving genotypes and `loc.metrics` out of sync. NA-depth loci are now
  removed cleanly (and counted in the `verbose >= 3` summary); datasets
  without NA read depths are unaffected. (2) Specifying `plot.file` with
  `plot.display = FALSE` crashed ("object 'p3' not found"); plots are now
  always built and the RDS save works without displaying. (3) The
  progress message claimed boundary loci are removed ("rdepth <= lower
  and >= upper") although they are retained; it now reads "< lower or
  > upper", matching the documented and actual behaviour. (4)
  `verbose = 0` now forces `plot.display = FALSE`, and the return is
  invisible, matching the other filter functions. (5) The indented
  `@family` tag leaked into the rendered help title
  ("... (read depth) @family matched filter"); the header has been
  rewritten and ?gl.filter.rdepth renders correctly.
* `gl.report.rdepth()`: four fixes. (1) Specifying `plot.file` with
  `plot.display = FALSE` crashed ("object 'p3' not found") because the
  plots were only built when displayed; plots are now always built and
  the RDS save works without displaying. (2) The summary statistics and
  quantile table printed at every verbosity level including
  `verbose = 0`; they are now gated at `verbose >= 1`, and `verbose = 0`
  forces `plot.display = FALSE` (matching the other report functions).
  (3) The "Retained" counts in the quantile table counted loci with NA
  read depth as retained; they now exclude NAs (tables change only for
  datasets with NA read-depth metrics). (4) The "3r quartile" typo.
* `gl.filter.secondaries()`: two behaviour fixes and cleanups.
  (1) `method = "best"` now actually retains the best SNP per sequence
  tag: previously the sort ran on the full `AlleleID` string, which is
  unique per locus, so the documented RepAvg/AvgPIC criteria could never
  engage and selection was effectively alphabetical (on platypus.gl, 5 of
  9 multi-SNP clones kept a lower-quality SNP, e.g. RepAvg 0.95 retained
  over 1.0). Loci retained under `method = "best"` will therefore change.
  (2) The output now preserves the input locus order; previously the
  object came back shuffled (method = "random") or AlleleID-sorted
  (method = "best") even when no secondaries were removed. Also: the
  invalid-`method` warning is gated at `verbose >= 1` (and `method` is
  actually set to "random"), the return is now invisible (matching the
  other filter functions), and an unreachable SilicoDArT check was
  removed.
* `gl.report.secondaries()`: four fixes. (1) The function crashed
  ("Subsetting resulted in zero loci") on any dataset with no secondaries
  -- an unused leftover subset (`x[, duplicated(b)]`) errored before the
  documented no-secondaries branch could run; that line is removed and
  such datasets now return the documented parameter data.frame. (2) The
  results block printed at every verbosity level including `verbose = 0`;
  it is now gated at `verbose >= 1`. (3) The 'TrimmedSequence not found'
  warning likewise printed at `verbose = 0`; now gated at `verbose >= 1`.
  (4) The raw per-iteration lambda estimates (several hundred lines on a
  typical dataset) printed at default `verbose = 2`; they now print only
  at `verbose >= 5`, with the existing "Converged on Lambda" summary
  retained at `verbose >= 2`.
* `gl.filter.allna()`: five fixes. (1) Removing all-NA individuals (with
  no all-NA loci) previously left the locus-metric flags stale
  (`CallRate` still TRUE though every locus's denominator changed) and
  recorded no history entry -- both were gated on locus-count change
  only; they now fire on any removal. (2) The return is invisible, so an
  unassigned call no longer prints the object summary. (3) The standard
  datatype check had been commented out; restored, so non-genlight input
  fails fast with a clear message. (4) `by.pop = TRUE` now records one
  history entry instead of leaking a second, internal `gl.drop.loc`
  entry. (5) The all-NA-individuals listing previously printed a literal
  "NULL" for every healthy individual; it now names only the affected
  individuals, with a count. Documentation also corrected (single
  correct `@family`, fixed title, tag order).
* `gl.filter.monomorphs()`: two behavioural fixes. (1) The return is now
  invisible, so an unassigned call no longer prints the full object
  summary (assigned use is unchanged). (2) Each call now records exactly
  one entry in `@other$history` -- previously the internal
  `gl.drop.loc()` delegation leaked a second, implementation-detail
  entry carrying the full list of removed locus names. Also: a
  works-by-coincidence `length(loc.list > 0)` corrected to
  `length(loc.list) > 0` (provably identical behaviour), and
  documentation tidied (Author(s) line, tag order, unused imports
  removed).
* `gl.report.allna()`: four fixes. (1) When individuals scored all-NA
  were present, the listing printed a literal "NULL" for every healthy
  individual, burying the real names -- it now names only the affected
  individuals, with a count. (2) The standard datatype check had been
  commented out; it is restored, so non-genlight input fails fast with
  the standard clear message. (3) The results are gated at
  `verbose >= 1` (previously printed at every verbosity level including
  0). (4) Documentation: real `@return` text (was the junk string
  "gl.report.allna"), the erroneous second `@family filter functions`
  tag removed, and the title corrected.
* `gl.rename.pop()`: three behaviour fixes. (1) Renaming a population
  that does not exist was a silent no-op that still recorded a history
  entry claiming the rename happened; it is now a fatal error listing
  the populations present. (2) Renaming TO an existing population name
  silently MERGED the two populations (an artefact of R's `levels<-`
  merge semantics); it is now a fatal error -- use `gl.recode.pop()` to
  amalgamate populations deliberately. (3) A pop-less object produced
  the cryptic "attempt to set an attribute on NULL"; now a clear fatal
  error. Also: invisible return, a redundant class check removed, and
  the header brought to the ratified template.
* `gl.report.monomorphs()`: the function previously returned the object
  with all monomorphic and all-NA loci silently REMOVED (an undocumented
  filtering side effect of deriving the counts), plus a phantom history
  entry -- so `gl <- gl.report.monomorphs(gl)` quietly deleted loci
  despite the documentation promising an unaltered return. It now
  returns the input untouched; use `gl.filter.monomorphs()` to actually
  remove them. The results block is also gated at `verbose >= 1` (it
  previously printed at every verbosity level including 0).
* `gl.filter.callrate()`: four fixes. (1) `method = 'pop'` previously
  returned an object whose locus metrics (`@other$loc.metrics`) were
  entirely NA -- the metadata was re-subset by locus name against row
  names that are not locus names; it is now subset positionally and the
  metrics come back intact. (2) An unset `loc.metrics.flags$monomorphs`
  no longer crashes with "argument is of length zero" (also guarded in
  the callee `utils.recalc.callrate()`). (3) Individuals with a call
  rate exactly equal to the threshold were retained but also printed in
  the "Individuals deleted" listing; the listing now shows only
  individuals actually removed. (4) After filtering individuals with
  `mono.rm = FALSE`, `loc.metrics.flags$monomorphs` is now set FALSE
  (monomorphic loci may have arisen), matching the behaviour of
  gl.drop.ind/gl.keep.ind. Also: verbose = 0 is now fully silent, and
  assorted documentation corrections (the input is a genlight, not a
  genind; no "summary" is returned; the threshold is a proportion, not
  an integer).
* `gl.report.hamming()` now reports the exact number of loci that
  `gl.filter.hamming()` would remove at candidate thresholds 0-10, by running
  the filter's own comparison engine in simulation (same worst-to-best
  call-rate ordering). Distances are reported as counts of mismatching bases
  over `min.length` bases, matching the filter's threshold units, and are
  computed in compiled code (the former O(n^2) R loop is gone). Arguments
  `tag.length` and `probar` are deprecated and ignored; a new `min.length`
  argument matches `gl.filter.hamming()`.

* `gl.filter.hamming()`: `threshold` is a count of mismatching bases (e.g. 3),
  no longer a proportion of sequence length as in earlier versions;
  proportion-style values (0 < threshold < 1) are rejected with an error.
  Comparable loci are ordered worst-to-best call rate before deduplication,
  so the retained locus of every duplicate pair is the one with the better
  call rate. The comparison engine is shared with `gl.report.hamming()` and
  compiled once per session.

* `gl.drop.ind()`: fixed a bug where locus-metric flags
  (`AvgPIC`, `OneRatioRef`, `OneRatioSnp`, `PICRef`, `PICSnp`, `CallRate`,
  `maf`, `FreqHets`, `FreqHomRef`, `FreqHomSnp`) were only reset to `FALSE`
  after dropping individuals when `verbose >= 2`. At `verbose = 0` or `1`
  with the default `recalc = FALSE`, these flags stayed `TRUE` (stale) even
  though the underlying statistics no longer reflected the retained
  individuals. They now reset correctly at every verbosity level, matching
  the behaviour `verbose >= 2` already had.

* `gl.keep.loc()`: three edge-case fixes. (1) `last` now defaults to the
  last locus in the dataset when omitted, as the documentation always
  stated -- previously `gl.keep.loc(gl, first = 100)` crashed with
  "argument is of length zero". (2) Calling with neither `loc.list` nor
  `first` now fails with a clear parameter error instead of
  "object 'flag' not found". (3) The out-of-range check on the locus range
  previously tested `first` but clamped `last`; an out-of-range `last` now
  warns and clamps as intended, and an out-of-range `first` is now a clear
  fatal error (previously it silently returned a single arbitrary locus).
  Zero-length `loc.list` input still returns the object unchanged with a
  warning, as before.
* `gl.report.bases()`: the results printout (sequence length, base
  frequencies, transitions/transversions) previously printed at every
  verbosity level, including `verbose = 0`. It is now gated at
  `verbose >= 1`, so fully-quiet calls are silent as documented. The
  default (`verbose = 2`) behaviour is unchanged.
* `gl.keep.ind()`: fixed the identical bug -- the same locus-metric flags
  were only reset to `FALSE` after retaining individuals when
  `verbose >= 2`. Now reset correctly at every verbosity level.

* `gl.drop.loc()`: four fixes. (1) The "not present in the dataset"
  warning previously named the wrong loci -- it indexed the dataset's
  locus names with positions from the user's `loc.list`, so a typo in one
  locus name produced a warning about a different, valid locus. It now
  names the loci the user actually listed. (2) `last` now defaults to the
  last locus when omitted, as documented -- previously
  `gl.drop.loc(gl, first = 100)` crashed. (3) The out-of-range check on
  the locus range tested `first` but clamped `last`; an out-of-range
  `last` now warns and clamps, and an out-of-range `first` is a clear
  fatal error (previously it silently dropped a single arbitrary locus).
  (4) The range-clamp warnings are now silent at `verbose = 0`.
* `gl.report.callrate()`: five fixes. (1) The results tables previously
  printed at every verbosity level including `verbose = 0`; they are now
  gated at `verbose >= 1` (default `verbose = 2` output unchanged).
  (2) The returned object is now truly unaltered, as the documentation
  always stated -- previously it came back with the `CallRate` locus
  metric recalculated in place; the report still uses freshly
  recalculated values internally. (3) An unknown `method` (anything other
  than "loc"/"ind") now warns and coerces to "loc" instead of silently
  producing no output at all. (4) Two documented examples calling a
  nonexistent `by.pop` argument (silently swallowed by `...`) were
  removed. (5) Cosmetic: a stray ")" printed after the individuals table,
  the "3r quartile" typo in both branches, `ind.to.list = 0` listing one
  individual instead of none, and the `@param bins` default documented as
  25 when the signature default is 50.
* `gl.alf()`: retained as the documented fast path for per-locus allele
  frequencies rather than deprecated. The inert commented-out deprecation
  notice in the function body and the "gl.alf is deprecated" developer
  note in `gl.allele.freq.r` are removed: the nominated replacement,
  `gl.allele.freq(x, simple = TRUE)`, rounds to 4 decimal places, costs
  5-19x more because it splits the object by population first, cannot be
  passed as a bare one-argument function to `lapply()` or
  `utils.jackknife()`, and rejects genlight objects not built by dartR --
  and after the companion SilicoDArT fix it will no longer return the same
  quantity for Tag P/A data. Fifteen call sites across dartR.base,
  dartR.captive, dartR.popgen and dartR.sim depend on the fast path.
  Two behaviour fixes. (1) Duplicate locus names silently discarded the
  locus keys: `data.frame()` substitutes row names `1:n` whenever the
  names it is given are not unique, so callers that read
  `rownames(gl.alf(x))` as locus names got integers --
  `gl.report.heterozygosity(method = 'pop')` then reported every locus as
  polymorphic and returned negative monomorphic-locus counts. Row names
  are now set from `locNames()` and disambiguated with `make.unique()`,
  keeping positional correspondence with the loci of the input. (2) The
  function now rejects SilicoDArT data with a fatal error; presence/
  absence data has ploidy 1 but was divided by the SNP ploidy, so `alf2`
  was half the presence frequency and `alf1` was not a complement of
  anything. Values for SNP data are unchanged. Documentation brought to
  standard: `@return` named columns "ref" and "alt" that the function has
  never produced (they are `alf1` and `alf2`), and the header gains a
  description, details stating the semantics and the relationship to
  `gl.allele.freq()`, an author/custodian line, and cross-links.
