# dartR.base 1.2.3 (development)

* `gl2genalex()`: SilicoDArT (presence/absence) objects were silently
  exported as meaningless codominant GenAlEx genotypes; the function now
  accepts SNP data only. A missing poppr installation raises a proper
  error condition instead of returning -1 (scripts relying on the -1
  return now see an error). Also: poppr's conversion chatter is
  suppressed below `verbose = 2`, all-NA loci dropped before export are
  reported at `verbose >= 1`, the NULL return is invisible, and the
  `overwrite = FALSE` doc now states the actual raise-on-exists
  behaviour.
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
  individually but never compared; now fatal. (5) Messages used
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
