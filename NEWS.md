# dartR.base 1.2.3 (development)

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
