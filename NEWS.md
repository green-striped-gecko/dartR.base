# dartR.base 1.2.3 (development)

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

* `gl.drop.ind()`: fixed a bug where locus-metric flags
  (`AvgPIC`, `OneRatioRef`, `OneRatioSnp`, `PICRef`, `PICSnp`, `CallRate`,
  `maf`, `FreqHets`, `FreqHomRef`, `FreqHomSnp`) were only reset to `FALSE`
  after dropping individuals when `verbose >= 2`. At `verbose = 0` or `1`
  with the default `recalc = FALSE`, these flags stayed `TRUE` (stale) even
  though the underlying statistics no longer reflected the retained
  individuals. They now reset correctly at every verbosity level, matching
  the behaviour `verbose >= 2` already had.

* `gl.keep.ind()`: fixed the identical bug -- the same locus-metric flags
  were only reset to `FALSE` after retaining individuals when
  `verbose >= 2`. Now reset correctly at every verbosity level.

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
