# dartR.base 1.2.3 (development)

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
