# dartR.base 1.2.3 (development)

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
