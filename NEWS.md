# dartR.base 1.2.3 (development)

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
