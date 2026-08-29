# dartR.base 1.2.3 (development)

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
