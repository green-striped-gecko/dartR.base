# dartR.base 1.2.3 (development)

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
