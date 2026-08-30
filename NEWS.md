# dartR.base 1.2.3 (development)

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
