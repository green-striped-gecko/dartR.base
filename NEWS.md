# dartR.base 1.2.3 (development)

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
