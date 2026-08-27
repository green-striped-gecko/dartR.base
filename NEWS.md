# dartR.base 1.2.3 (development)

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

* `gl.keep.ind()`: fixed the identical bug -- the same locus-metric flags
  were only reset to `FALSE` after retaining individuals when
  `verbose >= 2`. Now reset correctly at every verbosity level.
