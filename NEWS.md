# dartR.base 1.2.3 (development)

* `gl.drop.ind()`: fixed a bug where locus-metric flags
  (`AvgPIC`, `OneRatioRef`, `OneRatioSnp`, `PICRef`, `PICSnp`, `CallRate`,
  `maf`, `FreqHets`, `FreqHomRef`, `FreqHomSnp`) were only reset to `FALSE`
  after dropping individuals when `verbose >= 2`. At `verbose = 0` or `1`
  with the default `recalc = FALSE`, these flags stayed `TRUE` (stale) even
  though the underlying statistics no longer reflected the retained
  individuals. They now reset correctly at every verbosity level, matching
  the behaviour `verbose >= 2` already had.

* `gl.report.bases()`: the results printout (sequence length, base
  frequencies, transitions/transversions) previously printed at every
  verbosity level, including `verbose = 0`. It is now gated at
  `verbose >= 1`, so fully-quiet calls are silent as documented. The
  default (`verbose = 2`) behaviour is unchanged.
