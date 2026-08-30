# dartR.base 1.2.3 (development)

* `gl.load()`: the `compliance` parameter (documented [default FALSE]) was
  inert -- `gl.compliance.check()` ran unconditionally on every load,
  potentially modifying the loaded object and printing its full output
  even at `verbose = 0`. As directed by the module coordinator, the check
  now runs only when `compliance = TRUE` (with verbose passed through).
  The "Loaded object" and fbm-conversion messages are gated at
  `verbose >= 2` (they printed at verbose 0). A missing file and an RDS
  that does not contain a genlight object now give clear fatal errors
  (previously a raw connection error and a cryptic "no applicable method
  for `@`" failure). Docs corrected: `fbm` default is FALSE (was
  documented TRUE), `file` is the file to read (was "receive data"), an
  examples block added.

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
