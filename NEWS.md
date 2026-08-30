# dartR.base 1.2.3 (development)

* `gl.colors()`: three fixes. (1) Both invalid-type exits used
  `cat(error(...))` followed by `stop(-1)`, so the error condition an
  upstream `tryCatch()` received carried the message "-1" while the real
  message printed to stdout where even `try(silent = TRUE)` could not
  suppress it; they now use `stop(error(...))`. (2) The return is now
  visible -- the documented example `gl.colors(2)` previously displayed
  nothing because the result was returned `invisible()`. (3) Documentation:
  the description listed a `"pal"` category that the function has never
  accepted (it was a fatal error) -- removed; the implemented but
  undocumented `"structure"` type (35 discrete colors, used by
  `dartR.popgen::gl.plot.snmf`) is now documented, along with header
  conformance fixes and an accurate `@return` (the four palette types
  return a function, not a vector). Note: `gl.colors()` evaluated as a
  default argument in a caller's signature still prints its banner at the
  caller's `verbose = 0` unless the default passes `verbose = 0`; changing
  the default verbosity was considered and deliberately not adopted.

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
