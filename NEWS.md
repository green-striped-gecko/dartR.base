# dartR.base 1.2.3 (development)

* `gl.report.allelerich()`: plumbing fixes; the rarefaction calculation
  itself was verified against an independent recomputation and is
  unchanged. (1) The plot rendered at `verbose = 0` (missing
  `plot.display` guard) and the lazy signature default
  `gl.colors("dis")` printed a 3-line banner at every verbosity -- both
  silenced (the default now passes `verbose = 0`). (2) An unrecognized
  `error.bar` value crashed downstream with "object 'max_val' not
  found"; unknown values now coerce to "SD" with a gated warning, and the
  silent override of the user's error-bar choice when `nboots > 0` is
  announced at `verbose >= 2`. (3) The package check called
  `requireNamespace()` on a vector, silently checking only dplyr, and
  returned -1 through a cat(); each package is now checked individually
  with `stop(error(...))`, and boot/Rcpp (needed for bootstrapping) are
  checked before the bootstrap path. (4) Dead code removed (an unused
  first plot built on every call, a commented parallel block, duplicated
  global declarations); `boot.method` validated; header canon.

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
