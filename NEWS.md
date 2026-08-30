# dartR.base 1.2.3 (development)

* `gl.select.shapes()`: four fixes. (1) The range validation was a
  parenthesis slip (`min(select < 0 | max(select > 25))`) that only fired
  when every element was negative -- a partially-negative `select` such as
  `c(-1, 5)` passed straight through to `pch`. It now correctly rejects any
  value outside 0-25. (2) The documented `x=` genlight argument was
  non-functional: `nPop(x)` was computed and discarded, so a `select` of
  the wrong length passed silently and a NULL `select` returned all 26
  shapes regardless of the number of populations. Now (matching
  `gl.select.colors`) a length mismatch with `nPop(x)` is a fatal error, a
  NULL `select` with `x` returns one shape per population, and more than
  26 populations without an explicit `select` is a fatal error (only 26
  distinct shapes exist). (3) New `plot.display` parameter (default TRUE);
  the palette chart was previously drawn unconditionally and can now be
  suppressed, and is suppressed automatically at `verbose = 0` (which
  previously also leaked the datatype banner). (4) Cosmetic: the genlight
  argument `x` was shadowed by the plot x-coordinates mid-function, a
  "Requires shapes" typo, and header conformance.

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
