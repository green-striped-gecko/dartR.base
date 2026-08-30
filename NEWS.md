# dartR.base 1.2.3 (development)

* `gl.select.colors()`: five fixes. (1) An unrecognised `library` value
  silently returned base R's `colors()` FUNCTION as the "colour vector"
  (the internal variable was never assigned and lexical scoping found
  grDevices::colors); unknown libraries now coerce to the default
  scales/hue_pal with a warning at `verbose >= 1`. (2) Brewer requests
  are honoured: fewer than 3 colours are trimmed from the 3-colour pull
  (you get exactly what you asked for), and requests above the palette
  maximum return the maximum with a clear gated warning (previously 2
  requested delivered 3, and 12 requested from Blues silently delivered
  9). (3) Out-of-bounds `select` indices, which produced NA colours,
  are now a fatal error. (4) baseR palette='heat' passed validation but
  returned rainbow; it now dispatches to heat.colors. (5) The datatype
  banner no longer prints at `verbose = 0` when a genlight object is
  supplied, and missing packages raise a fatal error instead of
  returning -1. The internal default-colour contract used across the
  package (brewer, Blues, select=c(7,5)) is unchanged and
  regression-tested.

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
