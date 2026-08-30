# dartR.base 1.2.3 (development)

* `gl.filter.locmetric()`: major fixes. (1) `keep = "outside"` had
  never worked -- its condition (`metric <= lower AND >= upper`) is
  impossible whenever lower < upper, so every call crashed with
  "Subsetting resulted in zero loci". It now retains the exact
  complement of 'within' (metric < lower or > upper; 'within' keeps
  [lower, upper] inclusive). (2) An invalid `keep` value crashed with
  the cryptic "object 'x2' not found"; it now coerces to 'within' with
  a warning at `verbose >= 1`. A non-numeric metric produced factor
  warnings then a crash; it now stops with a clear fatal error (the
  check the report sibling already had). (3) Two irrelevant preambles
  removed: a full monomorphs scan run on every call solely to warn that
  monomorphic loci exist, and a block that silently stamped
  `pop = 'pop1'` onto objects without population assignments. (4) Loci
  with NA metric values (already removed, correctly) are now itemised
  in the `verbose >= 3` summary, and the return is invisible.

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
