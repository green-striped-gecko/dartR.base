# dartR.base 1.2.3 (development)

* `gl.filter.replicates()`: four fixes. (1) When the members of a
  replicate pair had tied missing-data rates -- the exact-duplicate case
  -- BOTH were removed (the doubled pair table evaluated the drop rule in
  both orientations); pairs are now canonicalised and deduplicated before
  the rule is applied, so one member per pair is removed (ties remove the
  alphabetically later individual), and this works with report tables
  generated both before and after the companion gl.report.replicates
  fix. (2) Re-thresholding to an empty set crashed via gl.drop.ind ("no
  individuals to drop"); the object is now returned unchanged with a
  gated message. (3) `replicates.report` was never validated -- the
  report's old no-pairs string return crashed with "$ operator is invalid
  for atomic vectors"; malformed input is now a clear fatal error, and an
  empty table is handled as nothing-to-drop. (4) The history of the
  returned object recorded gl.drop.ind's call instead of
  gl.filter.replicates' (gl.drop.ind now runs quietly, the dropped
  individuals are itemised at `verbose >= 2`, and a single
  gl.filter.replicates entry is appended); a datatype check was added and
  the object is returned invisibly per house standard.

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
