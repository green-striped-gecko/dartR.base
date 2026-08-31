# dartR.base 1.2.3 (development)

* `gl.join()`: five fixes. (1) A join by shared loci LOST the individual
  metrics entirely -- plain `rbind()` returns NULL metadata and the
  function never rebuilt it; the combined object now carries the
  row-bound ind.metrics of both inputs (with the id column re-synced when
  duplicate names are made unique). (2) A join by shared individuals
  CRASHED ("replacement has 0 rows") on objects whose loc.metrics.flags
  data.frame lacks the OneRatio/PIC columns -- which is standard
  SNP-report data; only flags present in both objects are now combined.
  (3) The same flags block was triplicated (each path set the flags and a
  third copy ran again for both), warnings printed at `verbose = 0`
  (method deprecation, missing metrics/flags), and two
  `cat(error()) + stop()` splits returned no condition message; the
  duplicate block is gone, warnings are gated at `verbose >= 2`, and the
  exits use `stop(error(...))`. (4) SNP and SilicoDArT objects with
  matching names could be joined silently -- the datatypes were checked
  individually but never compared; now fatal. (4a, amendment) The
  historical legacy values `method='end2end'` and `method='sidebyside'`
  -- accepted by the pre-refactor implementation and described in the
  documentation ever since, and still used by real callers (the
  dartR.popgen gl.assign functions) -- were fatal because the legacy shim
  mapped only join.by.loc/join.by.ind; they are now mapped to their
  historical meanings, and any explicitly requested join is validated
  against the data so a mismatch fails with a clear message rather than
  a cryptic cbind/rbind error. (5) Messages used
  `substitute()` inside `cat()`, which printed garbage at `verbose = 2`
  and crashed at `verbose >= 3` whenever the arguments were expressions
  rather than names (e.g. `gl.join(x[1:7, ], y)`); arguments are now
  deparsed once for display. Docs: the description claimed the history was
  cleared (it is carried from the first object and appended); @details
  described method='sidebyside'/'end2end' values that never existed;
  duplicate/incorrect end-of-run summary removed; typos.

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
