# dartR.base 1.2.3 (development)

* `gl.report.hwe()` / `gl.filter.hwe()`: the functionality of the
  deprecated `gl.report.excess.het()` / `gl.filter.excess.het()` has been
  migrated into this pair. Both functions gain `direction` ('both',
  'excess', 'deficit') to restrict attention to heterozygote excess or
  deficit (an expected-heterozygote column, Het.exp, is now included in
  the report output), and `min.hobs` to restrict testing to loci with
  observed heterozygosity at or above a threshold, applied before the
  multiple-comparison adjustment. The published Robledo-Ruiz et al.
  (2023) excess-heterozygosity workflow is reproduced with
  direction='excess', min.hobs=0.5, ChiSquare test and fdr adjustment;
  the deprecated functions remain as thin wrappers that emit a
  deprecation warning and produce the same flagged/removed loci as
  before (verified exactly on the LBP dataset), and will be removed in a
  future release. Migrating also retires a defect in the old
  gl.filter.excess.het, which computed its per-population genotype
  counts on the wrong individuals (a recycled per-population index).

  Review fixes applied to the pair in the same change: (1) the skipping
  of monomorphic and small-sample populations happened only at
  `verbose >= 2`, so the tested set -- and any multiple-comparison
  pool -- depended on verbosity (verbose 0 tested all 30 testset.gl
  populations, verbose 2 only 23); populations are now skipped at every
  verbosity. (2) A missing HardyWeinberg package now raises a fatal
  error instead of returning -1, and ggtern is required only when
  ternary plots are requested (plot.out=TRUE). (3) The out-of-range
  alpha warning is gated at `verbose >= 1` and no longer calls alpha an
  "integer". (4) gl.report.hwe now carries a @family tag, and
  gl.filter.hwe's duplicated/indented family tags are repaired.

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
