# dartR.base 1.2.3 (development)

* `gl.smearplot()`: six fixes. (1) The documented `plot.display` parameter
  was accepted and guarded but never used -- `print()` was unconditional,
  so the plot rendered even at `verbose = 0` or `plot.display = FALSE`;
  the print is now gated (saving via `plot.file` remains available either
  way). (2) The SilicoDArT legend has never shown its intended labels: a
  named assignment appended "Absence"/"Presence" to an unnamed c("0","1")
  vector instead of replacing, so the legend displayed "0" and "1"; it now
  reads Absence/Presence. (3) `het.only = TRUE` on SilicoDArT data warned
  "Set to FALSE" but rendered BOTH presence and absence in the het-only
  gray (#d3d3d3) because the palette was overridden before the datatype
  branch; the override now applies to SNP data only, and the warning is
  gated at `verbose >= 2` instead of printing at verbose 0. (4) The two
  package checks returned -1 after a cat() instead of raising an error;
  now `stop(error(...))`, and plotly (Suggests) is checked before
  `interactive = TRUE` use. (5) Dead "Missing data" label lines removed
  from both label blocks (NAs were stripped before the check, so the
  legend entry never appeared); documented that NA cells are colored
  without a legend entry. (6) Documentation and style: `group.pop` default
  documented as TRUE but is FALSE, stray duplicated doc line, grammar and
  quote fixes, `seq(1:n)` idiom, dead assignment in the dendrogram branch,
  a gated message when `den = TRUE` overrides `group.pop`, and the
  `loc.order` chromosome guard hardened against a zero-length chromosome
  slot (which would silently have dropped every locus).

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
