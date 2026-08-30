# dartR.base 1.2.3 (development)

* `gl.report.replicates()`: four fixes. (1) The pair table carried BOTH
  orderings of every pair, and the drop rule picked the opposite member in
  each ordering whenever missing-data rates tied -- which is precisely the
  exact-duplicate case -- so BOTH replicates landed in `ind.list.drop` and
  the histograms double-counted every pair. The table now holds one row
  per unordered pair (ind1 = the individual earlier in the object) and a
  tie deterministically drops ind2. (2) When no pairs passed the
  thresholds the function returned a bare character message instead of the
  documented 3-element list, crashing `gl.filter.replicates()` downstream;
  it now returns the documented structure with an empty table and a gated
  message. (3) `ind.list.rep` used `>= perc_geno` where every other output
  used `>` (aligned), and its NaN diagonal previously injected NA entries;
  a datatype check was added, Rcpp/RcppParallel are checked before
  compiling (RcppParallel added to Suggests), and the plot no longer
  renders at `verbose = 0`. (4) The result is returned invisibly (house
  standard) with the pair table printed at `verbose >= 3`; docs cleaned
  (@family, literal "##" headings, typos).

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
