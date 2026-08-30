# dartR.base 1.2.3 (development)

* `gl.fixed.diff()`: four fixes. (1) The documented `mono.rm` parameter
  ([default TRUE], remove monomorphic loci before computation) was never
  referenced in the body, and the flag logic standing in its place was
  inverted -- `gl.filter.monomorphs()` ran only when the flag certified
  the data already monomorph-free, so monomorphic loci were effectively
  never removed. `$nloc` was inflated and `$pcfd` denominators with it
  (roughly fivefold understatement of percent fixed differences on a test
  subset with 208/255 monomorphic loci). `mono.rm` is now honoured:
  TRUE (default) removes monomorphs, FALSE retains them with the warning
  and reproduces the previous numbers. Raw `$fd` counts are unaffected
  (a monomorphic locus cannot carry a fixed difference), so
  `gl.fdsim()`/`gl.collapse()` amalgamation decisions are unchanged;
  `gl.dist.pop(method = "fixed-diff")` distances change with the
  corrected `$pcfd`. (2) The "false positives can only be simulated for
  tloc=0" warning printed at `verbose = 0`; now gated at `verbose >= 2`.
  A dead, unreachable `tloc.hold` block inside the pairwise loop (with
  its own ungated output) was removed. (3) The `verbose >= 4` return
  listing misnamed `$pval` as `$prob`, described `$sdfpos` with the
  copy-pasted expected-count text, and called `$gl` the input object;
  matrix dimnames and diagonals are now set once instead of on every
  pass of the pairwise loop. (4) Documentation: the details promised a
  per-comparison warning at sample sizes below 5 that never existed
  (actual behaviour: a global minimum-sample warning at n < 10,
  `verbose >= 3`); the progress bar's `verbose >= 2` gate documented;
  header conformance.

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
