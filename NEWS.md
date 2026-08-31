# dartR.base 1.2.3 (development)

* `gl.report.monomorphs()`: the function previously returned the object
  with all monomorphic and all-NA loci silently REMOVED (an undocumented
  filtering side effect of deriving the counts), plus a phantom history
  entry -- so `gl <- gl.report.monomorphs(gl)` quietly deleted loci
  despite the documentation promising an unaltered return. It now
  returns the input untouched; use `gl.filter.monomorphs()` to actually
  remove them. The results block is also gated at `verbose >= 1` (it
  previously printed at every verbosity level including 0).
* `gl.filter.callrate()`: four fixes. (1) `method = 'pop'` previously
  returned an object whose locus metrics (`@other$loc.metrics`) were
  entirely NA -- the metadata was re-subset by locus name against row
  names that are not locus names; it is now subset positionally and the
  metrics come back intact. (2) An unset `loc.metrics.flags$monomorphs`
  no longer crashes with "argument is of length zero" (also guarded in
  the callee `utils.recalc.callrate()`). (3) Individuals with a call
  rate exactly equal to the threshold were retained but also printed in
  the "Individuals deleted" listing; the listing now shows only
  individuals actually removed. (4) After filtering individuals with
  `mono.rm = FALSE`, `loc.metrics.flags$monomorphs` is now set FALSE
  (monomorphic loci may have arisen), matching the behaviour of
  gl.drop.ind/gl.keep.ind. Also: verbose = 0 is now fully silent, and
  assorted documentation corrections (the input is a genlight, not a
  genind; no "summary" is returned; the threshold is a proportion, not
  an integer).
* `gl.report.hamming()` now reports the exact number of loci that
  `gl.filter.hamming()` would remove at candidate thresholds 0-10, by running
  the filter's own comparison engine in simulation (same worst-to-best
  call-rate ordering). Distances are reported as counts of mismatching bases
  over `min.length` bases, matching the filter's threshold units, and are
  computed in compiled code (the former O(n^2) R loop is gone). Arguments
  `tag.length` and `probar` are deprecated and ignored; a new `min.length`
  argument matches `gl.filter.hamming()`.

* `gl.filter.hamming()`: `threshold` is a count of mismatching bases (e.g. 3),
  no longer a proportion of sequence length as in earlier versions;
  proportion-style values (0 < threshold < 1) are rejected with an error.
  Comparable loci are ordered worst-to-best call rate before deduplication,
  so the retained locus of every duplicate pair is the one with the better
  call rate. The comparison engine is shared with `gl.report.hamming()` and
  compiled once per session.

* `gl.drop.ind()`: fixed a bug where locus-metric flags
  (`AvgPIC`, `OneRatioRef`, `OneRatioSnp`, `PICRef`, `PICSnp`, `CallRate`,
  `maf`, `FreqHets`, `FreqHomRef`, `FreqHomSnp`) were only reset to `FALSE`
  after dropping individuals when `verbose >= 2`. At `verbose = 0` or `1`
  with the default `recalc = FALSE`, these flags stayed `TRUE` (stale) even
  though the underlying statistics no longer reflected the retained
  individuals. They now reset correctly at every verbosity level, matching
  the behaviour `verbose >= 2` already had.

* `gl.keep.loc()`: three edge-case fixes. (1) `last` now defaults to the
  last locus in the dataset when omitted, as the documentation always
  stated -- previously `gl.keep.loc(gl, first = 100)` crashed with
  "argument is of length zero". (2) Calling with neither `loc.list` nor
  `first` now fails with a clear parameter error instead of
  "object 'flag' not found". (3) The out-of-range check on the locus range
  previously tested `first` but clamped `last`; an out-of-range `last` now
  warns and clamps as intended, and an out-of-range `first` is now a clear
  fatal error (previously it silently returned a single arbitrary locus).
  Zero-length `loc.list` input still returns the object unchanged with a
  warning, as before.
* `gl.report.bases()`: the results printout (sequence length, base
  frequencies, transitions/transversions) previously printed at every
  verbosity level, including `verbose = 0`. It is now gated at
  `verbose >= 1`, so fully-quiet calls are silent as documented. The
  default (`verbose = 2`) behaviour is unchanged.
* `gl.keep.ind()`: fixed the identical bug -- the same locus-metric flags
  were only reset to `FALSE` after retaining individuals when
  `verbose >= 2`. Now reset correctly at every verbosity level.

* `gl.drop.loc()`: four fixes. (1) The "not present in the dataset"
  warning previously named the wrong loci -- it indexed the dataset's
  locus names with positions from the user's `loc.list`, so a typo in one
  locus name produced a warning about a different, valid locus. It now
  names the loci the user actually listed. (2) `last` now defaults to the
  last locus when omitted, as documented -- previously
  `gl.drop.loc(gl, first = 100)` crashed. (3) The out-of-range check on
  the locus range tested `first` but clamped `last`; an out-of-range
  `last` now warns and clamps, and an out-of-range `first` is a clear
  fatal error (previously it silently dropped a single arbitrary locus).
  (4) The range-clamp warnings are now silent at `verbose = 0`.
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
