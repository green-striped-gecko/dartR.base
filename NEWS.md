# dartR.base 1.2.3 (development)

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
