# Audit: population colours on the wrong populations (2026-10-22)

Follow-up to the `gl.report.heterozygosity` fix (dartR.base PR #455). The question was which other dartRverse functions put per-population colours on the wrong populations. The code checked was `origin/dev` of each package: dartR.base cf29619, dartR.popgen 94b630c, dartR.spatial 0aed56d, dartR.captive f91737e, dartR.sim 6a70557, plus dartR.sexlinked. Nothing was changed.

**Test object:** `platypus.gl` with `pop(x) <- factor(pop, levels = c("SEVERN_BELOW", "TENTERFIELD", "SEVERN_ABOVE"))`. In this object, first-appearance order, level order and alphabetical order all differ. Two inputs were tried:
- named colours `c(TENTERFIELD = "green", SEVERN_ABOVE = "red", SEVERN_BELOW = "blue")`;
- unnamed colours `c("red", "blue", "green")`, where the expected result is SEVERN_BELOW red, TENTERFIELD blue, SEVERN_ABOVE green.

`gl.genleastcost` used possums A/B/C with levels B, C, A instead. The fixed `gl.report.heterozygosity` passed every case, which shows the test harness works.

Mechanisms:
- (a) colours assigned in first-appearance order;
- (b) user names overwritten or ignored;
- (c) colours assigned by rank;
- (d) unnamed colours passed by position while the plotted column is a character column (ggplot sorts it alphabetically) or a factor whose levels differ from `levels(pop(x))`.

## Tier 1: wrong with default or unnamed colours, whenever populations are not listed in level order

| Package | Function | Mechanism | Lines | Evidence |
|---|---|---|---|---|
| base | gl.report.polyploid_heterozygosity | a, b | :695, :704, :854, :865, :887 | Same code as heterozygosity before #455. With `n.invariant = 0`, unnamed colours gave SB blue, T red, SA green. With `n.invariant > 0`, named colours were also wrong. The subsample panels were wrong in both cases. |
| popgen | gl.LDNe | b, d | :514, :527, :529, :550, :552 | Colours are attached by row, while the bars use an alphabetical factor. **The x-axis labels are attached to the wrong bars**: the SEVERN_ABOVE bar is labelled "SEVERN_BELOW \| 15". The plot is only right when the levels are alphabetical and the colours are unnamed. |
| base | gl.report.pa (pairwise sankey) | a, b | :378–380 | `d3.scaleOrdinal()` has no domain, so colours go to populations in the order links are drawn. Emulation gave T red, SA blue, SB green. The plot was not rendered. |

## Tier 2: unnamed colours go in alphabetical order (wrong only when the levels are not alphabetical)

| Package | Function | Lines | Note |
|---|---|---|---|
| base | gl.report.allelerich | :583–586, :603 | `fill = pop` is a character column. Named colours are OK. |
| base | gl.report.ld.map | :308, :359–362 | Same mechanism. Populations dropped by `ind.limit` probably shift the colours too (not tested). |
| spatial | gl.spatial.autoCorr | :628, :653 | `rbindlist(idcol = "Population")` produces a character column. |
| popgen | gl.ld.distance | :120, :122, :159 | Unnamed colours go in alphabetical order, not in the population order of `ld.report`. |
| spatial | gl.run.eems | :459–466 | Same mechanism as autoCorr. Not run (needs the EEMS binary). |

## Tier 3: level order is right, but user names are overwritten or ignored

| Package | Function | Lines |
|---|---|---|
| base | gl.plot.heatmap | :235 |
| base | gl.report.contamination | :257 |
| base | gl.map.interactive | :180 |
| base | gl.report.pa (one2rest) | :507–511, :556. Also, "Rest" reuses colour 1 with the default palette. |
| spatial | gl.grm2 | :185 |
| spatial | gl.genleastcost | :269 |
| captive | gl.grm | :188 |
| captive | gl.grm.network | :373 |
| sim | gl.report.nall | :253 |
| popgen | gl.plot.snmf (`den = TRUE`) | :206–207. The `is.null(names())` guard is commented out. |

## OK

- **dartR.base:** gl.report.diversity, gl.pcoa.plot (2D; the 3D/plotly path was not run).
- **dartR.popgen:** gl.plot.structure, gl.plot.popcluster, gl.plot.faststructure, gl.map.structure and gl.map.snmf colour by cluster, not by population. gl.ld.haplotype and gl.run.bayescan colour by category. gl.run.mongrail does not colour populations.
- **dartR.captive:** gl.plot.network, gl.parentage.power, gl.report.gd.projection, gl.report.kinship, gl.run.EMIBD9 and gl.relatedness. gl.relatedness passes its colours on to gl.plot.heatmap.
- **dartR.sim:** gl.diagnostics.sim, gl.sim.emigration and gl.sim.ind.af.
- **dartR.sexlinked:** its two fixed `plot.colors` are not per population.

## Fix pattern

The fix is the #455 block. Keep the user's names when they cover every `popNames(x)`; otherwise name the colours by `popNames(x)`. Subset with `[popNames(x)]`, and then:
- index colours by population name, never by row position;
- pass the named vector to `scale_*_manual`, or make the plotted column `factor(..., levels = popNames(x))`.

For the d3 sankey, add `.domain([...])` built from the named vector.

## Separate issue found

In gl.plot.snmf and gl.plot.popcluster, the cluster labels "Pop_1"…"Pop_K" are text, so they sort alphabetically. With K ≥ 10, the second colour goes to Pop_10.
