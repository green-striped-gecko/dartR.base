# Review: gl.ld.haplotype (dartR.popgen)

- Family mode: analysis (LD computation plus plotting; returns a haplotype table)
- Date: 2026-09-18
- Reviewer: Claude (claude-fable-5-1), dartr-function-review v2.0.0
- Package commit: 31fb40c (dartR.popgen `dev_luis`, synced with `origin/dev`)
- Datasets: platypus.gl (dartR.data 1.2.5) as in the documented example
  (TENTERFIELD, first 15 individuals, call rate 1, chromosome and position
  slots filled from `Chrom_Platypus_Chrom_NCBIv1` /
  `ChromPos_Platypus_Chrom_NCBIv1`: 569 loci); testset.gl and testset.gs for
  the no-chromosome path
- Environment: R 4.4.2 (framework), dartR.base 1.2.3, snpStats 1.56.0,
  raster 3.6.32, sp 2.2.1, ggplot2 with `fortify(<SpatialPolygonsDataFrame>)`
  deprecated
- Baseline: `tests/testthat/test-gl.ld.haplotype.R` in dartR.popgen (new
  file; 29 assertions, all passing on the reviewed state). It pins the
  returned table for the default call and for `haplo_id = TRUE`, the files
  written, the non-silent `verbose = 0`, the `ind.limit` boundary and the
  error on an object without chromosome information.
- Verification method: a scratch copy of the function with the local
  environment captured at the end of each chromosome iteration
  (`as.list(environment())`), run beside the package function. No package
  file was edited.
- Community check: no thread on the dartR Google Group and no GitHub issue
  mentions `gl.ld.haplotype` (web search, 2026-09-18).

## Verdicts

**Standards: Needs work** -- the FS backbone is present, but the dependency
guards return `-1` instead of stopping, four Suggests packages are used
without a guard, nothing is gated at `verbose = 0`, `plot.save` is ignored,
and the documentation misdescribes several arguments.

**Spec: Rework** -- the haplotype table, the function's only return value, is
wrong: on chromosome 1 of the example one adjacent SNP pair reaches r2 0.72
and the function reports a single haplotype of 18 SNPs spanning the whole
chromosome. The default all-chromosome mode crashes on the documented
dataset, and with more than one chromosome the heterozygosity and SNP
position tracks are drawn from the wrong loci. The haplotype step needs to be
rebuilt from the LD matrix rather than from the rotated image.

What works well: the LD values themselves are computed by `snpStats::ld` on
correctly ordered, MAF-filtered, position-deduplicated loci, and the plotted
triangle geometry is consistent with them.

## Findings

**F1 [BLOCKER, confidence: high] -- haplotype boundaries do not follow the LD
threshold (Spec axis: result correctness; DOC5 (proposed rule))**
`R/gl.ld.haplotype.r:449-476`. Haplotypes are read off the first row of a
45-degree rotated copy of the LD matrix (`first_row_2`) whose columns do not
map one to one onto adjacent SNP pairs (30 columns for 21 pairs on chromosome
1). Runs above the threshold are then converted to start and end indices with
`start_haplo <- c(1, start_haplo)` prepended unconditionally, so whenever the
first pair is below the threshold the starts and ends are paired off by one,
the last start is discarded (`:474-476`), and the reported block begins at
SNP 1.
Failure scenario: chromosome 1, `haplo_id = TRUE, min_snps = 3,
ld_threshold_haplo = 0.5`. Adjacent-pair r2 recomputed from the function's
own `ld_snps` matrix: 21 pairs, one at or above 0.5 (0.72, SNPs 165892285 and
171211062), longest run 1 pair. Correct answer under the documented rule:
no haplotype of 3 or more SNPs. Reported: haplotype 1, start 10300000, end
166000000, 18 SNPs, i.e. the whole chromosome, with the one pair actually in
LD placed at its end.
Proposed change: identify blocks directly from the superdiagonal of the LD
matrix (runs of adjacent pairs with LD >= `ld_threshold_haplo`, block size =
run length + 1, kept when >= `min_snps`) and take start and end from
`ld_map_loci$loc_bp`; derive `start_ld_plot` / `end_ld_plot` from the SNP
index afterwards. The rotated image stays for drawing only.

**F2 [HIGH, confidence: high] -- `plot.save` is ignored; a PDF is written on
every call, and never when haplotypes are found (DOC5 (proposed rule), PLT2
(proposed rule))**
`R/gl.ld.haplotype.r:848-874`. In the no-haplotype branch `ggsave` runs
unconditionally (`:850`), then again inside `if (!is.null(plot.save))`
(`:859`), which is always true for a logical. The haplotype branch
(`:599-707`) has no save at all.
Failure scenario: `gl.ld.haplotype(x, chrom_name = chr1, plot.out = FALSE,
plot.save = FALSE)` writes `TENTERFIELD_<chr>.pdf` to `plot.dir` (observed);
`gl.ld.haplotype(x, chrom_name = chr1, haplo_id = TRUE, plot.save = TRUE)`
writes nothing (observed). The second save re-saves `last_plot()`, which was
verified identical to `p_temp`, so the file content is right; the file is just
written twice.
Proposed change: one save per plot, in both branches, executed only when
`plot.save` is `TRUE`.

**F3 [HIGH, confidence: high] -- the default all-chromosome mode crashes on
chromosomes with 2 or 3 informative SNPs (FS5)**
`R/gl.ld.haplotype.r:298-300` guards only `nrow(ld_map_loci) <= 1`. With 2
loci the rotation/`second_row` code fails with "incorrect number of
dimensions"; with 3 it fails with "replacement has 2 rows, data has 0"
(scratch runs on 2..8 loci: 2 and 3 crash, 4 and above run).
Failure scenario: `gl.ld.haplotype(x)` on the documented example
(`chrom_name = NULL`, the documented "all chromosomes" mode) works through
23 chromosomes and aborts on `NC_041750.1_chromosome_X2` (5 loci, 2 after the
MAF filter), losing all plots already produced and returning nothing. An
unknown `chrom_name` ("nope") is not detected: `gl.keep.loc` warns and returns
the object unchanged, so the call silently runs all chromosomes and hits the
same crash.
Proposed change: skip a chromosome with fewer than 4 informative SNPs with a
warning naming it; validate `chrom_name` and `pop_name` against the object
and `stop(error())` when none match.

**F4 [HIGH, confidence: high] -- heterozygosity track, SNP position track,
plot title and x-axis come from all loci of the population, not from the
chromosome being drawn (Spec axis: plot content; STY2)**
`R/gl.ld.haplotype.r:436-446, 669, 726, 754-774`. `utils.basic.stats(pop_ld)`
and `pop_ld$position` are evaluated on the whole population object inside the
chromosome loop.
Failure scenario: `chrom_name = c(chromosome_1, chromosome_2)`. The
chromosome 2 plot draws 40 heterozygosity points and 40 position ticks for
its 18 loci, is titled "40 SNPs", and its axis runs to 179 Mbp (the maximum
of chromosome 1) while chromosome 2 ends at 165 Mbp (observed). The
single-chromosome example masks this. Even with one chromosome the tracks
include loci removed as duplicated positions (`:293-297`), so the track has
more points than the triangle has columns.
Proposed change: subset `pop_ld` to the chromosome's retained loci before
building the tracks, title and axis; compute `utils.basic.stats` once per
population outside the chromosome loop.

**F5 [HIGH, confidence: high] -- PLINK intermediates are written to one
directory and read from another (FS7)**
`R/gl.ld.haplotype.r:254-268`. `gl2plink` is called without `outpath`, so it
writes to `gl.check.wd(NULL)`, which is `options()$dartR_wd` when
`gl.set.wd()` has been used; `utils.read.ped` always reads from `tempdir()`.
Failure scenario: fresh session, `gl.set.wd("<dir>")`, then the documented
example: "cannot open the connection ... gl_plink_TENTERFIELD.ped" (observed).
In a session where an earlier call left `gl_plink_TENTERFIELD.ped` in
`tempdir()`, the stale file is read silently and the LD is computed on the
earlier subset (observed: the same call "ran ok" after a previous run).
Proposed change: write and read from the same per-call path
(`tempfile("gl_plink_")` passed as `outpath`/`outfile` and to
`utils.read.ped`), and remove the files afterwards.

**F6 [MEDIUM, confidence: high] -- dependency guards return -1, guard a
package that is not used, and miss four that are (DEP1)**
`R/gl.ld.haplotype.r:143-172`. `snpStats`, `fields` and `sp` are checked with
`cat(error()); return(-1)`. `fields` is not called anywhere in the function
(`rotate.matrix` is local to `utils.ld.r`). `raster` (`:331-340`), `zoo`
(`:376-430`), `scales` (`:391-446`) and `viridis` (`:640, 722`) are Suggests
and are called unguarded.
Failure scenario: without `raster`, the call writes the PLINK files, computes
LD and then fails with "there is no package called 'raster'"; without
`snpStats` at `verbose = 0` the call returns `-1` and code that expects a data
frame fails later.
Proposed change: `stop(error())` guards for `snpStats`, `sp`, `raster`,
`zoo`, `scales`, `viridis`; drop the `fields` guard.

**F7 [MEDIUM, confidence: high] -- nothing is gated at `verbose = 0` (VRB3,
VRB4 (proposed rule))**
`R/gl.ld.haplotype.r:195, 230, 308, 503, 709, 886`. The coordinates note, the
population-skip warning, the `ld_max_pairwise` cap warning, the `min_snps`
warning, the no-haplotype warning and the final `print(haplo_table)` are
unconditional; `raster::rasterToPolygons` also emits "Regions defined for each
Polygons" on every chromosome.
Failure scenario: the documented example at `verbose = 0` prints 7 lines
(observed).
Proposed change: warnings that change the result (population skipped,
`ld_max_pairwise` raised) at `verbose >= 1`; other warnings at `>= 2`; the
table at `>= 3`; wrap the raster call in `suppressMessages()`.

**F8 [MEDIUM, confidence: high] -- "No haplotypes were identified" is printed
whenever `haplo_id = FALSE` (DOC5 (proposed rule))**
`R/gl.ld.haplotype.r:449-450, 709`. The `else` branch serves both "not asked
to identify" and "none found".
Failure scenario: every default call (`haplo_id = FALSE`) warns that no
haplotypes were identified on each chromosome, although identification was
never attempted (observed).
Proposed change: warn only when `haplo_id` is `TRUE`.

**F9 [MEDIUM, confidence: high] -- pairs with LD exactly 0 or negative are
dropped from the heatmap (Spec axis: plot content; same class as
gl.report.ld.map F1)**
`R/gl.ld.haplotype.r:327`: `ld_columns_2[-ld_columns_2$Freq < 0, ]` keeps
rows with `Freq > 0`. The sparse matrix from `snpStats::ld` holds 0 both for
pairs outside the depth band and for computed pairs with LD 0.
Failure scenario: chromosome 1, `R.squared`: 231 pairs computed within the
depth band, 136 of them exactly 0, 95 cells drawn; the zero pairs are holes,
not "no LD" cells. With `ld_stat = "R"` the 44 negative pairs are dropped as
well (observed).
Proposed change: select cells by band index (`abs(i - j) <= depth`) rather
than by value, so 0 and negative values are drawn.

**F10 [MEDIUM, confidence: high] -- reported haplotype coordinates are cut()
labels with three significant digits (Spec axis: result correctness)**
`R/gl.ld.haplotype.r:496-515`. Start and end are parsed back from the
interval labels of `cut()` (`dig.lab = 3` by default); `midpoint` and
`labels` inherit them.
Failure scenario: the SNP at 165892285 bp is reported as `end = 166000000`
(observed); above 100 Mbp the error reaches 0.5 Mbp, and two block edges within
the same 3-digit rounding produce an identical label.
Proposed change: carry the bp coordinates numerically (falls out of the F1
rebuild).

**F11 [LOW, confidence: high] -- `ind.limit` skips a population that has
exactly `ind.limit` individuals (DOC5 (proposed rule))**
`R/gl.ld.haplotype.r:229`: `nInd(pop_ld) <= ind.limit`. The roxygen text says
"Minimum number of individuals that a population should contain".
Failure scenario: 10 individuals, `ind.limit = 10`: skipped with the message
"less than 10 individuals" (observed).
Proposed change: use `<`. Consequence: populations with exactly `ind.limit`
individuals are now analysed.

**F12 [LOW, confidence: high] -- an object without chromosome information,
or a SilicoDArT object, fails with an internal error (FS5, DAT5)**
`R/gl.ld.haplotype.r:139, 214-223`. `utils.check.datatype` accepts
SilicoDArT; an empty `@chromosome` gives `chr_list` of length 0 and
`names(x) <- paste0("chr_", chr_list)` fails.
Failure scenario: `gl.ld.haplotype(testset.gl)` and
`gl.ld.haplotype(testset.gs)`: "'names' attribute [1] must be the same length
as the vector [0]" (observed).
Proposed change: `stop(error())` naming the `@chromosome`/`@position` slots
when they are empty; restrict the datatype check to SNP.

**F13 [LOW, confidence: high] -- documentation does not match behaviour
(DOC2, DOC5 (proposed rule), DOC7 (proposed rule))**
`R/gl.ld.haplotype.r:1-96`. `verbose` text is not the standard wording;
`plot.dir` says "[default = working directory]" but resolves to `tempdir()`
unless `gl.set.wd()` was used; `snp_pos`, `target.snp1-3` and `col.*` take
effect only when no haplotypes are drawn, which is not stated;
`@description` says the function identifies haplotypes but the default
`haplo_id = FALSE` never does; typos "Nme", "analised"; `@author` has a
custodian but no `Author(s):` line. `utils.flag.start(build = "Jody")` is the
outdated proforma argument.
Proposed change: docs only, plus drop `build =`.

**F14 [LOW, confidence: high] -- dead code, a misspelt argument and a
partial-match typo (STY1, STY3)**
`R/gl.ld.haplotype.r:216-223` builds a list `p` that is never used;
`mean_column`, `df_col`, `second_row_ver_2`, `real_distance_4`, `test_var_2`,
`enlarge_factor`, `correction_factor`, `reduce_factor` are computed and
unused; `fortify(polygon_haplo, polygon_haplo = "id")` (`:347`) passes a
non-existent argument and triggers "Arguments in `...` must be used" on every
call; `haplo_temp_a$end_ld_plo` (`:646`) works only through `$` partial
matching; the `pop_name` argument is overwritten in the loop (`:227`).
Failure scenario: two warnings on every call; a future ggplot2 that errors on
unused `...` breaks the function.
Proposed change: remove the unused objects, use `region = "id"` or drop the
argument, fix the typo, rename the loop variable.

**F15 [INFO, confidence: high] -- the plotting path rests on deprecated
spatial tooling (DEP2)**
`fortify(<SpatialPolygonsDataFrame>)` has been deprecated since ggplot2 3.4.4
and warns on every call; `raster` and `sp` are in maintenance mode and the
package already imports `terra`. No change proposed in this review; a
migration to `terra`/`sf` (or drawing the triangle directly with
`geom_polygon` from computed coordinates) is a follow-up.

## Proposed changes

1. Rebuild haplotype identification from the adjacent-pair LD of the LD
   matrix and report bp coordinates numerically (F1, F10).
   **Consequence: numerical output changes -- the haplotype table for any
   call with `haplo_id = TRUE` changes; on the example it goes from one
   whole-chromosome block to no block.**
2. Save the plot once, in both branches, only when `plot.save = TRUE` (F2).
   **Consequence: user-visible behaviour change -- default calls stop
   writing `<pop>_<chr>.pdf` to `plot.dir`; calls with `plot.save = TRUE`
   start writing it when haplotypes are drawn.**
3. Skip chromosomes with fewer than 4 informative SNPs with a warning;
   error on `chrom_name` / `pop_name` values absent from the object (F3).
4. Build heterozygosity and position tracks, title and axis from the
   chromosome's retained loci; compute `utils.basic.stats` once per
   population (F4).
5. Write and read the PLINK intermediates from one per-call temporary path
   and delete them afterwards (F5).
6. Replace the `return(-1)` guards with `stop(error())` for `snpStats`, `sp`,
   `raster`, `zoo`, `scales`, `viridis`; drop the `fields` guard (F6).
7. Gate all messages by `verbose` (result-affecting warnings at >= 1, others
   at >= 2, table at >= 3), suppress the raster message, and warn about
   missing haplotypes only when `haplo_id = TRUE` (F7, F8).
8. Select heatmap cells by depth band instead of by non-zero value so LD of
   0 and negative LD are drawn (F9).
9. Use `nInd(pop_ld) < ind.limit` (F11).
   **Consequence: populations with exactly `ind.limit` individuals are now
   analysed.**
10. Clear `stop(error())` for empty `@chromosome`/`@position` and for
    SilicoDArT input (F12).
11. Documentation: standard `verbose` text, correct `plot.dir` default,
    scope of `snp_pos`/`target.snp*`/`col.*`, `haplo_id` default stated in
    the description, typos, `Author(s):` line; drop `build =` (F13). Docs
    only.
12. Remove dead code, fix the `fortify` argument and the `end_ld_plo` typo,
    rename the shadowed loop variable (F14).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY -- run
- Spec: behaviour vs roxygen on the documented platypus example -- run
- Independent computation: adjacent-pair r2 taken from the function's own
  `snpStats::ld` matrix (superdiagonal) and compared with the reported
  blocks -- run (F1); zero/negative pair counts -- run (F9)
- Multi-chromosome and few-SNP paths: scratch-copy instrumentation -- run
  (F3, F4)
- `gl.set.wd()` interaction in a fresh session -- run (F5)
- SilicoDArT and no-chromosome input -- run (F12)
- FBM path (DAT6): SKIPPED -- no FBM fixture; the function densifies through
  `gl2plink` in any case
- Plot rendering correctness beyond the data fed to the layers: SKIPPED --
  visual inspection not automated; the pinned baseline covers the returned
  table and files, not the image
- `target.snp1-3` nearest-SNP matching: not exercised beyond code reading
- Downstream callers: dartr2shiny `shiny_fun/Fun_gl.ld.haplotype.R` calls
  the function with `pop_name`, `chrom_name`, `ld_max_pairwise`, `maf`,
  `ld_stat`, `ind.limit`, `haplo_id`, `min_snps`, `ld_threshold_haplo`,
  `plot_het`, `snp_pos`, `col.all`, `color_haplo`, `color_het` and captures
  the printed plot; it does not pass `plot.save`, `plot.dir` or `plot.out`.
  No sibling `dartR.*` package calls the function.

## Approval

All twelve changes approved by Luis on 2026-09-18 (AskUserQuestion selections;
changes 1, 2 and 9 approved with their stated consequences).

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | numerical output change accepted |
| 2 | approved | Luis | file-writing change accepted |
| 3 | approved | Luis |  |
| 4 | approved | Luis |  |
| 5 | approved | Luis |  |
| 6 | approved | Luis |  |
| 7 | approved | Luis |  |
| 8 | approved | Luis |  |
| 9 | approved | Luis | ind.limit boundary change accepted |
| 10 | approved | Luis |  |
| 11 | approved | Luis |  |
| 12 | approved | Luis |  |

## Outcome

Applied on dartR.popgen `dev_luis` (all twelve changes). Evidence:

- Characterization test `tests/testthat/test-gl.ld.haplotype.R`: the pre-change
  baseline (29 assertions) run against the revised function moved exactly the
  pins tagged in the file: PDF and PLINK files no longer present (changes 2, 5),
  `haplo_id = TRUE` table (change 1), `verbose = 0` output (change 7), the
  `ind.limit` boundary (change 9) and the error text for `testset.gl`
  (change 10). No other pin moved. The updated file (48 assertions) passes.
- Change 1, independent check: a `SnpMatrix` built directly from the genlight
  dosages of the 22 retained chromosome-1 SNPs, `snpStats::ld(depth = 1)`,
  gives adjacent-pair r2 identical to the function's values (max difference
  2.5e-16). At threshold 0.3 the four blocks (34408296-39659332,
  55819584-57787059, 104075140-111265950, 165892285-171211062) match the
  function's table exactly; at 0.5 with `min_snps = 2` the single block is
  165892285-171211062 (previously 10300000-166000000 with `min_snps = 3`);
  at 0.73 none; at 0 one block from the first to the last SNP.
- Change 2: `plot.save = FALSE` writes nothing in either branch;
  `plot.save = TRUE` writes `TENTERFIELD_<chr>.pdf` in both.
- Change 3: `gl.ld.haplotype(x)` (all chromosomes) completes on the example,
  skipping 9 chromosomes with fewer than 4 SNPs by name, and returns the
  0 x 10 table; `chrom_name = "nope"` and `pop_name = "nope"` stop with the
  name in the message.
- Change 4 (scratch instrumentation): with `chrom_name = c(chr1, chr2)` the
  chromosome 2 plot has 18 heterozygosity points and 18 position ticks for
  18 loci, title "18 SNPs", axis maximum 165010140 (its own last SNP);
  chromosome 1: 22/22/22, 179168634.
- Change 5: with `options(dartR_wd = <dir>)` in a fresh session the call
  runs, the PDF lands in `plot.dir`, and no `gl_plink_*` file remains in
  either directory.
- Change 7: `verbose = 0` prints 0 lines; `verbose = 3` prints the table or
  "No haplotypes in the results table".
- Change 8: chromosome 1, depth 5: 190 band cells expected, 190 drawn, 95 of
  them with r2 = 0 (previously 95 drawn in total); `ld_stat = "R"`: 44
  negative cells drawn (previously 0).
- Change 9: 10 individuals with `ind.limit = 10` are analysed; 9 are skipped.
- `devtools::document()` regenerated `man/gl.ld.haplotype.Rd`; NAMESPACE
  unchanged. `zoo` is no longer used by the function (the rotated-matrix
  code it served was removed with change 1), so its guard was not added.
- `devtools::check()`: 0 errors; 2 warnings and 4 notes pre-existing and
  unrelated (git-ignored binaries and scratch files in the working tree,
  `gl.plot.popcluster` globals, NEWS parsing).
- NEWS entry added under "Bug fixes".
- PR: green-striped-gecko/dartR.popgen#91, from branch `review-gl.ld.haplotype`
  (cut from `origin/dev`), commit 52eced6. A parallel session's commit on
  `dev_luis` (8eea717, gl.blast) had swept up the staged files; with Luis's
  approval it was rewritten to a4fd6ed (gl.blast only, PR #90) and the
  gl.ld.haplotype work moved to its own branch.

```json
{
  "function": "gl.ld.haplotype",
  "package": "dartR.popgen",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "31fb40c",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "BLOCKER", "confidence": "high", "rule": "spec-correctness/DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5/PLT2", "status": "approved", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "FS5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "HIGH", "confidence": "high", "rule": "spec-plot-content/STY2", "status": "approved", "change": 4},
    {"id": "F5", "severity": "HIGH", "confidence": "high", "rule": "FS7", "status": "approved", "change": 5},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "DEP1", "status": "approved", "change": 6},
    {"id": "F7", "severity": "MEDIUM", "confidence": "high", "rule": "VRB3/VRB4", "status": "approved", "change": 7},
    {"id": "F8", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 7},
    {"id": "F9", "severity": "MEDIUM", "confidence": "high", "rule": "spec-plot-content", "status": "approved", "change": 8},
    {"id": "F10", "severity": "MEDIUM", "confidence": "high", "rule": "spec-correctness", "status": "approved", "change": 1},
    {"id": "F11", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 9},
    {"id": "F12", "severity": "LOW", "confidence": "high", "rule": "FS5/DAT5", "status": "approved", "change": 10},
    {"id": "F13", "severity": "LOW", "confidence": "high", "rule": "DOC2/DOC5/DOC7", "status": "approved", "change": 11},
    {"id": "F14", "severity": "LOW", "confidence": "high", "rule": "STY1/STY3", "status": "approved", "change": 12},
    {"id": "F15", "severity": "INFO", "confidence": "high", "rule": "DEP2", "status": "no-change", "change": null}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "plot image inspection not automated", "target.snp matching not exercised"],
  "status": "pr-open",
  "pr": 91
}
```
