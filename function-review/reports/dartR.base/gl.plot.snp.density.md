# Review: gl.plot.snp.density (dartR.base)

## Provenance

- Model: Claude Fable 5.1 (claude-fable-5-1)
- Skill: dartr-function-review v2.0.0 (Phase A, read-only review)
- Package commit: 562befa (dev_luis after merging origin/dev fadab59;
  `R/gl.plot.snp.density.r` identical on both)
- Date: 2026-09-15
- Family mode: graphics (pure plotting; report-mode checks applied: input
  untouched, no history append)
- Datasets: platypus.gl with `@chromosome`/`@position` set from
  `Chrom_Platypus_Chrom_NCBIv1`/`ChromPos_Platypus_Chrom_NCBIv1` (104
  chromosomes, 921 positioned SNPs); testset.gl (no chromosome data);
  testset.gs (SilicoDArT); gl.gen2fbm(platypus)
- Baseline: tests/testthat/test-gl.plot.snp.density.R (new file; 11
  tests, 33 expectations, all passing at 562befa; defects pinned as-is and
  tagged with finding IDs; binned counts checked against an independent
  `tapply`/`table` computation: 25 chromosomes, 888 SNPs, 331 bins agree)
- Callers: none in dartR.base or sibling packages; dartr2shiny
  (`shiny_fun/Fun_gl.plot.snp.density.R`) passes `x`, `bin.size`,
  `min.snps`, `min.length`, `chr.info`, `plot.theme`, `color.palette` by
  name.
- Checks skipped: dartR Google Group not searched (no browser session);
  GitHub issue search across the org found nothing for this function.

## Verdicts

**Standards: Needs work** -- entry checks, dependency guards, verbosity
gating and the invisible ggplot return conform, and an FBM-backed object
plots without densification, but `verbose = 3` prints nothing beyond
`verbose = 2`, two argument-check messages state the wrong bound, the
plot cannot be suppressed and saving still uses the `save2tmp`/`tempdir()`
idiom, and the header lacks `@family` and the standard verbose text.

**Spec: Needs work** -- the binned counts are correct, but two claims in
`@details` are false: chromosomes are ordered alphabetically, not longest
to shortest, and bins with no SNPs are not drawn at all rather than in
the lowest palette colour; an object without chromosome data fails with
an R-internal message.

What works well: SNP counts per bin match an independent computation
exactly; the input object comes back identical; `verbose = 0` is silent;
the SilicoDArT gate works; the example runs in 0.07 s.

## Findings

**F1 [MEDIUM, confidence: high] -- documented chromosome order is not the
one drawn (DOC5 proposed rule; STY1)**
`R/gl.plot.snp.density.r:36-37, 134, 161-163` -- `@details` says
"ordered from longest (bottom) to shortest (top)". `chr_stats` is
arranged by size, but the factor levels are then set to
`sort(unique(chr_label), decreasing = TRUE)`, so the plot is alphabetical
and the size ordering is dead code. Commit d41107a8 ("order
alphabetically") is newer than the documentation, so the code is the
intent and the documentation is wrong. A plain string sort also
interleaves unpadded numbers: `chr1, chr10, chr11, chr2`.
Failure scenario: a user reading the help expects size order and reads
the top row as the shortest chromosome; on a reference with names
`chr1..chr22`, `chr10` sits between `chr1` and `chr2`.
Proposed change: keep alphabetical order but numeric-aware
(`stringr::str_sort(numeric = TRUE)`; stringr is already in Imports) so
`chr2` precedes `chr10`; drop the dead `arrange()`; document the order.
**Consequence: the row order changes for chromosome names with unpadded
numbers.**

**F2 [MEDIUM, confidence: high] -- empty bins are not drawn (DOC5
proposed rule; PLT1)**
`R/gl.plot.snp.density.r:37-38, 151-159` -- `@details` says "Bins
containing no SNPs are rendered in the lowest colour of the palette".
`plyr::count()` only produces rows for bins with at least one SNP, so an
empty bin is panel background, indistinguishable from the region past
the chromosome end (baseline test "bins with no SNPs are absent": minimum
`n_snps` is 1 over 331 bins).
Failure scenario: a 20 Mb SNP desert inside a chromosome looks like a
gap in the assembly, and the colour scale starts at 1 rather than 0.
Proposed change: complete every chromosome's bins from 0 to its last
occupied bin with `n_snps = 0` before plotting, so empty bins take the
lowest palette colour and the legend starts at 0.
**Consequence: the plot changes -- empty bins are coloured and the fill
scale starts at 0.**

**F3 [MEDIUM, confidence: high] -- object without chromosome data fails
with an R-internal error (FS5, DAT5)**
`R/gl.plot.snp.density.r:104-108` -- when `@chromosome` is empty (every
packaged dataset as shipped, including testset.gl), `data.frame()` stops
with "arguments imply differing number of rows: 0, 255".
Failure scenario: `gl.plot.snp.density(testset.gl)`.
Proposed change: stop with a message naming the slots to fill
(`x$chromosome`, `x$position`, as in the example) when either is empty
or not of length `nLoc(x)`.

**F4 [LOW, confidence: high] -- verbose 3 adds nothing (VRB1; DOC5
proposed rule)**
`R/gl.plot.snp.density.r:122-138` -- the only progress line is the
retained SNP count; nothing reports how many chromosomes survived
`min.snps`/`min.length` or how many were dropped, at any verbosity.
Failure scenario: 104 platypus scaffolds go in, 25 come out, and the
user has no way to see which filter removed 79 of them without
recomputing.
Proposed change: at `verbose >= 2` report chromosomes kept and dropped
by each criterion; at `verbose >= 3` print the per-chromosome table
(name, SNPs, last position, bins).

**F5 [LOW, confidence: high] -- argument-check messages state the wrong
bound (FS5)**
`R/gl.plot.snp.density.r:96-101` -- `min.snps < 1` stops with "must be
> 1" and `min.length < 1` with "must be > 1 bp"; both accept 1.
Proposed change: ">= 1".

**F6 [LOW, confidence: high] -- plot always printed; saving uses the
tempdir idiom (PLT2 proposed rule; PLT3)**
`R/gl.plot.snp.density.r:195-206` -- `print(p1)` is unconditional and
`save2tmp` writes an RDS to `tempdir()` for `gl.print.reports()`. 34
functions in the package use `plot.display`/`plot.file`/`plot.dir` with
`utils.plot.save()`; 3 still use `save2tmp`.
Failure scenario: in an R Markdown chunk or a loop that assigns the
result, the plot prints once by the function and again by the user.
Proposed change: add `plot.display = TRUE`, `plot.file = NULL`,
`plot.dir = NULL` (resolved by `gl.check.wd()`, saved by
`utils.plot.save()`) and remove `save2tmp`.
**Consequence: a call passing `save2tmp` stops with "unused argument".**
dartr2shiny does not pass it.

**F7 [LOW, confidence: high] -- documentation (DOC1, DOC2, DOC5 proposed
rule, DOC7 proposed rule)**
`R/gl.plot.snp.density.r:1-57`:
- no `@family` tag (DOC1).
- `verbose` text is not the DOC2 standard; `@author` has no `Author(s):`
  part (DOC7).
- "chromosome length (Mb)" in `@description`, `@param min.length` and
  the labels is the position of the last SNP on the chromosome, not the
  assembly length (`chr_size = max(pos)`).
- `@details` claims corrected by F1 and F2; `@param x` should say the
  slots are filled with `x$chromosome <- ...` as in the example.
Proposed change: docs-only rewrite; `devtools::document()`.

**F8 [INFO, confidence: high] -- argument naming and redundant guards
(PLT1, DEP1)**
`color.palette` differs from the house `plot.colors`; dartr2shiny passes
it by name, so a rename is an API2 change; no change proposed. The
`ggplot2` and `dplyr` guards test packages already in Depends; harmless.
Positions equal to 0 are dropped by `pos > 0`, which excludes a
0-based first coordinate; DArT positions are 1-based, so no change
proposed.

## Proposed changes

1. Complete each chromosome's bins with zero counts so empty bins are
   drawn in the lowest palette colour and the fill scale starts at 0
   (F2). **Consequence: the plot changes for every dataset with empty
   bins.**
2. Numeric-aware alphabetical chromosome order; dead size ordering
   removed; order documented (F1). **Consequence: row order changes for
   chromosome names with unpadded numbers.**
3. Clear error when `@chromosome`/`@position` are empty or the wrong
   length (F3).
4. Chromosome filter summary at `verbose >= 2`; per-chromosome table at
   `verbose >= 3` (F4).
5. Correct the `min.snps`/`min.length` messages (F5).
6. `plot.display`/`plot.file`/`plot.dir` with `utils.plot.save()`
   replace `save2tmp` (F6). **Consequence: `save2tmp` is removed from the
   signature.**
7. Documentation rewrite (F7), docs only.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY -- run
- Spec: behaviour vs roxygen on platypus.gl with chromosome data -- run;
  binned counts vs independent computation -- run (agree)
- Family checks (graphics/report mode): input untouched (run, identical);
  no history append (run, none)
- FBM path (DAT6): run on gl.gen2fbm(platypus); plot data identical to
  the in-memory object; only `@chromosome`/`@position` are read
- SilicoDArT: run; rejected by `utils.check.datatype(accept = "SNP")`
- Callers (API3): dartr2shiny read; named arguments only, none of them
  touched by the proposed changes
- Google Group: SKIPPED -- no browser session

## Report notes (other functions, not fixed here)

- Every packaged dataset ships with an empty `@chromosome` slot; the
  example fills it from `loc.metrics`. A helper that sets
  `@chromosome`/`@position` from named `loc.metrics` columns would serve
  this function, `gl2plink` and `gl2hapmap`, which each do it by hand.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | plot change approved as stated |
| 2 | approved | Luis | row-order change approved as stated |
| 3 | approved | Luis |  |
| 4 | approved | Luis |  |
| 5 | approved | Luis |  |
| 6 | approved | Luis | save2tmp removal approved as stated |
| 7 | approved | Luis |  |

## Outcome

Branch `review-gl.plot.snp.density` from origin/dev fadab59. PR #407 (-> dev).

- Change 1 (F2): every bin from the chromosome start to the bin of its
  last SNP is in the plot data with `n_snps = 0` where empty. Evidence:
  platypus at 5 Mb bins has 349 rows (sum over chromosomes of last bin
  index + 1; was 331), minimum count 0 (was 1), occupied bins identical
  to the independent count, per-chromosome bin sequences checked.
  NEWS entry added.
- Change 2 (F1): `rev(stringr::str_sort(chr, numeric = TRUE))` sets the
  levels; the dead `arrange(desc(chr_size))` removed. Evidence: platypus
  order unchanged; a relabelled object with `chr1..chr12` gives levels
  `chr12, chr11, ..., chr1` (was `chr9, chr8, ..., chr2, chr12, chr11,
  chr10, chr1`). NEWS entry added.
- Change 3 (F3): `length(@chromosome)`/`length(@position)` checked
  against `nLoc(x)`. Evidence: testset.gl stops with "x needs a
  chromosome name and a position for every locus; found 0 chromosome
  entries and 255 position entries for 255 loci ...".
- Change 4 (F4): verbose 2 line "Chromosomes with positioned SNPs: 46;
  dropped 21 with fewer than 10 SNPs and 0 with the last SNP below
  2,000,000 bp; retained 25"; verbose 3 prints the per-chromosome table
  top to bottom (chromosome, n_snps, last_snp_Mb, n_bins). Evidence:
  verbose-3 run in the test and end to end.
- Change 5 (F5): messages read ">= 1". Evidence: test "argument checks".
- Change 6 (F6): `plot.display`, `plot.file`, `plot.dir` added;
  `save2tmp` removed; `gl.check.wd()` resolves the directory,
  `utils.plot.save()` writes `<plot.file>.RDS`. Evidence: `plot.display =
  FALSE` writes no page to a png device, TRUE writes one; `dens.RDS` is
  written and reloads as a ggplot. NEWS entry added.
- Change 7 (F7): header rewritten (`@family graphics`, DOC2 verbose text,
  author block, last-SNP wording, new arguments); `devtools::document()`
  run. Adding the family tag regenerates the "Other graphics" cross-links
  in seven sibling Rd files (gl.colors, gl.map.interactive,
  gl.plot.heatmap, gl.select.colors, gl.select.shapes, gl.smearplot,
  gl.tree.nj); they are included in the commit.

Characterization test: 11 tests, 66 expectations, 0 failures. Every
flipped pin carries an `[approved diff, change n]` tag: changes 1-6. No
unexplained diff. Example runs clean via `tools::Rd2ex()` (0.07 s).
Function run end to end at verbose 3 on platypus with chromosome data.

Package check (`R CMD check --no-tests`): no finding from the package
code. The remaining output (non-portable names of local scratch files,
non-standard top-level files, "built under R 4.4.3" install warnings)
comes from the local checkout and the R installation, not from the
package.

```json
{
  "function": "gl.plot.snp.density",
  "package": "dartR.base",
  "family": "graphics",
  "skill_version": "2.0.0",
  "commit": "562befa",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5",
     "status": "approved", "change": 2},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5",
     "status": "approved", "change": 1},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "FS5",
     "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "VRB1",
     "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "FS5",
     "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "PLT2",
     "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC1",
     "status": "approved", "change": 7},
    {"id": "F8", "severity": "INFO", "confidence": "high", "rule": "PLT1",
     "status": "noted", "change": null}
  ],
  "coverage_skipped": ["Google Group: no browser session"],
  "status": "pr-open",
  "pr": 407
}
```
