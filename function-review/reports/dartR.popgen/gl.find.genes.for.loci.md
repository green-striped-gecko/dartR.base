# Review: gl.find.genes.for.loci (dartR.popgen)

- Family mode: analysis
- Date: 2026-09-04
- Reviewer: Claude (claude-fable-5-1), dartr-function-review v2.0.0
- Package commit: 34546bd (dev_luis, reviewed state before the fix)
- Datasets: testset.gl (dartR.data) with synthetic chromosome/position slots
  and a synthetic four-gene GFF3 written by the test; taro DArTseq genets
  (392 x 16,910, project `brad`) mapped to the Ahi taro assembly
  GCA_009445465.1 with its NCBI GFF (56,238 genes) as the real-data check
- Baseline: tests/testthat/test-gl.find.genes.for.loci.R. The tests
  directory held only test-gl.nhybrids.R. The new file asserts the
  documented behaviour rather than snapshotting the defective output; the
  pre-fix run is recorded below as the baseline evidence.

**Standards: Needs work** — structure follows the house order, but the
roxygen header has no author/custodian, a non-standard verbose text and an
example that cannot run, and the end-of-function message is missing.
**Spec: Needs work** — before the fix the nearest-gene branch, half of what
the function promises, returned no usable rows; with the fix applied the
remaining gaps are the unmet claims about `product` parsing and tie-breaks,
and silent NA output when sequence names do not match.

What works well: the overlap branch (`foverlaps` with the documented
midpoint/length/id tie-break) is correct, and the GFF attribute parsing
handles the NCBI `ID=gene-…;Name=…` style.

## Findings

**F1 [BLOCKER, confidence: high] — nearest-gene rows lose their locus
identifiers and collapse to one row (DOC5, DAT2 principle: results must
track the loci they describe)**
`R/gl.find.genes.for.loci.r:259-279` (pre-fix) — the roll join
`gene_pts[points, roll = "nearest"]` returns the join columns under the
names of `gene_pts` (`gene_seqid`, `join_coord`) and the other columns of
`points` under their own names (`pos`, `locus`). The block then reads
`i.locus`, `i.chrom` and `i.pos`, which do not exist as columns, so R falls
back to the NSE placeholders assigned at the top of the function
(`i.pos <- i.locus <- i.chrom <- NA`). Every non-overlapping locus therefore
gets `locus = NA`, `distance_bp = NA`, `nearest_side = NA`, and the final
`.SD[1L], by = locus` keeps one NA-locus row for all of them. A second
defect in the same block: `abs(i.pos - join_coord)` would be 0 even with
the identifiers fixed, because `join_coord` carries the locus position
after the join, not the gene edge.
Failure scenario: eight taro candidate SNPs, five outside genes. Output has
four rows: the three overlapping loci and one row with an empty locus
naming the gene nearest to whichever locus sorted first. The synthetic
test (loci at 300, 1000 and 1450 bp on chr1) gives `gene_name = NA` and
`distance_bp = NA` for all three.
Proposed change: build the output from the columns the join actually
returns (`locus`, `gene_seqid`, `pos`), take the distance from the gene
interval (`pmin(abs(pos - gene_start), abs(pos - gene_end))`), and give a
locus on a sequence without genes NA gene columns rather than an error.
Applied (approved by the custodian's request to fix the bug).

**F2 [MEDIUM, confidence: high] — silent NA output when locus sequence
names never match the GFF, and silent dropping of unmapped loci (VRB4,
proposed rule; DOC5)**
`R/gl.find.genes.for.loci.r:188-208` — loci with NA chromosome or position
are dropped without a message unless all of them are, and a chromosome
naming mismatch (`chr1` in the genlight, `NC_012345.1` in the GFF) yields a
table of NA genes with no hint of the cause. The `@return` text promises
"one row per input locus".
Failure scenario: a user whose genlight carries DArT chromosome labels
runs the function against an NCBI GFF; every row comes back NA and nothing
says why.
Proposed change: at `verbose >= 1`, report how many requested loci had no
position and how many locus sequence names are absent from the GFF, with an
example of each; amend `@return` to "one row per locus with a position".

**F3 [MEDIUM, confidence: high] — `gene_product` is empty for standard
GFF3 (DOC5)**
`R/gl.find.genes.for.loci.r:147, 180` — `product` is read from the gene
feature's own attributes. In NCBI and most Ensembl GFF3 files `product`
sits on the mRNA or CDS children, so `gene_product` is NA for every gene
even though `@details` says the function parses `product` "to provide
informative gene labels".
Failure scenario: the Ahi taro GFF (NCBI): 56,238 genes, `gene_product` NA
for all of them; the user gets `Taro_020035` and no description.
Proposed change: after selecting gene features, join the first non-NA
`product` of child features (`Parent == ID`, one level) into
`gene_product`; keep the gene's own value when present.

**F4 [LOW, confidence: high] — `nearest_side` is in coordinate space only;
gene strand is not reported (principle: output interpretable by the
biologist)**
`R/gl.find.genes.for.loci.r:263` — "left"/"right" says where the locus is
on the sequence, but the question a user asks is whether the SNP is
upstream (promoter side) or downstream of the gene, which needs the strand.
The GFF strand column is read and discarded.
Failure scenario: a SNP 884 bp "right" of a minus-strand gene is actually
884 bp upstream of its start; the user has to open the GFF to learn this.
Proposed change: add `gene_strand` and `relative_position`
("inside", "upstream", "downstream", NA when strand is unknown) as two
extra columns; existing columns unchanged.

**F5 [LOW, confidence: medium] — the documented tie-break is not applied in
the nearest branch (DOC5)**
`R/gl.find.genes.for.loci.r:259` — `@details` promises "closest to gene
midpoint, then shorter gene length, then lexicographic gene_id" for equally
close genes. In the nearest branch `roll = "nearest"` returns one edge per
locus by data.table's own rule, so a locus exactly midway between two
genes gets whichever edge data.table prefers; the final de-duplication
never sees the alternative.
Failure scenario: locus at 350 with gene A ending at 300 and gene B
starting at 400: the reported gene depends on data.table's tie handling,
not on the documented rule.
Proposed change: roll backward and forward (`roll = Inf` and `roll = -Inf`)
to get the nearest edge on each side, then choose by distance and the
documented tie-break. **Consequence: the reported gene can change for a
locus at exactly equal distance from two genes.**

**F6 [LOW, confidence: high] — roxygen gaps (DOC1, DOC2, DOC3, DOC7
proposed)**
`R/gl.find.genes.for.loci.r:1-56` — no `@author` tag at all (DOC1; DOC7
wants `Author(s):` and `Custodian:`); the `verbose` text is not the
standard house text (DOC2); the example refers to `species.gff3`, which
does not exist, and is wrapped in `\dontrun{}` (DOC3); `@param gff.file`
says "plain or with a .gz alongside" but a `.gz` path passed directly also
works (ape::read.gff reads gzip through R connections), which is worth
saying; `@return` says "one row per input locus" (see F2).
Proposed change: add the author/custodian block, the DOC2 verbose text, a
runnable example built on a small GFF shipped in `inst/extdata` (or
written with `writeLines` in the example) with `testset.gl` given
synthetic positions, and the two wording fixes.

**F7 [LOW, confidence: high] — no end-of-function message (FS9)**
`R/gl.find.genes.for.loci.r:316-318` — the function flags its start with
`utils.flag.start` but never prints `Completed: gl.find.genes.for.loci`.
Failure scenario: at `verbose = 1` the user sees "Starting" and nothing
else.
Proposed change: add `if (verbose > 0) cat(report("Completed:", funname, "\n"))`
before the return.

**F8 [LOW, confidence: high] — `loci` is required although "all loci" is
the common case (API1 proposed: additive default, no change for existing
callers)**
`R/gl.find.genes.for.loci.r:59, 90` — the user must pass `locNames(x)` to
annotate every locus.
Proposed change: default `loci = NULL` meaning all loci with a position.

**F9 [INFO, confidence: high] — redundant dependency guards and legacy
save idiom (DEP1, PLT2 proposed)**
`R/gl.find.genes.for.loci.r:81-85, 307-314` — `ape`, `data.table` and
`stringr` are in Imports, so the `requireNamespace` loop can never fail;
`save2tmp` follows the tempdir idiom the campaign is retiring.
Proposed change: drop the loop; leave `save2tmp` until the team ratifies
PLT2 (no change proposed now).

## Proposed changes

1. Fix the nearest-gene branch: keep locus identifiers, take the distance
   from the gene interval, NA gene columns for a locus on a sequence
   without genes (F1). **Consequence: results change for every locus that
   does not overlap a gene — from an empty collapsed row to one row per
   locus.** Applied.
2. Report unmapped loci and sequence-name mismatches at `verbose >= 1`;
   amend `@return` (F2).
3. Fill `gene_product` from child mRNA/CDS features when the gene line has
   none (F3).
4. Add `gene_strand` and `relative_position` output columns (F4).
   **Consequence: two new columns; existing columns unchanged.**
5. Apply the documented tie-break in the nearest branch (F5).
   **Consequence: the reported gene can change for a locus at exactly
   equal distance from two genes.**
6. Roxygen: author/custodian, standard verbose text, runnable example,
   `.gz` note, `@return` wording (F6). Run `devtools::document()`.
7. Add the FS9 completion message (F7).
8. Default `loci` to all loci with a position (F8). **Consequence: none for
   existing callers; a new default.**
9. Remove the redundant `requireNamespace` loop (F9).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: behaviour vs roxygen on synthetic GFF and on the Ahi taro GFF — run
- Real-data check: eight taro SNPs on the Ahi assembly, compared against an
  independent nearest-gene computation from the same GFF — run (matches
  after the fix)
- FBM path (DAT6): SKIPPED — the function reads only `chromosome`,
  `position` and `locNames`; no genotype access
- SilicoDArT dispatch: SKIPPED — `utils.check.datatype(accept = "SNP")`
  rejects it by design
- Google Group / GitHub issues search: SKIPPED — the function is new in the
  dev branch; no user reports expected

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | requested: "please fix the bug you found" |
| 2 | approved | Luis | |
| 3 | approved | Luis | |
| 4 | approved | Luis | |
| 5 | approved | Luis | recommended option |
| 6 | approved | Luis | |
| 7 | approved | Luis | |
| 8 | approved | Luis | |
| 9 | approved | Luis | |

## Outcome

- Changes 1 to 9 applied in `R/gl.find.genes.for.loci.r`; `man/gl.find.genes.for.loci.Rd`
  regenerated with `devtools::document()` (NAMESPACE unchanged); NEWS
  entries added for the bug fix and the new columns and default.
- Evidence, before (change 1): `tests/testthat/test-gl.find.genes.for.loci.R`
  at commit 34546bd fails 7 expectations in the first test (rows 2 to 5 have
  `locus`, `gene_name`, `distance_bp` and `nearest_side` all NA; 2 rows
  returned instead of 5).
- Evidence, after: 25 of 25 expectations pass across 8 tests (rows per
  locus, distances and sides; product from mRNA/CDS children; nested-gene
  tie-break; equal-distance tie-break; `loci = NULL` default; verbose 1
  messages for unmapped loci and unknown sequence names; `.gz` path;
  unknown locus names). `devtools::check_man()`: no issues. The roxygen
  example runs.
- Real data: the eight taro candidates on the Ahi assembly return 884, 0,
  610, 0, 0, 231, NA and 12,211 bp, identical to an independent
  computation from the same GFF; strand-aware positions upstream, inside,
  downstream, inside, inside, downstream, NA, downstream; `gene_product`
  now "hypothetical protein" from the mRNA rows instead of NA. All 16,910
  loci with `loci = NULL`: 14,074 rows (7,012 inside a gene, 6,150 near,
  912 on scaffolds without genes; 2,836 without a position reported at
  verbose 1) in 12.8 s.
- Not verified: `R CMD check` for the whole package (not run; the change is
  confined to one function, its Rd and its test).
- PR: green-striped-gecko/dartR.popgen #87 (dev_luis -> dev, commit 5fdf605).

```json
{
  "function": "gl.find.genes.for.loci",
  "package": "dartR.popgen",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "34546bd",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "BLOCKER", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "VRB4", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "principle: interpretable output", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "medium", "rule": "DOC5", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "FS9", "status": "approved", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "API1", "status": "approved", "change": 8},
    {"id": "F9", "severity": "INFO", "confidence": "high", "rule": "DEP1", "status": "approved", "change": 9}
  ],
  "coverage_skipped": ["DAT6: no genotype access", "SilicoDArT: rejected by design", "issues search: new function"],
  "status": "pr-open",
  "pr": 87
}
```
