# Review: gl.find.loci.in.genes (dartR.popgen)

- Family mode: analysis
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: f2c8b8a (origin/dev, reviewed state)
- Datasets: testset.gl and testset.gs (dartR.data) with synthetic
  chromosome/position slots and a synthetic GFF3 written by the test;
  platypus.gl (1,000 loci, 921 mapped) placed on the NCBI platypus assembly
  GCF_004115215.1 (mOrnAna1.p.v1) using `Chrom_Platypus_Chrom_NCBIv1` and
  `ChromPos_Platypus_Chrom_NCBIv1 + SnpPosition`, with the matching NCBI GFF
- Baseline: tests/testthat/test-gl.find.loci.in.genes.R (8 tests,
  snapshot captured pre-review; defects marked `BASELINE (F<n>)`)

**Standards: Needs work** — the function has none of the house entry and exit
steps: `verbose` is never read, output is uncoloured `cat()`, and it
validates no input.
**Spec: Needs work** — the overlap step is correct, but gene detection
misses genes whose matching text sits on transcript rows, and a mistyped
file name can silently use another annotation.

What works well: once the matching genes are known, the `foverlaps` join is
correct (closed intervals, gzip handled both ways) and on the platypus data
it agrees with an independent overlap for every pattern where the gene
detection agrees.

## Findings

**F1 [HIGH, confidence: high] — a missing GFF file silently falls back to
any object called `gff` in the workspace (FS5)**
`R/gl.find.loci.in.genes.r:84-91` — `gff` is only assigned when
`gff.file` or `gff.file.gz` exists. Otherwise `as.data.table(gff)` looks the
name up outside the function (package namespace, then the global
environment). If the user has a `gff` object there, it is used as the
annotation. If not, the error is `object 'gff' not found`, or, when `gff`
is a file path string, "`string` must be a vector, not a primitive
function".
Failure scenario: a user reads the annotation with
`gff <- ape::read.gff(...)` to inspect it, then calls the function with a
typo in `gff.file`. The function returns loci from the workspace table
without any message (reproduced; the baseline test asserts it).
Proposed change: stop with `error()` naming both paths tried when neither
exists (change 1).

**F2 [HIGH, confidence: high] — genes are found only through their own row
or their CDS rows; text on mRNA, lnc_RNA and other transcript rows is
ignored (DOC5)**
`R/gl.find.loci.in.genes.r:106-129` — the pattern is tested against gene
rows and CDS rows only, and CDS matches are linked to genes by the `gene=`
key. NCBI GFF3 files put `product` on transcript rows (mRNA, lnc_RNA,
transcript) and on CDS rows. The CDS product often differs from the
transcript product, and non-coding genes have no CDS at all. `@description`
promises genes "whose attributes match". The synthetic LOC1 gene, with
"major histocompatibility" only on its mRNA, is missed.
Failure scenario: platypus.gl on the NCBI annotation, compared with an
independent overlap that walks every matching feature up its `Parent` chain
to its gene:

| Pattern | Loci (independent) | Loci (function) |
|---|---|---|
| `(?i)kinase` | 24 | 24 |
| `(?i)receptor` | 38 | 33 |
| `(?i)zinc finger` | 16 | 14 |
| `(?i)uncharacterized` | 50 | 6 |

RGMB is one missed gene: its mRNA product contains "receptor" but its CDS
product reads "LOW QUALITY PROTEIN: RGM domain family member B". The 44
missed "uncharacterized" loci sit in lncRNA genes.
Proposed change: test the pattern against every feature row and map each
match to its gene by walking `Parent` links, as `gl.find.genes.for.loci`
does for products (change 2).

**F3 [MEDIUM, confidence: high] — `verbose` is ignored and the house
structure is missing (FS2, FS3, FS9, VRB1, VRB2, VRB3)**
`R/gl.find.loci.in.genes.r:79, 133, 166-167` — no `gl.check.verbosity()`,
no `utils.flag.start()`, no "Completed" line. Three uncoloured `cat()`
lines print at every verbosity. The no-match warning uses base `warning()`
and says "No MHC genes detected … Consider widening 'mhc_pat'". That text
is wrong for any other pattern, and `mhc_pat` is not an argument.
Failure scenario: `gl.find.loci.in.genes(x, gff, gene = "TAP", verbose = 0)`
prints "LOADING GFF FILE...", "MHC genes detected: 1" and "Loci overlapping
MHC intervals: 1".
Proposed change: add the standard verbosity, flag-start and flag-end
blocks, use `report()`/`warn()`, gate progress at 2 and the summary at 3,
and use generic wording ("genes matching the pattern") (change 3).

**F4 [MEDIUM, confidence: high] — no datatype or argument checks (FS4, FS5)**
`R/gl.find.loci.in.genes.r:57-61` — nothing checks `x`, `gene`, or the
presence of chromosome and position. `gene = c("MHC", "TAP")` fails
inside data.table with "Recycling of logical i is no longer allowed";
`gene = NA` fails inside stringr; an object without chromosome or
position returns `character(0)` silently.
Failure scenario: `testset.gl` as shipped (no mapping) returns an empty
vector with no message after loading the GFF.
Proposed change: `utils.check.datatype(x, accept = c("SNP",
"SilicoDArT"))`, then check that `gene` is a single non-NA string and that
`x$chromosome` and `x$position` exist with length `nLoc(x)` (change 1).

**F5 [MEDIUM, confidence: high] — sequence-name mismatches and unmapped loci
are dropped without a message (VRB4, proposed rule)**
`R/gl.find.loci.in.genes.r:147-158` — loci with NA chromosome or position
are dropped without a count. A chromosome name that never occurs in the
GFF produces no overlap and no message.
Failure scenario: platypus.gl stores chromosomes as
`NC_041728.1_chromosome_1`; the NCBI GFF uses `NC_041728.1`. The function
finds 7 MHC genes, returns 0 loci, and gives no sign that no locus
sequence name matched. `gl.find.genes.for.loci` already reports both cases
at `verbose >= 1`.
Proposed change: report the number of loci without a position and the
sequence names absent from the GFF at `verbose >= 1`, with an explicit
warning when none match (change 4).

**F6 [LOW, confidence: medium] — pseudogene features are never matched
(principle: parity with `gl.find.genes.for.loci`, which includes
`pseudogene` by default)**
`R/gl.find.loci.in.genes.r:120` — only `type == "gene"` rows are
intervals. NCBI annotates pseudogenes as type `pseudogene`, and MHC
regions carry many of them.
Failure scenario: the synthetic "MHC class I pseudogene" P holds locus 4,
which is not returned. In the platypus data, no locus falls in a matching
pseudogene, so there is no real-data effect there.
Proposed change: treat `gene` and `pseudogene` rows as gene intervals
(change 5).

**F7 [LOW, confidence: high] — `save2tmp` is documented but does nothing
(DOC5)**
`R/gl.find.loci.in.genes.r:18-19, 60` — the argument is never read.
Failure scenario: `save2tmp = TRUE`; `gl.list.reports()` lists nothing.
Proposed change: save the table of matching genes and hits to `tempdir()`
as the sibling does (change 6).

**F8 [LOW, confidence: high] — roxygen gaps (DOC1, DOC2, DOC3, DOC5, DOC7
(proposed rule))**
`R/gl.find.loci.in.genes.r:1-56` — `@author` has no Author(s)/Custodian
parts; the `verbose` text is non-standard; `@family` sits after
`@examples`; the only example is `\dontrun{}` on a file that does not
exist. `@details` does not say that:
- `x$position` must be a genome coordinate. dartR sets it to the SNP
  position within the tag (platypus.gl: 36, 13, 23, …), so a user who
  sets only the chromosome gets meaningless overlaps.
- the pattern is matched against the whole attribute string. `"TAP"`
  matches 13 platypus genes (WTAP, METAP1, STAP1, TAPBP …), and `"ID"`
  matches every gene.
- `gl.find.genes.for.loci()` on the result says which gene each locus
  falls in.
Failure scenario: a reader following the help page gets zero or wrong
loci and cannot see why.
Proposed change: rewrite the header in house order with a runnable
example on a small temporary GFF (change 7).

**F9 [INFO] — output order**
Loci come back in genome order (chromosome, position), not in
`locNames(x)` order. Not documented either way; no change proposed.

## Proposed changes

1. Input checks: stop with `error()` when neither `gff.file` nor its `.gz`
   companion exists; datatype check accepting SNP and SilicoDArT; `gene`
   must be one non-NA string; `x$chromosome`/`x$position` must exist with
   length `nLoc(x)` (F1, F4). Valid calls return the same result; calls
   that today use a workspace `gff` object now error.
2. Match the pattern on every feature row and map each hit to its gene
   through the `Parent` chain (F2).
   **Consequence: the returned loci change: loci in genes matched only
   through transcript rows are added (platypus: +5 for "receptor", +44 for
   "uncharacterized").**
3. Standard verbosity and structure: `gl.check.verbosity`, flag start/end,
   `report()`/`warn()`, progress at 2, summary at 3, generic warning text
   (F3).
4. Report loci without position and unmatched sequence names at
   `verbose >= 1`; warn when no sequence name matches the GFF (F5).
5. Treat `pseudogene` rows as gene intervals (F6).
   **Consequence: the returned loci change: loci inside matching
   pseudogenes are added.**
6. Implement `save2tmp` (F7).
7. Roxygen rewrite with a runnable example and the three `@details` points
   (F8).

No change alters the signature. The only external caller found is the
dartr2shiny generator copy (`input_generator/dartR.popgen/`), which takes
the roxygen header; it needs regenerating after change 7.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run. DAT1–DAT4 and
  FS8 not applicable (returns a character vector, does not modify `x`).
  PLT not applicable (no plot).
- Spec: behaviour vs roxygen on the synthetic GFF — run.
- Real-data check: platypus.gl on the NCBI mOrnAna1.p.v1 GFF, compared
  with an independent Parent-chain overlap for four patterns — run.
  MHC pattern: 7 genes, 0 loci in both computations (1,000 loci is too
  sparse to hit them).
- SilicoDArT: runs on testset.gs with synthetic positions — run.
- DEP1: `ape`, `data.table` and `stringr` are in Imports — no guard
  needed; checked.
- dartR Google Group / GitHub issues: SKIPPED — no search run in this
  session.
- FBM path (DAT6): not applicable — genotypes are never read.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis |  |
| 2 | approved | Luis | consequence approved: returned loci change |
| 3 | approved | Luis |  |
| 4 | approved | Luis |  |
| 5 | approved | Luis | consequence approved: returned loci change |
| 6 | approved | Luis |  |
| 7 | approved | Luis |  |

## Outcome

- Changes 1-7 applied in `R/gl.find.loci.in.genes.r` (branch
  `review-find-loci-in-genes`); `man/gl.find.loci.in.genes.Rd` regenerated;
  NEWS entry added.
- Snapshot diffs against the pre-review baseline: 6, all mapped. Locus 2
  (mRNA-only match) added by change 2; locus 4 (pseudogene) added by change
  5, in two tests; `verbose = 0` output and the base `warning()` gone by
  change 3; missing file now errors by change 1, including with a `gff`
  object in the global environment. Tests rewritten to assert the approved
  behaviour: 12 tests, 31 expectations, all pass.
- Real-data check (platypus.gl, NCBI mOrnAna1.p.v1 GFF): the function now
  equals the independent Parent-chain overlap for all four patterns
  (kinase 24/24, receptor 38/38, zinc finger 16/16, uncharacterized 50/50;
  before: 24, 33, 14, 6). With the DArT sequence names
  (`NC_041728.1_chromosome_1`) the run at `verbose = 1` now warns that no
  locus sequence name matches the GFF.
- `devtools::document()` also rewrote 7 unrelated Rd files (duplicated
  family link lists left by an earlier merge); reverted to keep this PR to
  one function. Worth a separate `document()` commit on `dev`.
- PR: dartR.popgen#109 (commit 9db88dd, branch `review-find-loci-in-genes`).

## Machine block

```json
{
  "function": "gl.find.loci.in.genes",
  "package": "dartR.popgen",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "f2c8b8a",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "FS5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "FS2", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS4", "status": "approved", "change": 1},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "VRB4", "status": "approved", "change": 4},
    {"id": "F6", "severity": "LOW", "confidence": "medium", "rule": "principle: sibling parity", "status": "approved", "change": 5},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 6},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 7},
    {"id": "F9", "severity": "INFO", "confidence": "high", "rule": "DOC5", "status": "no_change", "change": null}
  ],
  "coverage_skipped": ["forum/issue search: not run in this session", "DAT6: genotypes not read"],
  "status": "done",
  "pr": 109
}
```
