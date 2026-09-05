# Review: gl2eigenstrat (dartR.base)

- Family mode: io
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203
- Datasets: testset.gl, testset.gs, platypus.gl
- Baseline: tests/testthat/test-gl2eigenstrat.R (snapshot captured pre-review; 15 assertions passing)

## Verdicts

**Standards: Needs work** — preamble, `gl.check.wd`, `file.path`-based output and `invisible(NULL)` return all conform, but the datatype check admits SilicoDArT and metadata fields pass through an unsafe factor-to-numeric coercion.

**Spec: Rework** — the genotype file encodes the opposite of what EIGENSTRAT defines (variant-allele counts written where reference-allele counts belong, relative to the ref/var columns the function itself writes), the `@details` promise to remove illegal-chromosome SNPs is not implemented, and factor chromosome fields are exported as factor level codes. The three files are well-formed; their content misrepresents the data.

## Findings

**F1 [HIGH, confidence: high] — geno file counts the wrong allele (DOC5, proposed rule; target-format spec)**
`R/gl2eigenstrat.r:91-104,122-126` — the dartR score (0 = homozygous for the first `loc.all` allele; the score counts the second, variant allele) is written unchanged, while the EIGENSTRAT geno format defines each value as the number of copies of the reference allele — the allele this function writes in `.snp` column 5 (`ref_allele <- substring(x@loc.all, 1, 1)`).
Failure scenario: verified 2026-09-05 — testset.gl individuals 1-5 at locus 1 (`C/T`) score `2 2 NA 2 2` (homozygous T) and the file line starts `22922`, which EIGENSTRAT reads as homozygous C. Every genotype is allele-flipped relative to the declared ref/var; allele frequencies, f-statistics and any downstream analysis keyed to allele identity are inverted.
Proposed change: write `2 - score` (keeping 9 for missing), or swap the `.snp` ref/var columns — one of the two, consistently. **Consequence: numerical output changes for every exported dataset.**

**F2 [HIGH, confidence: high] — factor metadata fields exported as factor level codes (DAT5; no exact catalogue rule — factor-coercion principle)**
`R/gl2eigenstrat.r:110-120` — `as.numeric(unname(unlist(snp_temp[snp.chr])))` on a factor column returns the level codes, not the values. dartR loc.metrics fields imported from CSV are routinely factors.
Failure scenario: verified 2026-09-05 — the roxygen example itself (`platypus.gl`, `snp.chr = 'Chrom_Platypus_Chrom_NCBIv1'`, a factor) writes chromosome "numbers" 1-103 that are alphabetical level indices; chromosome X1 gets an arbitrary code, not the 23 the `@details` prescribes. The same coercion path applies to `snp.pos`, where a factor field would silently corrupt every position.
Proposed change: coerce via `as.numeric(as.character(...))`, then validate: map X/Y/mtDNA per the documented encoding or stop on non-numeric chromosome labels.

**F3 [HIGH, confidence: high] — documented removal of illegal-chromosome SNPs is not implemented (DOC5, proposed rule)**
`R/gl2eigenstrat.r:33-37` vs the function body — `@details` states "SNPs with illegal chromosome values, such as 0, will be removed"; no code removes anything.
Failure scenario: verified 2026-09-05 — the platypus run writes all 1000 loci regardless of chromosome value. A user relying on the documented clean-up ships an invalid file to EIGENSOFT.
Proposed change: implement the filter (with a VRB4-level warning stating how many loci were dropped) or delete the claim from the docs.

**F4 [MEDIUM, confidence: high] — SilicoDArT admitted; malformed .snp output (DAT7)**
`R/gl2eigenstrat.r:86` — permissive `utils.check.datatype` default.
Failure scenario: verified 2026-09-05 — `gl2eigenstrat(testset.gs, ...)` runs; presence/absence 0/1 is written as diploid dosage, and because SilicoDArT has no `loc.all` the ref/var columns vanish, leaving a 4-column `.snp` file where the layout is positional.
Proposed change: `accept = "SNP"` in the datatype call.

**F5 [MEDIUM, confidence: high] — parameter docs disagree with the signature (DOC5, proposed rule)**
`R/gl2eigenstrat.r:23-24,68` — `pos.cM` is documented "[default 1]" but the signature default is `0` (the file's genetic-position column is all `0`, verified); `snp.pos`/`snp.char` docs call the default `1` a "field name", when it is a numeric sentinel meaning "write the constant 1".
Failure scenario: a user reading the docs expects genetic position 1 and gets 0; the sentinel semantics are discoverable only from source.
Proposed change: align docs with the signature and describe the sentinel behaviour (constant column when left at default).

**F6 [LOW, confidence: medium] — sex.code/phen.value lengths never validated (FS5)**
`R/gl2eigenstrat.r:140-144` — `cbind` recycles short vectors; a sex vector whose length divides nInd recycles without any warning.
Failure scenario: a 2-value `sex.code` on an even-sized dataset alternates M/F down the `.ind` file silently.
Proposed change: check `length(sex.code) %in% c(1, nInd(x))` (same for `phen.value`), else stop.

**F7 [LOW, confidence: high] — author block names a custodian only (DOC7) (proposed rule)**
`R/gl2eigenstrat.r:39-40` — "Custodian: Luis Mijangos" with no `Author(s):` line. `@return` also precedes `@export` at the tail of the block rather than sitting in house position.
Proposed change: DOC7 author structure; house tag order.

## Proposed changes

1. Fix the geno-file allele orientation: write `2 - score` (NA stays 9) so geno counts the `.snp` reference allele (F1). **Consequence: numerical output changes for every exported dataset.**
2. Coerce metadata fields via `as.numeric(as.character(...))` and validate chromosome values against the documented numeric encoding (F2).
3. Implement the illegal-chromosome removal (warning at `verbose >= 1`) or delete the claim from `@details` (F3).
4. Add `accept = "SNP"` to the datatype check (F4).
5. Align `pos.cM`/`snp.pos`/`snp.chr` docs with actual defaults and sentinel semantics (F5).
6. Validate `sex.code`/`phen.value` lengths (F6).
7. Roxygen pass: DOC7 author structure, tag order (F7).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP (no Suggests packages used), PLT (n/a), STY — run
- Spec: all three files generated on testset.gl and inspected against the EIGENSTRAT
  format description (one geno line per SNP, one char per individual, 9 = missing,
  .snp and .ind column layouts) — run; genotype orientation cross-checked against
  `as.matrix(x)` and `loc.all` (F1)
- Roxygen example path (platypus.gl with named chrom/pos fields) — run (F2, F3)
- SilicoDArT admission probe on testset.gs — run (admitted; F4)
- verbose = 0 silence via `capture.output` — run (0 lines)
- Genome-position lens: the function never reads `@position`/`@chromosome` (conforms —
  coordinates come only from nominated loc.metrics fields); defaults emit constant
  chrom = 1 / physical position = 1, which EIGENSOFT tools accept for unmapped data
- FBM path (DAT6): SKIPPED — no FBM fixture
- EIGENSOFT round-trip (convertf/smartpca): SKIPPED — binaries unavailable; format
  checked against the EIGENSTRAT/convertf documentation

## Approval (Phase B)

All findings at every severity approved 2026-09-05 by Arthur Georges via the
formal approval boxes, as part of the blanket class approval for the io
converter batch: converted outputs change where they were wrong (explicitly
including change 1's numerical inversion for every exported dataset), and the
DAT7 class fix (fatal `accept = "SNP"` gate wherever SilicoDArT is silently
admitted) applies.

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges | geno = 2 - score; output changes for every dataset |
| 2 | Approved | Arthur Georges | as.character coercion + X/Y/MT/XY mapping; stop when nothing is encodable |
| 3 | Approved | Arthur Georges | Filter implemented (warning at verbose >= 1) |
| 4 | Approved | Arthur Georges | DAT7 fatal gate (class approval) |
| 5 | Approved | Arthur Georges | |
| 6 | Approved | Arthur Georges | |
| 7 | Approved | Arthur Georges | |

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-gl2eigenstrat` (base `upstream/dev` =
ddaed27). All seven changes applied. Composition note for changes 2/3: labels
'X', 'Y', 'MT'/'mtDNA', 'XY' map to 23/24/90/91; remaining illegal values
(non-numeric label, or < 1) are removed with a warning at `verbose >= 1`
(F3's filter option); when NO locus is encodable the function stops with a
fatal error (F2's stop option) rather than writing an empty file set -- the
platypus labels ('NC_041731.1_chromosome_4', ...) are all un-encodable, so
the roxygen example now nominates only the position field. Baseline
characterization test updated for the approved behaviour changes (F1
inversion 22922 -> 00900 on ind 1-5 x locus 1, F2/F3 fatal on platypus
chrom field, F4 SilicoDArT rejection, F6 length validation) and extended
with a synthetic mapping/removal test (8 loci with chrom labels
1,X,Y,0,MT,scaffold_7,XY,2 -> 6 rows written as 1,23,24,90,91,2, warning
names the 2 removed loci, geno file shrinks in step). 21 assertions pass.
End-to-end runs at `verbose = 3` on testset.gl (755 rows, `.snp` layout
unchanged) and on platypus.gl with the position field verified.
Sibling-clone grep across the 8 dartR-verse packages: no callers of
`gl2eigenstrat` -- all clear. NEWS.md entry added.
PR: green-striped-gecko/dartR.base#355.

```json
{
  "function": "gl2eigenstrat",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203482f4cc51dd2cef70f07493200e534ee7",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "principle:factor-coercion", "status": "applied", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DAT7", "status": "applied", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "medium", "rule": "FS5", "status": "applied", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "applied", "change": 7}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "EIGENSOFT round-trip: binaries unavailable"],
  "status": "pr-open",
  "pr": 355
}
```
