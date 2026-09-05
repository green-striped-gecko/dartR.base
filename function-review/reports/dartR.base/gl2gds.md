# Review: gl2gds (dartR.base)

- Family mode: io
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203
- Datasets: testset.gl, testset.gs, platypus.gl
- Baseline: tests/testthat/test-gl2gds.R (snapshot captured pre-review)

## Verdict

**Standards: Needs work** — structure and docs are largely in shape, but the
function is loudly non-silent at `verbose = 0` and admits SilicoDArT.
**Spec: Rework** — with `snp.pos`/`snp.chr` supplied, the written gds
assigns each SNP the position and chromosome of a *different* locus (7 of 8
records wrong in the platypus demonstration), and the stored dosage is
inverted relative to SNPRelate's allele convention. Both corruptions are
silent; every positional or frequency-aware downstream analysis is affected.

## Findings

**F1 [BLOCKER, confidence: high] — sorted output mixes fields from different loci (DAT2 principle, spec)**
`R/gl2gds.r:123-136` — after
`snp_order_temp <- snp_order_temp[order(chrom, snp.pos), ]`, the column
`snp_order_temp$order` holds the sort permutation `p` (original index of
each sorted row). The genotype matrix, `snp.id`, and `snp.allele` are then
subset with `order(snp_order_temp$order)` — the *inverse* permutation —
while `snp.position`/`snp.chromosome` take the forward-sorted values. The
id/allele/genotype triple stays internally consistent, but its record order
differs from the coordinate columns, so record k pairs locus
`inverse_p[k]`'s identity and genotypes with locus `p[k]`'s coordinates.
Failure scenario: `gl2gds(platypus.gl, snp.pos =
'ChromPos_Platypus_Chrom_NCBIv1', snp.chr = 'Chrom_Platypus_Chrom_NCBIv1')`
— on an 8-locus subset, 7 of 8 records carry another locus's position
(empirically confirmed; baseline test freezes this). Any SNPRelate use of
coordinates — LD pruning by distance, per-chromosome analysis,
`autosome.only` filtering, export to other tools — silently operates on
wrong coordinates. Defaults (`snp.chr = "0"`, `snp.pos = "0"`) escape only
because the sort is a no-op.
Proposed change: compute `ord <- order(chrom, snp.pos)` once on the
un-sorted table and subset *everything* (genotype rows, snp.id, snp.allele,
positions, chromosomes) by `ord`; delete the `order` bookkeeping column.

**F2 [HIGH, confidence: high] — dosage written in the opposite orientation to snp.allele (spec)**
`R/gl2gds.r:127-129,148-149` — the genlight dosage (count of the *second*
allele in `@loc.all`, "ref/alt") is copied into the gds unchanged, but
`snpgdsCreateGeno` defines genotype as the count of the *first* (A) allele
of `snp.allele`. Empirically: a locus with first-allele frequency 0 in the
genlight reports `AlleleFreq = 1.0` from `snpgdsSNPRateFreq` (baseline test
comment records this).
Failure scenario: allele frequencies, reference-allele reports, and any
export from the gds swap ref and alt for every locus; PCA/IBD distances are
unaffected, which makes the flip easy to miss.
Proposed change: write `2 - dosage` (keeping NA -> 3), or swap each
`snp.allele` string to "alt/ref"; the former preserves the documented
allele order.

**F3 [MEDIUM, confidence: high] — 16+ lines of ungated console output at verbose = 0 (VRB5)**
`R/gl2gds.r:177-179` — `cat(important("Structure of gds file\n\n"))`,
`SNPRelate::snpgdsSummary(genofile)` and `print(genofile)` run
unconditionally (confirmed at `verbose = 0`); SNPRelate's own
invalid-position warnings add stderr noise whenever the defaults leave
pos/chrom at 0.
Failure scenario: silent pipelines get a full gds structure dump per call.
Proposed change: gate the summary block at `verbose >= 3` (it is a results
summary) and the banner via `report()` at the same level.

**F4 [MEDIUM, confidence: high] — SilicoDArT admitted (DAT7)**
`R/gl2gds.r:85` — permissive datatype default; a presence/absence object
(ploidy 1, `@loc.all` NULL) converts without error into a "SNP" gds with
absence coded as homozygote and `snp.allele` omitted (confirmed).
Failure scenario: `gl2gds(testset.gs)` yields a structurally valid gds that
SNPRelate analyses as diploid SNP data.
Proposed change: `accept = "SNP"` in `utils.check.datatype()`.

**F5 [LOW, confidence: high] — verbose text and author block off-canon (DOC2, DOC7) (proposed rule)**
`R/gl2gds.r:21-23,49-50` — `verbose` default clause reads "[default 2 or as
specified using gl.set.verbosity]" (DOC2 canon differs); `@author` names a
custodian only, with no `Author(s):` part (DOC7).
Failure scenario: doc drift; authorship unstated.
Proposed change: canonical DOC2 clause; `Author(s): Luis Mijangos.
Custodian: Luis Mijangos -- Post to ...`.

Genome-position note (no finding): coordinates come from named
`loc.metrics` fields via `snp.pos`/`snp.chr`, not from
`@position`/`@chromosome` — consistent with the convention that the slots
hold genome coordinates only and may be NULL. With the defaults the gds
carries all-zero coordinates ("unknown" per SNPRelate), which the `@details`
text already explains.

## Proposed changes

1. Fix the sort so genotype/id/allele and position/chromosome share one
   permutation (F1). **Consequence: numerical output changes for every call
   that supplies snp.pos/snp.chr — the current misaligned files were wrong.**
2. Align dosage orientation with snp.allele (`2 - dosage`) (F2).
   **Consequence: stored genotype values change for all loci; allele
   frequencies from the gds become correct.**
3. Gate the summary/print block at `verbose >= 3` (F3).
4. `accept = "SNP"` (F4).
5. Canonical DOC2 verbose text and DOC7 author labels (F5).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP (SNPRelate is in Imports — no DEP1
  guard required), PLT (n/a), STY — run
- Spec: behaviour vs roxygen and the SNPRelate gds contract — run
  empirically (SNPRelate 1.x installed): default path, sorted path
  (platypus), dosage orientation via snpgdsSNPRateFreq
- SilicoDArT probe — run (F4)
- verbose = 0 silence probe — run (F3)
- chr.format = "numeric" branch: read statically only — no numeric-chromosome
  fixture was run
- FBM path (DAT6): SKIPPED — no FBM fixture; note `t(as.matrix(x))`
  densifies the whole object
- Genome-position lens: run — see note above

## Approval (Phase B)

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges (2026-09-05) | Blanket approval, incl. output-change consequence |
| 2 | Approved | Arthur Georges (2026-09-05) | Blanket approval, incl. output-change consequence |
| 3 | Approved | Arthur Georges (2026-09-05) | Blanket approval, all severities |
| 4 | Approved | Arthur Georges (2026-09-05) | Blanket approval, all severities |
| 5 | Approved | Arthur Georges (2026-09-05) | Blanket approval, all severities |

All findings approved 2026-09-05 via the formal approval boxes: blanket
class approval at every severity, explicitly including the consequence
that converted outputs change where they were wrong (F1 coordinates, F2
stored dosage), and the DAT7 fatal `accept = "SNP"` gate wherever
SilicoDArT was silently admitted.

## Outcome (Phase C)

All five changes applied on branch `review-gl2gds` (base `upstream/dev`,
ddaed27). One permutation (`ord <- order(chrom, snp.pos)`) now reorders
genotype rows, snp.id, snp.allele, positions and chromosomes together
(F1); the stored genotype is `2 - dosage` with NA -> 3 (F2); the
banner/summary/print block is gated at `verbose >= 3` (F3);
`accept = "SNP"` (F4); DOC2/DOC7 doc canon (F5). Baseline
characterization test updated for four approved diffs (F1 coordinate
mapping, F2 dosage orientation, F3 silence, F4 rejection); 16/16
assertions pass. Empirical SNPRelate verification on the full
`platypus.gl` (1000 loci) with `snp.pos`/`snp.chr` supplied: every
record's position and chromosome match the supplied `loc.metrics`
fields 1:1 (with the function's documented NA-position -> 0 recode),
and `snpgdsSNPRateFreq` allele frequencies equal the genlight
first-allele frequencies (`all.equal` TRUE; NaN only at all-missing
loci, on both sides). End-to-end run on `testset.gl` at `verbose = 3`
wrote a 274 x 755 gds with the structure summary printed at level 3.
Sibling-caller grep across the 8 dartRverse clones: no callers of
`gl2gds` — all clear. NEWS entry added.

```json
{
  "function": "gl2gds",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "BLOCKER", "confidence": "high", "rule": "DAT2", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DAT7", "status": "applied", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "DOC2", "status": "applied", "change": 5}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "chr.format='numeric' branch not executed"],
  "status": "phase-c",
  "pr": null
}
```
