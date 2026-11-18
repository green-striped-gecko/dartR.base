# Review: gl2gi (dartR.base)

- Family mode: io (reviewed as a pair with `gi2gl`; round-trip evidence in
  the `gi2gl` report)
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203
- Datasets: testset.gl, testset.gs, adegenet::glSim fixture
- Baseline: tests/testthat/test-gl2gi.R (snapshot captured pre-review)

## Verdict

**Standards: Needs work** — the preamble conforms, but locus metadata desyncs
from the genotypes on the packaged reference data, and SilicoDArT input is
admitted where only SNP data makes sense.
**Spec: Needs work** — genotype coding is verified correct for biallelic SNP
data (0 -> hom first allele, 1 -> het, 2 -> hom second allele, NA survives),
individual names, populations, and ploidy 2 all carry over; but the promise
of "a genind object, with all slots filled" fails on data containing all-NA
loci, and a single-locus object crashes.

## Findings

**F1 [HIGH, confidence: high] — loc.metrics desync when all-NA loci are dropped (DAT2)**
`R/gl2gi.r:88-99` — `df2genind` silently removes loci with no scored alleles
("Markers with no scored alleles have been removed"), then `gen@other <-
x@other` copies the full locus metadata wholesale.
Failure scenario: `gl2gi(testset.gl)` returns a 752-locus genind carrying a
755-row `loc.metrics` table. Every row past the first dropped locus describes
the wrong locus, and `gi2gl(gl2gi(testset.gl))` crashes inside
`gl.recalc.metrics` ("replacement has 752 rows, data has 755") — the round
trip fails outright on the packaged reference dataset.
Proposed change: remove all-NA loci explicitly before conversion (warning at
`verbose >= 1`, per the `gl2geno` precedent), then subset `loc.metrics` to
the loci actually present in the genind (DAT3 idiom) before assigning
`gen@other`.

**F2 [HIGH, confidence: high] — SilicoDArT admitted (DAT7)**
`R/gl2gi.r:37` — `utils.check.datatype(x, verbose = verbose)` uses the
default `accept`, so presence/absence data passes.
Failure scenario: `gl2gi(testset.gs)` runs to completion and returns a
ploidy-2 genind in which presence/absence scores are re-labelled as
homozygote/heterozygote calls at fabricated A/T loci — structurally valid,
biologically meaningless, and produced without any warning. Verified
empirically.
Proposed change: `accept = "SNP"`.

**F3 [MEDIUM, confidence: high] — single-locus genlight crashes (FS5)**
`R/gl2gi.r:40,90` — `as.matrix(x[,])` and `xx[,]` use default `drop = TRUE`,
so a one-locus object collapses to a vector before reaching `df2genind`.
Failure scenario: `gl2gi(x[, 1])` fails with the opaque adegenet error
"X is not a matrix". Verified empirically.
Proposed change: `xx[, , drop = FALSE]` (and the same for the initial
`as.matrix` subset).

**F4 [MEDIUM, confidence: high] — silent fabrication of allele labels (DAT5, DOC5 proposed rule)**
`R/gl2gi.r:51-54` — when `x@loc.all` is NULL the function invents `"A/T"`
for every locus (`"C/G"` for locus 1) with no message at any verbosity.
Failure scenario: a genlight not built by dartR (e.g. `glSim`, or the output
of `gi2gl`) converts cleanly, and the resulting genind reports concrete
nucleotide alleles that were never in the data; downstream allele-frequency
tables and exports present the fabrication as fact. Verified empirically.
Proposed change: warn at `verbose >= 1` that allele labels are placeholders,
or use neutral labels (`"1/2"`) that cannot be mistaken for sequence data.

**F5 [LOW, confidence: high] — progress message at the begin/end level (VRB1)**
`R/gl2gi.r:84-86` — "Matrix converted.. Prepare genind object..." prints at
`verbose >= 1`; level 1 is reserved for begin/end.
Failure scenario: a user at verbosity 1 receives mid-run progress chatter.
Proposed change: gate at `verbose >= 2`.

**F6 [LOW, confidence: high] — roxygen gaps and a wrong documented default (DOC1, DOC5 proposed rule)**
`R/gl2gi.r:1-21` — no `@description` tag; `probar` is documented
"[default TRUE]" but the signature default is `FALSE`; `@return` sits after
`@export`.
Failure scenario: a user reading the manual expects a progress bar by default
and does not get one.
Proposed change: add `@description`, correct the `probar` default text,
reorder tags to the ratified order, and run `devtools::document()` (DOC4).

**F7 [LOW, confidence: high] — author block lacks Custodian line (DOC7) (proposed rule)**
`R/gl2gi.r:12` — `@author` names Bernd Gruber with no `Custodian:` label.
Proposed change: add `Custodian: Bernd Gruber` per the DOC7 default.

**F8 [LOW, confidence: medium] — full densification plus elementwise double loop (DAT6, STY2) (proposed rule)**
`R/gl2gi.r:40,62-79` — the whole genlight is densified and then recoded one
cell at a time in nested loops.
Failure scenario: on a large or FBM-backed object the function materialises
the full matrix and spends O(nInd x nLoc) scalar iterations; a 100k-locus
object is slow and memory-hungry where a vectorised `matrix(hets[col],...)`
lookup would not be.
Proposed change: vectorise the recode (index `homs1/hets/homs2` by column in
one pass); defer FBM support to the team-level DAT6 decision.

## Proposed changes

1. Remove all-NA loci up front with a gated warning and subset `loc.metrics`
   to the surviving loci before assigning `gen@other` (F1).
2. Restrict datatype with `accept = "SNP"` (F2).
   **Consequence: SilicoDArT callers that previously received a fabricated
   genind now receive an error.**
3. Add `drop = FALSE` to the two matrix subsets (F3).
4. Warn (or use neutral labels) when allele labels are fabricated from a
   NULL `loc.all` (F4).
5. Re-gate the progress message to `verbose >= 2` (F5).
6. Roxygen: add `@description`, fix the `probar` default text, reorder tags,
   add the Custodian line, re-document (F6, F7).
7. Vectorise the genotype recode loop (F8).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run (no plot bundle; no
  file output so FS7 not applicable; adegenet is a Depends package so DEP1
  not applicable; FS8 not applicable — the return is a genind, not a
  modified genlight).
- Spec: genotype coding verified allele-by-allele against `tab(gi)` on a
  testset.gl subset; indNames/pop/ploidy/locNames survival verified on full
  testset.gl; NULL-`loc.all` path run on `glSim`; single-locus edge run;
  SilicoDArT probe run — run.
- `verbose = 0`: zero console lines verified. Note: the un-gated `df2genind`
  R-level warning still fires at verbose 0; treated as part of F1 rather
  than a separate VRB5 finding, since the fix removes the condition.
- Round trip gl -> gi -> gl: run; results reported in the `gi2gl` report.
- FBM path (DAT6): SKIPPED — no FBM fixture; hazard noted statically in F8.

## Approval (Phase B)

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges | 2026-09-05, via approval boxes |
| 2 | Approved | Arthur Georges | 2026-09-05, via approval boxes |
| 3 | Approved | Arthur Georges | 2026-09-05, via approval boxes |
| 4 | Approved | Arthur Georges | 2026-09-05, via approval boxes |
| 5 | Approved | Arthur Georges | 2026-09-05, via approval boxes |
| 6 | Approved | Arthur Georges | 2026-09-05, via approval boxes |
| 7 | Approved | Arthur Georges | 2026-09-05, via approval boxes |

All findings at every severity approved 2026-09-05 (blanket class
approval, explicitly acknowledging that converted outputs change where
they were wrong, and including the DAT7 fatal `accept = "SNP"` gate —
applied here as change 2).

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-gl2gi` (base `upstream/dev`,
ddaed27). PR: green-striped-gecko/dartR.base#338.

- F1: all-NA loci removed up front (warning at `verbose >= 1`);
  `loc.metrics` re-subset from the original object (DAT3), and a
  belt-and-braces re-register against `locNames(gen)` after `df2genind`.
  Implemented inline rather than via `gl.filter.allna`, which would also
  have removed all-NA individuals — beyond the approved scope.
- F2: `accept = "SNP"` — SilicoDArT input now errors.
- F3: `drop = FALSE` on the `df2genind` input plus a matrix guard on the
  `as.matrix` result; single-locus objects convert.
- F4: warning at `verbose >= 1` when placeholder alleles are fabricated
  (labels kept as A/T, C/G-first for back-compatibility).
- F5: progress message re-gated to `verbose >= 2`.
- F6/F7: `@description` added, `probar` default corrected to FALSE, tags
  reordered, DOC2 verbose clause, `Author(s)`/`Custodian` block;
  re-documented.
- F8: recode vectorised (three indexed assignments replace the
  elementwise double loop); `probar` retained for back-compatibility and
  documented as such.

Verification: 27 baseline assertions pass. Updated expectations, each
mapped to an approved finding: loc.metrics rows 755 -> 752 (F1);
SilicoDArT now `expect_error` (F2); single-locus now converts (F3).
End-to-end `gl2gi(testset.gl, verbose = 3)`: 752 loci / 752 loc.metrics
rows / 274 individuals / 274 ind.metrics rows, gated messages correct.
Sibling caller grep across the 8 clones: internal callers
`gl2genalex`, `gl2genepop`, `gl.report.ld` (dartR.base) and
`gl.outflank`, `gl.run.structure` (dartR.popgen),
`gl.diagnostics.sim` (dartR.sim), `gl.genleastcost` (dartR.spatial) —
all SNP contexts, none breaks; a SilicoDArT object routed through them
now errors in `gl2gi`, which is the approved DAT7 consequence. NEWS
entry added.

Discovered during verification (not fixed — outside approved scope,
recorded for a follow-up finding): `df2genind` also removes individuals
with no scored genotypes across the converted loci while `ind.metrics`
is copied wholesale — the individual-side analogue of F1 (DAT2).
Pre-existing; unreachable on the packaged multi-locus datasets, visible
on a single-locus conversion (274 -> 273 individuals with 274-row
`ind.metrics`). Characterized in the test file.

### Round-trip re-verification with both fixes (this branch + gi2gl PR #344)

Run 2026-09-05 with this branch's `gl2gi.r` and `review-gi2gl`'s
`gi2gl.R` loaded together:

- `gl.filter.allna(testset.gl)` (752 loci): dimensions, indNames,
  locNames, pop and NA pattern identical; `loc.metrics` 752 rows,
  `ind.metrics` 274 rows (in register); 596 of 752 loci come back as
  exact 0 <-> 2 complements and every locus is identical-or-complement
  (the documented reference-allele ambiguity, unchanged); `loc.all` is
  set-equal to the original pair on all 611 loci biallelic in the data.
- Raw `testset.gl` (755 loci incl. 3 all-NA): the round trip now
  SUCCEEDS — 752 loci with 752-row `loc.metrics` (previously crashed in
  `gl.recalc.metrics`).

```json
{
  "function": "gl2gi",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT2", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DAT7", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DAT5", "status": "applied", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "VRB1", "status": "applied", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "applied", "change": 6},
    {"id": "F8", "severity": "LOW", "confidence": "medium", "rule": "DAT6", "status": "applied", "change": 7}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture"],
  "status": "phase-c-complete",
  "pr": 338
}
```
