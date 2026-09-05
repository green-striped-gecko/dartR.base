# Review: utils.read.fasta (dartR.base)

- Family mode: io (internal engine for gl.read.fasta; not exported)
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203 (R/utils.read.fasta.r byte-identical to
  upstream/dev ddaed27)
- Datasets: synthetic FASTA fixtures (inventory below); locus_1.fas
  (dartR.data extdata) exercised via gl.read.fasta
- Baseline: tests/testthat/test-utils.read.fasta.R (29 assertions, snapshot
  captured pre-review, defective behaviour pinned with [BUG] tags)

## Verdict

**Standards: Needs work** — the signature's self-referential defaults, the
ungated console output and the cat()-as-return idiom all need fixing;
`merge_gl_fasta` (same file) shares the report.
**Spec: Rework** — on a clean, uppercase, equal-length, two-line FASTA the
happy path is correct (hand-verified: ambiguity hets R/Y/K/M -> 1, ref/alt
and locus naming right), but every departure from that narrow case corrupts
genotypes silently: missing data and gaps become heterozygotes, third
alleles become heterozygotes, lowercase input fabricates alleles, unequal
lengths recycle, and the documented reference-allele rule inverts when the
modal genotype is a het. The genotype-classification core needs redesign,
not patching.

What works well: the two-pass pool/classify structure keeps memory modest,
IUPAC ambiguity codes at genuinely biallelic uppercase columns convert
correctly, and locus names encode source file plus alignment column
faithfully (verified even after multiallelic drops).

## Synthetic fixtures

All written by the baseline test file itself (tempdir) and mirrored in the
review scratchpad; two-line FASTA records throughout except `wrapped.fas`.

| Fixture | Purpose |
|---|---|
| `clean.fas` | 6 individuals x 10 columns: A/G SNP with R het + N (col 3), C/T SNP with Y het (col 5), triallelic homozygote-classes A/C/G (col 6), triallelic with het-minority S (col 7, the only one the skip catches), C/T SNP with a gap (col 8); rest monomorphic. Every expected genotype hand-computed |
| `clean2.fas` | 5 individuals (ind1-4, ind7) x 6 columns: A/C SNP with M het, G/T SNP with K het; merge partner for gl.read.fasta |
| `lowercase.fas` | softmasked bases: case-mixed A/G column and case-mixed monomorphic-A column |
| `unequal.fas` | 10/8/10 bp sequences; recycling probe |
| `wrapped.fas` | sequences wrapped over two lines; format assumption probe |
| `dupnames.fas` | two records named ind1; duplicate-name handling |
| `mono.fas` | no polymorphism; NULL-return path |
| `refinv.fas` | column Y,Y,Y,T,T,C — modal genotype is the het; ref-polarity probe |

## Findings

**F1 [HIGH, confidence: high] — missing data is coded as heterozygous, never NA (DAT1, DOC5 proposed rule)**
`R/utils.read.fasta.r:181-194` — genotype assignment is a three-way
fall-through: hom-ref 0, hom-alt 2, everything else 1. `N`, `-`, and
V/H/D/B — which `gl.read.fasta`'s documentation explicitly lists as
"taken as missing data" — are excluded from the allele pool (lines
132-146) but the individuals carrying them fall into the else-branch.
Verified: an `N` and a `-` cell both return 1 where NA is expected; no NA
is ever produced within a file.
Failure scenario: any alignment with masked or missing sites — the norm —
inflates heterozygosity silently; downstream He/Ho, HWE and relatedness
estimates shift with no warning.
Proposed change: classify hom-ref/het/hom-alt explicitly and return NA for
anything outside the recognized genotype set.

**F2 [HIGH, confidence: high] — third alleles survive as fake heterozygotes; the documented skip rarely fires (DAT1, DOC5 proposed rule)**
`R/utils.read.fasta.r:148-167,181-194` — ref/alt come from the most and
least frequent GENOTYPE classes (`which.max`/`which.min` on the class
table). A truly triallelic column of three homozygote classes (A/A x3,
C/C x2, G/G x1) yields ref A, alt G — two letters — so the >2-alleles
detection never triggers; the column is kept, labelled `A/G`, and every
C/C individual is coded 1 by F1's fall-through (verified). The skip fires
only when the least-frequent class is a heterozygote spanning a third
letter (also verified: a column with minority S=C/G is correctly dropped).
Failure scenario: "SNPs with more than two alleles are skipped" is false
for the commonest triallelic configuration; plausible wrong genotypes under
a wrong allele pair.
Proposed change: detect multiallelism from the allele pool (which the
function already computes at lines 80-88), not from the two extreme
genotype classes.

**F3 [HIGH, confidence: high] — lowercase bases fabricate alleles and heterozygotes (DAT1)**
`R/utils.read.fasta.r:82,115-124` — `letterOK` admits `a/c/t/g`, but the
hom-expansion (`"A" -> "AA"` etc.) and the ambiguity substitutions are
uppercase-only, so a lowercase base survives as a one-character genotype
class that can never match `hom_ref`/`hom_alt`. Verified on a softmasked
fixture: a case-mixed monomorphic column returns as a locus with alleles
`a/a` and every individual heterozygous; a genuine A/G column with one
lowercase `a` reports alleles `A/a` and codes the G/G individual as het.
Which variant becomes ref additionally depends on locale-sensitive table
name ordering.
Failure scenario: softmasked FASTA (repeat-masker output, many aligners)
is standard; whole objects are fabricated without any message.
Proposed change: `toupper()` the sequence data once on read (or reject
lowercase input loudly).

**F4 [HIGH, confidence: high] — unequal-length sequences recycle into fabricated loci (DAT2)**
`R/utils.read.fasta.r:70-73` — `matrix(unlist(txt), byrow = TRUE, nrow =
length(txt))` assumes equal lengths; a short sequence shifts every
subsequent base and base-R recycling pads the tail. Verified: a 10/8/10 bp
trio emits only the base-R "data length [28] is not a sub-multiple..."
warning and returns 8 loci where the alignment contains exactly one true
SNP — the killer silent-recycling class. The second pass (line 105-112)
then indexes the short individual's vector beyond its length, yielding NAs
that F1 converts to fake hets.
Failure scenario: one truncated record in a FASTA silently garbles every
locus; the only signal is a generic matrix warning.
Proposed change: `stop(error(...))` when sequence lengths differ, naming
the offending records.

**F5 [MEDIUM, confidence: high] — reference allele follows the modal genotype, not allele frequency (DOC5, proposed rule)**
`R/utils.read.fasta.r:148-155` — "The allele with the highest frequency is
taken as the reference allele" (gl.read.fasta docs) is implemented as
`which.max` over genotype classes. When the modal class is the
heterozygote the polarity flips: verified on a column Y,Y,Y,T,T,C (allele
T = 7/12) where C becomes ref — T/T codes 2, C/C codes 0, `loc.all` says
`C/T`. The object is internally consistent but documented polarity is
inverted for maf-near-0.5 loci.
Proposed change: choose ref/alt by summed allele counts from the expanded
genotypes.

**F6 [MEDIUM, confidence: high] — self-referential defaults make the signature unusable and the docs false (FS5, DOC5 proposed rule)**
`R/utils.read.fasta.r:22-25` — `parallel = parallel, verbose = verbose`.
A call omitting either argument dies with "promise already under
evaluation" (verified); the roxygen claims `[default FALSE]` and
`[default 4]` (`n.cores`'s actual signature default is NULL). Internal
callers happen to pass everything, which is why it has survived.
Proposed change: `parallel = FALSE`, `verbose = NULL` with
`gl.check.verbosity()` resolution; align the `n.cores` doc with the code.

**F7 [MEDIUM, confidence: high] — the parallel path is broken on its documented default and on Windows (FS5, DEP2)**
`R/utils.read.fasta.r:54-60,101-107` — with `parallel = TRUE` and
`n.cores = NULL` (documented as "maximum number of cores"),
`mc.cores = NULL` errors with "argument is of length zero" on every
platform (verified); with `n.cores = 2` on Windows, fork-based
`parallel::mclapply` errors "'mc.cores' > 1 is not supported on Windows"
(verified). No NULL-resolution (`parallel::detectCores()`) and no platform
guard exist.
Failure scenario: every user who sets `parallel = TRUE` as documented
fails; Windows users fail regardless of `n.cores`.
Proposed change: resolve NULL to `parallel::detectCores()`, and on Windows
either fall back to serial with a gated warning or use a PSOCK cluster.

**F8 [MEDIUM, confidence: high] — results-affecting messages ignore verbosity (VRB3, VRB4 proposed, VRB5)**
`R/utils.read.fasta.r:90,163-166` — the no-polymorphism warning and the
"more than 2 alleles ... skipped" notice are raw `cat(warn())`/
`cat(important())` with no `verbose` gate; both print at `verbose = 0`
(verified through `gl.read.fasta` as well). Under VRB4 these belong at
`verbose >= 1` — they change what the output contains — but they must not
print at 0.
Proposed change: gate both at `verbose >= 1`.

**F9 [MEDIUM, confidence: high] — merge joins duplicate individual names many-to-one silently (DAT2)**
`R/utils.read.fasta.r:261-264` (merge_gl_fasta) — `merge(..., by =
"ind_names", all = TRUE)` treats equal names as the same sample: with a
file containing two records named `ind1`, the second file's single `ind1`
genotypes are copied onto both rows (verified). `gl.compliance.check`
de-duplicates names only afterwards, sealing the wrong join. Within a
single file duplicates pass through untouched (verified, baseline test).
Failure scenario: re-sequenced or sloppily named records silently share
genotypes; with duplicates in both files the join multiplies rows.
Proposed change: error (or warn and uniquify) on duplicate names per file
before merging.

**F10 [LOW, confidence: high] — no-polymorphism path returns cat()'s NULL (FS10)**
`R/utils.read.fasta.r:89-91` — `return(cat(warn(...)))` returns NULL as a
side effect of printing; `merge_gl_fasta` compensates by class-testing for
"NULL" strings. Fragile, and the reason an all-monomorphic single file
crashes downstream in `gl.read.fasta` (see that report, F3).
Proposed change: `return(NULL)` explicitly with the message gated (F8),
and let the caller handle the empty case.

**F11 [LOW, confidence: high] — roxygen gaps (DOC2, DOC7 proposed rules)**
`R/utils.read.fasta.r:8-17` — `@author` has no `Author(s):` line (DOC7);
the `verbose` text deviates from the DOC2 canon; `n.cores` documents a
default the signature does not have (folded into F6's fix); the internal-use
warning required for `utils.*` is present and correct.
Proposed change: align the header; re-document (DOC4).

### @position / @chromosome (PR #330 interaction)

Verified empirically: the constructed genlight carries `@position` and
`@chromosome` NULL. The alignment column of each SNP is encoded only in the
locus name (`<file>_<column>` — the "third semantic" noted in the
position-slot audit), and this is verified correct even after multiallelic
columns are dropped (columns 3,5,6,8 of a 10-column fixture come back as
`clean_3, clean_5, clean_6, clean_8`). No loc.metrics `SnpPosition` is
created, so neither the ddaed27 compliance fill nor the #330 clearing path
engages. Consistent with the genome-only convention; nothing to change
here. A future improvement (out of scope): mirror the alignment column into
a dedicated loc.metrics column rather than the name string.

## Proposed changes

1. Explicit genotype classification with NA for unrecognized cells (F1).
   **Consequence: genotypes at missing/masked sites change from 1 to NA;
   heterozygosity-based results change for any affected dataset.**
2. Pool-based multiallelic detection and skip (F2).
   **Consequence: triallelic columns previously kept are now dropped;
   locus counts change.**
3. `toupper()` normalization on read (F3).
4. Fail fast on unequal sequence lengths (F4).
   **Consequence: previously "successful" reads of ragged files now
   error.**
5. Allele-frequency-based ref/alt choice (F5). **Consequence: dosage
   polarity flips at loci whose modal genotype is the heterozygote.**
6. Real defaults (`parallel = FALSE`, `verbose = NULL`) and doc alignment
   (F6, F11).
7. Working `n.cores` resolution and a Windows-safe parallel path (F7).
8. Gate the two ungated messages at `verbose >= 1` (F8, F10).
9. Duplicate-name guard before the merge join (F9).

Changes 1-5 belong to one redesign of the classification core and should
be reviewed as a unit even though they are listed separately.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run (PLT n/a; FS2-FS4
  n/a by utils convention — no flag start/datatype check expected; stringr
  is in Imports, DEP1 n/a).
- Spec: clean alignment hand-verified per cell (hets R/Y/K/M, ref/alt,
  locus naming, dropped column); missing N and gap cells — run; true
  triallelic and het-minority triallelic — run; lowercase — run; unequal
  lengths — run; wrapped (multi-line) records — run (opaque failure;
  reported under gl.read.fasta F4 as the user-facing surface);
  no-polymorphism file — run; duplicate names within file and across merge
  — run; ref-inversion fixture — run; parallel path — run (both failure
  modes).
- W ambiguity code (A/T het): SKIPPED — same substitution table as the
  verified R/Y/K/M codes; asserted from code reading.
- `parallel = TRUE` on a fork-capable OS: SKIPPED — Windows host; the
  `n.cores = NULL` failure is platform-independent and verified.
- FBM path (DAT6): n/a — this utility always builds an in-memory genlight;
  FBM conversion happens in gl.read.fasta.

## Approval (Phase B)

All BLOCKER, HIGH and MEDIUM findings approved by Arthur Georges
2026-09-05 via the formal approval boxes, explicitly acknowledging the
stated consequences — in particular that utils.read.fasta stops
fabricating genotypes from missing, lowercase and recycled data, so
objects users previously read will differ where the current behaviour
is wrong. LOW findings were not approved this round and are deferred.

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur Georges, 2026-09-05 | F1 (HIGH); consequence acknowledged: genotypes at missing/masked sites change from 1 to NA |
| 2 | approved | Arthur Georges, 2026-09-05 | F2 (HIGH); consequence acknowledged: triallelic columns previously kept are dropped; locus counts change |
| 3 | approved | Arthur Georges, 2026-09-05 | F3 (HIGH) |
| 4 | approved | Arthur Georges, 2026-09-05 | F4 (HIGH); consequence acknowledged: previously "successful" reads of ragged files now error |
| 5 | approved | Arthur Georges, 2026-09-05 | F5 (MEDIUM); consequence acknowledged: dosage polarity flips at het-modal loci |
| 6 | approved | Arthur Georges, 2026-09-05 | F6 (MEDIUM); the F11 (LOW) roxygen items beyond the n.cores default remain deferred |
| 7 | approved | Arthur Georges, 2026-09-05 | F7 (MEDIUM) |
| 8 | approved | Arthur Georges, 2026-09-05 | F8 (MEDIUM); the F10 (LOW) portion is deferred — gating the message necessarily makes the NULL return explicit, but the return VALUE is NULL before and after (no behaviour change) |
| 9 | approved | Arthur Georges, 2026-09-05 | F9 (MEDIUM) |

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-utils.read.fasta` (base upstream/dev
ddaed27): changes 1-9 (F1-F4 HIGH, F5-F9 MEDIUM), with changes 1-5
implemented as the single classification-core redesign the report calls
for: genotypes classify against a per-column allele pool built from
expanded IUPAC codes (NA for anything unrecognized), multiallelism is
detected from that pool, sequences are upper-cased on read, ref/alt
follow summed allele counts (ties resolve alphabetically), and ragged
alignments stop with the offending records named. F10/F11 (LOW)
deferred.

Verification:

- Baseline characterization suite: 32/32 pass. Ten [BUG] pins flipped,
  each mapped to its approved finding ([approved F1]-[approved F9] tags
  in test-utils.read.fasta.R); no unexplained diffs. The dupnames
  single-file passthrough test (not a bug pin) is unchanged and passes.
- E2e engine call at verbose = 3 on clean.fas: the full 6 x 3 genotype
  matrix matches the hand-computed expectation cell for cell, including
  NA at the N and gap cells, loc.all = A/G, C/T, C/T, and skipped
  positions 6 7 reported by alignment column.
- Two-file merge through the (unfixed, this-branch) gl.read.fasta
  wrapper at verbose = 3: nInd 7, nLoc 5, NA fill for absent
  individuals correct.
- Duplicate-name merge: refused with a fatal error naming the records.
- parallel = TRUE with n.cores = NULL and n.cores = 2: identical
  matrices to the serial read (Windows fall-back path).

Integration probe (combined behaviour, this branch's engine + the
`review-gl.read.fasta` wrapper loaded together in a scratch state):

- test-gl.read.fasta.R: 31 pass, 0 fail, 0 skip — the verbose-0
  silence test, which self-gates on the unfixed engine, activates and
  passes with both patches.
- test-utils.read.fasta.R: 32 pass, 0 fail.
- End to end at default verbosity, clean.fas + clean2.fas: nInd 7,
  nLoc 5, genotypes present in @gen, N and gap cells NA, triallelic
  columns absent, absent individuals NA at the other file's loci.
- verbose = 0 through the full public surface: zero captured lines.

```json
{
  "function": "utils.read.fasta",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT1", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DAT1", "status": "applied", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "DAT1", "status": "applied", "change": 3},
    {"id": "F4", "severity": "HIGH", "confidence": "high", "rule": "DAT2", "status": "applied", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 5},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "applied", "change": 6},
    {"id": "F7", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "applied", "change": 7},
    {"id": "F8", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "applied", "change": 8},
    {"id": "F9", "severity": "MEDIUM", "confidence": "high", "rule": "DAT2", "status": "applied", "change": 9},
    {"id": "F10", "severity": "LOW", "confidence": "high", "rule": "FS10", "status": "deferred", "change": 8},
    {"id": "F11", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "deferred", "change": 6}
  ],
  "coverage_skipped": ["W ambiguity code: same table as verified codes", "fork-parallel on non-Windows: no host available"],
  "status": "phase-c-applied",
  "pr": null
}
```
