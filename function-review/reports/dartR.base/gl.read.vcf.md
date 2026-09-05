# Review: gl.read.vcf (dartR.base)

- Family mode: io (reviewed as a chain with `utils.vcfr2genlight.polyploid`,
  its genotype-conversion engine)
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ddaed27 (upstream/dev; file verified identical on local HEAD ed99203)
- Datasets: eight synthetic VCF 4.2 fixtures (see Fixtures), scrambled
  ind.metafile CSV
- Baseline: tests/testthat/test-gl.read.vcf.R (snapshot captured pre-review,
  bugs included)

## Verdict

**Standards: Needs work** — the preamble is present but the function is
never silent at `verbose = 0`, prints no FLAG SCRIPT END of its own, and
the header misses required tags.
**Spec: Needs work** — the core diploid biallelic path is correct
(genotypes, phasing, missing data, coordinates and QUAL/FILTER/INFO all
hand-verified), but dosage mode ships polyploid dosages stamped diploid,
INFO metadata can be silently misassigned, and FILTER is destroyed whenever
a multi-allelic record is present.

What works well: ALT-allele dosage, phased separators, and `./.` -> NA all
verified cell-by-cell; POS/CHROM land in `@position` (integer) and
`@chromosome` (factor) per the genome-only convention (PR #330); the
multi-allelic drop is applied consistently to genotypes, slots and
loc.metrics (probed for desync with an `ALT=.` monomorphic record and with
comma-ALT records — lengths stayed in step); the ind.metafile join matches
by id and realigns correctly.

## Findings

**F1 [HIGH, confidence: high] — dosage mode returns dosages up to ploidy-n under a forced ploidy of 2 (DAT1)**
`R/gl.read.vcf.r:125` — after the helper builds the matrix (correctly:
`1/1/1` -> 3, `0/0/1` -> 1, hand-verified), `ploidy(x) <- rep(2, nInd(x))`
stamps every individual diploid; the comment above ("allow varied ploidy
level") documents the opposite intent, and the helper's roxygen says dosage
mode assigns ploidy from copy number. Verified: the triploid fixture
returns a matrix containing 3s with `ploidy(x) == c(2,2,2)`, and
`gl.compliance.check` does not object.
Failure scenario: every downstream statistic that assumes dosage in 0..2
(maf, heterozygosity, distances) computes garbage on dosage-mode imports;
nothing warns.
Proposed change: keep the helper's inferred ploidy in dosage mode (or set
it from max copy number per the docs); force 2 only in genotype mode.
**Consequence: ploidy values change for dosage-mode callers.**

**F2 [HIGH, confidence: high] — INFO fields misassigned when key order or content varies across records (DAT2)**
`R/gl.read.vcf.r:74-92` — the INFO string is split on `=|;` per record and
`rbind`ed; column names are taken from record 1 and values are filled by
position. Only a count mismatch triggers the fallback. Verified: with
record 1 `DP=10;AC=2` and record 2 `AC=3;DP=99`, the returned loc.metrics
give record 2 `DP = 3, AC = 99` — values swapped, silently. Flag-type INFO
keys (no `=`) shift pairing the same way.
Failure scenario: any VCF whose INFO keys vary in order or presence across
records — common in merged or multi-caller VCFs — gets plausible but wrong
locus metrics; filtering on rdepth/AC then filters the wrong loci.
Proposed change: parse key=value pairs per record into named columns
(e.g. build a named list per record and fill by key), falling back to
QUAL/FILTER when keys are absent.

**F3 [MEDIUM, confidence: high] — never silent at verbose = 0 (VRB5)**
`R/gl.read.vcf.r:127,274` — `gl.compliance.check(x)` and
`gl.recalc.metrics(x)` are called without the user's `verbose`, so both run
at their own default (2). Verified: `verbose = 0` on the basic fixture
prints 48 lines. The ind.metafile mismatch warnings (lines 162-180) are
likewise ungated.
Failure scenario: batch pipelines running quietly still emit a full
compliance transcript per file.
Proposed change: pass `verbose = verbose` to both calls (VRB4 argues the
id-mismatch warnings should print from `verbose >= 1`).

**F4 [MEDIUM, confidence: high] — FILTER (and any character INFO column) destroyed when a multi-allelic record is present (DAT2)**
`R/gl.read.vcf.r:98-99` — in the AC branch, `info[] <- lapply(info,
as.numeric)` runs only when multi-allelic records were found, coercing
every column; FILTER "PASS" becomes NA. Verified: identical fixture with
and without a multi-allelic record returns FILTER "PASS" vs FILTER NA.
Failure scenario: presence of one unrelated multi-allelic record silently
erases the FILTER column (and any character INFO field) for every retained
locus.
Proposed change: numerify only columns that are numeric-parseable (or
QUAL/DP-style known-numeric keys), in both branches or neither.

**F5 [MEDIUM, confidence: high] — individuals absent from the ind.metafile are dropped from the returned object (DAT2, DOC5)**
`R/gl.read.vcf.r:183-209` — the id join subsets `x <- x[ord2, ]` to the
matched individuals. Verified: a metafile listing S3,S1 returns a 2-
individual object with S2 gone; the notice is an ungated `cat(warn(...))`
listing unmatched ids, and `@details` does not say individuals will be
removed.
Failure scenario: a metafile missing one animal silently shrinks the
dataset; downstream sample sizes drop without an explicit decision.
Proposed change: keep unmatched individuals with NA metadata (or stop),
and document whichever behaviour is chosen; at minimum print the drop at
`verbose >= 1` as a consequence-bearing warning (VRB4).

**F6 [MEDIUM, confidence: medium] — haploid and partial calls silently recoded as diploid (DAT1)**
`R/gl.read.vcf.r:4,125` — `@param vcffile` says "works only for diploid
data" but nothing enforces it. Verified: a haploid VCF returns ploidy-2
homozygotes (`1` -> 2); half-missing diploid calls are mis-coded by the
helper (`0/.` -> 1; see the chained review F1 there). No warning either way.
Failure scenario: mitochondrial/sex-chromosome or low-confidence VCFs
import as confident diploid genotypes.
Proposed change: detect non-diploid GT arity and stop (or warn at
`verbose >= 1`) per the documented diploid-only contract.

**F7 [LOW, confidence: high] — DEP1 guard returns -1; no FLAG SCRIPT END (DEP1, FS9)**
`R/gl.read.vcf.r:55-62` — missing vcfR prints via `cat(error(...))` and
`return(-1)` instead of `stop(error(...))`; and the function never prints
`Completed: gl.read.vcf` (verified: 49 lines at verbose 2, none its own
end-flag — the last "Completed" is `gl.recalc.metrics`'s).
Proposed change: `stop(error(...))` in the guard; add FS9 before `return`.

**F8 [LOW, confidence: high] — misleading tail messages (VRB2, DOC5)**
`R/gl.read.vcf.r:276-290` — at `verbose > 2` the function prints
"Genlight object does not have individual metrics" even when a metafile was
just ingested, and "Created an file-backed matrix (fbm) dartR object" even
when `fbm = FALSE` (the message sits outside the `if (fbm)` block).
Proposed change: gate the first on `is.null(ind.metafile)`, move the second
inside the fbm block.

**F9 [LOW, confidence: high] — header misses required tags; examples unrunnable (DOC1, DOC3, DOC7)**
`R/gl.read.vcf.r:1-35` — no `@name`/`@title`/`@family`; `@examples` is
`\dontrun{}` referencing `system.file('extdata/test.vcf', package='dartR')`
(wrong package, file not shipped); `@param mode`/`@details` text is garbled
("sites will be codes as", "will issue the user with an error");
`@author` names authors with no `Custodian:` label (DOC7, proposed rule);
`@param fbm` lacks a `[default FALSE]` tag.
Proposed change: complete the header, ship a small `inst/extdata` VCF so
the example runs, add the Custodian line; re-document (DOC4).

**F10 [INFO, confidence: high] — monomorphic (`ALT=.`) records are retained with `loc.all` "C/NA" (DOC5)**
`R/gl.read.vcf.r:69` — `paste0(REF, "/", ALT)` on a missing ALT yields the
literal string "C/NA"; the locus imports as all-0. Indels (`A/AT`) are
likewise retained silently and coded as ALT dosage. Neither behaviour is
documented. Recorded for the docs pass; no data corruption observed.

Compliance-check interaction (not a finding against this function): the
first loc.metrics column of VCF-imported objects is named
`array(NA, nLoc(x))` — an unnamed-expression cbind artefact arising inside
`gl.compliance.check`, noted here for the campaign but belonging to that
function's review (PR #330 territory adjacent; no hunks proposed).

## Proposed changes

1. Stop forcing ploidy 2 in dosage mode; set ploidy per the documented
   copy-number semantics (F1). **Consequence: ploidy (and all
   ploidy-dependent statistics) change for dosage-mode callers.**
2. Parse INFO by key, not position (F2). **Consequence: loc.metrics values
   change for VCFs with variable INFO layouts — to the correct values.**
3. Pass `verbose` through to `gl.compliance.check` and
   `gl.recalc.metrics`; gate the metafile warnings (F3).
4. Restrict the AC-branch numeric coercion to numeric-parseable columns
   (F4).
5. Decide and document the unmatched-individual policy; stop dropping
   silently (F5). **Consequence: returned dimensions change for callers
   with partial metafiles if retention is chosen.**
6. Enforce or document the diploid-only contract for haploid input (F6).
7. DEP1 stop, FS9 end-flag, tail-message gating (F7, F8).
8. Complete the roxygen header and make the example runnable (F9, F10).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run (no plot bundle;
  FS4 n/a — input is a file path; DAT5 satisfied via the internal
  `gl.compliance.check`).
- Spec: biallelic diploid path hand-verified cell-by-cell (incl. phased
  `|`, `./.` -> NA, QUAL/FILTER/DP/AC landing); multi-allelic and indel
  records — run; `ALT=.` monomorphic desync probe — run (no desync found;
  hypothesis refuted); INFO key-order variation — run (F2); dosage mode on
  triploid GT incl. `0/0/1` — run, hand-computed match (F1); haploid — run
  (F6); ind.metafile scramble/subset — run (F5); gzipped `.vcf.gz` — run
  (reads correctly; undocumented); `verbose = 0` silence — run (F3; text
  side only, no plot path exists).
- FBM path (DAT6): SKIPPED — `fbm = TRUE` exercises `gl.gen2fbm`, which the
  gl.read.PLINK review covers; no large-data test attempted here.
- Shipped fixtures: none found — no `inst/extdata` in dartR.base, none in
  dartR.data; all fixtures synthetic and documented below.

## Fixtures (all synthetic VCF 4.2, written by the probes/tests)

- `basic.vcf`: 3 samples x 3 biallelic SNPs; contig headers; DP/AC INFO;
  phased and unphased separators; one `./.`.
- `acmulti.vcf`: adds `ALT=G,T` multi-allelic and an `A/AT` indel.
- `mono.vcf`: `ALT=.` monomorphic record, no AC key (desync probe).
- `infoswap.vcf`: INFO key order swapped between records.
- `poly.vcf`: triploid GTs incl. `0/0/1`, `1/1/1`, `././.`.
- `haploid.vcf`: single-allele GTs.
- `partial.vcf`: `./1` and `0/.` half-missing calls (helper review).
- `basic.vcf.gz`: gzip of basic.vcf.
- `indmeta.csv`: scrambled two-row subset (S3, S1) with pop/age.

## Approval (Phase B)

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges, 2026-09-05 | Consequence acknowledged: ploidy (and ploidy-dependent statistics) change for dosage-mode callers — stamped correctly. |
| 2 | Approved | Arthur Georges, 2026-09-05 | Consequence acknowledged: loc.metrics values change for VCFs with variable INFO layouts — to the correct values. |
| 3 | Approved | Arthur Georges, 2026-09-05 | |
| 4 | Approved | Arthur Georges, 2026-09-05 | |
| 5 | Approved | Arthur Georges, 2026-09-05 | Consequence acknowledged: returned dimensions change for callers with partial metafiles — retention with NA metadata chosen. |
| 6 | Approved | Arthur Georges, 2026-09-05 | Warning at verbose >= 1; recoding documented. |
| 7 | Deferred | Arthur Georges, 2026-09-05 | LOW findings F7, F8 not approved this round. |
| 8 | Deferred (except the orientation sentence) | Arthur Georges, 2026-09-05 | The dosage-orientation documentation (chained gl.read.PLINK F4) IS applied; the rest of the F9 header repair and F10 are deferred/no action. |

All HIGH/MEDIUM findings approved 2026-09-05, with explicit
acknowledgement that objects users previously read will differ where
current behaviour is wrong (ploidy stamped correctly, INFO parsed by
key, id-matched joins retained). The MEDIUM dosage-orientation item from
the chained gl.read.PLINK review (its F4) is applied as DOCUMENTATION in
this function's roxygen. LOW findings (F7, F8, F9) deferred; F10 (INFO)
no action.

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-gl.read.vcf` (base upstream/dev
ddaed27):

- F1 (HIGH): ploidy 2 is forced only in genotype mode; dosage mode
  stamps the data's maximum ALT copy number (uniform across
  individuals — `gl.compliance.check` errors on mixed per-individual
  ploidy, discovered during verification). The triploid fixture now
  returns ploidy 3,3,3 (was 2,2,2 under dosages of 3).
- F2 (HIGH): INFO parsed per record into key=value pairs and filled by
  key; flag-type keys (no `=`) become "TRUE" columns; keys absent from
  a record are NA. The count-mismatch fallback that discarded all INFO
  is gone.
- F3 (MEDIUM): `verbose` passed to `gl.compliance.check` and
  `gl.recalc.metrics`; metafile id-mismatch warnings gated at
  `verbose >= 1`. `verbose = 0` is now fully silent (was 48 lines).
- F4 (MEDIUM): the AC-branch numeric coercion is restricted to
  numeric-parseable columns; FILTER survives multi-allelic files.
- F5 (MEDIUM): individuals absent from the ind.metafile are retained
  with NA metadata (id filled from the object); the retention is listed
  in a warning at `verbose >= 1` and documented in `@details`.
- F6 (MEDIUM): non-diploid GT arity in genotype mode triggers a warning
  at `verbose >= 1`; the diploid recoding is documented in `@details`.
- Chained gl.read.PLINK F4 (documentation): ALT-count orientation stated
  in `@details` with the cross-reference to `gl.read.PLINK`.
- F7, F8, F9 (LOW): deferred. F10 (INFO): no action.

Verification: baseline characterization test updated; all 39 assertions
pass; every diff from the pre-review snapshot is annotated
`[approved F1..F6]`. Unchanged (LOW-pinned) snapshots: no FLAG SCRIPT
END at verbose 2 (F7), indel retention, gzip path, biallelic dosage
matrix and coordinate slots. E2e hand-verified against the raw files:
basic (dosages, QUAL/FILTER/DP/AC), infoswap (DP=99/AC=3 by key),
flag-key INFO (DR2/AF/IMP incl. NA for the absent key), mono `ALT=.`
(no desync, loc.metrics in step), triploid dosage (matrix and ploidy
3,3,3), haploid (recoding warning at verbose >= 1). Half-missing calls
are unchanged on this branch — their fix is the chained
`utils.vcfr2genlight.polyploid` PR.

Caller grep (all 8 clones): one in-package caller, `gl.impute`
(beagle-imputed diploid VCF, `verbose = 0`, no metafile). Effect: its
read is now silent as requested; beagle's flag-type `IMP` INFO key
previously triggered the count-mismatch fallback (all INFO discarded) —
it now yields DR2/AF/IMP columns; genotype use is untouched. Not
breaking — all clear.

```json
{
  "function": "gl.read.vcf",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ddaed27",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT1", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DAT2", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DAT2", "status": "applied", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DAT2", "status": "applied", "change": 5},
    {"id": "F6", "severity": "MEDIUM", "confidence": "medium", "rule": "DAT1", "status": "applied", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DEP1", "status": "deferred", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "deferred", "change": 7},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "deferred", "change": 8},
    {"id": "F10", "severity": "INFO", "confidence": "high", "rule": "DOC5", "status": "no-action", "change": 8}
  ],
  "coverage_skipped": ["DAT6 FBM path: covered under gl.read.PLINK review; no large-data test"],
  "status": "phase-c-complete",
  "pr": null
}
```
