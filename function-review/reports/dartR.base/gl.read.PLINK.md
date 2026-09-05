# Review: gl.read.PLINK (dartR.base)

- Family mode: io (reviewed as a chain with `utils.read.ped`; note the chain
  is nominal — `gl.read.PLINK` reads via `snpStats::read.plink`, and
  `utils.read.ped`'s only caller is `gl.report.ld.map`)
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ddaed27 (upstream/dev; file verified identical on local HEAD ed99203)
- Datasets: synthetic PLINK fileset (4 ind x 4 SNP `toy.ped`/`toy.map`
  hand-built; `toy.bed`/`.bim`/`.fam` generated with PLINK v1.90b7.2),
  scrambled ind.metafile CSV, loc.metafile CSVs — see Fixtures
- Baseline: tests/testthat/test-gl.read.PLINK.R (snapshot captured pre-review;
  the .bed bytes are embedded so the tests need no plink executable)

## Verdict

**Standards: Needs work** — the preamble (FS2, FS3, DEP1), double
`gl.compliance.check`, and history entry are in place, but metadata joins
and error exits break DAT2/FS5, and the roxygen header has no examples.
**Spec: Rework** — at any verbosity <= 2, including the default, the
returned object contains no genotype data at all; the `fbm` argument is
dead code and verbosity decides the storage backend instead.

What works well: the `.bim` genome coordinates land in `@position`/
`@chromosome` exactly as the genome-only convention (PR #330) prescribes
(verified: positions 100,200,150,300 and chromosomes 1,1,2,2 on the
fixture), and a valid loc.metafile is realigned by `AlleleID` row name, so
scrambled locus metadata merges correctly.

## Findings

**F1 [BLOCKER, confidence: high] — returned object has no genotypes at default verbosity; `fbm` argument is dead (DAT1, VRB1)**
`R/gl.read.PLINK.r:316-321` — the closing block reads `if (fbm) {}` (empty
braces), then unconditionally `gl <- gl.gen2fbm(gl, verbose = verbose)`,
which moves `@gen` into `@fbm` and empties `@gen`; then
`if (verbose>2) {...} else gl@fbm <- NULL` discards the FBM copy whenever
verbose <= 2. Net behaviour, verified on the fixture: at `verbose = 0`, 1
or 2 (2 is the default) the object comes back with `length(@gen) == 0`,
`@fbm NULL`, `nInd() == 0` and `as.matrix()` erroring — total genotype
loss; at `verbose = 3` the object is FBM-backed even with `fbm = FALSE`.
Verbosity, not the `fbm` flag, decides the result, and no setting yields
the ordinary `@gen`-backed object the docs describe.
Failure scenario: every default-settings call — `gl.read.PLINK("data")` —
returns an empty shell; any downstream use fails or silently reports 0
individuals.
Proposed change: `if (fbm) gl <- gl.gen2fbm(gl, verbose = verbose)` and
delete the verbose-conditioned nulling; gate the "Created a file-backed
matrix" report on `fbm && verbose >= 2`.

**F2 [HIGH, confidence: high] — ind.metafile rows joined by position, not id; pop column never applied (DAT2)**
`R/gl.read.PLINK.r:282-291` — `gl@other$ind.metrics <- cbind(gl@other$ind.metrics, ind.metrics)`
binds the CSV rows in file order against individuals in `.fam` order. The
preceding check only tests membership of ids, not order. Verified: with a
CSV listing ind3,ind1,ind4,ind2, individual `ind1` receives ind3's `pop`
and `tag` values. Additionally the CSV's `pop` column lands only in
`ind.metrics`; `pop(gl)` keeps the hardcoded "A" for every individual.
Failure scenario: any metafile whose row order differs from the `.fam` —
the normal case for files assembled by hand or exported from a database —
silently assigns every metadata field, including population, to the wrong
animals.
Proposed change: reorder before binding
(`ind.metrics <- ind.metrics[match(fam$member, ind.metrics$id), ]`) and set
`pop(gl) <- ind.metrics$pop` when present.

**F3 [MEDIUM, confidence: high] — missing AlleleID column prints a fatal error but does not stop (FS5)**
`R/gl.read.PLINK.r:211-218` — the mandatory-column check is
`cat(error(...))` with no `stop()`, so execution continues into
`loc.metrics[, "AlleleID"]`, which dies with "undefined columns selected"
(verified). The compliant genotype checks at lines 181-195 use the sibling
idiom `cat(error(...)); stop()` — non-standard but at least terminal.
Failure scenario: a user with a headerless or misnamed locus metafile gets
an unrelated subscript error after the real diagnosis has scrolled past.
Proposed change: use `stop(error("..."))` in all three places.

**F4 [MEDIUM, confidence: high] — dosage orientation undocumented and opposite to gl.read.vcf (DOC5, proposed rule)**
`R/gl.read.PLINK.r:144-155` — the snpStats numeric coding counts
`allele.2` (`.bim` column 6, the major allele under PLINK 1.x defaults), so
2 = homozygous major and 0 = homozygous minor. Verified on the fixture:
ind1 `A A` at snp1 (bim `G A`) reads 2. `gl.read.vcf` counts the ALT
allele, so the two importers give complementary dosages for the same
variant, and neither roxygen header says which allele is counted.
Failure scenario: allele-frequency or association work mixing PLINK- and
VCF-imported datasets flips ref/alt semantics with no warning.
Proposed change: state the counted allele in `@details` (and in
`gl.read.vcf`'s); if the team wants ALT/minor-consistent dosage, that is an
API1 change. **Consequence: numerical output changes for every locus if
re-oriented.**

**F5 [MEDIUM, confidence: high] — .ped path spams the console at any verbosity and writes into the input directory (VRB5, FS7)**
`R/gl.read.PLINK.r:117-134` and `R/utils.plink.run.r:59-60` — the PLINK
conversion runs via `system(cmd)` with no output suppression; the child
process's progress chatter prints even at `verbose = 0` (verified — it
bypasses `sink()`), and the generated `toy.bed/.bim/.fam/.log` are written
into `dirname(filename)`, the user's data directory, although `dir.out <-
tempdir()` is created at line 88 with a comment saying transformations
should happen there.
Failure scenario: `verbose = 0` batch jobs still emit hundreds of progress
lines; read-only or shared input directories fail or accumulate artefacts.
Proposed change: convert in `dir.out` (copy the `.ped`/`.map` or point
plink's `--out` there) and call `system(cmd, ignore.stdout = verbose < 3)`.

**F6 [MEDIUM, confidence: high] — roxygen header gaps (DOC1, DOC2, DOC5, DOC7)**
`R/gl.read.PLINK.r:1-39` — no `@examples` at all (DOC1 required tag; TST1);
`@param fbm` says the matrix lands "in its @genome slot" but the slot is
`@fbm` and, per F1, no argument value currently produces it (DOC5);
`@param verbose` deviates from the DOC2 default clause; `@author` has a
`Custodian:` line but no `Author(s):` (DOC7, proposed rule); no `@family`.
Failure scenario: `example(gl.read.PLINK)` gives nothing; users hunt for a
nonexistent `@genome` slot.
Proposed change: add a `\donttest{}` example built on a shipped or
generated fileset, correct the `fbm` text, standardize verbose wording, add
`Author(s):`, re-document (DOC4).

**F7 [LOW, confidence: high] — verbose-2 report prints the wrong slot name and fires per column (VRB2)**
`R/gl.read.PLINK.r:241-249` — the loc.metrics update message says the names
were added "to the other$ind.metrics slot" and, being `paste()` over the
name vector inside one `cat(report(...))`, prints one line per column.
Proposed change: name the correct slot; collapse the vector.

**F8 [LOW, confidence: medium] — FLAG SCRIPT END precedes the FBM stage (FS9)**
`R/gl.read.PLINK.r:310-321` — "Completed:" prints before `gl.gen2fbm`
runs, so a failure in the FBM stage arrives after the completion flag.
Fixing F1 should move the fbm block above the flag.

## Proposed changes

1. Wire the `fbm` flag: convert only when `fbm = TRUE`, never null the
   result by verbosity, and return a `@gen`-backed object otherwise (F1, F8).
   **Consequence: numerical output changes for every current caller — the
   returned object regains its genotypes.**
2. Match ind.metafile rows to `fam$member` by id before binding, and apply
   the `pop` column to `pop(gl)` (F2). **Consequence: metadata values (and
   populations) change for callers whose CSV order differed from the .fam.**
3. Replace the three `cat(error(...))` exits with `stop(error(...))` (F3).
4. Document the counted allele in `@details`; put the PLINK-vs-VCF
   orientation question to the team (F4).
5. Run the .ped->.bed conversion in `tempdir()` and suppress plink stdout
   below verbose 3 (F5).
6. Roxygen repairs: examples, `fbm` slot name, DOC2 verbose text,
   `Author(s):` line, `@family`; re-document (F6, F7).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run (no plot bundle;
  DEP1 guard present for snpStats, which is in Imports so the guard is
  belt-and-braces; FS4 datatype check not applicable — input is a file
  path).
- Spec: .bed path on the synthetic fileset — run (genotypes hand-verified
  cell-by-cell at verbose = 3 via `@fbm`; missing `0 0` genotype comes
  through as NA); .ped path with a real PLINK 1.9 binary — run (same
  genotypes; console/side-effect behaviour under F5); ind.metafile scramble
  — run; loc.metafile scramble and missing-column — run; `verbose = 0`
  silence on the .bed path — verified silent (0 captured lines; plot check
  n/a, no plots).
- Genome-only @position convention (PR #330): verified — `.bim` bp
  positions and chromosomes land in the slots; no tag-offset material
  exists for PLINK input, so `loc.metrics$SnpPosition` is legitimately
  absent (created NA-filled by `gl.compliance.check`).
- FBM path (DAT6): partially run — the always-on `gl.gen2fbm` call was
  exercised (F1); no large-data densification test attempted.
- Shipped fixtures: none — dartR.base has no `inst/extdata` PLINK files
  and dartR.data ships only .rda genlights; all fixtures are synthetic and
  documented in the baseline test file.
- Duplicated-id inputs (duplicate `fam$member` / `snp.name` stops):
  SKIPPED — asserted from code reading only (lines 181-195).

## Fixtures (all synthetic, written by the probes/tests)

- `toy.ped`/`toy.map`: 4 individuals (2 families, both sexes, case/control
  phenotype), 4 SNPs incl. one missing genotype (`0 0`); full truth table
  in the test file header.
- `toy.bed`/`toy.bim`/`toy.fam`: produced by PLINK v1.90b7.2 from the
  above; the 7-byte `.bed` is embedded in the baseline test.
- `indmeta.csv` (ids deliberately out of .fam order), `locmeta.csv`
  (valid, scrambled) and a no-AlleleID variant.

## Approval (Phase B)

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved (F1; F8 deferred) | Arthur Georges, 2026-09-05 | Consequence acknowledged: default-settings callers regain genotypes; the FLAG SCRIPT END position (F8, LOW) is untouched this round. |
| 2 | Approved | Arthur Georges, 2026-09-05 | Consequence acknowledged: metadata values and populations change for callers whose CSV order differed from the .fam — to the correct values. |
| 3 | Approved | Arthur Georges, 2026-09-05 | |
| 4 | Approved as DOCUMENTATION | Arthur Georges, 2026-09-05 | Orientation recorded in both readers' roxygen; no dosage flip. |
| 5 | Approved | Arthur Georges, 2026-09-05 | |
| 6 | Approved (F6; F7 deferred) | Arthur Georges, 2026-09-05 | |

All BLOCKER/HIGH/MEDIUM findings approved 2026-09-05, with explicit
acknowledgement that objects users previously read will differ where
current behaviour is wrong. LOW findings (F7, F8) were not approved this
round — deferred, not applied.

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-gl.read.PLINK` (base upstream/dev
ddaed27):

- F1 (BLOCKER): the `fbm` flag now decides the backend —
  `if (fbm) gl <- gl.gen2fbm(...)`; the verbose-conditioned nulling is
  deleted. Default calls return a `@gen`-backed object with genotypes.
- F2 (HIGH): ind.metafile rows matched to `fam$member` by `id` before
  binding; `pop` column applied to `pop(gl)`.
- F3 (MEDIUM): the three `cat(error(...))` exits (duplicate individual
  labels, duplicate AlleleID, missing AlleleID column) replaced with
  `stop(error(...))`.
- F4 (MEDIUM, documentation only): counted allele (allele.2, the PLINK
  1.x major) stated in `@details`, with the cross-reference to
  `gl.read.vcf`'s ALT orientation.
- F5 (MEDIUM): .ped-to-.bed conversion copies the fileset to `tempdir()`
  and runs there; `utils.plink.run` gains
  `system(cmd, ignore.stdout = verbose < 3)`.
- F6 (MEDIUM): roxygen — `@family io`, runnable `\donttest{}` example
  (embedded 4x4 .bed fileset), corrected `fbm` slot name (`@fbm`), DOC2
  verbose text, `Author(s):` line; re-documented.
- F7, F8 (LOW): deferred — not applied.

Verification: baseline characterization test updated; all 26 assertions
pass; every diff from the pre-review snapshot is annotated
`[approved F1]`/`[approved F2]`/`[approved F3]`; LOW-pinned snapshots
(loc.metafile AlleleID merge, verbose-0 silence, dosage orientation,
coordinate slots) unchanged. E2e at verbose = 3: .bed path genotypes
hand-verified against the truth table; .ped path run with PLINK
v1.90b7.2 — zero files added to the input directory and zero console
lines at verbose = 0.

Caller grep (all 8 clones): no callers of `gl.read.PLINK` anywhere in
the dartRverse — all clear.

```json
{
  "function": "gl.read.PLINK",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ddaed27",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "BLOCKER", "confidence": "high", "rule": "DAT1", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DAT2", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "applied-as-documentation", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "applied", "change": 5},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "deferred", "change": 6},
    {"id": "F8", "severity": "LOW", "confidence": "medium", "rule": "FS9", "status": "deferred", "change": 1}
  ],
  "coverage_skipped": ["duplicated-id stop paths: asserted from code reading", "DAT6 large-data FBM: not attempted"],
  "status": "phase-c-complete",
  "pr": null
}
```
