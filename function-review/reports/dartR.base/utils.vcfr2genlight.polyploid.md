# Review: utils.vcfr2genlight.polyploid (dartR.base)

- Family mode: io (chained with `gl.read.vcf`, its only in-package caller)
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ddaed27 (upstream/dev; file verified identical on local HEAD ed99203)
- Datasets: synthetic VCF 4.2 fixtures read via `vcfR::read.vcfR` (see the
  gl.read.vcf report's fixture list; `poly.vcf` and `partial.vcf` are the
  load-bearing ones here)
- Baseline: tests/testthat/test-utils.vcfr2genlight.polyploid.R (snapshot
  captured pre-review, bugs included)

## Verdict

**Standards: Needs work** — the utils warning banner is present, but
dependency guards, error idioms and the details text need repair.
**Spec: Needs work** — the dosage arithmetic for well-formed polyploid
calls is correct (hand-verified, `0/0/1` -> 1 and `1/1/1` -> 3 included),
but genotype mode invents heterozygotes from half-missing calls, and the
ploidy the object carries is an artefact of per-individual max dosage.

What works well: multi-allelic records are dropped with an explicit
counting warning before conversion (verified), CHROM/POS/ID land in
`@chromosome`/`@position`/`locNames` (genome coordinates, PR #330
conform), and both documented modes produce the hand-computed matrix for
complete calls.

## Findings

**F1 [HIGH, confidence: high] — genotype mode codes half-missing calls as heterozygous (DAT1)**
`R/utils.vcfr2genlight.polyploid.R:54-57` — after separators are stripped,
the genotype-mode rule is "anything not all-0, not all-1, and longer than
one character becomes 1". A half-missing diploid call survives as `".1"` or
`"0."` and falls into that bucket. Verified: `0/.` — a call with no ALT
allele observed — returns 1 (heterozygous); `./1` also returns 1. Fully
missing `./.` is safe (vcfR delivers it as NA, and NA survives the
length-1 logical assignments). Dosage mode differs: `./1` -> 1 and `0/.`
-> 0, the latter only because `as.numeric("0.")` happens to parse as 0
when adegenet coerces the character matrix — accidental, not designed.
Failure scenario: low-coverage VCFs, where half-calls are routine, inflate
heterozygosity in genotype mode with no warning; the same file gives
different answers in the two modes.
Proposed change: convert tokens still containing "." to NA before the
mode rules (e.g. `x[grepl(".", x, fixed = TRUE)] <- NA`).

**F2 [MEDIUM, confidence: high] — ploidy is per-individual max dosage, not the data's ploidy (DAT1, DOC5)**
`R/utils.vcfr2genlight.polyploid.R:77` — `new("genlight", t(x))` lets
adegenet infer ploidy from each individual's maximum value. Verified on
the triploid fixture: dosage mode returns ploidy 2,3,3 (the first triploid
individual carries no `1/1/1` call, so it is stamped diploid); genotype
mode returns ploidy 1,2,2 (an individual whose calls are all-ref or coded
1 is stamped haploid). The `@details` wording ("assign ploidy levels as
maximum copy number of alternate alleles") describes the mechanism but not
its wrongness as a ploidy estimate; the caller then overwrites all of it
with 2 (see gl.read.vcf F1).
Failure scenario: mixed-ploidy stamps on uniform-ploidy data; every
ploidy-aware downstream computation misinterprets the object.
Proposed change: derive ploidy from GT arity (allele count per call,
e.g. from the pre-strip token length), which is available and exact.

**F3 [MEDIUM, confidence: high] — omitted mode2 fails with a closure-coercion error (FS5)**
`R/utils.vcfr2genlight.polyploid.R:34` — the signature default
`mode2 = mode` resolves to `base::mode` (a function). Verified:
`utils.vcfr2genlight.polyploid(x = v)` fails with "comparison (==) is
possible only for atomic and list types". The `@details` claim that a
missing mode "will issued the user with a error" documents this accident
as a feature.
Failure scenario: any direct caller omitting `mode2` gets an error naming
neither the argument nor the fix.
Proposed change: default `mode2 = "genotype"` (matching the caller) or a
`match.arg` check with an informative stop.

**F4 [LOW, confidence: high] — invalid mode2 stops with an empty condition message (VRB2)**
`R/utils.vcfr2genlight.polyploid.R:71-74` — `cat(error(...))` then bare
`stop()`. Verified: `conditionMessage` of the thrown error is `""`, so
`tryCatch` consumers and logs see a blank error.
Proposed change: `stop(error("  Please choose 'genotype' or 'dosage' mode\n"))`.

**F5 [LOW, confidence: high] — dependency guards absent or inverted (DEP1)**
`R/utils.vcfr2genlight.polyploid.R:35,76-80` — `vcfR::` and `stringr::`
are used with no DEP1 guard (vcfR is only in Suggests; the exported helper
can be called without the guarded caller), while adegenet — a hard
dependency — gets `if (requireNamespace('adegenet'))` whose else-branch
warns "adegenet not installed" and then proceeds to call
`adegenet::chromosome()` on the raw character matrix, guaranteeing a
downstream crash.
Proposed change: add the DEP1 stop for vcfR; drop the adegenet
conditional (or make it stop).

**F6 [LOW, confidence: high] — details text garbled; example not runnable (DOC5, DOC3)**
`R/utils.vcfr2genlight.polyploid.R:10-21` — "missing input of mode will
issued the user with a error", "codes were modified from 'vcfR2genlight'
in vcfR packge"; the `\dontrun` example references an undefined `vcfr`
object and assigns to a variable named `datatype`. `@author` has the
Custodian line (DOC7 satisfied on the custodian side; no `Author(s):`
label).
Proposed change: rewrite details, provide a self-contained example,
add `Author(s):`; re-document (DOC4).

## Proposed changes

1. Treat any residual "."-bearing token as NA before mode coding (F1).
   **Consequence: numerical output changes — half-missing calls become NA
   instead of 1 (genotype mode) or 0/1 (dosage mode).**
2. Derive ploidy from GT arity instead of max dosage (F2).
   **Consequence: ploidy values change for polyploid imports.**
3. Give `mode2` a real default or `match.arg` (F3).
4. `stop(error(...))` for the invalid-mode branch (F4).
5. Fix the dependency guards (F5).
6. Documentation repair (F6).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, STY — run (PLT n/a: no plots;
  FS2/FS3/FS9 absent as is usual for utils helpers — accepted; the
  function prints nothing except its error branch, so verbosity gating is
  n/a).
- Spec: triploid dosage and genotype matrices hand-verified per cell
  (incl. the assigned `0/0/1` -> 1 killer case); ploidy inference — run
  (F2); half-missing calls both modes — run (F1); fully-missing `./.` and
  `././.` — run (NA preserved); multi-allelic drop — run; haploid GTs —
  run via the caller (coded 0/2); mode2 omitted/invalid — run (F3, F4).
- Phased polyploid separators (`0|0|1`): SKIPPED — asserted from code
  reading (the `gsub` strips `|` identically to `/`); diploid phased calls
  were verified empirically via the caller.
- n.cores > 1 path: SKIPPED — single-core only; the argument is passed
  straight to `new("genlight", ...)`.

## Approval (Phase B)

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges, 2026-09-05 | Consequence acknowledged: half-missing calls become NA instead of 1 (genotype mode) or 0/1 (dosage mode) — handled correctly. |
| 2 | Approved | Arthur Georges, 2026-09-05 | Consequence acknowledged: ploidy values change for polyploid imports — derived exactly from GT arity. |
| 3 | Approved | Arthur Georges, 2026-09-05 | |
| 4 | Deferred | Arthur Georges, 2026-09-05 | LOW — not approved this round. |
| 5 | Deferred | Arthur Georges, 2026-09-05 | LOW — not approved this round. |
| 6 | Deferred | Arthur Georges, 2026-09-05 | LOW — not approved this round. |

All HIGH/MEDIUM findings approved 2026-09-05, with explicit
acknowledgement that objects users previously read will differ where
current behaviour is wrong (half-missing calls handled correctly, ploidy
stamped from the data). LOW findings (F4, F5, F6) were not approved this
round — deferred, not applied.

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-utils.vcfr2genlight.polyploid`
(base upstream/dev ddaed27):

- F1 (HIGH): any GT token still carrying "." after separator stripping
  is converted to NA before the mode rules — half-missing calls are
  missing data in both modes.
- F2 (MEDIUM): per-individual ploidy is computed from the allele count
  (arity) of the raw GT calls before separator stripping and stamped
  over adegenet's max-dosage inference; individuals with no called
  genotypes take the modal ploidy (2 if none at all).
- F3 (MEDIUM): the signature default is `mode2 = "genotype"` (matching
  the caller) instead of the accidental `mode2 = mode` closure; the
  `@details`/`@param` text is corrected accordingly.
- F4, F5, F6 (LOW): deferred — the empty stop() message, the dependency
  guards and the remaining documentation gaps are untouched.

Verification: baseline characterization test updated; all 21 assertions
pass; every diff from the pre-review snapshot is annotated
`[approved F1]`/`[approved F2]`/`[approved F3]`. Unchanged (LOW-pinned)
snapshots: multi-allelic drop warning, invalid-mode empty condition
message (F4), CHROM/POS/ID slot landing, dosage and genotype matrices
for complete calls.

Caller effect (`gl.read.vcf`, the only in-package caller — checked
across all 8 clones): completes against the fixed helper on mixed
diploid data (hand-verified: `./1` and `0/.` come through as NA, all
other calls unchanged) and on the triploid dosage fixture; on this
branch the caller still forces ploidy 2 — its companion fix is
PR #364. Note: a degenerate VCF in which every locus becomes
monomorphic after the half-missing NA conversion trips
`gl.compliance.check`'s pre-existing all-monomorphic crash
("Subsetting resulted in zero loci") inside the caller — a known
weakness of that function (same class as the gl.report.hwe monomorph
issue in NEWS), recorded for its own review, not a regression of this
change on realistic data. All clear.

```json
{
  "function": "utils.vcfr2genlight.polyploid",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ddaed27",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT1", "status": "applied", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DAT1", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "deferred", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "DEP1", "status": "deferred", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "deferred", "change": 6}
  ],
  "coverage_skipped": ["phased polyploid separators: asserted from code reading", "n.cores > 1: not exercised"],
  "status": "phase-c-complete",
  "pr": null
}
```
