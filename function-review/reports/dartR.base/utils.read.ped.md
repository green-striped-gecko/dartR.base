# Review: utils.read.ped (dartR.base)

- Family mode: io (chained with `gl.read.PLINK` in the manifest, but the
  chain is nominal: `gl.read.PLINK` reads via `snpStats::read.plink`; the
  only in-package caller of `utils.read.ped` is `gl.report.ld.map`, which
  passes `show_warnings = F` and `na.strings = NA`)
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ddaed27 (upstream/dev; file verified identical on local HEAD ed99203)
- Datasets: synthetic `.ped`/`.map` filesets (see Fixtures)
- Baseline: tests/testthat/test-utils.read.ped.R (snapshot captured
  pre-review, bugs included)

Provenance note: the function is a vendored, lightly modified copy of
`snpStats::read.pedfile`. Three of the defects below (F1, F3, F4) are
present verbatim in the installed snpStats original — they are inherited,
not introduced — but they ship in dartR.base and act on dartR data, so they
are findings here regardless; F2 is introduced by the dartR modification.

## Verdict

**Standards: Needs work** — the roxygen header is a placeholder skeleton
("Custodian to provide" throughout) with a wrong `@return`, and the file
carries an unrelated dead function.
**Spec: Needs work** — the default parse path is correct (hand-verified
against the fixture, including `na.strings` missing-data handling), but
`lex.order = TRUE` silently inverts dosage semantics, `show_warnings =
FALSE` disables a data-cleaning step, and two argument/detection paths
corrupt silently.

What works well: the core two-pass parse reproduces snpStats semantics
faithfully — fam columns, first-seen/second-seen allele assignment, and
`na.strings = "0"` missing handling all verified cell-by-cell.

## Findings

**F1 [HIGH, confidence: high] — lex.order=TRUE swaps map alleles but not genotypes (DAT2)**
`R/utils.read.ped.r:176-183` — `snpStats::switch.alleles(result, swa)` is
called for value, and the value is discarded; only the `a1`/`a2` bookkeeping
is swapped. Verified: on the fixture, `lex.order = TRUE` returns a map with
snp3/snp4 alleles swapped to lexicographic order while the genotype raw
bytes are `identical()` to the `lex.order = FALSE` run — the dosages at
those loci now count the wrong allele relative to the returned map.
(Inherited: the installed `snpStats::read.pedfile` has the same unassigned
call.)
Failure scenario: any caller using `lex.order = TRUE` to normalise allele
order gets silently complemented genotypes at every locus whose first-seen
allele was lexicographically later.
Proposed change: `result <- snpStats::switch.alleles(result, swa)`.

**F2 [HIGH, confidence: high] — show_warnings=FALSE also disables the multi-allelic NA reset (VRB3, DAT1)**
`R/utils.read.ped.r:155-158` — the dartR modification wraps the original
unconditional cleanup `result[, mallelic] <- r0` inside
`if (any(mallelic & show_warnings==T))`. Verified: on a fixture with a
flagged multi-allelic locus, `show_warnings = TRUE` returns the whole locus
as NA (snpStats-documented behaviour); `show_warnings = FALSE` returns
retained genotypes for the same locus (0, 1, NA). A flag named for message
suppression changes the data. The in-package caller `gl.report.ld.map`
passes `show_warnings = F`, so it is on the affected path (its inputs come
from `gl2plink` on a genlight and are biallelic by construction, which is
why this has not surfaced).
Failure scenario: quiet callers keep genotypes at loci the function itself
has classified as unreliable, with no message and no NA reset.
Proposed change: perform the reset unconditionally; gate only the
`warning()` calls on `show_warnings`.

**F3 [HIGH, confidence: high] — multi-allelic detection is masked when a novel allele pairs with a known one (DAT1)**
`R/utils.read.ped.r:125-140` — `mallelic <- mallelic | !(akm | one | two)`
uses `one`/`two`, which are set per individual across both allele columns,
so a genotype like `A T` (A = a1, T = a third allele) sets `one` on column
1 and thereby exempts column 2's novel T from detection. Verified: an
`A T` carrier is coded raw 1 = homozygous A/A, no warning, locus not
flagged; order matters — `T A` at the same locus would be flagged.
(Inherited from snpStats.)
Failure scenario: real tri-allelic loci pass undetected whenever novel
alleles happen to co-occur with a1/a2, with het-with-novel calls silently
promoted to homozygotes.
Proposed change: track per-column matches (`br2`/`br5` results) instead of
the accumulated `one`/`two` when computing `mallelic`.

**F4 [MEDIUM, confidence: high] — the split argument is ignored on data lines (API2 principle; inherited)**
`R/utils.read.ped.r:114` — the locus-count probe (line 93) honours `split`,
but the per-line loop hardcodes `strsplit(line, "\t| +")`. Verified: a
comma-separated ped read with `split = ","` returns without error —
individual ids equal to the entire raw line, sex/affected NA, and all
genotypes NA.
Failure scenario: a caller with a non-default separator receives an all-NA
SnpMatrix and mangled ids instead of an error.
Proposed change: `strsplit(line, split)` in the loop.

**F5 [MEDIUM, confidence: high] — placeholder documentation with a wrong @return (DOC1, DOC5)**
`R/utils.read.ped.r:1-26` — `@title` is "An internal script [Custodian to
provide a title]", every `@param` reads "Custodian to provide", and
`@return` claims "The resultant genlight object" while the function returns
a `list(genotypes = SnpMatrix, fam, map)`. `@export` is commented out, so
the roxygen block documents nothing that is exposed (the style guide
expects `utils.*` exported with the internal-use warning; the warning line
is present).
Failure scenario: a developer trusting `@return` treats the result as a
genlight.
Proposed change: write the real title/params/details; correct `@return`;
decide export status per the utils convention; re-document (DOC4).

**F6 [LOW, confidence: high] — unrelated dead function in the same file (FS1)**
`R/utils.read.ped.r:197-218` — `rotate.matrix()`, an image-rotation helper,
shares the file; it is unexported, undocumented, and has no callers
anywhere in `R/` (grep). Proposed change: delete it (or move it beside its
user if one exists elsewhere in the verse).

**F7 [LOW, confidence: high] — DEP1 guard returns -1 instead of stopping (DEP1)**
`R/utils.read.ped.r:39-47` — missing snpStats triggers `cat(error(...))`
and `return(-1)`, so callers receive `-1` where a list is expected and fail
later. snpStats is in Imports, so the guard is unreachable in a working
install; if kept, it should `stop(error(...))`.

## Proposed changes

1. Assign the `switch.alleles` result (F1). **Consequence: numerical output
   changes for `lex.order = TRUE` callers — genotypes at swapped loci are
   complemented to match the returned map.**
2. Decouple the multi-allelic NA reset from `show_warnings` (F2).
   **Consequence: numerical output changes for `show_warnings = FALSE`
   callers with flagged multi-allelic loci (including `gl.report.ld.map`
   inputs, if any ever contain them) — those loci become NA.**
3. Fix the per-column multi-allelic detection (F3). **Consequence: loci and
   genotypes previously mis-coded as homozygous become flagged/NA.**
4. Honour `split` on data lines (F4).
5. Real documentation and correct `@return` (F5); delete `rotate.matrix`
   (F6); make the DEP1 guard stop (F7).
6. Offer F1/F3/F4 upstream to snpStats (they reproduce in
   `snpStats::read.pedfile` as installed) — team decision on whether to
   carry local fixes meanwhile.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, STY — run (PLT n/a: no plots;
  FS2/FS3/FS9 n/a by design for this vendored helper — it has no
  verbosity/flag preamble at all, accepted as-is for a vendored copy but
  noted for the promote decision).
- Spec: default parse, na.strings, fam/map contents — run and
  hand-verified; lex.order — run (F1); show_warnings T/F on flagged loci —
  run (F2); novel-allele masking — run (F3); split override — run (F4);
  monomorphic and no-data warnings — observed incidentally.
- gzip input path (`gzfile`): SKIPPED — plain-text connections only were
  exercised; `gzfile()` transparently reads both, asserted from code
  reading.
- Caller interaction (`gl.report.ld.map` passing `na.strings = NA`, which
  stops "0" being treated as missing): noted, out of scope here — flag for
  that function's own review.

## Fixtures (all synthetic, written by the probes/tests)

- `toy.ped`/`toy.map`: as for the gl.read.PLINK review (4x4, one `0 0`
  missing genotype).
- `multi.ped` (`A T` het-with-novel carrier), `multi3.ped` (`T T` novel
  homozygote), `comma.ped` (comma-separated), `miss.ped` (`0 0` pair).

## Approval (Phase B)

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges, 2026-09-05 | Consequence acknowledged: `lex.order = TRUE` genotypes at swapped loci are complemented to match the returned map. |
| 2 | Approved | Arthur Georges, 2026-09-05 | Consequence acknowledged: `show_warnings = FALSE` callers with flagged multi-allelic loci now get those loci as NA — as documented. |
| 3 | Approved | Arthur Georges, 2026-09-05 | Consequence acknowledged: loci previously mis-coded as homozygous become flagged/NA. |
| 4 | Approved | Arthur Georges, 2026-09-05 | |
| 5 | Approved (F5 only; F6, F7 deferred) | Arthur Georges, 2026-09-05 | Fixes applied in the vendored copy, not snpStats upstream. |
| 6 | Noted — team decision pending | Arthur Georges, 2026-09-05 | Offering F1/F3/F4 upstream to snpStats remains open; local fixes carried meanwhile. |

All HIGH/MEDIUM findings approved 2026-09-05, with explicit
acknowledgement that objects users previously read will differ where
current behaviour is wrong (`lex.order` and `show_warnings` now behave
as documented). LOW findings (F6 `rotate.matrix`, F7 DEP1 guard) were
not approved this round — deferred, not applied.

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-utils.read.ped` (base upstream/dev
ddaed27), in the vendored copy as scoped (snpStats upstream untouched):

- F1 (HIGH): `result <- snpStats::switch.alleles(result, swa)` — the
  switched genotype matrix is now assigned.
- F2 (HIGH): the multi-allelic NA reset runs unconditionally; only the
  `warning()` is gated on `show_warnings`.
- F3 (HIGH): multi-allelic detection now tests each allele column's
  match (`br2`/`br5`) instead of the accumulated `one`/`two`, so a novel
  allele paired with a known one no longer escapes detection.
- F4 (MEDIUM): the per-line loop now splits on `split` rather than the
  hardcoded `"\t| +"`.
- F5 (MEDIUM): real roxygen documentation (title, params, details,
  corrected `@return` — the function returns
  `list(genotypes, fam, map)`, not a genlight); exported with
  `@keywords internal` and the internal-use warning per the utils
  convention (NAMESPACE gains `export(utils.read.ped)`); re-documented.
- F6, F7 (LOW): deferred — `rotate.matrix` and the `return(-1)` DEP1
  guard are untouched.

Verification: baseline characterization test updated; all 28 assertions
pass; every diff from the pre-review snapshot is annotated
`[approved F1]`/`[approved F2]`/`[approved F3]`/`[approved F4]`.
Unchanged snapshots: default parse (fam/map/genotype coding,
`na.strings = "0"` missing handling) — verified identical.

Caller effect (`gl.report.ld.map`, the only in-package caller; passes
`show_warnings = F`, `na.strings = NA`): with the F2 fix the
multi-allelic NA cleanup now runs on its path regardless of the flag.
Its inputs are written by `gl2plink` from a genlight and are biallelic
by construction, so the cleanup cannot trigger on them — no behaviour
change. Verified empirically: `gl.report.ld.map` run end-to-end on a
mapped 30-ind x 60-loci `testset.gl` subset against the fixed helper —
completes and returns the LD table. It uses `lex.order` default FALSE
(F1 path not taken) and the default `split`, which matches `gl2plink`'s
space-separated output (F4 no-op). dartR.popgen's `gl.ld.haplotype`
calls its own vendored copy in `dartR.popgen/R/utils.ld.r`, not this
one — unaffected. All clear.

```json
{
  "function": "utils.read.ped",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ddaed27",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT2", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "VRB3", "status": "applied", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "DAT1", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "API2", "status": "applied", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "FS1", "status": "deferred", "change": 5},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DEP1", "status": "deferred", "change": 5}
  ],
  "coverage_skipped": ["gzip input path: asserted from code reading", "caller-side na.strings interaction: deferred to gl.report.ld.map review"],
  "status": "phase-c-complete",
  "pr": 361
}
```
