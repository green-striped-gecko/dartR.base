# Review: gl.write.csv (dartR.base)
- Family mode: io
- Date: 2026-09-14
- Reviewer: Claude (Claude Opus 4.8), dartr-function-review v1.0.0
- Package commit: f02bd34 (upstream/dev; working copy verified identical by `git diff upstream/dev -- R/gl.write.csv.r`)
- Datasets: testset.gl (SNP), testset.gs (SilicoDArT) — dartR.data 1.2.5
- Baseline: tests/testthat/test-gl.write.csv.R (new file, snapshot captured pre-review; 15 assertions, all passing)

## Verdicts

**Standards: Needs work** — the FS backbone is present and the written file
is structurally correct, but the function returns a visible NULL (printing
"NULL" at the console on every un-assigned call, even at `verbose = 0`),
skips the family-standard `gl.check.wd()` on `outpath` (so a bad path throws
an opaque connection error), and the documentation has gaps.

**Spec: Ready** — the output matches the documentation exactly: header row
of loc.metrics column names then specimen ids; a filler/pop row (`*` under
the metadata, population under each specimen); one row per locus carrying
its metadata and the genotype scores. Confirmed on `testset.gl` (257 rows =
header + filler + 255 loci; 271 columns = 21 metadata + 250 individuals) and
`testset.gs` works too. No spec defect found.

## Blast radius

`gl.write.csv` has **no internal callers** anywhere in the dartRverse — it
is a user-facing io function. (The eBook, p48, currently shows base
`write.csv(gl, "genotypes.csv")` where it means `gl.write.csv` — an eBook
error recorded in the eBook audit, not a code issue here.)

## Independent verification (spec axis)

The write path was traced and read back (test 1): `cbind.data.frame(pop(x),
as.matrix(x))` then transpose puts specimens as columns and loci as rows;
`rbind(filler, loc.metrics)` aligns a `*` row with the population row; the
`cbind` places metadata to the left of the genotypes. The header carries the
loc.metrics names and the specimen ids; the first data row carries `*` and
the populations; each subsequent row carries a locus's metadata and its
0/1/2/NA scores (0/1/NA for SilicoDArT). Every element matched the
documentation.

## Findings

**F1 [MEDIUM, confidence: high] — returns a visible NULL (FS10, VRB5)**
`R/gl.write.csv.r:85` — `return(NULL)` is visible, so an un-assigned call
(the normal use of a side-effect writer) prints `NULL` at the console. It
prints at every verbosity, including `verbose = 0`, so the function is not
silent even when asked to be.
Failure scenario: `gl.write.csv(x, verbose = 0)` at the console prints a
bare `NULL` line; in a script sourced at `verbose = 0` the autoprint still
fires. Confirmed: `withVisible()` reports `visible = TRUE` (test 3).
Proposed change: `invisible(NULL)` — the house idiom for a side-effect
function (FS10).

**F2 [MEDIUM, confidence: high] — outpath is not validated via gl.check.wd (FS7)**
`R/gl.write.csv.r:53` — `outfilespec <- file.path(outpath, outfile)` is
built directly; the function does not call `gl.check.wd(outpath)`. Every
other io / output function in dartR.base resolves its output directory with
`outpath <- gl.check.wd(outpath, verbose = 0)` first.
Failure scenario: `gl.write.csv(x, outpath = "typo/dir")` fails inside
`write.table` with "cannot open the connection" — an opaque error naming a
file, not the missing directory — after the (potentially large) matrix
assembly has already run. Confirmed (test 4).
Proposed change: add `outpath <- gl.check.wd(outpath, verbose = 0)` after
the datatype check, matching the io family.
**Consequence: a non-existent outpath changes behaviour — it currently
errors; after the fix it falls back to tempdir() (with a warning at
verbose >= 1, per the reviewed gl.check.wd), as every sibling io function
already does.**

**F3 [INFO, confidence: high] — densifies the full genotype matrix (DAT6 (proposed rule))**
`R/gl.write.csv.r:58` — `as.matrix(x)` materializes the whole
individual-by-locus matrix. This is largely inherent to writing every
genotype to a flat CSV, but for a very large FBM-backed object it
materializes the entire matrix in memory (the example even wraps the input
in `gl.gen2fbm`). No action proposed; recorded so the cost is on record.

**F4 [LOW, confidence: high] — the coding note documents only SNP data (DOC5 (proposed rule))**
`R/gl.write.csv.r:11-13` — "0 = reference homozygous, 2 = alternate
homozygous, 1 = heterozygous, and NA = missing" describes SNP data only.
The function also accepts SilicoDArT (the second `@examples` line uses
`testset.gs`), where the values written are 0/1 presence-absence and NA.
Proposed change: note the SilicoDArT coding (0 absent, 1 present, NA
missing) alongside the SNP coding.

**F5 [LOW, confidence: high] — verbose param text predates DOC2 (DOC2)**
`R/gl.write.csv.r:21-23` — the `verbose` text is the pre-DOC2 form
("[default 2 or as specified using gl.set.verbosity]").
Proposed change: adopt the standard DOC2 verbose text.

## Proposed changes

1. Return `invisible(NULL)` instead of a visible NULL (F1).
2. Resolve `outpath` through `gl.check.wd(outpath, verbose = 0)`, matching
   the io family (F2).
   **Consequence: a non-existent outpath falls back to tempdir() with a
   warning at verbose >= 1 instead of throwing a connection error.**
3. Documentation: add the SilicoDArT coding note and adopt the standard
   verbose text (F4, F5).

## Coverage

- Standards walk: FS, DOC, VRB, DAT (DAT6 noted), DEP (n/a), PLT (n/a),
  STY — run
- Spec: output structure vs docs traced and read back on SNP and
  SilicoDArT; input untouched; verbose messaging — run
- FBM path (DAT6): densification identified (F3); not exercised on a real
  FBM fixture
- Caller survey (API3): no internal callers; user-facing io only
- GitHub issues / Google Group: SKIPPED — not queried this session

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur | |
| 2 | approved | Arthur | family-standard gl.check.wd fallback (consequence approved) |
| 3 | approved | Arthur | |

## Outcome

(pending Phase C)

```json
{
  "function": "gl.write.csv",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "f02bd34",
  "verdict_standards": "needs_work",
  "verdict_spec": "ready",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "FS10", "status": "proposed", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "FS7", "status": "proposed", "change": 2},
    {"id": "F3", "severity": "INFO", "confidence": "high", "rule": "DAT6", "status": "proposed", "change": null},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "proposed", "change": 3},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "DOC2", "status": "proposed", "change": 3}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture exercised", "GitHub issues not queried", "Google Group not queried"],
  "status": "awaiting-approval",
  "pr": null
}
```
