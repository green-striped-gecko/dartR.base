# Review: gl2gapit (dartR.base)

- Family mode: io
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203
- Datasets: testset.gl, testset.gs, platypus.gl
- Baseline: tests/testthat/test-gl2gapit.R (snapshot captured pre-review)

## Verdict

**Standards: Rework** — the entire roxygen header belongs to `gl2geno`
(`@name gl2geno`, LEA geno/lfmm description), so the exported function has no
man page of its own and its documentation describes a different converter.
**Spec: Rework** — despite `outfile`/`outpath` parameters and a progress
message naming output files, the function writes nothing to disk; the
returned hapmap data.frame is the only product, and the docs say the return
is `NULL`. What the function does and what any of its surfaces claim have to
be renegotiated, not patched.

The hapmap object itself is sound: header-row duplication matches the GAPIT
`header = FALSE` reading convention, genotypes double the allele letters, and
"00" is within GAPIT's accepted missing codes.

## Findings

**F1 [HIGH, confidence: high] — roxygen header is gl2geno's, so gl2gapit has no man page (DOC1, DOC4)**
`R/gl2gapit.r:1-32` — the block opens `@name gl2geno` with the LEA geno/lfmm
title and description. roxygen therefore never generates `man/gl2gapit.Rd`
(confirmed: absent from `man/`), while the function is exported in
NAMESPACE.
Failure scenario: `?gl2gapit` finds nothing; a user reading the header is
told the function writes LEA geno/lfmm files.
Proposed change: write a genuine header (`@name gl2gapit`, GAPIT hapmap
description, correct `@return`, working `@examples` — the current example
block also carries `gl2geno` remnants and mangled lines such as
`# assigning chromosomet1`) and run `devtools::document()`.

**F2 [HIGH, confidence: high] — no files are written, yet parameters and a message promise them (DOC5)**
`R/gl2gapit.r:34-36,113-117` — `outfile` and `outpath` are accepted and
resolved via `gl.check.wd()`, and at `verbose > 0` the function prints
`Output files: <outfile>.geno.lfmm.` (another gl2geno remnant). No write
call exists anywhere in the body; confirmed empirically that the output
directory is untouched.
Failure scenario: a user points `outpath` at a directory, sees the "Output
files" message, and finds nothing; also `@return "returns no value (i.e.
NULL)"` while the function actually returns the hapmap data.frame
invisibly.
Proposed change: decide the contract with the custodian — either write the
hapmap object to `file.path(outpath, outfile)` (and keep the message,
corrected), or drop `outfile`/`outpath` and the message and document the
invisible data.frame return. Either way, align `@return`.

**F3 [MEDIUM, confidence: high] — chromosome names are replaced by factor level codes (DOC5)**
`R/gl2gapit.r:93` — `x$chromosome <- as.factor(as.numeric(x$chromosome))`
converts the chromosome factor to its level indices. The mapping depends on
the object's level set, not on the names: platypus chromosome
`NC_041731.1_chromosome_4` became `5`, `..._chromosome_X1` became `23`
(confirmed). Subsetting or re-factoring the same data yields different
codes for the same chromosome.
Failure scenario: GAPIT Manhattan plots and per-chromosome output are
labelled with arbitrary, unstable integers; two exports of overlapping data
disagree on chromosome numbering. Nothing documents the recode.
Proposed change: document and stabilise the mapping (e.g. code chromosomes
as `match(chrom, sort(unique(chrom)))` with the name-to-code table reported
at `verbose >= 2`), or pass names through when GAPIT accepts them.

**F4 [MEDIUM, confidence: high] — slot warnings are not verbosity-gated (VRB3, VRB5)**
`R/gl2gapit.r:86,90` — the two `cat(warn(...))` calls for empty
chromosome/position slots print unconditionally; confirmed at
`verbose = 0`.
Failure scenario: a pipeline running silent still gets two warning lines per
call — and under the current dartRverse convention (`@position`/
`@chromosome` NULL until genome coordinates are assigned) the empty-slot
path is the normal case, so this fires on nearly every object.
Proposed change: gate with `if (verbose >= 2)` (or `>= 1` under VRB4 logic,
since fabricated coordinates affect the output).

**F5 [MEDIUM, confidence: high] — SilicoDArT admitted, then opaque crash (DAT7)**
`R/gl2gapit.r:52,57-59` — the datatype check keeps the permissive default;
SilicoDArT objects have `@loc.all` NULL, so `homs1`/`hets`/`homs2` are
zero-length and the fill loop dies with "number of items to replace is not a
multiple of replacement length" (confirmed).
Failure scenario: `gl2gapit(testset.gs)` fails with an error that names no
cause.
Proposed change: `accept = "SNP"` in `utils.check.datatype()`.

**F6 [LOW, confidence: high] — hardcoded `assembly = "Oilpalm"` (STY3)**
`R/gl2gapit.r:100` — every exported record carries the assembly string
"Oilpalm", a leftover from whatever example the code was adapted from.
Failure scenario: cosmetic but public: hapmap files from any species declare
an oil palm assembly.
Proposed change: set `assembly = NA` or expose it as a parameter.

**F7 [LOW, confidence: medium] — per-cell double loop densifies and crawls (STY2)**
`R/gl2gapit.r:61-75` — the nested `for` over individuals x loci performs
nInd x nLoc scalar assignments where three vectorised `ifelse`/`matrix`
operations on `x_mat` would do.
Failure scenario: large objects (100k loci) take minutes for a conversion
that should be seconds.
Proposed change: vectorise the genotype-to-letter mapping (index
`homs1`/`hets`/`homs2` by column with `col(x_mat)`).

## Proposed changes

1. Rewrite the roxygen header for GAPIT: name, title, description, `@return`
   documenting the invisible hapmap data.frame, runnable example;
   `devtools::document()` (F1).
2. Resolve the file-writing contract with the custodian; make code, message,
   parameters, and `@return` agree (F2). **Consequence: if writing is added,
   new files appear at outpath for existing callers.**
3. Stabilise and document (or drop) the chromosome recode (F3).
4. Gate the two slot warnings by verbosity (F4).
5. `accept = "SNP"` (F5).
6. Parameterise or NA the assembly field (F6).
7. Vectorise the genotype mapping loop (F7).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP (all deps in Imports), PLT (n/a),
  STY — run
- Spec: behaviour vs roxygen vs GAPIT hapmap conventions on testset.gl — run;
  file-output claim tested against the filesystem — run
- Genome-position lens: run — @position/@chromosome NULL fills dummies
  (chrom "1", pos 1..nLoc) with an ungated warning (F4); populated slots pass
  positions through but recode chromosome names to level indices (F3)
- SilicoDArT probe — run (F5)
- verbose = 0 silence probe — run (F4)
- "00" missing-code acceptance by GAPIT: checked against GAPIT
  numericalisation conventions statically, not by running GAPIT — GAPIT
  itself was not executed
- FBM path (DAT6): SKIPPED — no FBM fixture exercised for this function

## Approval (Phase B)

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges (2026-09-05) | Blanket approval, all severities |
| 2 | Approved | Arthur Georges (2026-09-05) | File-writing added: behaviour now matches docs |
| 3 | Approved | Arthur Georges (2026-09-05) | Blanket approval, all severities |
| 4 | Approved | Arthur Georges (2026-09-05) | Blanket approval, all severities |
| 5 | Approved | Arthur Georges (2026-09-05) | Blanket approval, all severities |
| 6 | Approved | Arthur Georges (2026-09-05) | Blanket approval, all severities |
| 7 | Approved | Arthur Georges (2026-09-05) | Blanket approval, all severities |

All findings approved 2026-09-05 via the formal approval boxes: blanket
class approval at every severity, explicitly including the consequence
that converted outputs change where they were wrong, and the DAT7 fatal
`accept = "SNP"` gate wherever SilicoDArT was silently admitted. For F2
the writing contract was resolved as: write the hapmap table to
`file.path(outpath, outfile)` so behaviour matches the documented
parameters and message.

## Outcome (Phase C)

All seven changes applied on branch `review-gl2gapit` (base
`upstream/dev`, ddaed27). New roxygen header generates
`man/gl2gapit.Rd` for the first time (F1); the hapmap table is now
written tab-delimited to `outpath/outfile` (F2); chromosome codes are
assigned by alphabetical order of the distinct names with the mapping
reported at `verbose >= 2` (F3); slot warnings gated (F4);
`accept = "SNP"` (F5); `assembly = NA` (F6); the per-cell double loop
replaced by a vectorised mapping (F7). Baseline characterization test
updated for five approved diffs (F2 file now written, F3 comment, F4
warnings gated, F5 clean datatype error, F6 assembly NA); 19/19
assertions pass — SNP-path genotype snapshots (doubled letters, "00"
missing, duplicated header row, dims) unchanged. End-to-end run on
`testset.gl` at `verbose = 3` returned a 756 x 285 data.frame and wrote
756 lines to tempdir. Sibling-caller grep across the 8 dartRverse
clones: no callers of `gl2gapit` — all clear. NEWS entry added.
PR: green-striped-gecko/dartR.base#337.

```json
{
  "function": "gl2gapit",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203",
  "verdict_standards": "rework",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "applied", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DAT7", "status": "applied", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "STY3", "status": "applied", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "medium", "rule": "STY2", "status": "applied", "change": 7}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "GAPIT not executed; missing-code acceptance checked statically"],
  "status": "phase-c",
  "pr": 337
}
```
