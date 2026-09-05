# Review: gl.read.fasta (dartR.base)

- Family mode: io (wrapper over utils.read.fasta + merge_gl_fasta —
  reviewed alongside; genotype-calling findings live in the
  utils.read.fasta report)
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203 (R/gl.read.fasta.r byte-identical to upstream/dev
  ddaed27)
- Datasets: locus_1.fas (dartR.data extdata, the documented example);
  synthetic FASTA fixtures (inventory below)
- Baseline: tests/testthat/test-gl.read.fasta.R (28 assertions, snapshot
  captured pre-review, defective behaviour pinned with [BUG] tags)

## Verdict

**Standards: Needs work** — preamble, compliance pass and history conform;
the closing fbm block violates the phase order and the verbosity contract.
**Spec: Rework** — at the default verbosity the function returns an object
with no genotype data at all (nInd = 0, `@gen` and `@fbm` both empty), for
every input including the packaged documentation example; the `fbm`
argument has no effect in any direction. What the user receives is decided
by `verbose`, not by the data or the arguments.

What works well: everything up to the final block. The multi-file merge is
correct — individuals absent from a file get NA at that file's loci,
ploidy 2 throughout, locus metrics registered 1:1 — verified on the
`verbose = 3` path where the data survive.

## Synthetic fixtures

Shared with the utils.read.fasta review (full inventory and per-column
design in that report): `clean.fas` (6 x 10, four verifiable SNP columns),
`clean2.fas` (5 x 6, overlapping individuals ind1-4 plus ind7 — merge/NA
fixture), `mono.fas` (no polymorphism), `wrapped.fas` (multi-line
records), `dupnames.fas` (duplicate record names).

## Findings

**F1 [BLOCKER, confidence: high] — the fbm block destroys or mistypes the returned object; the fbm argument is dead (FS6, DAT2; principle: the returned object must carry the data read — no catalogue rule covers return-payload integrity, worth reporting to the skill maintainer)**
`R/gl.read.fasta.r:114-119` — the block reads:

```r
if (fbm) {}
x <- gl.gen2fbm(x, verbose = verbose)
if (verbose>2) {
  cat(report(" Created an  file-backed matrix (fbm) dartR object\n"))
} else x@fbm <- NULL
```

The `if (fbm)` guard is an empty statement, so `gl.gen2fbm()` always runs
and clears `@gen` (its documented XOR). The `else` then belongs to the
*verbosity* test: at `verbose <= 2` it wipes `@fbm` as well. All four
combinations verified empirically:

- `verbose = 2` (default), `fbm = FALSE`: `@gen` empty, `@fbm` NULL,
  `nInd(x) = 0`, `as.matrix()` errors — all genotype data gone.
- `verbose = 2`, `fbm = TRUE`: identical — data gone despite the request.
- `verbose = 3`, `fbm = FALSE`: FBM-backed object returned although fbm
  was declined (genotypes intact on this path, hand-verified).
- The packaged example (`locus_1.fas`, default verbosity) returns a
  0-individual shell.

Failure scenario: every default-settings user since this block landed gets
an empty object; scripts that only check `nLoc` or locus names proceed on
a data-less shell.
Proposed change: rewrite as
`if (fbm) x <- gl.gen2fbm(x, verbose = verbose)` with the report message
gated inside, and move the block above FLAG SCRIPT END (F5).
**Consequence: numerical output changes for every caller at
verbose <= 2 — from an empty object to the actual data.**

**F2 [HIGH, confidence: high] — the roxygen missing-data and multiallelic claims are false of the delegate (DOC5, proposed rule)**
`R/gl.read.fasta.r:16-30` — this header documents the read semantics, but
the engine does otherwise (utils.read.fasta report, F1/F2): V/H/D/B/N and
gaps are "taken as missing data" only in the allele pool — the carrying
individuals are coded heterozygous, never NA; "SNPs with more than two
alleles are skipped" fails for the commonest triallelic configuration; and
"the allele with the highest frequency is taken as the reference" inverts
when the modal genotype is a het. All verified through this function's
public surface at `verbose = 3`.
Failure scenario: a user relies on the documented semantics; heterozygosity
is silently inflated and polarity silently flipped.
Proposed change: fix the engine (utils.read.fasta changes 1-5); this
header is then accurate as written.

**F3 [MEDIUM, confidence: high] — an all-monomorphic input fails with an opaque error (FS5)**
`R/gl.read.fasta.r:86-92` with `R/utils.read.fasta.r:233-243` — a file
without polymorphism makes utils.read.fasta return NULL; merge_gl_fasta
filters the NULLs and, left with an empty list, takes the multi-file
branch and dies in `Reduce`/`gl.compliance.check` with "argument of length
0" (verified).
Failure scenario: a legitimate but invariant alignment produces an error
that names no cause; the only explanation printed is the ungated warn
line, which a `verbose = 0` logger would attribute to nothing.
Proposed change: after the NULL filter, stop with a clear message when no
file yielded SNPs.

**F4 [MEDIUM, confidence: high] — wrapped (multi-line) FASTA is unsupported and fails opaquely (FS5)**
`R/utils.read.fasta.r:49-53` (surfaces here) — records are split assuming
exactly two lines per individual; standard line-wrapped FASTA mis-groups
and dies with "missing value where TRUE/FALSE needed" plus split warnings
(verified). Neither header documents the two-line requirement.
Failure scenario: FASTA from most aligners/exporters wraps at 60-80
columns; users get an error naming no cause — or, for some record counts,
a mis-parse instead of an error.
Proposed change: either concatenate lines between `>` headers on read
(correct fix) or stop with "sequences must be on a single line".

**F5 [MEDIUM, confidence: high] — verbose = 0 is not silent (VRB5)**
`R/gl.read.fasta.r` (whole surface) — the delegate's ungated
multiallelic-skip and no-polymorphism messages print at `verbose = 0`
(verified: two lines captured). Same root as utils.read.fasta F8; listed
here because VRB5 is a property of the exported surface.
Proposed change: covered by utils.read.fasta change 8.

**F6 [LOW, confidence: high] — FLAG SCRIPT END precedes the final phase, and the reader returns invisibly (FS9, FS10)**
`R/gl.read.fasta.r:108-122` — "Completed: gl.read.fasta" prints before the
fbm conversion runs, so gl.gen2fbm's own start/end banners appear after
the wrapper's completion flag (visible in every verbose run); and the
result is returned `invisible(x)` although a reader's value is always
assigned — `gl.read.csv` returns visibly. Folded into F1's rewrite for the
ordering; the return style is a one-word change.
Proposed change: move FLAG SCRIPT END last; `return(x)`.

**F7 [LOW, confidence: high] — roxygen details (DOC5 proposed rule, DOC1, DOC7 proposed rule)**
`R/gl.read.fasta.r:39-46,54-55` — the `fbm` param says the matrix lands
"in its @genome slot" (the slot is `@fbm`); the `verbose` default clause
deviates from DOC2; `@author` has no `Author(s):` line (DOC7); tag order
places `@details` after `@examples`' neighbours per the old convention.
Proposed change: align header, re-document (DOC4).

### @position / @chromosome (PR #330 interaction)

Verified empirically on the surviving (`verbose = 3`) path, single- and
multi-file: `@position` and `@chromosome` are NULL on return. Alignment
coordinates exist only inside locus names (`<file>_<column>`); no
`SnpPosition` loc.metrics column is created, so neither the ddaed27
compliance fill nor the #330 clearing logic engages. Conforms to the
genome-only convention; the docs make no coordinate claim. Nothing to
re-propose.

## Proposed changes

1. Rewrite the fbm block: convert only when `fbm = TRUE`, gate the report
   message, never null `@fbm` on a converted object; move it before FLAG
   SCRIPT END and return visibly (F1, F6).
   **Consequence: numerical output changes for every caller at
   verbose <= 2 — from an empty object to the actual data; verbose >= 3
   callers now get a plain (non-FBM) object unless they ask for fbm.**
2. Stop with a clear message when no input file yields SNPs (F3).
3. Support or explicitly reject wrapped FASTA (F4).
4. Engine fixes for missing data, multiallelic skip and ref polarity are
   utils.read.fasta changes 1-5 (F2, F5 here); adopt that report's
   consequence lines.
5. Roxygen alignment and re-document (F7).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run (PLT n/a; FS4
  n/a, file-path input; DEP: bigsnpr/bigstatsr in Imports so the
  unconditional gl.gen2fbm call at least resolves; FS8 verified — history
  carries the gl.read.fasta call, plus gl.recalc.metrics' own entry).
- Spec: single file, default and verbose = 3 — run; fbm = TRUE/FALSE at
  both verbosity bands — run (all four combinations); two-file merge with
  disjoint and overlapping individuals, NA fill hand-verified — run;
  packaged example locus_1.fas — run (empty object confirmed); duplicate
  names through compliance (renamed ind1_1) — run; all-monomorphic input —
  run; wrapped FASTA — run; verbose = 0 silence — run (fails, F5).
- DAT1/DAT2 on the surviving path: verified (ploidy 2, loc.metrics rows =
  nLoc, merge NA semantics correct).
- FBM genotype fidelity: verified at verbose = 3 (decoded matrix matches
  the hand-computed genotypes).
- `parallel = TRUE` pass-through: SKIPPED here — failure modes verified
  and reported at the utils.read.fasta level (F7 there).

## Approval (Phase B)

All BLOCKER, HIGH and MEDIUM findings approved by Arthur Georges
2026-09-05 via the formal approval boxes, explicitly acknowledging the
stated consequences — in particular that gl.read.fasta finally returns
genotypes at default verbosity, so objects users previously read will
differ where the current behaviour is wrong. LOW findings were not
approved this round and are deferred.

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur Georges, 2026-09-05 | F1 (BLOCKER); consequence acknowledged: output changes for every caller at verbose <= 2 (empty shell -> actual data), verbose >= 3 callers get a plain object unless fbm is requested. The F6 (LOW) portions of this change (moving FLAG SCRIPT END, visible return) are NOT applied — F6 deferred |
| 2 | approved | Arthur Georges, 2026-09-05 | F3 (MEDIUM) |
| 3 | approved | Arthur Georges, 2026-09-05 | F4 (MEDIUM); rejection option implemented (clear error + docs), keeping the parse fix out of the engine file to avoid cross-branch conflicts |
| 4 | approved | Arthur Georges, 2026-09-05 | F2 (HIGH) + F5 (MEDIUM); the fix is the utils.read.fasta engine changes, applied on that function's own branch — no hunks here |
| 5 | deferred | Arthur Georges, 2026-09-05 | F7 (LOW) not approved this round |

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-gl.read.fasta` (base upstream/dev
ddaed27): change 1 (F1 BLOCKER, fbm-guard portion only), change 2 (F3
MEDIUM), change 3 (F4 MEDIUM, reject option). Change 4 (F2 HIGH, F5
MEDIUM) is approved but carries no hunks on this branch — the engine fix
lands on `review-utils.read.fasta`; the branches are independent, each
based on upstream/dev. Change 5 (F7 LOW) and the F6 (LOW) portions of
change 1 deferred.

Verification (against the unfixed engine, as the baseline was
snapshotted):

- Baseline characterization suite: 30 pass, 0 fail, 1 skip. The skip is
  the verbose = 0 silence test ([approved F5]), self-gated at run time on
  whether the loaded engine still prints ungated — it activates
  automatically once the utils.read.fasta branch merges.
  Engine-dependent expectations (locus sets, locus counts) are derived
  from the loaded engine at run time so the file holds on this branch
  alone and after the engine branch merges. All [BUG] pins for F1, F3,
  F4 flipped with [approved F<n>] tags; no unexplained diffs.
- E2e at verbose = 3, single file: nInd 6, 6 @gen slots, @fbm NULL,
  genotype matrix intact (clean_5 = 0,2,1,0,2,0 hand-verified).
- Packaged documentation example at default verbosity: nInd 117,
  nLoc 9, genotypes present, as.matrix 117 x 9 — previously a
  0-individual shell.
- fbm = TRUE at default verbosity: @gen empty, @fbm populated, nInd 6 —
  previously the data were lost.
- Wrapped FASTA: clear "must be on a single line" fatal error;
  all-monomorphic input: clear "no SNPs found" fatal error.

Integration probe (combined behaviour: this branch's wrapper + the
`review-utils.read.fasta` engine, PR #362, loaded together in a scratch
state, 2026-09-05):

- test-gl.read.fasta.R: 31 pass, 0 fail, 0 skip -- the verbose-0
  silence test activates against the gated engine and passes.
- test-utils.read.fasta.R: 32 pass, 0 fail.
- End to end at default verbosity, clean.fas + clean2.fas: nInd 7,
  nLoc 5, genotypes present in @gen (@fbm NULL), N and gap cells read
  as NA, triallelic columns absent, absent individuals NA at the other
  file's loci -- the clean/clean2 fixtures read correctly end to end.
- verbose = 0 through the full public surface: zero captured lines.

```json
{
  "function": "gl.read.fasta",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "BLOCKER", "confidence": "high", "rule": "FS6", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "applied-in-utils.read.fasta", "change": 4},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "applied", "change": 2},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "applied", "change": 3},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "applied-in-utils.read.fasta", "change": 4},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "FS9", "status": "deferred", "change": 1},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "deferred", "change": 5}
  ],
  "coverage_skipped": ["parallel pass-through: verified at utils.read.fasta level"],
  "status": "phase-c-applied",
  "pr": 359
}
```
