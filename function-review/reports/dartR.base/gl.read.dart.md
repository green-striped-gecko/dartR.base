# Review: gl.read.dart (dartR.base)

- Family mode: io (entry point of the DArT read chain; delegates to `utils.read.dart` and `utils.dart2genlight`)
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ddaed27 (upstream/dev); `R/gl.read.dart.r` identical on the empirical checkout ed99203
- Datasets: `testset_SNPs_2Row.csv` + `testset_metadata.csv` (dartR.data); synthetic fixtures (see `utils.read.dart.md`)
- Baseline: `tests/testthat/test-gl.read.dart.R` (new; snapshot captured pre-review, 39 expectations passing)

## Chain state (context, not findings)

On the reviewed chain the returned object is class `dartR`, ploidy 2 on all
individuals, `loc.metrics` rows register 1:1 with the 255 loci, `ind.metrics`
rows 1:1 with the 250 individuals, and the chain ends in
`gl.compliance.check` (DAT5 satisfied). With PR #330 applied (as on the
checkout), `@position` and `@chromosome` are NULL and the tag offset lives in
`loc.metrics$SnpPosition` — recorded here as context; on pre-#330 upstream/dev
the compliance check still copies SnpPosition into `@position`. Genotype
encoding was verified against raw csv cells (see `utils.dart2genlight.md`).

## Verdicts

**Standards: Needs work** — FS2/FS3 present but preceded by executable code;
verbose = 0 leaks a full downstream log; one dropped return value.
**Spec: Needs work** — the documented happy path works and was verified
cell-by-cell; robustness claims break on files missing optional DArT columns,
and two inherited chain defects (service/plate garbage, metafile filtering)
surface in this function's output.

## Findings

**F1 [MEDIUM, confidence: high] — verbose = 0 leaks the full gl.compliance.check log (VRB5)**
`R/gl.read.dart.r:282` — `gl.compliance.check(glout)` is called without
`verbose`, so it resolves its own default (2) regardless of what the user
passed.
Failure scenario (verified): `gl.read.dart(..., verbose = 0)` prints 18 lines
of compliance-check progress.
Proposed change: `gl.compliance.check(glout, verbose = verbose)`.

**F2 [MEDIUM, confidence: high] — utils.recalc.maf result is discarded (DAT4, DOC5)**
`R/gl.read.dart.r:215` — `utils.recalc.maf(glout, verbose = 0)` without
assignment; `glout` is unchanged, yet the next lines print "MAF calculated
and added to the locus metrics" and `loc.metrics.flags$maf` was already set
TRUE at line 163.
Failure scenario: between this point and the compliance check the object
carries a TRUE maf flag with no maf column; the final object only has maf
because `gl.compliance.check` independently recalculates metrics (verified:
maf present and non-NA in the final object). Any refactor of the compliance
check, or any future caller of the pre-compliance object, inherits a silent
lie.
Proposed change: `glout <- utils.recalc.maf(glout, verbose = 0)`.

**F3 [MEDIUM, confidence: high] — lastmetric autodetection runs before the preamble and fails opaquely (FS2, FS5)**
`R/gl.read.dart.r:87-98` — `getLastMarkerMetaDataField()` executes before SET
VERBOSITY/FLAG SCRIPT START, reads the file with no existence check, and
indexes `last(which(top[,1]=="*"))+1` with no guard.
Failure scenario (verified): a valid headerless csv read with
`topskip = 0` dies with "argument of length 0" before any dartR banner; a
mistyped filename errors from inside the helper. The helper also duplicates
the `*`-row detection already in `utils.read.dart` and reads the file an
extra time.
Proposed change: move the call after the preamble, guard the no-`*`-rows case
with the same clear error `utils.read.dart` uses, and skip the helper
entirely when the user supplied `lastmetric`.

**F4 [MEDIUM, confidence: high] — rdepth calculation crashes when AvgCount columns are absent (FS5)**
`R/gl.read.dart.r:189-211` — the loop indexes
`loc.metrics$CallRate/AvgCountRef/AvgCountSnp` without presence checks; a
missing column yields NULL and the assignment dies with "replacement has
length zero".
Failure scenario (verified, `synth_1row_nocounts.csv`): a report without
`AvgCountRef`/`AvgCountSnp` — the roxygen explicitly says metrics are
computed "because sometimes they are not reported" — fails with the opaque
replacement error.
Proposed change: guard the block on all three columns existing; otherwise set
rdepth NA with a gated warning.

**F5 [LOW, confidence: high] — roxygen defects (DOC2, DOC5, DOC7 proposed, DOC3)**
`R/gl.read.dart.r:13, 20-21, 29-31, 53, 55-59` — `covfilename` doc has typo
"sse" and no standard default bracket; `lastmetric` is described as
"Deprecated" but is an active override (and when NULL triggers F3's helper,
so `utils.read.dart`'s documented "RepAvg" default is unreachable from this
wrapper); the `verbose` text deviates from the DOC2 canon ("default 2, or as
set by gl.set.verbose()" — the function is `gl.set.verbosity`); `@author` has
a Custodian but no Author(s) line (DOC7, proposed); the `@examples` run two
full imports including `fbm = TRUE` (file-backed conversion) unwrapped by
`\donttest{}` (DOC3).
Proposed change: correct the six fragments; `devtools::document()` (DOC4).

**F6 [INFO, confidence: high] — ordering and hygiene**
`R/gl.read.dart.r:284-295` — FLAG SCRIPT END prints before the fbm
conversion, so "Completed:" precedes real work (FS9);
`:166` `glout@other$verbose <- 2` hardcodes the object's verbosity default
regardless of the call; `:244` `flags$monomorphs <- FALSE` is set before the
optional `mono.rm` filtering rather than reflecting its outcome. No change
proposed; touch alongside approved changes only.

## Proposed changes

1. Pass `verbose` through to `gl.compliance.check` (F1).
2. Assign the `utils.recalc.maf` result (F2).
3. Relocate and guard the lastmetric helper; honour user-supplied
   `lastmetric` without reading the file twice (F3).
4. Guard the rdepth block against missing CallRate/AvgCount columns (F4).
5. Roxygen corrections + `devtools::document()` (F5).

Inherited-defect note for the apply phase: fixes to service/plate garbage and
metafile filtering live in `utils.read.dart` (change 3 there) and
`utils.dart2genlight` (change 1 there); this function needs no code change
for them but its baseline test snapshots their current surface behaviour.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run (PLT: no plots,
  n/a)
- Spec: behaviour vs roxygen on real 2-row fixture + synthetic 1-row,
  headerless, and missing-column variants — run
- 1-row vs 2-row autodetection, topskip autodetect/explicit, service/plate
  rows, NA code "-", duplicate ids, recalc/mono.rm flags — run
- `fbm = TRUE` path (DAT6): SKIPPED — bigstatsr backing-file setup is
  environment-dependent; reviewed by reading only (`glout@fbm <- NULL` on the
  non-fbm path is fine on the dartR class)
- `covfilename` deprecated path: run by reading only (trivial alias)
- recalc=FALSE / mono.rm=TRUE combinations: SKIPPED — behaviour dominated by
  gl.filter.monomorphs/gl.filter.allna, reviewed separately in the campaign

## Approval (Phase B)

Approved 2026-09-05 by Arthur Georges via the formal approval boxes, covering
all BLOCKER/HIGH/MEDIUM findings for the read-chain review round, with
explicit acknowledgment of the consequence that objects users previously read
will differ where current behaviour is wrong. LOW findings were not approved
this round and are deferred; INFO items carry no action.

| Change | Decision | By | Note |
|---|---|---|---|
| 1 (F1, MEDIUM) | approved | Arthur Georges | verbose = 0 now fully silent |
| 2 (F2, MEDIUM) | approved | Arthur Georges | |
| 3 (F3, MEDIUM) | approved | Arthur Georges | clear errors replace opaque pre-banner failures |
| 4 (F4, MEDIUM) | approved | Arthur Georges | reports without AvgCount columns now load, rdepth NA |
| 5 (F5, LOW) | deferred | Arthur Georges | LOW findings not approved this round |
| — (F6, INFO) | no action | | |

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-gl.read.dart` (base upstream/dev
ddaed27): changes 1-4 (F1-F4). Change 5 (F5, LOW) deferred - roxygen
untouched, so `man/gl.read.dart.Rd` and NAMESPACE are unchanged
(`devtools::document()` run to confirm).

Verification: baseline characterization test run against the branch -
40 expectations, 0 failures. Diffs from the Phase A baseline all map to
approved findings and are marked `# [approved F<n>]` in the test file:
verbose = 0 now prints nothing (F1), the no-`*`-rows and mistyped-filename
cases now fail with clear errors (F3), and a report missing
AvgCountRef/AvgCountSnp now loads with rdepth NA (F4). The @position
expectations were made #330-state-agnostic (pre-#330 base: @position carries
SnpPosition; post-#330: NULL) - a base-commit difference, not a behaviour
change of this branch. End-to-end run on `testset_SNPs_2Row.csv` +
`testset_metadata.csv` at verbose = 3: class dartR, 250 individuals,
255 loci, 30 populations, ploidy 2, maf and rdepth complete.

Caller grep across the 8 dartR-verse clones: no sibling package calls
`gl.read.dart` in code (doc examples only, happy path unchanged).

```json
{
  "function": "gl.read.dart",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ddaed27",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "applied", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DAT4", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "FS2", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "applied", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "deferred", "change": 5},
    {"id": "F6", "severity": "INFO", "confidence": "high", "rule": "FS9", "status": "no-action", "change": null}
  ],
  "coverage_skipped": [
    "fbm=TRUE path: environment-dependent, read-only review",
    "recalc/mono.rm combinations: dominated by separately-reviewed filters"
  ],
  "status": "pr-open",
  "pr": null
}
```
