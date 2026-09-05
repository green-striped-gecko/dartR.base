# Review: gl.read.silicodart (dartR.base)

- Family mode: io (SilicoDArT sibling of the SNP read chain; self-contained — does not use `utils.read.dart`/`utils.dart2genlight`)
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ddaed27 (upstream/dev); `R/gl.read.silicodart.r` identical on the empirical checkout ed99203
- Datasets: `testset_SilicoDArT.csv` + `testset_metadata_silicodart.csv` (dartR.data); synthetic fixtures `synth_silico_dupclone.csv`, `synth_silico_allpresent.csv`, no-pop metafile (layouts in the test file)
- Baseline: `tests/testthat/test-gl.read.silicodart.R` (new; snapshot captured pre-review, 17 expectations passing)

## Chain state (context, not findings)

On the testset the function returns a class `dartR` object (cast by the
closing `gl.compliance.check`), ploidy 1 on all individuals, genotypes
strictly 0/1/NA, `loc.metrics` rows 1:1 with the 255 loci (including
`Reproducibility`, `OneRatio`, `PIC`), `ind.metrics` rows 1:1 with
individuals, `@position`/`@chromosome` NULL. DAT1 (ploidy 1 for silico) and
DAT5 (ends compliant) hold. The 2-row question is moot: the function
hardcodes `nrows <- 1` ("there is no two row SilicoFormat??") and no 2-row
silico format is known.

## Verdicts

**Standards: Needs work** — preamble correct; verbosity contract broken
throughout (31 lines print at `verbose = 0`); one unused parameter.
**Spec: Needs work** — the happy path parses correctly, but a documented
metafile scenario always crashes on a variable-name typo, and duplicate
non-numeric CloneIDs silently corrupt locus names.

## Findings

**F1 [HIGH, confidence: high] — `pop(out)` typo: metafile without a pop column always errors (FS5; spec)**
`R/gl.read.silicodart.r:310` — the no-pop branch assigns
`pop(out) <- array(NA, nInd(glout))`; no object `out` exists.
Failure scenario (verified): the canonical silico fixture with a metafile
lacking `pop` dies with "object 'out' not found" immediately after printing
"Please note: there is no pop column". The branch has never worked. The SNP
path (`utils.dart2genlight`) handles the same case by assigning "pop1".
Proposed change: assign to `glout` and mirror the SNP path ("pop1" default
rather than NA, which `gl.compliance.check` would flag anyway); align the
warning text.

**F2 [HIGH, confidence: high] — duplicate non-numeric CloneIDs become NA locus names (DAT2)**
`R/gl.read.silicodart.r:199-214` — the uniquification loop assigns pasted
strings (`CLONEA_1`) into `covmetrics$CloneID`, which `read.csv(...,
stringsAsFactors = TRUE)` made a factor; new strings are not levels, so the
assignment produces NA.
Failure scenario (verified, `synth_silico_dupclone.csv`): both copies of a
duplicated alphanumeric CloneID end up NA in `locNames` and in
`loc.metrics$CloneID`; the object builds and passes compliance with NA locus
names. Numeric CloneIDs (as in the testset) escape only because their column
is not a factor.
Proposed change: convert `CloneID` with `as.character()` before the loop (or
read with `stringsAsFactors = FALSE` and re-factor selected columns).

**F3 [MEDIUM, confidence: high] — verbose = 0 is not silent (VRB5, VRB2)**
`R/gl.read.silicodart.r:105, 141-151, 158-159, 202-207, 255-346, 375` — most
progress and warning messages are ungated `cat()`/`cat(report())` calls
(line 105 and 158-159 are raw `cat`/`paste`, not even colour helpers), and
the closing `gl.compliance.check(glout)` is called without `verbose`.
Failure scenario (verified): `verbose = 0` prints 31 lines on the testset.
Proposed change: gate every message per VRB1-VRB4, route raw `cat()` through
the helpers, and pass `verbose` to `gl.compliance.check`.

**F4 [MEDIUM, confidence: high] — individuals absent from the metafile are dropped, and the message calls them loci (DOC5)**
`R/gl.read.silicodart.r:285-300` — same defect class as
`utils.dart2genlight` F1: `glout <- glout[ord2,]` filters the object to
metafile-matched individuals. The progress message says "Subsetting loci
now!".
Failure scenario (verified): 250 individuals in the csv, 218 in the metafile
— the returned object has 218; 32 genotype columns discarded with no
statement that removal happened.
Proposed change: same contract decision as the SNP path (join vs loud drop) —
apply the two consistently; fix the message wording either way.
**Consequence: if join semantics are adopted, output dimensions change for
partial-metafile calls (API1).**

**F5 [LOW, confidence: high] — `probar` parameter is documented but never used (DOC5)**
`R/gl.read.silicodart.r:23, 70` — `probar = TRUE` is accepted and documented
"Show progress bar"; no progress bar exists in the body.
Failure scenario: users pass `probar = FALSE` expecting effect (the
`@examples` in sibling docs do); silent no-op.
Proposed change: remove the parameter (loudly, per API2) or implement the
bar as in `utils.dart2genlight`.

**F6 [LOW, confidence: medium] — all-present (or all-absent) matrices are rejected as fatal (spec)**
`R/gl.read.silicodart.r:181-184` — the coding gate requires `max == 1 &&
min == 0`; a legitimate report in which every retained tag is present in
every individual (or a monomorphic subset) cannot be read.
Failure scenario (verified, `synth_silico_allpresent.csv`): "Fatal Error: Tag
P/A data must be 0 or 1!" on a valid 0/1-coded file containing only 1s.
Proposed change: reject only values outside {0, 1, NA}
(`!all(unlist(datas) %in% c(0, 1, NA))`) rather than requiring both values
to occur.

**F7 [LOW, confidence: high] — roxygen and message defects (DOC5, DOC2, DOC7 proposed)**
`R/gl.read.silicodart.r:2, 48, 55-58, 147` — title typo "\{agegenet\}";
`@details` claims the data are "combined into a genind object" (it builds a
genlight); the `@examples` call
`system.file('extdata', ind.metafile = 'testset_metadata_silicodart.csv', ...)`
passes `ind.metafile=` as a named argument into `system.file`, which works
only because `file.path()` drops names — the example misleads readers;
the duplicate-individual warning says "Rendering locus names unique" for
individuals; `verbose` doc deviates from DOC2 canon; `@author` has Custodian
only (DOC7, proposed).
Proposed change: correct the five fragments; `devtools::document()` (DOC4).

**F8 [INFO, confidence: high] — accidental-truthy length idioms**
`R/gl.read.silicodart.r:201, 335` — `if (length(index > 0))` and
`if (length(other.col > 0))` test the length of a logical vector; both work
only because non-zero length is truthy. Touch when neighbouring changes are
applied (STY1). Also noted: no service/plate extraction parity with the SNP
path (the silico header rows carry the same information) — parity is a
possible enhancement, not a defect.

## Proposed changes

1. Fix the `pop(out)` typo; default missing pop to "pop1" as in the SNP path
   (F1). **Consequence: previously-crashing calls now succeed — behaviour
   change for any workaround scripts (API1).**
2. Uniquify CloneID on character, not factor (F2).
3. Gate all messaging per VRB1-VRB4; pass `verbose` to `gl.compliance.check`
   (F3).
4. Align metafile semantics with the SNP path per the custodian's F1 decision
   on `utils.dart2genlight`; fix "Subsetting loci" wording (F4).
5. Remove or implement `probar` (F5).
6. Relax the 0/1 coding gate to set-membership (F6).
7. Roxygen corrections + `devtools::document()` (F7).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run (PLT: no plots,
  n/a)
- Spec: behaviour vs roxygen on the real fixture + 3 synthetic variants — run
- Ploidy-1 construction, presence/absence coding, Reproducibility metric
  carriage, topskip autodetection — run
- Explicit `topskip`/numeric `lastmetric` arguments: run by reading only
  (identical logic to `utils.read.dart`, exercised there)
- 2-row silicodart: SKIPPED — format does not exist (`nrows` hardcoded 1)
- FBM path (DAT6): SKIPPED — no fbm option on this function
- Duplicate individual names path (lines 141-151): run by reading only — the
  uniquified names ARE applied here (unlike `utils.read.dart`), no fixture
  run

## Approval (Phase B)

Approved 2026-09-05 by Arthur Georges via the formal approval boxes, covering
all BLOCKER/HIGH/MEDIUM findings for the read-chain review round, with
explicit acknowledgment of the consequence that objects users previously read
will differ where current behaviour is wrong (including the F1 consequence:
previously-crashing no-pop calls now succeed, API1). LOW findings were not
approved this round and are deferred; INFO items carry no action.

| Change | Decision | By | Note |
|---|---|---|---|
| 1 (F1, HIGH) | approved | Arthur Georges | consequence acknowledged: no-pop metafile calls now succeed with 'pop1' |
| 2 (F2, HIGH) | approved | Arthur Georges | duplicate CloneIDs get real _n names, not NA |
| 3 (F3, MEDIUM) | approved | Arthur Georges | verbose = 0 now fully silent |
| 4 (F4, MEDIUM) | approved | Arthur Georges | aligned with the SNP path per the utils.dart2genlight F1 decision (explicit loud drop); "Subsetting loci" wording removed |
| 5 (F5, LOW) | deferred | Arthur Georges | LOW findings not approved this round (probar untouched) |
| 6 (F6, LOW) | deferred | Arthur Georges | LOW findings not approved this round (0/1 coding gate untouched) |
| 7 (F7, LOW) | deferred | Arthur Georges | LOW findings not approved this round (roxygen untouched) |
| — (F8, INFO) | no action | | length-idiom lines left verbatim |

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-gl.read.silicodart` (base upstream/dev
ddaed27); PR #365 to `dev`
(https://github.com/green-striped-gecko/dartR.base/pull/365): changes 1-4 (F1-F4). Changes 5-7 (F5-F7, LOW) deferred - roxygen
untouched, so `man/gl.read.silicodart.Rd` and NAMESPACE are unchanged
(`devtools::document()` run to confirm). Change 4 follows the custodian
contract applied on the SNP path (`utils.dart2genlight` change 1): the drop
is retained, announced with a "REMOVED" warning at `verbose >= 1`, and the
"Subsetting loci now!" wording is gone; the deferred-LOW message texts
(F7's "Rendering locus names unique" for individuals) were left verbatim
inside the newly gated blocks.

Addendum (fixture correction, below-HIGH): the F4 failure scenario states
"250 individuals in the csv, 218 in the metafile". On dartR.data 1.2.5 the
`testset_SilicoDArT.csv` fixture itself holds 218 individuals and the
metafile lists the same 218, so no drop occurs on the canonical pair and
the Phase A "verified 250 -> 218" figure cannot have come from this
fixture. The defect class is real and the fix is exercised in the baseline
test with a 100-row partial metafile (118 individuals dropped, loudly).

Verification: baseline characterization test run against the branch -
22 expectations, 0 failures. Diffs from the Phase A baseline all map to
approved findings and are marked `# [approved F<n>]` in the test file:
the no-pop metafile now yields a 'pop1' object instead of "object 'out'
not found" (F1), duplicate alphanumeric CloneIDs come back as CLONEA_1/
CLONEA_2 instead of NA locus names (F2), and `verbose = 0` prints nothing
(F3). LOW-pinned snapshots unchanged: the all-present matrix is still
rejected ("must be 0 or 1", F6 deferred) and `probar` remains inert (F5
deferred). End-to-end run on `testset_SilicoDArT.csv` +
`testset_metadata_silicodart.csv` at verbose = 3: class dartR,
218 individuals, 255 loci, 29 populations, ploidy 1, genotypes all in
{0, 1, NA}, loc.metrics/ind.metrics registered 1:1.

Caller grep across the 8 dartR-verse clones: no sibling package calls
`gl.read.silicodart` in code.

```json
{
  "function": "gl.read.silicodart",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ddaed27",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "FS5", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DAT2", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "deferred", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "medium", "rule": "none:input-validation", "status": "deferred", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "deferred", "change": 7},
    {"id": "F8", "severity": "INFO", "confidence": "high", "rule": "STY1", "status": "no-action", "change": null}
  ],
  "coverage_skipped": [
    "2-row silicodart: format does not exist",
    "DAT6 FBM: no fbm option",
    "duplicate individual names: read-only, uniquification is applied here"
  ],
  "status": "pr-open",
  "pr": 365
}
```
