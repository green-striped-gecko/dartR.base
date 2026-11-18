# Review: utils.dart2genlight (dartR.base)

- Family mode: io (DArT read chain, called by `gl.read.dart`)
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ddaed27 (upstream/dev) reviewed by reading. Empirical runs on checkout ed99203 (integration-local), which carries PR #330's hunks for this file — see "PR #330 interaction" below.
- Datasets: `testset_SNPs_2Row.csv` + `testset_metadata.csv` (dartR.data); synthetic 2-row fixtures (see `utils.read.dart.md` for layouts)
- Baseline: `tests/testthat/test-utils.dart2genlight.R` (new; snapshot captured pre-review, 23 expectations passing)

## PR #330 interaction (context, not findings)

Open PR #330 (branch `position-genome-only`) removes the
`pos <- sraw$SnpPosition[esl]` / `position = pos` pair from this file
(upstream/dev lines 118 and 171) so that `@position`/`@chromosome` hold genome
coordinates only, with the tag offset kept in `loc.metrics$SnpPosition`. The
review checkout already carries those hunks; empirically the chain now
delivers `@position = NULL`, `@chromosome = NULL`, and
`loc.metrics$SnpPosition` populated (verified on the testset:
first five SnpPosition values 19, 38, 24, 33, 21). No findings are raised
against those lines. Any approved change touching lines 114-172 must be
applied on top of #330 (conflicts-with #330 if applied to upstream/dev
directly).

## Verdicts

**Standards: Needs work** — FS2/FS3 preamble present; ploidy (DAT1) and
loc.metrics/ind.metrics registration (DAT2) verified correct on the testset;
verbosity gating and one error idiom fail.
**Spec: Needs work** — genotype encoding is correct for both formats (verified
cell-by-cell against the raw csv), but the ind.metafile behaves as a filter,
not a metadata join: individuals without a metafile row are removed, which the
roxygen never states.

## Encoding verification (context)

Hand-read cells from `testset_SNPs_2Row.csv` match the object exactly:
`AA010915` and `UC_00126`, loci 1-4, raw pairs `0/1`, `-/-`, `1/0`, `1/0`
convert to `2, NA, 0, 0`. Synthetic 1-row raws `0/1/2` convert to `0/2/1`
(hom ref / hom SNP / het). Hom-ref = 0 and hom-alt = 2 agree between the
1-row and 2-row paths; ploidy 2 in both. The historic inverted-coding trap is
not present at this commit.

## Findings

**F1 [HIGH, confidence: high] — individuals absent from ind.metafile are silently dropped (DOC5; spec)**
`R/utils.dart2genlight.r:271-293` (upstream numbering) — `ord` keeps only
DArT ids found in the metafile and `gout <- gout[ord2, ]` subsets the object
to them. The roxygen sells `ind.metafile` as metadata ("Optional file ... with
metadata for each individual"); nothing states that genotypes are discarded.
Failure scenario (verified): the 250-individual testset read with a
100-row metafile returns 100 individuals; 150 genotype columns vanish with
only a warning ending "Maybe this is fine if a subset matches."
Proposed change: either keep unmatched individuals with NA metadata (join
semantics) or make the drop explicit — a `verbose >= 1` warning stating "N
individuals in the DArT file had no metafile row and were REMOVED", plus a
`@details` sentence. Custodian to choose the contract.
**Consequence: if join semantics are adopted, output dimensions change for
partial-metafile calls (API1).**

**F2 [MEDIUM, confidence: high] — missing SNP/Variant column fails with an opaque closure error (FS5)**
`R/utils.dart2genlight.r:125-135` — `alleles` is only assigned when a `SNP`
or `Variant` column exists. When neither does, the name resolves to
adegenet's `alleles()` function and `substr()` dies with "cannot coerce type
'closure' to vector of type 'character'". The explicit header check that used
to catch this is commented out at lines 86-90.
Failure scenario (verified): a report with the SNP column renamed produces
the closure error with no hint about headers.
Proposed change: reinstate a fail-fast check — if neither `SNP` nor
`Variant` is present, `stop(error(...))` naming the expected columns.

**F3 [MEDIUM, confidence: high] — 2-row translation coerces unlisted patterns to NA via warnings (DAT1)**
`R/utils.dart2genlight.r:138-150` — only `0/1`, `1/0`, `1/1`, `NA/NA` are
mapped; anything else (`0/0` double-null, partial `NA/1`/`1/NA`, or any
unexpected token in a malformed file) falls through to
`as.numeric()`, becoming NA with a raw "NAs introduced by coercion" warning
per affected individual.
Failure scenario (verified): a `0/0` pair converts to NA — defensible as a
missing call, but indistinguishable from data corruption; a file with a stray
"2" in one row would likewise dissolve into NAs instead of erroring.
Proposed change: map `0/0` and partial-NA pairs explicitly to NA (silencing
the spurious warning), and count-and-report any residual unmatched patterns
as a gated warning (or error above a threshold).

**F4 [MEDIUM, confidence: high] — duplicate metafile ids abort with an empty error message (VRB2, FS5)**
`R/utils.dart2genlight.r:240-245` — `cat(error("Individual names are not
unique..."))` followed by a bare `stop()`. The condition message the user's
tooling sees is empty; `tryCatch` handlers and logs capture nothing.
Failure scenario: any metafile with a repeated id; batch pipelines record an
error with no message.
Proposed change: `stop(error("Fatal Error: individual names in the
ind.metafile are not unique\n"))` in the house idiom.

**F5 [MEDIUM, confidence: high] — verbose = 0 is not silent (VRB5)**
`R/utils.dart2genlight.r:192-209, 249-268` — the TrimmedSequence fallback
warnings and the id-mismatch warnings are ungated `cat(warn(...))`/`print()`
calls.
Failure scenario (verified): at `verbose = 0` the TrimmedSequence warning
prints (no trailing newline, so it concatenates into the next line of
output), and a subset metafile prints the full unmatched-id listing — 150
lines on the test run.
Proposed change: gate result-affecting warnings at `verbose >= 1`, purely
informational ones at `verbose >= 2`; add the missing `\n` to the
TrimmedSequence warnings.

**F6 [LOW, confidence: high] — roxygen defects (DOC5, DOC7 proposed)**
`R/utils.dart2genlight.r:11, 27, 30-31` — `covfilename` documented
"Depreciated" (spelling) with no default bracket; `@author` uses
"Maintainer:" where the house structure is "Author(s): ... Custodian: ..."
(DOC7, proposed rule); `@return` is garbled ("Including all available slots
are filled.") and does not mention `loc.metrics`/`ind.metrics`/`service`/
`plate_location`, which the function does attach.
Proposed change: rewrite the three doc fragments; run `devtools::document()`
(DOC4).

**F7 [INFO, confidence: high] — accidental-but-working idioms**
`R/utils.dart2genlight.r:284` `length(ord == nind)` equals `length(ord)` —
the message is right by accident; `:179` `x <- NULL` is dead; `:37` `probar`
default TRUE here but FALSE in the `gl.read.dart` wrapper (harmless
inconsistency). No change proposed beyond touching these lines when F5 is
applied.

## Proposed changes

1. Decide and implement the ind.metafile contract: join semantics (keep
   unmatched individuals, NA metadata) or explicit loud drop + documentation
   (F1). **Consequence: output dimensions change for partial-metafile calls
   if join semantics are chosen (API1).**
2. Reinstate the SNP/Variant header check with a clear fatal error (F2).
3. Make the 2-row translation table exhaustive; report unmatched patterns
   (F3).
4. Replace `cat(error(...)); stop()` with `stop(error(...))` (F4).
5. Gate ungated warnings per VRB3/VRB4; fix missing newlines (F5, F7 lines).
6. Roxygen corrections + `devtools::document()` (F6).

All changes to this file must be sequenced after PR #330 merges
(conflicts-with #330 in the 114-172 region for change 3).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, STY — run
- Spec: behaviour vs roxygen on real 2-row fixture, synthetic 1-row and
  malformed variants — run
- Genotype encoding vs raw csv by hand: run (both formats)
- DAT1 ploidy: run — ploidy 2 confirmed on every path exercised
- DAT2 registration: run — loc.metrics rows == nLoc, ind.metrics rows == nInd
  after every metafile scenario
- `Variant`/`TrimmedSequenceSnp`/`AlleleSequenceSnp` fallback columns:
  SKIPPED — no recent-format fixture with those headers; branches read but
  not executed
- FBM path (DAT6): SKIPPED — conversion happens in `gl.read.dart`
- PLT: not applicable

## Approval (Phase B)

Approved 2026-09-05 by Arthur Georges via the formal approval boxes, covering
all BLOCKER/HIGH/MEDIUM findings for the read-chain review round, with
explicit acknowledgment of the consequence that objects users previously read
will differ where current behaviour is wrong. LOW findings were not approved
this round and are deferred; INFO items carry no action.

| Change | Decision | By | Note |
|---|---|---|---|
| 1 (F1, HIGH) | approved | Arthur Georges | applied as the explicit-loud-drop branch of the proposed change (see Outcome); the join-semantics alternative remains open as a future custodian decision |
| 2 (F2, MEDIUM) | approved | Arthur Georges | clear header error replaces the closure error |
| 3 (F3, MEDIUM) | approved | Arthur Georges | exhaustive translation table; unmatched patterns reported |
| 4 (F4, MEDIUM) | approved | Arthur Georges | house stop(error(...)) idiom |
| 5 (F5, MEDIUM) | approved | Arthur Georges | verbose = 0 now fully silent; F7 lines touched as scoped |
| 6 (F6, LOW) | deferred | Arthur Georges | LOW findings not approved this round |
| — (F7, INFO) | folded into change 5 | | length(ord == nind), dead x <- NULL, missing newlines |

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-utils.dart2genlight` (base upstream/dev
ddaed27); PR #363 to `dev`
(https://github.com/green-striped-gecko/dartR.base/pull/363): changes 1-5 (F1-F5, with F7's line touches folded into change 5 as
scoped). Change 6 (F6, LOW) deferred.

Contract decision on change 1: the proposed change offered join semantics
(keep unmatched individuals, NA metadata - an API1 dimension change) or an
explicit loud drop. The blanket approval did not record a choice of the
join alternative, so the conservative branch was applied: the drop is
retained (dimensions unchanged from baseline), announced with a
"REMOVED" warning at `verbose >= 1`, and documented in `@details`
(`man/utils.dart2genlight.Rd` regenerated). Adopting join semantics
remains open as a follow-up custodian decision.

PR #330 disposition: the branch is plain upstream/dev - no #330 hunks were
carried. Checked against the actual #330 diff: its two hunks on this file
(old lines 114-121 and 165-171) share no lines with this branch's hunks
(nearest approach: two unchanged lines between #330's `position = pos`
removal and the change-5 region), so both merge cleanly in either order.
The baseline test's @position expectation is written state-agnostically
(pre-#330: @position mirrors `loc.metrics$SnpPosition`; post-#330: NULL),
so it passes on this branch and after #330 lands.

Verification: baseline characterization test run against the branch -
25 expectations, 0 failures. Diffs from the Phase A baseline all map to
approved findings and are marked `# [approved F<n>]` in the test file:
the partial-metafile drop is unchanged in effect but now loud (F1) and
silent at verbose = 0 (F5), the missing-SNP-column error is now the clear
header message (F2), and 0/0 pairs convert to NA without coercion warnings
(F3). End-to-end run on `testset_SNPs_2Row.csv` + `testset_metadata.csv`
at verbose = 3: class dartR, 250 individuals, 255 loci, 30 populations,
ploidy 2, loc.metrics 255 rows, ind.metrics 250 rows, hand-checked cells
2,NA,0,0 unchanged.

Caller check: `utils.dart2genlight` is called only by `gl.read.dart`
(grep across the 8 dartR-verse clones); the full chain was re-run with the
fixed util (250 x 255, ploidy 2). `gl.read.silicodart` is self-contained
and does not call this util (verified by grep).

```json
{
  "function": "utils.dart2genlight",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ddaed27",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DAT1", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "VRB2", "status": "applied", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "applied", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "deferred", "change": 6},
    {"id": "F7", "severity": "INFO", "confidence": "high", "rule": "STY1", "status": "applied", "change": 5}
  ],
  "coverage_skipped": [
    "Variant/TrimmedSequenceSnp fallback branches: no recent-format fixture",
    "DAT6 FBM: out of scope for this function"
  ],
  "interactions": {"pr330": "no shared lines with #330 hunks (checked against the PR diff); merges cleanly in either order; no #330 hunks carried on the branch"},
  "status": "pr-open",
  "pr": 363
}
```
