# Review: gl2bpp (dartR.base)

- Family mode: io
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203
- Datasets: testset.gl, testset.gs (rejection probe), derived subsets
- Baseline: tests/testthat/test-gl2bpp.R (snapshot captured pre-review; 18 assertions passing)

## Verdicts

**Standards: Needs work** — the function gets the datatype restriction right (`accept = "SNP"`, SilicoDArT rejected with a clear message — the only converter in this batch that does) and uses the tag-offset `SnpPosition` correctly, but it guards a dependency it never uses, validates neither `method` nor the presence of `loc.metrics.flags`, and redirects the console without protection.

**Spec: Needs work** — on dartR-class SNP input, methods 1 and 2 produce a correct phylip-per-locus alignment (headers, ambiguity codes, all-N rows for missing genotypes all verified), but a plain genlight input silently yields genetically wrong sequences, and the `merge.secondaries` block as written cannot produce a valid merged file.

## Findings

**F1 [HIGH, confidence: high] — plain genlight input silently produces wrong sequences (DAT2/DAT3, DAT5)**
`R/gl2bpp.r:166` — `x <- x[, order(x@other$loc.metrics$SnpPosition)]` relies on the `[` method co-reordering `loc.metrics`. The dartR-class method does (verified); adegenet's genlight method does not.
Failure scenario: verified 2026-09-05 — the same 30-locus object run once as dartR class and once as plain genlight (the class any adegenet/vcfR-built object has) differs on 8220 of 8250 output lines; locus names are paired with the wrong `TrimmedSequence`, so the alignment looks valid and is genetically wrong. No warning is emitted.
Proposed change: re-subset the metadata explicitly from the pre-sort object per DAT3 (`ord <- order(...); x <- x[, ord]; x@other$loc.metrics <- old.metrics[ord, , drop = FALSE]` where needed), or require/coerce dartR class up front.

**F2 [HIGH, confidence: medium] — `merge.secondaries` produces a structurally invalid file (DOC5, proposed rule; code trace)**
`R/gl2bpp.r:366-453` — four defects in the merge block:
(a) line 433 replaces sequences but the retained block header (written earlier) still carries the pre-merge sequence length;
(b) line 437 deletes only `min(dell)-1` — one header line — so a clone with three or more secondaries leaves orphan `"<n> <l>"` header lines with no sequences;
(c) line 427 appends `substr(seq, snppos[ii-1]+1, nchar(seq))` after the loop already appended up to `snppos[ii]`, duplicating the segment between the last two SNP positions;
(d) line 430 reads names from `bppf` (the original file) at indices computed against `b2`, which shrinks after each clone's deletions — labels for the second and later merged clones come from the wrong lines.
Failure scenario: any dataset with retained secondaries and `merge.secondaries = TRUE` gets a file BPP either refuses or misparses (header/sequence length disagreement).
Confidence is medium because the path could not be executed: testset.gl contains no clone with more than one retained locus (checked), so only the no-secondaries no-op path was run (verified harmless). The defects are from code trace.
Proposed change: rewrite the merge block — recompute each merged block's header from the merged length, delete all superseded headers, take the tail from `snppos[ii]+1`, and index names from the current `b2`.

**F3 [MEDIUM, confidence: high] — `method` never validated; invalid value silently writes no alignment (FS5)**
`R/gl2bpp.r:76,172,278` — only `method == 1` and `method == 2` write output; anything else skips both branches.
Failure scenario: verified 2026-09-05 — `method = 3` runs to "Completed", writes the Imap only, and no alignment file exists. A typo costs the user a silent empty result.
Proposed change: fail fast (`stop(error(...))`) or clamp to 1 with a warning, per the house error-handling pattern.

**F4 [MEDIUM, confidence: high] — missing `loc.metrics.flags` crashes with an opaque error (DAT5)**
`R/gl2bpp.r:116` — `if (x@other$loc.metrics.flags$monomorphs == FALSE)` on an object without flags evaluates `if (logical(0))`.
Failure scenario: verified 2026-09-05 — a genlight lacking `loc.metrics.flags` (any non-dartR-built object) fails with "argument is of length zero" instead of a usable message.
Proposed change: guard with `is.null()` (treat absent flags as "monomorphs not confirmed removed") or require `gl.compliance.check`.

**F5 [MEDIUM, confidence: high] — dependency guard for a package the function never uses, with a non-stopping idiom (DEP1)**
`R/gl2bpp.r:105-113` — `seqinr` is checked but no `seqinr::` call exists anywhere in the function; on failure the code does `cat(error(...)); return(-1)` instead of `stop()`.
Failure scenario: a user without seqinr is refused for no reason, and a script gets `-1` back and continues as if a file had been written.
Proposed change: delete the guard; if a guard were ever needed, use the DEP1 `stop(error(...))` idiom.

**F6 [LOW, confidence: medium] — vectorized `cat` puts leading/trailing spaces on sequence and Imap lines (STY1)**
`R/gl2bpp.r:268,356,459` — `cat(paste0(...vector..., "\n"))` joins elements with a space, so every line after the first in a block begins with a space, and Imap lines carry trailing blanks (verified in output). BPP tokenizes on whitespace, so this likely parses, but the files are noisy and line-diff unfriendly.
Proposed change: `cat(..., sep = "")` or `writeLines`.

**F7 [LOW, confidence: high] — visible NULL return (FS10)**
`R/gl2bpp.r:468` — `return(NULL)` prints `NULL` on an unassigned call (verified with `withVisible`).
Proposed change: `return(invisible(NULL))`.

**F8 [LOW, confidence: high] — roxygen order, dead locals (DOC1, STY1)**
`R/gl2bpp.r:1-74,174-176` — `@details` precedes `@author` but `@return` sits last after `@export` (old order); `@author` names a custodian only (DOC7, proposed rule); `allnames` (line 174) and the method-1 `trimmed` vector (line 176) are assigned and never used.
Proposed change: roxygen pass; delete the dead assignments.

## Proposed changes

1. Make the internal locus reorder metadata-safe for plain genlight input (F1). **Consequence: output changes (becomes correct) for non-dartR genlight callers.**
2. Rewrite the `merge.secondaries` block: correct headers, full deletion of superseded blocks, non-overlapping segment concatenation, current-file name indexing (F2). **Consequence: merged output changes for any dataset with secondaries.**
3. Validate `method` at entry; error or clamp with a warning (F3).
4. Null-guard `loc.metrics.flags$monomorphs` (F4).
5. Remove the vestigial seqinr guard (F5).
6. Emit lines via `writeLines`/`cat(sep = "")` to drop stray spaces (F6).
7. Return `invisible(NULL)`; roxygen pass; delete dead locals (F7, F8).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT (n/a — no plot), STY — run
- Spec method 1: run on the roxygen-example pipeline (callrate 1, monomorphs removed,
  50 loci, seed 42) — 45 blocks of 274+1 lines after the internal overshoot filter;
  header, ambiguity coding, and ^-tag naming vs the BPP/Imap contract verified
- Spec method 2: run (seeded) — same shape verified
- All-missing locus behaviour: verified — genuinely all-NA loci yield all-N rows for
  every individual (correct), confirmed against the genotype matrix
- SNP vs SilicoDArT dispatch: run — SilicoDArT rejected with a clear fatal error
- Genome-position lens: conforms — positions come from `loc.metrics$SnpPosition`
  (tag offset) as the convention requires; `@position`/`@chromosome` untouched
- verbose = 0 silence via `capture.output` — run (0 lines)
- `merge.secondaries` with real secondaries: SKIPPED — testset.gl retains no multi-locus
  clones; only the no-secondaries no-op path was executed (F2 is code-trace)
- FBM path (DAT6): SKIPPED — no FBM fixture
- BPP executable round-trip: SKIPPED — BPP binary unavailable

## Approval (Phase B)

All findings at every severity approved 2026-09-05 by Arthur Georges via the
formal approval boxes, as part of the blanket class approval for the io
converter batch: converted outputs change where they were wrong, and the DAT7
class fix applies (already satisfied here -- gl2bpp had `accept = "SNP"`).

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges | Output changes (becomes correct) for plain genlight callers |
| 2 | Approved | Arthur Georges | Merged output changes for datasets with secondaries |
| 3 | Approved | Arthur Georges | Fatal error chosen (not clamp) |
| 4 | Approved | Arthur Georges | Absent flags treated as monomorphs-not-confirmed |
| 5 | Approved | Arthur Georges | |
| 6 | Approved | Arthur Georges | |
| 7 | Approved | Arthur Georges | |

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-gl2bpp` (base `upstream/dev` = ddaed27).
All seven changes applied. Base note: the review was performed at ed99203,
which includes the position-genome-only audit (PR #330) moving gl2bpp's
position source from `@position` to `loc.metrics$SnpPosition`; the
upstream/dev version cannot run at all on dartR.data >= 1.2.3 (`testset.gl`
has an empty `@position` slot), so this branch carries the same gl2bpp hunks
as PR #330 -- they become a no-op merge once #330 lands. The merge.secondaries
rewrite is block-wise: block j of the file maps to locus j of the
position-sorted object, each merged block's header is recomputed from the
merged length, all superseded blocks are deleted, segments concatenate
without overlap, and labels come from the clone's first block. Baseline
characterization test updated for the four approved behaviour changes (F1
class-identical output, F3 method validation, F4 flag guard, F6 clean lines)
and extended with a synthetic-secondaries property test (testset.gl carries
no multi-locus clones): 3 loci in, 2 blocks out, header `4 40`, merged
sequence exactly seq1[1..6] + seq2[7..21] + seq2[22..40] per individual,
singleton block untouched. 25 assertions pass; SNP output on the roxygen
example pipeline is unchanged (12375 alignment lines, header `274 69`).
End-to-end runs at `verbose = 3` for methods 1 and 2 verified. Sibling-clone
grep across the 8 dartR-verse packages: no callers of `gl2bpp` -- all clear.
NEWS.md entry added.
PR: green-striped-gecko/dartR.base#347.

```json
{
  "function": "gl2bpp",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203482f4cc51dd2cef70f07493200e534ee7",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT3", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "medium", "rule": "DOC5", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DAT5", "status": "applied", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DEP1", "status": "applied", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "medium", "rule": "STY1", "status": "applied", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "FS10", "status": "applied", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 7}
  ],
  "coverage_skipped": ["merge.secondaries with real secondaries: no multi-locus clones in testset.gl", "DAT6: no FBM fixture", "BPP executable round-trip: binary unavailable"],
  "status": "pr-open",
  "pr": 347
}
```
