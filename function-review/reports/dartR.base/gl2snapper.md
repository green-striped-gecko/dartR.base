# Review: gl2snapper (dartR.base)

- Family mode: io
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203
- Datasets: testset.gl, testset.gs
- Baseline: tests/testthat/test-gl2snapper.R (snapshot captured pre-review; 22 assertions passing)

## Verdicts

**Standards: Needs work** — the preamble and datatype rejection are correct
(`accept = c("genlight","SNP")` rejects SilicoDArT as documented), but the
preprocessing paths leak console output at `verbose = 0` and the write is not
sink-safe.
**Spec: Needs work** — the nexus output itself is well formed and verified
against the matrix (ntax/nchar, `symbols="012"`, `missing=?`, pop_ind taxon
labels, name de-spacing/de-duplication all confirmed), but a documented
default is wrong and the autapomorphy classifier misclassifies near-fixed
polymorphic populations.

What works well: the nexus body reproduces `as.matrix(x)` exactly, and the
name-sanitisation branch (spaces, duplicates) produces labels BEAUti can
parse.

## Findings

**F1 [MEDIUM, confidence: high] — preprocessing chatter defeats verbose = 0 (VRB5, VRB3)**
`R/gl2snapper.r:151,156,167,172-173,177-179` — the space/duplicate warnings
(`cat(warn(...))` at 151, 156) and the two `cat(report(...))` progress lines
(173, 179) are ungated, and `gl.allele.freq` (167) and `gl.drop.loc` (172)
are called without `verbose`, so they run at the global default of 2.
Failure scenario: `gl2snapper(x, rm.autapomorphies = TRUE, nloc = 100,
verbose = 0)` prints 13 lines (measured), including two "Completed:"
banners for other functions.
Proposed change: gate the four `cat()` calls at `verbose >= 2`, pass
`verbose = 0` (or the resolved `verbose`) to `gl.allele.freq` and
`gl.drop.loc`.

**F2 [MEDIUM, confidence: high] — documented default for rm.autapomorphies is wrong (DOC5, proposed rule)**
`R/gl2snapper.r:14-16` vs `:127` — roxygen says `[default TRUE]`; the
signature says `FALSE`.
Failure scenario: a user relying on the help page believes autapomorphies
are pruned by default and reports branch lengths from a dataset that still
contains them.
Proposed change: correct the roxygen to `[default FALSE]` (or flip the
default if pruning-by-default was the intent — custodian call).

**F3 [MEDIUM, confidence: medium] — autapomorphy classifier drops loci polymorphic in more than one population (DOC5, proposed rule)**
`R/gl2snapper.r:167-172` — a population counts as polymorphic only when its
allele frequency lies strictly inside (0.0001, 0.99) on the proportion
scale, computed from values `gl.allele.freq` has rounded to 2 dp on the
percent scale.
Failure scenario: a locus at 199/200 = 0.995 in population A and 0.5 in
population B is polymorphic in both, but A fails the `< 0.99` cut, the
locus tallies one population, and `rm.autapomorphies = TRUE` silently drops
a phylogenetically informative locus. Symmetrically, a frequency rounding
to 0.01% fails the `> 0.0001` cut.
Proposed change: classify on allele counts (`0 < alt < 2 * nobs`) rather
than thresholded rounded frequencies, or document the tolerance bounds in
`@details`.

**F4 [MEDIUM, confidence: high] — unguarded sink can hijack the console (no catalogue rule; gl2fasta F3 precedent)**
`R/gl2snapper.r:208-224` — `sink(outfilespec)` has no `on.exit(sink())`.
Failure scenario: any error or user interrupt between lines 208 and 224
(full disk, permission error) leaves all subsequent console output diverted
to the half-written nexus file for the rest of the session.
Proposed change: `on.exit(sink(), add = TRUE)` immediately after the
`sink()` call (idiom already applied to `gl2fasta` in this campaign).

**F5 [LOW, confidence: high] — NULL return is visible (FS10, VRB5)**
`R/gl2snapper.r:238` — `return(NULL)` prints `NULL` at the console on every
unassigned call, including at `verbose = 0` (verified via `withVisible`).
Proposed change: `invisible(NULL)` (as `gl2vcf` already does).

**F6 [LOW, confidence: high] — roxygen block drift (DOC1; DOC2, DOC6, DOC7 proposed rules)**
`R/gl2snapper.r:10-123` — `@param` order (`outfile, nloc,
rm.autapomorphies`) differs from the signature (`outfile,
rm.autapomorphies, nloc`); the `verbose` default clause is the outdated
"[default 2 or as specified using gl.set.verbosity]" wording; the `@author`
block names only a custodian (DOC7); `@return` sits after `@export`; typos
"suitable for for", "phylogentically".
Failure scenario: help-page readers get parameters in the wrong order and
the wrong verbosity contract.
Proposed change: reorder `@param` to match the signature, adopt the DOC2
default clause, add `Author(s): Arthur Georges.`, move `@return` before
`@author`, fix typos; run `devtools::document()` (DOC4).

## Proposed changes

1. Gate the four ungated `cat()` calls and pass `verbose` through to
   `gl.allele.freq` and `gl.drop.loc` (F1).
2. Correct the `rm.autapomorphies` documented default to FALSE, or flip the
   signature default with custodian approval (F2).
   **Consequence if the default flips: output loci change for default
   callers.**
3. Reimplement the autapomorphy test on allele counts instead of rounded
   frequency thresholds (F3).
   **Consequence: numerical output changes for datasets with near-fixed
   polymorphic populations under rm.autapomorphies = TRUE.**
4. Add `on.exit(sink(), add = TRUE)` around the nexus write (F4).
5. Replace `return(NULL)` with `invisible(NULL)` (F5).
6. Roxygen conformance pass + `devtools::document()` (F6).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT (n/a — no plot), STY — run
- Spec: nexus structure vs matrix on testset.gl subsets; SilicoDArT
  rejection; name sanitisation; nloc subsampling; rm.autapomorphies path —
  run empirically
- Loading the nexus into BEAUti/Snapper: SKIPPED — no BEAST2 install on the
  review machine; format checked structurally against the snapper template
  expectations only
- FBM path (DAT6): SKIPPED — no FBM fixture for converters

## Approval (Phase B)

Blanket approval by Arthur Georges, 2026-09-05, via the formal approval
boxes: all findings at every severity approved, explicitly acknowledging
that converted outputs change where they were wrong; the DAT7 class fix
(fatal `accept = "SNP"` gate wherever SilicoDArT is silently admitted) is
approved campaign-wide (gl2snapper already gates correctly, so no DAT7
change was needed here).

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges | 2026-09-05 |
| 2 | Approved | Arthur Georges | 2026-09-05; doc corrected to FALSE, signature unchanged |
| 3 | Approved | Arthur Georges | 2026-09-05; output-change consequence acknowledged |
| 4 | Approved | Arthur Georges | 2026-09-05 |
| 5 | Approved | Arthur Georges | 2026-09-05 |
| 6 | Approved | Arthur Georges | 2026-09-05 |

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-gl2snapper` (base upstream/dev,
ddaed27). All six changes applied:

1. F1 — name-mangling warnings and the two preprocessing `cat(report())`
   lines gated at `verbose >= 2`; `verbose = 0` passed to
   `gl.allele.freq` and `gl.drop.loc`.
2. F2 — roxygen `rm.autapomorphies` default corrected to
   `[default FALSE]`; signature unchanged.
3. F3 — autapomorphy classifier now scores a population polymorphic on
   allele counts (`0 < sum < 2 * nobs` from the popxloc table) instead
   of thresholded rounded frequencies.
4. F4 — sink guarded with a depth-recorded
   `on.exit(while (sink.number() > sink.depth) sink(), add = TRUE)`
   (a plain `on.exit(sink())` would warn "no sink to remove" on the
   success path, which closes the diversion explicitly).
5. F5 — `invisible(NULL)`.
6. F6 — roxygen reordered to house order with the DOC2 verbose clause
   and DOC7 author labels; `devtools::document()` run (gl2snapper.Rd
   only).

Verification: baseline characterization test passes (21 assertions; the
flipped baseline-bug expectations map to F1/F2/F5, and the two
ungated-warning assertions merged into one silence check under F1).
End-to-end run on testset.gl (`rm.autapomorphies = TRUE`, `nloc = 50`,
`verbose = 3`) clean: 54 autapomorphic loci removed, 50 retained, valid
nexus written. Sibling grep across the 8 clones: no callers of
gl2snapper outside dartR.base — all clear. NEWS.md entry added.

```json
{
  "function": "gl2snapper",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203482f4cc51dd2cef70f07493200e534ee7",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "applied", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "medium", "rule": "DOC5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "none:sink-safety", "status": "applied", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "FS10", "status": "applied", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 6}
  ],
  "coverage_skipped": ["BEAUti load test: no BEAST2 install", "DAT6: no FBM fixture"],
  "status": "pr-open",
  "pr": null
}
```
