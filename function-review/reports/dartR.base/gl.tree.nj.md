# Review: gl.tree.nj (dartR.base)
- Family mode: analysis
- Date: 2026-09-06
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ddaed27 (upstream/dev; working-copy `R/gl.tree.nj.r` verified identical by `git diff upstream/dev -- R/gl.tree.nj.r`, so `devtools::load_all()` exercised the reviewed state)
- Datasets: testset.gl (SNP), testset.gs (SilicoDArT) — dartR.data 1.2.5
- Baseline: tests/testthat/test-gl.tree.nj.R (new file, snapshot captured pre-review; 30 assertions, all passing)

## Verdicts

**Standards: Needs work** — the FS backbone (verbosity, flag start, datatype
check, flag end, explicit return) is present and correctly ordered, and the
input object comes back untouched; but plotting is unconditional and coupled
to the return path (PLT3, VRB5), one warning prints ungated at `verbose = 0`
(VRB3), and the roxygen block has a batch of accuracy gaps.

**Spec: Needs work** — the core computation matches its own formula exactly
(independent recomputation: RF distance 0, identical edge sums), but three
documented behaviours fail in practice: the default distance silently discards
23.8% of the frequency matrix through missing-data handling inconsistent with
the cross-referenced `gl.dist.pop`; `by.pop = FALSE` errors on any ordinary
dataset; and the UPGMA method is only reachable through a misspelling, with
the correct spelling silently returning an NJ tree.

## Independent verification (spec axis)

The function computes: per-population allele frequencies
`mean(genotypes)/2` per locus (no `na.rm`), `stats::dist` (Euclidean) on the
population-by-locus frequency matrix, rounded to 4 dp, then `ape::nj`
(or `stats::hclust(method = "average")` + `as.phylo` for "ugpma").
Recomputing that pipeline independently on `testset.gl` and running
`ape::nj` on it reproduces the returned tree exactly: RF distance 0,
edge-length sums equal (14.33359), tip labels = `popNames`. Same result on
`testset.gs`. So the code does what the code says; the findings below are
where that diverges from what the documentation and cross-references say.

## Findings

**F1 [HIGH, confidence: high] — default distance discards data at any locus with a missing genotype (DOC5 (proposed rule), API1 (proposed rule))**
`R/gl.tree.nj.r:120-121` — frequencies are computed with
`mean(e) / 2` and no `na.rm`, so one missing genotype in a population
nullifies that population's frequency at that locus; `stats::dist` then
excludes the NA'd locus pairwise and rescales.
Failure scenario: on `testset.gl`, 23.8% of the 30x255 frequency matrix is NA
(only 77 of 255 loci are complete in every population). Rebuilding the matrix
with `na.rm = TRUE` (the policy `gl.allele.freq`/`gl.dist.pop` use — the
function the `@details` section points users to for alternative matrices)
changes mean pairwise distance from 1.62 to 1.95 and the topology by RF
distance 32. `gl.tree.nj(x)` and
`gl.tree.nj(x, dist.matrix = gl.dist.pop(x, method = "euclidean"))` therefore
give materially different trees with no hint to the user.
Proposed change: `mean(e, na.rm = TRUE) / 2`.
**Consequence: numerical output changes for any dataset with missing data —
that is, essentially all callers using the default distance.**

**F2 [HIGH, confidence: high] — `by.pop = FALSE` errors on any ordinary dataset (DOC5 (proposed rule))**
`R/gl.tree.nj.r:90` — `popNames(x) <- indNames(x)` assigns `nInd` names to
the `nPop` levels of the population factor.
Failure scenario: `gl.tree.nj(testset.gl, by.pop = FALSE)` stops with
adegenet's "Vector length does no match number of populations" — the
documented individual-tree mode is unusable whenever `nInd != nPop`, i.e.
always in practice.
Proposed change: `pop(x) <- indNames(x)`.

**F3 [HIGH, confidence: high] — UPGMA reachable only via the misspelling "ugpma"; correct spelling silently returns an NJ tree (DOC5 (proposed rule), VRB3)**
`R/gl.tree.nj.r:82-85, 134` — the accepted value is `"ugpma"` (docs spell it
"UGPMA" too); anything else, including the correct `"UPGMA"`/`"upgma"`,
triggers a fallback to `"nj"`.
Failure scenario: `gl.tree.nj(testset.gl, method = "UPGMA")` returns the NJ
tree — which differs from the actual UPGMA tree by RF distance 24 on
`testset.gl` — and the fallback warning is not verbosity-gated, so it prints
even at `verbose = 0` (the only console output the function emits at 0). A
user who scripts with the correct spelling and `verbose = 0` in a pipeline
gets the wrong algorithm's tree with the warning buried.
Proposed change: accept `"upgma"` (keep `"ugpma"` as a silent back-compat
alias), fix the spelling in the roxygen, and gate the fallback warning at
`verbose >= 1` per VRB4 since it changes the result.

**F4 [MEDIUM, confidence: high] — plotting is unconditional and its failure destroys the computed tree (PLT3, VRB5)**
`R/gl.tree.nj.r:143-167` — `ape::plot.phylo` runs on every call; there is no
`plot.display`/`plot.out` parameter and no `verbose == 0` gate, and `type` is
validated only inside ape after the tree is built.
Failure scenario: (a) `verbose = 0` still renders a plot, violating the
silence contract; (b) `gl.tree.nj(x, type = "banana")` errors inside
`plot.phylo` after the tree exists, so the user loses the result; (c) any
headless environment without a device gets plot-side failures on a
computation that succeeded.
Proposed change: add `plot.display = TRUE` gated off at `verbose == 0` (house
preamble idiom), validate `type` with `match.arg` during error checking, and
plot inside a non-fatal wrapper so the tree is returned regardless.

**F5 [MEDIUM, confidence: high] — outgroup rooting returns a tree that is not rooted in ape's sense (DOC5 (proposed rule))**
`R/gl.tree.nj.r:146` — `ape::root(tree, outgroup)` with the default
`resolve.root = FALSE` leaves a basal trichotomy.
Failure scenario: `tr <- gl.tree.nj(testset.gl, outgroup = "EmmacBrisWive")`
gives `ape::is.rooted(tr) == FALSE`; downstream tools that require a rooted
tree (and any `treefile` written in the same call) reject or misread it.
Proposed change: `ape::root(tree, outgroup, resolve.root = TRUE)`.

**F6 [LOW, confidence: high] — `as.pop` error message points at the wrong slot (VRB2)**
`R/gl.tree.nj.r:108` — the fatal error for an unknown `as.pop` says "Check
names(gl@other$loc.metrics)" but the lookup is in `ind.metrics`.
Failure scenario: a user mistyping an individual metric is sent to inspect
locus metrics, where the metric can never be.
Proposed change: reword to `names(x@other$ind.metrics)`.

**F7 [LOW, confidence: high] — SilicoDArT admitted with SNP-ploidy arithmetic (DAT7)**
`R/gl.tree.nj.r:70, 121` — the default `accept` admits SilicoDArT (and the
`@examples` run `testset.gs`), but `mean(e) / 2` halves presence/absence
frequencies (ploidy 1).
Failure scenario: the scaling is uniform, so topology is unaffected
(verified: RF 0 against the same-formula reference on `testset.gs`), but all
edge lengths are half the true P/A-frequency Euclidean distances — wrong
scale for anyone reading branch lengths off a SilicoDArT tree.
Proposed change: divide by `ploidy(x)[1]` (or by 1 for SilicoDArT), or
document the scale explicitly.

**F8 [LOW, confidence: medium] — a plain-matrix `dist.matrix` breaks the ugpma path with an opaque error (DAT5)**
`R/gl.tree.nj.r:131, 136` — `dist.matrix` is documented as "Distance matrix";
`ape::nj` accepts a matrix but `stats::hclust` requires a `dist` object.
Failure scenario:
`gl.tree.nj(x, dist.matrix = as.matrix(d), method = "ugpma")` fails with
"missing value where TRUE/FALSE needed" from inside `hclust` — no hint the
input class is the problem.
Proposed change: coerce with `d <- as.dist(dist.matrix)` (or validate class)
when a user matrix is supplied.

**F9 [LOW, confidence: high] — roxygen accuracy batch (DOC1, DOC2, DOC5 (proposed rule), DOC7 (proposed rule), DOC4)**
`R/gl.tree.nj.r:1-48` — (a) `@return` claims "A tree file of class phylo";
the return is a `phylo` object, a file only when `treefile` is set;
(b) `@param` order swaps `type`/`outgroup` relative to the signature;
(c) the `verbose` text does not match the DOC2 canon; (d) tag order is the
outdated one (`@details` after `@param`, `@return` after `@export`) — per
DOC1's note, flag the file; (e) `@author` has a Custodian but no
`Author(s):` line (DOC7, proposed rule); (f) `@param method` carries the
"UGPMA" misspelling (see F3) and `@param type` omits ape's accepted
"radial"/"tidy"; (g) `@description` says "hclust ... applied to Euclidean
distances" — hclust applies only under `method = "ugpma"`.
Failure scenario: users act on the documented contract (e.g. spell UPGMA
correctly, pass `type` values from ape's docs, expect a file back) and hit
F3/F4-class surprises.
Proposed change: one documentation pass fixing a-g, then
`devtools::document()` (DOC4).

**F10 [INFO, confidence: high] — minor code-quality notes (STY1, FS10, VRB3, FS5)**
`R/gl.tree.nj.r` — (a) line 126: distances are rounded to 4 dp before tree
building — deliberate or leftover display code, worth a comment either way;
(b) line 120: local `t` shadows `base::t`; (c) line 187: visible
`return(tree)` prints the phylo summary when the documented example calls the
function unassigned — house style is `invisible()` for that case;
(d) lines 74-78: the dartR-cast warning is gated at `verbose > 2` where VRB3
warnings sit at `>= 2`; (e) lines 148-166: the `par` setup block is
duplicated verbatim in both branches; (f) fewer than 3 populations surfaces
raw ape's "cannot build an NJ tree with less than 3 observations" rather than
a dartR-styled FS5 guard; (g) `@importFrom stringr str_pad` and
`@importFrom graphics hist` are unused in this file.
No failure scenario rises above nuisance level; bundled as notes.

Out-of-scope note: `gl.dist.pop`'s unscaled Euclidean
(`sqrt(sum(sq))` over pairwise-complete loci, no rescaling by locus count)
and `stats::dist`'s proportional NA rescale are different estimators under
missing data — worth attention in `gl.dist.pop`'s own review, not here.

## Proposed changes

1. Compute frequencies with `mean(e, na.rm = TRUE) / 2` (F1).
   **Consequence: numerical output changes for any dataset with missing
   data under the default distance.**
2. Replace `popNames(x) <- indNames(x)` with `pop(x) <- indNames(x)` so
   `by.pop = FALSE` works (F2).
3. Accept `"upgma"` (retaining `"ugpma"` as an alias), correct the spelling
   in the roxygen, and gate the method-fallback warning at `verbose >= 1`
   (F3, part of F9).
4. Add `plot.display = TRUE` (forced FALSE at `verbose == 0`), validate
   `type` up front with `match.arg`, and decouple plotting from the return
   so a plot failure cannot lose the tree (F4).
5. Root with `resolve.root = TRUE` when `outgroup` is supplied (F5).
   **Consequence: returned tree topology gains a resolved root node for
   outgroup callers.**
6. Correct the `as.pop` error message to reference `ind.metrics` (F6).
7. Divide frequencies by ploidy (or document the halved SilicoDArT edge
   scale) (F7).
8. Coerce `dist.matrix` via `as.dist()` before use (F8).
9. Documentation pass per F9 (a-g) + `devtools::document()` (F9, F10g).
10. Optional tidy per F10 a-f.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run (DEP1 n/a: ape and
  stringr are in Imports, `import(ape)` in NAMESPACE covers `as.phylo`)
- Spec: independent recomputation of the distance + `ape::nj`/`hclust` on
  testset.gl and testset.gs — run (RF 0, edge sums equal)
- NA handling: quantified on testset.gl (23.8% NA cells; RF 32 vs na.rm) — run
- Dispatch SNP vs SilicoDArT — run (both admitted; F7)
- Edge cases: 2 pops (ape error), 3 pops (OK), `as.pop`, bad `as.pop`,
  outgroup, bad outgroup, `treefile`, invalid `type`, user `dist.matrix`
  (dist and matrix), `method` all spellings, `by.pop = FALSE`,
  `verbose = 0` silence via `capture.output` — run
- Population of n=1: present in testset.gl (2 such pops) and exercised
  implicitly by the baseline; not isolated as its own fixture
- Populations with identical allele frequencies (zero distances): SKIPPED —
  no natural fixture; `ape::nj` on tied distances is deterministic but the
  case was not synthesized
- `labelsize` effect on the rendered plot: code-inspection only (passed as
  `cex` to `plot.phylo`); pixel-level verification impracticable headless
- VRB5 plot side at `verbose = 0`: code-inspection (plot call is
  unconditional, no display gate exists to test empirically)
- FBM path (DAT6): SKIPPED — no FBM fixture; `as.matrix(x)` densifies
  regardless (noted, not counted as a finding for a 30x255 test object)

## Approval (Phase B)

All decisions recorded 2026-09-06 via the formal approval boxes, each with
its stated consequence acknowledged.

| Change | Decision | By | Note |
|---|---|---|---|
| 1 (F1) | approved | Arthur Georges | Consequence acknowledged: the default tree changes on any dataset with missing data (RF 32 on testset.gl); na.rm = TRUE matches gl.dist.pop's policy |
| 2 (F2) | approved | Arthur Georges | Repair via `pop(x) <- indNames(x)`; documented by.pop = FALSE mode becomes functional |
| 3 (F3) | approved | Arthur Georges | Accept "UPGMA"/"upgma"; KEEP "ugpma" as a silent back-compat alias; unknown-method fallback becomes a gated message |
| 4 (F4) | approved | Arthur Georges | plot.display argument + failure-proofing; plot decoupled from the result |
| 5 (F5) | approved | Arthur Georges | resolve.root = TRUE; outgroup trees rooted in ape's sense |
| 6 (F6) | approved | Arthur Georges | |
| 7 (F7) | approved | Arthur Georges | Divide by ploidy; SilicoDArT branch lengths change scale (doubled) |
| 8 (F8) | approved | Arthur Georges | |
| 9 (F9) | approved | Arthur Georges | Includes F10g (unused @importFrom removal) per the change list |
| 10 (F10) | no action | Arthur Georges | INFO items a-f left as-is |

## Outcome (Phase C)

Applied 2026-09-06 on branch `review-gl.tree.nj` (base ddaed27,
upstream/dev). Findings F1-F9 implemented as scoped; F10 no action
(F10g folded into the change-9 documentation pass). `devtools::document()`
run; only `man/gl.tree.nj.Rd` and the two gl.tree.nj-only `importFrom`
removals committed from `NAMESPACE`. NEWS entry added with the F1
topology change stated prominently.

Caller grep: no callers of `gl.tree.nj` outside its own file, Rd and test
in any of the 8 dartRverse clones under `D:\workspace\R`
(dartR.base, dartR.captive, dartR.popgen, dartR.sim, dartR.spatial,
dartR.sexlinked, dartRstartup, dartRverse); dartr2shiny has no local
clone; archive-only hits ignored. All clear.

Verification (characterization test updated pin-by-pin, every diff maps
to an approved finding; 40 assertions, 0 failures):

- F1: default tree matches an independent na.rm = TRUE rebuild — RF 0,
  edge sums equal (18.23083 on testset.gl, dartR.data 1.2.5); RF 32 from
  the old no-na.rm tree, as predicted.
- F2: `by.pop = FALSE` returns a phylo with 250 tips labelled
  `indNames(testset.gl)`.
- F3: `method = "UPGMA"` and `"ugpma"` produce identical UPGMA trees
  (RF 0 to `as.phylo(hclust(d, "average"))`); unknown method falls back
  to nj with the warning silent at verbose 0, printed at verbose 1.
- F4: with `ape::plot.phylo` mocked to fail, the tree is still returned
  (gated warning at verbose >= 1); `type` validated up front; at
  verbose = 0 the plot is gated off entirely.
- F5: outgroup rooting yields `ape::is.rooted() == TRUE`.
- VRB5: fully silent at verbose = 0 (zero captured lines) including the
  former ungated fallback warning.
- F7: testset.gs tree matches the ploidy-1 reference — RF 0, edge sums
  equal (47.33798).
- F8: plain-matrix `dist.matrix` on the upgma path now works (RF 0 to
  the dist-object result).
- End-to-end run at verbose = 3 with outgroup + treefile: clean.

PR: green-striped-gecko/dartR.base #370
(https://github.com/green-striped-gecko/dartR.base/pull/370).

```json
{
  "function": "gl.tree.nj",
  "package": "dartR.base",
  "family": "analysis",
  "skill_version": "1.0.0",
  "commit": "ddaed27",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "PLT3", "status": "applied", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "applied", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DAT7", "status": "applied", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "medium", "rule": "DAT5", "status": "applied", "change": 8},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 9},
    {"id": "F10", "severity": "INFO", "confidence": "high", "rule": "STY1", "status": "no-action", "change": 10}
  ],
  "coverage_skipped": [
    "DAT6: no FBM fixture",
    "zero-distance population pair: no fixture synthesized",
    "labelsize/VRB5 plot rendering: code-inspection only (headless)"
  ],
  "status": "pr-open",
  "pr": 370
}
```
