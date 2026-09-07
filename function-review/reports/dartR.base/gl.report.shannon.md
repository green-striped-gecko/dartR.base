# Review: gl.report.shannon (dartR.base)

- Family mode: report (with analysis-grade numerical verification — the function computes diversity statistics)
- Date: 2026-09-07
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v2.0.0
- Package commit: ed99203 (working tree; `git diff upstream/dev -- R/gl.report.shannon.r` vs ddaed27 is empty — the reviewed file is byte-identical to upstream/dev, so `load_all()` exercised the reviewed code)
- Datasets: testset.gl (250 x 255), testset.gs (218 x 255, rejection path), plus a constructed 3 x 4 fixture with an all-reference individual
- Baseline: `tests/testthat/test-gl.report.shannon.R` (new file; 24 assertions, all passing pre-review)

## What the function actually computes

Not Shannon entropy per locus or per population. For **each individual**, it takes the row of the dosage matrix (`as.matrix(x)`, 0/1/2), drops NAs and zeros, normalises the surviving dosages to proportions `p`, and computes **Hill numbers (diversity of order q)** for `q = 0 .. order-1`:

- `q = 0`: count of loci with non-zero dosage (richness);
- `q = 1`: `exp(-sum(p log p))` — the exponential of Shannon entropy in nats (Hill numbers are log-base free);
- `q >= 2`: `(sum(p^q))^(1/(1-q))` (q = 2 is inverse Simpson).

This is the individual-level diversity profile of Ma, Li & Zhang (2020); the internal `d.chao()` is their published OSI code. There is no per-locus, per-population, or hierarchical aggregation: population assignments are ignored entirely (the pop-check block is commented out), and each row of the returned data frame is one individual. Monomorphic loci contribute nothing beyond dosage counts; `0 * log(0)` never arises because zeros are dropped before the entropy sum. NA policy: missing genotypes are dropped per individual before normalisation.

## Verdicts

**Standards: Needs work** — the skeleton is recognisably house style (flag start, datatype check, plot save via `utils.plot.save`), but the preamble deviates on three confirmed rules: hardcoded `verbose = 2` default, a malformed and pointless dependency guard, and no verbose-0 gating of output.

**Spec: Needs work** — the `alpha` path is numerically exact against independent hand computation, and the report contract (input untouched, plot-decoupled, invisible return) holds; but two of the three documented `level` options return degenerate values, and the roxygen block is placeholder text that describes none of this.

What works well: the core `alpha` computation reproduces the cited OSI algorithm exactly, and the read-only contract is honoured on every path tested.

## Findings

**F1 [HIGH, confidence: high] — `level = "beta"` and `"gamma"` are degenerate (DOC5, proposed rule)**
`R/gl.report.shannon.r:90-125,145` — `d.chao()` implements alpha/beta/gamma for an abundance *matrix* (individuals x loci within a group), but the function only ever passes it a single individual's *vector*. For a vector, `cA <- A` and `N <- 1`, so gamma collapses to alpha and beta = gamma/alpha = 1.
Failure scenario (verified on testset.gl): `level = "gamma"` returns values byte-identical to `level = "alpha"` for all 250 individuals and all orders; `level = "beta"` returns exactly 1.0 everywhere. The docs say `level` "also accept 'beta', 'gamma'" — a user partitioning diversity gets meaningless numbers with no warning.
Proposed change: either restrict `level` to `"alpha"` (validate and error on the rest, docs corrected), or implement the Ma et al. matrix form — group individuals (e.g. by population) into abundance matrices so beta/gamma carry their published meaning. The first is the bounded fix; the second is a design decision for the custodian.

**F2 [HIGH, confidence: high] — global verbosity silently ignored (FS2, DOC2)**
`R/gl.report.shannon.r:43,45` — the signature defaults `verbose = 2` instead of `NULL`. `gl.check.verbosity(2)` returns 2 without consulting `options(dartR_verbose)`.
Failure scenario (verified): with `gl.set.verbosity(0)` in effect and no `verbose` argument, the function still prints 6 lines ("Starting gl.report.shannon", processing messages, "Completed:"). The `@param verbose` text "default 2, unless specified using gl.set.verbosity" claims the opposite of what happens.
Proposed change: `verbose = NULL` in the signature and the DOC2 canonical param text.

**F3 [MEDIUM, confidence: high] — verbose = 0 is not silent (VRB5)**
`R/gl.report.shannon.r:154,192-193` — two leaks. (a) `reshape2::melt(div_mat)` with no `id.vars` emits the message "Using ID as id variables" at every verbosity (verified: 1 message line at `verbose = 0`). (b) There is no `if (verbose == 0) plot.display <- FALSE` gate, so the default `plot.display = TRUE` prints the plot at `verbose = 0`.
Failure scenario: a user running silent (scripts, loops) gets console noise and a plot window per call.
Proposed change: `reshape2::melt(div_mat, id.vars = "ID")` and add the standard verbose-0 plot gate.

**F4 [MEDIUM, confidence: high] — no validation of `level` or `order` (FS5)**
`R/gl.report.shannon.r:41-42` — neither parameter is checked.
Failure scenario (verified): `level = "banana"` fails with `object 'D.value' not found` (no `d.chao` branch matches, from deep inside the loop); `order = 0` fails with `subscript out of bounds` (the inner loop `0:list_order` becomes `0:-1` and indexes column 0 of a zero-column matrix). Both errors are opaque and fire after the preamble has already printed.
Proposed change: fail fast with `stop(error(...))` — `level` in `c("alpha","beta","gamma")` (or just `"alpha"` per F1), `order` a positive integer.

**F5 [MEDIUM, confidence: high] — malformed, dead dependency guard (DEP1)**
`R/gl.report.shannon.r:47-55` — `requireNamespace(c("dplyr", "tidyr"), quietly = TRUE)` only ever checks the first element (verified: `requireNamespace(c("dplyr", "nonexistentpkg999"))` returns TRUE), the failure path uses `cat(error(...)); return(-1)` instead of `stop(error(...))`, and neither dplyr nor tidyr is used anywhere in the function — the real dependency is `reshape2::melt`, which is in Imports and needs no guard. The `@description` compounds this by naming a third package, adegenet, as the requirement.
Failure scenario: latent — if the guard ever fired, a report function would return `-1` to an unsuspecting caller instead of stopping. Today it is dead code that misleads the reader.
Proposed change: delete the guard block.

**F6 [MEDIUM, confidence: high] — all-reference / all-NA individuals return a 0/1/Inf row with no warning (VRB4, proposed rule)**
`R/gl.report.shannon.r:140-147` — an individual whose retained dosages are empty (all homozygous-reference, or all missing) feeds `d.chao(numeric(0))`, yielding q0 = 0, q1 = 1, q2+ = Inf (verified on the 3 x 4 fixture; the all-NA case is identical).
Failure scenario: `Inf` propagates into the returned data frame and the plot silently; a downstream `colMeans` or sort is poisoned without any indication which individuals were uninformative.
Proposed change: warn at `verbose >= 1` naming the affected individuals, and return `NA` (or document the convention) for their rows.

**F7 [MEDIUM, confidence: high] — roxygen block is placeholder text (DOC1, DOC5)**
`R/gl.report.shannon.r:1-34` — the title spans two lines, ends with a citation and full stop; there is no `@name`, `@title`, `@family`; the `@description` slot is occupied by the false sentence "This function needs package adegenet, please install it."; `@details` is literally "details" plus a one-item list. Nothing documents what is computed (Hill numbers, orders q0..q(order-1)), the dosage-as-abundance semantics, the NA/zero policy, or that `q1` is the exponential of Shannon entropy — the only connection to the function's name.
Failure scenario: a user cannot determine from `?gl.report.shannon` what any returned column means.
Proposed change: rewrite the header to the DOC1 tag order with a real description and details section stating the formulas and conventions above.

**F8 [LOW, confidence: high] — example wrapped in `\dontrun` (DOC3)**
`R/gl.report.shannon.r:31-34` — the only example never runs under any check, so nothing exercises the function on CRAN.
Proposed change: unwrap (it runs in ~1 s on `possums.gl[1:30,]`) or use `\donttest`.

**F9 [LOW, confidence: high] — `@author` lacks a Custodian line (DOC7, proposed rule)**
`R/gl.report.shannon.r:23-24` — "Ching Ching Lau (Post to ...)" names an author only.
Proposed change: `Author(s): Ching Ching Lau. Custodian: Ching Ching Lau -- Post to \url{...}`.

**F10 [LOW, confidence: high] — output phase after FLAG SCRIPT END; plot bundle incomplete (FS9, PLT1)**
`R/gl.report.shannon.r:186-193` — the plot prints after "Completed:" (phase order inverted relative to the family pattern), and the plot bundle omits `plot.colors` (the bar fill ignores house palettes).
Proposed change: move the print block before FLAG SCRIPT END; add `plot.colors` if the custodian wants palette control, otherwise note the omission as accepted.

### Notes (not findings)

- The function name promises Shannon but delivers full Hill-number profiles; only `q1` is Shannon-derived. Covered by the F7 rewrite rather than a rename (API3 territory).
- Dead code: the commented-out population-assignment block at lines 70-84 should go when the file is next touched.
- DAT6 (proposed): `as.matrix(x)` densifies the whole object; the example even suggests `gl.gen2fbm` first, after which the densification defeats the purpose. Row-wise access would serve. Noted, not counted.
- Other-function note: `tests/testthat/test-gl.report.basics.R` pins testset.gl at 755 loci x 274 individuals "dartR.data 1.2.5", but the installed restored-classics testset.gl is 250 x 255 — that baseline's pins belong to a different data version. One line, not this review's scope.

## Proposed changes

1. Validate `level`, restricting it to `"alpha"` with an informative error for `"beta"`/`"gamma"` until a matrix-form implementation exists; correct the `level` docs (F1, F4-level part). **Consequence: `level = "beta"`/`"gamma"` calls that currently return degenerate values will error instead.**
2. Change the signature default to `verbose = NULL` and adopt the DOC2 param text (F2). **Consequence: users relying on the global `gl.set.verbosity` setting will see it honoured; effective default remains 2.**
3. Silence the melt message (`id.vars = "ID"`) and gate `plot.display` off at `verbose = 0` (F3).
4. Validate `order` as a positive integer with `stop(error(...))` (F4, order part).
5. Delete the dead dplyr/tidyr guard block (F5).
6. Warn at `verbose >= 1` for individuals with no non-zero, non-NA dosages and return `NA` for their rows (F6). **Consequence: numerical output changes for all-reference/all-NA individuals (0/1/Inf becomes NA).**
7. Rewrite the roxygen header: DOC1 tag order, real title/description/details with formulas, log-base note, NA/zero policy, `\dontrun` unwrapped, Custodian line (F7, F8, F9). Docs only.
8. Move the plot-print phase before FLAG SCRIPT END (F10). Docs/structure only, no behaviour change.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: independent hand computation of Hill numbers q0..q4 on testset.gl (250 individuals) — run, exact agreement (max abs diff 0) for `alpha`
- Degenerate-level check: gamma vs alpha byte-compare, beta == 1 — run, confirmed
- Report contract: input serialize-byte-compare before/after — run, untouched; no history append (FS8) — confirmed by source, no append on any path; invisible return — run, confirmed
- PLT3 plot decoupling: results identical `plot.display` on/off — run, identical
- Population-independence: shuffled pop labels — run, identical results
- VRB5: `capture.output` at verbose 0, stdout and message streams separately — run (stdout clean, message leak found); plot gate — confirmed absent in source
- Edge cases: all-reference individual, all-NA individual, invalid `level`, `order` 0 and 1, single individual, SilicoDArT rejection — run
- Non-alphabetical pop levels lens: NOT APPLICABLE — populations are ignored by the computation (verified)
- FBM path (DAT6): SKIPPED — no FBM fixture in the test environment
- Google Group / GitHub issue search: SKIPPED — no web access in this session
- Small-sample bias correction check: NOT APPLICABLE — no corrected estimator is claimed or implemented; raw plug-in Hill numbers only

## Approval

Approved by Arthur Georges via the formal approval boxes, 2026-09-07. All 10 findings approved; none rejected or deferred.

- F1 — approved as the **implement** variant of the proposed change: build the real population-level alpha/beta/gamma partition (population abundance matrices feeding the existing Hill-number machinery), not the restrict-to-alpha bounded fix. Consequence acknowledged: `level = "beta"`/`"gamma"` outputs change from meaningless constants (beta identically 1, gamma identical to alpha) to real diversities, one row per population; the partition identity gamma = alpha x beta to be verified numerically at every order on the fixtures.
- F2 — approved: `verbose = NULL` default + `gl.check.verbosity()`; `gl.set.verbosity()` honoured; the 6-line leak at global verbosity 0 stops.
- F3, F4, F5, F6, F7 (all MEDIUM) — approved, including F6's consequence (0/1/Inf rows for degenerate individuals become NA with a gated warning) and F7's full roxygen rewrite to describe what the function actually computes.
- F8, F9, F10 (all LOW) — approved.

## Outcome

All 8 approved changes applied on branch `review-gl.report.shannon` (from `upstream/dev` ddaed27).

- F1 implemented: for `level = "beta"`/`"gamma"`, each population's individuals form an abundance matrix (individuals x loci, NA as zero, degenerate individuals excluded) passed to the unchanged `d.chao()` matrix forms; the return is one row per population (`pop`, q0 ... q(order-1)). Objects without population assignments are assigned to a single population 'pop1' with a gated warning. `level = "alpha"` output is byte-identical to the pre-review implementation (serialize-compare on testset.gl and fixtures at orders 3, 5, 7).
- Partition identity: max |gamma - alpha x beta| = 1.4e-14 across 30 testset.gl populations at orders q0..q6 (alpha recomputed independently of the function); on the hand-computable 2-pop toy, gamma matched pooled-frequency gamma exactly (diff 0) at every order and the q0/q1/q2 values matched closed-form hand calculations.
- Verbosity: global `gl.set.verbosity(0)` with no `verbose` argument → 0 stdout lines, 0 messages, 0 graphics devices opened (plot gate verified empirically without a null device); `verbose = 3` runs clean end to end for alpha and beta.
- Validation: `level = "banana"` → "Fatal Error: level must be one of 'alpha', 'beta' or 'gamma'"; `order = 0` and `order = 2.5` → "Fatal Error: order must be a positive whole number ...". Degenerate-individual warning fires at `verbose >= 1` naming the individuals, silent at 0.
- Characterization test rerun: all assertions pass; every flip from the pinned baseline is confined to an approved finding and tagged `[approved Fn]` in the test file; partition-identity assertions added.
- Report contract re-verified: input serialize-byte-identical before/after; plot on/off results identical for alpha and beta.
- Caller grep across all 8 dartRverse clones: no callers of `gl.report.shannon` outside its own file; no dartr2shiny clone present locally. All clear.
- NEWS.md entry added, leading with the beta/gamma consequence.
- `devtools::document()` run; `man/gl.report.shannon.Rd` regenerated (NAMESPACE unchanged).

PR: green-striped-gecko/dartR.base #(see JSON `pr` field), branch `review-gl.report.shannon` → `dev`.

```json
{
  "function": "gl.report.shannon",
  "package": "dartR.base",
  "family": "report",
  "skill_version": "2.0.0",
  "commit": "ed99203",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "FS2", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "applied", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DEP1", "status": "applied", "change": 5},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "VRB4", "status": "applied", "change": 6},
    {"id": "F7", "severity": "MEDIUM", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "DOC3", "status": "applied", "change": 7},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "applied", "change": 7},
    {"id": "F10", "severity": "LOW", "confidence": "high", "rule": "FS9", "status": "applied", "change": 8}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "issue-tracker search: no web access"],
  "status": "pr-open",
  "pr": null
}
```
