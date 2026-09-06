# Review: gl.dist.phylo (dartR.base)

- Family mode: analysis
- Date: 2026-09-06
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v2.0.0
- Package commit: ddaed27 (upstream/dev; local HEAD ed99203 — `git diff upstream/dev -- R/gl.dist.phylo.r` is empty, so `load_all()` exercised the reviewed code exactly)
- Datasets: platypus.gl (all fixtures; dartR.data 1.2.5), testset.gs (SilicoDArT rejection path). testset.gl checked for prerequisites (`TrimmedSequence`, `SnpPosition` present) but not run through the function.
- Baseline: `tests/testthat/test-gl.dist.phylo.R` (new file; 40 assertions, all pass at the reviewed state; defect pins tagged `[pins defect]`)

What the function is: a thin pipeline, not a distance engine of its own. It (1) optionally filters on tag length (`gl.filter.taglength`), (2) delegates sequence assembly to `gl2fasta(method = 1)` — one concatenated pseudo-sequence per individual, built by substituting the SNP state into `TrimmedSequence` at the 0-based `loc.metrics$SnpPosition`, heterozygotes as IUPAC ambiguity codes, NA genotypes as an all-N tag — then (3) reads the FASTA back with `ape::read.dna` and computes all distances with `ape::dist.dna(model = subst.model, pairwise.deletion = pairwise.missing)`, and (4) for `by.pop = TRUE` (the default) averages the individual-pair distances between populations (mean over all between-population pairs; within-population diagonal never computed). The substitution models are ape's, applied to the assembled tag alignment. The pipeline uses `loc.metrics$SnpPosition` throughout (`gl2fasta` was migrated off `@position`; only commented-out lines remain) — consistent with the PR #330 invariant, no stale `@position` reliance.

## Verdict

**Standards: Needs work** — the FS preamble conforms and messaging is gated, but the SilicoDArT branch is dead code, the working directory is stranded on any mid-pipeline failure, and the roxygen block deviates from the house canon in several places.

**Spec: Needs work** — the numerics are exactly right end-to-end (independent verification below), but the numbers do not mean what a reader of the docs would assume: heterozygous sites are silently excluded from every distance, the promised gamma corrections cannot be requested, and four edge paths fail opaquely.

What works well: the sequence assembly, both substitution-model distances checked, the population averaging, and the label handling are all exactly correct, including under non-alphabetical population factor levels.

## Independent verification (spec axis, analysis mode)

Fixture: 12 platypus individuals (4 per population, 3 populations), 35 polymorphic loci, callrate 1 (no missing data), 141 heterozygous genotype cells, no secondaries, every SNP position strictly inside its tag.

- Sequence assembly: an independent reimplementation (character-vector substitution at `SnpPosition + 1`, IUPAC lookup, N-tag for NA) reproduces the `gl2fasta` output character-for-character for all 12 individuals (~2,400 bp each); FASTA labels are `indName_pop` in individual order.
- Distances, hand formulas from the assembled sequences (sites where both bases in {A,C,G,T}): K80 (`-log(1-2P-Q)/2 - log(1-2Q)/4`) max |Δ| 3.6e-17; F81 (`-b log(1 - p/b)`, `b = 1 - Σπ²`, π from whole-alignment counts) exact 0; raw exact 0 — against `gl.dist.phylo(by.pop = FALSE)`.
- Population averaging: hand-averaged individual distances by name match `by.pop = TRUE` output exactly (max |Δ| 0).
- Missing data: hand recomputation with pairwise and with complete deletion both match exactly on a 15-NA fixture.
- Label ordering (the `gl.dist.pop` F1 BLOCKER class): with population levels reversed to a non-alphabetical order, labels follow `popNames()` and every value agrees under name lookup with the forward fixture — `avg.dist` maps individuals by population name, not position. No defect.

## Findings

**F1 [HIGH, confidence: high] — heterozygous sites contribute nothing to any distance, undocumented (DOC5)**
`R/gl.dist.phylo.r:247-263` — assembly is hard-wired to `gl2fasta` method 1 (heterozygote → IUPAC ambiguity code), and `ape::dist.dna` treats ambiguity codes as missing data.
Failure scenario (reproduced): an individual heterozygous at 10 loci sits at raw distance exactly 0 from one homozygous-reference at every locus; with `pairwise.missing = FALSE` the effect is global — any site at which any individual is heterozygous is deleted for everyone, and on the review fixture every variable site vanished, collapsing all distances to 0. Distances are driven by homozygous differences only; between-population distances are biased downward wherever heterozygosity is appreciable. The `@details` describe twelve substitution models at length and never mention this.
Proposed change: state the ambiguity-as-missing behaviour in `@details` (including the `pairwise.missing = FALSE` amplification), and consider exposing `gl2fasta`'s method 2 (random allele assignment, which makes heterozygous sites count) as a het-handling argument, default unchanged.

**F2 [HIGH, confidence: high] — working directory stranded in tempdir() on any mid-pipeline failure (FS7; same class as gl.tree.fitch F5)**
`R/gl.dist.phylo.r:239-254` — `hold <- getwd(); setwd(tempdir())` with the restore only on the success path; no `on.exit`.
Failure scenario (reproduced): `gl.dist.phylo(testset.gs)` fails inside `gl2fasta` and leaves the session in `C:/Users/.../Temp/Rtmp...`; every later relative path in the user's session lands in a directory that vanishes at session end.
Proposed change: drop the `setwd` dance entirely — `gl2fasta` already takes `outpath` and `ape::read.dna` accepts a full path (`file.path(tempdir(), "tmp.fas")`); the two `setwd` calls do no work. Failing that, `on.exit(setwd(hold), add = TRUE)`.

**F3 [MEDIUM, confidence: high] — the SilicoDArT rejection branch is dead code (DAT7, FS5; VRB2 for the idiom)**
`R/gl.dist.phylo.r:160-167` — `datatype == "silicodart"` can never be true: `utils.check.datatype` returns the string `"SilicoDArT"` (verified). The branch also uses the empty-condition idiom `cat(error(...)); stop()`.
Failure scenario (reproduced): a SilicoDArT object passes the guard, wastes a full `gl.filter.monomorphs` pass, then dies inside `gl2fasta` with the generic "found SilicoDArT expecting SNP" — the bespoke message directing users to `gl.dist.pop`/`gl.dist.ind` never prints, and the working directory is already stranded (F2). The gl.tree.fitch report's coverage note that this function "has a silico branch" is true only of the source text.
Proposed change: pass `accept = "SNP"` to `utils.check.datatype` and delete the dead branch (the standard rejection fires before any work); if the bespoke redirect message is wanted, fix the comparison to `"SilicoDArT"` and use `stop(error(...))`.

**F4 [MEDIUM, confidence: high] — a single population with by.pop = TRUE fails with "subscript out of bounds" (FS5)**
`R/gl.dist.phylo.r:207` — `for (i in 1:(n_pop - 1))` with `n_pop = 1` runs `i = 1, j = 2` against a 1×1 matrix.
Failure scenario (reproduced): `gl.dist.phylo(one_pop_object)` errors with the bare subscript message, nothing pointing at the population count.
Proposed change: `stop(error(...))` when `by.pop = TRUE` and `nPop(x) < 2`, naming the requirement.

**F5 [MEDIUM, confidence: high] — all-monomorphic input crashes inside the advisory monomorph probe (FS5)**
`R/gl.dist.phylo.r:187` — `tmp <- gl.filter.monomorphs(x, verbose = 0)` exists only to warn, but when every locus is monomorphic the filter's internal subsetting errors first.
Failure scenario (reproduced): the call dies with "Subsetting resulted in zero loci." from `gl.drop.loc` internals, at every verbosity, before any message from this function.
Proposed change: make the probe safe (count monomorphs without filtering, or wrap in `tryCatch`) and issue an informative fatal error when no polymorphic loci remain.

**F6 [MEDIUM, confidence: medium] — documented gamma corrections cannot be requested (DOC5)**
`R/gl.dist.phylo.r:41-122` — the `@details` for JC69, K80 and TN93 state a gamma correction is available and its "parameter must be given by the user", and the raw entry discusses the `variance` option; the signature (`x, subst.model, min.tag.len, pairwise.missing, by.pop, verbose`, verified via `formals()`) has no `gamma`, `variance` or `base.freq` argument, and the `dist.dna` call passes none.
Failure scenario: a user following the docs has no way to supply the gamma parameter; the text describes ape's interface, not this function's.
Proposed change: either forward `gamma`/`variance`/`base.freq` (or `...`) to `ape::dist.dna`, or cut the unfulfillable claims from `@details`.

**F7 [LOW, confidence: high] — BH87 returns an asymmetric matrix, not the documented dist; by.pop averaging silently symmetrises it (DOC5)**
`R/gl.dist.phylo.r:263, 137` — ape returns a full asymmetric matrix for BH87; with `by.pop = FALSE` that matrix is returned as is (reproduced: class `matrix`), contradicting `@return` ("object of class dist"); with `by.pop = TRUE`, `avg.dist` averages only the `[inds_i, inds_j]` direction block and assigns it to both triangles.
Failure scenario: downstream code expecting a `dist` (e.g. `gl.tree.fitch`, which validates its `D`) breaks or silently mishandles the object.
Proposed change: document the BH87 exception (and the direction used in averaging), or exclude asymmetric models from `subst.model`.

**F8 [LOW, confidence: high] — by.pop = FALSE labels are undocumented `indName_pop` composites (DOC5)**
`R/gl2fasta.r:337` (consumed here via `read.dna`) — labels come from the FASTA headers, e.g. `T27_TENTERFIELD`, not `indNames(x)`.
Failure scenario: matching the matrix back to the genlight by `indNames()` fails; parsing the composite is ambiguous when individual or population names contain underscores.
Proposed change: document the label format in `@return`, or relabel with `indNames(x)` after `dist.dna` (relabelling is an output change for existing consumers — **Consequence: labels change for `by.pop = FALSE` callers.**).

**F9 [LOW, confidence: high] — roxygen deviations (DOC1, DOC2, DOC5, DOC6)**
`R/gl.dist.phylo.r:11-137` — `@param` order (`min.tag.len` before `subst.model`) disagrees with the signature; the `verbose` text deviates from the DOC2 canon; `pairwise.missing` is documented as the last `\item` inside the `subst.model` itemize list, as if it were a model; `@details` names the parameter `min.tag.length` (the argument is `min.tag.len`); the prerequisites (`TrimmedSequence`, `SnpPosition`, `loc.all`) are never stated — a VCF-derived genlight fails with `gl2fasta`'s message referencing a function the user did not call; curly quotes in `@details` ("saturation plots", "Kimura's 2-parameters distance") violate DOC6 (proposed rule).
Proposed change: documentation pass fixing the six items; `devtools::document()` in the same change (DOC4).

**F10 [LOW, confidence: high, proposed rule] — `@author` names a custodian only, no `Author(s):` line (DOC7)**
`R/gl.dist.phylo.r:123-124`.
Proposed change: `Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to \url{https://groups.google.com/d/forum/dartr}`.

**F11 [INFO, confidence: high] — minor style and gating items (STY1, STY2, VRB1)**
`R/gl.dist.phylo.r:187` — the advisory `gl.filter.monomorphs` pass runs unconditionally, including at `verbose = 0` where its warning can never print; `:171` — the adegenet-coercion warning is gated `verbose > 2`, so `verbose = 2` (the warning level) does not see it (reproduced: absent at 2, present at 3); `:217-221, 272-289` — three blocks of commented-out dead code (Phylip infile writing, matrix conversion). The visible `return(D)` versus the style guide's analysis-reporter `invisible()` convention is left as an open question for the custodian — the `@examples` call the function unassigned, suggesting the visible print is intended.

**F12 [INFO, confidence: high] — subst.model validation is delegated to ape; behaviour acceptable (FS5, no action)**
`R/gl.dist.phylo.r:263` — an invalid model errors inside `dist.dna` with `'model' must be one of: "RAW" "JC69" ...` (reproduced). Unlike the `gl.tree.nj` "ugpma" precedent there is no silent fallback, and the message lists the valid values, so no change is required; noted so the check is on record.

### Notes on other functions (scope rule — one line each, not findings here)

- `gl2fasta`: the sink guard `on.exit(if (sink.number() > 0) sink())` removes the *caller's* sink after its own explicit `sink()` has already closed the file — reproduced in isolation (a pre-existing sink is gone after a successful `gl2fasta` call) and via `capture.output(gl.dist.phylo(...))`, which loses every line printed after the `gl2fasta` return; the guard needs to remember its own sink depth.
- `gl2fasta` owns the heterozygote policy (`method` 1–4) that F1 turns on; any het-handling argument added here should pass through to it.
- `gl.tree.fitch` (PR #371): its forwarding of `subst.model`/`pairwise.missing`/`min.tag.len` into replicate `gl.dist.phylo` calls matches this function's actual signature — compatible; its report's "silico branch" coverage note should read "dead silico branch" (see F3).
- PR #330 interaction: none adverse — the pipeline reads `loc.metrics$SnpPosition` only; the commented `@position` lines in `gl2fasta` are already-migrated residue.

## Proposed changes

1. Document heterozygote handling (ambiguity codes → missing in ape; `pairwise.missing = FALSE` deletes any site with a heterozygote anywhere) and optionally add a het-handling argument exposing `gl2fasta` method 2, default unchanged (F1). **Consequence: none by default; numerical output changes only if a user selects the new option.**
2. Remove the `setwd` dance (full paths into `gl2fasta`/`read.dna`), or add `on.exit(setwd(hold), add = TRUE)` (F2).
3. Reject SilicoDArT via `accept = "SNP"` and delete the dead branch; if the redirect message is kept, fix the case comparison and use `stop(error(...))` (F3).
4. Fatal, informative error when `by.pop = TRUE` and `nPop(x) < 2` (F4).
5. Make the monomorph probe crash-proof and error informatively when no polymorphic loci remain (F5).
6. Forward `gamma`/`variance`/`base.freq` to `ape::dist.dna`, or delete the gamma/variance claims from `@details` (F6).
7. Document the BH87 asymmetric-matrix return and the averaging direction, or exclude asymmetric models (F7; exclusion is API2). **Consequence if excluded: `subst.model = "BH87"` calls error instead of returning a matrix.**
8. Document the `indName_pop` label format for `by.pop = FALSE`, or relabel with `indNames(x)` (F8). **Consequence if relabelled: labels change for `by.pop = FALSE` callers.**
9. Documentation pass: param order, DOC2 verbose text, `pairwise.missing` out of the model list, `min.tag.len` name, prerequisites stated, ASCII quotes, `Author(s):` line (F9, F10); `devtools::document()` in the same change.
10. Style tidy: skip the monomorph probe below `verbose 2`, gate the coercion warning at `verbose >= 2`, delete the commented dead code (F11).

## Coverage

- Standards walk (FS, DOC, VRB, DAT, DEP, PLT, STY): run. PLT n/a — the function has no plot, so plot decoupling and the plot bundle do not apply. `ape` is in Imports (DEP satisfied); `seqinr` is guarded inside `gl2fasta`.
- Independent verification (analysis mode): run — assembly character-identical to an independent reimplementation; K80/F81/raw hand formulas exact; population averaging exact; ape cross-check on the same FASTA agrees.
- Label ordering with non-alphabetical levels: run — no defect.
- Missing data: run — NA → all-N tag; pairwise vs complete deletion both hand-verified; complete-deletion collapse pinned.
- Heterozygote handling: run — exclusion reproduced and pinned.
- Edge cases: run — SilicoDArT, invalid model, 1 population, 2 populations, all-monomorphic input, missing `TrimmedSequence`, `min.tag.len` path, BH87, `verbose = 0` silence (zero lines), return visibility.
- Secondaries: assessed, not deep-run — the fixture had none; real data carry them (14/1000 loci in platypus.gl, 4/255 in testset.gl share a tag), and each such locus contributes a full duplicate copy of the tag with only its own SNP substituted — a documentation caveat (recommend noting `gl.filter.secondaries`), folded into change 9's prerequisites text rather than a separate finding.
- FBM path (DAT6): SKIPPED — no FBM fixture; `gl2fasta` densifies via `as.matrix` regardless.
- testset.gl end-to-end: SKIPPED — prerequisites verified present; platypus fixtures already cover the SNP path.
- gamma/variance behaviour: NOT EXERCISABLE — no parameter exists to pass (F6).
- dartR Google Group / GitHub issue sweep: SKIPPED — not searched in this pass.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges, 2026-09-06 | Document + het-load warning at `verbose >= 1` (fraction of het cells), stronger warning when `pairwise.missing = FALSE`; numerics unchanged; no het-handling argument |
| 2 | Approved | Arthur Georges, 2026-09-06 | Remove/guard the setwd dance; no output changes |
| 3 | Approved | Arthur Georges, 2026-09-06 | Early SilicoDArT rejection with the correct message |
| 4 | Approved | Arthur Georges, 2026-09-06 | Informative fatal error for single-population `by.pop = TRUE` |
| 5 | Approved | Arthur Georges, 2026-09-06 | Crash-proof monomorph probe; informative all-monomorphic error |
| 6 | Approved | Arthur Georges, 2026-09-06 | Expose the parameter or correct the docs to what `dist.dna` is asked |
| 7 | Approved | Arthur Georges, 2026-09-06 | Document (LOW batch) |
| 8 | Approved | Arthur Georges, 2026-09-06 | Document the `by.pop = FALSE` label format |
| 9 | Approved | Arthur Georges, 2026-09-06 | Docs pass incl. prerequisites and DOC7 author line |
| 10 | No action | Arthur Georges, 2026-09-06 | INFO items left as they stand |

## Outcome

Applied 2026-09-06 on branch `review-gl.dist.phylo` (base `upstream/dev`
ddaed27); PR #377.

- Changes 1-9 applied to `R/gl.dist.phylo.r`; change 10 not applied (INFO,
  no action per approval).
- Change 6 disposition: parameters exposed, not docs cut — new `gamma`
  (default FALSE) and `variance` (default FALSE) arguments forwarded to
  `ape::dist.dna()`, placed after `by.pop` and before `verbose` so all
  existing named and positional-through-`by.pop` calls are unaffected;
  `@details` reworded to point at the `gamma` argument; a gated
  (`verbose >= 2`) note reports that the variance attribute is not
  preserved by population averaging. `base.freq` not exposed (not claimed
  by the details text).
- Change 3 disposition: comparison fixed to `"SilicoDArT"` and the branch
  converted to `stop(error(...))` with the original redirect message, so
  the rejection fires before any processing (the `accept = "SNP"` variant
  would have replaced the bespoke message with the generic one).
- Change 5 disposition: monomorphs counted from the genotype matrix
  (unique-genotype check treating any heterozygote as polymorphic);
  `gl.filter.monomorphs()` no longer called; all-monomorphic input gets
  "no polymorphic loci" as a fatal error.

Verification (all rerun at the applied state):

- (a) All 12 pinned fixture matrices byte-identical to the pre-change
  snapshot at default arguments. Two explained exceptions, values exact
  in both: individual-level `dist` objects differ only in ape's recorded
  `call` attribute (now names the forwarded `variance`/`gamma`
  arguments — a metadata consequence of change 6); the BH87 matrix
  diagonal is nondeterministic in ape itself (uninitialized memory:
  7 of 20 repeated *pre-change-style* `dist.dna` calls differed,
  diagonal cells only, values from 1e-312 to 8e+298 — pre-existing ape
  defect, off-diagonal distances stable and exact).
- (b) Het-load warning fires at `verbose >= 1` with the exact count and
  fraction (fixture: "141 of 420 genotype calls (33.6%)"); silent at
  `verbose = 0`; the stronger global-deletion warning fires exactly when
  `pairwise.missing = FALSE`.
- (c) Forced mid-pipeline failure (missing `TrimmedSequence`, dies inside
  `gl2fasta`) leaves `getwd()` unchanged.
- (d) SilicoDArT rejected early with the redirect message; single-pop
  `by.pop = TRUE` and all-monomorphic inputs give the new informative
  fatal errors.
- (e) `gamma = 0.5` inflates all non-zero K80 distances (zero distances
  stay zero); `variance = TRUE` attaches ape's variance attribute at
  `by.pop = FALSE`; defaults reproduce the previous numerics exactly
  (see (a)).
- (f) Baseline characterization file: 48 assertions, 0 failures. Flips
  confined to approved findings, each tagged `# [approved Fn]` in the
  test file (F2 wd preserved, F3 message, F4 message, F5 message); new
  assertions added for F1 warnings and F6 pass-through.
- (g) `verbose = 3` end-to-end on `platypus.gl` clean (het warning, model
  log, Completed); `verbose = 0` silent — after filtering one known
  third-party leak: `gl.filter.overshoot` at ddaed27 prints "There were
  no loci with SNP falling outside the trimmed sequence" ungated (its
  own defect, already fixed on the integration branch by the pending
  overshoot review; commented in the test).
- Caller grep across all 8 clones: only `gl.tree.fitch` (dartR.base)
  calls `gl.dist.phylo` — all-named arguments on both dev and the PR
  #371 branch, compatible with the appended parameters. All clear.

```json
{
  "function": "gl.dist.phylo",
  "package": "dartR.base",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "ddaed27",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "FS7", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DAT7", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "applied", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "applied", "change": 5},
    {"id": "F6", "severity": "MEDIUM", "confidence": "medium", "rule": "DOC5", "status": "applied", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 8},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 9},
    {"id": "F10", "severity": "LOW", "confidence": "high", "rule": "DOC7 (proposed rule)", "status": "applied", "change": 9},
    {"id": "F11", "severity": "INFO", "confidence": "high", "rule": "STY2", "status": "no-action", "change": 10},
    {"id": "F12", "severity": "INFO", "confidence": "high", "rule": "FS5", "status": "no-action", "change": null}
  ],
  "coverage_skipped": [
    "DAT6/FBM: no fixture",
    "testset.gl end-to-end: platypus fixtures cover the SNP path",
    "gamma/variance: no parameter exists to exercise",
    "Google Group / GitHub issue sweep: not searched this pass"
  ],
  "status": "pr-open",
  "pr": 377
}
```
