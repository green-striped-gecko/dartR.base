# Review: gl.propShared (dartR.base)

## Provenance

- Model: Claude Fable 5 (claude-fable-5, Claude Code) via dartr-dev agent;
  Skill: dartr-function-review v2.0.0; Base: upstream/dev at ddaed27
  (`git diff upstream/dev -- R/gl.propShared.r` empty — the loaded code is
  the reviewed code); working branch integration-local at ed99203.
- Family mode: analysis (pairwise individual similarity kernel).
- Datasets: testset.gl (full and `[1:8, 1:60]`, `[1:20, 1:200]`,
  `[1:100, ]`), testset.gs `[1:6, 1:60]`, bandicoot.gl (the documented
  example), constructed missing-data fixtures (6x10 with an all-NA
  individual and a disjoint-call pair; 12x40 with 25 per cent missing),
  600x1000 random dosage fixture for benchmarking. dartR.data 1.2.5,
  R 4.4.2, Rcpp 1.0.13, Rtools44 present.
- Baseline: `tests/testthat/test-gl.propShared.R` (13 blocks, 35
  assertions, all pass at the reviewed state).
- Checks skipped: Google Group not searched (no browser session); the
  missing-Rcpp path (F2) reviewed statically — Rcpp cannot be removed from
  this machine's library without breaking the harness.

## Verdicts

**Standards: Rework** — no part of the house preamble is present (no
`verbose` argument at all, no `gl.check.verbosity`, no
`utils.flag.start`, no `utils.check.datatype`, no completion flag), the
Suggests guard returns the sentinel `-1` instead of erroring, and the
roxygen header names a man page that does not match the function and
loses its `@family` tag. Bringing the file to standard is a rebuild
against the skeleton, not a patch list.
**Spec: Needs work** — for SNP dosages the statistic is correct and the
returned object is a well-formed full symmetric similarity matrix (the
half-empty-matrix hypothesis is disproved below); the defects are the
missing SNP-only gate, an undocumented NaN policy, and a title that calls
one result both a similarity and a distance.

What works well: `res <- as.matrix(as.dist(res)); diag(res) <- 1` is a
correct and cheap symmetrisation of the half-filled C++ output, and the
statistic agrees with hand computation to machine precision.

## What the returned matrix actually contains

Verified on `testset.gl[1:8, 1:60]` and the missing-data fixtures:

| Property | Value |
|---|---|
| Class / shape | `matrix`, `nInd(x)` x `nInd(x)` |
| Symmetry | `isSymmetric()` TRUE — full matrix, both triangles populated |
| Upper triangle | populated (zero cells: 0 of 190 on the 20-individual case) |
| Diagonal | forced to 1 |
| Dimnames | `indNames(x)` on both margins |
| Orientation | similarity in [0, 1]; 1 = identical at all compared loci |
| Missing pairs | `NaN` where no locus is called in both individuals |

The C++ kernel does fill only `out(j, i)` for `j > i`, leaving the upper
triangle and diagonal at zero (verified directly against the raw kernel
output). Line 51, `res <- as.matrix(as.dist(res))`, then reads exactly
that lower triangle and mirrors it, so the half-empty intermediate never
reaches the caller. This is not a repeat of the `gl.dist.pop`
`type = 'matrix'` defect.

The statistic is `1 - mean(|g_i - g_j|) / 2` over loci called in both
individuals — pairwise-complete, verified equal to a hand-coded R
reference on `testset.gl[1:8, 1:60]` (exact) and confirmed as the
standard proportion of shared alleles for biallelic dosages: at a locus,
`2 - |g_i - g_j|` is the number of shared alleles (0, 1 or 2), so the
mean of `|g_i - g_j|/2` is the proportion of unshared alleles.

## Findings

**F1 [HIGH, confidence: high] — no datatype gate; SilicoDArT halves the dissimilarity (DAT7)**
`R/gl.propShared.r:24-36` — the function never calls
`utils.check.datatype()`. `gl.propShared(testset.gs)` runs to completion
and returns plausible-looking numbers.
Failure scenario: presence-absence data is scored 0/1 with ploidy 1, so
`|g_i - g_j|` is already the full mismatch; dividing by 2 halves it. On
`testset.gs[1:6, 1:60]` the returned similarity is exactly
`(1 + simple matching) / 2` (verified: `(1 - propShared) * 2` equals
`1 - simple matching` to machine precision), so the values are
compressed into [0.5, 1] and can never reach 0. A user comparing SNP and
SilicoDArT similarity matrices from the same function reads them on
different scales with no warning.
Proposed change: `datatype <- utils.check.datatype(x, accept = "SNP",
verbose = verbose)` in the CHECK DATATYPE slot. This is the recurring
DAT7 defect class already recorded against `gl2structure`, `gl2vcf` and
the `utils.recalc.*` siblings.
**Consequence: `gl.propShared(<SilicoDArT>)` errors instead of returning
a halved similarity.**

**F2 [HIGH, confidence: high] — Suggests guard returns -1 instead of erroring (DEP1, VRB2)**
`R/gl.propShared.r:27-34` — when Rcpp is absent the function prints with
`cat(error(...))` and executes `return(-1)`.
Failure scenario: the caller receives a length-1 numeric, not a matrix.
`dartR.spatial::gl.ibd` then evaluates `as.dist(1 - gl.propShared(x))`,
which is `as.dist(2)` — a valid, empty `dist` object — and the analysis
proceeds on nothing rather than stopping. `gl.spatial.autoCorr` behaves
the same way. The printed message is also ungated console output, which
`stop()` would not be.
Proposed change: replace the `cat()`/`return(-1)` pair with the DEP1
idiom, `stop(error("Package Rcpp needed for this function to work.
Please install it.\n"))`.
**Consequence: callers that currently receive -1 now get an error.**

**F3 [MEDIUM, confidence: high] — the C++ kernel is recompiled in every R session and is slower than the engine that already exists (DEP3 (proposed rule), STY2)**
`R/gl.propShared.r:38-49` — `Rcpp::cppFunction()` runs inside the
function body on every call. Rcpp's source cache is per-session, so the
first call in each session pays the full compile.
Measured on this machine (R 4.4.2, Rcpp 1.0.13, Rtools44):

| Case | Time |
|---|---|
| First call in a fresh session, 5x20 object | 3.55 s (compile only) |
| Second call, same session | 0.12 s |
| `testset.gl` (250 x 255), warm | 0.15 s |
| Same, plain R double loop | 0.75 s |
| 600 x 1000 fixture, first call | 5.81 s |
| 600 x 1000 fixture, warm | 2.75 s |
| Same, plain R double loop | 8.41 s |
| Same, `gl.dist.ind(method = "simple")` warm | 1.24 s |

The kernel is about 3x faster than plain R, but 2.2x *slower* than the
already-compiled `utils.dist.ind.snp` kernel that computes the same
quantity: the sugar expression `1 - mean(abs(na_omit(x(i,_) - x(j,_))) / 2)`
allocates three temporary vectors per pair, where `dist_mod` runs a
single scalar loop. A session that uses both functions pays two separate
~3 s compiles.
Failure scenario: `gl.spatial.autoCorr` calls `gl.propShared` once per
distance class per bootstrap; each new R session adds ~3.5 s before any
result, and on a machine without Rtools the call fails at compile time
with a toolchain error rather than a dartR message.
Proposed change: see Overlap analysis — delegate to
`gl.dist.ind(method = "simple")` and delete the kernel. If the kernel is
kept, hoist it out of the function body (compile once at first use and
cache in the package environment).

**F4 [MEDIUM, confidence: high] — no standard preamble and no verbosity control (FS2, FS3, FS4, FS9, VRB5)**
`R/gl.propShared.r:24-55` — the signature is `function(x)`. There is no
`verbose` argument, no `gl.check.verbosity()`, no `utils.flag.start()`,
no `utils.check.datatype()` and no completion flag.
Failure scenario: `gl.propShared(as.matrix(testset.gl[1:3, 1:10]))`
converts the input, compiles the kernel (about 3.5 s in a fresh
session), then fails at line 53 with
`unable to find an inherited method for function 'indNames' for
signature 'x = "matrix"'` — an S4 dispatch message that names neither
the function nor the problem. A user running the 250-individual case has
no way to ask for progress, and no way to suppress the F2 message.
The absence of `verbose` is not defensible for this function under
FS2/FS4: it is an exported user-facing `gl.*` that takes seconds to
minutes, has a dependency guard to report, and is called from other
packages' pipelines that pass `verbose = 0` to everything else. It is
silent by accident, not by contract.
Proposed change: add `verbose = NULL` to the signature and fit the
standard preamble (SET VERBOSITY, FLAG SCRIPT START, CHECK DATATYPE with
`accept = "SNP"`, FUNCTION SPECIFIC ERROR CHECKING, FLAG SCRIPT END).
FS8 does not apply — the return value is a matrix, not a genlight.
**Consequence: adds a `verbose` argument; at `verbose > 0` the function
prints start and end lines where it previously printed nothing.**

**F5 [MEDIUM, confidence: high] — roxygen header does not compile to the man page it appears to promise (DOC1, DOC4)**
`R/gl.propShared.r:1-4` — `@name gl.prop.shared` does not match the
function `gl.propShared`, and `@family distance` on line 4 is indented
under the wrapped `@title`, so roxygen parses it as title text, not a
tag. Verified in the generated `man/gl.prop.shared.Rd`:

- the file is named `gl.prop.shared.Rd`, not `gl.propShared.Rd`;
- `\title{}` ends with the literal string `@family distance`, which is
  what the CRAN PDF manual and the HTML index display;
- there is no `\seealso{}` block, and `gl.propShared` is absent from the
  `Other distance:` lists in `gl.dist.ind.Rd`, `gl.dist.pop.Rd`,
  `gl.fdsim.Rd` and `utils.dist.ind.snp.Rd`.

`?gl.propShared` does still work — roxygen writes `\alias{gl.propShared}`
alongside `\alias{gl.prop.shared}`.
Failure scenario: the function is invisible in the distance family
cross-references, and its manual entry is titled with a raw roxygen tag.
Proposed change: set `@name gl.propShared`, put `@family distance` at
column 0 as its own tag, then run `devtools::document()` in the same
change (DOC4) so `man/gl.prop.shared.Rd` is replaced by
`man/gl.propShared.Rd`.

**F6 [MEDIUM, confidence: high] — documentation contradicts itself and omits the return contract (DOC5 (proposed rule), DOC7)**
`R/gl.propShared.r:2-3,6-8,10,12,22` — several separate gaps:

- the title calls the result a "similarity (distance) matrix"; it is a
  similarity, and the two callers in dartR.spatial disagree about which
  it is (see Overlap analysis);
- `@description` calls it "an individual based distance matrix";
- `@return` is "A similarity matrix" with no statement of the statistic,
  the [0, 1] range, the forced diagonal, the dimnames, or the NaN policy
  (F7);
- `@param x` says "the genlight containing the SNP genotypes" but nothing
  enforces it (F1);
- there is no `@details` section and no `@references`;
- `@author` is `Bernd Gruber (Post to ...)` — no `Author(s):` label and
  no `Custodian:` label (DOC7; this file is the example cited in the
  rule text).

Failure scenario: a user reading "distance matrix" writes
`as.dist(gl.propShared(x))` and gets a similarity treated as a distance
— which is exactly what `gl.spatial.autoCorr` does.
Proposed change: docs-only. Title and description say "similarity";
`@details` states the formula and its equivalence to
`1 - gl.dist.ind(method = "simple")`; `@return` states shape, range,
diagonal, dimnames and NaN; `@author` gains the two labels.

**F7 [LOW, confidence: high] — undocumented NaN and a self-similarity of 1 for an all-NA individual (DOC5 (proposed rule))**
`R/gl.propShared.r:44,52` — for a pair with no locus called in both
individuals, `na_omit()` leaves an empty vector and `mean()` returns
`NaN`; `diag(res) <- 1` then reports self-similarity 1 for an individual
with no calls at all.
Failure scenario: verified on a 6x10 fixture — an all-NA individual
yields an entire `NaN` row with a 1 on the diagonal. Downstream,
`as.dist(1 - gl.propShared(x))` in `gl.ibd` carries the `NaN` into the
Mantel test, where `vegan::mantel` errors or returns `NA` depending on
the path, with no indication of which individual caused it.
Proposed change: document the policy in `@return`; optionally warn at
`verbose >= 1` (VRB4 — the result is affected) naming the individuals
with no overlapping calls.

**F8 [LOW, confidence: high] — duplicate individual names produce duplicate dimnames**
`R/gl.propShared.r:53` — `colnames(res) <- rownames(res) <- indNames(x)`
copies duplicates through.
Failure scenario: verified — with `indNames` `c("A", "A", "B", "C")`,
`res["A", "B"]` silently returns the row for the first `A` only. Nothing
is corrupted, but name-based subsetting of the result is not reliable.
Proposed change: none required beyond a note; if `gl.compliance.check()`
is invoked (DAT5) the names are made unique upstream.

**F9 [INFO, confidence: high] — FBM-backed objects are fully densified (DAT6 (proposed rule))**
`R/gl.propShared.r:36` — `as.matrix(x)`. Verified working on an FBM
object built with `gl.gen2fbm()` (result identical to the in-memory
path), so the documented example is correct, but the whole genotype
matrix is materialised. For the object sizes FBM exists to handle, that
defeats the purpose.
Proposed change: none in this review; the same densification is present
in `utils.dist.ind.snp`, so the fix belongs to whichever kernel survives
the overlap decision.

**F10 [INFO, confidence: high] — the C++ string cannot carry leading comments (STY1)**
`R/gl.propShared.r:39` — the source string begins directly with the
function signature, so the hazard recorded in the `gl2fasta` and PR #315
work is not triggered as written. Confirmed still live in Rcpp 1.0.13: a
leading comment containing parentheses or a pipe makes `cppFunction()`'s
signature scanner derive the wrapper name and argument types from the
comment text, and compilation fails
(`'dosage' was not declared in this scope`).
Proposed change: none now — a note for anyone applying F3 or adding
explanatory comments; put commentary inside the function body, as PR
#315 did.

## Overlap analysis

### Is `gl.propShared` the same statistic as `gl.dist.ind`?

Yes, exactly, for three of its methods.

`utils.dist.ind.snp` (the engine behind `gl.dist.ind`) computes, for
`method = "manhattan"` and `method = "czekanowski"`,
`sum(|g_i - g_j|) / (2 * L)` over pairwise-complete loci — the exact
complement of `gl.propShared`. After PR #315, `method = "simple"` is
`1 - mean(2 - |g_i - g_j|) / 2`, which is arithmetically the same thing
(recorded in that PR's own docs).

Verified cell for cell:

| Comparison | Result |
|---|---|
| `1 - gl.dist.ind(method = "simple")` vs `gl.propShared`, `testset.gl[1:20, 1:200]` | identical (max abs deviation 1.11e-16) |
| same, `method = "czekanowski"` | identical |
| same, `method = "manhattan"` | identical |
| same, 12x40 fixture with 25 per cent missing | identical |
| same, 600x1000 fixture | identical |
| `method = "absolute"` | different statistic (mean relative difference 0.0067) |
| `method = "euclidean"` | different statistic |

At **ddaed27 as released**, `method = "manhattan"`/`"czekanowski"` is
already the exact complement (the ddaed27 kernel's Manhattan branch is
`sumabs / (2 * L)`; PR #315 left it untouched and re-verified it against
manual recomputation). `method = "simple"` at ddaed27 is *not*
equivalent — its allele-sharing scoring was reference-allele asymmetric
until #315; reproducing the ddaed27 Simple branch in R and comparing on
the 12x40 fixture gives a maximum absolute deviation of 0.2368. So the
equivalence claim is: `gl.propShared == 1 - gl.dist.ind(method =
"manhattan")` at every state, and `== 1 - gl.dist.ind(method =
"simple")` from #315 onward.

Interface differences:

| | `gl.propShared` | `gl.dist.ind(method = "simple")` |
|---|---|---|
| Orientation | similarity | distance (`1 - x`) |
| Return class | `matrix` always | `dist` by default, `matrix` with `type = "matrix"` |
| Triangles | full symmetric | full symmetric (`type = "matrix"`) |
| Diagonal | forced to 1 | 0 |
| Dimnames | `indNames(x)` | `indNames(x)` |
| No-overlap pairs | `NaN` | `NA_real_` (both are missing values; `is.na()` TRUE either way) |
| Datatype gate | none (F1) | `accept = c("SNP", "SilicoDArT")` in the wrapper, `accept = "SNP"` in the engine |
| Verbosity | none (F4) | full standard preamble |
| Plot | none | histogram of distances, `plot.display`/`plot.file` |

A drop-in replacement body was tested:

```r
res <- 1 - as.matrix(gl.dist.ind(x, method = "simple", type = "matrix",
                                 plot.display = FALSE, verbose = 0))
diag(res) <- 1
colnames(res) <- rownames(res) <- indNames(x)
```

On the 6x10 all-NA/disjoint fixture this reproduces the current output
exactly (`all.equal` TRUE); the only difference is `NA` where the
current code puts `NaN`.

### Other pairwise-individual measures

- `utils.dist.ind.snp` — the engine above. Same quantity under three
  method names. Open PR #315.
- `gl.dist.ind` — the wrapper. Open review, PR pending; its F1 is a
  missing `sorensen` entry, unrelated to this overlap.
- `gl.dist.pop` — population-level, not individual-level. No overlap.
- `gl.dist.phylo` — model-based nucleotide distances from sequence tags.
  No overlap.
- `dartR.captive::gl.grm` — `rrBLUP::A.mat(as.matrix(x) - 1)`, the
  VanRaden genomic relationship matrix. Allele-frequency-centred and
  scaled, admits negative values, diagonal `1 + F`. A different quantity
  (expected additive relatedness, not observed allele sharing) and not
  redundant with `gl.propShared`.
- `dartR.spatial::gl.grm2` — the same GRM family, in the spatial package.
- `dartR.base::gl.fdsim`, `gl.kosman` (dartR.spatial) — different
  statistics.

A grep of `R/` across all eight live clones found no other implementation
of a shared-allele similarity from dosages.

### Caller inventory (eight live clones; dated backup copies excluded)

| Caller | Line | Call | Orientation verdict |
|---|---|---|---|
| `dartR.spatial::gl.ibd` | `R/gl.ibd.r:312` | `Dgen <- as.dist(1 - gl.propShared(x))` | **Correct.** Similarity converted to a distance; the forced diagonal of 1 maps to 0, which `as.dist` discards anyway. |
| `dartR.spatial::gl.spatial.autoCorr` | `R/gl.spatial.autoCorr.r:357` | `Dgen <- as.dist(gl.propShared(x_temp))` | **Sign-flipped.** See below. |
| `dartR.spatial::gl.genleastcost` | `R/gl.genleastcost.r:311` | `as.matrix(as.dist(propShared(xx)))` | Different function — `PopGenReport::propShared` on a genind, not `gl.propShared`. Same similarity-as-distance orientation question; out of scope here. |

No caller in dartR.base, dartR.captive, dartR.data, dartR.popgen,
dartR.sexlinked, dartR.sim, dartRstartup or dartRverse.

The `gl.spatial.autoCorr` case, in full: line 357 assigns the similarity
matrix to `Dgen`, and lines 365-374 then apply `Dgen <- 1 - Dgen` to a
branch whose condition includes `propShared`, with the comment "Reverse
genetic distance matrix so that correlated values indicated more similar
individuals". For the `Simple` and `Absolute` branches that reversal
turns a distance into a similarity, which is the stated intent. For the
`propShared` branch the input is *already* a similarity, so the same
line turns it into a distance — the opposite of every sibling branch.
The autocorrelation coefficient for `Dgen_method = "propShared"`
therefore carries the opposite sign to the other SNP methods. One-line
note only: different package, different job, and the same block also has
an operator-precedence problem (`&` binds tighter than `|`, so the
`dt == "SNP"` guard applies only to the `propShared` term).

### Recommendation

**Deprecate the kernel, keep the name: reimplement `gl.propShared` as a
thin wrapper over `gl.dist.ind(method = "simple")`.**

Reasoning: the two implementations compute the same statistic to machine
precision, the surviving one is 2.2x faster and already carries the
standard preamble, datatype gate and verbosity contract, and the wrapper
form reproduces the current return object exactly — so both dartR.spatial
callers keep working unchanged and no user-visible numbers move (on
dev, post-#315). Deleting the runtime `cppFunction()` call removes a
~3.5 s per-session compile and a hard Rtools dependency from a function
that gains nothing from having its own kernel.

Sequencing: this change must land after PR #315 is merged, because at
ddaed27 `method = "simple"` is not yet the equivalent statistic
(`method = "manhattan"` is, and could be used instead if the change is
wanted first). Keeping `gl.propShared` as an exported name — rather than
telling users to call `gl.dist.ind` — costs one small file and preserves
the two dartR.spatial call sites, the documented `bandicoot.gl` example,
and any user scripts.

The alternative, "keep but fix", means applying F1-F7 to a second
implementation of a statistic dartR.base already computes faster
elsewhere, and leaves two kernels to keep in step.

## Proposed changes

1. Replace the function body with a wrapper over
   `gl.dist.ind(x, method = "simple", type = "matrix",
   plot.display = FALSE, verbose = 0)`, preserving the current return
   shape (full symmetric matrix, diagonal 1, `indNames` dimnames), and
   delete the `Rcpp::cppFunction()` call and its dummy stub (F3, F9,
   F10). Depends on PR #315 being merged.
   **Consequence: no numerical change on dev post-#315 (verified exact);
   no-overlap pairs report `NA` instead of `NaN`.**
2. Add the standard preamble: `verbose = NULL` in the signature, SET
   VERBOSITY, FLAG SCRIPT START, CHECK DATATYPE with `accept = "SNP"`,
   FLAG SCRIPT END (F1, F4).
   **Consequence: adds a `verbose` argument; SilicoDArT input errors
   instead of returning a halved similarity; the function prints start
   and end lines at `verbose > 0`.**
3. Replace the `cat(error(...))` / `return(-1)` pair with the DEP1
   `stop(error(...))` idiom — needed only if change 1 is rejected and the
   Rcpp path is retained (F2).
   **Consequence: callers that currently receive -1 now get an error.**
4. Fix the roxygen header: `@name gl.propShared`, `@family distance` as
   its own unindented tag, and run `devtools::document()` in the same
   change so `man/gl.prop.shared.Rd` is replaced by
   `man/gl.propShared.Rd` (F5).
5. Rewrite the documentation prose: title and description say
   "similarity"; add `@details` with the formula and the `gl.dist.ind`
   equivalence; `@return` states shape, range, forced diagonal, dimnames
   and the missing-value policy; `@param x` notes SNP only; `@author`
   gains `Author(s):` and `Custodian:` labels (F6, F7, F8).

Notes carried elsewhere, not findings against this file:

- `dartR.spatial::gl.spatial.autoCorr:357` uses `gl.propShared` output as
  a distance and then inverts it, giving the `propShared` branch the
  opposite sign to its `Simple`/`Absolute` siblings; the same `if` has an
  `&`/`|` precedence problem. Different package, different job.
- `dartR.spatial::gl.ibd:312` is correct as written and stays correct
  under all five proposed changes.
- `dartR.spatial::gl.genleastcost:311` calls `PopGenReport::propShared`,
  not this function.
- `utils.dist.ind.snp` and `gl.dist.ind` defects are covered by their own
  reports (PR #315 applied; `gl.dist.ind` awaiting approval).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, STY — run
- PLT — not applicable (no plot bundle)
- FS8 history — not applicable (returns a matrix, not a genlight)
- Spec: behaviour vs roxygen on testset.gl, testset.gs, bandicoot.gl and
  constructed fixtures — run
- Returned-matrix structure (triangles, diagonal, dimnames, symmetry),
  including the raw pre-`as.dist` kernel output — run
- Statistic verified against hand computation — run
- NA policy (all-NA individual, disjoint-call pair) — run
- Edge cases: 1 individual, 2 individuals, duplicate `indNames`,
  non-genlight input — run
- Zero-individual case — SKIPPED: dartR's `[` refuses the subset before
  the function is reached
- Datatype gate on SilicoDArT — run
- FBM path (DAT6) — run, via `gl.gen2fbm(testset.gl[1:10, 1:100])`
- Rcpp leading-comment hazard — run, reproduced on Rcpp 1.0.13
- Timing: cold/warm, vs plain R, vs `gl.dist.ind` — run
- Overlap equivalence vs `gl.dist.ind` at the post-#315 state — run
  empirically; at ddaed27, `method = "simple"` compared against an R
  reproduction of the ddaed27 branch, `method = "manhattan"` from the
  #315 report's verification
- Missing-Rcpp path (F2) — SKIPPED: reviewed statically; Rcpp cannot be
  removed from this machine's library
- Google Group — SKIPPED: no browser session
- dartr2shiny — SKIPPED: not present in the local workspace

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | | | |
| 2 | | | |
| 3 | | | |
| 4 | | | |
| 5 | | | |

## Outcome

Phase C not started.

```json
{"function": "gl.propShared", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "rework", "verdict_spec": "needs_work",
 "findings": [
  {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT7", "loc": "R/gl.propShared.r:24-36", "status": "proposed", "change": 2},
  {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DEP1", "loc": "R/gl.propShared.r:27-34", "status": "proposed", "change": 3},
  {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DEP3", "proposed_rule": true, "loc": "R/gl.propShared.r:38-49", "status": "proposed", "change": 1},
  {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS2,FS3,FS4,FS9", "loc": "R/gl.propShared.r:24-55", "status": "proposed", "change": 2},
  {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DOC1,DOC4", "loc": "R/gl.propShared.r:1-4", "status": "proposed", "change": 4},
  {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5,DOC7", "proposed_rule": true, "loc": "R/gl.propShared.r:2-22", "status": "proposed", "change": 5},
  {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC5", "proposed_rule": true, "loc": "R/gl.propShared.r:44,52", "status": "proposed", "change": 5},
  {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "DAT5", "loc": "R/gl.propShared.r:53", "status": "proposed", "change": 5},
  {"id": "F9", "severity": "INFO", "confidence": "high", "rule": "DAT6", "proposed_rule": true, "loc": "R/gl.propShared.r:36", "status": "note"},
  {"id": "F10", "severity": "INFO", "confidence": "high", "rule": "STY1", "loc": "R/gl.propShared.r:39", "status": "note"}],
 "overlap": {"duplicates": "utils.dist.ind.snp via gl.dist.ind(method='simple'/'manhattan'/'czekanowski')",
  "equivalence": "exact, max abs deviation 1.11e-16, verified on testset.gl and missing-data fixtures",
  "equivalence_at_ddaed27": "manhattan/czekanowski exact; simple differs (max abs deviation 0.2368) until PR #315",
  "recommendation": "deprecate kernel, reimplement as wrapper over gl.dist.ind(method='simple') after PR #315",
  "callers": ["dartR.spatial::gl.ibd:312 (orientation correct)",
              "dartR.spatial::gl.spatial.autoCorr:357 (orientation sign-flipped)"],
  "timing_s": {"compile_per_session": 3.5, "propShared_600x1000_warm": 2.75,
               "gl.dist.ind_simple_600x1000_warm": 1.24, "plain_R_600x1000": 8.41}},
 "datasets": ["testset.gl", "testset.gs", "bandicoot.gl", "constructed"],
 "baseline_test": "tests/testthat/test-gl.propShared.R",
 "status": "awaiting-approval", "pr": null}
```


## Relocation addendum (2026-09-08)

gl.propShared moved to dartR.spatial at the custodian's direction (green-striped-gecko/dartR.spatial#34), where its only two callers live (gl.ibd, gl.spatial.autoCorr). The custodian directed a move with NO other changes, so none of this report's 10 findings were applied; they remain recorded here against the relocated function and are for dartR.spatial's own review to take up. The overlap finding stands: gl.propShared equals 1 - gl.dist.ind(method = 'simple') to 1.11e-16 once dartR.base PR #315 merges.
