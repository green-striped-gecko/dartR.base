# Review: gl.sim.genotypes (dartR.base)

- Family mode: analysis (simulation constructing a new genlight)
- Date: 2026-09-09
- Reviewer: Claude (claude-opus-5), dartr-function-review v2.2.0
- Package commit: `ddaed27` (upstream/dev; `git diff upstream/dev -- R/gl.sim.genotypes.r` is empty, so the reviewed state is upstream's)
- Working checkout: `D:\workspace\R\dartR.base` at `ed99203` (integration-local; uncommitted test edits present, none touching `R/`)
- Datasets: possums.gl, testset.gl, testset2.gl, platypus.gl, testset.gs (dartR.data 1.2.5)
- Baseline: `tests/testthat/test-gl.sim.genotypes.R` (new file; 22 blocks, all passing against the reviewed state)

## Verdict

**Standards: Needs work** — the FS preamble is present and correct, but the
returned object records the wrong history, the `n.ind` warning prints at
`verbose = 0`, `n.ind` is never validated, and the roxygen block is missing
`@details`, `@examples` and an `Author(s):` line.

**Spec: Rework** — the sampling engine itself is sound (verified: two
independent draws per genotype, Hardy-Weinberg recovered, no recycling), but
three of the function's own claims are untrue: the announced `n.ind` cap is
never applied, "a single population" is not enforced or checked so
multi-population objects are silently pooled, and the function aborts on four
of the five packaged datasets.

What works well: the core loop is the right algorithm. Each locus draws two
independent haplotype vectors of length `n.ind`, so the joint structure is
clean — this function does **not** carry the `gl.sim.crosses` recycling
defect.

## Findings

**F1 [BLOCKER, confidence: high] — announced `n.ind` cap is never applied (DOC5 proposed / numerical-correctness principle)**
`R/gl.sim.genotypes.r:44-47` — when `n.ind > n.loc` the function prints
`Setting n.ind to <n.loc>` and then simulates the full requested `n.ind`.
There is no assignment: the `if` block contains two `cat()` calls and nothing
else.
Failure scenario: `gl.sim.genotypes(possums.gl, n.ind = 300)` (200 loci)
prints "Setting n.ind to 200" and returns 300 individuals. A user who reads
the message believes the object was capped and reports a PCA on 300
individuals x 200 loci as if it were 200 x 200. Verified: requested 300,
`nInd()` 300.
This is the same class as `gl.sim.crosses`'s `n` — a stated behaviour with no
implementing line.
Proposed change: either apply the cap (`n.ind <- n.loc`) or reword the message
to a warning that states no cap was applied. Applying the cap changes returned
dimensions, so it crosses the escalation gate.

**F2 [BLOCKER, confidence: high] — an all-NA locus aborts the run (FS5, DAT5)**
`R/gl.sim.genotypes.r:55,57` — `gl.allele.freq()` returns `NA` for a locus
with no calls; `sample(prob = c(NA, NA))` then stops with
`NA in probability vector`, from inside `sample.int`, naming neither the locus
nor the function.
Failure scenario: `gl.sim.genotypes(testset.gl)` fails outright. Verified NA
frequency counts: testset.gl 3/255, testset2.gl 3/755, platypus.gl 6/1000,
testset.gs 3/255; possums.gl 0/200. The function runs on exactly one of the
five packaged SNP/SilicoDArT datasets without pre-filtering.
Proposed change: check `is.na(df$frequency)` after the frequency call and
`stop(error(...))` naming the offending loci and pointing at
`gl.filter.allna()`; or drop those loci with a gated warning.

**F3 [HIGH, confidence: high] — "a single population" is neither enforced nor checked (DOC5 proposed)**
`R/gl.sim.genotypes.r:40` — `@description` says the genotypes are drawn from
"the allele frequencies from that population" for "a single population", but
`gl.allele.freq(x, by = "loc")` pools every individual regardless of `pop()`.
A structured object is collapsed to one panmictic gene pool with no check and
no message.
Failure scenario: possums.gl has 10 populations. Source mean Ho = 0.347;
simulated mean Ho = 0.472, which is the pooled He (0.471). The Wahlund
deficit in the source is erased, so the simulated object is not a null model
for the data the user supplied. On a genuine single-population source the
match is faithful (source Ho 0.328, simulated 0.323), which is why the defect
is invisible unless the user tests a structured object.
Proposed change: `stop(error(...))` or a `verbose >= 1` warning when
`nPop(x) > 1`, and state in `@details` that frequencies are pooled across all
individuals in `x`.

**F4 [HIGH, confidence: high] — SilicoDArT input is accepted and returned as SNP (DAT7)**
`R/gl.sim.genotypes.r:36` — `utils.check.datatype()` is called without
`accept`, so its default admits SilicoDArT. The algorithm is diploid dosage
arithmetic (`v1 + v2` giving 0/1/2) with `ploidy = 2`, which has no meaning
for presence/absence data.
Failure scenario: `gl.sim.genotypes(gl.filter.allna(testset.gs))` returns an
object that `utils.check.datatype()` reports as `SNP`, with ploidy 2 and
genotype values 0, 1, 2 — from a source whose only values are 0 and 1. The
datatype has been silently changed.
Proposed change: `accept = "SNP"` on the `utils.check.datatype()` call.

**F5 [HIGH, confidence: high] — history omits the call that created the object (FS8, FS12 proposed)**
`R/gl.sim.genotypes.r:63-72` — the function never appends `match.call()`. The
returned object's `@other$history` holds exactly two entries, both internal:
`gl.recalc.metrics(x = x, verbose = 0)` and
`gl.compliance.check(x = gl, verbose = 0)`, left behind by the compliance
call at line 72.
Failure scenario: a user who saves a simulated dataset and later runs
`gl.print.history()` sees a compliance check and a metrics recalculation, with
no record that the data are simulated or of the `n.ind` used. The provenance
of a simulated object is exactly the case where history matters most.
Proposed change: reset `gl@other$history` to a single entry holding
`match.call()` after the compliance check, so the returned object's history is
the call that made it.

**F6 [MEDIUM, confidence: high] — the `n.ind` warning prints at `verbose = 0` (VRB3, VRB5)**
`R/gl.sim.genotypes.r:45-46` — two bare `cat(warn(...))` calls with no
verbosity gate.
Failure scenario: `capture.output(gl.sim.genotypes(possums.gl, n.ind = 300,
verbose = 0))` returns 2 lines. Verified: an ordinary `verbose = 0` call
returns 0 lines, so this branch is the only leak.
Proposed change: wrap in `if (verbose >= 1)` — the message reports a
result-affecting condition, so `>= 1` rather than `>= 2` (VRB4).

**F7 [MEDIUM, confidence: high] — `n.ind` is never validated (FS5)**
`R/gl.sim.genotypes.r:44-50` — no parameter checking between the datatype
check and the work. Every invalid value dies inside base R.
Failure scenario, verified: `n.ind = 0` gives
`Error: subscript out of bounds`; `n.ind = -5` and `n.ind = "a"` give
`Error: negative length vectors are not allowed`; `n.ind = NA` gives
`Error: missing value where TRUE/FALSE needed`; `n.ind = 2.7` silently
returns 2 individuals with no message. None of these name `n.ind` or the
function.
Proposed change: add an FS5 block requiring a single finite numeric
`n.ind >= 1`, rounding non-integers with a gated warning.

**F8 [MEDIUM, confidence: high] — `ploidy` is built with the wrong length (DAT1)**
`R/gl.sim.genotypes.r:68` — `ploidy = rep(2, n.loc)`. In a genlight, `ploidy`
is per individual, so this should be `rep(2, n.ind)`. It is correct by
accident only because every element is 2; `new("genlight", ...)` recycles or
truncates the vector to `n.ind` without complaint.
Failure scenario: no wrong output today. The line is a latent trap — any
future change to a non-constant ploidy vector, or a stricter adegenet
constructor, produces silently mis-assigned ploidy. Verified: with `n.ind =
50` and `n.loc = 200`, `length(ploidy(out))` is 50 and all values are 2.
Proposed change: `ploidy = rep(2, n.ind)`.

**F9 [MEDIUM, confidence: high] — the returned metadata is incomplete and carries a junk column (DAT2)**
`R/gl.sim.genotypes.r:63-72` — the new genlight is built with no `pop`, no
`ind.metrics` and no `loc.metrics`, and `gl.compliance.check()` is left to
invent them.
Failure scenario, verified: `@other$ind.metrics` has one column, `id` — the
standard `pop` column is missing, so any downstream code that rebuilds
populations from `ind.metrics$pop` finds nothing. `@other$loc.metrics`
gains a column named literally `array(NA, nLoc(x))`, and
`@other$loc.metrics.flags` a matching `array(NA, 1) = NA` entry. Source
population names are discarded: a single-population source named `A` returns
`pop1`.
Root cause of the junk column is in `gl.compliance.check` (scope note below),
but the object this function returns is where it surfaces.
Proposed change: build `pop`, `ind.metrics` (with `id` and `pop`) and an
empty-but-named `loc.metrics` before the compliance call; carry the source
population name through when `nPop(x) == 1`.

**F10 [MEDIUM, confidence: high] — a plain genlight fails with an opaque error (DAT5)**
`R/gl.sim.genotypes.r:36-40` — the input is never passed through
`gl.compliance.check()`, so an object not built by dartR reaches
`gl.allele.freq()` and dies there.
Failure scenario, verified: `new("genlight", gen = m, ploidy = rep(2, 10))`
gives `Error: argument is of length zero`, from a
`x@other$loc.metrics.flags$monomorphs` test inside `gl.allele.freq`.
Proposed change: run `x <- gl.compliance.check(x, verbose = 0)` on the input
before extracting frequencies.

**F11 [MEDIUM, confidence: high] — roxygen block is incomplete (DOC1, DOC2, DOC3, DOC7)**
`R/gl.sim.genotypes.r:1-19` — no `@details`; no `@examples` (DOC3 requires
one and `R CMD check` flags an exported function without them); `@return`
sits after `@export`, against the DOC1 house order; the `verbose` text uses
the old `[default 2 or as specified using gl.set.verbosity]` clause instead
of the DOC2 wording; `@author` gives `Custodian:` only, with no `Author(s):`
line (DOC7); the `@param n.ind` text repeats the "should be less than the
number of loci" advice that F1 shows is not enforced.
Failure scenario: users have no worked example, and the PDF manual reports a
verbosity default that does not describe the `gl.check.verbosity()` cascade.
Proposed change: add `@details` stating the model (per-locus HWE draws from
pooled frequencies, no missing data simulated), add a runnable `@examples`
on `possums.gl`, reorder to the DOC1 sequence, adopt the DOC2 verbose text,
and add `Author(s): Arthur Georges`. Run `devtools::document()` in the same
change (DOC4).

**F12 [LOW, confidence: high] — `verbose = 3` promises a results summary and prints none (VRB1)**
`R/gl.sim.genotypes.r:76-78` — the only output above `verbose = 2` comes from
`utils.flag.start`. Verified: the `verbose = 2` and `verbose = 3` transcripts
are byte-identical (3 lines each).
Failure scenario: a user raising verbosity to see what was simulated learns
nothing.
Proposed change: at `verbose >= 3`, report `n.ind`, `n.loc`, the number of
source individuals the frequencies came from, and the mean simulated
heterozygosity.

**F13 [LOW, confidence: high] — outdated `build =` argument to `utils.flag.start` (FS3)**
`R/gl.sim.genotypes.r:32` — `build = "v.2023.3"` is the 2021-proforma form;
FS3 records it as outdated, and it surfaces in the `verbose = 5` banner as
`Build = v.2023.3` on a package at version 1.2.3.
Proposed change: drop the `build` argument.

**F14 [INFO, confidence: high] — commented-out parameter left in the signature**
`R/gl.sim.genotypes.r:23` — `#error.check = TRUE,` sits in the argument list.
Proposed change: delete, or implement it as the F7 validation switch.

### Scope notes (other functions — recorded, not fixed)

- `gl.compliance.check` is the source of the `array(NA, nLoc(x))` column and
  the `array(NA, 1)` flag when it is handed an object with no `loc.metrics`.
  Any function that builds a bare genlight and hands it to the compliance
  check inherits them.
- `gl.allele.freq` (already at `awaiting-approval`) computes `m$frequency` at
  `R/gl.allele.freq.r:194` with `colMeans(as.matrix(x))/2` after reordering
  `m` by `loc_order` at line 192; the two orderings are assumed to agree.
  Not re-examined here.

## Proposed changes

1. Apply the `n.ind > n.loc` cap, or replace the message with one that does
   not claim a cap (F1). **Consequence: `nInd()` of the returned object
   changes for every call with `n.ind > nLoc(x)`.**
2. Error with a named message when any locus frequency is `NA`, pointing at
   `gl.filter.allna()` (F2). **Consequence: the current opaque
   `sample.int` error is replaced; calls that fail today still fail.**
3. Refuse or warn on a multi-population input, and document that frequencies
   are pooled (F3). **Consequence: multi-population calls that silently
   succeed today either stop or print a warning.**
4. Set `accept = "SNP"` on `utils.check.datatype()` (F4).
   **Consequence: SilicoDArT input is rejected instead of returning a SNP
   object.**
5. Record `match.call()` as the returned object's history, replacing the two
   internal entries (F5).
6. Gate the `n.ind` warning at `verbose >= 1` (F6).
7. Add an FS5 validation block for `n.ind` (F7).
8. Change `ploidy = rep(2, n.loc)` to `rep(2, n.ind)` (F8).
9. Construct `pop`, `ind.metrics` and `loc.metrics` before the compliance
   check; carry a single source population's name through (F9).
10. Run `gl.compliance.check()` on the input object (F10).
11. Rewrite the roxygen block to DOC1/DOC2/DOC3/DOC7 and re-document (F11).
12. Add a `verbose >= 3` results summary (F12).
13. Drop `build =` from `utils.flag.start()` and delete the commented-out
    `error.check` argument (F13, F14).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run.
- Spec: behaviour vs roxygen on possums.gl, testset.gl, testset2.gl,
  platypus.gl, testset.gs — run.
- Simulation model: per-locus HWE chi-square over 200 loci x 2000
  individuals — run. Against the simulated frequency (1 df): mean chi-square
  0.819 (expectation 1), proportion p < 0.05 = 0.030 (expectation 0.05).
  Against the input frequency (2 df): mean 1.782 (expectation 2), proportion
  p < 0.05 = 0.045. Allele-frequency recovery: `cor(input, simulated)` =
  0.99799, max absolute deviation 0.0228 against a 4-SE bound of 0.0316.
- Joint structure (recycling lens, three probes) — run, all clean:
  (i) 83 locus pairs sharing an allele frequency, n = 1000 individuals:
  mean |r| = 0.027, max |r| = 0.094, zero pairs above 0.3;
  (ii) 500 random pairs of individuals, congruent-position match rate
  observed 0.3962 against a chance rate of 0.3935 (z = 1.79);
  (iii) draw accounting: the RNG stream advances by exactly `2 * n.loc`
  vectorised draws of length `n.ind`, i.e. 16,000 independent draws for
  8,000 cells at `n.ind = 40`, and the output is reproduced exactly by an
  explicit two-haplotype reimplementation.
- Seeded reproducibility (TST4 proposed): same seed gives identical
  genotypes and identical `loc.metrics`; a different seed differs — run.
- Missing data: none simulated (0 NAs at every `n.ind` tested); the function
  makes no claim about missingness, so there is nothing to contradict —
  recorded as a documentation gap under F11, not a separate finding.
- Monomorphic loci (frequency exactly 0 or 1): correct — the simulated locus
  stays fixed at 0 or 2 respectively — run on a synthetic fixture; possums.gl
  contains no monomorphic locus.
- Every documented parameter exercised individually: `x` (SNP, SilicoDArT,
  plain genlight, single-locus, single-individual), `n.ind` (0, 1, 2.7, -5,
  "a", NA, `n.loc`, `n.loc + 100`), `verbose` (0, 1, 2, 3, 5) — run.
- Downstream usability: `gl.compliance.check`, `gl.filter.callrate`,
  `gl.recalc.metrics` and `gl.pcoa` all accept the returned object — run.
- Input object untouched — run, `identical()` before and after.
- PLT rules: not applicable — the function draws no plot.
- DEP rules: not applicable — no Suggests package is used.
- FBM path (DAT6): SKIPPED — this checkout has no FBM machinery
  (`gl.genlight2fbm` and `gen2fbm` are both absent), so no FBM fixture can
  be built.
- dartR Google Group / GitHub issue search: SKIPPED — no network access in
  this session.

## Approval (Phase B)

All 14 findings approved by arthur on 2026-09-09, with two design choices
recorded against the escalation-gate items.

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | arthur | Apply the cap: `n.ind <- n.loc`, message unchanged |
| 2 | Approved | arthur | Design choice: simulate all-NA loci AS all-NA rather than error. nLoc and the source's missingness structure are preserved |
| 3 | Approved | arthur | Design choice: simulate PER POPULATION rather than refuse. `n.ind` becomes per-population and is documented as such. Consequence acknowledged: multi-population input now returns structured data where it returned pooled data |
| 4 | Approved | arthur | `accept = "SNP"` |
| 5 | Approved | arthur | |
| 6 | Approved | arthur | |
| 7 | Approved | arthur | |
| 8 | Approved | arthur | |
| 9 | Approved | arthur | |
| 10 | Approved | arthur | |
| 11 | Approved | arthur | |
| 12 | Approved | arthur | |
| 13 | Approved | arthur | Covers F13 and F14 |

## Outcome (Phase C)

Applied 2026-09-09 by Claude (claude-opus-5) on branch
`review-gl.sim.genotypes` from `upstream/dev` (`ddaed27`). All 14 findings
applied; 0 rejected, 0 deferred.

### Engine preservation

The sampling engine is unchanged: two independent full-length haplotype
draws per locus, now made once per population. Verified on a single
population of `platypus.gl` (1000 loci) at `n.ind = 1000`, so the new
`n.ind <= nLoc` cap does not bind and the probes keep the sample size
Phase A used:

- Per-locus HWE, 487 polymorphic loci. Against the simulated frequency
  (1 df): mean chi-square 1.010 (expectation 1), proportion p < 0.05 =
  0.045. Against the input frequency (2 df): mean 1.999 (expectation 2),
  proportion p < 0.05 = 0.043.
- Allele-frequency recovery over 991 loci: `cor(input, simulated)` =
  0.99981, max absolute deviation 0.0267 against a 4-SE bound of 0.0447;
  zero loci fall outside their own 4-SE bound.
- Joint structure, all three probes clean. (i) 4305 equal-frequency locus
  pairs: mean |r| = 0.0257, max |r| = 0.1556, zero pairs above 0.3 (SD of
  r under independence 0.0316). (ii) 500 random individual pairs,
  congruent-position match rate 0.8062 against a chance rate of 0.8058
  (z = 1.01). (iii) Draw accounting on `possums.gl` at `n.ind = 40`: the
  RNG stream advances by exactly `2 * nLoc * nPop` = 4000 vectorised
  draws of length 40 for 80,000 cells, identical to a manual stream of
  the same length, and the output is reproduced exactly by an explicit
  two-haplotype reimplementation drawing from each population's own
  frequencies.
- Seeded reproducibility: the same seed gives identical genotypes and
  identical `loc.metrics`; a different seed differs.

### Finding verification

- **F1** — `gl.sim.genotypes(<200-locus single population>, n.ind = 300)`
  returns 200 individuals; the message is unchanged ("Setting n.ind to
  200") and now describes what happened. On `possums.gl` (200 loci, 10
  populations) the same call returns 2000 individuals, 200 in each
  population.
- **F2** — all five packaged datasets. `possums.gl`, `testset.gl`,
  `testset2.gl` and `platypus.gl` run; `testset.gs` stops under F4.
  `nLoc` is preserved in every case (200, 255, 755, 1000), and the set of
  all-NA loci in the result is identical to the set in the source (0, 3, 3
  and 6 respectively). Missingness is reproduced population by
  population: a locus is missing in the simulated individuals of
  population p exactly when it has no calls in population p of the
  source.
- **F3** — per-population Ho tracking on `possums.gl` at `n.ind = 1000`,
  source against simulated: A 0.3282/0.3221, B 0.4198/0.3962,
  C 0.4280/0.4035, D 0.2753/0.2692, E 0.2945/0.2836, F 0.2820/0.2635,
  G 0.3482/0.3386, H 0.2917/0.2801, I 0.4075/0.3969, J 0.3952/0.3832.
  Overall source Ho 0.3470, simulated Ho 0.3337; the pooled He the old
  code produced is 0.4711. Single-population case: source Ho 0.3282,
  simulated 0.3217.
- **F4** — `gl.sim.genotypes(testset.gs)` stops with "Fatal Error:
  inappropriate object passed to function, found SilicoDArT expecting
  SNP".
- **F5** — `@other$history` holds one entry, the `gl.sim.genotypes(...)`
  call, with the arguments used.
- **F6/F7** — `verbose = 0` returns zero lines both for an ordinary call
  and for a call that trips the `n.ind` cap. `n.ind` values 0, -5, "a",
  NA, `c(10, 20)`, NULL and Inf each stop with a message naming `n.ind`;
  2.7 is rounded to 3 with a `verbose >= 1` warning.
- **F8/F9** — `length(ploidy(out)) == nInd(out)`, all 2.
  `@other$ind.metrics` holds `id` and `pop`; no `array(NA, nLoc(x))`
  column and no `array(NA, 1)` flag; source population names are carried
  through.
- **F10** — a plain `new("genlight", ...)` is accepted and returns a
  simulated object.
- **F11/F12/F13/F14** — roxygen rewritten to the DOC1 order with
  `@details`, runnable `@examples` (2.3 s) and an `Author(s):` line;
  `devtools::document()` run and `man/gl.sim.genotypes.Rd` regenerated;
  `verbose = 3` prints the results summary (7 lines against 3 at
  `verbose = 2`); `build =` and the commented-out `error.check` argument
  removed, so `verbose = 5` no longer reports "Build = v.2023.3".
- Downstream: the returned object passes `gl.compliance.check()` and
  survives `gl.filter.callrate()` and `gl.recalc.metrics()`, including
  when it carries all-NA loci. The input object is unmodified
  (`identical()` before and after).
- Baseline test rerun: 22 blocks, all passing. Every flipped expectation
  is annotated `[approved Fn]` in
  `tests/testthat/test-gl.sim.genotypes.R`; no unexplained diff. The
  model probes were moved onto a 1000-locus single population so the new
  `n.ind` cap does not reduce the sample size they rely on.
- Caller grep across all eight dartRverse clones: no package calls
  `gl.sim.genotypes`; the only references are `@family` cross-links in
  `dartR.base/man/*.Rd`. No breaking caller.

### Scope note discovered during application

Combining F2 and F3 means a locus that is called overall but has no calls
within one population is simulated as missing for that population's
individuals. On `testset.gl` (30 populations, 5 individuals each) this
produces substantially more missing data than the source's three global
all-NA loci would suggest. This is the intended reading of "preserve the
source's missingness structure" under per-population simulation, and it is
stated in `@details`; it is recorded here because it is not obvious from
the finding text.

```json
{
  "function": "gl.sim.genotypes",
  "package": "dartR.base",
  "family": "analysis",
  "skill_version": "2.2.0",
  "commit": "ddaed27",
  "model": "claude-opus-5",
  "model_phase_c": "claude-opus-5",
  "approved_by": "arthur",
  "approved_date": "2026-09-09",
  "design_choices": {
    "F2": "simulate all-NA loci as all-NA, preserving nLoc and the source missingness structure",
    "F3": "simulate per population; n.ind is per population and documented as such"
  },
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1",  "severity": "BLOCKER", "confidence": "high", "rule": "DOC5",  "status": "applied", "change": 1},
    {"id": "F2",  "severity": "BLOCKER", "confidence": "high", "rule": "FS5",   "status": "applied", "change": 2},
    {"id": "F3",  "severity": "HIGH",    "confidence": "high", "rule": "DOC5",  "status": "applied", "change": 3},
    {"id": "F4",  "severity": "HIGH",    "confidence": "high", "rule": "DAT7",  "status": "applied", "change": 4},
    {"id": "F5",  "severity": "HIGH",    "confidence": "high", "rule": "FS8",   "status": "applied", "change": 5},
    {"id": "F6",  "severity": "MEDIUM",  "confidence": "high", "rule": "VRB5",  "status": "applied", "change": 6},
    {"id": "F7",  "severity": "MEDIUM",  "confidence": "high", "rule": "FS5",   "status": "applied", "change": 7},
    {"id": "F8",  "severity": "MEDIUM",  "confidence": "high", "rule": "DAT1",  "status": "applied", "change": 8},
    {"id": "F9",  "severity": "MEDIUM",  "confidence": "high", "rule": "DAT2",  "status": "applied", "change": 9},
    {"id": "F10", "severity": "MEDIUM",  "confidence": "high", "rule": "DAT5",  "status": "applied", "change": 10},
    {"id": "F11", "severity": "MEDIUM",  "confidence": "high", "rule": "DOC1",  "status": "applied", "change": 11},
    {"id": "F12", "severity": "LOW",     "confidence": "high", "rule": "VRB1",  "status": "applied", "change": 12},
    {"id": "F13", "severity": "LOW",     "confidence": "high", "rule": "FS3",   "status": "applied", "change": 13},
    {"id": "F14", "severity": "INFO",    "confidence": "high", "rule": "STY1",  "status": "applied", "change": 13}
  ],
  "coverage_skipped": [
    "DAT6: no FBM machinery in this checkout, no fixture buildable",
    "issue-tracker search: no network access in this session"
  ],
  "joint_structure_probes": {
    "equal_frequency_locus_pairs": {"n": 83, "mean_abs_r": 0.027, "max_abs_r": 0.094, "verdict": "independent"},
    "congruent_position_match": {"expected": 0.3935, "observed": 0.3962, "z": 1.79, "verdict": "chance rate"},
    "draw_accounting": {"cells": 8000, "draws": 16000, "ratio": 2, "verdict": "no recycling"}
  },
  "hwe_verification": {
    "chisq_1df_mean": 0.819, "chisq_1df_prop_sig": 0.030,
    "chisq_2df_vs_input_mean": 1.782, "chisq_2df_prop_sig": 0.045,
    "allele_freq_cor": 0.99799, "allele_freq_max_dev": 0.0228,
    "verdict": "Hardy-Weinberg model confirmed"
  },
  "joint_structure_probes_phase_c": {
    "equal_frequency_locus_pairs": {"n": 4305, "mean_abs_r": 0.0257, "max_abs_r": 0.1556, "verdict": "independent"},
    "congruent_position_match": {"expected": 0.8058, "observed": 0.8062, "z": 1.01, "verdict": "chance rate"},
    "draw_accounting": {"cells": 80000, "draws": 4000, "draw_length": 40, "ratio": 2, "verdict": "no recycling"}
  },
  "hwe_verification_phase_c": {
    "chisq_1df_mean": 1.010, "chisq_1df_prop_sig": 0.045,
    "chisq_2df_vs_input_mean": 1.999, "chisq_2df_prop_sig": 0.043,
    "allele_freq_cor": 0.99981, "allele_freq_max_dev": 0.0267,
    "allele_freq_4se_bound": 0.0447,
    "verdict": "Hardy-Weinberg model preserved"
  },
  "fixtures_candidate": true,
  "fixtures_candidate_finding": "F1",
  "status": "pr-open",
  "pr": null
}
```
