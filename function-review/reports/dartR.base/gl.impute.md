# Review: gl.impute (dartR.base)

- Family mode: modify
- Date: 2026-09-06
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ddaed27 (upstream/dev)
- Datasets: testset.gl, testset.gs (dartR.data 1.2.5), synthetic fixtures
- Baseline: tests/testthat/test-gl.impute.R (48 assertions, all passing, captured pre-review)

State handling: `git diff upstream/dev -- R/gl.impute.r` is empty, so `load_all()`
exercised the reviewed code. Three helpers differ in the working tree:
`utils.impute.R` (doc conversion plus the `parallel = TRUE` fix, approved diff I7 —
sampler semantics identical at ddaed27), `gl.compliance.check.r` (position-slot fix,
metric recalculation unchanged), and `gl2vcf.r` (PR #352). The beagle path is
therefore reviewed statically against ddaed27's copies; all empirical results for
the four native methods are valid for ddaed27.

## Verdict

**Standards: Needs work** — the FS anatomy, history append, metadata resync and
verbose-0 silence all check out empirically for the SNP paths, but `method` is never
validated, the beagle guards break DEP1/VRB5, and two branches report diagnostics
for the wrong population.

**Spec: Rework** — the documented contrast between `frequency` and `HW` does not
exist (both are the same Hardy–Weinberg draw), the claimed presence–absence support
corrupts or crashes, and the FBM path silently converts residual NAs to genotype 0.

What works well: population-wise imputation, individual-order restoration after
`seppop()`/`rbind()` (the 9897d48 fix), locus-metric recalculation via the internal
compliance pass, and seeded reproducibility are all verified correct.

## Method-by-method verification

| Method | Docs claim | Verified behaviour | Match |
|---|---|---|---|
| frequency | "imputed using the average allele frequencies at that locus in the population" | Random draw of two Bernoulli(q) alleles per cell (`s_alleles`), per population. At q=0.5, imputed 0/1/2 at 28/47/25% (n=1249) — a HW draw, not an average fill | **No (F1)** |
| HW | "sampling at random assuming Hardy-Weinberg equilibrium" | Trinomial draw from (p², 2pq, q²) per population — 28/49/23% at q=0.5. Distributionally identical to `frequency` | Yes, but ≡ frequency (F1) |
| neighbour | nearest Euclidean neighbour, then next nearest, until filled | Confirmed: copies from the nearest non-NA carrier, falls through to the second nearest, deterministic across runs; ties broken by individual index (stable `order()`); self excluded; distance matrix global, computed once from the un-imputed data | Yes |
| random | "random values (0, 1 or 2)" | Uniform over 0:2 (33/33/33%), ignoring frequencies | Yes (but invalid for SilicoDArT, F3) |
| beagle | BEAGLE via VCF export | Static review only; guard defects and chromosome mutation (F7, F8) | Skipped (no jar) |
| fill.residual | residuals "set to 0, 1, 2 at random, taking into account global allele frequencies" | Confirmed global HW draw; but globally-all-NA loci remain NA (F9) | Partial |

Cross-cutting, verified on `testset.gl[1:40, 1:80]` (288 NAs, 2 all-NA loci):
originally scored cells byte-identical after every method; dims, `indNames`,
`locNames`, `pop`, ploidy unchanged; `ind.metrics$id` still aligned; `CallRate`
recalculated (1 for imputable loci, 0 for all-NA loci — DAT2/DAT3/DAT4 satisfied via
the internal `gl.compliance.check`); history appended once (FS8); `verbose = 0`
fully silent for all four native methods (VRB5); `set.seed()` reproduces stochastic
results exactly. No imputation path subsets a genlight with repeated indices, so the
known SNPbin duplicate-index defect is not in play here.

## Findings

**F1 [HIGH, confidence: high] — `frequency` does not do what it documents (DOC5)**
`R/gl.impute.r:164-169` with `R/utils.impute.R` (`s_alleles`) — the "frequency"
method draws two alleles Bernoulli(q), which is exactly a Hardy–Weinberg draw;
`sample_genotype` (the "HW" method) draws from (p², 2pq, q²). The two documented
methods are distributionally identical (verified: 350/588/311 vs 349/615/285 over
1249 imputed cells at q=0.5).
Failure scenario: a user selects `frequency` expecting the documented deterministic
average-frequency fill (expected value 2q, or the most frequent genotype) and gets a
stochastic draw; results differ run to run and method choice between `frequency` and
`HW` changes nothing.
Proposed change: either implement a genuine frequency fill for `method="frequency"`
or re-document both methods as the same HW draw and deprecate one.

**F2 [HIGH, confidence: high] — FBM write-back silently converts residual NA to genotype 0 (DAT1)**
`R/gl.impute.r:328` (also 171, 206, 367, 457) — `x3@fbm[] <- x_matrix` coerces
NA/NaN to 0 in the FBM.code256 backing (bigstatsr warns "nan -> 0"). Verified: the
FBM object holds all 288 NAs before imputation; after `method="neighbour"` the 80
cells that legitimately remain NA in the dense run (two all-NA loci) come back as
genotype 0 — homozygous reference — with no NA left. All other cells agree with the
dense run exactly.
Failure scenario: any FBM-backed object with all-NA loci (or `fill.residual=FALSE`)
gets fabricated 0 genotypes; `gl.pcoa`'s FBM path (PR #369) depends on
`gl.impute(method="neighbour")` and would ingest them silently.
Proposed change: write residual NAs through the code256 NA code (as `gl.gen2fbm`
does), or stop with an informative error when residual NAs remain in FBM mode.

**F3 [HIGH, confidence: high] — SilicoDArT admitted but corrupted or crashed (DAT7, DAT1, DOC5)**
`R/gl.impute.r:116` — `utils.check.datatype` is called with the permissive default,
and the description promises "imputes the state for presence-absence data", but the
samplers are diploid. Verified on `testset.gs[1:30, 1:60]`: `random` writes 16
genotype-2 values into 0/1 data; `frequency` writes 23 twos, which pushes `glMean`
of the part-imputed object to 1.172, so the `fill.residual` pass stops with the
opaque error "negative probability" (`HW` fails the same way). Only `neighbour`
preserves the 0/1 domain.
Failure scenario: a presence-absence user gets either corrupt data (`random`) or an
uninterpretable crash (`frequency`, `HW`).
Proposed change: dispatch on `datatype` — restrict `frequency`/`HW`/`beagle` to SNP
data with a clear error, make `random` sample 0/1 for ploidy-1 data, keep
`neighbour` as is.

**F4 [MEDIUM, confidence: high] — `method` never validated (FS5)**
`R/gl.impute.r:118-122` — the only check is the "neighbor" alias; any other typo
falls through every branch and dies at `glMean(x3)` with "object 'x3' not found".
Failure scenario: `gl.impute(gl, method="freq")` produces an internal-variable error
with no hint at the cause.
Proposed change: validate `method` against the five allowed values and
`stop(error(...))` with the list.

**F5 [MEDIUM, confidence: high] — all-missing detection triggers at >50% missing (VRB3, DOC5)**
`R/gl.impute.r:135, 247, 337, 395` — `loci_all_nas <- sum(glNA(y) > nInd(y))`;
`glNA` counts alleles (2 per missing diploid genotype), so the condition is true for
any locus with more than 50% missing, not only all-missing loci (verified: a
67%-missing locus produces "Population A has 1 loci with all missing values"). The
neighbour branch's `>= nInd(yy)` variant fires at ≥50%. The derived
"values to be imputed" count (`nas_number - loci_all_nas * nInd`) is wrong whenever
the miscount is non-zero, and `sum(glNA(y))/2` halves the true NA count for
ploidy-1 data.
Failure scenario: verbose output warns of all-missing loci that are not, and reports
negative or understated imputation counts.
Proposed change: test `glNA(y) == ploidy * nInd(y)` (ploidy-aware) and derive
`nas_number` from `glNA(y) / ploidy`.

**F6 [MEDIUM, confidence: high] — random/beagle diagnostics report only the last population (STY1, VRB3)**
`R/gl.impute.r:336-359` and `394-417` — the per-population loop computes
`loci_all_nas` but the warning block sits after the loop, so only the final
population's numbers are ever printed (contrast the frequency/HW branches, where the
block is inside the loop).
Failure scenario: a population with all-missing loci earlier in the pop list is
never reported at any verbosity.
Proposed change: move the warning block inside the loop, as in the frequency branch.

**F7 [MEDIUM, confidence: medium] — beagle path guards and silence (DEP1, VRB5, FS5)**
`R/gl.impute.r:373-381, 421-432` — the R.utils guard uses `cat(error(...))` plus
`return(-1)` instead of `stop(error(...))`, so a missing package returns `-1` in
place of a genlight; there is no existence check for the beagle jar, the PLINK
binary, or java (verified: the failure surfaces as the raw
`'C:\...\Temp\...\plink' not found`); and `gl2vcf` is called without `verbose`, so
the path prints "Starting gl2vcf ..." even at `verbose = 0`.
#352 interaction (known, not re-found): the call sets `x_tmp@position <- 1:nLoc`
and passes no `snp.pos`/`snp.chr`, which under #352's explicit-arg precedence still
resolves to the slots — compatible; note the empirical probe above exercised the
working-tree (post-#352) `gl2vcf`, not ddaed27's.
Proposed change: DEP1 stop idiom; `file.exists()` checks on the jar and PLINK paths
before running; pass `verbose = 0`-consistent verbosity through to `gl2vcf`.

**F8 [MEDIUM, confidence: medium] — beagle branch permanently blanks singleton-scaffold chromosome names (DAT2)**
`R/gl.impute.r:384-389` with `:464` — the singleton-scaffold rename (to `""`) is
applied to `x` itself, and `x3$chromosome <- x@chromosome` copies the mutated factor
into the returned object. Static finding (path not executed — no jar).
Failure scenario: after beagle imputation the user's chromosome assignments for
single-SNP scaffolds are silently erased in the returned data.
Proposed change: apply the rename to the throwaway `x_tmp` only and copy chromosome
names from `x_hold`.

**F9 [LOW, confidence: high] — residual-NA exceptions undocumented (DOC5)**
Roxygen `@return` ("missing data imputed") and the `fill.residual` text promise
complete filling, but globally-all-NA loci always remain NA (`glMean` is NaN, and
`s_alleles` returns NA) — verified with `fill.residual = TRUE`. The details section
warns only about per-population all-NA loci. Related interaction, one line per the
scope rule: at ddaed27 `parallel = TRUE` crashes in `matrix2gen` ("object 'i' not
found") — already pinned and approved as diff I7 in `test-utils.impute.R`.
Proposed change: document that loci missing across all individuals are never imputed
and recommend `gl.filter.allna` (no `by.pop`) as the companion filter.

**F10 [LOW, confidence: medium] — full densification on FBM-backed objects (DAT6, proposed rule)**
`R/gl.impute.r:165, 199, 277, 363, 452` — every path materialises the whole matrix
with `as.matrix(x)` even when FBM-backed, and frequency/HW additionally
`seppop()`/`rbind()` per population.
Failure scenario: a genuinely large FBM object is densified in RAM, defeating the
backing.
Proposed change: none at function level yet — flag for the team's FBM strategy.

**F11 [LOW, confidence: high] — documentation defects (DOC1, DOC2, DOC7, DOC3)**
`R/gl.impute.r:1-94` — roxygen order departs from the ratified sequence (`@param`
before `@details`; `@return` after `@export`); the `verbose` text uses the outdated
default clause ("default 2 or as specified using gl.set.verbosity"); `@author` names
only a custodian (no `Author(s):` line; proposed rule DOC7); `beagle.bin.path` and
`plink.bin.path` defaults read "[default getwd())]" (stray parenthesis); and the
presence-absence example is broken — `gs <- gl.filter.callrate(testset.gs, ...)` is
followed by `gl <- gl.filter.allna(gl)` (filters the SNP object again) and the
actual `gl.impute(gs, ...)` line is commented out, so the P/A example never runs the
function. The accepted "neighbor" spelling alias is undocumented.
Proposed change: one documentation pass fixing order, verbose text, author line,
default typos, the example, and documenting the alias.

**F12 [INFO, confidence: high] — minor style and edge notes**
`utils.flag.start(build = "v.2023.3")` uses the argument the conventions mark
outdated (FS3). The frequency/HW banner prints once per population rather than once.
The neighbour path prints one `important()` line per un-imputable individual (40
lines on the test fixture) where a single summary would serve. An imputed result
that is entirely monomorphic crashes inside the mandatory compliance pass
("Subsetting resulted in zero loci") — pathological input, noted under other
functions below.

Other-function notes (scope rule, one line each):
- `gl.filter.monomorphs`/`gl.drop.loc` cannot drop all loci — "Subsetting resulted
  in zero loci" surfaces through `gl.impute`'s compliance pass on fully-monomorphic
  imputed data.
- `matrix2gen(parallel = TRUE)` crashes at ddaed27 — already covered by approved
  diff I7 (`test-utils.impute.R`).
- `gl2vcf` explicit-arg precedence and REF handling — known, PR #352; not re-found.

## Proposed changes

1. Resolve the `frequency`/`HW` duplication: implement a genuine
   average-frequency fill for `method = "frequency"`, or re-document both as HW
   draws and deprecate one (F1). **Consequence: numerical output changes for
   `method = "frequency"` if reimplemented.**
2. Preserve NA in the FBM write-back (code256 NA code), or error when residual NAs
   remain in FBM mode (F2). **Consequence: FBM results change where residual NAs
   existed (0 becomes NA).**
3. Datatype dispatch for SilicoDArT: block `frequency`/`HW`/`beagle` with a clear
   error, ploidy-aware `random` (F3). **Consequence: presence-absence callers of
   the blocked methods now error instead of receiving corrupt output.**
4. Validate `method` against the allowed set with an informative stop (F4).
5. Ploidy-aware all-missing detection and imputation counts (F5).
6. Move the per-population warning block inside the loop in the random and beagle
   branches (F6).
7. Beagle guards: DEP1 stop idiom, existence checks for the jar/PLINK/java, pass
   verbosity through to `gl2vcf` (F7).
8. Confine the singleton-scaffold chromosome rename to the temporary export object
   (F8).
9. Documentation pass: residual-NA exceptions, roxygen order, verbose text, author
   line, default typos, working P/A example, "neighbor" alias (F9, F11).
10. Team-level: FBM-aware (column-wise) processing strategy for imputation (F10,
    proposed rule).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, STY — run. PLT — trivially satisfied (no
  plotting in this function).
- Spec, empirical per method: `frequency`, `HW`, `neighbour`, `random` on synthetic
  fixtures (known q per population, known nearest neighbours) and
  `testset.gl[1:40, 1:80]` — run.
- Post-conditions: untouched-cell identity, dims/names, residual-NA accounting,
  metadata resync, history, `verbose = 0` silence, seeded reproducibility — run.
- Edge cases: pop-all-NA locus, globally-all-NA locus, all-NA individual,
  single-individual populations, no-missing no-op, invalid method — run.
- SilicoDArT: `random`, `frequency`, `HW`, `neighbour` on `testset.gs[1:30, 1:60]`
  — run.
- FBM (DAT6): `neighbour`, `frequency` via `gl.gen2fbm` — run.
- beagle: SKIPPED empirically — no beagle jar/PLINK/java fixture; an empirical run
  would also exercise the working-tree (post-#352) `gl2vcf` rather than ddaed27's.
  Static review done (F7, F8).
- SNPbin duplicate-index hazard: static — no imputation path subsets a genlight
  with repeated indices (the only fancy indexing is a `match()` permutation at
  `R/gl.impute.r:230`).

## Approval (Phase B)

Approved 2026-09-06 by Arthur Georges via the formal approval boxes, with
the stated consequences acknowledged: change 1 changes numerical output for
every past `method = "frequency"` user (they were getting HW draws); change
2 changes FBM results where residual NAs existed (0 becomes NA); change 3
re-specifies frequency/HW/random for the 0/1 domain rather than blocking
them (documented presence-absence support becomes real; beagle blocked).

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur Georges | option: implement the documented deterministic fill (not re-document); statistic recorded in Outcome |
| 2 | approved | Arthur Georges | option: NA-preserving write-back (code256 NA code), not stop-on-residual |
| 3 | approved | Arthur Georges | option: domain-aware methods (band-frequency Bernoulli for HW-analogue, deterministic band fill for frequency, uniform 0:1 for random, neighbour unchanged); beagle blocked for SilicoDArT |
| 4 | approved | Arthur Georges | |
| 5 | approved | Arthur Georges | |
| 6 | approved | Arthur Georges | |
| 7 | approved | Arthur Georges | stop(error()) idiom, gated output, jar/java/PLINK presence checks |
| 8 | approved | Arthur Georges | |
| 9 | approved | Arthur Georges | includes documenting globally-all-NA residual behaviour |
| 10 | approved as recorded | Arthur Georges | no function-level action; flagged for the team FBM strategy |

F12 (INFO): no action, per the approval.

## Outcome (Phase C)

Applied 2026-09-06 on branch `review-gl.impute` (base upstream/dev ddaed27).

Changes 1-9 applied (F1-F9, F11); change 10 (F10) no function-level action,
per the approved proposal; F12 no action.

Option choices recorded precisely:

- **F1 statistic implemented (SNP):** imputed genotype = the expected
  dosage 2q (q = alternate-allele frequency at the locus in the
  individual's population, from `glMean`) mapped to the nearest valid
  genotype: 2q < 0.5 -> 0; 0.5 <= 2q <= 1.5 -> 1; 2q > 1.5 -> 2. The
  boundary cases 2q = 0.5 and 2q = 1.5 are exact ties between two valid
  genotypes and resolve to the heterozygote (1). Per population;
  deterministic; NA frequency (all-NA locus) stays NA. Recorded in the
  roxygen details in these words.
- **F1/F3 statistic implemented (SilicoDArT):** imputed state = the
  majority band state in the population (1 if band frequency >= 0.5 else
  0); the tie at exactly 0.5 resolves to presence (1).
- **F2:** residual NAs are written to the FBM backing as raw code 3, which
  decodes to NA under the house `bigsnpr::CODE_012` mapping (the same
  convention `gl.gen2fbm` writes); dense and FBM-backed runs now return
  identical genotypes including residual NAs.
- **F3:** HW-analogue for SilicoDArT = Bernoulli(band frequency) per
  population; `random` draws uniformly from 0:1; `fill.residual` draws
  Bernoulli from the global band frequency; `neighbour` unchanged;
  `beagle` stops with an informative error for SilicoDArT.
- **F5:** detection is `glNA(y) == ploidy * nInd(y)` with
  `nas_number = sum(glNA(y))/ploidy`; also applied to the same defective
  idiom in the `verbose >= 3` end-of-run summary counts (same finding
  class, same fix).
- **F7:** R.utils guard now `stop(error(...))`; `file.exists()` checks for
  `beagle.27Feb25.75f.jar` and `plink`/`plink.exe`; `Sys.which("java")`
  check; `verbose` passed through to `gl2vcf`; the beagle/java `system()`
  call output suppressed below verbose 2.

Verification (test-gl.impute.R, 75 assertions, all passing):

- Every flipped assertion is annotated `# [approved Fn]`; unflipped
  baseline assertions (untouched-cell identity, dims/names/pop/ploidy,
  CallRate recalculation, history +1, seeded reproducibility, verbose-0
  silence) pass unchanged.
- F1: same-seed and different-seed `frequency` runs identical; imputed
  cells match the hand-computed statistic for every cell on both the
  200-locus q-fixture and a two-population fixture covering both tie
  cases; now distributionally distinct from `HW` (which still scatters
  0/1/2 at q = 0.5).
- F2: dense vs FBM `neighbour` and `frequency` runs on
  `testset.gl[1:40, 1:80]` cell-identical including the 80 residual NAs
  (NA positions `identical`; values equal -- the FBM decode returns
  doubles where the dense path returns integers).
- F3: all four native methods on `testset.gs[1:30, 1:60]` stay in 0/1,
  originally scored cells untouched, ploidy 1 preserved; `frequency`
  deterministic across seeds (with `fill.residual = FALSE`, since the
  residual fill is a documented draw); majority-band fill verified
  cell-level including the 0.5 tie.
- F4-F6: unknown method -> "method must be one of ..."; a 67%-missing
  locus no longer warns; a genuinely all-missing locus warns with the
  right count; the random branch names the affected population (popA)
  rather than only the last.
- Beagle: static + gated tests only (no jar/PLINK/java fixture) --
  SilicoDArT block and missing-jar guard exercised; the jar-present path
  remains statically reviewed.
- Verbose-3 end-to-end runs on `testset.gl` (all four native methods) and
  `testset.gs` (`frequency`) clean; verbose-0 silent apart from a known
  sibling leak (below).

Caller grep (all 8 clones): `gl.pcoa` and `gl.mahal.assign` (dartR.base),
`gl.assign.mahal`/`gl.assign.mahalanobis`/`gl.assign.pca` (dartR.captive),
`gl.assign.mahalanobis`/`gl.assign.pca` (dartR.popgen) all call the default
`method = "neighbour"` on dense objects -- semantics byte-identical, no
caller breaks. `gl.sfs` (dartR.popgen) mentions gl.impute in a warning
string only. dartr2shiny: not available in this workspace (not checked).
`gl.pcoa`'s FBM path fails at ddaed27 with "'k' must satisfy 0 < k <
nrow(A)" before and regardless of gl.impute (verified identical with the
unmodified function) -- pre-existing, out of scope.

NEWS.md entry added, leading with the `frequency` semantics change.

Addendum (found during apply, not applied):

- **F13 [MEDIUM, confidence: high] -- FBM path imputes the caller's input
  object in place.** `x3 <- x` shares the FBM backing file, so
  `x3@fbm[] <- ...` mutates the input object (verified: the input's NA
  count drops from 288 to 80 after `gl.impute(fb, method = "neighbour")`).
  Pre-existing at ddaed27 (the baseline behaviour wrote the fabricated
  zeros into the input too); after F2 the in-place copy at least equals
  the returned object. A fix needs a fresh-backing-file policy and belongs
  with the team FBM strategy (change 10).
- Baseline-test note: the verbose-0 silence assertion filters one known
  sibling leak -- at ddaed27 `utils.dist.ind.snp` prints its banner
  ungated at any verbosity (fix pending on `review-utils.dist.ind.snp`);
  gl.impute itself prints nothing at verbose 0.

```json
{
  "function": "gl.impute",
  "package": "dartR.base",
  "family": "modify",
  "skill_version": "1.0.0",
  "commit": "ddaed27",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DAT1", "status": "applied", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "DAT7", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "applied", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "VRB3", "status": "applied", "change": 5},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "STY1", "status": "applied", "change": 6},
    {"id": "F7", "severity": "MEDIUM", "confidence": "medium", "rule": "DEP1", "status": "applied", "change": 7},
    {"id": "F8", "severity": "MEDIUM", "confidence": "medium", "rule": "DAT2", "status": "applied", "change": 8},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 9},
    {"id": "F10", "severity": "LOW", "confidence": "medium", "rule": "DAT6", "status": "approved-no-action", "change": 10},
    {"id": "F11", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 9},
    {"id": "F12", "severity": "INFO", "confidence": "high", "rule": "FS3", "status": "no-action", "change": null},
    {"id": "F13", "severity": "MEDIUM", "confidence": "high", "rule": "DAT6", "status": "addendum-open", "change": null}
  ],
  "coverage_skipped": [
    "beagle empirical: no beagle jar/PLINK/java fixture; working-tree gl2vcf post-dates reviewed SHA (#352)"
  ],
  "approved_by": "Arthur Georges",
  "approved_date": "2026-09-06",
  "status": "pr-open",
  "pr": 373
}
```
