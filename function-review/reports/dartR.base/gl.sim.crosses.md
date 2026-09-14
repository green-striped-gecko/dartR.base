# Review: gl.sim.crosses (dartR.base)

- Family mode: analysis/simulation (returns a newly constructed genlight, so the DAT1 ploidy,
  DAT2/DAT3 metadata and FS8 history checks apply to the offspring object)
- Date: 2026-09-09
- Reviewer: Claude Opus 5 (`claude-opus-5`), dartr-function-review v2.1.0
- Package commit: `ddaed27` (upstream/dev; `git diff upstream/dev -- R/gl.sim.crosses.r` is empty,
  so the reviewed state is upstream/dev unmodified. Working tree is `integration-local` at
  `ed99203` with unrelated uncommitted test edits.)
- Dependency versions: dartR.data 1.2.5, adegenet 2.1.11, R 4.4.2
- Datasets: testset.gl, testset.gs, testset2.gl (inspected for pedigree columns only), plus nine
  constructed fixtures (all six single-pair cross types, two engineered gamete-collision panels,
  a known-NA panel, a 20 x 30 random panel with 10 per cent missing, monomorphic parents, an
  all-NA locus, a single-locus object, and mismatched locus panels)
- Baseline: `tests/testthat/test-gl.sim.crosses.R` (new file, 85 assertions, all passing against
  the reviewed state)

## Verdict

**Standards: Needs work** — the entire standard preamble (verbosity, `utils.flag.start`,
`utils.check.datatype`, argument validation) and the history append sit inside an optional
`if (error.check)` block, `verbose = 0` still prints, the datatype check has no
`accept = "SNP"`, and the roxygen block is documented under the wrong `@name` with its
`@examples` commented out.

**Spec: Rework** — three of the function's documented contracts do not hold. Gametes are not
drawn independently: one random vector of length `mhet` is recycled across the whole genotype
matrix, so heterozygous calls at congruent positions always transmit the same allele, even in
different parents. Marginal Mendelian ratios are exactly right, which is why the defect survives
a doc-versus-output comparison; the joint structure of the simulated cohort is fabricated. The
documented `n` parameter still has no effect, and `error.check = FALSE` — which `@details`
instructs the user to set for simulations — aborts.

What works: the dosage arithmetic `(ova + sperm) / 2` is a correct and neat encoding of Mendelian
transmission; marginal segregation ratios match expectation at every cross type; missing data
propagates exactly (either parent missing gives a missing offspring call, with no spurious calls);
seeded runs are reproducible; the parent objects come back untouched.

## Historical-defect status

**Still present.** The `n` parameter — cited in this skill's own `SKILL.md` as the example of a
documented parameter that is silently never used — is unfixed at `ddaed27`.

Evidence:

- Static: `n` appears twice in `R/gl.sim.crosses.r` outside roxygen — once in the signature
  (line 77) and once in the comparison `if (noff < n)` at line 112, which only decides whether a
  warning prints. No code path subsets the brood.
- Empirical: 3 mothers, `broodsize = 5`, `n = 2` returns 15 individuals.
- Provenance: only three commits have ever touched the file (`e4ea42c` and `8cfbf3a`, both
  2023-09-13; `79ed9c9`, 2023-10-02). All three are ancestors of `upstream/dev`
  (`git merge-base --is-ancestor`), so nothing is stranded. `git log --all --oneline --
  R/gl.sim.crosses.r` lists no other commit on any local or remote branch — the
  local-only-`dev` pattern found for `gl.sample` (`11e8deb`) does not recur here.

The review did not stop at the flagged defect. It found two further BLOCKERs, one of them
(F1, non-independent gamete draws) more consequential than the flagged one, plus five HIGHs.

## Findings

**F1 [BLOCKER, confidence: high] — gamete draws are not independent (spec axis: sampling
correctness; no catalogue rule fits — see the skill-maintainer note)**

`R/gl.sim.crosses.r:119-126` and `:130-137` —
`ova <- ifelse(mmat == 1, sample(c(0, 2), mhet, replace = T), mmat)` draws exactly `mhet` random
values, where `mhet <- sum(mmat == 1, na.rm = TRUE)` is the total number of heterozygous calls in
the whole parent matrix. `ifelse()` then evaluates `rep(yes, length.out = length(test))[ypos]`, so
the het call at column-major position `k` receives draw number `((k - 1) %% mhet) + 1`. Two het
calls whose positions are congruent modulo `mhet` therefore always transmit the same allele — in
the same brood replicate, regardless of which locus or which parent they belong to.

Failure scenario, engineered: a mother scored `1,0,1,2` is het at L1 and L3, giving `mhet = 2`;
both het cells take draw 1. Over 500 offspring the L1 and L3 genotypes are identical in every
single individual (agreement rate 1.000; 0.500 expected). Across parents: three mothers scored
`1,2` / `0,2` / `1,2` give `mhet = 2` with the two het cells at positions 1 and 3, so mothers M1
and M3 — unrelated individuals — transmit the same L1 allele in all 300 brood replicates.

Failure scenario, real data: a 10-mother x 255-locus cohort drawn from `testset.gl` carries 29
heterozygous calls but consumes only 14 distinct random draws; 24 of the 29 het calls share a
draw with at least one other. Empirically, mother 7 at locus 5 and mother 2 at locus 194
transmitted the same maternal allele in 200 of 200 brood replicates. The recycling factor for a
single parent is severe: one `testset.gl` female has 2 het calls in a 255-cell row, so a
two-element draw vector is recycled 127 times.

The consequence is that per-locus segregation ratios are correct while the simulated cohort
carries fabricated linkage between loci and fabricated allele sharing between unrelated parents.
Any downstream use that depends on the joint structure — relatedness or kinship estimation,
parentage assignment, LD, power analysis, Fst between simulated cohorts — is invalid, and the
error is invisible in a per-locus summary.

Proposed change: draw one value per heterozygous call rather than one vector per matrix. Replace
the `ifelse()` idiom with an explicit index assignment, for example

```r
gam <- mmat
het <- which(mmat == 1)                       # NA-safe: which() drops NA
gam[het] <- sample(c(0, 2), length(het), replace = TRUE)
```

inside each brood iteration. Applying this changes every simulated genotype under any given seed.
**Consequence: numerical output changes for every call.**

---

**F2 [BLOCKER, confidence: high] — the documented `n` parameter has no effect (spec axis:
documented parameter is a no-op; the closest catalogue rule, DOC5, is `[proposed]` and cannot
carry BLOCKER, so the principle is named directly)**

`R/gl.sim.crosses.r:77` and `:110-114` — `@param n Number of offspring to retain [default 1000 or
mothers*broodsize whichever is the lesser]`, and `@return A genlight object with n offspring of
both sexes`. The value is read only by `if (noff < n)`, which decides whether a warning prints.
The function always returns `nInd(mothers) * broodsize` individuals.

Failure scenario: a simulation that budgets a fixed cohort size calls
`gl.sim.crosses(fa, mo, broodsize = 3, n = 10)` on 10 mothers and receives 30 individuals. Nothing
in the output records that the request was ignored, and at `verbose = 0` the only clue is a
warning that fires in the opposite case (`noff < n`). The signature default is a flat `1000`, not
the documented `min(1000, mothers * broodsize)`, so the default path is also wrong whenever the
brood total exceeds 1000: 200 mothers with `broodsize = 10` returns 2000 individuals, not 1000.

Proposed change: after building `offmat`, retain a random subset of `min(n, nrow(offmat))` rows,
and subset the sex vector and individual names with it. Set the signature default so that it
resolves as documented, or restate the documentation to match a plain default. Guard `n < 1`.
**Consequence: the returned object shrinks for every caller who passes `n`, and for default
calls where the brood total exceeds 1000.**

---

**F3 [BLOCKER, confidence: high] — `error.check = FALSE`, the setting `@details` recommends for
simulations, aborts (FS6, DOC5)**

`R/gl.sim.crosses.r:111` and `:146` — `noff <- nInd(mothers) * broodsize` is computed inside the
`if (error.check)` block, but the object-construction step uses it for
`ind.names = paste0("Po_", 1:noff)`. `@details` states "Set check.error to FALSE if using this
script in simulations".

Failure scenario: `gl.sim.crosses(fa, mo, broodsize = 5, error.check = FALSE)` fails with
`object 'noff' not found`. The documented fast path for repeated simulation use has never worked.

Proposed change: hoist `noff` (and the verbosity resolution needed by the completion message) out
of the `if (error.check)` block, leaving only the validation and messaging inside it. Better,
derive individual names from `nrow(offmat)` so the count cannot disagree with the matrix — which
also fixes F6's blank names.

---

**F4 [HIGH, confidence: high] — SilicoDArT data is admitted and yields impossible scores
(DAT7, DAT1)**

`R/gl.sim.crosses.r:96,98` — `utils.check.datatype()` is called without `accept = "SNP"`, so its
default admits presence/absence data. The algorithm reads `mmat == 1` as "heterozygote", but in
SilicoDArT 1 means "present" and the valid alphabet is 0/1/NA at ploidy 1.

Failure scenario: `gl.sim.crosses(fa, mo, broodsize = 2)` on two 5-individual cohorts taken from
`testset.gs` returns an object with `ploidy = 2` whose scores include the value 2 — impossible for
presence/absence — in 6.5 per cent of cells. `utils.check.datatype()` subsequently reports the
result as SNP data. Nothing warns. The two `datatype.dad` / `datatype.mum` values the function
does compute are assigned and never used.

Proposed change: pass `accept = "SNP"` to both datatype calls so SilicoDArT input fails with the
standard message, and use the returned values (or drop the assignments).

---

**F5 [HIGH, confidence: high] — `compliance.check = FALSE` returns an object dartR cannot consume
(DAT2, DAT5)**

`R/gl.sim.crosses.r:152-157` — `gl2@other$ind.metrics$sex <- sr` assigns into a `NULL`, because
the freshly constructed genlight has an empty `@other`. The result is a bare list, not a
`data.frame`. When `compliance.check = FALSE` the object is returned with no `loc.metrics`, no
`loc.metrics.flags` and no population assignment.

Failure scenario: `gl.filter.callrate()` on such an object fails with
`incorrect number of dimensions`. Since `compliance.check = FALSE` is the documented way to skip
an expensive step in simulation loops, the output of a fast simulation cannot be filtered,
reported on or recalculated without a manual repair.

Proposed change: build `ind.metrics` as a `data.frame` explicitly
(`data.frame(id = ind.names, sex = sr, stringsAsFactors = FALSE)`), populate
`@other$loc.metrics` (at minimum a `data.frame` with `nLoc` rows) and set `pop()` before the
optional compliance call, so the object is structurally valid whether or not
`gl.compliance.check()` runs.

---

**F6 [HIGH, confidence: high] — non-positive or fractional `broodsize` silently returns extra
offspring with blank names (FS5, DAT2)**

`R/gl.sim.crosses.r:101-103,122,133,145` — the guard prints
`"Error: Brood Size must be a positive interger\n Set to 10"` but never assigns 10, and
`for (i in 1:broodsize)` counts down when `broodsize <= 0`. `noff` is then computed from the
unvalidated value and disagrees with `nrow(offmat)`, so `paste0("Po_", 1:noff)` produces a vector
of the wrong length and adegenet fills individual names with empty strings.

Failure scenario, on 2 mothers: `broodsize = 0` returns 4 individuals, `broodsize = -1` returns 6,
`broodsize = 2.7` returns 4 — in all three cases with every `indNames()` entry an empty string, so
individual names are neither unique nor informative and `gl.compliance.check()` cannot repair them
sensibly. No error is raised.

Proposed change: validate `broodsize` as a positive integer in an FS5 block, either erroring via
`stop(error(...))` or actually assigning the documented fallback of 10 with a gated warning; derive
individual names from `nrow(offmat)`.

---

**F7 [HIGH, confidence: high] — parent cohorts are neither size-checked nor locus-checked
(FS5, DAT5)**

`R/gl.sim.crosses.r:140,147` — `offmat <- (ova + sperm) / 2` requires the two gamete arrays to have
identical dimensions, which requires `nInd(fathers) == nInd(mothers)` and `nLoc(fathers) ==
nLoc(mothers)`. Neither is checked. Locus identity is not checked at all; the offspring adopt
`locNames(mothers)`.

Failure scenario A: 3 fathers with 5 mothers fails with the base R message
`non-conformable arrays`, with no indication of which argument is wrong. Users following the
`@details` recipe hit this immediately — see F8.

Failure scenario B: two cohorts with the same number of loci but different loci are crossed
silently. Fathers scored on panel `Z1..Z4` crossed with mothers scored on `L1..L4` return offspring
labelled `L1..L4` whose paternal contribution comes from unrelated markers. Nothing warns.

Proposed change: add an FS5 block erroring with `stop(error("Fatal Error: ..."))` when
`nInd(fathers) != nInd(mothers)` or `!identical(locNames(fathers), locNames(mothers))`.

---

**F8 [HIGH, confidence: high] — the workflow in `@details` cannot be executed as written
(DOC5) (proposed rule)**

`R/gl.sim.crosses.r:22-56` — the five-step recipe fails at two points on the packaged data it
implies.

Failure scenario: step (a) prescribes `gl.keep.pop(x, pop.list = "male", as.pop = "sex")`. The
`sex` levels in `testset.gl` are `Female` / `Male` / `Unknown`, so the lower-case call fails with
`Fatal Error: no populations listed to keep!`. Step (c) prescribes `gl.subsample.ind()` to build
the two cohorts; because that function subsamples per population, a request for `n = 10` on the
`testset.gl` sex-split objects returns 28 fathers and 368 mothers, and `gl.sim.crosses()` then
fails with `non-conformable arrays`. The polygyny and polyandry scenarios described in `@details`
are unreachable for the further reason that pairing is positional (F11) — sampling mothers with
replacement produces repeated mothers, but each still crosses exactly one father.

Proposed change: rewrite `@details` as a runnable recipe against a packaged dataset, with matched
cohort sizes, and state the equal-cohort requirement. Add a working `@examples` block (F16) that
executes the recipe.

---

**F9 [MEDIUM, confidence: high] — `verbose = 0` is not silent (VRB5, VRB3)**

`R/gl.sim.crosses.r:95,98,113` — `cat(report("father --"))` and `cat(report("mother --"))` are
ungated, and the brood-total warning is emitted with no `if (verbose >= 2)` guard.

Failure scenario: `capture.output(gl.sim.crosses(fa, mo, broodsize = 5, n = 1, verbose = 0))`
returns one line, `"father --mother --  Error: Sum of broods less than specified number of
offspring to return, returning 10 "`. A simulation loop calling the function thousands of times at
`verbose = 0` floods the console. At `verbose = 1` the ungated prefix runs into the completion
message: `"father --mother --Completed: gl.sim.crosses"`. The warning text is also mislabelled
`Error:` for a condition that is not fatal, and the message is emitted through `warn()` while
reading as an error.

Proposed change: gate both `cat(report(...))` calls and the warning behind
`if (verbose >= 2)`, and relabel the warning `Warning:`. Add `\n` where the prefix currently runs
into following output.

---

**F10 [MEDIUM, confidence: high] — history records the internal helper calls, and only when
`error.check = TRUE` (FS8)**

`R/gl.sim.crosses.r:159-169` — the history append is inside the `if (error.check)` block, so the
documented simulation path would produce an object with no record of its own construction (it
currently errors first, F3). When the append does run, `@other$history` already contains the
`gl.compliance.check()` and nested `gl.recalc.metrics()` calls made inside this function.

Failure scenario: `off@other$history` on a normal call has three entries —
`gl.recalc.metrics(x = x, verbose = 0)`, `gl.compliance.check(x = gl2, verbose = 0)` and the
`gl.sim.crosses(...)` call. A user replaying the history cannot distinguish the user-level
provenance from this function's internals, and the first two calls reference local variable names
(`x`, `gl2`) that do not exist in the user's session.

Proposed change: append history unconditionally on the returned object, and reset
`@other$history` to a single entry — this call — before appending, since the object is newly
constructed and has no prior provenance. This is the case the proposed FS12 history-discipline
rule is intended to cover.

---

**F11 [MEDIUM, confidence: high] — crosses are positional, not random, and parentage is not
recorded (DOC5) (proposed rule); plus a metadata gap (DAT2)**

`R/gl.sim.crosses.r:1-8,140,152` — `@title` says "Generates random crosses between fathers and
mothers" and `@description` says "Generates random crosses ... then randomly selects a specified
number of offspring to retain". Because the two gamete arrays are stacked in the same order and
added element-wise, mother *i* is crossed with father *i*, in every brood replicate. The only
randomness in the pairing comes from the caller's prior subsampling. The offspring object records
no parentage at all: `ind.metrics` holds a single `sex` column, individual names are `Po_1..Po_n`,
and all offspring are placed in one population, `pop1`.

Failure scenario: a user cannot check the simulated pedigree, cannot compute realised relatedness
against a truth set, and cannot reproduce the parent assignment. `testset2.gl` ships `sire`, `dam`,
`father`, `mother` and `cohort` columns for exactly this kind of validation; nothing comparable can
be built from this function's output. Row order is the only link back to the parents, and it is
undocumented.

Proposed change: state in `@description` and `@details` that pairing is by position (father *i* x
mother *i*) and that randomisation is the caller's responsibility; add `mother` and `father`
columns to `ind.metrics` carrying the parental `indNames()`, and consider a `cohort` or brood index
column.

---

**F12 [MEDIUM, confidence: high] — `sexratio` direction is undocumented and out-of-range values
are accepted (DOC5) (proposed rule), (FS5)**

`R/gl.sim.crosses.r:13,106-108,152` — `@param sexratio Sex ratio of simulated offspring
[default 0.5]` does not say which sex the number refers to. The code assigns `"female"` when
`runif() < sexratio`, so it is the expected proportion of females — the opposite of the common
convention of quoting the proportion or number of males. The range guard prints "Set to 0.5"
without assigning it.

Failure scenario: `sexratio = 0.9` returns 901 females and 99 males out of 1000. A user intending a
male-biased cohort gets the reverse and nothing warns. `sexratio = 1.7` is accepted after the
misleading warning and returns an all-female cohort whose `sex` factor has a single level, which
breaks any downstream sex-split. `sexratio = 0` and `1` likewise return single-level factors.

Proposed change: rename the documentation to "expected proportion of female offspring", or rename
the parameter; validate the range in an FS5 block and either error or actually assign 0.5; give
`sr` both levels explicitly via `factor(..., levels = c("female", "male"))`.

---

**F13 [MEDIUM, confidence: high] — parental locus metrics are discarded (DAT2, DAT3)**

`R/gl.sim.crosses.r:142-149,155-157` — the offspring object is built from a bare matrix, so
`@other$loc.metrics` is created from scratch by `gl.compliance.check()`. Every parental locus
metric is lost.

Failure scenario: crossing two cohorts taken from `testset.gl` returns an object whose
`loc.metrics` has the right row count (255) but none of `TrimmedSequence`, `AlleleID`, `SNP` or
`RepAvg`. `gl.filter.taglength()`, `gl.filter.secondaries()` and the `gl2fasta` family then fail
or silently mis-handle the object, even though the loci are the same loci.

Proposed change: carry `mothers@other$loc.metrics` (and the sequence-level columns) onto the
offspring object before the compliance call, then let `gl.compliance.check()` recalculate the
frequency-derived columns. This is the case the proposed DAT8 derived-column-refresh rule
addresses: the inherited columns that are still valid (sequence, position) must be kept and the
ones that are not (`CallRate`, `maf`, `FreqHets`) must be recalculated.

---

**F14 [LOW, confidence: high] — the roxygen `@name` does not match the function (FS1, DOC1)**

`R/gl.sim.crosses.r:1` — `@name gl.sim.cross`, but the function is `gl.sim.crosses`, so the
generated help topic is `man/gl.sim.cross.Rd`.

Failure scenario: the manual indexes the function under a name that does not exist. `?gl.sim.crosses`
resolves only because roxygen adds the function name as a second alias. A search of the PDF manual
or the pkgdown index for `gl.sim.crosses` lands on a topic titled `gl.sim.cross`.

Proposed change: set `@name gl.sim.crosses` and re-run `devtools::document()`, removing the stale
`man/gl.sim.cross.Rd`.

---

**F15 [LOW, confidence: high] — the function ships with no examples (DOC1, DOC3)**

`R/gl.sim.crosses.r:62-67` — the `@examples` block is commented with `#` rather than `#'`, so it is
invisible to roxygen. The intended example uses `glSim()`, which returns a plain genlight with no
dartR metadata.

Failure scenario: `?gl.sim.crosses` shows no examples and `R CMD check` never exercises the
function, which is how F2, F3 and F6 have survived since 2023.

Proposed change: add a runnable `@examples` block using `testset.gl` and matched cohorts, per
DOC3/TST1.

---

**F16 [LOW, confidence: high] — roxygen block order and author structure (DOC1, DOC2, DOC7)**

`R/gl.sim.crosses.r:1-71` — `@return` sits after `@export` rather than before `@author`; there is
no `@seealso` linking `gl.subsample.ind` or the `gl.sim.*` family; the `@author` block reads
`Custodian: Bernd Gruber (Post to ...)` with no `Author(s):` part and with the URL in parentheses
rather than after `--`; the `verbose` parameter text ends `[default NULL, unless specified using
gl.set.verbosity]` instead of the DOC2 canonical clause. `@param error.check` and
`@param compliance.check` omit the `[default TRUE]` bracket punctuation used elsewhere (they carry
the value but no full stop).

Failure scenario: the rendered manual page differs in structure from its siblings, and the
custodian's own name is not recorded as an author.

Proposed change: reorder to the DOC1 house order, add
`Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to \url{...}`, and adopt the DOC2
`verbose` wording.

---

**F17 [LOW, confidence: high] — dead code, non-namespaced idioms and typos (STY1, STY3)**

`R/gl.sim.crosses.r:96,98,120,123,131,134,102` — `ova <- array(data = NA, dim = dim(mmat))` and the
matching `sperm` initialisation are overwritten on the first loop iteration and never read;
`replace = T` uses the abbreviation rather than `TRUE`; `datatype.dad` and `datatype.mum` are
assigned and never used; "interger" is misspelt; the two validation messages are labelled
`Error:` although they are warnings and nothing stops.

Failure scenario: `T` is a bindable name — a user or package that has assigned `T <- 0` in the
global environment changes the sampling behaviour of this function. The remainder is a readability
cost only.

Proposed change: delete the dead initialisations, use `TRUE`, drop or use the datatype values, fix
the spelling and the message labels.

---

**F18 [LOW, confidence: medium] — full densification and a brood-sized copy (DAT6)
(proposed rule)**

`R/gl.sim.crosses.r:118,129,124,135` — both parent objects are densified with `as.matrix()`, and
the `rbind(hold, ova)` idiom grows the gamete array by reallocation on every brood iteration, so
peak memory holds roughly two brood-sized double matrices plus `offmat`.

Failure scenario: 50 parents x 5000 loci at `broodsize = 20` allocates about 40 MB per gamete
array, and the quadratic `rbind` pattern reallocates 20 times. FBM-backed parents were tested and
do work — `gl.gen2fbm()` inputs returned a correct 12-individual offspring object — but only by
full densification, which is the DAT6 concern rather than a failure.

Proposed change: preallocate the gamete array at `broodsize * nInd` rows and fill by block rather
than `rbind`; consider looping over brood replicates writing into the preallocated matrix.

## Proposed changes

1. Draw one random value per heterozygous call instead of recycling a single `mhet`-length vector
   (F1). **Consequence: every simulated genotype changes under any given seed; simulated cohorts
   lose the fabricated between-locus and between-parent allele sharing.**
2. Implement `n`: retain a random subset of `min(n, brood total)` offspring, subsetting genotypes,
   names and sex together; make the signature default resolve as documented and guard `n < 1`
   (F2). **Consequence: the returned object shrinks for any caller passing `n`, and for default
   calls whose brood total exceeds 1000.**
3. Hoist `noff` and the verbosity resolution out of the `if (error.check)` block, and derive
   individual names from `nrow(offmat)` so the count cannot disagree with the matrix (F3, F6).
4. Pass `accept = "SNP"` to both `utils.check.datatype()` calls, and either use or drop
   `datatype.dad` / `datatype.mum` (F4). **Consequence: SilicoDArT input now errors instead of
   returning a ploidy-2 object.**
5. Build `ind.metrics` as a `data.frame`, populate `loc.metrics` and set `pop()` so the returned
   object is structurally valid with `compliance.check = FALSE` (F5).
6. Add an FS5 validation block: positive-integer `broodsize` (assign the documented fallback of 10
   or error), `sexratio` in [0, 1], `nInd(fathers) == nInd(mothers)`, and identical locus panels —
   all via `stop(error("Fatal Error: ..."))` where fatal (F6, F7, F12).
7. Gate `cat(report("father --"))`, `cat(report("mother --"))` and the brood-total warning behind
   `if (verbose >= 2)`, relabel the warning `Warning:`, and add the missing newlines (F9).
8. Append history unconditionally, resetting `@other$history` to this single call rather than
   inheriting the internal `gl.compliance.check` / `gl.recalc.metrics` entries (F10).
9. Add `mother` and `father` columns to the offspring `ind.metrics`, carrying parental
   `indNames()` (F11).
10. Carry the parental sequence-level locus metrics onto the offspring object before the
    compliance call, letting the frequency-derived columns be recalculated (F13).
11. Rewrite `@details` as a runnable recipe with matched cohort sizes and correct `sex` labels;
    state that pairing is positional and that `sexratio` is the expected proportion of females;
    add a working `@examples` block (F8, F11, F12, F15).
12. Fix `@name` to `gl.sim.crosses`, reorder the roxygen tags to the DOC1 house order, add
    `Author(s):` alongside `Custodian:`, adopt the DOC2 `verbose` wording, add `@seealso`, and
    re-run `devtools::document()`, removing `man/gl.sim.cross.Rd` (F14, F16).
13. Remove the dead `array()` initialisations, replace `replace = T` with `replace = TRUE`, and fix
    the "interger" spelling and the mislabelled `Error:` messages (F17).
14. Preallocate the gamete arrays instead of growing them with `rbind()` in the brood loop (F18).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run. PLT: not applicable, the function
  produces no plot and has no plot parameter bundle. DEP: not applicable, no Suggests package is
  used; `stats::runif` is declared via `@importFrom`.
- Spec: behaviour versus roxygen — run against nine constructed fixtures plus `testset.gl` and
  `testset.gs`. Every documented parameter was exercised individually: `fathers`, `mothers`,
  `broodsize` (1, 0, -1, 2.7, 7, 200, 500, 4000), `sexratio` (0, 0.5, 0.9, 1, 1.7), `n` (1, 2, 10,
  default 1000), `error.check` (TRUE, FALSE), `compliance.check` (TRUE, FALSE), `verbose`
  (0, 1, 2, 3, 5).
- Mendelian correctness (spec axis 1): run. All six single-pair cross types verified over 4000
  offspring — hom x hom fixed, hom x hom (opposite) all het, het x het 0.251/0.4905/0.2585
  (chi-square p = 0.388 against 1:2:1), hom x het 2009/1991 (binomial p = 0.788). Joint structure
  tested separately and found wrong (F1).
- NA integrity / SNPbin duplicate-index hazard: run and **clear**. The function reads its inputs
  through `as.matrix()` and constructs the offspring object with `new("genlight", gen = ...)`, so
  it never invokes the SNPbin `[` subsetter and cannot hit the duplicate-index NA corruption
  recorded against `utils.dartR.class.def.r`. Verified empirically on a 20 x 30 panel with 10 per
  cent missing: the offspring NA pattern is byte-identical to the parental NA pattern, missing
  calls in either parent always yield a missing offspring call, and no offspring call appears
  where a parent was missing. Repeated parents (the same genotypes supplied as both fathers and
  mothers, and duplicated individuals within a cohort from `gl.subsample.ind(replace = TRUE)`)
  produce no corruption — the positional algorithm indexes each row once per brood replicate.
- Determinism (proposed TST4): run. Two `set.seed(555)` runs give identical genotypes and
  identical sex assignment.
- Input immutability: run. `expect_identical` on both parent objects before and after.
- FBM path (DAT6): run. `gl.gen2fbm()` parents accepted; the offspring object is correct but is
  produced by full densification.
- Downstream consumption: run. On the `compliance.check = TRUE` output,
  `gl.report.callrate`, `gl.filter.callrate`, `gl.recalc.metrics`, `gl.filter.monomorphs`,
  `gl.report.heterozygosity` and `gl.dist.ind` all succeed; a re-run of
  `gl.compliance.check(verbose = 3)` reports only the expected monomorphic and all-NA locus
  notices. On the `compliance.check = FALSE` output, `gl.filter.callrate` errors (F5).
- Caller grep across siblings: run. No `dartR.*` package or report references `gl.sim.crosses`
  except this campaign's own reports; the only in-repo mentions are `man/gl.sim.cross.Rd` and the
  `@details` text.
- Pedigree validation against `testset2.gl` sire/dam truth: SKIPPED — the function records no
  parentage in its output (F11), so simulated relatedness cannot be matched to a truth set without
  first reconstructing pairs from row order. `testset2.gl` was inspected to confirm the columns
  exist (`id, pop, lat, lon, sex, maturity, sire, dam, cohort, offspring, father, mother`).
- `possums.gl` / `platypus.gl`: SKIPPED — neither adds a code path beyond `testset.gl`; the
  function reads no locus metric, so `platypus.gl`'s `TrimmedSequence` is only relevant to F13,
  which was demonstrated on `testset.gl`.
- Plot/results coupling (PLT3): not applicable.
- dartR Google Group / GitHub issue search: SKIPPED — no network access in this session.

## Notes on other functions (scope rule: reported, not fixed)

- `gl.compliance.check` creates a locus-metrics column literally named `array(NA, nLoc(x))` and a
  matching entry in `loc.metrics.flags`. Reproduced on a bare
  `new("genlight", gen = matrix(...))` with no involvement of `gl.sim.crosses`, so the defect is
  in `gl.compliance.check`, not here. Every object this function returns with
  `compliance.check = TRUE` carries the junk column.
- `gl.subsample.ind` subsamples per population, so `n = 10` on a sex-split `testset.gl` returns 28
  and 368 individuals. Whether that is the intended contract is a question for that function's
  review; it is the reason the `@details` recipe here cannot work (F8).
- `testset.gl`'s `sex` levels are `Female` / `Male` / `Unknown`. Several functions' documentation
  assumes lower-case sex labels.

## Note for the skill maintainer

- **New rule needed: documented parameters must affect behaviour.** F2 is the campaign's headline
  fixture for this class, and the only catalogue rule that reaches it is `[proposed]` DOC5, which
  by contract cannot carry BLOCKER severity — so the finding had to name the principle directly.
  Recommend a `[confirmed]` rule (SPEC1, or API4) reading: "Every parameter in the signature must
  influence the returned value or the emitted output on at least one code path. A parameter read
  only by a validation or messaging branch is a defect regardless of documentation." A mechanical
  check is available: for each formal, grep the body for uses outside `if`/`cat`/`stop` blocks.
- **DOC5 leaned on again** — this is the seventh review to cite `[proposed]` DOC5 (F8, F11, F12
  here). Recommend ratifying it as `[confirmed]`; at seven citations across the campaign its
  utility is established, and its `[proposed]` status is now actively distorting severities.
- **Proposed FS12 (history discipline) applies** (F10): a function that constructs a new object
  should reset history to its own call rather than inheriting the entries left by internal helper
  calls. This is the third review to want that rule.
- **Proposed DAT8 (derived-column refresh) applies** (F13): sequence-level locus metrics must be
  carried over and frequency-derived ones recalculated when an object is rebuilt from a matrix.
- **Proposed TST4 (seed reproducibility) applies and passes** here, unlike `gl.fst.pop`.
- **Fixtures candidate.** F1, F2 and F3 all qualify for `references/fixtures.md`. F1 is the
  stronger fixture of the three, because it is the one a doc-versus-output comparison cannot find:
  every marginal segregation ratio is correct. Suggested rows:

  | Fixture | Pre-fix state | Defect the review must find | Expected class | Min severity |
  |---|---|---|---|---|
  | `gl.sim.crosses` | `ddaed27` (see report `gl.sim.crosses.md`) | `ifelse(mmat == 1, sample(c(0, 2), mhet, replace = TRUE), mmat)` draws only `mhet` values and lets `ifelse` recycle them, so het calls at congruent column-major positions always transmit the same allele -- within and between parents. Marginal Mendelian ratios stay correct | Spec axis: sampling correctness (vector recycling in `ifelse`) | BLOCKER |
  | `gl.sim.crosses` (same state) | `ddaed27` | the documented `n` ("offspring to retain") is read only by a warning test; the full brood is always returned | Spec axis: documented parameter is a no-op | BLOCKER |
  | `gl.sim.crosses` (same state) | `ddaed27` | `noff` is computed inside `if (error.check)` and used by the constructor, so the `error.check = FALSE` path `@details` recommends fails with `object 'noff' not found` | FS6 scoping / spec axis | BLOCKER |

## Approval (Phase B)

All 18 findings approved by Arthur Georges on 2026-09-09, with the F1
consequence acknowledged explicitly: every simulated dataset this function
produces changes, because the current output carries fabricated correlation
between loci and between parents.

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | arthur | F1. Consequence acknowledged: all simulated output changes |
| 2 | Approved | arthur | F2 |
| 3 | Approved | arthur | F3, F6 |
| 4 | Approved | arthur | F4 |
| 5 | Approved | arthur | F5 |
| 6 | Approved | arthur | F6, F7, F12 |
| 7 | Approved | arthur | F9 |
| 8 | Approved | arthur | F10 |
| 9 | Approved | arthur | F11 |
| 10 | Approved | arthur | F13 |
| 11 | Approved | arthur | F8, F11, F12, F15 |
| 12 | Approved | arthur | F14, F16 |
| 13 | Approved | arthur | F17 |
| 14 | Approved | arthur | F18 |

## Outcome (Phase C)

Applied 2026-09-09 by Claude Opus 5 (`claude-opus-5`), on branch
`review-gl.sim.crosses` off `upstream/dev` (`ddaed27`). All 18 findings applied
(3 BLOCKER, 5 HIGH, 5 MEDIUM, 5 LOW); none deferred.

### F1 verification gate

| Condition | Before | After |
|---|---|---|
| (i) mother scored `1,0,1,2`: agreement of L1 and L3 over 500 offspring | 1.000 (500/500) | 0.496 |
| (ii) two unrelated mothers (`1,2` / `0,2` / `1,2`) at L1, 300 brood replicates | 1.000 (300/300) | 0.503 |
| (iii) distinct random draws consumed by the 29 het calls in the `testset.gl` 10-mother panel | 14 | 29 |
| (iii, real data) the two het cells the recycling forced to share a draw (mother 7 locus 5, mother 6 locus 8), 200 replicates, homozygous-reference fathers | 1.000 (200/200) | 0.460 |
| (iv) per-locus Mendelian ratios, 4000 offspring | correct | correct |

Condition (iv) in full: `0x0` all 0, `2x2` all 2, `0x2` all 1, `1x1`
0.2535/0.4998/0.2468 (chi-square p = 0.833 against 1:2:1), `0x1` 1958/2042
(binomial p = 0.189), `1x2` 2012/1988 (p = 0.716); only 0/1/2 produced. The
marginals were right before and stay right; only the joint structure changes.
Both cells of the real-data pair still segregate 1:1 individually (0.500 and
0.460).

### Other verification

- **F2**: 3 mothers x `broodsize = 5` with `n = 2` returns 2 (was 15); `n = 6`
  returns 6. Genotypes, individual names and `ind.metrics` are subset together.
  The signature default is now `n = NULL`, resolving to the documented lesser
  of 1000 and the brood total: `broodsize = 5` returns 15, `broodsize = 400`
  returns 1000. `n = 0` and `n = 2.5` raise `Fatal Error: Number of offspring
  to retain (n) must be a positive integer`. An `n` above the brood total warns
  and returns the brood total.
- **F3**: `gl.sim.crosses(fa, mo, broodsize = 5, error.check = FALSE)` returns
  10 individuals of class `dartR` with non-empty names. It previously failed
  with `object 'noff' not found`.
- **F4**: two 5-individual cohorts from `testset.gs` now stop with
  `Fatal Error: inappropriate object passed to function, found SilicoDArT
  expecting SNP`.
- **F5**: with `compliance.check = FALSE`, `ind.metrics` is a `data.frame` with
  columns `id, sex, mother, father`, `loc.metrics` has `nLoc` rows,
  `loc.metrics.flags` is present and `pop` is `pop1`.
  `gl.filter.callrate(out, threshold = 0.5)` completes and returns a `dartR`
  object with the same individual count. On a `testset.gl`-derived cross,
  `gl.filter.callrate(threshold = 0.9)` returns a `dartR` object of 30
  individuals and 192 of the 255 loci.
- **F6**: `broodsize` of 0, -1 and 2.7 all fall back to 10 (20 offspring from
  2 mothers) with the warning `Brood size must be a positive integer. Set to
  10`; names are `Po_1..Po_20`, never blank. Previously 4, 6 and 4 offspring
  with every name an empty string.
- **F7**: 3 fathers with 5 mothers raises `Fatal Error: Pairing is positional,
  so the two parent cohorts must hold the same number of individuals. Found 3
  fathers and 5 mothers` (was `non-conformable arrays`). Equal-length but
  different locus panels raise `Fatal Error: The two parent cohorts have the
  same number of loci but not the same loci` (was silently adopting the
  mother's names). Unequal locus counts are fatal regardless of `error.check`.
- **F8**: the corrected `@details` recipe was run end to end and produces 20
  offspring x 255 loci of class `dartR`. The two corrections are the
  capitalised `sex` labels (`"Male"` / `"Female"`) and assigning a single
  population before `gl.subsample.ind()`, which subsamples per population.
- **F12**: `sexratio = 1.7` assigns 0.5 with a warning, and the `sex` factor
  keeps both levels even at `sexratio = 1`.
- **Preamble and history**: verbosity, `utils.flag.start`, the datatype checks
  and all parameter validation now run unconditionally. History is appended
  unconditionally and reset to this call alone (length 1, was 3 entries
  including the internal `gl.compliance.check` and `gl.recalc.metrics` calls).
  `error.check` now governs only the `identical(locNames(...))` comparison,
  which scales with the number of loci.
- **Verbosity**: `verbose = 0` produces zero lines with `error.check` and
  `compliance.check` in either state; it previously printed
  `father --mother --  Error: Sum of broods ...`. `verbose = 2` prints the
  start banner, both datatype lines and the completion banner.
- **NA integrity**: unchanged. Missing in either parent gives a missing
  offspring call, with no spurious calls; the offspring NA pattern on a
  20 x 30 panel with 10 per cent missing is identical to the parental pattern.
  The SNPbin duplicate-index negative control (repeated parents supplied in
  both cohorts) is exact.
- **Compliance and downstream**: `gl.compliance.check(verbose = 3)` on the
  offspring reports no errors. `gl.report.callrate`, `gl.filter.callrate`,
  `gl.recalc.metrics`, `gl.filter.monomorphs`, `gl.report.heterozygosity` and
  `gl.dist.ind` all complete. The junk `array(NA, nLoc(x))` column noted under
  "Notes on other functions" no longer appears on this function's output,
  because `loc.metrics` is populated before the compliance call.
- **Determinism**: two `set.seed(555)` runs give identical genotypes and
  identical sex assignment.
- **Test file**: 29 blocks, 114 assertions, all passing. Every changed
  expectation carries an `[approved Fn]` marker. The 85-assertion baseline
  produced 27 failures across 13 blocks against the applied code; all 13 map
  to approved findings. Twelve were DEFECT PINs; the thirteenth, the marginal
  Mendelian block, failed only because its 4000-offspring brood is now capped
  by the default `n`, so it passes `n = 4000` explicitly.
- **Caller grep**: no sibling package references `gl.sim.crosses`
  (`dartR.captive`, `dartR.popgen`, `dartR.sim`, `dartR.spatial`,
  `dartR.sexlinked`, `dartRstartup`, `dartRverse` all clean). Within
  `dartR.base` only `NAMESPACE`, the source file, its test file and the
  now-removed `man/gl.sim.cross.Rd` mention it.
- **Docs**: `devtools::document()` run; `man/gl.sim.cross.Rd` removed and
  `man/gl.sim.crosses.Rd` added. `NAMESPACE` is unchanged. Unrelated roxygen
  drift in 34 other `man/*.Rd` files was discarded. The `@examples` block runs
  in 0.6 s.

### Not applied

None. All 14 proposed changes were approved and applied.

```json
{
  "function": "gl.sim.crosses",
  "package": "dartR.base",
  "family": "analysis",
  "skill_version": "2.1.0",
  "model": "claude-opus-5",
  "commit": "ddaed27",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "historical_defect_status": "still_present",
  "historical_defect_note": "the documented `n` parameter is a no-op at ddaed27; all three commits touching the file are ancestors of upstream/dev, no stranded fix on any branch",
  "datasets": ["testset.gl", "testset.gs", "testset2.gl", "constructed fixtures (9)"],
  "baseline_test": "tests/testthat/test-gl.sim.crosses.R",
  "baseline_assertions": 85,
  "findings": [
    {"id": "F1", "severity": "BLOCKER", "confidence": "high", "rule": "spec:sampling-correctness", "status": "proposed", "change": 1},
    {"id": "F2", "severity": "BLOCKER", "confidence": "high", "rule": "spec:parameter-no-op", "status": "proposed", "change": 2},
    {"id": "F3", "severity": "BLOCKER", "confidence": "high", "rule": "FS6", "status": "proposed", "change": 3},
    {"id": "F4", "severity": "HIGH", "confidence": "high", "rule": "DAT7", "status": "proposed", "change": 4},
    {"id": "F5", "severity": "HIGH", "confidence": "high", "rule": "DAT2", "status": "proposed", "change": 5},
    {"id": "F6", "severity": "HIGH", "confidence": "high", "rule": "FS5", "status": "proposed", "change": 3},
    {"id": "F7", "severity": "HIGH", "confidence": "high", "rule": "FS5", "status": "proposed", "change": 6},
    {"id": "F8", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "proposed", "change": 11},
    {"id": "F9", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "proposed", "change": 7},
    {"id": "F10", "severity": "MEDIUM", "confidence": "high", "rule": "FS8", "status": "proposed", "change": 8},
    {"id": "F11", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "proposed", "change": 9},
    {"id": "F12", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "proposed", "change": 6},
    {"id": "F13", "severity": "MEDIUM", "confidence": "high", "rule": "DAT2", "status": "proposed", "change": 10},
    {"id": "F14", "severity": "LOW", "confidence": "high", "rule": "FS1", "status": "proposed", "change": 12},
    {"id": "F15", "severity": "LOW", "confidence": "high", "rule": "DOC3", "status": "proposed", "change": 11},
    {"id": "F16", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "proposed", "change": 12},
    {"id": "F17", "severity": "LOW", "confidence": "high", "rule": "STY3", "status": "proposed", "change": 13},
    {"id": "F18", "severity": "LOW", "confidence": "medium", "rule": "DAT6", "status": "proposed", "change": 14}
  ],
  "coverage_skipped": [
    "pedigree validation against testset2.gl sire/dam: the function records no parentage in its output (F11)",
    "possums.gl / platypus.gl: no additional code path",
    "PLT3 plot/results coupling: not applicable, no plot",
    "DEP1 dependency guards: not applicable, no Suggests package used",
    "dartR Google Group / GitHub issue search: no network access in this session"
  ],
  "status": "pr-open",
  "approved": "all",
  "approved_by": "arthur",
  "approved_date": "2026-09-09",
  "applied_by": "claude-opus-5",
  "applied_date": "2026-09-09",
  "branch": "review-gl.sim.crosses",
  "test_assertions_after": 114,
  "pr": 388
}
```
