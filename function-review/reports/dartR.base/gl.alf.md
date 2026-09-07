# Review: gl.alf (dartR.base)

## Provenance

- Model: Claude Fable 5 (claude-fable-5, Claude Code) via dartr-dev agent;
  Skill: dartr-function-review v2.0.0; Base: upstream/dev at ddaed27
  (`git diff upstream/dev -- R/gl.alf.r` empty — the loaded code is the
  reviewed code); working branch integration-local at ed99203.
- Reviewed with a family-level redundancy analysis of the
  allele-frequency functions, requested by the custodian. Companion
  report: `gl.allele.freq.md` (same wave, awaiting approval).
- Datasets: testset.gl, testset.gs, possums.gl, platypus.gl
  (dartR.data 1.2.5), constructed plain-genlight and duplicate-name
  fixtures, FBM conversion of possums.gl.
- Family mode: analysis (pure per-locus accessor; `@family utilities`).
- Baseline: tests/testthat/test-gl.alf.R (20 tests, 64 assertions, all
  pass at the reviewed state).
- Checks skipped: Google Group not searched (not available: no browser
  session). Zero-locus and zero-individual input not reachable — the
  dartR `[` method errors first ("Subsetting resulted in zero loci").

## Verdicts

**Standards: Needs work** — the preamble-free accessor style is
defensible and is not reported as a finding, following the gl.Ho/gl.He
ruling (PR #273: "the pure-function style is judged appropriate for an
accessor whose vector IS the product"). What remains is a roxygen header
below template: no `@description`, no `@details`, no `@seealso`, a
DOC7-noncompliant `@author`, and a `@return` naming columns the function
does not produce.

**Spec: Needs work** — the arithmetic is exact. `alf2` equals
`colMeans(as.matrix(x), na.rm = TRUE)/2` cell for cell on testset.gl,
possums.gl and platypus.gl, `alf1 + alf2` is 1 to 1e-12 wherever
defined, and all-NA loci give NaN in both columns. Two behaviours break
the contract: duplicate locus names silently replace the row keys with
integers (and one caller then produces negative locus counts), and
SilicoDArT is admitted but divided by the SNP ploidy, returning half the
presence frequency.

What works well: the NA policy is the same per-locus `na.rm = TRUE`
policy that gl.allele.freq applies within a population x locus cell, so
the two functions agree to 4 dp on SNP data; and the function is the
only member of the family that tolerates a genlight not built by dartR.

## Findings

**F1 [HIGH, confidence: high] — duplicate locus names silently drop the locus keys (DAT3)**
`R/gl.alf.r:24-25` — `colMeans()` returns a named vector; `data.frame()`
accepts those names as row names only when they are unique, and
otherwise falls back to `1:n` without warning. Two live callers read
`rownames(gl.alf(x))` as locus names.
Failure scenario: with one duplicated locus name in a 30-locus subset of
testset.gl, `gl.report.heterozygosity(method = "pop")` returns
`polyLoc = 30` for every population (all loci) and `monoLoc` values of
`-1` and `-2` — negative counts, no error, no warning (verified against
the same object with unique names, which returns `polyLoc` 3-4 and
`monoLoc` 24-26). `gl.report.heterozygosity.r:509` builds `loc.list`
from those row names and passes it to `gl.drop.loc`, which matches
nothing. `dartR.popgen/R/gl.select.panel.R:141,148` has the same
dependency through `names(sort(rowSums(...)))`. Duplicate locus names
are permitted by genlight and `utils.dart2genlight.r` contains no
uniqueness check (verified: zero occurrences of `duplicated`).
Proposed change: build the frame with explicit, order-preserving row
keys — `data.frame(alf1 = 1 - alf, alf2 = alf, row.names = NULL)` plus
`rownames(out) <- make.unique(locNames(x))`, or return the locus name as
a column. Either keeps positional correspondence and stops the silent
fallback.

**F2 [MEDIUM, confidence: high] — no datatype gate: SilicoDArT returns half the presence frequency (DAT7)**
`R/gl.alf.r:22-26` — there is no `utils.check.datatype()` call, so Tag
P/A data is accepted and divided by 2 despite `ploidy(testset.gs) == 1`.
Failure scenario: `gl.alf(testset.gs)$alf2` spans 0 to 0.5 and equals
`colMeans/2` (verified); `alf1` never falls below 0.5. A user reading
`alf2` as the presence frequency is out by a factor of 2, and `alf1` is
neither the absence frequency nor a complement of anything meaningful.
This is the same defect class as gl.Ho F1 (SilicoDArT silently returning
meaningless heterozygosities) and gl.allele.freq F2.
Proposed change: `utils.check.datatype(x, accept = "SNP", verbose = 0)`
as the first line — SilicoDArT input becomes a fatal error.
**Consequence: SilicoDArT callers of gl.alf get an error instead of a
number.** The caller inventory (below) shows no SilicoDArT caller across
the eight packages.

**F3 [LOW, confidence: high] — roxygen header below template (DOC1, DOC7)**
`R/gl.alf.r:1-20` — `@return` reads "A simple data.frame with ref
(reference allele), alt (alternate allele)"; the columns are named
`alf1` and `alf2` (verified). There is no `@description` (the second
line of `@title` is doing that job), no `@details`, no `@seealso`, and
`@author` reads "Bernd Gruber (bugs? Post to ...)" with neither the
`Author(s):` nor the `Custodian:` label DOC7 requires.
Failure scenario: a user following the manual writes `out$ref` and gets
`NULL`.
Proposed change: docs-only — split `@title`/`@description`, correct
`@return` to name `alf1`/`alf2` and state the NaN-for-all-NA-loci
behaviour and the SNP-only restriction, add
`Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to ...`, add
`@seealso gl.allele.freq, gl.Ho, gl.He`.

**F4 [LOW, confidence: high] — the deprecation notice is commented out and inert (API1 (proposed rule))**
`R/gl.alf.r:23` — `#cat(warn("Deprecated: Please use
gl.allele.freq(x,simple=TRUE)\n"))` is disabled, while
`gl.allele.freq.r:42` states "gl.alf is deprecated" in a developer note
and `gl.alf`'s own `@examples` demonstrate the replacement side by side.
Failure scenario: the status is ambiguous in three places at once. A
user reading the source sees "deprecated"; a user running the function
sees nothing; a maintainer reading `gl.allele.freq.r` sees a settled
decision that was never enacted. Fifteen call sites across four packages
continue to depend on it.
Proposed change: settle the status — see the Redundancy analysis
recommendation below. Either remove the commented line and document
gl.alf as the supported exact fast path, or enable a real deprecation
with a migration plan for the fifteen call sites.

**F5 [INFO, confidence: high] — full densification on every call (DAT6 (proposed rule))**
`R/gl.alf.r:24` — `as.matrix(x)` materialises the whole genotype matrix.
Verified working on an FBM-backed possums.gl (identical result, row
names preserved), so the function does not refuse FBM input; it pays the
full memory cost instead. Not proposed for change on its own: every
member of this family densifies the same way, and a column-chunked
`colMeans` is a family-level decision, not a per-function one.

## Redundancy analysis

The custodian asked whether gl.alf duplicates the rest of the
allele-frequency family. It does duplicate `gl.allele.freq(simple=TRUE)`
today — but not exactly, and the two are about to diverge.

### (a) Equivalence matrix

Verified empirically on testset.gl (255 loci, 250 individuals, 30
populations), testset.gs, possums.gl and platypus.gl.

| Function / call | Statistic returned | NA policy | Agrees with `gl.alf$alf2`? |
|---|---|---|---|
| `gl.alf(x)` | `colMeans(m, na.rm=TRUE)/2` per locus, whole object; `alf1 = 1 - alf2` | per locus, `na.rm = TRUE`; all-NA locus -> NaN | reference |
| `gl.allele.freq(x, simple=TRUE)` | same quantity, then `round(., 4)`; forces `by='loc'`, `percent=FALSE` | identical — same 3 NaN loci, same positions | **yes to 4 dp, not identical**: 80/255 loci differ, max abs diff 5.0e-05; `alf2_simple == round(alf2_alf, 4)` exactly |
| `gl.allele.freq(x, by='loc')$frequency` | same column before the alf1/alf2 split | identical | yes, same 4 dp relation |
| `gl.allele.freq(x, by='popxloc')$frequency` | per population x locus, percentage, 2 dp | per cell, `na.rm=TRUE`; empty cell -> NA | **different aggregation**: mean of the per-population cells differs from the global frequency by up to 0.095 on testset.gl (unequal population sizes). Matches `gl.alf` applied per `seppop` element to 1e-4 (max diff 4.5e-05) |
| `gl.allele.freq(x, by='pop')$frequency` | per population, averaged over loci | cell NaNs dropped in the aggregate | different statistic |
| `utils.recalc.freqhomref` | `colMeans(m == 0, na.rm=TRUE)` | per locus, `na.rm=TRUE` | genotype-class frequency, not allele frequency |
| `utils.recalc.freqhets` | `colMeans(m == 1, na.rm=TRUE)` | same | same |
| `utils.recalc.freqhomsnp` | `colMeans(m == 2, na.rm=TRUE)` | same | same |
| the three together | `FreqHomRef + FreqHets + FreqHomSnp == 1` (verified, 1e-12) | same | **`FreqHomSnp + FreqHets/2 == gl.alf$alf2` to 1.1e-16** — an exact, unrounded reparameterisation |
| `utils.recalc.maf` | `ifelse(alf > 0.5, 1 - alf, alf)` | inherits gl.alf's; NaN becomes NA through `ifelse` | **calls `gl.alf(x)[, 2]` directly** — same numbers by construction |
| `gl.report.allelerich` | rarefied allelic richness per population x site (hypergeometric, Rcpp) | own `na.rm` handling inside the kernel | **genuinely distinct** — an expected-distinct-alleles count at a standardised sample size, not a frequency |

Missing-data fixture check: on a 5-locus subset with locus 3 blanked,
`gl.alf` returns NaN in both columns for that locus and exact values for
the rest; `gl.allele.freq(simple=TRUE)` returns NaN in the same
position. The two agree cell for cell up to the 4 dp rounding on every
dataset tested, SNP and SilicoDArT alike.

Two shapes where they do **not** agree, both in gl.alf's favour:

- Plain genlight not built by dartR: `gl.alf` returns the frame;
  `gl.allele.freq(simple=TRUE)` fails with "argument is of length zero"
  (its unguarded `loc.metrics.flags$monomorphs` access — F3 in that
  report).
- Duplicate locus names: `gl.alf` degrades silently to integer row keys
  (gl.alf F1 above); `gl.allele.freq(simple=TRUE)` errors with
  "replacement has 5 rows, data has 4" (its `aggregate(. ~ locus)`
  collapses the duplicate).

### (b) What is genuinely distinct

- **`gl.allele.freq(by='pop'|'popxloc')`** — different aggregation unit.
  `by='popxloc'` is the only one of these that carries `sum`, `nobs`,
  `nmissing` and `n` alongside the frequency, and it is what
  `gl.dist.pop.r:184` consumes.
- **The `utils.recalc.freq*` trio** — genotype-class frequencies
  (`m == 0`, `== 1`, `== 2`), used to repopulate DArT's `loc.metrics`
  columns after individuals are removed. They are not duplicate
  frequency machinery in the sense of computing the same number twice;
  they compute a finer decomposition from which the allele frequency
  follows exactly. The redundancy is one level up: `utils.recalc.maf`
  recomputes all three and then calls `gl.alf`, densifying
  `as.matrix(x)` four separate times for one object.
- **`gl.report.allelerich`** — related but distinct: richness, not
  frequency.

### (c) Consumer consistency

Three consumers, three different frequency sources:

| Consumer | Frequency source | NA policy |
|---|---|---|
| `gl.dist.pop.r:184` | `gl.allele.freq(percent=TRUE, by='popxloc')` | per cell `na.rm=TRUE`, 2 dp on the percentage scale |
| `gl.report.heterozygosity.r:508` | `gl.alf` per `seppop` element | per locus `na.rm=TRUE`, exact |
| `gl.tree.nj.r:120` | `apply(as.matrix(x), 2, tapply, pop(x), function(e) mean(e)/2)` | **no `na.rm`** — NA propagates |

The first two agree to 1e-4 (verified on a 3-population subset: 42 NaN
cells either way, max abs diff 4.5e-05). `gl.tree.nj` disagrees: on the
same subset it produces 192 NA cells out of 765 against 42 for the other
two, because a single missing genotype voids the whole population x
locus cell. That is a defect in `gl.tree.nj`, already covered by PR #370
and its report — recorded here as the inconsistency requested under (d),
not raised as a new finding under the scope rule.

### (d) Cost

`gl.allele.freq(simple=TRUE)` runs `seppop()` and builds the full
population x locus table before collapsing it, so its cost scales with
population count while gl.alf's does not.

| Dataset | `gl.alf` | `gl.allele.freq(simple=TRUE)` | ratio |
|---|---|---|---|
| testset.gl (250 x 255, 30 pops) | 19.5 ms | 161.5 ms | 8.3x |
| possums.gl (300 x 200, 10 pops) | 15.0 ms | 83.0 ms | 5.5x |
| platypus.gl (81 x 1000, 3 pops) | 6.0 ms | 113.0 ms | 18.8x |
| `lapply` over 30 seppop elements | 0.05 s | 0.67 s | 13.4x |

`gl.alf` is also exact where the replacement rounds, and silent where
the replacement prints 7 lines at default verbosity.

### (e) Caller inventory

Live clones under `D:\workspace\R\` only; dated backup copies
(`dartR.base _2026-04-13` and siblings) excluded.

| Package | Site | Shape | If gl.alf were removed |
|---|---|---|---|
| dartR.base | `utils.recalc.maf.r:91` | `gl.alf(x)[, 2]` | `gl.allele.freq(x, simple=TRUE, verbose=0)[, 2]`. Would introduce 4 dp rounding into the `maf` locus metric and add a `seppop()` pass; `maf` is consumed by `gl.filter.maf` and `gl.report.maf` |
| dartR.base | `gl.report.heterozygosity.r:508` | `gl.alf(y_temp)` on `seppop` elements, then reads `rownames` | direct swap works numerically; inherits gl.alf F1's row-name hazard either way, and adds 30x the `seppop` overhead in the population loop |
| dartR.base | `utils.jackknife.R:43` (`@examples`) | `FUN = "gl.alf"` — function passed by name | needs a named wrapper; `gl.allele.freq` cannot be passed as a bare one-argument function |
| dartR.captive | `utils.assignment.r:68` and `_2`, `_3`, `_4` variants, all line 68 | `gl.alf(y)` inside `lapply` over `seppop`, columns read positionally | 4 sites; wrapper `function(p) gl.allele.freq(p, simple=TRUE, verbose=0)` |
| dartR.popgen | `gl.select.panel.R:141,148,189,211,224,236` | 6 sites; `gl.alf(x)[,1]`, `rowSums(cbind(gl.alf(p1p) * !gl.alf(p2p)))`, and `names()` of the result | needs the wrapper at every site, and `names(sort(rowSums(...)))` depends on locus row names surviving |
| dartR.sim | `gl.sim.WF.run.r:307,311` | `lapply(pop_list_freq_temp, gl.alf)` — bare function reference in a simulation loop | needs a wrapper; simulated genlights may lack `loc.metrics.flags`, which `gl.allele.freq` requires (it errors on plain genlight), and the 13x cost lands inside the generation loop |
| dartR.data, dartR.sexlinked, dartR.spatial, dartRstartup, dartRverse | none | — | — |

Total: 15 call sites across 4 packages, plus `@seealso` cross-links from
`gl.He.r:32` and `gl.Ho.r:32` and the `utils.jackknife` example.

### Recommendation

**Keep gl.alf as a documented, guarded fast path; do not deprecate.**
Retire the commented-out deprecation line and the developer note in
`gl.allele.freq.r:42` instead.

Reasoning, in order of weight:

1. **The two functions are about to stop being equivalent.** PR #374's
   F2 fix makes `gl.allele.freq(simple=TRUE)` return the presence
   frequency for SilicoDArT, where it currently returns half. Simulated
   post-fix values are exactly `2 x gl.alf$alf2` (verified). After #374,
   "use `gl.allele.freq(x, simple=TRUE)` instead" is no longer a true
   statement for Tag P/A data, so a deprecation notice pointing there
   would be wrong for one of the two datatypes. Applying gl.alf F2
   (`accept = "SNP"`) resolves this cleanly: gl.alf becomes SNP-only and
   exactly equivalent to the fixed replacement on its whole domain,
   modulo rounding.
2. **The replacement is lossier and 5-19x more expensive**, and the cost
   is in the wrong place — inside `dartR.sim`'s generation loop and
   `dartR.captive`'s per-population `lapply`. `utils.recalc.maf` would
   also inherit 4 dp rounding in a locus metric that `gl.filter.maf`
   thresholds against.
3. **15 call sites across 4 packages** would each need a wrapper,
   because `gl.allele.freq` cannot be passed as a bare one-argument
   function to `lapply` or `utils.jackknife`. That is a four-repository
   coordinated change for no numerical gain.
4. **gl.alf is the only member of the family that accepts a genlight not
   built by dartR** — which is what `dartR.sim` produces.

The alternative worth considering, if the family is to shrink: invert
the dependency. Make `gl.allele.freq`'s `by='loc'` branch call `gl.alf`
for its frequency column (it already computes the identical expression
inline at `gl.allele.freq.r:194`), so there is one implementation of the
allele-frequency formula in the package rather than three. That is a
change to `gl.allele.freq`, outside this review's scope; recorded here
as a note for that function's Phase C.

## Proposed changes

1. Give `gl.alf` order-preserving, unique row keys so duplicate locus
   names cannot silently become integers (F1).
   **Consequence: the returned row names change from `1:n` to
   `make.unique(locNames(x))` for objects with duplicate locus names;
   objects with unique names are unaffected.**
2. Add `utils.check.datatype(x, accept = "SNP", verbose = 0)` (F2).
   **Consequence: SilicoDArT input errors instead of returning half the
   presence frequency. No SilicoDArT caller exists across the eight
   packages.**
3. Docs-only: split `@title`/`@description`, correct `@return` to name
   `alf1`/`alf2` and state the NaN and SNP-only behaviour, bring
   `@author` to DOC7 form, add `@seealso` (F3).
4. Settle the deprecation status per the recommendation above — remove
   the commented-out deprecation line in `gl.alf.r:23` and the "gl.alf is
   deprecated" developer note in `gl.allele.freq.r:42`, and add a
   `@details` paragraph stating gl.alf's role as the exact SNP fast path
   relative to `gl.allele.freq(simple=TRUE)` (F4).
   *(The second half touches `gl.allele.freq.r`; under the scope rule it
   is recorded here and applied only if the custodian folds it into that
   function's Phase C.)*

## Report notes

Defects observed in other functions during this review. Recorded only —
not fixed, per the scope rule.

- `gl.report.heterozygosity.r:508-511` reads `rownames(gl.alf(...))` as
  locus names with no uniqueness check of its own, and returns negative
  `monoLoc` counts when they are not unique. Fixing gl.alf F1 removes
  the trigger but the arithmetic
  (`nLoc(hold) - nLoc(y_temp) - length(loc.list_NA)`) is unguarded
  against overlap between the two removal lists.
- `utils.recalc.maf.r:87-92` densifies `as.matrix(x)` four times for one
  object (once in each of the three `utils.recalc.freq*` calls, once in
  `gl.alf`). One densification would serve all four.
- `gl.tree.nj.r:120` uses a stricter NA policy than the rest of the
  family (no `na.rm`), producing 4.6x as many missing frequency cells.
  Already covered by PR #370.
- `gl.He.r:32` and `gl.Ho.r:32` cross-link `gl.alf` in `@seealso`. Both
  are still present in `dartR.base/R/` at ed99203 (PR #273 open); if
  they move to dartR.sim, those links become cross-package.

## Coverage

- Standards walk: FS, DOC, VRB (n/a — no `verbose` parameter; function
  is silent by design, verified `capture.output` length 0), DAT, DEP
  (n/a — base R only), PLT (n/a — no plot), STY — run.
- Spec, value-level: hand computation `colMeans(as.matrix(x), na.rm =
  TRUE)/2` against every locus on testset.gl, possums.gl and
  platypus.gl; `alf1 + alf2 == 1`; pinned anchors on testset.gl and
  possums.gl — run.
- Edge cases: all-NA loci (NaN, verified in a constructed fixture and in
  testset.gl's three natural cases); single locus; single individual;
  duplicate locus names; plain non-dartR genlight — run.
- Zero loci / zero individuals: NOT REACHABLE — the dartR `[` method
  errors before gl.alf is called ("Subsetting resulted in zero loci" /
  "zero individuals"). A directly constructed empty genlight fails
  inside `colMeans` with "'x' must be numeric" — an uninformative but
  unreachable-in-practice path; not raised as a finding.
- FBM path (DAT6): run — `gl.gen2fbm(possums.gl)` then `gl.alf`, result
  identical to the in-memory object, row names preserved.
- Cross-function equivalence: gl.alf vs `gl.allele.freq(simple=TRUE)`,
  `by='loc'`, `by='popxloc'`, `by='pop'`; the three `utils.recalc.freq*`
  decomposition; `utils.recalc.maf`; consumer NA policies for
  `gl.dist.pop`, `gl.report.heterozygosity`, `gl.tree.nj` — run on
  testset.gl and testset.gs.
- Timing: 20 repetitions per function on testset.gl and possums.gl, 10
  on platypus.gl, plus a 30-element `lapply` — run.
- Google Group search: SKIPPED — not available, no browser session.

## Approval

Approved 2026-09-07 by Arthur Georges, including the redundancy verdict.

**Redundancy verdict: keep gl.alf as a documented fast path; do not
deprecate.** The retirement of the deprecation intent is a decision, not
an omission: after PR #374 the nominated replacement returns twice
gl.alf's value on SilicoDArT, it is 5-19x slower, it rounds to 4 dp, and
it rejects inputs gl.alf accepts. The commented-out deprecation line in
`gl.alf.r` and the matching developer note in `gl.allele.freq.r` are
removed.

| Change | Decision | By | Note |
|---|---|---|---|
| 1 (F1 HIGH) | approved | Arthur Georges 2026-09-07 | Preserve the original locus names on the returned rows so duplicate names flow through honestly instead of silently reindexing to 1:n |
| 2 (F2 MEDIUM) | approved | Arthur Georges 2026-09-07 | Add the `accept = "SNP"` gate; SilicoDArT becomes a fatal error. All 15 callers verified SNP-context |
| 3 (F3 LOW) | approved | Arthur Georges 2026-09-07 | `@description`/`@details`/`@seealso`, DOC7 author line, `@return` corrected to alf1/alf2 |
| 4 (F4 LOW) | approved | Arthur Georges 2026-09-07 | Deprecation intent retired per the verdict above; both the commented line and the developer note deleted, `@details` documents the fast-path role and cross-references gl.allele.freq |
| F5 INFO | noted | — | Family-level decision; not actioned |

## Outcome

Applied on `review-gl.alf` from upstream/dev at ddaed27. All four
approved changes are in; F5 remains noted only.

**Scope note.** The gate is a deliberate departure from the
preamble-free accessor style ratified for gl.Ho/gl.He (PR #273). It is
the datatype check alone, in the gl.Ho/gl.He form
(`utils.check.datatype(x, accept = "SNP", verbose = 0)`); no
`gl.check.verbosity`, no `utils.flag.start`, no history append. The
function stays silent and returns visibly.

### Verification

(a) **Numbers unchanged for SNP input.** `gl.alf` output is identical
cell for cell, NaN positions included, on testset.gl (255 loci),
possums.gl (200) and platypus.gl (1000): maximum absolute difference 0,
row names identical. `alf1 + alf2 - 1` is exactly 0 on all three. The
all-NA loci still return NaN in both columns (3 loci in testset.gl, 6 in
platypus.gl, 0 in possums.gl, unchanged).

(b) **F1 — duplicate locus names.** On the 30-locus testset.gl fixture
with `locNames(d)[5] <- locNames(d)[4]`, the returned row names are now
the original locus names with the duplicate disambiguated
(`100049805-56-T/A`, `100049805-56-T/A.1`), not `1:30`; the values are
unchanged and positional. `gl.report.heterozygosity(method = "pop")` on
that fixture returned `polyLoc = 30` for all 30 populations and
`monoLoc` of -1 or -2 for 27 of them before the change. It now returns
`monoLoc` 23-27 and `polyLoc` 3-6, against 22-26 and 3-7 for the same
object with unique names — no negative counts remain.

(c) **F2 — SilicoDArT.** `gl.alf(testset.gs)` now stops with "Fatal
Error: inappropriate object passed to function, found SilicoDArT
expecting SNP". Previously it returned `alf2` spanning 0 to 0.5.

(d) **Timing — fast-path advantage retained.** Median ms per call, 20
reps (10 for platypus), same machine and session protocol as Phase A:

| Dataset | gl.alf before | gl.alf after | gl.allele.freq(simple=TRUE) | ratio after |
|---|---|---|---|---|
| testset.gl (250 x 255, 30 pops) | 21.0 | 19.5 | 158.0 | 8.1x |
| possums.gl (300 x 200, 10 pops) | 15.0 | 14.5 | 80.5 | 5.6x |
| platypus.gl (81 x 1000, 3 pops) | 8.0 | 8.0 | 113.0 | 14.1x |
| `lapply` over 30 seppop elements | 36.0 | 47.0 | 613.3 | 12.8x |

The gate costs 0.35 ms per call measured on its own: 1.8% of a
whole-object call, but 15% of a 8-individual per-population call, which
is the 36 -> 47 ms movement in the `lapply` row. The advantage over the
replacement is retained everywhere (5.6x to 14.1x per call, 12.8x over
the 30-population loop). Recorded, not treated as material erosion.

(e) **FBM path.** `gl.gen2fbm(possums.gl)` then `gl.alf`: values
identical to the before state and to the in-memory object, row names
preserved as the locus names.

(f) **Baseline rerun.** 20 tests / 64 assertions passed at the reviewed
state before the change. After: 21 tests / 68 assertions, 0 failures.
Every flip carries an `# [approved Fn]` comment — the duplicate-name row
keys and the gl.report.heterozygosity consequence (F1), the SilicoDArT
acceptance and the SilicoDArT half of the cross-function equivalence
pin (F2), and the `@return` column-name comment (F3). No unexplained
diff. The three sibling test files that touch this code path
(test-gl.Ho.He, test-gl.filter.heterozygosity,
test-gl.report.heterozygosity) pass unchanged: 13, 19 and 21
assertions.

**Caller check.** All 15 sites re-run in their exact call shapes against
the patched function, on the live clones only: dartR.base
`utils.recalc.maf.r:91` (maf still the exact `ifelse` transform of
alf2), `gl.report.heterozygosity.r:508` (both methods), the
`utils.jackknife` example with `FUN = "gl.alf"`; dartR.captive
`utils.assignment{,_2,_3,_4}.r:68` (the `lapply(seppop, ...)` frame,
columns read positionally); dartR.popgen `gl.select.panel.R:141,148`
(`names(sort(rowSums(...)))` returns locus names), `:189`, `:211`,
`:224`, `:236`; dartR.sim `gl.sim.WF.run.r:307,311` (bare function
reference in `lapply`, including the `gl.keep.loc` subset path and a
plain genlight with no locus metrics). No site breaks. Passing `gl.alf`
bare still works because the gate takes only the object.

```json
{
  "function": "gl.alf",
  "package": "dartR.base",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "ddaed27",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT3", "status": "applied", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DAT7", "status": "applied", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "API1", "status": "applied", "change": 4},
    {"id": "F5", "severity": "INFO", "confidence": "high", "rule": "DAT6", "status": "noted", "change": null}
  ],
  "redundancy": {
    "equivalent_to": "gl.allele.freq(x, simple=TRUE) to 4 dp on SNP data; 80/255 loci differ on testset.gl, max abs diff 5.0e-05",
    "exact_reparameterisation": "FreqHomSnp + FreqHets/2 == alf2 to 1.1e-16",
    "distinct": ["gl.allele.freq by='pop'", "gl.allele.freq by='popxloc'", "utils.recalc.freqhomref/freqhets/freqhomsnp", "gl.report.allelerich"],
    "speed_ratio_range": [5.5, 18.8],
    "callers": 15,
    "caller_packages": ["dartR.base", "dartR.captive", "dartR.popgen", "dartR.sim"],
    "recommendation": "keep_as_documented_fast_path",
    "verdict": "kept_as_documented_fast_path",
    "approved_by": "Arthur Georges",
    "approved_date": "2026-09-07"
  },
  "coverage_skipped": ["Google Group: no browser session", "zero loci / zero individuals: not reachable through the dartR [ method"],
  "baseline_test": "tests/testthat/test-gl.alf.R",
  "baseline_before": {"tests": 20, "assertions": 64, "failures": 0},
  "baseline_after": {"tests": 21, "assertions": 68, "failures": 0},
  "timing_after_ms": {"testset.gl": 19.5, "possums.gl": 14.5, "platypus.gl": 8.0, "lapply30": 47.0, "gate_overhead": 0.35},
  "status": "pr-open",
  "pr": null
}
```
