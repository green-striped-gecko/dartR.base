# Review: gl.recalc.metrics (dartR.base)

## Provenance

- Model: claude-opus-5 (Claude Code); Skill: dartr-function-review v2.0.0;
  Base: `upstream/dev` at `ddaed27`. `git diff upstream/dev -- R/gl.recalc.metrics.r`
  is empty, so the reviewed file is the upstream/dev state verbatim.
- Phase C (apply) run 2026-09-08 by claude-opus-5 (Claude Code) on branch
  `review-gl.recalc.metrics`, base `ddaed27`.
- Runtime state: the working tree is `integration-local` at `ed99203`, which
  already carries the merged recalc-battery fixes (PRs #302-#308) in the six
  helper files this function calls. Behaviour recorded here is therefore
  post-helper-fix behaviour; where a helper defect is still visible it is
  marked `propagates-from #<PR>`.
- Family mode: MODIFY (recalculates locus metrics, returns a changed object).
- Datasets: `testset.gl` (250x255, plus a 159-individual 10-population-dropped
  subset), `testset.gs` (218x255, plus a 5-population-dropped subset),
  `platypus.gl` (81x1000), `testset2.gl` (274x755), `possums.gl` (column set
  inspected only). dartR.data 1.2.5.
- Baseline: `tests/testthat/test-gl.recalc.metrics.R` (new file, 87 assertions,
  all passing against the reviewed state).
- Custodian: Luis Mijangos.
- Checks skipped: zero-locus fixture (the `[` method refuses to produce one:
  "Subsetting resulted in zero loci"); Google Group not searched (no browser
  session); dartr2shiny not present in the workspace.

## Verdicts

**Standards: Needs work** — no precondition guard on the `loc.metrics` slot, so
a non-conforming object is silently corrupted rather than rejected; `mono.rm`
is unvalidated; `verbose` is passed unfiltered to six helpers; `@family` is
swallowed into `@title`.

**Spec: Needs work** — every metric the function computes is numerically exact
against independent hand computation on all four SNP fixtures and on
SilicoDArT, and datatype dispatch is correct. The defect is completeness: the
function refreshes `OneRatioRef`/`OneRatioSnp` but leaves `rdepth` — which is
defined from them — at its pre-subset value, so the object it returns is
internally inconsistent for 897 of 1000 loci on `platypus.gl`.

What works well: all ten SNP metrics and all three SilicoDArT metrics match a
hand computation from `as.matrix(x)` to machine precision; genotypes,
individuals, populations and ploidy are byte-identical in and out; `verbose = 0`
is silent on both streams and in both `mono.rm` branches; the return is by
value; repeated calls are idempotent in the metrics.

### Metric correctness matrix

Recalculated and verified exact (max absolute difference 0.000e+00, no NA
mismatches) on the `testset.gl` subset, `platypus.gl` and `testset2.gl`:

| Metric | Datatype | Recalculated | Verified against |
|---|---|---|---|
| `CallRate` | SNP + silico | yes | `signif(1 - colSums(is.na(t))/nInd, 6)` |
| `OneRatioRef` | SNP | yes | `(c0+c1)/ctot` |
| `OneRatioSnp` | SNP | yes | `(c1+c2)/ctot` |
| `PICRef` | SNP | yes | `1-(p^2+(1-p)^2)` |
| `PICSnp` | SNP | yes | `1-(p^2+(1-p)^2)` |
| `AvgPIC` | SNP | yes | `(PICRef+PICSnp)/2` |
| `FreqHomRef` | SNP | yes | `c0/ctot` |
| `FreqHomSnp` | SNP | yes | `c2/ctot` |
| `FreqHets` | SNP | yes | `c1/ctot` |
| `maf` | SNP | yes (undocumented) | `pmin(alf, 1-alf)` |
| `OneRatio` | silico | yes | `colMeans(t == 1, na.rm = TRUE)` |
| `PIC` | silico | yes | `1-(OneRatio^2+(1-OneRatio)^2)` |

Not recalculated:

| Metric | Datatype | Individual-dependent | Status |
|---|---|---|---|
| `rdepth` | SNP | yes | **stale and inconsistent (F1)** |
| `AvgCountRef`, `AvgCountSnp` | SNP | yes | stale; not recoverable from genotypes, but undocumented (F9) |
| `AvgReadDepth`, `StDevReadDepth` | silico | yes | **stale (F1)** |
| `RepAvg`, `TrimmedSequence`, `AlleleSequence`, `SNP`, `SnpPosition`, `clone`, `uid` | SNP | no | correctly untouched |
| `Qpmr`, `Reproducibility` | silico | no | correctly untouched |
| `OneRatio`, `PIC`, `monomorphs` columns in a SNP object | SNP | n/a | never refreshed; spurious columns, propagates-from #308 |

## Findings

**F1 [HIGH, confidence: high] — `rdepth` is left inconsistent with the metrics
the function has just refreshed (DAT4)**

`R/gl.recalc.metrics.r:58-66` — the SNP branch calls `utils.recalc.avgpic`,
`utils.recalc.callrate` and `utils.recalc.maf`; nothing touches `rdepth`.
`rdepth` is defined in `gl.read.dart` (`R/gl.read.dart.r:195-210`) as
`OneRatioRef * AvgCountRef + OneRatioSnp * AvgCountSnp`, rounded to 1 decimal.
That formula reproduces the shipped `testset.gl` `rdepth` exactly (maximum
absolute difference 0 over all 255 loci), so the dependency is definitional,
not incidental. Because the function refreshes `OneRatioRef` and
`OneRatioSnp` and not `rdepth`, the returned object contradicts itself.

Measured: on the 159-individual subset of `testset.gl`, 135 of 252 scored loci
carry an `rdepth` that disagrees with its own refreshed metrics (median
absolute drift 0.1 reads, maximum 5.5, 20 loci off by more than 1 read). On
`platypus.gl`, 897 of 1000 loci disagree. The SilicoDArT branch has the same
shape: `AvgReadDepth` is byte-identical to the pre-subset value.

Failure scenario: `gl.drop.pop()` -> `gl.recalc.metrics()` ->
`gl.filter.rdepth(lower = 10, upper = 100)` retains 125 loci using the stale
`rdepth` and 124 using the consistent one, with 5 loci differing between the
two locus sets. The user has explicitly called the function whose purpose is to
prevent exactly this, and gets a silently wrong filter.

Proposed change: recalculate `rdepth` from the refreshed `OneRatio*` and the
surviving `AvgCount*` columns (a `utils.recalc.rdepth` helper mirroring the
`gl.read.dart` formula), and do the same for SilicoDArT `AvgReadDepth`. Where a
metric cannot be recovered from the genotypes (`AvgCountRef`, `AvgCountSnp`,
which need raw read counts), set its flag FALSE and say so in `@details`.

---

**F2 [HIGH, confidence: high] — a missing `loc.metrics` slot silently
fabricates a bogus metrics table (DAT2, DAT5)**

`R/gl.recalc.metrics.r:54-66` — there is no FUNCTION SPECIFIC ERROR CHECKING
section, so an object without `@other$loc.metrics` reaches the helpers
directly. In the helpers, `x@other$loc.metrics` is read with `$`, which
partial-matches to `@other$loc.metrics.flags` when `loc.metrics` is absent
(verified: `identical(e@other$loc.metrics, e@other$loc.metrics.flags)` is
`TRUE`). Two outcomes, neither safe:

- `nLoc(x) > 1`: an opaque base R error, `"replacement has 255 rows, data has
  1"`, raised from `` `$<-.data.frame` `` inside `utils.recalc.avgpic`. No
  dartR message, no indication of the real cause.
- `nLoc(x) == 1`: **no error at all.** The metric vectors are written into the
  flags data frame, which is then stored back as `@other$loc.metrics`. The
  returned object carries a 14-column table whose columns are the flag names
  (`AvgPIC`, `OneRatioRef`, ..., `monomorphs`, `OneRatio`, `PIC`, `allna`), and
  every real DArT metadata column (`AlleleID`, `TrimmedSequence`, `rdepth`,
  `AvgCountRef`, `AvgCountSnp`, `RepAvg`) is gone.

A plain non-dartR genlight is a third variant: `loc.metrics` is created as a
**list**, not a data frame, so `loc.metrics[keep, , drop = FALSE]` in every
downstream `gl.filter.*` cannot subset it row-wise.

Failure scenario: a user reconstructs a genlight by hand, or strips
`loc.metrics` while keeping flags, calls `gl.recalc.metrics` to "fix" the
object, and receives an object that looks repaired and has lost its sequence
and read-depth metadata without a single message.

Proposed change: add an FS5 guard — require `@other$loc.metrics` to be a data
frame with `nrow == nLoc(x)`, otherwise route through `gl.compliance.check()`
or `stop(error(...))` with a clear message (DAT5). Independently, change the
helpers to read `x@other[["loc.metrics"]]` so partial matching cannot fire.

---

**F3 [MEDIUM, confidence: high] — history entries multiply; one call can
append two (FS8)**

`R/gl.recalc.metrics.r:72-81` — the append at line 80-81 is unconditional, and
`gl.filter.monomorphs` (line 73) appends its own entry first. Two consequences:

- `mono.rm = TRUE` adds **two** entries for one user call:
  `gl.filter.monomorphs(x = x, verbose = 0)` then
  `gl.recalc.metrics(x = x, mono.rm = TRUE, verbose = 0)`. The first records an
  internal implementation step the user never invoked.
- The function is called internally at 26 sites across dartR.base, all of which
  inherit an entry with no way to suppress it. `gl.compliance.check` adds two
  entries per call (`gl.recalc.metrics` then `gl.compliance.check`); measured
  1 -> 3 -> 5 over two successive compliance checks on `testset.gl`. This is the
  same leak reported as F9 of the `gl.compliance.check` review (PR #367); the
  root cause is here.

Failure scenario: a pipeline that calls `gl.compliance.check` at each step
produces a history in which half the entries are `gl.recalc.metrics(x = x,
verbose = 0)`, and the user's actual analysis steps are buried. History is the
package's provenance record; it becomes unusable for reconstruction.

Proposed change: (a) restore `x@other$history` to its pre-filter length after
the `gl.filter.monomorphs` call, so `gl.recalc.metrics` contributes exactly one
entry (`gl.filter.monomorphs` already uses this `hold` idiom internally);
(b) give the function a documented way for internal callers to skip the append,
and use it at the 26 internal call sites. See the skill-maintainer note below —
the conventions catalogue has no rule for this.

---

**F4 [MEDIUM, confidence: medium] — the `monomorphs` flag is never re-examined
when `mono.rm = FALSE` (DAT4)**

`R/gl.recalc.metrics.r:72-77` — the flag is set only as a side effect of
`gl.filter.monomorphs`. On the default path the flag is passed through
unchanged. Verified: an object whose flag reads TRUE, and from which
individuals have been dropped, still reads TRUE after the call, with
monomorphic loci present in the data.

Failure scenario: every `utils.recalc.*` helper gates its "Dataset contains
monomorphic loci which will be included in the calculations" warning on
`!isTRUE(flags$monomorphs)`. A stale TRUE suppresses that warning for the rest
of the session, so downstream heterozygosity and diversity reports are computed
over undetected monomorphs with no notice.

Proposed change: on the `mono.rm = FALSE` path, set
`loc.metrics.flags$monomorphs <- FALSE` (the honest "unknown" state after a
recalculation), or run the detection and set the flag to its true value.

---

**F5 [MEDIUM, confidence: high] — `mono.rm = TRUE` on all-monomorphic data
fails with an internal message**

`R/gl.recalc.metrics.r:73` — when every locus is monomorphic,
`gl.filter.monomorphs` removes them all and the `[` method raises
`"Subsetting resulted in zero loci."`. Verified on a 5-locus constant-genotype
subset of `testset.gl`.

Failure scenario: a single-population subset, or a heavily filtered dataset, is
passed with `mono.rm = TRUE`; the user gets an error naming neither
`gl.recalc.metrics` nor monomorphism.

Proposed change: detect the all-monomorphic case before filtering and either
`stop(error(...))` with a dartR message naming the condition, or return the
object unfiltered with a gated warning.

---

**F6 [LOW, confidence: high] — `verbose` is passed unfiltered to six helpers
(VRB1, FS9)**

`R/gl.recalc.metrics.r:59-65` — each helper runs its own `utils.flag.start` and
`Completed:` banner and its own datatype and monomorphs messages. Measured on
the `testset.gl` subset: `verbose = 1` produces 14 lines (seven
Starting/Completed pairs) where VRB1's "1, begin and end" means two;
`verbose = 2` produces 49 lines, including the all-NA warning five times and
the monomorphs warning four times.

Failure scenario: a user at `verbose = 1` cannot tell which function they
called; a user at `verbose = 2` sees the same warning repeated and reasonably
concludes something is looping.

Proposed change: call the helpers with `verbose = 0` (the idiom
`gl.compliance.check` already uses when calling this function), or with
`verbose - 1`, and emit one consolidated message from `gl.recalc.metrics`.

---

**F7 [LOW, confidence: high] — `@family` is swallowed into `@title` (DOC1,
DOC4)**

`R/gl.recalc.metrics.r:4` — `#'  @family environment` is indented by two spaces,
so roxygen treats it as a continuation of `@title` rather than a tag. The
generated `man/gl.recalc.metrics.Rd:4-7` shows the title ending
`"... object @family environment"`, and the function appears in no family
index. `environment` is a valid family in this package (six other files use it).

Failure scenario: the PDF manual and the help page show a title with a stray
roxygen tag; the function is missing from its `@family` cross-reference block.

Proposed change: unindent the tag to `#' @family environment` and re-run
`devtools::document()` in the same change (DOC4).

---

**F8 [LOW, confidence: high] — documentation does not match behaviour (DOC5,
proposed rule; DOC2; DOC7)**

`R/gl.recalc.metrics.r:1-37`:

- `maf` is recalculated (via `utils.recalc.maf`) but is not in the
  `@description`'s list of recalculated metrics. The nine listed metrics are
  all genuinely recalculated.
- "Metrics that remain unaltered are RepAvg and TrimmedSeq as they are
  unaffected by the removal of individuals" is false as stated: `rdepth`,
  `AvgCountRef` and `AvgCountSnp` are also unaltered, and `rdepth` *is*
  affected by the removal of individuals (F1). For SilicoDArT the untouched set
  also includes `AvgReadDepth`, `StDevReadDepth`, `Qpmr` and `Reproducibility`.
- `@param x` says "containing SNP genotypes"; the function also handles
  SilicoDArT, and does so correctly.
- `@param verbose` reads "[default 2 or as specified using gl.set.verbosity]"
  where DOC2 canon is "[default NULL, adopting the global verbosity set by
  gl.set.verbosity(), or 2 if no global is set]".
- `@author` gives `Custodian:` only, with no `Author(s):` line (DOC7, proposed
  rule).
- There is no `@details` section (DOC1).

Failure scenario: a user reads the description, believes `rdepth` is unaffected
by dropping individuals, and filters on it after subsetting.

Proposed change: list `maf`; replace the "unaltered" sentence with the accurate
split (unaffected by individual removal vs. affected but not recoverable);
correct `@param x`; adopt the DOC2 verbose text; add `Author(s):`; add
`@details`.

---

**F9 [LOW, confidence: high] — `mono.rm` is not validated (FS5)**

`R/gl.recalc.metrics.r:39-52` — there is no parameter check. Measured:
`mono.rm = "yes"` gives `"argument is not interpretable as logical"`,
`mono.rm = NA` gives `"missing value where TRUE/FALSE needed"`,
`mono.rm = NULL` gives `"argument is of length zero"`, and `mono.rm = 1` is
silently accepted as TRUE.

Failure scenario: a scripted call passing a character flag fails at line 72
with a base R message that names neither the parameter nor the function.

Proposed change: add the standard logical check in a FUNCTION SPECIFIC ERROR
CHECKING section, erroring with `stop(error(...))`.

---

**F10 [LOW, confidence: medium] — `utils.flag.start` carries the outdated
`build=` argument (FS3)**

`R/gl.recalc.metrics.r:47-49` — `build = "v.2023.2"`. FS3 records `build=` and
`v=` as outdated. Failure scenario: none at run time; the string is stale
metadata that drifts further from reality with each release.

Proposed change: drop `build =` when the file is next touched.

---

**F11 [INFO, confidence: high] — six full densifications of the genotype matrix
per SNP call (DAT6, proposed rule; STY2)**

`R/gl.recalc.metrics.r:59-61` — `utils.recalc.avgpic`,
`utils.recalc.callrate`, `utils.recalc.freqhets`, `utils.recalc.freqhomref`,
`utils.recalc.freqhomsnp` and `gl.alf` each contain exactly one `as.matrix(x)`
call, so one `gl.recalc.metrics` call materialises the full genotype matrix six
times. All six compute functions of the same three column counts
(`c0`, `c1`, `c2`).

FBM handling was checked and works: a `gl.gen2fbm()`-backed object
recalculates without error and its `CallRate` matches the dense result exactly.
The cost is the concern, not correctness. Timing on `platypus.gl` (81x1000):
0.24 s elapsed.

Proposed change: compute `c0`/`c1`/`c2` and the NA counts once in a shared
internal and have the battery consume them. This is a battery-level refactor
touching six files, so it belongs in a separate change from the rest of this
list.

---

### Propagation notes (helper defects, already reported — not re-found here)

- Spurious `monomorphs`, `OneRatio` and `PIC` **columns** in SNP `loc.metrics`
  — propagates-from #308 (`utils.reset.flags`). Visible in this function's
  output: `platypus.gl`, `testset2.gl` and any `gl.drop.pop`-subsetted
  `testset.gl` carry `OneRatio` and `PIC` columns that `gl.recalc.metrics`
  never refreshes in a SNP object, so they stay stale permanently.
- NULL-unsafe `monomorphs` flag check — propagates-from #307, #303, #304, #305,
  #306. Reachable through `gl.recalc.metrics`, which passes objects straight
  through to the helpers.
- SilicoDArT accepted by the `freq*` and `maf` helpers — propagates-from #303,
  #304, #305, #306. **Not** reachable through `gl.recalc.metrics`: dispatch here
  is correct, sending SilicoDArT to `utils.recalc.avgpic` and
  `utils.recalc.callrate` only. Verified — no SNP-only column is created on
  `testset.gs` and the SNP-only flags stay FALSE.
- All-NA loci yield `NaN` metrics (`AvgPIC`, `FreqHets`, `OneRatio*`, `PIC*`)
  and `NA` `maf`, with `CallRate = 0` — propagates-from #307. Present in the
  `testset.gl` subset (3 all-NA loci).

### Verdicts on the three converging campaign leads

1. **History leak from `gl.compliance.check` (PR #367 F9)** — **confirmed, and
   the root cause is here.** `gl.recalc.metrics` appends unconditionally at
   `R/gl.recalc.metrics.r:80-81` with no suppression mechanism, and
   `gl.compliance.check` calls it at `R/gl.compliance.check.r:186`. Each
   compliance call therefore adds two entries (measured 1 -> 3 -> 5). The direct
   call is correct on its own terms — exactly one entry, recording the actual
   call (F3 covers the two-entry `mono.rm = TRUE` variant, which is a second,
   separate leak from `gl.filter.monomorphs`).
2. **`gl.read.dart` discarding a return value (PR #356)** — **return contract
   confirmed value-returning, not by-reference.** Discarding the result of
   `gl.recalc.metrics(x, ...)` leaves the caller's object byte-identical,
   including its history length. Note for the #356 report: the discarded call
   in `gl.read.dart` is `utils.recalc.maf(glout, verbose = 0)` at
   `R/gl.read.dart.r:215`, not `gl.recalc.metrics` — the same defect class, and
   the same by-value semantics make it a real defect.
3. **Recalc-battery helper defects** — reviewed and recorded above under
   Propagation notes. One propagates into this function's output as a permanent
   staleness (the spurious `OneRatio`/`PIC` columns, #308); the SilicoDArT
   admission defects do not, because dispatch here is correct.

## Proposed changes

1. Recalculate `rdepth` (SNP) and `AvgReadDepth` (SilicoDArT) from the
   refreshed metrics, and set FALSE flags for metrics that cannot be recovered
   from genotypes (`AvgCountRef`, `AvgCountSnp`) (F1).
   **Consequence: numerical output changes for `rdepth` on any object whose
   individuals have changed since import; `gl.filter.rdepth` results change
   accordingly.**
2. Add a FUNCTION SPECIFIC ERROR CHECKING section that requires
   `@other$loc.metrics` to be a data frame with one row per locus, and change
   the helpers to index it with `[["loc.metrics"]]` so `$` partial matching
   cannot reach `loc.metrics.flags` (F2).
3. Restore the pre-filter history length after the `gl.filter.monomorphs` call
   so one user call appends one entry, and add a documented suppression path
   for internal callers, applying it at the 26 in-package call sites (F3).
   **Consequence: history contents change for every function that calls
   `gl.recalc.metrics` internally, including `gl.compliance.check`.**
4. Set `loc.metrics.flags$monomorphs <- FALSE` on the `mono.rm = FALSE` path so
   the flag stops asserting a check that was not performed (F4).
5. Guard the all-monomorphic case before `gl.filter.monomorphs` and fail (or
   pass through) with a dartR message (F5).
6. Call the six helpers with `verbose = 0` and emit one consolidated progress
   message from `gl.recalc.metrics` (F6).
7. Unindent `@family` so roxygen reads it as a tag, and re-run
   `devtools::document()` (F7).
8. Correct the roxygen header: list `maf`; replace the "unaltered metrics"
   sentence with the accurate split; widen `@param x` to both datatypes; adopt
   the DOC2 verbose wording; add an `Author(s):` line and an `@details`
   section (F8).
9. Validate `mono.rm` as a length-1 logical with `stop(error(...))` (F9).
10. Drop the `build =` argument from `utils.flag.start` (F10).
11. (Separate, battery-level) Compute the genotype-count summaries once and
    share them across the six helpers instead of densifying six times (F11).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run. PLT not applicable
  (no plot bundle). DEP not applicable (no Suggests dependency).
- Metric correctness against independent hand computation from `as.matrix(x)`:
  run on the `testset.gl` 10-population-dropped subset, `platypus.gl`,
  `testset2.gl` (SNP) and the `testset.gs` 5-population-dropped subset
  (SilicoDArT). All 12 recalculated metrics exact.
- Completeness / staleness matrix after a subsetting operation: run;
  `rdepth`, `AvgCountRef`, `AvgCountSnp`, `AvgReadDepth`, `StDevReadDepth`
  identified as not refreshed (F1).
- `mono.rm` path, including DAT2 row tracking after locus removal and the
  documented removal of all-NA loci: run. The all-NA removal claim holds
  (3 all-NA loci -> 0).
- Flags (DAT4) before and after, both datatypes: run.
- History (FS8): direct call, `mono.rm = TRUE`, nested in
  `gl.compliance.check`, repeated calls, and call-recording accuracy: run.
- Dispatch (DAT7): SNP vs SilicoDArT on `testset.gs`: run. Dispatch correct.
- Verbosity: `verbose = 0` silence on stdout and stderr in both `mono.rm`
  branches and both datatypes; line counts at 1, 2, 3, 5: run. No plot bundle,
  so the VRB5 graphics half is not applicable.
- Input/output identity: genotypes, individual names, populations, ploidy,
  `ind.metrics`: run, all byte-identical.
- Return contract (by value): run.
- Edge cases: absent `loc.metrics`, absent `loc.metrics.flags`, missing
  individual metric columns, all-NA loci, one-locus object, monomorphic-only
  data, plain non-dartR genlight, invalid `mono.rm`, out-of-range `verbose`:
  run.
- FBM path (DAT6): run — `gl.gen2fbm()` object recalculates and matches the
  dense result; densification count taken from source (F11).
- Zero-locus object: SKIPPED — the `[` method refuses to construct one
  ("Subsetting resulted in zero loci"), so the fixture does not exist.
- Google Group search for user reports: SKIPPED — no browser session.
- dartr2shiny signature sweep (API3): SKIPPED — not present in the workspace.
- `possums.gl`: column set inspected only; not used for numeric pinning
  (simulated data with no `rdepth` or `AvgCount*` columns).

## Approval (Phase B)

| Change | Findings | Decision | By | Note |
|---|---|---|---|---|
| 1 | F1 | **deferred (custodian, 2026-09-08)** | Arthur Georges, 2026-09-08 | `rdepth` and silico `AvgReadDepth` stay stale. No recomputation, no warning, no flag change -- deferred means untouched. The evidence stays in this report (897/1000 loci on `platypus.gl`; the `gl.filter.rdepth` 125-vs-124 divergence) so the finding can be revisited. |
| 2 | F2 | approved | Arthur Georges, 2026-09-08 | Guard the missing `loc.metrics` slot so `$` can never partial-match to `loc.metrics.flags`. Implemented as documented construction of a proper data-frame `loc.metrics` (see Outcome), with a fatal error for a table whose rows do not track the loci. |
| 3 | F3 | approved | Arthur Georges, 2026-09-08 | History discipline: suppress the append on internal calls; one entry for a direct call; fix the `mono.rm = TRUE` double entry. Resolves the deferred F9 of the `gl.compliance.check` review (PR #367) at its root. |
| 4 | F4 | approved | Arthur Georges, 2026-09-08 | Re-examine the `monomorphs` flag when `mono.rm = FALSE`. |
| 5 | F5 | approved | Arthur Georges, 2026-09-08 | All-monomorphic input with `mono.rm = TRUE` must complete or fail informatively; reuse the PR #367 guard shape. |
| 6 | F6 | approved | Arthur Georges, 2026-09-08 | Filter `verbose` before passing it to the six helpers. |
| 7 | F7 | approved | Arthur Georges, 2026-09-08 | Separate `@family` from `@title`. |
| 8 | F8 | approved | Arthur Georges, 2026-09-08 | Documentation. |
| 9 | F9 | approved | Arthur Georges, 2026-09-08 | Validate `mono.rm`. |
| 10 | F10 | approved | Arthur Georges, 2026-09-08 | Drop `build =`. |
| 11 | F11 (INFO) | no action | Arthur Georges, 2026-09-08 | Battery-level refactor, out of scope here. |

## Outcome (Phase C)

Applied 2026-09-08 on branch `review-gl.recalc.metrics` (base `ddaed27`).
Nine findings applied (F2-F10); F1 deferred and untouched; F11 no action.
Files touched: `R/gl.recalc.metrics.r`, `man/gl.recalc.metrics.Rd` and the
five `environment`-family `.Rd` files that gain the cross-reference F7
restores, `NEWS.md`, `tests/testthat/test-gl.recalc.metrics.R`, this report.

- **F2** — a FUNCTION SPECIFIC ERROR CHECKING section now reads the slot as
  `x@other[["loc.metrics"]]`, so exact matching decides. If the slot is absent
  or is not a data frame, a data frame with one row per locus is created to
  receive the metrics and a `verbose >= 1` message states that the DArT
  metadata is not present and cannot be reconstructed; the helpers then write
  into that table instead of into the flags table. If the slot is a data frame
  whose row count does not match `nLoc(x)`, the object is not repairable here
  and the function stops with a message naming the condition. A second
  precondition on the same path is seeded: the helpers read
  `loc.metrics.flags$monomorphs` before testing it, so a missing flags list is
  given `monomorphs = FALSE` (the "not checked" state) before they run.
- **F3** — the append is skipped when any frame above this one belongs to a
  `dartR*` namespace, so a call made as an implementation step of another
  dartRverse function records nothing. The `mono.rm = TRUE` branch restores the
  pre-filter history, so `gl.filter.monomorphs`' entry is not inherited. A
  direct call still appends exactly one entry, recording the user's call.
  **This resolves the deferred F9 of the `gl.compliance.check` review
  (PR #367) at its root**: no change to `gl.compliance.check` is needed.
- **F4** — on the `mono.rm = FALSE` path the flag is set from a check made on
  the metrics just recalculated: a locus is monomorphic if `CallRate == 0`, or
  `FreqHomRef == 1`, or `FreqHomSnp == 1` (SNP), or `OneRatio` is 0 or 1
  (SilicoDArT). That is the `gl.filter.monomorphs` definition read from the
  metrics rather than from a second pass over the genotypes.
- **F5** — the `gl.filter.monomorphs` call is wrapped in `tryCatch`, the same
  shape PR #367 uses; the "Subsetting resulted in zero loci" error is read as
  the all-monomorphic case, the flag is set FALSE, a `verbose >= 1` warning
  names the condition, and the function completes.
- **F6** — the six helpers are called with `verbose = 0`; this function emits
  one consolidated set of messages, including the all-NA count and the
  monomorph count that the helpers previously repeated four and five times.
- **F7** — `@family environment` unindented; documented in the same change.
- **F8** — `maf` listed; the "unaltered metrics" sentence replaced with the
  accurate split (unaffected by individual removal vs. affected but not
  recoverable from genotypes, which names `rdepth`, `AvgCountRef`,
  `AvgCountSnp`, `AvgReadDepth`, `StDevReadDepth`); `@param x` widened to both
  datatypes; DOC2 `verbose` wording; `Author(s):` line; `@details` added.
  Stating that `rdepth` is not recalculated and can be inconsistent is a
  documentation correction, not an F1 fix: no behaviour, warning or flag
  attached to `rdepth` was changed.
- **F9** — `mono.rm` must be a single non-NA logical; `stop(error(...))`.
- **F10** — `build = "v.2023.2"` dropped from `utils.flag.start`.

Verification (R 4.4.2, `pdf(NULL)`, dartR.data 1.2.5), run in two
environments: the branch base `ddaed27`, and an integration preview carrying
the merged helper battery at `ed99203` (the runtime tree the review was
recorded against). Results identical in both unless stated.

- **(a) Metric values unchanged.** Every metric the function computes still
  matches an independent hand computation from `as.matrix(x)`:
  max |difference| = 0.000e+00 over the ten SNP metrics on the `testset.gl`
  159-individual subset (255 loci), `platypus.gl` (1000 loci) and
  `testset2.gl` (755 loci), and over the three SilicoDArT metrics on the
  `testset.gs` subset (255 loci). The SilicoDArT column set is unchanged, so
  no SNP-only column is invented.
- **(b) F2.** `nLoc == 1` with no `loc.metrics`: returns a 1-row data frame
  carrying the ten recalculated metrics, not the flags table
  (`identical(loc.metrics, loc.metrics.flags)` is FALSE, `monomorphs` is not a
  column); the DArT metadata is absent because the input did not carry it, and
  the `verbose >= 1` message says so. `nLoc > 1` with no `loc.metrics`:
  completes with a 255-row table and the message "no locus metrics data frame
  found; creating one with 255 rows", in place of "replacement has 255 rows,
  data has 1". A plain `new("genlight", ...)`: `loc.metrics` is a data frame
  with `nrow == nLoc`, row-subsettable as every `gl.filter.*` requires (it was
  a list before). A 10-row table on a 255-locus object: fatal, "the locus
  metrics data frame has 10 rows for 255 loci ... must track loci one for
  one".
- **(c) F3.** Direct call: history 2 -> 3, one entry,
  `gl.recalc.metrics(x = sub_snp, verbose = 0)`. `mono.rm = TRUE`: 2 -> 3 (was
  2 -> 4), and `gl.filter.monomorphs` appears nowhere in the history. Nested:
  three successive `gl.compliance.check` calls on `testset.gl` now give
  1 -> 2 -> 3 -> 4, where the measured baseline was 1 -> 3 -> 5 -> 7; the
  entries read `gl.read.dart | gl.compliance.check | gl.compliance.check |
  gl.compliance.check`. `gl.drop.pop(recalc = TRUE)`: 1 -> 2, entries
  `gl.read.dart | gl.drop.pop`.
- **(d) F4.** An object whose flag reads TRUE and which holds 160 monomorphic
  or all-NA loci now returns the flag FALSE; the same object after
  `gl.filter.monomorphs` returns TRUE. The metric-derived count agrees exactly
  with the `gl.filter.monomorphs` definition applied to the genotype matrix on
  all four fixtures: 160/160 (`testset.gl` subset), 359/359 (`platypus.gl`),
  144/144 (`testset2.gl`), 61/61 (`testset.gs` subset).
- **(e) F5.** A 5-locus all-monomorphic object with `mono.rm = TRUE`
  completes, keeps its 5 loci, sets `monomorphs = FALSE`, prints "Warning: all
  5 loci are monomorphic or scored all NA; none removed" at `verbose >= 1` and
  nothing at `verbose = 0`. No "Subsetting resulted in zero loci".
- **(f) F1 untouched.** `rdepth` is byte-identical in and out; 135 of 252
  scored loci on the subset and 897 of 994 on `platypus.gl` still disagree
  with their own refreshed metrics; silico `AvgReadDepth` is byte-identical;
  `gl.filter.rdepth(lower = 10, upper = 100)` still retains 125 loci on the
  subset. Both `[pins defect F1]` assertions pass unchanged.
- **(g) Identity, silence, baseline.** Genotypes, individual names,
  populations, ploidy and `ind.metrics` are identical in and out.
  `verbose = 0` emits zero lines on stdout and stderr in every branch tested,
  including the new guard, all-monomorphic and `mono.rm = TRUE` paths;
  `verbose = 1` emits 2 lines (was 14); `verbose = 2/3/5` emit 8. Metrics are
  idempotent across repeated calls. The FBM path still works
  (`gl.gen2fbm` object recalculates, `CallRate` matches the dense result).
  Baseline rerun: 26 tests / 105 assertions, 0 failures, 0 errors, in both
  environments. Every flip is annotated `[approved Fn]` and maps to an
  approved finding: F4 (1 pin), F3 (2 pins), F5 (1), F9 (1), F2 (2), F6 (1).
  One test that pinned no defect, "a missing flags slot is tolerated and
  rebuilt", changes from error to pass: at `ddaed27` it was blocked by the
  NULL-unsafe `monomorphs` read in `utils.recalc.avgpic` (open PR #307), which
  the F2 flags precondition now avoids. It passed already in the integration
  preview.
- **Callers.** `gl.recalc.metrics` is called at 26 sites in dartR.base and at
  2 sites in dartR.popgen (`gl.ld.haplotype.r:250`, `gl.select.panel.R:248`);
  no other dartRverse clone calls it. All 28 are internal implementation
  steps that consume the returned metrics, so the only change they see is the
  suppressed history entry. Full dartR.base test suite before and after the
  patch: of 53 test files, the only one whose result changes is
  `test-gl.recalc.metrics.R`. The other 52 are identical line for line, with
  the same 9 failures and 2 errors before and after (pre-existing at
  `ddaed27` in `test-gl.filter.callrate`, `test-gl.filter.hwe`,
  `test-gl.fixed.diff`, `test-gl.report.allelerich`, `test-gl.report.basics`,
  `test-gl.report.callrate` and `test-gl.report.hwe`).
- **Merge-order interaction with PR #367.** The `gl.compliance.check` baseline
  on the `review-gl.compliance.check` branch contains
  `test_that("BUG(F9): one call appends two history entries (internal leak)")`,
  which asserts `added == 2` and that `gl.recalc.metrics` appears in the
  history. That pin is the defect this PR removes at its root, so whichever
  branch merges second must flip it to `added == 1` with no
  `gl.recalc.metrics` entry, and re-tag it as approved rather than pinned.
- **Checks skipped.** `R CMD check` was not run (the branch base has
  pre-existing failures in unrelated test files); the zero-locus fixture is
  still not constructible; Google Group and dartr2shiny sweeps unchanged from
  Phase A.

## Skill-maintainer notes

Two catalogue gaps hit in this review, both already recurring in the campaign:

1. **No rule covers history discipline for functions called internally.** FS8
   says to append when returning a modified object, but says nothing about a
   function that is itself an implementation step of another exported function.
   The campaign has now hit this twice — `gl.compliance.check` F9 (PR #367) and
   `gl.recalc.metrics` F3 — with only FS8 to cite, and FS8 arguably endorses
   the current behaviour. Proposed **FS12**: a function called internally by
   other exported functions must either offer a documented way to suppress its
   history append or have its callers restore the history length; the
   user-visible history records the call the user made, not the implementation
   steps beneath it.
2. **No rule covers completeness of a repair function.** DAT4 obliges the
   function that *invalidates* a metric to flag or recalculate it, but nothing
   obliges the designated repair function (`gl.recalc.metrics`) to leave every
   derived locus metric either refreshed or flagged FALSE. F1 had to be argued
   from DAT4 by extension. Proposed **DAT8**: a function whose contract is to
   restore metric validity must, for every column in `loc.metrics`, either
   recalculate it or set its flag FALSE; a column left both stale and flagged
   TRUE is a finding.

Also, once again: "implemented != documented" (F8) has only the `[proposed]`
DOC5 to cite. This is the fourth review in the campaign to lean on it; it looks
ready for ratification.

```json
{
  "function": "gl.recalc.metrics",
  "package": "dartR.base",
  "family": "modify",
  "skill_version": "2.0.0",
  "commit": "ddaed27",
  "runtime_tree": "ed99203",
  "model": "claude-opus-5",
  "phase_c_model": "claude-opus-5",
  "phase_c_date": "2026-09-08",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1",  "severity": "HIGH",   "confidence": "high",   "rule": "DAT4",  "status": "deferred (custodian, 2026-09-08)", "change": 1},
    {"id": "F2",  "severity": "HIGH",   "confidence": "high",   "rule": "DAT2",  "status": "applied", "change": 2},
    {"id": "F3",  "severity": "MEDIUM", "confidence": "high",   "rule": "FS8",   "status": "applied", "change": 3},
    {"id": "F4",  "severity": "MEDIUM", "confidence": "medium", "rule": "DAT4",  "status": "applied", "change": 4},
    {"id": "F5",  "severity": "MEDIUM", "confidence": "high",   "rule": "FS5",   "status": "applied", "change": 5},
    {"id": "F6",  "severity": "LOW",    "confidence": "high",   "rule": "VRB1",  "status": "applied", "change": 6},
    {"id": "F7",  "severity": "LOW",    "confidence": "high",   "rule": "DOC1",  "status": "applied", "change": 7},
    {"id": "F8",  "severity": "LOW",    "confidence": "high",   "rule": "DOC5",  "status": "applied", "change": 8},
    {"id": "F9",  "severity": "LOW",    "confidence": "high",   "rule": "FS5",   "status": "applied", "change": 9},
    {"id": "F10", "severity": "LOW",    "confidence": "medium", "rule": "FS3",   "status": "applied", "change": 10},
    {"id": "F11", "severity": "INFO",   "confidence": "high",   "rule": "DAT6",  "status": "no-action", "change": 11}
  ],
  "resolves": [
    {"pr": 367, "note": "deferred F9 (history leak from gl.compliance.check) fixed at its root here; #367's BUG(F9) pin must flip in whichever branch merges second"}
  ],
  "propagates_from": [
    {"pr": 308, "note": "spurious monomorphs/OneRatio/PIC columns in SNP loc.metrics, never refreshed here"},
    {"pr": 307, "note": "NULL-unsafe monomorphs check; NaN metrics on all-NA loci"},
    {"pr": 303, "note": "NULL-unsafe monomorphs check (silico admission not reachable via this function)"},
    {"pr": 304, "note": "NULL-unsafe monomorphs check (silico admission not reachable via this function)"},
    {"pr": 305, "note": "NULL-unsafe monomorphs check (silico admission not reachable via this function)"},
    {"pr": 306, "note": "NULL-unsafe monomorphs check (silico admission not reachable via this function)"}
  ],
  "baseline_test": "tests/testthat/test-gl.recalc.metrics.R (26 tests, 105 assertions after the approved flips)",
  "coverage_skipped": [
    "zero-locus fixture: the [ method refuses to construct one",
    "Google Group search: no browser session",
    "dartr2shiny signature sweep: not present in the workspace"
  ],
  "status": "phase-c-complete",
  "pr": null
}
```
