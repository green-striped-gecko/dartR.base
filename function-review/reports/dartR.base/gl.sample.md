# Review: gl.sample (dartR.base)

- Family mode: modify
- Date: 2026-09-08 (Phase A), 2026-09-09 (Phase C)
- Reviewer: Claude Opus 5 (`claude-opus-5`), dartr-function-review v2.0.0
- Phase C applied by: Claude Opus 5 (`claude-opus-5`), on branch `review-gl.sample` cut from
  `upstream/dev` (`ddaed27`) in worktree `D:/workspace/R/wt-base-AA`
- Package commit: `ddaed27` (upstream/dev; `git diff upstream/dev -- R/gl.sample.r` is empty, so
  the reviewed state is upstream/dev unmodified. Working tree is `integration-local` at `ed99203`
  with unrelated uncommitted test edits.)
- Dependency versions: adegenet 2.1.10, dartR.data 1.2.5, R 4.4.2
- Datasets: testset.gl, testset.gs, platypus.gl, plus constructed fixtures (all-NA individual,
  object with an unused population factor level, FBM-backed platypus.gl and testset.gl)
- Baseline: `tests/testthat/test-gl.sample.R` (new file, 80 assertions, all passing against the
  reviewed state)

## Verdict

**Standards: Rework** — the assembly step rebuilds the result with `do.call(rbind, ...)`, which
discards the entire `@other` list. Locus metrics, individual metrics, coordinates, metric flags
and history are all gone from the returned object, so DAT2, DAT3, DAT4 and FS8 fail together.
Meeting the metadata contract means replacing the assembly approach, not patching it. There is no
FS5 validation block, no FS9 completion message, and the roxygen block carries duplicate `@details`
and `@family` tags.

**Spec: Rework** — the function does not deliver what `@return` promises. Any population with
exactly one member is silently dropped from the result and replaced by unrelated individuals drawn
from elsewhere in the object. On `testset.gl` this happens on every call, including the default
call, and the returned individual count is still correct, so nothing signals the substitution.

What works well: seeded reproducibility holds exactly; `verbose = 0` is fully silent; SilicoDArT
dispatch preserves ploidy 1; the input object is not mutated; and resampling individuals with
replacement does *not* trigger the known adegenet SNPbin duplicate-index defect (see Coverage).

## Fixture status (for the skill maintainer)

The task framing described `gl.sample` as a registered golden fixture of this skill, with the
documented defect "`sample(x, n)` returning random draws from `1:x` instead of resampling `x`
itself whenever a group had exactly one member".

Three separate facts, kept apart:

1. **`gl.sample` is not in `fixtures.md`.** The table at
   `C:\Users\arthu\.claude\skills\dartr-function-review\references\fixtures.md` holds five rows —
   `gl.filter.overshoot`, `utils.het.report.r`, `gl.test.heterozygosity`, `gl.read.csv`,
   `gl.filter.hamming`. There is no `gl.sample` row. The fixture entry described in the task does
   not exist in the file.

2. **The defect is STILL PRESENT at the reviewed state.** `R/gl.sample.r:68` at `ddaed27` reads:

   ```r
   samps <- unlist(lapply(ss, function(x) sample(x, nsample, replace=replace)))
   ```

   Confirmed empirically, not by reading — see F1. Line 68 is the only `sample()` call in the
   file, so there is no second instance of the idiom to check.

3. **A correct fix exists but is stranded.** Commit `11e8deb` "Fix singleton-population corruption
   in gl.sample" (Arthur Georges, 2026-08-19) changes line 68 to the safe idiom
   `x[sample.int(length(x), nsample, replace=replace)]`. That commit is on the local-only `dev`
   branch (ahead of `origin/dev` by 110 commits). It is **not** an ancestor of `upstream/dev`
   (`ddaed27`) and **not** an ancestor of the working tree (`ed99203`), verified with
   `git merge-base --is-ancestor`. It has no manifest row and no review report, so it never went
   through the campaign workflow.

**Verdict: still defective.** The fixture is live, not spent — a Phase A review of `ddaed27` must
find this, and this one did. Maintenance items for the skill owner:

- Add the `gl.sample` row to `fixtures.md` if it is meant to be a registered fixture. Suggested
  entry: pre-fix state `ddaed27` (or "parent of `11e8deb`"), defect "`sample(v, n)` on a length-1
  index vector draws from `1:v`, so single-individual populations are replaced by unrelated
  individuals", expected class "Spec axis: sampling correctness", min severity BLOCKER.
- Decide what to do with `11e8deb`. It fixes F1 correctly and completely, but leaves F2 through
  F10 untouched, so merging it alone would close the loudest defect while the returned object
  still arrives with no metadata.
  *(Resolved 2026-09-09: Phase C reused its idiom and credits it in the commit message. The branch
  itself was neither merged nor cherry-picked, so it can be discarded once this PR lands.)*

## Findings

**F1 [BLOCKER, confidence: high] — single-individual populations are silently replaced (spec axis:
sampling correctness; no rule ID — see skill-maintainer note)**

`R/gl.sample.r:68` — the anonymous function's parameter is named `x`, shadowing the genlight `x`,
and receives one population's vector of individual indices. When that vector has length 1,
`sample(x, nsample, replace = replace)` follows R's documented length-1 convention and draws from
`1:x` — the index treated as a range — instead of returning the single index repeated.

Failure scenario: `testset.gl` has two single-member populations, `EmmacNormLeic` (index 185) and
`EmmacNormSalt` (index 178). `gl.sample(testset.gl, nsample = 3, replace = TRUE)` returns 90
individuals spanning 28 populations, not 30. Both singleton populations are absent. Their six
slots are filled by individuals drawn uniformly from `1:185` and `1:178`, which arrive carrying
their own population labels, so five unrelated populations come back over-represented (counts of
4 or 5 where 3 was requested). Over 20 consecutive runs the singleton populations appeared zero
times — this is deterministic, not a rare race. The total individual count is always correct,
which is precisely why the substitution is silent.

The default call is affected too: `nsample` defaults to `min(table(pop(x)))`, which is 1 on
`testset.gl` *because* of those singleton populations, so `gl.sample(testset.gl)` corrupts as
well. Any bootstrap loop over a dataset with a rare population — the documented use case in
`@details` — produces replicates that never contain that population and are quietly enriched for
others.

In isolation: `set.seed(1); sample(c(185), 3, replace = TRUE)` returns `68 167 129`.

Proposed change: `x[sample.int(length(x), nsample, replace = replace)]`, which samples from the
vector regardless of its length. Rename the lambda parameter so it does not shadow `x`.

---

**F2 [BLOCKER, confidence: high] — the returned object has no `@other` metadata at all (DAT2,
DAT3, DAT4, FS8)**

`R/gl.sample.r:72-76` — the result is assembled by splitting the sampled indices into chunks,
subsetting each, and calling `do.call(rbind, px)`. The genlight `rbind` method's SNPbin fallback
drops the `@other` list wholesale.

Failure scenario: `gl.sample(platypus.gl, nsample = 4, replace = TRUE)` returns an object where
`names(x@other)` is empty. `loc.metrics`, `ind.metrics`, `latlon`, `loc.metrics.flags` and
`history` are all `NULL`, where the input carried all five.

Concrete downstream break, run against the reviewed state:

```
gl.filter.callrate(out, threshold = 0.5, plot.display = FALSE, verbose = 0)
#> Error: incorrect number of dimensions
```

`gl.compliance.check()` only partially repairs this. It recomputes `loc.metrics` from the
genotypes, after which `gl.filter.callrate` runs — but `ind.metrics` comes back as a single `id`
column, and the 16 other columns present in `platypus.gl` (`lat`, `lon`, `Sex`, `AgeClass`,
`Weight`, `Microchip`, `PlateID`, and the rest) are gone permanently, along with `@other$latlon`.
Any spatial or sex-linked analysis on a `gl.sample` result loses its coordinates and its sex
assignments with no message.

The mechanism was isolated: `platypus.gl[c(1,2,3), ]` preserves `@other` intact with
`ind.metrics` correctly subset to 3 rows. `rbind()` on that same result empties `@other`. The `[`
method is not at fault; `rbind` is.

Proposed change: drop the split/rbind machinery and subset once with a positive index vector,
`xx <- x[samps, ]`, which already tracks `ind.metrics`, `pop` and `ploidy` correctly and leaves
`loc.metrics` untouched. This is also the fix for F10.

---

**F3 [HIGH, confidence: high] — the same call returns a structurally different object depending on
`options(dartR_fbm)` (DAT6, API1 — proposed rule)**

`R/gl.sample.r:75-76` — the `rbind` metadata loss in F2 occurs only on the SNPbin path. On the
FBM-backed path, `rbind` takes a different branch and `@other` survives with all six slots and
`ind.metrics` correctly in sync at `nrow == nInd`.

Failure scenario: identical code, identical seed, identical dataset. With
`options(dartR_fbm = FALSE)` the result has no metadata; with `options(dartR_fbm = TRUE)` it has
complete metadata. A script that works under one global option fails under the other, and the
option is a session-level setting the caller may not have set themselves. Verified side by side
on `platypus.gl` with `set.seed(1), nsample = 4, replace = TRUE`.

Proposed change: fixing F2 makes both paths take the same subsetting route, which closes this
divergence as a side effect. Add a test asserting both paths return the same `@other` structure.

---

**F4 [HIGH, confidence: high] — locus metrics are left flagged valid after the individuals change
(DAT4)**

`R/gl.sample.r:62-80` — nothing touches `@other$loc.metrics.flags`. Resampling individuals
invalidates every per-locus statistic computed across individuals: `CallRate`, `maf`, `FreqHets`,
`FreqHomRef`, `FreqHomSnp`, `AvgPIC`, `monomorphs`.

Failure scenario: on the FBM path (where the flags survive at all, per F3),
`gl.gen2fbm(testset.gl)` carries `loc.metrics.flags$CallRate == TRUE`. After
`gl.sample(xf, nsample = 6, replace = TRUE, onepop = TRUE)` the flag is still `TRUE`, but the
stored `CallRate` disagrees with the recomputed value for **248 of 255 loci**, with a maximum
absolute difference of 0.9855. Locus 2 is stored as 0.449 against an actual 0.167; locus 1 is
stored as 0.984 against an actual 1. A downstream `gl.report.callrate` or `gl.filter.callrate`
trusts the flag, skips recalculation, and filters on numbers describing the original object rather
than the sample.

Proposed change: set the affected entries of `loc.metrics.flags` to `FALSE` before returning, or
call `gl.recalc.metrics()`. Setting the flags is the cheaper and more conventional option here,
since a bootstrap loop does not want the recalculation cost on every replicate.

---

**F5 [HIGH, confidence: high] — the call is never appended to history (FS8)**

`R/gl.sample.r:77-80` — the function returns a modified genlight, so FS8 requires
`nh <- length(x@other$history); x@other$history[[nh + 1]] <- match.call()`. No such block exists.

Failure scenario: on the SNPbin path the history is destroyed outright (`NULL`, per F2). On the
FBM path it survives but is never extended — `platypus.gl` goes in with 4 entries and comes out
with 4. Either way the provenance record of a bootstrap replicate does not say it was sampled, so
a result cannot be traced back to the resampling step that produced it.

Proposed change: append `match.call()` to the history of the object actually returned (`xx`), on
the single return path, after the metadata is restored per F2.

---

**F6 [MEDIUM, confidence: high] — renamed individuals desync from `ind.metrics$id` (DAT2)**

`R/gl.sample.r:77-79` — `indNames(xx)` is rewritten with a zero-padded ordinal prefix, but
`@other$ind.metrics$id` is not updated to match.

Failure scenario: on the FBM path, `indNames` come back as `1_AA013220, 2_AA004861, ...` while
`ind.metrics$id` still reads `AA013220, AA004861, ...`. dartR convention, enforced by
`gl.compliance.check()`, is that these agree. Code that joins individual metadata by `id` against
`indNames` silently matches nothing. On the SNPbin path the question does not arise only because
`ind.metrics` has been destroyed (F2), so fixing F2 without fixing this converts a hidden problem
into a visible desync.

Proposed change: assign the new names to `@other$ind.metrics$id` in the same step, and note the
renaming in the `@return` text.

---

**F7 [MEDIUM, confidence: high] — no parameter validation; `nsample` failures surface as raw
`sample()` errors and one silent truncation (FS5, API1 — proposed rule)**

`R/gl.sample.r:60-61` — the FUNCTION SPECIFIC ERROR CHECKING section is present as a comment and
empty. `nsample`, `replace` and `onepop` are never validated. Observed on `platypus.gl`:

| Call | Current result |
|---|---|
| `nsample = 0` | `Error: Subsetting resulted in zero individuals.` — does not name `nsample` |
| `nsample = 0.5` | same error; a caller assuming a proportion gets no hint that `nsample` is a count |
| `nsample = 2.7` | **succeeds silently**, returns 2 per population |
| `nsample = 25`, `replace = FALSE` | `Error: cannot take a sample larger than the population when 'replace = FALSE'` — does not name which population (the smallest is 17) |
| `nsample = -1` | `Error: invalid 'size' argument` |
| `nsample = NA` | `Error: vector size cannot be NA/NaN` |

The `nsample = 2.7` row is the one that changes results rather than stopping: a computed
`nsample` such as `nInd(x)/3` truncates without a word.

Proposed change: add an FS5 block validating `nsample` as a positive whole number; error with
`stop(error(...))` naming `nsample` for zero, negative, `NA` and non-integer values; and when
`replace = FALSE`, check `nsample` against the smallest population and name that population in the
error.

---

**F8 [MEDIUM, confidence: high] — no completion message at any verbosity (FS9, VRB1)**

`R/gl.sample.r:80` — `return(xx)` is reached without
`if (verbose > 0) cat(report("Completed:", funname, "\n"))`.

Failure scenario: measured with the result assigned so the object is not printed —
`verbose = 1` emits 1 line ("Starting gl.sample"), `verbose = 2` and `3` emit 4, `verbose = 5`
emits 5, and none contains "Completed". At `verbose = 1`, documented as "begin and end", only the
beginning prints, so a user cannot tell a completed call from one that stopped inside. `verbose`
is otherwise unused by the function's own body — no progress or summary reporting at any level,
which also makes the `verbose` scale claim in `@param` untrue (DOC5).

Proposed change: add the FS9 line; add a `verbose >= 2` message stating individuals drawn per
population and whether replacement was used.

---

**F9 [MEDIUM, confidence: medium] — `onepop = TRUE` changes what `nsample` means and what its
default resolves to, undocumented (DOC5 — proposed rule; API1 — proposed rule)**

`R/gl.sample.r:46, 64, 68` — `nsample` is a **per-population** count when `onepop = FALSE` and a
**whole-object total** when `onepop = TRUE`. Neither `@param nsample` nor `@return` says so.

Separately, the default `nsample = min(table(pop(x)))` is a lazily-evaluated promise. It is not
forced until line 68, which is *after* line 64 has overwritten `pop(x)` with a single level `"A"`.
So under `onepop = TRUE` the default silently resolves to `nInd(x)`.

Failure scenario: on `platypus.gl` (81 individuals, 3 populations, smallest 17),
`gl.sample(x)` returns 51 individuals, while `gl.sample(x, onepop = TRUE)` returns 81 — the whole
object. A caller who adds `onepop = TRUE` to an existing default call expecting "same sample size,
ignore structure" gets a full-size resample instead. The behaviour is arguably the useful one for
`onepop`, but it arrives by accident of evaluation order rather than by design, and would break if
the default were ever forced earlier.

Proposed change: force `nsample` explicitly before line 64, or compute the two defaults
deliberately with a comment; document both meanings in `@param nsample` and `@return`.

---

**F10 [LOW, confidence: high] — the split/rbind chunking serves no purpose and reorders rows on
one branch (STY1, STY3)**

`R/gl.sample.r:72-76` — `ns <- ceiling(length(samps)/nInd(x))` is 1 whenever the total draw does
not exceed `nInd(x)`, which is the common case, so `split()` produces a single chunk and the
machinery is a no-op wrapped around `x[samps, ]`. It does not deduplicate repeated indices (the
maximum duplicate count within a chunk was 5 in a test draw), so whatever it was written to avoid,
it does not avoid it.

When the total draw does exceed `nInd(x)`, `ff <- rep(1:ns, nInd(x))` assigns chunk labels
cyclically, so `split()` interleaves and `rbind()` concatenates the interleaved groups. Verified on
`platypus.gl` with `nsample = 30` (90 draws, 81 individuals, `ns = 2`): the returned row order does
not match the draw order, though the multiset of rows is correct. Draws are exchangeable within a
population, so this is not a statistical error, but the row order becomes a function of whether
`nsample * nPop` happens to exceed `nInd`, and the ordinal prefixes in `indNames` are assigned to
the scrambled order.

Proposed change: delete lines 72-76 and subset once — `xx <- x[samps, ]`. Same fix as F2.

---

**F11 [MEDIUM, confidence: high] — roxygen block has duplicate tags, missing defaults, and a
non-standard `verbose` description (DOC1, DOC2, DOC7 — proposed rule, DOC5 — proposed rule)**

`R/gl.sample.r:1-43`:

- Two `@details` blocks (lines 12 and 15-21) and two `@family` tags (line 3 `data manipulation`,
  line 42 `base dartR`).
- `@param nsample` has no `[default ...]` clause; `@param replace` says "(default)" instead of
  `[default TRUE]`; `@param verbose` reads "set verbosity" instead of the DOC2 standard text.
- `@param onepop`'s `[default FALSE]` is followed by a trailing space.
- `@author` reads `Bernd Gruber (Post to \url{...})` with no `Author(s):` or `Custodian:` labels,
  against DOC7.
- `@return` says "returns a genlight object with nsample samples from each populations" — untrue
  under `onepop = TRUE` (F9), untrue for singleton populations (F1), and silent on the loss of
  `@other` (F2) and the renaming of individuals (F6). `@title` says "Samples individuals from
  populations" while `@description` says "subsample individuals", two terms for one action.
- Tag order does not follow DOC1 (`@name`, `@title`, `@family`, `@description`, `@details`,
  `@param`, `@return`, `@author`, `@examples`, `@export`).

Failure scenario: `?gl.sample` renders duplicated details, lists the function under two families,
and states a return contract the code does not honour. A user reading `@param replace` cannot tell
what the default is without reading the source.

Proposed change: rewrite the block to DOC1 order with a single `@details` and single `@family`,
DOC2 `verbose` text, `[default ...]` on every parameter, DOC7 `@author` structure, and a `@return`
that matches post-fix behaviour. Run `devtools::document()` in the same change (DOC4).

---

**F12 [LOW, confidence: high] — outdated `build=` argument on `utils.flag.start` (FS3)**

`R/gl.sample.r:57` — `utils.flag.start(func = funname, build = "Jody", verbose = verbose)`.
FS3 records `build=` as outdated. The parameter still exists in the signature
(`function(func = NULL, build = NULL, verbose = NULL)`) so nothing breaks, but the value is a
stale release codename.

Failure scenario: no runtime effect. The finding is drift — the argument survives only in
un-updated files and misleads anyone copying this file as a template.

Proposed change: `utils.flag.start(func = funname, verbose = verbose)`.

---

**F13 [LOW, confidence: high] — variable shadowing and commented-out dead code (STY1, STY3)**

`R/gl.sample.r:68` — the lambda parameter is named `x`, shadowing the genlight `x` inside
`lapply`. This is the direct enabler of F1: the shadowing is what makes `sample(x, nsample, ...)`
read as plausible code at a glance.

`R/gl.sample.r:50-52` — three commented-out lines that would null out `loc.metrics` and
`ind.metrics` "to speed up", under the heading "remove metadata". Given F2 they read as an early
draft of the metadata handling and should not be left in a file whose metadata behaviour is being
corrected.

Indentation is inconsistent (lines 56, 72-76 are off the surrounding level).

Proposed change: rename the lambda parameter to `idx`; delete the commented-out block; normalise
indentation on lines touched.

---

**F14 [INFO, confidence: high] — the documented example works around F2**

`R/gl.sample.r:35-36` — the `@examples` block calls `gl.compliance.check(dummy)` immediately after
every `gl.sample()` call. That step is only needed because the returned object is non-compliant
(F2). The workaround is baked into the documentation rather than reported as a defect.

Proposed change: once F2 is fixed, drop the `gl.compliance.check()` line from the example.

## Proposed changes

1. Replace `sample(x, nsample, replace = replace)` with
   `idx[sample.int(length(idx), nsample, replace = replace)]` and rename the shadowing lambda
   parameter (F1, F13). **Consequence: results change for any object containing a population with
   exactly one member — such populations now appear in the output instead of being replaced by
   unrelated individuals. Seeded output changes for those datasets.** Commit `11e8deb` on the
   local `dev` branch already implements this fix.
2. Replace the split/rbind assembly with a single positive-index subset `xx <- x[samps, ]`,
   restoring `@other` on the SNPbin path, aligning the SNPbin and FBM paths, and removing the
   row reordering (F2, F3, F10). **Consequence: the returned object gains `loc.metrics`,
   `ind.metrics`, `latlon`, `loc.metrics.flags` and `history`, which callers previously had to
   rebuild with `gl.compliance.check()`; row order changes when `nsample * nPop > nInd`.**
3. Set the individual-dependent entries of `@other$loc.metrics.flags` to `FALSE` before returning
   (F4). **Consequence: downstream report and filter functions recalculate locus metrics instead
   of using stale stored values, so their numerical output changes.**
4. Append `match.call()` to `@other$history` on the return path (F5).
5. Write the new `indNames` through to `@other$ind.metrics$id` so the two agree (F6).
6. Add an FS5 validation block for `nsample`: reject zero, negative, `NA` and non-integer values
   with `stop(error(...))` naming the argument, and when `replace = FALSE` check against the
   smallest population and name it in the error (F7). **Consequence: `nsample = 2.7` now errors
   instead of silently truncating to 2.**
7. Add the FS9 completion message and a `verbose >= 2` progress line reporting draws per
   population and the replacement setting (F8).
8. Force `nsample` before `pop(x)` is overwritten, and document that `nsample` is per-population
   when `onepop = FALSE` and a whole-object total when `onepop = TRUE` (F9).
9. Rewrite the roxygen block to DOC1 order: single `@details`, single `@family`, DOC2 `verbose`
   text, `[default ...]` on every parameter, DOC7 `@author` structure, and a `@return` describing
   post-fix behaviour including the individual renaming. Drop the `gl.compliance.check()` line
   from `@examples`. Run `devtools::document()` in the same change (F11, F14, DOC4).
10. Drop the `build = "Jody"` argument, delete the commented-out metadata-removal block, and
    normalise indentation on touched lines (F12, F13).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run. DEP: no Suggests packages are
  used, so DEP1 does not apply. PLT: the function produces no plots and writes no files, so
  PLT1-PLT3 and FS7 do not apply.
- Spec: behaviour vs roxygen on `platypus.gl`, `testset.gl`, `testset.gs` — run.
- Sampling verification matrix (all empirical, all on the reviewed state):

  | Mode | Verified | Result |
  |---|---|---|
  | per-population, `replace = TRUE` | counts, pool, population labels | correct except singleton populations (F1) |
  | per-population, `replace = FALSE` | counts, no duplicate draws | correct; errors when `nsample` exceeds the smallest population |
  | `onepop = TRUE` | total count, labels restored | correct; `nsample` semantics undocumented (F9) |
  | default `nsample` | resolves to `min(table(pop(x)))` | correct for `onepop = FALSE`; resolves to `nInd(x)` for `onepop = TRUE` (F9) |
  | seeded reproducibility | `indNames` and genotype matrix over repeated seeds | identical, no defect |
  | `nsample` as a proportion | not a documented mode | not supported; fails opaquely (F7) |
  | `nsample` = 0, 1, negative, `NA`, fractional | error paths | see F7 table |
  | `nsample` = smallest population size, `replace = FALSE` | boundary | correct |
  | unused population factor level | draws skip the empty level | correct, no error |
  | all-NA individual | draw succeeds | correct |
  | SNP / SilicoDArT dispatch | ploidy, allowed genotype values | correct on both |
  | FBM-backed object | runs, `@other` retained | runs, but diverges from the SNPbin path (F3) |
- **SNPbin duplicate-index hazard (DAT2, known package-wide defect): tested, NOT triggered.** The
  adegenet `SNPbin[]` defect corrupts `NA` when **loci** (columns) are subset with repeated
  indices, because the dartR `[` method delegates per-individual locus subsetting to
  `SNPbin`'s own `[` (`x@gen <- lapply(x@gen, function(e) e[jj])`,
  `R/utils.dartR.class.def.r:457`). `gl.sample` duplicates **individuals** (rows), which takes the
  branch at `R/utils.dartR.class.def.r:410` (`x@gen <- x@gen[i]`) and selects whole `SNPbin`
  objects from the list without ever calling `SNPbin`'s subsetter. Verified on `testset.gl`:
  `x[c(5,5,5,7,7,9), ]` reproduces the source rows exactly — 219 `NA` expected, 219 observed, zero
  NA-pattern mismatches, zero value mismatches. End to end,
  `gl.sample(testset.gl, nsample = 5, replace = TRUE)` returned 5261 `NA` against 5261 expected
  from the corresponding source rows, with zero mismatching cells. Pinned as a negative-control
  test in the baseline so a future change to the assembly path cannot introduce the corruption
  unnoticed.
- Input purity: `nInd`, `pop` and `names(@other)` of the input are unchanged after a call with
  `onepop = TRUE` — run, no defect.
- `verbose = 0` silence (VRB5): text side verified with `capture.output()` — zero lines. Plot side
  not applicable; the function has no plot bundle.
- Cross-package signature impact (API3): SKIPPED — not run for this review. `gl.sample` is
  exported and the proposed changes add validation and metadata without altering the signature,
  but sibling packages were not grepped.
- `R CMD check` / full package test suite: SKIPPED — out of scope for Phase A; the baseline test
  file was run in isolation (80 assertions, all passing).
- Related occurrence outside the review target, not investigated: `R/gl.report.hamming.r:225-226`
  uses `sample(idx, max.pairs, replace = TRUE)`, the same idiom as F1, and would exhibit the same
  length-1 behaviour if `idx` could have length 1. Flagged for whoever reviews that function.

### Note for the skill maintainer

Rule gaps this review had to work around, in addition to the `fixtures.md` items above:

- **No rule covers sampling or algorithmic correctness.** F1 is the most severe finding in this
  review and cites no rule ID, because the catalogue has no entry for "the function must draw what
  it says it draws". The `fixtures.md` table already uses an informal "Spec axis: ..." class for
  the same gap on `gl.filter.overshoot`. Proposal: a `SPEC` section, or at minimum a rule stating
  that stochastic functions must be verified empirically against every documented mode.
- **No rule requires stochastic output to be seed-reproducible.** `gl.sample` passes this check,
  but only by luck of implementation; nothing in the catalogue would have flagged it had it
  failed. Proposal: `TST4` — functions whose output depends on the RNG must produce identical
  results under a fixed seed, and reviews must verify it.
- **DOC5 leaned on again** — F9 and F11 both rest on it. That is the fifth and sixth use of a
  `[proposed]` rule across the campaign; DOC5 looks ready for ratification.
- **The proposed FS12 (history discipline) would have carried F5** more precisely than FS8, since
  the issue here is a missing append on a modify-family function combined with an assembly step
  that destroys the existing history.
- **The proposed DAT8 (derived-column refresh) is the natural home for F4** — stale locus metrics
  left flagged valid after the individuals change.

## Approval

All findings approved via the formal approval boxes on 2026-09-09: both BLOCKERs (F1, F2), all
three HIGHs (F3, F4, F5), all five MEDIUMs (F6, F7, F8, F9, F11) and all three LOWs (F10, F12,
F13). F14 is INFO and was applied as part of change 9. Consequences acknowledged at approval:
sampled objects change, the returned object retains `@other`, flags and history become honest, and
the FBM and SNPbin paths converge on one contract.

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges | 2026-09-09; sampled objects change for datasets with a singleton population |
| 2 | Approved | Arthur Georges | 2026-09-09; row order changes when nsample * nPop > nInd |
| 3 | Approved | Arthur Georges | 2026-09-09; downstream numerical output changes because metrics are recalculated |
| 4 | Approved | Arthur Georges | 2026-09-09 |
| 5 | Approved | Arthur Georges | 2026-09-09 |
| 6 | Approved | Arthur Georges | 2026-09-09; `nsample = 2.7` now errors |
| 7 | Approved | Arthur Georges | 2026-09-09 |
| 8 | Approved | Arthur Georges | 2026-09-09 |
| 9 | Approved | Arthur Georges | 2026-09-09; carries F14 |
| 10 | Approved | Arthur Georges | 2026-09-09 |

## Outcome

All ten changes applied on branch `review-gl.sample` from `upstream/dev` (`ddaed27`). Commit
`11e8deb`, the stranded fix on the local-only `dev` branch, supplied the F1 idiom and is credited
in the commit message; that branch was neither merged nor cherry-picked.

Implementation notes where the applied change differs from the proposal:

- **Change 8 (F9)** took the second of the two options F9 offered. `nsample` moves to a signature
  default of `NULL` and is resolved deliberately in the body: the size of the smallest population
  when `onepop = FALSE`, `nInd(x)` when `onepop = TRUE`. Forcing the old promise before `pop(x)`
  was overwritten would have changed the `onepop = TRUE` default from `nInd(x)` to the smallest
  population size, a result change that proposal 8 declared no consequence for. Sample sizes are
  therefore unchanged; only the mechanism is. Populations with no members (unused factor levels)
  are excluded from the default and from the draw.
- **Change 3 (F4)** sets the flags FALSE rather than recalculating, per the proposal's preferred
  option. All 13 entries of `loc.metrics.flags` are set FALSE, not a subset: every locus metric in
  that table is computed across individuals, so resampling invalidates all of them.
- **Change 5 (F6)** writes the new `indNames` through to `ind.metrics$id`. Duplicate individual
  names cannot arise: the zero-padded ordinal prefix is applied to every returned individual, so
  repeated draws of one individual get distinct names by construction. `make.unique()` is not
  needed and is not used. This is stated in `@details`.

### Discovered during Phase C, outside the review target

**The dartR `[` method cannot return more individuals than the object holds.**
`R/utils.dartR.class.def.r` assigns ploidy through the accessor (`ploidy(x) <- ploidy(x)[ii]`)
*after* `x@ind.names` has already been replaced, and the `ploidy` accessor names its result from
`@ind.names`. So `platypus.gl[c(1:81, 1:9), ]` fails with
`'names' attribute [90] must be the same length as the vector [81]`, on **both** the SNPbin branch
(line 410 region) and the FBM branch (line 380 region). Suggested fix, in the function that owns
the method rather than here: assign the slot directly, `x@ploidy <- x@ploidy[ii]`, before
`@ind.names` changes.

This falsifies F10's claim that the split/rbind chunking "serves no purpose": the chunking existed
to keep each `[` call within `nInd(x)`. It also means change 2 could not be applied as a single
unconditional `x[samps, ]`. Applied instead: a draw within `nInd(x)` — every draw of
`nsample <= the smallest population`, so every ordinary bootstrap — goes through `[` unchanged; a
larger draw is taken as the distinct individuals drawn and then expanded to the full draw by
repeating rows, with the individual-indexed elements of `@other` expanded by the same index. Both
routes were verified to produce identical structure and genotypes byte-for-byte against the source
rows, on both `options(dartR_fbm)` settings.

**The oversized draw returned desynchronised metadata on the FBM path at the reviewed state.**
`do.call(rbind, ...)` on FBM-backed chunks returned 90 individuals carrying an `ind.metrics` of 81
rows, silently. This is now correct, but it is a second reason not to leave the rbind assembly in
place, and it was not visible in Phase A because Phase A tested the oversized case only on the
SNPbin path.

### Verification

Run with `pdf(NULL)`, R 4.4.2, dartR.data 1.2.5, `devtools::load_all()` on the worktree.

- **F1** — `testset.gl`, `nsample = 3, replace = TRUE`: 90 individuals across all 30 populations,
  every population at exactly 3, `EmmacNormLeic` and `EmmacNormSalt` both present (they appeared
  in none of the 20 pre-fix runs). Across seeds 1-20 every population appears at its requested
  count, and each singleton population returns its own member three times, never a stranger. Every
  drawn individual's label matches its source population. The default call
  (`nsample = min(table(pop)) = 1`) returns 30 individuals, one per population, singletons
  included.
- **F2** — `platypus.gl`, `nsample = 4, replace = TRUE`: `@other` holds `loc.metrics`,
  `ind.metrics`, `latlon`, `loc.metrics.flags`, `verbose`, `history`. `ind.metrics` is 12 rows by
  17 columns (input width preserved), `latlon` 12 rows, `loc.metrics` 1000 rows. Every non-`id`
  column of `ind.metrics` and every `latlon` cell equals the input row for the drawn individual.
  Re-run with `nsample = 20, onepop = TRUE` (3 duplicated draws): metadata rows repeat with the
  duplicated individuals and `indNames` stay unique.
  `gl.filter.callrate(out, threshold = 0.5, plot.display = FALSE, verbose = 0)` completes and
  returns a `dartR` object of 12 individuals by 953 loci (from 1000); before the fix the same call
  raised `Error: incorrect number of dimensions`.
- **F3** — `platypus.gl`, `set.seed(1), nsample = 4, replace = TRUE`, run under
  `options(dartR_fbm = FALSE)` and `options(dartR_fbm = TRUE)`: `names(@other)`, `nInd`, `nLoc`,
  `dim(ind.metrics)`, `dim(loc.metrics)`, `dim(latlon)`, every flag, history length, `indNames`,
  `pop`, `ploidy` and `ind.metrics$id` all identical; genotypes equal (FBM returns doubles where
  SNPbin returns integers, values identical). Repeated for the oversized draw
  (`nsample = 30`, 90 draws): identical across paths on all the same fields.
- **F4** — flagged FALSE, not recalculated. On the FBM `testset.gl` case from the finding, the
  input carries `CallRate == TRUE` and the output now carries `FALSE`, with the stored `CallRate`
  still disagreeing with the recomputed value for 248 of 255 loci (max absolute difference
  0.9855) — that is the point of the flag. `gl.filter.callrate` on the output recalculates and
  keeps 227 of 255 loci, and its own output carries `CallRate == TRUE` again. All 13 flags are
  FALSE on both paths.
- **F5** — `platypus.gl` history goes from 4 entries to 5; the appended entry deparses to
  `gl.sample(x = x, nsample = ..., replace = TRUE, ...)`. On the FBM path history goes from 1 to
  2. Exactly one entry per call on both paths.
- **MEDIUMs** — `nsample = 2.7` and `nsample = 0.5` now stop with "nsample must be a whole number,
  it is a count of individuals and not a proportion"; `nsample = 0` and `nsample = -1` with
  "nsample must be 1 or more"; `nsample = NA` with "nsample must be a single positive whole
  number"; `nsample = 25, replace = FALSE` with "nsample (25) exceeds the size of the smallest
  population, SEVERN_BELOW with 17 individuals". The `onepop` default matches the documentation:
  `gl.sample(platypus.gl, onepop = TRUE)` returns 81 (`nInd(x)`), `gl.sample(platypus.gl)` returns
  51 (17 x 3). With `nsample = 30` (90 draws over 81 individuals) row order equals draw order for
  `pop`, `indNames` and the genotype matrix, on both paths.
- **SNPbin negative control** — unchanged and still passing. `testset.gl[c(5,5,5,7,7,9), ]`
  reproduces the source rows exactly: 219 `NA` expected, 219 observed, zero cell mismatches. End
  to end, `gl.sample(testset.gl, nsample = 5, replace = TRUE)` at `set.seed(1234)` returns 5249
  `NA` against 5249 expected from its source rows, with zero cell mismatches. The count differs
  from the Phase A figure of 5261 only because the F1 fix changes which rows that seed selects;
  the property under test — observed equals expected, zero mismatches — holds unchanged.
- **Genotype fidelity** — `as.matrix(out)` is identical to `as.matrix(x)[samps, ]` for every case
  tested: in-bounds and oversized draws, SNPbin and FBM, with and without duplicates.
- **Verbosity** — `verbose = 0` produces zero lines of output. `verbose = 1` produces 2 lines
  (Starting, Completed) and no progress line; `verbose = 2, 3, 5` produce 7 lines including
  exactly one "Completed" and the "Drew N individuals" progress line. The function has no plot
  bundle, so VRB5's plot side does not apply.
- **F14** — the documented example runs to completion with the `gl.compliance.check()` line
  removed: `gl.report.pa()` on the sampled object returns fixed-allele counts (53, 51, 56, 13, 5,
  11, 4, 4, 5 over the reduced grid).
- **Characterization suite** — `tests/testthat/test-gl.sample.R` rewritten from the 80-assertion
  baseline to 186 assertions, all passing. Every changed assertion carries a `# [approved Fn]`
  marker; no assertion changed without one.
- **Caller grep (API3, previously skipped)** — all eight clones under `D:\workspace\R\` plus
  `dartR.data`: no call site to `gl.sample()` anywhere outside its own `@examples`. The only
  references are `@seealso` cross-links in `man/*.Rd` of the "data manipulation" family, which the
  signature change does not touch. No breaking caller.
- **Still skipped** — `R CMD check` and the full package suite were not run; the change is
  confined to `R/gl.sample.r`, `man/gl.sample.Rd`, `NEWS.md` and
  `tests/testthat/test-gl.sample.R`, and no sibling function calls `gl.sample`.
  `R/gl.report.hamming.r:225-226` still carries the same `sample(idx, n, replace = TRUE)` idiom and
  remains for whoever reviews that function.

```json
{
  "function": "gl.sample",
  "package": "dartR.base",
  "family": "modify",
  "skill_version": "2.0.0",
  "model": "claude-opus-5",
  "model_phase_c": "claude-opus-5",
  "commit": "ddaed27",
  "verdict_standards": "rework",
  "verdict_spec": "rework",
  "fixture_status": "still_defective",
  "fixture_registered_in_fixtures_md": false,
  "fixture_note": "gl.sample has no row in fixtures.md; the sample(v,n) length-1 defect is present at ddaed27; a correct fix exists as commit 11e8deb on the local-only dev branch, not an ancestor of upstream/dev or of integration-local HEAD",
  "findings": [
    {"id": "F1", "severity": "BLOCKER", "confidence": "high", "rule": "SPEC (no rule ID)", "status": "applied", "change": 1},
    {"id": "F2", "severity": "BLOCKER", "confidence": "high", "rule": "DAT2", "status": "applied", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "DAT6", "status": "applied", "change": 2},
    {"id": "F4", "severity": "HIGH", "confidence": "high", "rule": "DAT4", "status": "applied", "change": 3},
    {"id": "F5", "severity": "HIGH", "confidence": "high", "rule": "FS8", "status": "applied", "change": 4},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "DAT2", "status": "applied", "change": 5},
    {"id": "F7", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "applied", "change": 6},
    {"id": "F8", "severity": "MEDIUM", "confidence": "high", "rule": "FS9", "status": "applied", "change": 7},
    {"id": "F9", "severity": "MEDIUM", "confidence": "medium", "rule": "DOC5", "status": "applied", "change": 8},
    {"id": "F10", "severity": "LOW", "confidence": "high", "rule": "STY3", "status": "applied", "change": 2},
    {"id": "F11", "severity": "MEDIUM", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 9},
    {"id": "F12", "severity": "LOW", "confidence": "high", "rule": "FS3", "status": "applied", "change": 10},
    {"id": "F13", "severity": "LOW", "confidence": "high", "rule": "STY1", "status": "applied", "change": 10},
    {"id": "F14", "severity": "INFO", "confidence": "high", "rule": "DOC3", "status": "applied", "change": 9}
  ],
  "approved": ["F1","F2","F3","F4","F5","F6","F7","F8","F9","F10","F11","F12","F13","F14"],
  "approved_by": "Arthur Georges",
  "approved_date": "2026-09-09",
  "rejected": [],
  "deferred": [],
  "coverage_skipped": [
    "R CMD check and full package suite: not run; change confined to gl.sample and its Rd, NEWS and test file",
    "gl.report.hamming.r:225-226 same idiom: outside review target, not investigated"
  ],
  "coverage_closed_in_phase_c": [
    "API3: all 8 clones plus dartR.data grepped; no call site outside gl.sample's own examples"
  ],
  "snpbin_duplicate_index_hazard": "not_triggered",
  "phase_c_discoveries": [
    "dartR '[' method rejects an index vector longer than nInd(x) on both the SNPbin and FBM branches: ploidy is assigned through the accessor after @ind.names has been replaced (R/utils.dartR.class.def.r). Falsifies F10's claim that the split/rbind chunking served no purpose. Not fixed here; belongs to the class-definition file.",
    "At the reviewed state, an oversized draw on the FBM path returned nInd individuals with a stale ind.metrics of the original row count, silently. Fixed as a side effect of change 2."
  ],
  "status": "pr-open",
  "pr": null
}
```
