# Review: gl.compliance.check (dartR.base)

- Family mode: modify (repairs the object: flags, metric recalculation, history)
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ddaed27 (upstream/dev). Empirics ran against a clean
  extraction of that tree; every pinned behaviour was re-verified identical on
  the local integration branch except `@position` handling, which follows
  PR #330 there (see Coverage).
- Datasets: testset.gl, testset.gs (dartR.data; clone at 1.2.4 for empirics),
  plus 20+ deliberately non-compliant synthetic genlights (out-of-range
  genotypes, NULL/mixed ploidy, desynced loc.metrics/ind.metrics, duplicate
  and NA names, all-NA rows/columns, all-monomorphic, bare adegenet objects)
- Baseline: tests/testthat/test-gl.compliance.check.R (new; 21 tests,
  51 assertions, all passing; known defects pinned as BUG(Fn))
- Context: this is the most heavily depended-on function in the package —
  every `gl.read.*` ends in it. PR #330 (position-genome-only) amends its
  `@position` hunk; that hunk is not re-reviewed here.

## Verdict

**Standards: Needs work** — the preamble, history append and return shape
conform, and `verbose = 0` is genuinely silent on the happy path; but two
warnings are ungated, one call appends two history entries, status lines
print at `verbose = 1`, and the roxygen header uses the outdated tag order
and verbose text.
**Spec: Rework** — the function's contract is "check compliance and rectify".
Against a 20-dimension matrix it detects and repairs the easy structural
gaps well, but its headline genotype-coding check is vacuous whenever any NA
is present (it prints "confirmed" over a genotype of 5), metadata desyncs
either crash opaquely or pass silently, the NULL-ploidy repair crashes on
exactly the non-uniform case it exists to handle, all-monomorphic input
crashes it, and one of its repairs (individual renaming) manufactures a new
desync. The checking core needs redesign, not spot fixes.

What works well: the structural repairs are solid — a bare adegenet genlight
comes out a compliant dartR object (pop1, ind.metrics, locus names, correct
recalculated CallRate); the function is idempotent (second run changes
nothing but history, verified with `all.equal`); and `verbose = 0` produces
zero lines on every non-crashing path tested.

## The detect/repair/pass/crash matrix

Empirical behaviour at ddaed27, one row per compliance dimension. "warn@n"
means the message appears only at `verbose >= n`.

| Dimension | Detect | Repair | Outcome |
|---|---|---|---|
| SNP coding outside {0,1,2,NA}, any NA present | no — `max(mat)` is NA, `NA %in% c(0,1,2,NA)` is TRUE | no | **prints "scored NA, 0, 1 or 2 confirmed"; a 5 survives** (F1) |
| SNP coding outside {0,1,2,NA}, no NA anywhere | max-only, warn@1 | no | non-fatal message, object returned; silent at v0 (F1) |
| SilicoDArT coding outside {0,1,NA} | same vacuous max check | no | false "confirmed" with NAs present (F1) |
| Negative genotypes | n/a | n/a | unconstructible — adegenet SNPbin refuses at construction |
| All loci monomorphic | — | — | **crash: "Subsetting resulted in zero loci."** (F2) |
| Some monomorphic loci | yes, warn@1 | flag `monomorphs=FALSE` (retained, by design) | ok |
| All-NA loci | yes, warn@1 (typo "dat") | flag `allna=FALSE` (retained) | ok |
| All-NA individuals | no — only `nLoc` compared | no | **silent pass; `allna` flag reads TRUE** (F7) |
| `@ploidy` NULL, uniform data | yes | ploidy stamped | ok |
| `@ploidy` NULL, heterogeneous computed ploidy | — | — | **crash: "the condition has length > 1"** (F3) |
| `@ploidy` set but non-uniform (e.g. 3,3,2,...) | no | no | silent pass, zero mentions at v3 (F8) |
| `loc.metrics` missing entirely | yes | created + recalculated | ok (junk first column name, F18) |
| `loc.metrics` rows != nLoc | designed warn unreachable | no | **crash: "replacement has 20 rows, data has 10"** in utils.reset.flags (F5) |
| `loc.metrics.flags` missing | yes | created, metrics recalculated | ok |
| `ind.metrics` missing | yes, warn@1 | id data frame created | ok |
| `ind.metrics` rows != nInd | no | no | silent pass (F6) |
| Duplicate indNames | yes, warn@2 | `make.unique` | **repairs names but not `ind.metrics$id` — creates desync** (F4) |
| Duplicate locNames | no | no | silent pass (F11) |
| NULL locNames | yes | `Loc1..n` | ok (silent at all levels) |
| NA locus name | no | no | **crash: "row names contain missing values"** in utils.recalc.maf/gl.alf (F10) |
| Missing pop | yes, warn@1 | pop1 assigned | ok |
| Missing `@loc.all` (SNP) | yes, warn@2 | dummy "A/C" | ok |
| latlong/long misspelling | yes | renamed latlon/lon | ok |
| Plain genlight (not dartR class) | yes | class converted | ok (message gated at v>2, F14) |
| `@position` | filled from SnpPosition at ddaed27 | — | superseded by PR #330 — not re-reviewed |
| History | — | own call appended | **plus a leaked internal `gl.recalc.metrics` entry** (F9) |

## Findings

**F1 [HIGH, confidence: high] — the genotype-coding check is vacuous with NAs present, non-fatal without them, and never repairs (DAT1, VRB3)**
`R/gl.compliance.check.r:108` (SNP) and `:128` (SilicoDArT) — the test is
`max(mat) %in% scores` with `scores <- c(0, 1, 2, NA)`. `max()` is called
without `na.rm`, so any missing value makes it NA, and `NA %in% scores` is
TRUE because `scores` contains NA. Verified: a matrix containing a genotype
of 5 plus one NA prints "SNP data scored NA, 0, 1 or 2 confirmed" and the 5
is returned intact. With zero NAs the check sees only the maximum (a 5 is
caught, but only because it is the max), prints an `error()`-styled message
gated at `verbose >= 1`, does not `stop()`, does not repair, and returns
the object; at `verbose = 0` the violation is invisible (0 lines captured).
Failure scenario: every real dataset has missing values, so the check has
never fired on realistic data — this is how the genotype-5 object from the
`gl.read.csv` review passed its reader's final compliance step with a green
"confirmed".
Proposed change: test the value set exactly
(`all(mat %in% scores | is.na(mat))` or a `setdiff` on observed values),
and on violation either `stop(error(...))` ungated or set the offending
cells to NA with a `verbose >= 1` warning — team decision on fatal vs
repair. **Consequence: objects that previously passed now error (or are
altered).**

**F2 [HIGH, confidence: high] — all-monomorphic input crashes the monomorph check (DAT5)**
`R/gl.compliance.check.r:192` — monomorphs are detected by filtering:
`x2 <- gl.filter.monomorphs(x, verbose = 0)` then comparing `nLoc`. When
every locus is monomorphic, `gl.drop.loc` subsets to zero loci and the
dartR `[` method stops with "Subsetting resulted in zero loci." (traceback:
gl.filter.monomorphs -> gl.drop.loc -> `[`). Reproduced with a 3-locus
all-monomorphic object; this is the crash surfaced when a degenerate VCF's
loci all went monomorphic after NA conversion.
Failure scenario: any degenerate-but-constructible dataset kills every
reader at its final step with a message that names neither the function nor
the cause.
Proposed change: count monomorphic loci directly (the loop in
`gl.filter.monomorphs` without the drop, or a vectorised scan), or guard
the `nLoc(x2) == 0` case; report "all loci monomorphic" and return with the
flag set FALSE. The same filter-to-count idiom at `:209`
(`gl.filter.allna`) crashes the same way if every locus is all-NA.

**F3 [HIGH, confidence: high] — the NULL-ploidy repair crashes on the heterogeneous case it exists for (DAT1, FS5)**
`R/gl.compliance.check.r:59` — inside `if (is.null(x@ploidy))`, the branch
tests `if (unique(ploidy(x)) == 2)`. With a NULL slot, adegenet's
`ploidy()` computes per-individual maxima, so mixed data (individuals whose
max genotype is 1 next to individuals with a 2) yields a vector — verified
computed ploidy `1,1,2,2` — and `if` fails with "the condition has
length > 1". Only uniformly-diploid or uniformly-haploid data survive the
branch. The polyploid `else` at `:66` is near-unreachable, its
`cat(warn(...))` is ungated (prints at `verbose = 0`), and
`ploidy(x) <- ploidy(x)` is a no-op.
Failure scenario: a reader hands over an object with NULL ploidy and any
mix of per-individual maxima — the user gets an R-internals error that
names no dartR function; this is the "condition has length > 1" failure
that forced the `gl.read.vcf` fix to stamp uniform ploidy itself.
Proposed change: derive one ploidy for the whole object
(`max(ploidy(x))` or from the global max genotype), assign it uniformly,
warn at `verbose >= 1`; gate the polyploid warning; drop the no-op.

**F4 [HIGH, confidence: high] — the duplicate-name repair desynchronises ind.metrics (DAT2)**
`R/gl.compliance.check.r:241-251` — duplicated `indNames` are repaired with
`make.unique(indNames(x), sep = '_')`, but `x@other$ind.metrics$id` (and
any other id-keyed column) is left holding the old duplicated names.
Verified: input ids `dup, dup, ...` come back as `indNames = dup, dup_1`
with `ind.metrics$id = dup, dup` — the repair manufactures exactly the
genotype-metadata desync the function exists to catch, and the run reports
"Individual metrics confirmed".
Failure scenario: any join or lookup on `ind.metrics$id` downstream
matches the wrong individual or drops rows.
Proposed change: apply the same `make.unique` result to
`x@other$ind.metrics$id` when that column exists (the ind.metrics block at
`:258` runs later, so the order also needs care for objects that get the id
column created from already-repaired names — that path is currently
correct).

**F5 [HIGH, confidence: high] — loc.metrics row desync crashes opaquely; the designed diagnostic is unreachable and repairs nothing (DAT2)**
`R/gl.compliance.check.r:224-233` — the row-count comparison and its
"major problem... Trace back to identify the cause" warning sit AFTER
`utils.reset.flags` (`:180`) and `gl.recalc.metrics` (`:186`), both of
which assign length-`nLoc(x)` vectors into the desynced data frame first.
Verified: 20 loci with 10 metric rows dies inside `utils.reset.flags` with
"replacement has 20 rows, data has 10" (traceback: `$<-.data.frame` on the
`monomorphs` column) — the designed warning never printed in any tested
desync. When it could fire it is also ungated (leaks at `verbose = 0`) and
is warn-only: the object is returned desynced.
Failure scenario: the one metadata fault the details section singles out as
"potentially a major problem" produces a bare base-R error with no dartR
context, at every verbosity.
Proposed change: move the row-count check before `utils.reset.flags` and
make it `stop(error(...))` with the counts in the message (no safe
automatic repair exists — the row-to-locus mapping is unknown). Check
`ind.metrics` the same way (F6).

**F6 [MEDIUM, confidence: high] — ind.metrics row count is never checked (DAT2)**
`R/gl.compliance.check.r:253-269` — the individual-metrics block only tests
for NULL. Verified: `ind.metrics` with 5 rows against 10 individuals passes
in silence and returns still desynced, while the loc.metrics side at least
attempts a warning.
Failure scenario: individual metadata off by any count survives the
package's compliance gate and misaligns every ind.metrics consumer.
Proposed change: `nrow(x@other$ind.metrics) != nInd(x)` -> fatal error
(same rationale as F5).

**F7 [MEDIUM, confidence: high] — all-NA individuals pass undetected and the allna flag claims clean (DAT2, DOC5 proposed rule)**
`R/gl.compliance.check.r:209-220` — `gl.filter.allna` drops all-NA
individuals as well as loci, but the check compares only `nLoc(x2)`.
Verified: an object with one fully-NA individual reports "No loci with all
missing data detected", keeps the individual, and returns
`allna flag = TRUE`.
Failure scenario: an individual with zero genotypes flows into distance,
kinship and assignment analyses that assume at least some data per
individual.
Proposed change: also compare `nInd(x2)` and report/flag all-NA
individuals separately.

**F8 [MEDIUM, confidence: high] — non-uniform stamped ploidy passes without a word (DAT1)**
`R/gl.compliance.check.r:56-72` — when `@ploidy` is non-NULL nothing
inspects it. Verified: ploidy `3,3,2,2,...` runs to completion with zero
ploidy-related output at `verbose = 3` and the vector preserved;
`utils.check.datatype` resolves it to "SNP" through its
`is.null(ploidy(x))==F` fallback. dartR's data model (DAT1) admits only
uniform 2 (SNP) or uniform 1 (SilicoDArT).
Failure scenario: the mixed-ploidy objects `gl.read.vcf` used to produce
are certified compliant; dosage-based maths downstream is silently wrong.
Proposed change: validate `length(unique(ploidy(x))) == 1` and the value is
in {1, 2}; warn (or error — team decision) otherwise.

**F9 [MEDIUM, confidence: high] — one call appends two history entries (FS8)**
`R/gl.compliance.check.r:186,317-319` — `gl.recalc.metrics` appends its own
internal `match.call()` (`gl.recalc.metrics(x = x, verbose = 0)`) before
the function appends its own. Verified: history grows by 2 per call;
shipped `testset.gl` already carries the stray entry. `gl.filter.monomorphs`
(`:131-135` of its file) documents the house idiom: discard internal
delegation entries, one entry per user call. Because every reader ends in
this function, every freshly read object starts life with a polluted
history.
Proposed change: snapshot `x@other$history` before the delegation block and
restore it before appending the single `gl.compliance.check` entry.

**F10 [MEDIUM, confidence: medium] — an NA locus name crashes metric recalculation opaquely (DAT2)**
`R/gl.compliance.check.r:95-97` — only `is.null(locNames(x))` is repaired.
An NA among the names survives to `gl.recalc.metrics` ->
`utils.recalc.maf` -> `gl.alf`, which dies with "row names contain missing
values" (verified traceback).
Failure scenario: one corrupt locus label in an imported dataset turns the
compliance gate into an unexplained base-R error.
Proposed change: repair NA (and empty-string) names in the same block:
`ln[is.na(ln)] <- paste0("Loc", which(is.na(ln)))`, warn at
`verbose >= 2`.

**F11 [MEDIUM, confidence: high] — duplicate locus names pass silently while duplicate individual names are repaired (DAT3)**
`R/gl.compliance.check.r:241-251` repairs indNames; no counterpart exists
for locNames. Verified: duplicated locus names survive the run. DAT3
records why this matters: name-based locus dropping is unsafe under
duplicates, and this very function delegates to `gl.filter.monomorphs` ->
`gl.drop.loc`, which drops by name — a duplicated name can drop a
polymorphic locus that shares its label with a monomorphic one (in the
delegate's discarded copy today, in real subsets downstream).
Proposed change: `make.unique` on locNames with a `verbose >= 2` warning,
mirroring the individual-name repair.

**F12 [MEDIUM, confidence: high] — three readers call it without verbose, defeating their own verbose=0 (VRB5, API3)**
The function itself is clean: `verbose = 0` produced zero lines on every
non-crashing fixture (verified with `capture.output`), and the default
resolves to 2 -> 20 lines. The leak reported across the reader reviews is
caller-side: `gl.read.dart.r:282`, `gl.read.silicodart.r:375` and
`gl.read.vcf.r:127` call `gl.compliance.check(glout)` with no `verbose`
argument (confirmed at ddaed27), so a silent read still prints the full
compliance log. All other callers forward `verbose` or pass 0 (table
below).
Failure scenario: `gl.read.dart(..., verbose = 0)` emits ~20 lines.
Proposed change: pass `verbose = verbose` at the three call sites
(caller-side change, one line each).

**F13 [MEDIUM, confidence: high] (proposed rule) — the documentation promises detect-and-rectify that the code does not deliver, and omits real side effects (DOC5)**
`R/gl.compliance.check.r:6-31` — "@description ... will check to see that
the genlight object conforms ... and if it does not, will rectify it."
Against the matrix above: coding violations are neither reliably detected
nor ever rectified (F1); metadata desyncs crash or pass (F5, F6); all-NA
individuals pass (F7). Meanwhile the function does things the header never
mentions: it recalculates all locus metrics unconditionally at `:186`
(overwriting supplied values, e.g. DArT-computed CallRate), converts the
object to class `dartR`, resets metric flags and stamps `@other$verbose`.
Failure scenario: a custodian or user reads the header and trusts the gate;
the trust is misplaced in both directions.
Proposed change: rewrite `@description`/`@details` to enumerate the actual
checks, repairs and side effects — after the code-side findings above are
decided, so docs describe the agreed behaviour.

**F14 [LOW, confidence: high] — message-layer gating and text defects (VRB1, VRB3)**
Four small defects, one change: (a) `:66` polyploid `cat(warn(...))` and
`:225` desync warning are ungated — they print at `verbose = 0`; (b) the
genlight->dartR conversion notice at `:78` is gated `verbose > 2`, so it
appears at 3 but not at 2 where progress belongs; (c) result-status lines
("SNP data scored... confirmed", "Population assignments confirmed", etc.)
are gated `verbose >= 1`, producing 7 lines at level 1 where VRB1 promises
begin/end only; (d) `:217` typo "all missing dat", and several messages
embed source-indentation linebreaks that render ragged (visible in the v2
transcript).
Proposed change: gate (a) at `verbose >= 2`/`>= 1`, move (b) to
`verbose >= 2`, move (c) to `verbose >= 2`, fix the string literals.

**F15 [LOW, confidence: high] — roxygen header: outdated tag order, nonstandard verbose text, custodian-only author (DOC1, DOC2, DOC7 proposed)**
`R/gl.compliance.check.r:1-41` — `@param` precedes `@details` and `@return`
sits last after `@export` (the pre-2026-08-27 order; DOC1 says flag the
file); the `verbose` text reads "[default 2 or as specified using
gl.set.verbosity]" instead of the DOC2 canon stating the NULL ->
`options(dartR_verbose)` -> 2 cascade; `@author` names only a custodian
with no `Author(s):` part (DOC7, proposed rule).
Proposed change: reorder tags to the house order, adopt the DOC2 verbose
text, add `Author(s):` (defaulting to the custodian per DOC7); run
`devtools::document()` in the same change (DOC4).

**F16 [LOW, confidence: medium] (proposed rule) — three full densifications per run (DAT6, STY2)**
`:103`/`:123` densify the whole object merely to bound-check values, and
the `gl.filter.monomorphs` and `gl.filter.allna` delegates each call
`as.matrix(x)` again on the full object. On an FBM-backed dartR object this
materialises the matrix three times per compliance run — and every read
pipeline ends here.
Proposed change: when the checks are redesigned (F1, F2), scan column-wise
or via accessors; at minimum reuse one `as.matrix` result for the in-file
checks.

**F17 [INFO, confidence: high] — dead code and misordered banners (STY1)**
`:156` `is.na(length(pop(x)))` is always FALSE (`length()` never returns
NA); `:67` `ploidy(x) <- ploidy(x)` is a no-op (F3 removes it); the
`# DO THE JOB` banner at `:54` precedes `# CHECKS DATATYPE` at `:71` — the
ploidy pre-check genuinely must run before `utils.check.datatype` (which
stops on NULL ploidy by telling the user to run this very function), but it
deserves its own banner so the FS4->FS6 order reads intentionally.
Proposed change: fold into whichever of F1-F5 touches each region.

**F18 [LOW, confidence: medium] — repairing a metrics-less object leaves a junk column named `array(NA, nLoc(x))` (DAT2)**
`utils.reset.flags.r` (upstream `:91`) creates a missing loc.metrics as
`as.data.frame(array(NA, nLoc(x)))`, whose auto-generated column name is
the literal expression text. Verified: a bare genlight comes back with
loc.metrics columns `array(NA, nLoc(x)), AvgPIC, ...` — the junk column
persists in the returned object. The defect is in the delegate but surfaces
through this function's contract.
Proposed change: in `utils.reset.flags`, create the frame empty
(`data.frame(row.names = seq_len(nLoc(x)))`) or name the seed column.

## Caller inventory

`gl.compliance.check` is called from 15 sites, all inside dartR.base; the
other seven clones (captive, popgen, sim, spatial, sexlinked, dartRstartup,
dartRverse) contain no calls. Every caller uses the return value. Inventory
taken on the working tree and cross-checked at ddaed27 for the readers.

| Caller | Site | verbose passed | Note |
|---|---|---|---|
| gl.read.dart | R/gl.read.dart.r:282 | **none (defaults to 2)** | F12 |
| gl.read.silicodart | R/gl.read.silicodart.r:375 | **none** | F12 |
| gl.read.vcf | R/gl.read.vcf.r:127 | **none** | F12 |
| gl.read.csv | R/gl.read.csv.r:252, 355 | `verbose = verbose` | |
| gl.read.fasta | R/gl.read.fasta.r:94 | `verbose = verbose` | |
| gl.read.PLINK | R/gl.read.PLINK.r:202, 304 | `verbose = verbose` | |
| gl.load | R/gl.load.r:78 | `verbose = verbose` | |
| gi2gl | R/gi2gl.R:64 | `verbose = verbose` | |
| gl.dist.pop | R/gl.dist.pop.r:125 | `verbose = 0` | |
| gl.drop.pop / gl.keep.pop | R/gl.drop.pop.r:93, R/gl.keep.pop.r:94 | `verbose = 0` | |
| gl.impute | R/gl.impute.r:473 | `verbose = 0` | |
| gl.report.maf | R/gl.report.maf.r:134 | `verbose = 0` | |
| gl.sim.crosses / gl.sim.genotypes | R/gl.sim.crosses.r:156, R/gl.sim.genotypes.r:72 | `verbose = 0` | |

Implication: the function's verbosity behaviour is correct; the three
readers are the leak. Its crash modes (F2, F3, F5, F10), however, propagate
to all 15 callers, which is why the crash findings carry HIGH severity.

## Proposed changes

1. Replace the genotype-coding check with an exact value-set test, ungated;
   on violation stop with `error()` (or repair offending cells to NA with a
   level-1 warning — team decision) (F1, part of F16).
   **Consequence: previously "confirmed" corrupt objects now error.**
2. Detect monomorphic and all-NA loci by counting, not by filtering-and-
   comparing; handle the all-monomorphic and all-all-NA cases; also compare
   `nInd` so all-NA individuals are reported and flagged (F2, F7, part of
   F16).
3. Rewrite the NULL-ploidy branch to derive a single object-wide ploidy,
   gate its warning, drop the no-op; add a validation that stamped ploidy is
   uniform and in {1, 2} (F3, F8).
   **Consequence: mixed-ploidy objects now warn or error instead of passing.**
4. Mirror the `make.unique` individual renaming into `ind.metrics$id`
   (F4).
5. Check `loc.metrics` and `ind.metrics` row counts up front, before
   `utils.reset.flags`/`gl.recalc.metrics`, failing fatally with counts in
   the message (F5, F6).
   **Consequence: the opaque "replacement has N rows" crash becomes an
   informative fatal error; silently-desynced ind.metrics now error.**
6. Repair NA/empty locus names and de-duplicate locus names alongside the
   existing NULL repair (F10, F11).
7. Restore the pre-delegation history and append exactly one entry per call
   (F9).
8. Pass `verbose = verbose` at the three reader call sites
   (`gl.read.dart:282`, `gl.read.silicodart:375`, `gl.read.vcf:127`) (F12).
9. Fix message gating, the v>2 conversion notice, the "dat" typo and the
   ragged embedded-newline strings (F14).
10. Rewrite the roxygen header: actual checks/repairs/side effects, house
    tag order, DOC2 verbose text, DOC7 author line; re-document (F13, F15).
11. Housekeeping: remove dead conditions, re-banner the ploidy pre-check
    (F17); fix the junk seed column in `utils.reset.flags` (F18).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, STY — run. DEP: no Suggests
  packages are used; nothing to guard. PLT: no plots; not applicable.
- Spec: behaviour vs roxygen exercised empirically on a 20+ dimension
  fixture matrix (testset.gl, testset.gs, synthetic non-compliant objects);
  crash sites confirmed by traceback.
- Idempotence: run-on-own-output verified identical bar history on
  testset.gl.
- Verbosity: v0/v1/v2/v3 transcripts captured; v0 empirically silent on all
  non-crashing paths (VRB5 text side; no plot side exists).
- Caller inventory: all 8 clones grepped; readers cross-checked at ddaed27.
- PR #330 interaction: the `@position` hunk (`:293-315` on the branch) was
  read first and is not re-reviewed; no finding above touches those lines —
  no conflicts-with-#330. The `@chromosome` slot is untouched by both the
  function and this review (genome-only convention leaves it NULL).
- FBM path (DAT6): SKIPPED empirically — no FBM fixture wired for this
  function; F16 is from reading the code.
- Zero-locus input: SKIPPED — not constructible through dartR subsetting
  (the `[` method refuses), and direct construction of a 0-locus genlight
  was not pursued.
- platypus.gl: not used — testset.* plus synthetic fixtures covered every
  dimension; documented here per the fixtures rule.
- Empirics ran on a `git archive` extraction of ddaed27 (the working tree
  carries campaign fixes to helpers); every pinned behaviour was then
  re-verified on the working tree — identical except `@position` (post-#330
  there, and the baseline test pins the post-#330 behaviour since that is
  what the repo's tests run against).
- dartR.data: empirics loaded the local clone (1.2.4); the baseline test
  run used the installed dartR.data via `load_all` of dartR.base.

## Approval (Phase B)

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved, minimal variant | Arthur Georges, 2026-09-05 | F1 only: fix the NA-vacuity (exact value-set test). Keep gating, warn-not-stop, no repair, no strengthened detection. The F16 densification rework bundled here is deferred. |
| 2 | approved in part, minimal variant | Arthur Georges, 2026-09-05 | F2 crash guards only: the all-monomorphic crash and the same idiom's all-all-NA crash must complete flag-only. No counting redesign of the filter-to-count idiom; the F7 nInd extension is deferred. |
| 3 | deferred | Arthur Georges, 2026-09-05 | F3, F8 |
| 4 | deferred | Arthur Georges, 2026-09-05 | F4 |
| 5 | deferred | Arthur Georges, 2026-09-05 | F5, F6 |
| 6 | deferred | Arthur Georges, 2026-09-05 | F10, F11 |
| 7 | deferred | Arthur Georges, 2026-09-05 | F9 |
| 8 | deferred | Arthur Georges, 2026-09-05 | F12 |
| 9 | deferred | Arthur Georges, 2026-09-05 | F14 |
| 10 | deferred | Arthur Georges, 2026-09-05 | F13, F15 |
| 11 | deferred; F17 (INFO) no action | Arthur Georges, 2026-09-05 | F18 deferred |

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-gl.compliance.check` (base ddaed27),
at the deliberately minimal scope decided by the member. Only
`R/gl.compliance.check.r` was touched in the R tree; the roxygen header is
unchanged, so no re-document was needed.

- **F1 (minimal)**: both coding checks (SNP `:108`, SilicoDArT `:128`) now
  test the value set exactly — `all(mat %in% scores)` in place of
  `max(mat) %in% scores` (`NA %in% scores` matches the NA entry, so missing
  values pass; any out-of-range value fails). Gating (`verbose >= 1`),
  warn-not-stop and no-repair semantics retained. The genotype-5-with-NA
  fixture now prints "Error: SNP data must be scored NA, 0 or 1 or 2"
  instead of "confirmed"; the no-NA path is behaviourally unchanged; clean
  testset.gl/testset.gs still print their "confirmed" lines.
- **F2 (minimal)**: the two filter-to-count sites (`:192`
  `gl.filter.monomorphs`, `:209` `gl.filter.allna`) are wrapped in
  `tryCatch`; the "Subsetting resulted in zero loci" error is read as the
  all-monomorphic / all-all-NA case — the relevant flag is set FALSE, a
  `verbose >= 1` warn prints ("All loci are monomorphic" / "All loci have
  all missing data"), and the function completes. No redesign of the idiom.

Verification (baseline characterization + empirics, R 4.4.2):

- Test file after pin updates: 22 tests / 58 assertions. Only the BUG(F1)
  and BUG(F2) pins flipped (annotated `[approved F1 minimal]` /
  `[approved F2 minimal]`); every other pin — repairs, F1's retained
  warn-not-stop path, and BUG(F3)-BUG(F11) — passes unchanged before and
  after the patch.
- The two `@position` assertions fail identically pre- and post-patch:
  the branch base ddaed27 predates PR #330, whose behaviour the baseline
  pins (documented in Coverage). Not a behaviour change of this PR.
- Empirics: genotype-5-with-NA now reports the violation at `verbose = 2`
  and the 5 survives (no repair); all-monomorphic completes silent at
  `verbose = 0` with `monomorphs = FALSE`; all-all-NA completes at
  `verbose = 0/1` with `allna = FALSE` and `monomorphs = FALSE`; clean
  testset.gl `verbose = 3` transcript unchanged, `verbose = 0` emits zero
  lines, idempotent bar history; testset.gs confirm line unchanged.

Residual, recorded not applied: with every locus all-NA and
`verbose >= 2`, `utils.check.datatype()`'s courtesy check
(`utils.check.datatype.R:87-99`) calls `gl.filter.allna(x, verbose = 0)`
unguarded and errors "subscript out of bounds" before this function's body
proceeds — verified identical at unpatched ddaed27, so pre-existing and
untouched by this PR. Same root cause as F2 (`gl.filter.allna` cannot
return a zero-locus object); belongs to a review of `utils.check.datatype`.
The all-all-NA completion guarantee therefore holds at `verbose` 0 and 1.

Deferred by the member's 2026-09-05 decision: F3, F4, F5 (HIGH), F6-F13
(MEDIUM), F14, F15, F16, F18 (LOW) — 15 findings recorded here, not
applied. F17 (INFO): no action.

```json
{
  "function": "gl.compliance.check",
  "package": "dartR.base",
  "family": "modify",
  "skill_version": "1.0.0",
  "commit": "ddaed27",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT1", "status": "applied-minimal", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DAT5", "status": "applied-minimal", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "DAT1", "status": "deferred", "change": 3},
    {"id": "F4", "severity": "HIGH", "confidence": "high", "rule": "DAT2", "status": "deferred", "change": 4},
    {"id": "F5", "severity": "HIGH", "confidence": "high", "rule": "DAT2", "status": "deferred", "change": 5},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "DAT2", "status": "deferred", "change": 5},
    {"id": "F7", "severity": "MEDIUM", "confidence": "high", "rule": "DAT2", "status": "deferred", "change": 2},
    {"id": "F8", "severity": "MEDIUM", "confidence": "high", "rule": "DAT1", "status": "deferred", "change": 3},
    {"id": "F9", "severity": "MEDIUM", "confidence": "high", "rule": "FS8", "status": "deferred", "change": 7},
    {"id": "F10", "severity": "MEDIUM", "confidence": "medium", "rule": "DAT2", "status": "deferred", "change": 6},
    {"id": "F11", "severity": "MEDIUM", "confidence": "high", "rule": "DAT3", "status": "deferred", "change": 6},
    {"id": "F12", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "deferred", "change": 8},
    {"id": "F13", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "proposed_rule": true, "status": "deferred", "change": 10},
    {"id": "F14", "severity": "LOW", "confidence": "high", "rule": "VRB3", "status": "deferred", "change": 9},
    {"id": "F15", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "deferred", "change": 10},
    {"id": "F16", "severity": "LOW", "confidence": "medium", "rule": "DAT6", "proposed_rule": true, "status": "deferred", "change": 1},
    {"id": "F17", "severity": "INFO", "confidence": "high", "rule": "STY1", "status": "no-action", "change": 11},
    {"id": "F18", "severity": "LOW", "confidence": "medium", "rule": "DAT2", "status": "deferred", "change": 11}
  ],
  "coverage_skipped": [
    "DAT6: no FBM fixture (static reading only)",
    "zero-locus input: not constructible via dartR subsetting",
    "PR #330 position hunk: excluded by campaign instruction",
    "platypus.gl: not needed, testset.* + synthetic fixtures covered the matrix"
  ],
  "status": "phase-c-complete",
  "pr": 367
}
```
