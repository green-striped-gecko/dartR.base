# Review: gl.read.csv (dartR.base)

- Family mode: io
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203 (R/gl.read.csv.r byte-identical to upstream/dev
  ddaed27; the local tree's gl.compliance.check carries the PR #330 hunks,
  so @position observations note the ddaed27 difference where it matters)
- Datasets: platy_test.csv + platy_ind.csv (dartR.data extdata); synthetic
  fixtures (inventory below)
- Baseline: tests/testthat/test-gl.read.csv.R (29 assertions, snapshot
  captured pre-review, defective behaviour pinned with [BUG] tags)

## Verdict

**Standards: Needs work** — preamble, history, double compliance pass and
`verbose = 0` silence all conform; the error-exit idiom and several guard
gaps need bounded fixes.
**Spec: Needs work** — both documented genotype encodings convert correctly
(hand-verified cell by cell, hets in both orders, ref = most frequent
allele, `transpose` round-trips), but the documented `loc.metafile` feature
crashes on every use, an id column anywhere but first is spuriously
rejected, and out-of-range numeric genotypes are admitted silently.

What works well: the core conversion is sound — numeric 0/1/2/NA passes
through cell-perfect, character data recodes correctly against the
most-frequent-allele reference exactly as documented, and ploidy 2 (DAT1),
metadata row-registration (DAT2) and the compliance pass at return (DAT5)
all check out empirically.

## Synthetic fixtures

All written by the baseline test file itself (tempdir) and mirrored in the
review scratchpad; 8 individuals x 6 loci unless stated.

| Fixture | Purpose |
|---|---|
| `num_basic.csv` | numeric 0/1/2 with NA cells; cell-level fidelity |
| `num_bad.csv` | num_basic with genotype 5 at i3/L2; validation probe |
| `char_basic.csv` | 5 loci A/C, G/T (ref = T by frequency), monomorphic A, C/T with `-/-`, C/G; hets in both orders |
| `char_tri.csv` | char_basic with third allele G/G at CL1/i8; biallelic guard |
| `char_halfmiss5.csv` (scratchpad probe) | `A/-` half-missing cell |
| `char_transposed.csv` | transpose of char_basic; `transpose = TRUE` |
| `small_num.csv` | 3 x 3; the 5x5 sniff crash |
| `ind_meta_good.csv` | id,pop,sex in order |
| `ind_meta_idsecond.csv` | same rows, id as second column |
| `ind_meta_reorder.csv` (scratchpad probe) | rows shuffled; order guard |
| `ind_meta_nopop.csv` | id,sex only; pop fallback |
| `loc_meta_good.csv` | AlleleID,SnpPosition,RepAvg in locus order |
| `loc_meta_noAID.csv` (scratchpad probe) | AlleleID column absent |
| `loc_meta_short.csv` (scratchpad probe) | 4 rows for 6 loci |

## Findings

**F1 [HIGH, confidence: high] — any loc.metafile crashes the read (FS5, DAT2)**
`R/gl.read.csv.r:273-282` — the AlleleID cross-check loop compares against
`gl@other$loc.metrics$AlleleID`, but the loc.metrics built by
`gl.compliance.check`/`gl.recalc.metrics` for a csv read has no `AlleleID`
column, so the guard reads `NULL[i]` and `if()` fails with "argument is of
length zero". Reproduced with a fully valid metafile (AlleleID first, rows
in order): every invocation with `loc.metafile` set dies before the metrics
are attached. Two secondary defects sit in the same block: the
missing-AlleleID guard at line 265 prints "Fatal Error" but has no `stop()`,
so execution continues into the same crash; and the row comparison uses
`loc.metrics[i, 1]` — the first column — rather than the `AlleleID` column
the documentation mandates.
Failure scenario: any user following the documented workflow for locus
metadata gets an uninterpretable length-zero error; the feature is unusable.
Proposed change: validate the incoming file's own `AlleleID` column against
`locNames(gl)` (position-wise, per the documented same-order contract),
`stop()` when the column is absent, and only then assign to
`gl@other$loc.metrics`.

**F2 [HIGH, confidence: high] — out-of-range numeric genotypes admitted silently (DAT1)**
`R/gl.read.csv.r:230` — the numeric-branch guard is
`any(names(tmp) %in% c("0","1","2","NA")) == FALSE`: it errors only when NO
value is a legal code. A file containing 0, 1, 2 and a stray 5 passes, and
the 5 survives into the returned object (verified: cell i3/L2 == 5 in the
final genlight; `gl.compliance.check` even prints "SNP data scored NA, 0, 1
or 2 confirmed" over it — a compliance gap worth flagging to that
function's custodian, DAT5).
Failure scenario: a typo or a 0/1 presence-absence file mixed with dosage
codes yields a nominally diploid SNP object whose downstream allele
frequencies are silently wrong.
Proposed change: `all(...)` instead of `any(...)`, and report the offending
values in the error message.

**F3 [MEDIUM, confidence: high] — files with fewer than 5 loci or 5 individuals crash (FS5)**
`R/gl.read.csv.r:113-121,151-154` — the type sniff indexes `data[1:5, 1:5]`
and the confirmation print indexes `df0[1, 2:6]` / `df0[2:6, 1]`
unconditionally. A 3-individual, 3-locus file dies with "subscript out of
bounds" at `verbose = 0` and "undefined columns selected" at
`verbose = 2` (both reproduced).
Failure scenario: small teaching or toy datasets — precisely the audience
for a csv reader — crash with errors that do not name the real problem.
Proposed change: clamp the sniff window to
`1:min(5, nrow(data))` / `1:min(5, ncol(data))` and the confirmation print
likewise.

**F4 [MEDIUM, confidence: high] — id column must be first, contrary to the docs (DOC5, proposed rule)**
`R/gl.read.csv.r:316` — the order check compares `ind.metrics[i, 1]`
(first column) instead of `ind.metrics$id`. A metafile with columns
`pop,id,sex` — valid per the docs, which require only "a mandatory id
column ... in the same order" — is rejected with "Fatal Error: id in the
individual metrics file does not correspond..." (reproduced). The same
first-column pattern is what F1's loop uses on the locus side.
Failure scenario: a correct metafile is refused with a message that sends
the user hunting for a nonexistent ordering problem.
Proposed change: compare `ind.metrics$id` (and the locus file's
`AlleleID`) by name, not by position.

**F5 [MEDIUM, confidence: high] — real allele pairs computed then discarded; loc.all is a fabricated placeholder (DAT2, DOC5 proposed rule)**
`R/gl.read.csv.r:193-215,242-249` — the character branch derives the true
ref/alt pair per locus (`homRef`, `homAlt`) but never stores it; the
genlight is built without allele labels and the compliance placeholder is
returned: `loc.all` reads `A/C` for every locus (verified on a fixture
whose loci are A/C, G/T, C/G — all reported as `A/C`).
Failure scenario: any downstream export that prints alleles (gl2vcf,
gl2fasta, gl2plink) emits wrong REF/ALT spellings under correct locus
names — plausible wrong data.
Proposed change: collect `paste0(names(tmp)[1], "/", names(tmp)[2])` per
locus in the conversion loop and assign to `alleles(gl)` before the
compliance pass.

**F6 [LOW, confidence: high] — error exits use cat(error()) + bare stop() (VRB2, FS5)**
`R/gl.read.csv.r:134-147,169-191,229-237,307-325` — fatal paths print the
red message with `cat(error(...))` then call `stop()` with no message; the
thrown condition is empty (probes show `ERROR:` followed by nothing).
`tryCatch`/testthat callers and logs capture nothing.
Failure scenario: scripted pipelines cannot distinguish which of the six
fatal exits fired.
Proposed change: the house idiom `stop(error("Fatal Error: ...\n"))`
throughout.

**F7 [LOW, confidence: high] — misleading and malformed progress messages (DOC5 proposed rule, VRB2)**
`R/gl.read.csv.r:287-295` — the `verbose >= 2` report after
`gl.recalc.metrics` prints one line per locus-metrics column but says
"...to the other$ind.metrics slot" (it is loc.metrics), and the first line
leaks a literal column named `array(NA, nLoc(x))` — that column name is
manufactured by `gl.compliance.check`'s loc.metrics constructor (flag to
that custodian) but this reader displays it. The two metafile advisories at
lines 72-87 are labelled "Warning:" yet issued through `report()` (green)
rather than `warn()`.
Proposed change: correct the slot name, collapse the column list into one
line, and use `warn()` for warnings.

**F8 [LOW, confidence: medium] — inconsistent file-encoding handling across the three inputs (FS5)**
`R/gl.read.csv.r:95-98,258-264,299-306` — only the ind.metafile read
passes `fileEncoding = "UTF-8-BOM"`. A BOM on the genotype file is
harmless by accident (the affected header cell is discarded), but a BOM on
the locus metafile corrupts the `AlleleID` header name, which the
mandatory-column check would then miss. Code-read; not empirically probed
(the loc.metafile path currently crashes earlier, F1).
Proposed change: apply the same `fileEncoding` to all three `read.csv`
calls when F1 is fixed.

**F9 [LOW, confidence: high] — roxygen deviations (DOC2, DOC7 proposed rules; DOC1)**
`R/gl.read.csv.r:38-45` — the `verbose` text ends
"[default 2 or as specified using gl.set.verbosity]" instead of the
canonical DOC2 default clause; the `fbm` param carries no
`[default FALSE]` marker; `@author` names only a custodian with no
`Author(s):` line (DOC7); there is no `@details` tag (the detail
paragraphs sit inside `@description`).
Proposed change: align the header with the DOC1/DOC2/DOC7 templates and
re-document (DOC4).

**F10 [INFO, confidence: high] — half-missing genotypes become NA by coercion accident (DOC5, proposed rule)**
`R/gl.read.csv.r:200-215` — a cell like `A/-` matches none of the four
substitution patterns and falls through `as.numeric()` to NA with a raw
"NAs introduced by coercion" warning (reproduced). The outcome is
defensible; the route is undocumented and the warning uninformative.
Proposed change: document `-/-` as the only missing form and convert
half-missing calls explicitly (to NA, with a gated warning).

### @position / @chromosome (PR #330 interaction)

Verified empirically: `gl.read.csv` returns `@position` and `@chromosome`
NULL — the reader itself never writes either slot, on synthetic fixtures
and on the packaged platy example alike. Correct under the genome-only
convention. On ddaed27 (pre-#330 compliance), a loc.metafile supplying a
`SnpPosition` column would have been copied into `@position` by
`gl.compliance.check` — but the loc.metafile path crashes first (F1), so
no live interaction exists on either commit. When F1 is fixed, the #330
compliance behaviour (no fill, clear stale copies) is the one the fix will
land on; nothing to re-propose.

## Proposed changes

1. Rebuild the loc.metafile validation: name-based `AlleleID` lookup,
   `stop()` on absence, position-wise check against `locNames(gl)`, then
   attach (F1, F4-locus side).
2. Numeric-code guard `all()` with offending values in the message (F2).
   **Consequence: files previously (wrongly) accepted with out-of-range
   codes now error.**
3. Clamp the 5x5 type sniff and the confirmation print to actual
   dimensions (F3).
4. Name-based `id` comparison for ind.metafile (F4).
   **Consequence: metafiles with id not in column 1 now load instead of
   erroring.**
5. Store the true allele pairs from the character branch in `alleles(gl)`
   (F5). **Consequence: loc.all changes from the uniform placeholder to
   real allele spellings; downstream exports change accordingly.**
6. Replace cat+stop() with `stop(error(...))` on all six fatal exits (F6).
7. Fix the loc.metrics report wording and warning helpers (F7).
8. Uniform `fileEncoding` across the three reads (F8).
9. Roxygen alignment: DOC2 verbose text, fbm default tag, DOC7 author
   block, `@details`; re-document (F9).
10. Handle half-missing character calls explicitly and document them (F10).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run (PLT n/a, no plot
  bundle; FS4 n/a, input is a file path not a genlight; DEP n/a for the
  read path — bigsnpr sits in Imports for `fbm = TRUE`; FS8 verified,
  history length 1 on return).
- Spec: both documented encodings run and hand-verified per cell (numeric
  incl. NA; character incl. reversed hets, missing `-/-`, monomorphic
  locus); `transpose = TRUE` round-trip — run; ind.metafile joined/reordered/
  id-misplaced/pop-absent — run; loc.metafile valid/absent-column/short —
  run (all crash, F1); packaged platy example — run; `verbose = 0` text
  silence — verified (zero captured lines; no plot surface exists).
- fbm = TRUE: run (`@fbm` populated, `@gen` cleared, ploidy 2).
- Allele-frequency tie in the character branch (which allele becomes ref
  on an exact tie): SKIPPED — outcome is order-of-table-names dependent;
  noted, not fixtured.
- BOM-corrupted metafile probe (F8): SKIPPED — path unreachable until F1
  is fixed; asserted from code reading.

## Approval (Phase B)

All BLOCKER, HIGH and MEDIUM findings approved by Arthur Georges
2026-09-05 via the formal approval boxes, explicitly acknowledging the
stated consequences (objects users previously read will differ where the
current behaviour is wrong). LOW findings were not approved this round
and are deferred; INFO requires no action.

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur Georges, 2026-09-05 | F1 (HIGH) |
| 2 | approved | Arthur Georges, 2026-09-05 | F2 (HIGH); consequence acknowledged: files previously (wrongly) accepted with out-of-range codes now error |
| 3 | approved | Arthur Georges, 2026-09-05 | F3 (MEDIUM) |
| 4 | approved | Arthur Georges, 2026-09-05 | F4 (MEDIUM); consequence acknowledged: metafiles with id not in column 1 now load |
| 5 | approved | Arthur Georges, 2026-09-05 | F5 (MEDIUM); consequence acknowledged: loc.all changes from the placeholder to real allele spellings |
| 6 | deferred | Arthur Georges, 2026-09-05 | F6 (LOW) not approved this round |
| 7 | deferred | Arthur Georges, 2026-09-05 | F7 (LOW) not approved this round |
| 8 | deferred | Arthur Georges, 2026-09-05 | F8 (LOW) not approved this round |
| 9 | deferred | Arthur Georges, 2026-09-05 | F9 (LOW) not approved this round |
| 10 | no action | Arthur Georges, 2026-09-05 | F10 (INFO) |

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-gl.read.csv` (base upstream/dev
ddaed27): changes 1-5 (F1, F2 HIGH; F3, F4, F5 MEDIUM). Changes 6-10
untouched (LOW deferred / INFO no action).

Verification:

- Baseline characterization suite: 37/37 pass. Five [BUG] pins flipped,
  each mapped to its approved finding ([approved F1]-[approved F5]
  comments in test-gl.read.csv.R); no unexplained diffs. LOW-pinned
  behaviour (cat+stop error exits, message wording, encoding, roxygen)
  unchanged.
- E2e at verbose = 3: numeric fixture with a valid loc.metafile and an
  id-second ind.metafile loads, attaches AlleleID/SnpPosition/RepAvg and
  sets pops B,B,B,C,C,C,C,C; character fixture returns
  loc.all = A/C, T/G, A/A, C/T, G/C (real pairs); packaged platy example
  runs clean and now reports real allele spellings
  (A/T G/C A/T A/A G/C C/T).
- verbose = 0: zero captured lines.
- Position-slot note: with a loc.metafile supplying SnpPosition, the
  ddaed27 compliance pass copies it into `@position` (integer, length
  nLoc) — the pre-#330 behaviour this branch is based on; once PR #330
  merges, compliance leaves the slot NULL per the genome-only
  convention. Nothing to change here, as anticipated in the PR #330
  note above.

```json
{
  "function": "gl.read.csv",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "FS5", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DAT1", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DAT2", "status": "applied", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "deferred", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "deferred", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "medium", "rule": "FS5", "status": "deferred", "change": 8},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "DOC2", "status": "deferred", "change": 9},
    {"id": "F10", "severity": "INFO", "confidence": "high", "rule": "DOC5", "status": "no-action", "change": 10}
  ],
  "coverage_skipped": ["allele-frequency tie in character branch: not fixtured", "BOM metafile probe: unreachable until F1 fixed"],
  "status": "phase-c-applied",
  "pr": null
}
```
