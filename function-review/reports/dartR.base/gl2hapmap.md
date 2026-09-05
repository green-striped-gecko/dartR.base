# Review: gl2hapmap (dartR.base)

- Family mode: io
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203
- Datasets: platypus.gl, testset.gl, testset.gs, adegenet::glSim fixture
- Baseline: tests/testthat/test-gl2hapmap.R (snapshot captured pre-review)

## Verdict

**Standards: Needs work** — the output message is not verbosity-gated, the
datatype rejection uses a broken stop idiom, and a genlight without
`loc.all` dies with an opaque error.
**Spec: Needs work** — the file itself is well-formed HapMap (verified: the
standard 11 metadata columns headed `rs#`, `A/B` allele definitions,
two-letter genotype pairs with `NN` for missing, one row per locus, one
column per sample), but the documented `pos` and `chrom` arguments are
silently ignored whenever the corresponding slot is already populated, and
the position column's fallback behaviour is undocumented.

## Findings

**F1 [MEDIUM, confidence: high] — documented pos/chrom arguments silently ignored (DOC5, proposed rule)**
`R/gl2hapmap.r:95-117,120-146` — the `pos` branch runs only when
`x$position` is NULL or mis-sized, and the `chrom` branch only when
`x$chromosome` is NULL; an explicit argument is otherwise discarded without
any message.
Failure scenario: verified — with `@position` populated,
`gl2hapmap(x, pos = "SnpPosition")` writes the slot values and ignores the
argument entirely; same for `chrom` with `@chromosome` set. A user
nominating a genome-coordinate field gets a file whose positions came from
somewhere else, with no indication. Notably, `platypus.gl` (the function's
own example) ships with `@position` pre-populated with within-tag offsets,
so the documented example emits tag offsets as genomic positions against
`chrom = 0` (see report `audit-position-slot.md`).
Proposed change: let an explicitly supplied argument take precedence over
the slot (or error on the conflict); document the precedence in `@param`.

**F2 [MEDIUM, confidence: high] — output message and behaviour not silent at verbose = 0 (VRB5)**
`R/gl2hapmap.r:230-233` — `cat(report("  The hapmap file is saved as: ...`
is not gated on `verbose`.
Failure scenario: verified — `capture.output()` around a `verbose = 0` call
returns the saved-file line (plus "NULL", see F5); pipelines that rely on
silent execution get console noise.
Proposed change: gate at `verbose >= 2`.

**F3 [MEDIUM, confidence: high] — opaque crash on a genlight without loc.all (DAT5)**
`R/gl2hapmap.r:159-179` — `homs1`/`hets`/`homs2` are built from `x@loc.all`
with no NULL check, then indexed per column.
Failure scenario: verified — a `glSim` genlight (or any object not built by
dartR's readers) fails with "number of items to replace is not a multiple
of replacement length", far from the actual cause.
Proposed change: fail fast in FUNCTION SPECIFIC ERROR CHECKING with a clear
message when `x@loc.all` is NULL (allele definitions are mandatory for
HapMap output).

**F4 [MEDIUM, confidence: medium] — position fallback undocumented; silent all-zero pos column (DOC5, proposed rule)**
`R/gl2hapmap.r:95-117` — when `@position` is NULL and `pos` is not given,
every SNP gets position 0 with no message at any verbosity (the analogous
chromosome fallback does print at `verbose >= 2`).
Failure scenario: verified on testset.gl — the emitted file has `pos = 0`
and `chrom = 0` throughout; the roxygen description promises "the standard
columns for ... position" and discloses the 0 convention only for
chromosome, so a consumer (e.g. TASSEL) receives coincident zero coordinates
for every marker without the docs warning anyone.
Proposed change: mirror the chromosome path (message at `verbose >= 2`) and
state the zero-fill convention for `pos` in `@description`/`@details`.

**F5 [LOW, confidence: high] — visible NULL return (FS10)**
`R/gl2hapmap.r:242` — `return(NULL)` is visible, so every top-level call
prints `NULL` after the output message.
Failure scenario: verified in the captured output.
Proposed change: `invisible(NULL)` (matches `gl2geno`).

**F6 [LOW, confidence: high] — broken stop idiom on the SilicoDArT rejection (FS5, VRB2)**
`R/gl2hapmap.r:80-85` — `cat(error(...)); stop()` prints the message but
raises a condition whose message is empty.
Failure scenario: verified — `tryCatch` callers receive
`conditionMessage(e) == ""`; logs and test harnesses record an error with no
text.
Proposed change: `stop(error("..."))`; equivalently, pass
`accept = "SNP"` to `utils.check.datatype` and drop the manual block.

**F7 [LOW, confidence: high] — message() used instead of the house cat(report()) idiom (VRB2)**
`R/gl2hapmap.r:109-112,125-129,138-141` — progress notes go through
`message(report(...))`, sending crayon-coloured text to stderr, unlike every
other message in the function.
Failure scenario: `suppressMessages()` hides these notes but not the
`cat()`-based ones; output routing is inconsistent within one function.
Proposed change: use `cat(report(...))` gated at `verbose >= 2`.

**F8 [INFO, confidence: medium] — lexicographic chromosome sort; non-ASCII hyphens (STY1)**
`R/gl2hapmap.r:215` — rows sort by the chromosome factor's alphabetical
levels, so "chr10" precedes "chr2" for character-named chromosomes; lines
96 and 155 contain U+2010 hyphens in comments.
Failure scenario: cosmetic ordering only; HapMap consumers do not require a
particular row order.
Proposed change: none required; note for a future touch (natural-order sort,
ASCII hyphens).

## Proposed changes

1. Give explicit `pos`/`chrom` arguments precedence over populated slots and
   document the precedence (F1).
   **Consequence: output positions change for callers who passed `pos`/
   `chrom` against an object with populated slots — previously the argument
   was ignored.**
2. Gate the saved-file message at `verbose >= 2` and return
   `invisible(NULL)` (F2, F5).
3. Fail fast with a clear message when `x@loc.all` is NULL (F3).
4. Document and message the zero-fill position fallback (F4).
5. Replace `cat(error()); stop()` with `stop(error(...))` or
   `accept = "SNP"` (F6).
6. Convert `message(report())` calls to gated `cat(report())` (F7).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run (FS7 satisfied via
  `gl.check.wd`; no Suggests dependency; no plot bundle; SNP-only dispatch
  present — explicit rejection satisfies the io-family DAT7 check, flawed
  idiom reported as F6).
- Spec: file layout, header, allele column, genotype alphabet, and missing
  code verified on platypus.gl; `@position`-NULL fallback verified on
  testset.gl; slot-vs-argument precedence verified both ways; loc.all-NULL
  crash verified on glSim; SilicoDArT rejection verified — run.
- `verbose = 0` silence: run — FAILS (F2, F5).
- Round trip / re-import into a HapMap consumer (TASSEL): SKIPPED — no
  external toolchain on the review machine; column set checked against the
  HapMap column standard.
- FBM path (DAT6): SKIPPED — no FBM fixture.

## Approval (Phase B)

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges | 2026-09-05, via approval boxes |
| 2 | Approved | Arthur Georges | 2026-09-05, via approval boxes |
| 3 | Approved | Arthur Georges | 2026-09-05, via approval boxes |
| 4 | Approved | Arthur Georges | 2026-09-05, via approval boxes |
| 5 | Approved | Arthur Georges | 2026-09-05, via approval boxes |
| 6 | Approved | Arthur Georges | 2026-09-05, via approval boxes |

All findings at every severity approved 2026-09-05 (blanket class
approval, explicitly acknowledging that converted outputs change where
they were wrong, and including the DAT7 fatal `accept = "SNP"` gate —
applied here as the F6 resolution, replacing the empty-message manual
rejection). F8 (INFO) requires no change, as scoped.

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-gl2hapmap` (base `upstream/dev`,
ddaed27). PR: green-striped-gecko/dartR.base#351.

- F1: explicitly nominated `pos`/`chrom` loc.metrics fields now take
  precedence over populated slots; precedence documented in `@param`.
- F2/F5: saved-file message gated at `verbose >= 2`; return is
  `invisible(NULL)` — `verbose = 0` is now fully silent (verified,
  zero captured lines).
- F3: fatal error with a clear message when `loc.all` is NULL.
- F4: the zero-fill position convention is documented in `@description`
  and messaged at `verbose >= 2`, mirroring the chromosome path.
- F6: SilicoDArT rejection now goes through
  `utils.check.datatype(accept = "SNP")`; the condition carries an
  informative message (was empty).
- F7: `message(report())` calls converted to gated `cat(report())`.
- F8 (INFO): no change, per the report.

Verification: all 21 baseline assertions pass. Updated expectations,
each mapped to a finding: pos/chrom argument precedence (F1),
verbose-0 silence and invisible return (F2, F5), loc.all fail-fast
message (F3), SilicoDArT condition message (F6). One further updated
expectation is upstream baseline drift, not a review change: the Phase A
baseline (ed99203) emitted platypus.gl's populated tag-offset
`@position` verbatim; the PR base (ddaed27, via PR #293) already
zero-fills positions whose maximum is < 1000 — verified on the pristine
base before applying any change. This branch matches the base for that
scenario and adds the F4 documentation/message for it.
End-to-end `gl2hapmap(testset.gl, verbose = 3)` completes with gated
zero-fill messages. Sibling caller grep across the 8 clones: no callers
of `gl2hapmap` — all-clear. NEWS entry added.

```json
{
  "function": "gl2hapmap",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DAT5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "medium", "rule": "DOC5", "status": "applied", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "FS10", "status": "applied", "change": 2},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "FS5", "status": "applied", "change": 5},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "applied", "change": 6},
    {"id": "F8", "severity": "INFO", "confidence": "medium", "rule": "STY1", "status": "no-change", "change": null}
  ],
  "coverage_skipped": ["TASSEL re-import: no external toolchain", "DAT6: no FBM fixture"],
  "status": "phase-c-complete",
  "pr": 351
}
```
