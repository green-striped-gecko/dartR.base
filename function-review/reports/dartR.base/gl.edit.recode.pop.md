# Review: gl.edit.recode.pop (dartR.base)
- Family mode: modify (interactively recodes population labels; deletes flagged pops; optional recalc/mono.rm)
- Date: 2026-09-15
- Reviewer: Claude (Claude Opus 4.8), dartr-function-review v1.0.0
- Package commit: f02bd34 (upstream/dev; working copy verified identical by `git diff upstream/dev -- R/gl.edit.recode.pop.r`)
- Datasets: testset.gl (SNP), testset.gs (SilicoDArT) — dartR.data 1.2.5
- Baseline: tests/testthat/test-gl.edit.recode.pop.R (new file, snapshot captured pre-review; 12 assertions, all passing)

## Verdicts

**Standards: Needs work** — the FS backbone is present (and mono.rm/recalc
are correctly gated to SNP data here, unlike the individual sibling), but
the metric flags depend on the verbosity level and a missing metrics flag
crashes the run.

**Spec: Rework** — the `pop.recode` parameter, documented (and used
everywhere else in the recode family) as the input recode table, is neither
read nor applied; instead the function **overwrites** it with a freshly
generated identity table, destroying a user's existing recode file. The
documented output parameter `out.recode.file` is dead, and `outpath` is
ignored. The parameter contract has to be redefined, not patched.

## Blast radius

`gl.edit.recode.pop` has **no internal callers** anywhere in the
dartRverse — a user-facing interactive function (it calls `edit()`).

## Independent verification (spec axis)

`edit()` no-ops non-interactively, so the surrounding logic was exercised
directly (tests below). The `pop.recode` semantics were checked against the
non-interactive sibling `gl.recode.pop`, which reads `pop.recode` with
`read.csv(pop.recode)` and applies it — confirming the family-standard
meaning of the parameter (input recode table).

## Findings

**F1 [HIGH, confidence: high] — pop.recode input file is silently overwritten; documented input/output params do not work (DAT2, DOC5 (proposed rule); data-loss)**
`R/gl.edit.recode.pop.r:102,121,128-147` — `pop.recode` is documented as
"Path to recode file" (the input table to edit), and the whole recode
family reads it as input. Here it is never read: the recode table is always
regenerated from `levels(pop(x))` (line 121). Worse, the write branch
(line 128-147) writes the freshly generated table to `file = pop.recode`,
**overwriting** the user's existing recode file with an identity mapping.
Meanwhile `out.recode.file` (the documented output name) is only used to
build `outfilespec` (line 102), which is never used, and `outpath` is
therefore ignored.
Failure scenario: a user follows the docs and passes their carefully
prepared recode table via `pop.recode`; the function ignores its contents
(populations are not recoded as intended) and overwrites the file with
`old,old` identity pairs — the recode table is destroyed. Confirmed
(test 4): the existing table's "MERGED" mapping is not applied and the file
is left containing identity pairs; `out.recode.file` in `outpath` is never
written.
Proposed change: restore the documented, family-consistent contract —
read and apply `pop.recode` as the input table when supplied (else generate
from `x`), and write the edited table to `outfilespec` (from
`out.recode.file` + `outpath`). This is a parameter-semantics change (see
Consequence).
**Consequence: `pop.recode` becomes a read-only input (no longer
overwritten) and `out.recode.file` becomes the output written to `outpath`.
Callers who currently pass `pop.recode` to receive the output file must
switch to `out.recode.file`.**

**F2 [HIGH, confidence: high] — metric-validity flags depend on the verbosity level (DAT4; verbosity-dependent-results class)**
`R/gl.edit.recode.pop.r:208-213` — for SNP data with `recalc = FALSE` (the
default), `x <- utils.reset.flags(x, verbose = 0)` — which marks the metrics
stale after populations are deleted — sits inside the `if (verbose >= 2)`
block, so the flags are reset only at `verbose >= 2`.
Failure scenario: `gl.edit.recode.pop(x, recalc = FALSE)` at `verbose = 0`
returns an object whose `loc.metrics.flags` still read TRUE (metrics claim
current) while the same call at `verbose = 2` returns them reset to FALSE.
Confirmed (test 1). Identical to the individual-recode sibling (PR #399).
Proposed change: move `utils.reset.flags()` out of the verbose gate.

**F3 [MEDIUM, confidence: high] — a missing monomorphs flag crashes the run (DAT5)**
`R/gl.edit.recode.pop.r:192` —
`if (x@other$loc.metrics.flags$monomorphs == FALSE)` evaluates
`if (logical(0))` when the flag is absent. Confirmed (test 2).
Proposed change: `if (!isTRUE(x@other$loc.metrics.flags$monomorphs))`.

**F4 [MEDIUM, confidence: high] — documented defaults contradict the signature (DOC5 (proposed rule))**
`R/gl.edit.recode.pop.r:19-21` document `recalc`/`mono.rm` as
`[default TRUE]`; the signature (lines 83-84) sets both `FALSE`. Confirmed
via `formals()` (test 3).
Proposed change: correct the documentation to `[default FALSE]` (signature
default preserved).

**F5 [LOW, confidence: high] — roxygen copy-paste errors (DOC1, DOC2)**
`R/gl.edit.recode.pop.r:15` — `@param out.recode.file "Name of the file to
output the new individual labels"` (should be population assignments);
`:17` `@param outpath "Directory to save the plot RDS files"` (this writes a
recode CSV); `:22-24` the `verbose` text is the pre-DOC2 form; `:26` the
`@details` tag is indented (`#'  @details`) so roxygen may not parse it.
Proposed change: correct the parameter descriptions, verbose text and the
`@details` indentation.

## Proposed changes

1. Restore the `pop.recode` (input) / `out.recode.file` (output, via
   `outfilespec`/`outpath`) contract: read and apply `pop.recode` when
   supplied instead of overwriting it; write the edited table to
   `outfilespec` (F1).
   **Consequence: `pop.recode` is read-only input; `out.recode.file` is the
   output. Callers passing `pop.recode` to obtain the output must move to
   `out.recode.file`.**
2. Move `utils.reset.flags()` out of the `if (verbose >= 2)` block (F2).
3. Null-safe monomorphs-flag check (F3).
4. Correct the documented `recalc`/`mono.rm` defaults to FALSE (F4).
5. Roxygen fixes: parameter descriptions, verbose text, `@details`
   indentation (F5).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP (n/a), PLT (n/a), STY — run
- Spec: pop.recode input/overwrite behaviour, out.recode.file/outpath
  handling, verbose-dependent flag state, missing-flag crash, defaults,
  no-population error, history — run; `pop.recode` semantics cross-checked
  against `gl.recode.pop`
- Interactive `edit()` path: SKIPPED for unit testing — no-ops
  non-interactively; the recode/delete application loop inspected
- Caller survey (API3): no internal callers; user-facing interactive
  function
- GitHub issues / Google Group: SKIPPED — not queried this session

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur | restore documented contract: pop.recode = read-only input, out.recode.file = output to outpath. Consequence (callers passing pop.recode for output must switch) approved. |
| 2 | approved | Arthur | |
| 3 | approved | Arthur | |
| 4 | approved | Arthur | docs corrected to FALSE (defaults unchanged) |
| 5 | approved | Arthur | |

## Outcome

(pending Phase C)

```json
{
  "function": "gl.edit.recode.pop",
  "package": "dartR.base",
  "family": "modify",
  "skill_version": "1.0.0",
  "commit": "f02bd34",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT2", "status": "proposed", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DAT4", "status": "proposed", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DAT5", "status": "proposed", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "proposed", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "proposed", "change": 5}
  ],
  "coverage_skipped": ["interactive edit() path not unit-tested", "GitHub issues not queried", "Google Group not queried"],
  "status": "awaiting-approval",
  "pr": null
}
```
