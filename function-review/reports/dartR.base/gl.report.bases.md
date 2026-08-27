# Review: gl.report.bases (dartR.base)

- Family mode: report
- Date: 2026-08-27
- Reviewer: Claude (claude-fable-5), dartr-function-review v1.6.0
- Package commit: 9344ef4
- Datasets: testset.gl, testset.gs
- Baseline: tests/testthat/test-gl.report.bases.R (snapshot captured pre-review)
- Special scope: this function is the nominated exemplar for the package;
  this review also cross-checks the conventions catalogue against it
  (meta-findings M1-M3, reported separately from the function findings).

## Verdict

**Standards: Needs work** — one HIGH verbosity-gating defect (the results
table prints even at `verbose = 0`), plus three LOW cosmetic items. The
report-mode core contracts all hold: input returned `identical()`, no
history append, results independent of plotting, plot suppressed at
`verbose = 0`.

**Spec: Ready** — base frequencies verified against an independent
computation; ts/tv reported for SNP and correctly withheld for SilicoDArT;
docs match behaviour.

## Findings (function)

**F1 [HIGH, confidence: high] — results table prints unconditionally, ignoring `verbose` (VRB5/VRB3)**
`R/gl.report.bases.r:189-200, 215-219` — the "Printing outputs" block
(sequence-length line, sequence count, base-frequency table) and the SNP
ts/tv lines are not gated on `verbose` at all. Confirmed empirically: at
`verbose = 0`, 11 lines print. The graphics half is handled correctly
(`verbose == 0` forces `plot.display <- FALSE`; confirmed no device is
touched), so only the text side violates VRB5. This is the same defect
class the skill's trial run found on `gl.report.callrate`.
Failure scenario: any scripted/pipeline call at `verbose = 0` gets 11
lines of unwanted console output per invocation.
Proposed change: gate both blocks on `verbose`. Two candidate levels:
`>= 3` (strict VRB1 — 3 is "results summary"; consequence: the default
`verbose = 2` shows progress but no results table, though the plot still
displays) or `>= 1` (minimal VRB5 compliance; default behaviour unchanged,
only fully-quiet calls go silent). Member's choice.

**F2 [LOW, confidence: high] — `@author` punctuation: `Author(s);` not `Author(s):` (DOC7)**
`R/gl.report.bases.r:59` — reads `Author(s); Arthur Georges. Custodian:
...` — a semicolon where DOC7's format has a colon. Both parts are present
(the exemplar supports DOC7); only the punctuation is off.
Proposed change: `Author(s):`.

**F3 [LOW, confidence: high] — one unqualified `str_count` call (DEP2)**
`R/gl.report.bases.r:159` — `sum(str_count(state.change, "A>C"))` is the
only one of 12 `str_count` calls not qualified as `stringr::str_count`. It
works (NAMESPACE has `import(stringr)`) but is inconsistent with the other
11 calls in the same expression.
Proposed change: add the `stringr::` prefix.

**F4 [LOW, confidence: high] — three messages embed source-indentation whitespace (VRB2-adjacent)**
`R/gl.report.bases.r:112-115, 176-182, 206-212` — the TrimmedSequence
fatal error and two warnings are multi-line string literals, so the
printed message contains a line break plus ~20 spaces of source
indentation mid-sentence.
Proposed change: reflow each onto a single-line string.

**F5 [INFO] — variable named `T` shadows the TRUE shorthand (STY3)**
`R/gl.report.bases.r:128` — `T <- sum(...)` masks the base-R `T` alias
within the function. No concrete failure exists (no bare `T`-as-TRUE is
used afterwards), so per the report contract this is a note, not a
finding. No change proposed.

## Meta-findings (conventions catalogue vs the exemplar)

The member nominated `gl.report.bases` as the model function and asked
that the skill check conformity with it. Walking the 50-rule catalogue
against the exemplar:

**M1 — DOC1's tag order contradicts the exemplar on two points.**
DOC1 [confirmed] prescribes `... @description, @details, @param, @return,
@author, @examples, @export`. The exemplar's actual order is `@name,
@title, @family, @description, @param, @details, @return-nowhere-early...
@author, @examples, @export, @return` — that is: `@param` BEFORE
`@details`, and `@return` LAST, after `@export`. The sandpit style guide
independently derived from this same exemplar agrees with the exemplar on
both points. DOC1's claimed order matches neither.
**Consequence needing a decision:** four earlier reviews (gl.drop.pop,
gl.drop.ind, gl.keep.pop, gl.keep.ind — PRs #237 merged, #238/#239/#240
open) each applied a "move @return to follow @param" fix citing DOC1.
Those fixes conform to DOC1 but CONTRADICT the exemplar. Either DOC1 is
amended to the exemplar's order and those four reorders should be
reverted in their open PRs (and noted for the merged one), or the team
ratifies DOC1's order and the exemplar itself is out of date. Not the
reviewer's call.

**M2 — DOC2's canonical `verbose` text matches neither the exemplar nor any observed variant.**
DOC2 [confirmed] canon: "2, brief progress messages ... [default 2,
unless specified using gl.set.verbosity]". The exemplar: "2, progress log
... [default NULL, unless specified using gl.set.verbosity]". The
drop/keep siblings: "2, progress but not results ... [default 2 or as
specified using gl.set.verbosity]". At least three wordings exist in the
codebase and DOC2's canon matches none of the surveyed ones. Proposed:
adopt the exemplar's scale wording ("2, progress log") as canonical; the
default clause is a separate choice ("[default 2, ...]" is factually what
the NULL resolves to, while the exemplar literally says "[default NULL,
...]").

**M3 — the skill has no exemplar-conformity provision.**
Nothing in SKILL.md or the catalogue names the exemplar or requires
checking rules against it. Proposed: record `gl.report.bases` as the
nominated exemplar in the conventions preamble; the standards walk checks
rule-vs-exemplar conformity; where a rule and the exemplar diverge, the
divergence is raised to the team (flag the rule, per the existing
false-positive discipline) rather than silently enforcing either — noting
the exemplar can itself be defective (F1 here is a live example, since
the team's own VRB5 mandate outranks the exemplar's current behaviour).

## Cleared during verification (checked, not a defect)

- **Report-mode contracts**: returned object `identical()` to input on
  both datatypes; no history append; results computed and printed
  independent of `plot.display` (PLT3); `utils.plot.save()` decoupled
  from display.
- **VRB5 graphics half**: `verbose = 0` forces `plot.display <- FALSE`;
  confirmed no graphics device is touched (no Rplots.pdf in a fresh
  Rscript session).
- **Spec numerical check**: reported A% matches an independent
  base-counting computation on `TrimmedSequence` exactly; ts/tv printed
  for SNP, correctly absent for SilicoDArT.
- **DAT7**: `accept = c("genlight", "SNP", "SilicoDArT")` is explicit and
  correct — the function genuinely serves both datatypes, and the
  TrimmedSequence guard covers the real precondition.
- **PLT1/PLT2**: full modern plot bundle (plot.theme, plot.colors,
  patchwork, plot.file/plot.dir + utils.plot.save) — the exemplar is the
  PLT2-recommended idiom.
- **FS walk**: preamble order, `gl.check.wd`, flag-start, explicit accept,
  `invisible(x)` return — all conform.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: behaviour vs roxygen on testset.gl / testset.gs, base-frequency
  numbers verified against an independent computation — run
- Exemplar conformity cross-check (catalogue vs this function) — run
  (M1-M3)
- FBM path (DAT6): SKIPPED — no FBM-backed fixture readily available
- Sibling `dartR.*` / dartr2shiny caller grep (API1-3): required for F1
  if approved (console-output change) — pending approval

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| F1 | approved, gate at `verbose >= 1` | Arthur | Escalation gate: console-change consequence approved; the stricter `>= 3` option (which would silence the default) was declined |
| F2 | approved | Arthur | |
| F3 | approved | Arthur | |
| F4 | approved | Arthur | |
| F5 | note only | — | not offered |
| M1 | resolved: **DOC1 wins** | Arthur | Catalogue order ratified as the standard; the four applied @return reorders stand; the exemplar's own tag order becomes an approved fix in this review (recorded as F6 below) |
| M2 | approved | Arthur | DOC2 canon becomes the exemplar's scale wording ("2, progress log") with "[default 2, unless specified using gl.set.verbosity]" |
| M3 | rejected | Arthur | No exemplar-conformity rule added — consistent with M1's ruling that the catalogue, not the exemplar, is the standard of record |

**F6 (from M1's ruling)** — reorder the exemplar's own roxygen tags to the
ratified DOC1 order: `@details` before `@param`, `@return` after the
`@param` block (currently last, after `@export`). Docs-only.

## Outcome

- **F1 (gate results printout at `verbose >= 1`)**: applied at
  `R/gl.report.bases.r` (the "Printing outputs" block and the SNP ts/tv
  lines). Verified: `verbose = 0` with default arguments now prints 0
  lines and touches no graphics device; `verbose = 2` (default) output
  and plot unchanged; `verbose = 3` end-to-end clean. Escalation gate:
  full-text grep of all 7 sibling checkouts found zero external callers
  of `gl.report.bases()` (only in-package examples and Rd cross-links);
  no local dartr2shiny (fourth consecutive negative). `NEWS.md` entry
  added.
- **F2 (`Author(s):` punctuation)**: applied.
- **F3 (qualify `str_count`)**: applied — all 12 calls now
  `stringr::`-qualified.
- **F4 (reflow 3 broken strings)**: applied — fatal error and both
  warnings now single-line.
- **F6 (from M1: reorder exemplar's tags to ratified DOC1 order)**:
  applied — `@details` now precedes `@param`, `@return` follows the
  `@param` block. `man/gl.report.bases.Rd` regenerated (the `\author{}`
  punctuation is the visible change; Rd section order is fixed by
  roxygen2 regardless of source order).
- **M1/M2 (skill)**: applied in dartrverse-skills v1.6.1 (local commit
  `7f05b17`); DOC1 carries the ratification note, DOC2's canon is now
  the exemplar's scale wording. Notified upstream as
  mijangos81/dartrverse-skills#3. M3 rejected — no change.
- **M2 addendum (v1.6.2, commit `650961c`)**: Arthur queried v1.6.1's
  "[default 2, ...]" clause. Verified against `gl.check.verbosity()`
  source and five empirical probes: the true resolution is NULL ->
  global (`options(dartR_verbose)` via `gl.set.verbosity()`) -> 2 if no
  global; explicit argument always overrides. DOC2's default clause now
  states this cascade, and the exemplar's `@param verbose` was updated
  to the corrected canon in this same review (its old text's "unless"
  clause attached misleadingly). Issue #3 comment records the
  correction. The end-to-end probe also confirmed the F1 fix holds
  through the global-verbosity path (global=0, `verbose` omitted ->
  fully silent).
- Characterization test (`tests/testthat/test-gl.report.bases.R`): 10/10
  pass. The two output-parsing tests moved from `verbose = 0` to
  `verbose = 1` and the VRB5 characterization test now asserts silence
  at 0 (approved behaviour change, not an unexplained diff).
- PR: not yet opened (pending pre-push confirmation).

```json
{
  "function": "gl.report.bases",
  "package": "dartR.base",
  "family": "report",
  "skill_version": "1.6.0",
  "commit": "9344ef4",
  "verdict_standards": "needs_work",
  "verdict_spec": "ready",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "VRB5", "status": "approved"},
    {"id": "F2", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "approved"},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "DEP2", "status": "approved"},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "approved"},
    {"id": "F5", "severity": "INFO", "confidence": "high", "rule": "STY3", "status": "note_only"},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved"},
    {"id": "M1", "severity": "HIGH", "confidence": "high", "rule": "DOC1", "status": "resolved_doc1_ratified", "meta": true},
    {"id": "M2", "severity": "MEDIUM", "confidence": "high", "rule": "DOC2", "status": "approved", "meta": true},
    {"id": "M3", "severity": "MEDIUM", "confidence": "high", "rule": "SKILL", "status": "rejected", "meta": true}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture available"],
  "status": "pr-open",
  "pr": null
}
```
