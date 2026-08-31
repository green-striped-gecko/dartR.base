# Review: gl.report.allna (dartR.base)

- Family mode: report
- Date: 2026-08-29
- Reviewer: Claude (claude-fable-5), dartr-function-review v2.0.0
- Package commit: 2739215
- Datasets: testset.gl, testset.gs, plus a crafted 3x4 genlight with one
  all-NA individual (testset objects contain all-NA loci but no all-NA
  individuals, so that branch needed a synthetic fixture)
- Baseline: tests/testthat/test-gl.report.allna.R (snapshot captured pre-review)

## Verdict

**Standards: Needs work** — the datatype check is commented out entirely
(FS4), the results print unconditionally (2 lines at `verbose = 0`; 32
lines with `by.pop = TRUE`), and the docs carry a junk `@return`, a
duplicated-and-wrong second `@family` tag, and a garbled title.

**Spec: Needs work** — the individuals listing is corrupted whenever an
all-NA individual actually exists: it prints a literal "NULL" for every
healthy individual (confirmed: "Individuals that are missing (NA) across
all loci: NULL, IND_ALLNA, NULL"). The loci and by.pop counts are correct
(verified independently), and the report-mode return contract holds
(input returned `identical()`, no history append).

## Findings

**F1 [MEDIUM, confidence: high] — individuals listing prints literal "NULL"s (STY3/Spec)**
`R/gl.report.allna.r:115, 121, 129` — the individuals branch initializes
`ind.list <- vector("list", nI)` (a list of NULLs) where the loci branch
correctly uses `array(NA, nL)`. `ind.list[!is.na(ind.list)]` removes
nothing from a list of NULLs, so `paste()` renders every empty slot as
"NULL". Confirmed on a crafted object with one all-NA individual among
three.
Failure scenario: on any dataset with an all-NA individual, the listing
is unreadable noise — one "NULL" per healthy individual, with the real
name buried among them.
Proposed change: `ind.list <- array(NA, nI)`, mirroring the loci branch,
and print the count like the loci branch does.

**F2 [MEDIUM, confidence: high] — datatype check commented out (FS4/DAT5)**
`R/gl.report.allna.r:61` — the standard
`utils.check.datatype()` line exists but is commented out; nothing
validates the input. A non-genlight input fails downstream with an
obscure error instead of the standard clear fatal message.
Proposed change: restore the check as written in the comment.

**F3 [MEDIUM, confidence: high] — results print unconditionally (VRB5/VRB3)**
`R/gl.report.allna.r:95-104, 127-134, 153, 158-159` — confirmed: 2 lines
at `verbose = 0` (default branch), 32 lines with `by.pop = TRUE`. Same
defect class fixed in #244/#247/#249, gated at `verbose >= 1`.
Proposed change: gate all results output at `verbose >= 1`.

**F4 [LOW, confidence: high] — documentation (DOC1/DOC5)**
`@return` is the junk string "gl.report.allna"; a second `@family filter
functions` tag (line 43) wrongly cross-lists this report among the
filters (the source of the odd Rd cross-links); the `@title` says
"populations with all NA across loci" where it means individuals; the
detail paragraphs sit inside `@description`; tags out of ratified DOC1
order; `@import utils patchwork` imports nothing this function uses.
Proposed change: correct all (docs only).

**F5 [notes, no fix proposed]**
(a) `verbose` doc clause is the "[default 2, unless...]" variant (DOC2
canon note, as all siblings).
(b) The loci listing prints a count, the individuals listing does not —
made symmetric as part of F1.

## Cleared during verification (checked, not a defect)

- Report-mode contract holds: input returned `identical()` on both
  datatypes, no history append.
- All-NA loci listing on testset.gl names exactly the loci an
  independent computation finds; by.pop=TRUE union count matches the
  independent per-population computation.
- No plot bundle — graphics half of VRB5 not applicable.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: counts/listings vs independent computation on testset.gl /
  testset.gs; all-NA-individual branch exercised via a crafted fixture —
  run
- FBM path (DAT6): SKIPPED — no FBM-backed fixture readily available
- Caller grep (for F2/F3 behaviour changes): run in Phase C, results
  recorded in Outcome

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur | listing-content consequence approved |
| 2 | approved | Arthur | fail-fast consequence approved |
| 3 | approved | Arthur | verbose=0 silence consequence approved |
| 4 | approved | Arthur | |

F5 notes only, not offered.

## Outcome

- **F1**: applied — `array(NA, nI)` mirroring the loci branch; listing
  now prints a count and only the affected names. Verified on the
  crafted fixture: one line, names IND_ALLNA, zero "NULL" strings.
- **F2**: applied — datatype check restored. Verified: data.frame input
  fails with the standard "inappropriate object passed" message.
- **F3**: applied — all results gated at `verbose >= 1`. Verified: 0
  lines at `verbose = 0` on both branches (was 2 / 32).
- **F4**: applied — real `@return`, wrong filter-family tag removed,
  title corrected, details moved out of description, unused imports
  dropped, ratified tag order. `man/gl.report.allna.Rd` regenerated.
- Caller grep: zero production call sites in dartR.base or any sibling
  (gl.filter.allna recomputes independently and does not delegate); no
  dartr2shiny. NEWS.md entry added.
- Characterization test: 14/14 pass. `verbose = 3` end-to-end clean on
  both datatypes and both branches.
- PR: [green-striped-gecko/dartR.base#250](https://github.com/green-striped-gecko/dartR.base/pull/250)

```json
{
  "function": "gl.report.allna",
  "package": "dartR.base",
  "family": "report",
  "skill_version": "2.0.0",
  "commit": "2739215",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "STY3", "status": "approved"},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "FS4", "status": "approved"},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "approved"},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved"},
    {"id": "F5", "severity": "LOW", "confidence": "medium", "rule": "DOC2", "status": "note_only"}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture available"],
  "status": "pr-open",
  "pr": 250
}
```
