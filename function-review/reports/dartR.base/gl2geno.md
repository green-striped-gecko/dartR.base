# Review: gl2geno (dartR.base)

- Family mode: io
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203
- Datasets: testset.gl, testset.gs
- Baseline: tests/testthat/test-gl2geno.R (snapshot captured pre-review)

## Verdict

**Standards: Ready** — structure conforms (verbosity, flag start/end, datatype
check, `gl.check.wd` for output); the fixes below are trivial deletions and
message corrections.
**Spec: Ready** — both output formats verified correct against the LEA layout:
`.geno` is one row per locus with the 0/1/2 alphabet and 9 as the missing
code, `.lfmm` is one space-separated row per individual, and values match
`as.matrix(x)` with `NA -> 9`. SilicoDArT writes the 0/1/9 alphabet as
documented. All-NA loci are removed before writing, with a warning gated at
`verbose >= 1` (VRB4-compliant).

## Findings

**F1 [LOW, confidence: high] — dead code block (STY1)**
`R/gl2geno.r:65-156` — the former conversion algorithm survives as ~90 lines
of commented-out code between the live preamble and the live 2-line
implementation.
Failure scenario: a maintainer reads or edits the dead block believing it is
live; the dead block also documents behaviour (monomorph removal) the live
code no longer performs, inviting wrong conclusions about what the function
does.
Proposed change: delete the commented-out block; git history preserves it.

**F2 [LOW, confidence: high] — garbled output-file message (DOC5, proposed rule)**
`R/gl2geno.r:177-181` — at `verbose >= 1` the function prints
`Output files: gl_geno.geno.lfmm.` — the two extensions concatenated onto one
phantom filename.
Failure scenario: a user looks for a file named `gl_geno.geno.lfmm` and
concludes the export failed; neither actual filename nor the output path is
shown.
Proposed change: print the two real paths, e.g.
`paste0(outfilespec, ".geno")` and `paste0(outfilespec, ".lfmm")`.

**F3 [LOW, confidence: high] — roxygen deviations (DOC1, DOC2)**
`R/gl2geno.r:1-29` — no `@details` tag; `@return` sits after `@export`
(pre-2026-08-27 order); the `verbose` default clause reads
"[default 2, unless specified using gl.set.verbosity]" rather than the
canonical cascade wording.
Failure scenario: none at runtime; the manual page omits the house-standard
default explanation and the header drifts from the ratified tag order.
Proposed change: add a short `@details` (note the all-NA removal and the
missing-data code), move `@return` before `@author`, adopt the DOC2 clause.

**F4 [LOW, confidence: high] — author block lacks Author(s) line (DOC7) (proposed rule)**
`R/gl2geno.r:19-20` — `@author` names only a custodian.
Failure scenario: authorship is unattributable from the manual page.
Proposed change: add `Author(s): Luis Mijangos.` per the DOC7 default.

## Proposed changes

1. Delete the commented-out former implementation, lines 65-156 (F1).
2. Replace the output message with the two real file paths (F2).
3. Roxygen tidy-up: add `@details`, reorder `@return`, canonical DOC2
   `verbose` wording, add `Author(s):` line (F3, F4).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run (no plot bundle;
  no Suggests dependency, so DEP1 not applicable; FS7 satisfied via
  `gl.check.wd`).
- Spec: `.geno`/`.lfmm` layout, alphabet, missing code, and value-for-value
  agreement with `as.matrix()` verified on testset.gl; SilicoDArT path run
  on testset.gs (0/1/9) — run.
- `verbose = 0` silence: verified empirically (zero captured lines) — run.
- Parse-back through `LEA::read.geno`/`read.lfmm`: SKIPPED — LEA not
  installed on the review machine; layout checked against the LEA format
  description instead.
- FBM path (DAT6): SKIPPED — no FBM fixture for this function.
- Note (not a finding): the dead code used to strip monomorphic loci
  generated during conversion; the live code retains monomorphic loci. No
  documented claim is violated.

## Approval (Phase B)

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges | 2026-09-05, via approval boxes |
| 2 | Approved | Arthur Georges | 2026-09-05, via approval boxes |
| 3 | Approved | Arthur Georges | 2026-09-05, via approval boxes |

All findings at every severity approved 2026-09-05 (blanket class
approval, including the acknowledged consequence that converted outputs
change where they were wrong, and the DAT7 fatal `accept = "SNP"` gate
wherever SilicoDArT is silently admitted — not applicable here:
`gl2geno` legitimately supports both datatypes).

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-gl2geno` (base `upstream/dev`,
ddaed27). PR: (recorded below after opening).

- F1: commented-out former implementation (~90 lines) deleted.
- F2: verbose >= 1 message now prints the two real output paths
  (`<outfilespec>.geno`, `<outfilespec>.lfmm`).
- F3/F4: `@details` added (all-NA removal, missing code 9), `@return`
  moved to the house position, DOC2 verbose clause adopted, `Author(s):`
  line added; re-documented.

Verification: all 13 baseline assertions pass; the only expectation
updated is the output-file message test, per F2. End-to-end run on
`testset.gl` at `verbose = 3` prints the corrected paths and completes.
Sibling caller grep across the 8 clones: one caller
(`dartR.popgen::gl.run.snmf`) — API and file output unchanged, no
impact. NEWS entry added for the message fix.

Note for a future review: `R/gl2gapit.r` declares `@name gl2geno`, so
`man/gl2geno.Rd` is a merged page (pre-existing; same defect class as
`gi2gl` F3). Out of scope for this branch.

```json
{
  "function": "gl2geno",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203",
  "verdict_standards": "ready",
  "verdict_spec": "ready",
  "findings": [
    {"id": "F1", "severity": "LOW", "confidence": "high", "rule": "STY1", "status": "applied", "change": 1},
    {"id": "F2", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "applied", "change": 3}
  ],
  "coverage_skipped": ["LEA parse-back: LEA not installed", "DAT6: no FBM fixture"],
  "status": "phase-c-complete",
  "pr": null
}
```
