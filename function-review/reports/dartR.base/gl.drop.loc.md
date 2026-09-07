# Review: gl.drop.loc (dartR.base)

- Family mode: modify
- Date: 2026-08-27
- Reviewer: Claude (claude-fable-5), dartr-function-review v1.6.2
- Package commit: b04c5bb
- Datasets: testset.gl, testset.gs
- Baseline: tests/testthat/test-gl.drop.loc.R (snapshot captured pre-review)

## Verdict

**Standards: Needs work** — a misleading-output bug unique to this
function (the not-present warning names the wrong loci), the two
range-path defects shared with the just-reviewed sibling `gl.keep.loc`
(PR #245), a VRB5 leak (clamp warnings print at `verbose = 0`), and the
family's routine cosmetic items.

**Spec: Needs work** — `@description` says the function "deletes
individuals and their associated metadata" (it deletes loci), the
documented `last` default is unimplemented (confirmed crash), and a
progress message says "loci to keep" on the drop function.

## Findings

**F1 [HIGH, confidence: high] — not-present warning names the wrong loci (Spec / STY)**
`R/gl.drop.loc.r:89-98` — `tmp2` holds indices into `loc.list` of entries
not found in the dataset, but the warning prints `locNames(x)[tmp2]` —
indexing the DATASET's locus names with `loc.list` positions. Confirmed
empirically: `gl.drop.loc(gl, loc.list = c("BOGUS_LOCUS_XYZ", <real>))`
warns that `100049687-12-C/T` — the dataset's locus #1, which IS present
and was never mentioned by the user — is "not present in the dataset",
while never naming `BOGUS_LOCUS_XYZ`. The removal itself
(`loc.list[-tmp2]`) is correct, so the data outcome is right and only the
message lies.
Failure scenario: a user with a typo in one locus name is told a
completely different, valid locus is missing — pointing any debugging at
the wrong locus.
Proposed change: print `loc.list[tmp2]`.

**F2 [HIGH, confidence: high] — `first` without `last` crashes; documented default never implemented (DOC5)**
`R/gl.drop.loc.r:143-148` — identical to `gl.keep.loc` F1 (PR #245):
docs promise `last` "[if not specified, last locus in the dataset]" but
`first:NULL` crashes with "argument is of length zero". Confirmed with
`first = 250`.
Proposed change: `if (is.null(last)) last <- nLoc(x)` at the top of the
range handling, mirroring the sibling fix.

**F3 [MEDIUM, confidence: high] — range clamp tests `first` but clamps `last` (DOC5)**
`R/gl.drop.loc.r:114-123` — identical block to `gl.keep.loc` F3.
Confirmed: `first = nLoc+50, last = nLoc+200` silently drops exactly one
arbitrary locus (254 of 255 remain); only-`last`-out-of-range goes
unclamped and unwarned (correct by accident via NA-name matching).
Proposed change: mirror the sibling fix — clamp on `last > nLoc(x)`,
fatal error for `first > nLoc(x)`.

**F4 [MEDIUM, confidence: high] — missing `invisible()` on return (FS10)**
`R/gl.drop.loc.r:182` — confirmed: unassigned call prints the 22-line
object summary.
Proposed change: `return(invisible(x2))`.

**F5 [MEDIUM, confidence: high] — description and messages say the wrong thing (DOC5)**
Three copy-paste survivals: `@description` (line 7) reads "This function
deletes individuals and their associated metadata" — it deletes loci; the
both-mode progress message (line 66) says "list of loci to keep has been
specified" on the drop function; the DO THE JOB banner comment (line 135)
reads "# Remove individuals ------".
Proposed change: correct all three to loci/drop wording. Docs/messages
only.

**F6 [LOW, confidence: high] — `@return` out of ratified house order (DOC1)**
`R/gl.drop.loc.r:39` — last, after `@export`. Per the DOC1 ratification
(gl.report.bases review), move to directly follow the `@param` block.

**F7 [LOW, confidence: high] — console hygiene: clamp warnings ungated (VRB5/VRB4) + stray `)` + raw `cat()` summary**
`R/gl.drop.loc.r:108-123, 164-170` — the two clamp warnings print at
every verbosity including 0 (confirmed: 2 lines at `verbose = 0`,
including a stray `)` on its own line from the string typo shared with
the sibling); the verbose>=3 summary block uses raw `cat()`.
Proposed change: gate the clamp warnings at `verbose >= 1` (they affect
results, so VRB4's level rather than VRB3's `>= 2`), remove both stray
parens, wrap the summary block in `report()`.

**F8 [LOW, confidence: high] — `@author` states only Custodian (DOC7)**
`R/gl.drop.loc.r:24` — same gap as siblings.
Proposed change: add `Author(s): Arthur Georges.`

**F9 [LOW, confidence: medium] — awareness notes, no fix proposed**
(a) `verbose` param wording is the "progress but not results" variant
(DOC2 canon since v1.6.2: "progress log" + the NULL->global->2 cascade
clause) — same widespread-deviation note as all siblings.
(b) `loc.metrics.flags` not invalidated after locus removal — same
DAT4-adjacent note as `gl.keep.loc` F8; callers (e.g. the
gl.filter.hamming fix earlier in this campaign) treat flag-resetting as
their responsibility.
(c) The just-fixed `gl.keep.loc` (PR #245) has the same ungated clamp
warnings as F7 — worth a one-line follow-up on that PR rather than a fix
here.

## Cleared during verification (checked, not a defect)

- **No-args path**: unlike the sibling pre-fix, this function already
  fails fast with a clear `stop(error(...))` — `gl.keep.loc`'s F2 class
  is NOT present here.
- **DAT2/DAT3**: `loc.metrics` re-subset from the original object with
  the same mask; row count matches `nLoc` in the characterization test.
- **DAT7**: unrestricted `accept=` correct — name/index-based locus
  removal, both datatypes pass.
- **FS8**: history appended on `x2`, the actual returned object, on both
  the normal and empty-list paths.
- **Empty `loc.list` path**: warns (gated `>= 1`, unlike the sibling's
  ungated version) and returns the object unchanged.
- Ploidy, individual counts preserved on both datatypes.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: behaviour vs roxygen on testset.gl / testset.gs, range-parameter
  edge probes — run
- FBM path (DAT6): SKIPPED — no FBM-backed fixture readily available
- Sibling `dartR.*` / dartr2shiny caller grep (API1-3): required for
  F2/F3/F7 if approved (edge-case/console behaviour changes) — pending
  approval

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur | Escalation gate: warning-content change approved |
| 2 | approved | Arthur | Escalation gate: crash -> documented default approved |
| 3 | approved | Arthur | Escalation gate: silent-wrong -> clamped/error approved |
| 4 | approved | Arthur | |
| 5 | approved | Arthur | |
| 6 | approved | Arthur | |
| 7 | approved | Arthur | |
| 8 | approved | Arthur | |
| F9(c) follow-up | approved | Arthur | gl.keep.loc's identical ungated clamp warnings to be gated at verbose >= 1 via a follow-up commit on the open PR #245 branch |

F9(a)/(b) carried no proposed fix (awareness only).

## Outcome

- **F1 (wrong-loci warning)**: applied — warning now prints
  `loc.list[tmp2]`. Verified: the warning names `BOGUS_LOCUS_XYZ` and no
  longer names the dataset's locus #1.
- **F2 (last default)**: applied. Verified: `first = 250` drops loci
  250-255 (6 loci) instead of crashing.
- **F3 (range clamp)**: applied — clamp tests/clamps `last`; `first >
  nLoc` is fatal. Verified both branches.
- **F4 (invisible)**: applied. Verified: unassigned call prints 0 lines
  (was 22).
- **F5 (loci-not-individuals wording)**: applied to `@description`, the
  both-mode progress message, and the banner comment.
- **F6 (@return order)**: applied per ratified DOC1.
- **F7 (console hygiene)**: applied — clamp warnings gated at
  `verbose >= 1`, stray parens removed, summary block `report()`-wrapped.
  Verified: range path at `verbose = 0` prints 0 lines (was 2).
- **F8 (Author(s) line)**: applied per DOC7.
- Escalation gate: caller audit found 14 production call sites (10 in
  dartR.base incl. the gl.filter.* delegates, 4 in dartR.popgen), all
  `loc.list`-based at `verbose = 0` — none pass `first`/`last`, none
  parse console output; zero blast radius. No local dartr2shiny (sixth
  check). NEWS.md entry added. Side-find for a future cleanup:
  `gl.filter.factorloadings.r:135,139` uses a stray in-argument
  assignment (`loclist<-tmp$locus`) that works only by positional
  accident.
- Characterization test: 23/23 pass; the two CHARACTERIZATION probes
  were rewritten as assertions of the fixed behaviour (approved changes,
  not unexplained diffs).
- PR: [green-striped-gecko/dartR.base#246](https://github.com/green-striped-gecko/dartR.base/pull/246)
  - open, based directly on the current dev (post-#238/#239/#240 merges).

```json
{
  "function": "gl.drop.loc",
  "package": "dartR.base",
  "family": "modify",
  "skill_version": "1.6.2",
  "commit": "b04c5bb",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "STY3", "status": "approved"},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved"},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved"},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS10", "status": "approved"},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved"},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved"},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "VRB5", "status": "approved"},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "approved"},
    {"id": "F9", "severity": "LOW", "confidence": "medium", "rule": "DOC2", "status": "note_only"}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture available"],
  "status": "pr-open",
  "pr": 246
}
```
