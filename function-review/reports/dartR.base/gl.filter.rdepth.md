# Review: gl.filter.rdepth (dartR.base)

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.filter.rdepth
- Datasets: testset.gl (SNP, rdepth), testset.gs (SilicoDArT, AvgReadDepth)
- Family mode: modify
- Checks skipped: FBM path not exercised (not available: no FBM fixture);
  Google Group not searched (not available: no browser session — GitHub
  issues showed no open complaint naming this function).

## Verdicts

- **Standards: FAIL** — visible return, missing `verbose == 0 →
  plot.display = FALSE` guard, `@family` tag indented so it leaks into the
  rendered help title instead of being parsed, header breaks
  DOC1/DOC2/DOC7.
- **Spec: FAIL** — loci with NA read depth silently desync genotypes from
  the locus metadata (127 loci vs 137 metric rows with 10 NAs injected),
  and `plot.file` with `plot.display = FALSE` crashes. On clean data the
  filter is numerically correct: retained counts match independent
  recomputation on both datatypes and boundary loci are retained inclusive
  as documented — though the progress message states the opposite
  (`<=`/`>=` where the behaviour is `<`/`>`).

## Findings

### F1 — NA read depths silently desync genotypes and metadata (HIGH, confidence: certain)

Rule: DAT2/DAT3 (genotype–metadata sync). Location:
`R/gl.filter.rdepth.r:128-132`.

`index <- (rdepth >= lower & rdepth <= upper)` contains NA wherever rdepth
is NA. The genlight subset `x[, index]` drops NA-index loci, but the
data.frame subset `x@other$loc.metrics[index, ]` returns an all-NA row for
each NA — so the output has fewer loci than metric rows and carries
garbage rows. Verified: with 10 NAs injected into testset.gl's rdepth,
nLoc = 127 but loc.metrics has 137 rows including 10 all-NA rows. Every
downstream locus-metric operation on such an object is corrupted, with no
warning.

Proposed change:
`index <- !is.na(rdepth) & rdepth >= lower & rdepth <= upper` — loci whose
read depth cannot be evaluated are removed cleanly (with the count visible
in the verbose >= 3 summary), and genotypes/metadata stay in sync.
Escalation-gate class: on datasets with NA read-depth metrics the output
changes from a corrupt object to a clean one with NA-depth loci removed.

### F2 — plot.file + plot.display=FALSE crashes (MEDIUM, confidence: certain)

Rule: spec axis; PLT3. Location: `R/gl.filter.rdepth.r:136-195`.

Same defect class as gl.report.rdepth (PR #255): `p1`/`p2`/`p3` are built
only inside `if (plot.display)`, but the save block references `p3`
unconditionally. Verified: "object 'p3' not found".

Proposed change: build the plots unconditionally; gate only `print(p3)`.

### F3 — Progress message states the wrong boundary logic (LOW, confidence: certain)

Rule: spec axis (message vs behaviour). Location:
`R/gl.filter.rdepth.r:118-126`.

The verbose >= 2 message says "Removing loci with rdepth <= lower and >=
upper" but boundary loci are retained (the docs and code agree:
removed strictly below lower / above upper). Verified: loci with
rdepth == 8 are retained at lower = 8 while the message claims they are
removed.

Proposed change: reword to "Removing loci with rdepth < lower or > upper".

### F4 — Missing `verbose == 0 → plot.display <- FALSE` guard (LOW, confidence: certain)

Rule: VRB2/PLT (ratified pattern). Location: `R/gl.filter.rdepth.r:72-74`.
Proposed change: add the standard guard after SET VERBOSITY.

### F5 — Visible return (LOW, confidence: certain)

Rule: FS/VRB5 (campaign precedent PRs #251/#252/#254). Location:
`R/gl.filter.rdepth.r:219`. Proposed change: `return(invisible(x2))`.

### F6 — @family tag leaks into the help title; header conformance (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC2, DOC7. Location: `R/gl.filter.rdepth.r:1-61`.

- Line 4 `#'  @family matched filter` is indented after the title
  continuation, so roxygen treats it as title text: the rendered
  `?gl.filter.rdepth` title ends with "… (read depth) @family matched
  filter" and the Rd has no \concept — verified in
  man/gl.filter.rdepth.Rd (0 concept entries).
- `@seealso` links to itself; should link gl.report.rdepth.
- Description calls the function its own "companion script" (text copied
  from the report's header).
- `@param upper` "[default infinite=1000]" — misleading; it is 1000.
- DOC2 verbose wording; DOC7 author missing "Author(s):"; tag order
  (`@details` after `@param`, `@return` after `@export`).

Proposed change: rewrite the header to the ratified template (docs-only;
fixes the visible title corruption).

### F7 — Structure: section order (LOW, confidence: certain)

Rule: FS4. Location: `R/gl.filter.rdepth.r:75-87`. SET WORKING DIRECTORY
and SET COLOURS precede FLAG SCRIPT START. Proposed change: move FLAG
SCRIPT START to directly follow SET VERBOSITY.

## Report notes (other functions / not fixed here)

- No locus-metric flag reset is needed here: removing loci leaves
  per-locus metrics valid for the retained loci (consistent with
  gl.filter.secondaries).
- Default plot.colors via gl.select.colors(...) returns the documented
  hexes (verified for gl.report.rdepth); left as-is.
- gl.report.rdepth's fixed version (PR #255) is the matched report; its
  quantile thresholds use the same inclusive `>=` convention as this
  filter's retention rule.

## Coverage

Characterization baseline: `tests/testthat/test-gl.filter.rdepth.R` — 15
assertions: retained counts vs independent recomputation on both
datatypes, inclusive boundaries, metadata sync on clean data, history +1,
input untouched, verbose-0 silence (display off), NA desync 127/137
(baseline), p3 crash (baseline), visible return (baseline), wrong message
wording (baseline). All 15 pass on the pre-fix code.

## Approval

Decisions recorded 2026-08-29 via approval boxes — all seven findings
**approved**, including both escalation-gate consequences (F1: NA-depth
loci removed cleanly instead of corrupting the object; F2: plot.file crash
becomes a working save). A follow-up decision during the Phase C caller
grep: the F4 guard is **kept** despite the
dartR.captive::gl.report.parent.offspring interaction (see Outcome).

## Outcome

Applied F1–F7. Verification: all 16 characterization assertions pass
post-fix; every diff from baseline maps to an approved finding (NA case:
127 loci = 127 metric rows, zero NA rows [F1]; RDS saved without display
[F2]; message "rdepth < 8 or > 50" [F3]; invisible return [F5]). The
verbose >= 3 summary now itemises NA-depth removals. End-to-end at
verbose = 3 on the NA-injected SNP set and on testset.gs; SilicoDArT
counts unchanged (220 = independent recomputation). The @family concept
now renders (man/gl.filter.rdepth.Rd has \concept{matched filter}); the
ten sibling matched-filter Rd files each gain the corresponding
cross-reference line (mechanical roxygen consequence of F6).

Caller grep across dartR.base + 7 siblings: two code callers in
dartR.captive, both assigning the return —
gl.filter.parent.offspring (plot.display = FALSE explicit; unaffected) and
gl.report.parent.offspring (verbose = 0 + plot.display = plot.filters).
For the latter, the kept F4 guard means plot.filters = TRUE no longer
displays the rdepth filter plots — decision: keep the convention guard;
**note for the future dartR.captive review: gl.report.parent.offspring
should pass verbose >= 1 or display its filter plots itself if
plot.filters = TRUE is to keep showing the rdepth histograms.**
dartr2shiny: not present in the workspace. NEWS entry added.

```json
{
  "function": "gl.filter.rdepth",
  "package": "dartR.base",
  "family_mode": "modify",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "HIGH", "rules": ["DAT2", "DAT3"], "loc": "R/gl.filter.rdepth.r:128-132", "status": "proposed"},
    {"id": "F2", "severity": "MEDIUM", "rules": ["spec", "PLT3"], "loc": "R/gl.filter.rdepth.r:136-195", "status": "proposed"},
    {"id": "F3", "severity": "LOW", "rules": ["spec"], "loc": "R/gl.filter.rdepth.r:118-126", "status": "proposed"},
    {"id": "F4", "severity": "LOW", "rules": ["VRB2", "PLT"], "loc": "R/gl.filter.rdepth.r:72-74", "status": "proposed"},
    {"id": "F5", "severity": "LOW", "rules": ["VRB5"], "loc": "R/gl.filter.rdepth.r:219", "status": "proposed"},
    {"id": "F6", "severity": "LOW", "rules": ["DOC1", "DOC2", "DOC7"], "loc": "R/gl.filter.rdepth.r:1-61", "status": "proposed"},
    {"id": "F7", "severity": "LOW", "rules": ["FS4"], "loc": "R/gl.filter.rdepth.r:75-87", "status": "proposed"}
  ],
  "datasets": ["testset.gl", "testset.gs"],
  "baseline_test": "tests/testthat/test-gl.filter.rdepth.R",
  "pr": 256
}
```
