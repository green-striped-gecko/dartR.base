# Review: gl.smearplot (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.smearplot
- Datasets: testset.gl, testset.gs
- Family mode: report (graphics; read-only — input untouched verified,
  no history append, correct)
- Checks skipped: visual appearance not assessed beyond scale
  colors/labels (headless); Google Group not searched (not available:
  no browser session).

## Verdicts

- **Standards: FAIL** — the documented `plot.display` parameter is
  accepted and guarded (`verbose==0 → plot.display <- FALSE`) but never
  used: `print(p3)` is unconditional, so the plot renders at verbose 0
  and with plot.display = FALSE alike (verified: full page drawn on a
  pdf device); both package checks use `cat(error()) + return(-1)`;
  the het.only-on-SilicoDArT warning prints at verbose 0; plotly is
  used for `interactive = TRUE` but never checked (it is in Suggests);
  header drift.
- **Spec: FAIL** — the SilicoDArT legend has never shown its intended
  labels: `labels_silicodart["0"] <- "Absence"` assigns by NAME into an
  unnamed `c("0","1")` vector, appending instead of replacing, so
  ggplot receives `c("0","1","Absence","Presence")` and displays "0"
  and "1" (verified). And `het.only = TRUE` with SilicoDArT data warns
  "Set to FALSE" but renders BOTH presence and absence as #d3d3d3 —
  the palette override at the top of the function happens before the
  datatype branch, and the "restore" snapshot was taken after the
  override (verified: `palette(2)` returns gray, gray).

## Findings

### F1 — plot.display never used; print(p3) unconditional (MEDIUM, confidence: certain)

Rule: PLT/VRB (the p3-gate class). Location: PRINTING OUTPUTS block.

`plot.display = FALSE` (or `verbose = 0`, which sets it) still renders
the plot — verified by page size on a pdf device. Proposed change:
`if (plot.display) print(p3)`. Saving via plot.file remains
unconditional (save-without-display stays possible).

### F2 — SilicoDArT legend labels never applied (MEDIUM, confidence: certain)

Rule: spec axis (verified against the description's promise of
Absence/Presence coding). Location: labels_silicodart block.

Named assignment into an unnamed vector appends: the labels vector
reaching ggplot is length 4 and the legend shows "0"/"1". Proposed
change: `labels_silicodart <- c("Absence", "Presence")` (direct), with
the dead `which(is.na())` line removed.

### F3 — het.only + SilicoDArT renders an all-gray plot (MEDIUM, confidence: certain)

Rule: spec axis / DAT dispatch. Location: het.only override (top) and
SilicoDArT branch.

The het.only palette override is applied before the datatype dispatch,
and `plot.colors.hold` snapshots the already-overridden palette, so the
"restored" colors in the SilicoDArT branch are the het-only grays:
presence and absence both #d3d3d3 despite the warning claiming the
option was disabled. The warning also prints at verbose 0. Proposed
change: apply the override only when `datatype == "SNP"` (datatype is
already known at that point); keep the SilicoDArT warning, gated at
verbose >= 2.

### F4 — return(-1) package checks; plotly unchecked (MEDIUM, confidence: certain)

Rule: DEP/STY (the return(-1) class). Location: package-check block and
interactive branch.

Both checks print via cat and `return(-1)` — callers get -1 instead of
an error condition. plotly (Suggests) is used when `interactive = TRUE`
with no check at all — a user without plotly gets a raw namespace
error. Proposed change: `stop(error(...))` for both existing checks;
add the same check for plotly inside the interactive branch. (reshape2
is in Imports, so its check is vestigial but harmless — retained.)

### F5 — Missing-data legend entry never appears; dead lines (LOW, confidence: certain)

Rule: spec/DOC. Location: both label blocks.

`labels_genotype[which(is.na(labels_genotype))] <- "Missing data"` is
dead in both blocks — NAs were removed from the vector beforehand — so
the legend never carries a Missing data entry even though the
description promises NA is color coded. NA cells ARE colored (via
na.value); only the legend entry is absent. Proposed change: remove the
dead lines; docs note that missing data are shown in the NA color
without a legend entry. (Adding a true NA legend entry would change the
rendered legend for every user — not proposed.)

### F6 — Header and style conformance (LOW, confidence: certain)

Rules: DOC1, DOC2, STY. `@param group.pop` documented "[default TRUE]"
but the signature default is FALSE; a stray duplicated "[default
FALSE]." line under `den`; the plot.dir doc line glued to a trailing
comment marker; description grammar ("cluttered if ind.labels If there
are too many"); curly quotes in the legend param; `@return` after
`@export`; `seq(1:nInd(x))` idiom; dead `individuals <-
df.matrix3$Label` assignment in the den branch; den silently forces
`group.pop <- FALSE` (message added at verbose >= 2); the loc.order
chromosome guard hardened to `length()` (a zero-length chromosome slot
would pass the `is.null` check and reorder by an empty vector, dropping
every locus).

## Report notes (other functions / not fixed here)

- `gl.dist.ind` is called in the den path with `dist()` applied to its
  result — distance-of-distances; the dendrogram is therefore built on
  a transformed metric. Behaviour is stable and the path works;
  flagged for consideration when gl.dist.ind or the den feature is
  reviewed on its own.
- Caller check for F1: `gl.random.snp` and `gl.randomize.snps` call
  gl.smearplot at verbose = 0 solely to build p1/p2 for their own
  composed plot — pre-fix, the unconditional print drew spurious
  intermediate pages inside them; F1 improves both, no breakage.
  Separately, both compose `p3 <- p1 / p2` OUTSIDE their plot.display
  gate, so with plot.display = FALSE they crash on 'p1' not found —
  their own defect (the p3-not-found class), for their own reviews.

## Coverage

Characterization baseline: `tests/testthat/test-gl.smearplot.R` — 15
assertions: SNP plot class/visibility/silence/input-untouched/labels,
page-rendered-despite-plot.display-FALSE (baseline), SilicoDArT
4-element label vector (baseline), het.only+SilicoDArT gray palette and
verbose-0 warning (baseline), SNP het.only palette, no-hom-alt subset
labels, loc.order warning + locus count, den path. All 15 pass on the
pre-fix code. Visual appearance beyond scales/labels not asserted
(headless).

## Approval

All six findings approved via the approval boxes (2026-08-30).

## Outcome

All six findings applied. Characterization suite: 15/15 pass; every
flipped assertion maps to an approved finding (no page rendered at
plot.display = FALSE [F1]; SilicoDArT legend labels now the 2-element
Absence/Presence vector [F2]; het.only + SilicoDArT palette #0000FF /
#FF0000 with the warning gated [F3]). One test-harness artifact fixed
during verification: capture.output() auto-prints a visibly returned
ggplot, so the no-page assertion needed the result assigned inside the
capture — a test correction, not a code diff. End-to-end at verbose 3
clean on SNP, SilicoDArT + het.only, and den + group.pop paths (the
new gated messages display as intended). Caller grep: no functional
callers of gl.smearplot in the live family packages. NEWS entry added.

```json
{
  "function": "gl.smearplot",
  "package": "dartR.base",
  "family_mode": "report",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "rules": ["PLT", "VRB"], "loc": "R/gl.smearplot.r print", "status": "applied"},
    {"id": "F2", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl.smearplot.r silico labels", "status": "applied"},
    {"id": "F3", "severity": "MEDIUM", "rules": ["spec", "DAT"], "loc": "R/gl.smearplot.r het.only", "status": "applied"},
    {"id": "F4", "severity": "MEDIUM", "rules": ["DEP", "STY"], "loc": "R/gl.smearplot.r pkg checks", "status": "applied"},
    {"id": "F5", "severity": "LOW", "rules": ["spec", "DOC"], "loc": "R/gl.smearplot.r label blocks", "status": "applied"},
    {"id": "F6", "severity": "LOW", "rules": ["DOC1", "DOC2", "STY"], "loc": "R/gl.smearplot.r header/body", "status": "applied"}
  ],
  "datasets": ["testset.gl", "testset.gs"],
  "baseline_test": "tests/testthat/test-gl.smearplot.R",
  "pr": 277
}
```
