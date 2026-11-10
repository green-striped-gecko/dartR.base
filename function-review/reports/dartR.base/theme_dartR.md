# Review: theme_dartR (dartR.base)

- Family mode: analysis (graphics utility; PLT checks applied)
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 2bd61c5 (`R/theme_dartR.r` identical to origin/dev)
- Environment: ggplot2 4.0.2
- Datasets: none needed (theme object; test plot on a small data frame)
- Baseline: `tests/testthat/test-theme_dartR.R` (12 expectations, captured pre-review, all pass)
- Author: not stated in the file; history shows Bernd Gruber (2023-06) and Arthur Georges (2023-07)

## Verdict

**Standards: Ready** — a single `ggplot2::theme()` call; trivial fixes only (return value, header).
**Spec: Ready** — builds plots without warnings under ggplot2 4.0.2; one element ignores `base_size`.

## Findings

**F1 [LOW, confidence: high] — facet strip text ignores `base_size` (PLT1)**
`R/theme_dartR.r:153` — `strip.text.x = element_text(size = 14)` is an absolute size, while every other text element is relative to `base_size`. With `base_size = 8`, facet labels are 14 pt against a 9.6 pt title; with `base_size = 20`, they are smaller than the axis text.
Failure scenario: a user shrinks or enlarges a figure through `theme_dartR(base_size = )`; facet labels stay at 14 pt.
Proposed change: `size = rel(14 / 11)`, which gives the same 14 pt at the default `base_size = 11` and scales with other values.

**F2 [LOW, confidence: high] — theme returned invisibly (FS10)**
`R/theme_dartR.r:41` — the last expression is the assignment `t <- theme(...)`, so the theme is returned invisibly; `theme_dartR()` typed at the console prints nothing. `+ theme_dartR()` is unaffected.
Failure scenario: a user inspects the theme at the console and sees no output.
Proposed change: return the `theme()` call directly.

**F3 [LOW, confidence: high] — documentation (DOC1, DOC7)**
`@param` entries lack defaults; `@return` reads "a the standard dartR theme"; no `@author` with custodian; the header comment refers to `half_size`, which the code calls `half_line`; stale "Version v.2023.2" comment.
Proposed change: rewrite the header. Docs only.

## Proposed changes

1. Make facet strip text relative to `base_size` (`rel(14 / 11)`); unchanged at the default size (F1).
2. Return the theme visibly (F2).
3. Rewrite the header; add author/custodian (F3). Docs only. Custodian to be confirmed.

## Coverage

- Standards walk: FS, DOC, PLT, STY — run. VRB, DAT, DEP: not applicable (no genlight input, no verbosity, ggplot2 in Depends).
- Plot build with facets, title, caption and tag under ggplot2 4.0.2 — run, no warnings.
- Complete-theme gaps: the theme is marked `complete = TRUE` but leaves 78 of `theme_grey()`'s elements unset (for example `axis.minor.ticks.*`, `geom`). ggplot2 fills missing elements of a complete theme from its own default theme, so plots build unchanged; not a finding.
- Deprecated elements `legend.text.align` / `legend.title.align` are set to `NULL` and raise no warning under 4.0.2.
- Callers: 37 in dartR.base, 34 files in sibling packages; changes 1-2 do not alter plots at the default `base_size`.
- Known complaints: none found.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | |
| 2 | approved | Luis | |
| 3 | approved | Luis | |

## Outcome

- Changes 1-3 applied on branch `review-theme_dartR` (from `origin/dev`).
- Characterization test: 12 expectations pass. Diffs from baseline map to approved changes only: strip text resolves to 14 pt at `base_size = 11` and 28 pt at 22 (1); `theme_dartR()` is visible (2).
- Old vs new theme at the default size: `strip.text.x` is the only element that differs, and it resolves to the same 14 pt in a built plot.
- No caller in dartR.base or the sibling packages passes a non-default `base_size`.
- `devtools::document()` run. Custodian set to Bernd Gruber (earliest author in history) — assumption, to be confirmed.
- Full `R CMD check` not run.
- PR: pending.

## Machine block

```json
{
  "function": "theme_dartR",
  "package": "dartR.base",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "2bd61c5",
  "verdict_standards": "ready",
  "verdict_spec": "ready",
  "findings": [
    {"id": "F1", "severity": "LOW", "confidence": "high", "rule": "PLT1", "status": "approved", "change": 1},
    {"id": "F2", "severity": "LOW", "confidence": "high", "rule": "FS10", "status": "approved", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "approved", "change": 3}
  ],
  "coverage_skipped": [],
  "status": "pr-open",
  "pr": null
}
```
