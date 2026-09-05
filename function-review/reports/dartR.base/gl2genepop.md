# Review: gl2genepop (dartR.base)

- Family mode: io
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203
- Datasets: testset.gl, testset.gs
- Baseline: tests/testthat/test-gl2genepop.R (snapshot captured pre-review)

## Verdict

**Standards: Needs work** — SNP-only rejection is present (a bright spot in
this batch), but it stops with an empty condition message, the save-path
message is ungated and malformed, and the conversion loop is quadratic in
loci.
**Spec: Needs work** — the ordinary path writes valid genepop (title line,
comma-separated locus line, Pop blocks, 2/4- and 3/6-digit coding, 0000
missing — all verified), but an all-monomorphic object yields a silently
malformed file, and `pop.order` silently drops populations.

## Findings

**F1 [HIGH, confidence: high] — all-monomorphic input writes a malformed file (spec)**
`R/gl2genepop.r:148-152` — for a monomorphic locus the genind carries one
allele column, so `loc.fac` has no duplicated entries. When *no* locus has
two columns, `loci_names_l[-which(duplicated(loci_names_l))]` negative-
indexes with `integer(0)`, which empties the vector: `n.loci` becomes 0 and
`for (i in 1:n.loci)` iterates over `c(1, 0)` — two passes with `loc = NA`,
each appending a column of "0000".
Failure scenario: `gl2genepop` on a 5-locus all-monomorphic subset writes a
header claiming 5 loci and data rows containing exactly two "0000" fields
(confirmed; frozen in the baseline test). Downstream genepop parsers error
or misread. Arises naturally after subsetting individuals.
Proposed change: `loci_names <- unique(as.character(loci_names_l))`,
`n.loci <- length(loci_names)`, loop with `seq_len(n.loci)`. Mixed
monomorphic/polymorphic data is unaffected (single-column homozygotes
already code correctly).

**F2 [HIGH, confidence: high] — pop.order silently drops populations (spec)**
`R/gl2genepop.r:98-103` — `tmp[match(pop.order, names(tmp))]` keeps only the
named populations; names absent from `pop.order`, and misspelled names
(match returns NA -> NULL element), are dropped by `Reduce(rbind, ...)`
without any message (both confirmed).
Failure scenario: `pop.order = c("C", "A")` on a three-population object
exports two populations; a typo (`"Bee"` for `"B"`) exports without pop B.
The user's genepop file is missing data with no diagnostic.
Proposed change: validate `setequal(pop.order, popNames(x))` and
`stop(error(...))` on mismatch, naming the offending entries.

**F3 [MEDIUM, confidence: high] — ungated, malformed save-path message (VRB5)**
`R/gl2genepop.r:208-211` — `cat(report("The genepop file is saved as:",
file.path(outpath, outfile, "\n")))` prints at every verbosity (confirmed at
`verbose = 0`), and passes `"\n"` to `file.path()` as a path component, so
the printed path gains a spurious trailing separator.
Failure scenario: silent pipelines emit one line per call showing a path
that, copy-pasted, ends in `/`.
Proposed change: gate at `verbose >= 2` and move `"\n"` outside
`file.path()`.

**F4 [MEDIUM, confidence: high] — quadratic conversion loop (STY2)**
`R/gl2genepop.r:152-183` — `col_loc <- which(loc_all[, "locus"] == loc)` is
recomputed inside the per-individual loop although it depends only on the
locus, giving O(nInd x nLoc^2) string comparisons overall.
Failure scenario: testset.gl-sized objects take noticeable seconds; a 50k-
locus dataset takes hours for a linear-time task.
Proposed change: hoist `col_loc` above the `j` loop (one-line fix); a
vectorised rewrite per locus is the fuller cleanup.

**F5 [LOW, confidence: high] — SilicoDArT rejection raises an empty error (VRB2)**
`R/gl2genepop.r:70-75` — the check prints via `cat(error(...))` then calls
bare `stop()`, so the raised condition has an empty message.
Failure scenario: `tryCatch` handlers and logs see `Error:` with no text;
with output sinks the printed explanation can be separated from the error.
Proposed change: `stop(error("Only SNP (diploid) data can be transformed
into genepop format!\n"))`.

**F6 [LOW, confidence: high] — single-individual object fails with an opaque error (spec)**
`R/gl2genepop.r:106-107` — a one-individual genlight dies inside the
genind conversion path with "X is not a matrix" (confirmed).
Failure scenario: exporting a single sample gives an error naming neither
the function nor the cause.
Proposed change: either support it (drop = FALSE discipline through the
matrix operations) or fail fast with a clear message.

**F7 [LOW, confidence: high] — all-NA loci silently dropped; verbose text off-canon (VRB4, DOC2) (proposed rule)**
`R/gl2genepop.r:92` — `gl.filter.allna(x, verbose = 0)` silently reduces the
locus set at any verbosity; the `verbose` doc default clause also deviates
from the DOC2 canon.
Failure scenario: file locus count differs from the object with no message.
Proposed change: report non-zero drops at `verbose >= 1`; canonical DOC2
clause.

## Proposed changes

1. Replace the duplicated()-negative-index idiom with `unique()` plus
   `seq_len()` (F1).
2. Validate `pop.order` against `popNames(x)`; error on mismatch (F2).
   **Consequence: calls that today silently export a subset will error.**
3. Gate and fix the save-path message (F3).
4. Hoist `col_loc` out of the individual loop (F4).
5. `stop(error(...))` for the SilicoDArT rejection (F5).
6. Clear failure (or support) for single-individual objects (F6).
7. Report all-NA locus drops; DOC2 canonical verbose text (F7).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP (stringr/tidyr in Imports), PLT
  (n/a), STY — run
- Spec: written file vs genepop format on testset.gl — run (title/locus
  lines, Pop blocks, alphabetic ordering, 2-digit and 3-digit coding,
  missing codes, ID-comma layout verified in the baseline test)
- Edge cases: all-monomorphic (F1), pop.order subset/typo (F2), single
  individual (F6), SilicoDArT rejection (F5) — run
- Heterozygote/homozygote coding cross-checked against dosage matrix — run
- Locus names containing "." (would trip the locus.allele split guard):
  not exercised — testset locus names contain no dots; the guard fails
  loudly, so no silent risk
- FBM path (DAT6): SKIPPED — no FBM fixture exercised for this function
- Genome-position lens: N/A — the format carries no coordinate fields

## Approval (Phase B)

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges (2026-09-05) | Blanket approval, incl. output-change consequence |
| 2 | Approved | Arthur Georges (2026-09-05) | Silently subsetting calls now error |
| 3 | Approved | Arthur Georges (2026-09-05) | Blanket approval, all severities |
| 4 | Approved | Arthur Georges (2026-09-05) | Blanket approval, all severities |
| 5 | Approved | Arthur Georges (2026-09-05) | Blanket approval, all severities |
| 6 | Approved | Arthur Georges (2026-09-05) | Fail-fast option adopted |
| 7 | Approved | Arthur Georges (2026-09-05) | Blanket approval, all severities |

All findings approved 2026-09-05 via the formal approval boxes: blanket
class approval at every severity, explicitly including the consequence
that converted outputs change where they were wrong. (The function
already rejected SilicoDArT; F5 upgrades the rejection to a proper
error condition.)

## Outcome (Phase C)

All seven changes applied on branch `review-gl2genepop` (base
`upstream/dev`, ddaed27): locus list built with `unique()` +
`seq_len()` so all-monomorphic objects write one correctly coded
genotype per claimed locus (F1); `pop.order` validated against
`popNames(x)` with a fatal error naming unlisted, misspelled and
duplicated entries (F2); save-path message gated at `verbose >= 2` with
the stray `"\n"` moved out of `file.path()` (F3); `col_loc` hoisted out
of the per-individual loop, removing the O(nInd x nLoc^2) recomputation
(F4); the SilicoDArT rejection raises `stop(error(...))` with its
message in the condition (F5); single-individual objects fail fast with
a clear message (F6); non-zero all-NA locus drops reported at
`verbose >= 1` and the DOC2 canonical verbose clause adopted (F7).
Baseline characterization test updated for five approved diffs (F1, F2,
F3, F5, F6); 21/21 assertions pass — the ordinary-path snapshots (title
and locus lines, Pop blocks, alphabetic ordering, 2-digit and 3-digit
coding, exact data rows) are byte-identical. Old-vs-new check on the
all-monomorphic fixture: upstream code reproduces the frozen malformed
`"A_AA010915, 0000 0000"` row; the fix writes 5 correctly coded fields.
End-to-end run on `testset.gl` at `verbose = 3` wrote 307 lines (2
header + 31 Pop lines + 274 individuals), the F7 warning reporting this
checkout's 3 all-NA loci. Sibling-caller grep across the 8 dartRverse
clones: one caller, `dartR.popgen/R/gl.LDNe.r:190`, which calls with
default `pop.order` on SNP data — unaffected. NEWS entry added.

```json
{
  "function": "gl2genepop",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "spec", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "spec", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "STY2", "status": "applied", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "applied", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "spec", "status": "applied", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "VRB4", "status": "applied", "change": 7}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "dotted locus names not exercised", "genome-position lens: format has no coordinate fields"],
  "status": "phase-c",
  "pr": null
}
```
