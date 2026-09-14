# Review: gl.filter.ld (dartR.base)
- Family mode: modify
- Date: 2026-09-10
- Reviewer: Claude (Claude Fable 5), dartr-function-review v1.0.0
- Package commit: f5e7b72 (upstream/dev; working copy verified identical by `git diff upstream/dev -- R/gl.filter.ld.r`)
- Datasets: platypus.gl (mapped, 383 loci after callrate/monomorph filtering), testset.gl (30-population subset) — dartR.data 1.2.5
- Baseline: tests/testthat/test-gl.filter.ld.R (new file, snapshot captured pre-review; 17 assertions, all passing)

## Verdicts

**Standards: Needs work** — the FS backbone is present, genotype–metadata
sync holds on every path tested, and the filter drops loci from the original
(un-subset) object as it should; but one message prints at `verbose = 0`,
each call writes two history entries, and objects without dartR flags die
with an opaque length-zero error.

**Spec: Needs work** — the documented tie-break rule ("if the SNP is already
in the list, the other SNP will be kept") is not what the code does: a locus
is still removed when its only LD partner was already removed, so chained LD
clusters are over-filtered. The `pop.limit` default works only by an
accident of lazy evaluation.

## Independent verification (spec axis)

On the platypus LD report (496 pairs), reconstructing the *documented*
algorithm and the *implemented* algorithm gives identical drop sets (34 loci
at `pop.limit = 1`) — real data did not happen to exercise the divergence.
A crafted three-locus chain does: pairs (L1,L2) and (L2,L3) with
`stat.keep` L1 > L2 > L3. The code drops L2 (correct) and then also L3,
although L3's only LD partner L2 is already gone; the documentation says L3
is kept. Confirmed empirically (test
"gl.filter.ld drops the partner of an already-dropped locus").

## Findings

**F1 [HIGH, confidence: high] — a locus is removed although its only LD partner was already removed (DOC5 (proposed rule))**
`R/gl.filter.ld.r:104-116` — the membership test `loci_tmp %in%
loci_list[[i]]` covers only the locus about to be dropped, not both members
of the pair, so a pair whose better member is already in the drop list still
nominates its worse member.
Failure scenario: loci A–B–C where A>B>C on `stat.keep`, pairs (A,B) and
(B,C) both above threshold, C in LD with nothing else: B is dropped for A,
then C is dropped for B — but B is no longer in the dataset, so C is
pseudo-replicated with nothing and is lost for no reason. Confirmed on a
crafted report; on strongly chained LD clusters the filter removes up to all
but one locus of material that the documented rule would keep.
Proposed change: skip the pair when either member is already listed
(`if (a %in% list || b %in% list) next`), which is exactly the documented
rule.
**Consequence: numerical output changes — fewer loci are removed wherever
LD clusters chain (the platypus baseline happens to be unaffected).**

**F2 [MEDIUM, confidence: high] — "No pair of loci" message prints at verbose = 0 (VRB5, VRB3)**
`R/gl.filter.ld.r:84-86` — the no-pairs branch `cat(report(...))` is not
gated by `verbose`.
Failure scenario: a pipeline run at `verbose = 0` still prints three lines
whenever the threshold removes nothing.
Proposed change: gate at `verbose >= 1` (the user should still learn the
filter was a no-op — VRB4).

**F3 [MEDIUM, confidence: high] — two history entries per call, one leaking the internal implementation (FS8)**
`R/gl.filter.ld.r:123,136-137` — `gl.drop.loc` appends its own
`match.call()` (`gl.drop.loc(x_hold, loc.list = loci_names, verbose = 0)`,
exposing internal variable names) and `gl.filter.ld` then appends a second
entry.
Failure scenario: `length(x@other$history)` grows by 2 per call (confirmed:
6 → 8); replaying a history that contains `x_hold`/`loci_names` fails.
Proposed change: reset the history to the input object's before appending
the single `gl.filter.ld` entry.

**F4 [MEDIUM, confidence: high] — missing loc.metrics.flags kills the call with an opaque error (DAT5)**
`R/gl.filter.ld.r:69` — `if (x@other$loc.metrics.flags$monomorphs == FALSE)`
on an object without the flag evaluates `if (logical(0))`.
Failure scenario: a genlight not built by dartR (flags absent) errors with
"argument is of length zero" before any useful message. Confirmed.
Proposed change: `if (!isTRUE(x@other$loc.metrics.flags$monomorphs))` so
both FALSE and absent produce the existing warning instead of a crash.

**F5 [MEDIUM, confidence: high] — ld.report is never validated (FS5)**
`R/gl.filter.ld.r:80-95` — the function indexes
`ld.report$pop`, `$ld.stat`, `$locus_a.stat.keep`, ... directly; any other
data frame (or a genlight passed by mistake) produces obscure downstream
errors.
Failure scenario: `gl.filter.ld(x, ld.report = x)` or a hand-rolled data
frame missing one column fails deep inside `gl.keep.pop` or the split loop
with no hint that `ld.report` must be the output of `gl.report.ld.map`.
Proposed change: fail fast when `ld.report` lacks any required column, with
a message naming `gl.report.ld.map`.

**F6 [LOW, confidence: high] — pop.limit default works by lazy-evaluation accident (STY1, DOC5 (proposed rule))**
`R/gl.filter.ld.r:49,80` — the default `ceiling(nPop(x)/2)` is evaluated at
first use, which is *after* `x <- gl.keep.pop(x, pops-in-report)`, so it
actually means "half the populations represented in the report", not "half
of the populations" of the input object. On testset.gl (30 populations, 3
in the report) the default resolves to 2, not 15 — which is the sensible
value, but only because of evaluation order that the next edit could break.
Failure scenario: none today; reordering the first lines of the body (or
using `pop.limit` before the reassignment) silently changes the default to
`ceiling(30/2) = 15` and turns the filter into a permanent no-op for
datasets with many small populations.
Proposed change: compute the default explicitly from the report's
populations at the top of the body and document it as such.

**F7 [LOW, confidence: high] — documentation gaps (DOC1, DOC3)**
`R/gl.filter.ld.r:15-16` — `threshold` says "above which" but the code uses
`>=`; there is no `@details`; no `@return` explanation of the invisible
return; the two `@examples` blocks overlap (the second, unconditional one
reruns the full report+filter and is the heavier of the two — CRAN runs it).
Failure scenario: CRAN incoming check time; user confusion at the boundary
value.
Proposed change: "at or above which"; add a short `@details`; wrap the
second example in the same `requireNamespace` guard but `\donttest{}` (or
merge the two blocks).

## Proposed changes

1. Implement the documented pair-resolution rule: skip a pair when either
   member is already listed (F1).
   **Consequence: numerical output changes — fewer loci removed for chained
   LD clusters; unchanged on the platypus baseline.**
2. Gate the no-pairs message at `verbose >= 1` (F2).
3. Append a single history entry on the returned object (F3).
4. Null-safe monomorphs-flag check (F4).
5. Validate `ld.report` columns at entry with an informative error (F5).
6. Evaluate the `pop.limit` default explicitly and document its actual
   meaning (F6).
7. Documentation fixes and example consolidation (F7).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT (n/a — no plot), STY — run
- Spec: documented vs implemented drop rule — run (real data + crafted
  chain); genotype–metadata sync after filtering — run; history — run
- GitHub issues: #210 ("gl.filter.ld breaks on subsets") checked — NOT
  reproducible on dev f5e7b72 (the example's 50-locus subset runs clean:
  report 1056x11, filter 50 → 46 loci); candidate for closing
- FBM path (DAT6): SKIPPED — no FBM fixture for the report+filter chain
- Google Group search: SKIPPED — not queried this session

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | rejected | Arthur | keep the implemented sequential rule; the documentation is corrected to describe it (folded into change 7) |
| 2 | approved | Arthur | |
| 3 | approved | Arthur | |
| 4 | approved | Arthur | |
| 5 | approved | Arthur | |
| 6 | approved | Arthur | |
| 7 | approved | Arthur | includes replacing the false pair-resolution sentence with the actual rule |

Cross-package caller grep (API3): no calls to `gl.filter.ld` in the local
dartR.* clones outside dartR.base. All clear.

## Outcome

Changes 2-7 applied on branch review-gl.filter.ld (commit a5df83d), PR
green-striped-gecko/dartR.base#391.

- Characterization suite green (20 assertions); every diff from the
  pre-review baseline maps to an approved change: silence at verbose 0
  (F2), single history entry (F3), no-crash on missing flags (F4), the new
  fail-fast validation error (F5).
- Loci removed unchanged: platypus baseline 383 -> 380 (defaults) and 349
  (pop.limit = 1); default-pop.limit equivalence with the explicit
  computation asserted on the 30-population testset fixture.
- End-to-end run at verbose = 3 on platypus.gl clean.
- The pop.limit signature default changed from `ceiling(nPop(x) / 2)` to
  `NULL` (resolved to the same value); recorded in NEWS and the PR's API
  impact section.

```json
{
  "function": "gl.filter.ld",
  "package": "dartR.base",
  "family": "modify",
  "skill_version": "1.0.0",
  "commit": "f5e7b72",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "rejected", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "FS8", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DAT5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "STY1", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 7}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "Google Group: not queried"],
  "status": "pr-open",
  "pr": 391
}
```
