# Review: utils.read.dart (dartR.base)

- Family mode: io (reviewed as part of the DArT read chain: `gl.read.dart` -> `utils.read.dart` -> `utils.dart2genlight`; sibling `gl.read.silicodart`)
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ddaed27 (upstream/dev). Empirical runs on checkout ed99203 (integration-local); `R/utils.read.dart.r` is identical between the two.
- Datasets: `testset_SNPs_2Row.csv` + `testset_metadata.csv` (dartR.data); synthetic fixtures (layouts below)
- Baseline: `tests/testthat/test-utils.read.dart.R` (new; snapshot captured pre-review, 20 expectations passing)

## Synthetic fixtures

No 1-row SNP fixture exists in dartR.data, so parsing branches were exercised with
minimal csvs built in the scratchpad (each: 5 `*` header rows carrying service,
barcode, plate, row, column; header `AlleleID,SNP,SnpPosition,CallRate,AvgCountRef,AvgCountSnp,RepAvg,ind1..3`):

- `synth_1row.csv` — 3 loci, codes 0/1/2 (valid 1-row)
- `synth_1row_nohet.csv` — 4 loci, codes 0/1 only (no het calls)
- `synth_2row_triplet.csv` — 3 loci in 2-row form, one uid present 3 times
- `synth_2row_quad.csv` — 3 loci, one uid present 4 times (true duplicate)
- `synth_2row_dupind.csv` — 3 loci, two individual columns both named `indA`

The baseline test file reconstructs these inline via `tempfile()`, so the tests
are self-contained.

## Verdicts

**Standards: Needs work** — preamble (FS2/FS3), colour helpers, and topskip
autodetection are sound; verbosity gating fails at `verbose = 0` and one
warning block is dead code.
**Spec: Needs work** — well-formed 1-row and 2-row files parse correctly
(verified against hand-read cells), but two malformed-input paths destroy data:
an odd-count clone ID wipes every locus, and a het-free 1-row file is silently
read as 2-row.

## Findings

**F1 [HIGH, confidence: high] — duplicate-ID cleanup removes every locus when any uid has an odd count (DAT2)**
`R/utils.read.dart.r:217-237` — when `table(table(uid))` has more than one
entry and the second entry is not "4", the else branch removes every uid with
count > 1. In 2-row data every valid locus has count 2, so all loci are
removed.
Failure scenario (verified, `synth_2row_triplet.csv`): one uid occurring 3
times (a hand-edited report with one allele row deleted or duplicated) strips
all 7 rows down to 0/1 orphans; the read then dies with the misleading error
"The DArT format must be either 1row or 2row". A dataset-scale loss triggered
by a single malformed row.
Proposed change: in the non-quad branch, remove only uids whose count differs
from the expected `nrows` count (or the modal count), and report which uids
were removed; never remove uids whose count equals the modal count.

**F2 [HIGH, confidence: high] — 1-row files without heterozygote calls are silently misread as 2-row (DAT1)**
`R/utils.read.dart.r:246` — format detection is `gnrows <- 3 - max(datas)`.
A genuine 1-row file whose called genotypes happen to be only 0/1 (small
subsets, inbred or haploidised panels, heavy filtering) has max 1 and is
detected as 2-row.
Failure scenario (verified, `synth_1row_nohet.csv`): 4 loci x 3 individuals
came back as 2 "loci" with genotypes fabricated by pairing unrelated adjacent
rows. Only an ungated console warning ("no. rows per Clone does not fit")
marks it; the function proceeds and returns the corrupt list.
Proposed change: cross-check `max(datas)` against the uid count structure
(`tt`) — the two independent signals already computed — and stop with a clear
message when they disagree, instead of warning and continuing.

**F3 [MEDIUM, confidence: high] — service/plate extraction reads past the header block (no rule fits; input-bounds principle)**
`R/utils.read.dart.r:112-124` — `service.row` (1), `plate.row` (3), and rows
`plate.row + 1/2` (4, 5) are indexed into `tdummy` with no check against
`topskip`. When the file has fewer header rows than `plate.row + 2`, the
column-header row and first data row are silently pasted into
`plate_location`.
Failure scenario (verified on the canonical `testset_SNPs_2Row.csv`, which has
3 header rows): `service` = "A5" (actually the plate-well row) and
`plate_location` = "UC_1-AA0109150" (plate name + individual name + first
genotype) flow into `ind.metrics` of every object built from this fixture.
Proposed change: validate `service.row` and `plate.row + 2` against `topskip`;
set service/plate to NA with a gated warning when out of range.

**F4 [MEDIUM, confidence: high] — individual-name uniquification is dead code and its warning misstates the result (DOC5)**
`R/utils.read.dart.r:126-137` — `ind.names` is trimmed and passed through
`make.unique(sep = "_")`, but the result is never applied: `gendata` keeps the
columns as `read.csv` delivered them, and `utils.dart2genlight` uses
`colnames(sdata)`. The printed warning claims "_n" suffixes are added.
Failure scenario (verified, `synth_2row_dupind.csv`): duplicate column `indA`
arrives as `indA.1` (read.csv's own suffix), not the promised `indA_1`; a user
whose ind.metafile lists `indA_1` per the warning will fail to match, and the
whole block is misleading to maintainers.
Proposed change: either assign the uniquified names to `colnames(datas)` (and
keep the `_n` claim true), or delete the block and document read.csv's
behaviour.

**F5 [MEDIUM, confidence: high] — verbose = 0 is not silent (VRB5, VRB3)**
`R/utils.read.dart.r:130, 212, 231, 266` — four `cat(warn(...))` calls are
ungated.
Failure scenario (verified): `verbose = 0` still prints the duplicate-name and
row-count warnings; scripted pipelines cannot silence the function.
Proposed change: gate result-affecting warnings at `verbose >= 1` (VRB4,
proposed rule) and the rest at `verbose >= 2`; the duplicate-removal and
format-mismatch warnings affect results, so `>= 1`.

**F6 [LOW, confidence: high] — user-supplied topskip silently disables service/plate extraction (DOC5)**
`R/utils.read.dart.r:110-124` — service/plate come from `tdummy`, which only
exists when topskip was auto-detected. Passing `topskip` explicitly yields all
NA service/plate with no message (verified: `topskip = 3` on the testset).
Proposed change: read the header rows for service/plate regardless of whether
topskip was supplied, or document the limitation in `@param topskip`.

**F7 [LOW, confidence: high] — @return documents 5 list elements; the function returns 7 (DOC5)**
`R/utils.read.dart.r:29-31` — the return list gained `service` and
`plate_location` without a doc update.
Proposed change: rewrite `@return` to name all 7 elements.

**F8 [INFO, confidence: medium] — minor robustness and wording**
`R/utils.read.dart.r:66` `sum(tdummy[, 1] == "*")` has no `na.rm` (an NA in
column 1 of the header block makes the `if` fail uninformatively);
`R/utils.read.dart.r:213` the duplicate-removal message reports `tt[2]` (count
of affected uids in one count class) as a SNP count; `R/utils.read.dart.r:273-283`
prints the pre-removal `tt`, so after a quad removal the log still shows "4".
No proposed change beyond folding into F1's rewrite.

## Proposed changes

1. Rewrite the duplicate-uid removal to target only uids whose count deviates
   from the modal/expected count, with an accurate report of removed uids (F1,
   F8). **Consequence: files previously wiped now load with only the bad uid
   removed — numerical output changes for malformed inputs (API1).**
2. Make format detection stop (clear error) when genotype-range and uid-count
   signals disagree, rather than warn-and-continue (F2). **Consequence:
   het-free 1-row files that previously "loaded" (corrupt) now error.**
3. Bounds-check `service.row`/`plate.row` against `topskip`; NA + gated
   warning when out of range (F3).
4. Apply (or remove) the individual-name uniquification; make the warning
   text match reality (F4).
5. Gate all warnings per VRB3/VRB4 (F5).
6. Extract service/plate independently of topskip auto-detection, or document
   the gap (F6).
7. Correct `@return` to the actual 7-element list (F7).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, STY — run
- Spec: behaviour vs roxygen on real 2-row fixture + 5 synthetic variants — run
- Genotype encoding hand-check: done downstream in `utils.dart2genlight`
  review (this function returns raw rows)
- 2-row silicodart variant: SKIPPED — handled by `gl.read.silicodart`, which
  hardcodes 1-row; no 2-row silico format is known to exist
- Recent-format header variants (`MarkerName`/`Variant` columns): partially
  run — `MarkerName` branch (line 163) not exercised; no fixture with that
  header was available and the branch mirrors the `AlleleID` branch
- FBM path (DAT6): SKIPPED — fbm conversion happens in `gl.read.dart`, out of
  this function's scope
- PLT: not applicable (no plots)

## Approval (Phase B)

Approved 2026-09-05 by Arthur Georges via the formal approval boxes, covering
all BLOCKER/HIGH/MEDIUM findings for the read-chain review round, with
explicit acknowledgment of the consequence that objects users previously read
will differ where current behaviour is wrong (including the API1
consequences called out on changes 1 and 2). LOW findings were not approved
this round and are deferred; INFO items carry no action.

| Change | Decision | By | Note |
|---|---|---|---|
| 1 (F1, HIGH) | approved | Arthur Georges | consequence acknowledged: malformed files previously wiped now load with only the bad uid removed |
| 2 (F2, HIGH) | approved | Arthur Georges | consequence acknowledged: het-free 1-row files that previously "loaded" (corrupt) now error |
| 3 (F3, MEDIUM) | approved | Arthur Georges | consequence: canonical testset now gets plate_location NA instead of header/data garbage |
| 4 (F4, MEDIUM) | approved | Arthur Georges | '_n' suffixes now actually applied |
| 5 (F5, MEDIUM) | approved | Arthur Georges | verbose = 0 now fully silent |
| 6 (F6, LOW) | deferred | Arthur Georges | LOW findings not approved this round |
| 7 (F7, LOW) | deferred | Arthur Georges | LOW findings not approved this round |
| — (F8, INFO) | folded into change 1 | | count-table recomputation shipped with the F1 rewrite |

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-utils.read.dart` (base upstream/dev
ddaed27): changes 1-5 (F1-F5, with F8's count-table recomputation folded
into change 1 as scoped). Changes 6-7 (F6-F7, LOW) deferred - roxygen
untouched, so `man/utils.read.dart.Rd` and NAMESPACE are unchanged
(`devtools::document()` run to confirm).

Verification: baseline characterization test run against the branch -
21 expectations, 0 failures. Diffs from the Phase A baseline all map to
approved findings and are marked `# [approved F<n>]` in the test file:
the triplet-uid file now loads with only the malformed uid removed (F1),
the het-free 1-row file now stops with the format-contradiction error (F2),
`plate_location` on the canonical 3-header-row testset is now NA instead of
pasted header/data garbage while the in-range `service` read is unchanged
(F3), duplicate individuals get the promised `_n` suffix (F4), and
`verbose = 0` prints nothing (F5). End-to-end run on `testset_SNPs_2Row.csv`
at verbose = 3: 2-row format, 250 individuals, 255 SNPs, 510 covmetrics
rows; full `gl.read.dart` chain re-run with the fixed util: 250 x 255,
ploidy 2, `plate_location` NA throughout (approved F3 surface change).

Caller grep across the 8 dartR-verse clones: `utils.read.dart` is called
only by `gl.read.dart` (dartR.base); the chain was re-verified above.
`gl.read.silicodart` does not use this util.

```json
{
  "function": "utils.read.dart",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ddaed27",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT2", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DAT1", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "none:input-bounds", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "applied", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "deferred", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "deferred", "change": 7},
    {"id": "F8", "severity": "INFO", "confidence": "medium", "rule": "STY1", "status": "applied", "change": 1}
  ],
  "coverage_skipped": [
    "2-row silicodart: format does not exist",
    "MarkerName header branch: no fixture",
    "DAT6 FBM: out of scope for this function"
  ],
  "status": "pr-open",
  "pr": null
}
```
