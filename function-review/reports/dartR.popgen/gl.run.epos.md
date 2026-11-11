# Review: gl.run.epos (dartR.popgen)
- Family mode: analysis (wrapper around EPOS, epos2plot and bootSfs)
- Custodian: Bernd Gruber (STY5: this report is the discussion record;
  changes approved by Luis)
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 1e90a1e (origin/dev)
- Datasets: possums.gl[1:30, 1:200] after gl.filter.monomorphs (30
  individuals, 176 loci), and without the filter (24 monomorphic loci) for
  the zero-class check; testset.gs[1:20, 1:60]
- External software: EPOS, epos2plot (`~/programs`) and bootSfs from the
  dartRverse `binaries/epos_mac.zip` (same `epos` binary); EPOS example data
  from github.com/EvolBioInf/epos (`testF.dat`, `fig2bFolded.sfs`,
  `fig2bUnfolded.sfs`) for the input format; macOS
- Baseline: tests/testthat/test-gl.run.epos.R (11 tests, 28 expectations;
  9 tests need EPOS_DIR, the bootstrap test also needs bootSfs)

## Verdict

**Standards: Needs work** — failures stop with an empty message or an
unrelated R error, SilicoDArT data is accepted, "Output written" prints at
`verbose = 0`, and there is no way to set a seed.
**Spec: Rework** — every setting other than the default `minbinsize = 1,
folded = TRUE` sends EPOS a wrongly labelled SFS: `folded = FALSE` runs as
a folded SFS of twice the sample size, `minbinsize = 2` shifts every class
down by one, and `minbinsize = 0` (the documented alternative to `L`)
always fails. `upper` and `lower` are never passed to epos2plot.

What works well: with the defaults the SFS matches the EPOS example format
(classes 1..n/2 for n sequences) and the diagnostics (sites, likelihood,
d^2, observed/expected SFS) are parsed correctly; the run folder is a
`tempfile()` and `cleanup` removes only that folder.

## Findings

**F1 [HIGH, confidence: high] — `folded = FALSE` is run as a folded SFS with twice the sample size (DOC5)**
`R/gl.run.epos.R:167-178` — the unfolded SFS from `gl.sfs` (classes
1..2n) is sent without `-U`, the option the documentation names.
Failure scenario (baseline test 7): 30 individuals give classes 1..60;
EPOS reads them as a folded SFS of 120 sequences and returns Ne 15,000.
Sent correctly (classes 1..59, `-U`), EPOS returns Ne 9,440.
Proposed change: add `-U` and send classes 1..(2n - 1) (the class 2n holds
loci fixed for the derived allele).
**Consequence: estimates change for every call with `folded = FALSE`.**

**F2 [HIGH, confidence: high] — `minbinsize >= 2` shifts the SFS by one class per dropped class (DOC5)**
`R/gl.run.epos.R:168` — the SFS is written as `r = 1:length(sfs)`, so
after `gl.sfs` drops classes below `minbinsize` the counts are relabelled
from 1.
Failure scenario (baseline test 6): with `minbinsize = 2` doubletons are
sent as singletons (EPOS sees classes 1..29): Ne 29,000 from 77,300
generations. Sending all classes and excluding singletons with EPOS's own
`-x 1` gives Ne 31,100 from 83,100 generations.
Proposed change: write each count with its true class number and exclude
classes below `minbinsize` with `-x`.
**Consequence: estimates change for every call with `minbinsize >= 2`.**

**F3 [HIGH, confidence: high] — `minbinsize = 0` (zero class instead of `L`) always fails (DOC5)**
`R/gl.run.epos.R:168,174` — the zero class is written as class 1 and `-l`
is dropped, so EPOS finds neither a zero class nor a sequence length.
Failure scenario (baseline test 5): possums.gl[1:30, 1:200] (24
monomorphic loci) stops with "ERROR[epos]: Please include either the
zero-clas in the SFS or enter the sequence length", then R's "'names'
attribute [4] must be the same length as the vector [1]". Written as
class 0, EPOS runs. Note: the zero class counts only the monomorphic loci
kept in the genlight (24 here against L = 1e5), so it rarely stands in for
`L` with DArT data; the documentation should say so.
Proposed change: write the zero class as class 0 (fixed together with F2
by writing true class numbers); document when it is meaningful.
**Consequence: `minbinsize = 0` returns results instead of an error.**

**F4 [MEDIUM, confidence: high] — `upper` and `lower` are ignored (DOC5)**
`R/gl.run.epos.R:192` — `epos2plot` is called without `-u`/`-l`.
Failure scenario: `boot = 100, upper = 0.75, lower = 0.25` returns the
95% limits (epos2plot defaults 0.975/0.025), labelled by the user as 50%
limits. dartR Shiny exposes both arguments.
Proposed change: pass `-u upper -l lower` to epos2plot.
**Consequence: `low`/`high` change for calls with non-default `upper` or
`lower`.**

**F5 [MEDIUM, confidence: high] — command built with empty options (FS5, DOC5)**
`R/gl.run.epos.R:174-178` — `L = NULL` gives `-l ` with no value,
`u = NULL` gives `-u ` with no value, and `other.options` is pasted with
no leading space.
Failure scenario (baseline test 8): each of the three stops with EPOS
"efopen(...) failed", then R's "'names' attribute [4] must be the same
length as the vector [1]". The documentation says that without `u` the
EPOS default (5e-9) is used; the example option style ("-m, -x") fails.
Proposed change: omit `-u` when `u` is NULL; stop with a message when `L`
is NULL and the SFS has no zero class; add a space before `other.options`.

**F6 [MEDIUM, confidence: high] — failures stop with empty or unrelated messages (FS5, VRB2)**
`R/gl.run.epos.R:146-161,186-200` — missing binaries print with `cat()`
then `stop()` with no message; exit statuses of bootSfs/epos/epos2plot
are never checked, so any EPOS error surfaces as the `names` error above.
`chmod 777` runs on Linux only.
Failure scenario (baseline test 1): a wrong `epos.path` gives `Error:`
with an empty message (the reason is on stdout, lost in dartR Shiny).
Proposed change: `stop(error(...))` naming the missing files and
`gl.download.binary("epos")`; check each exit status and stop with EPOS's
message; `Sys.chmod(0755)` on Linux and macOS.

**F7 [MEDIUM, confidence: high] — SilicoDArT data is accepted (FS4)**
`R/gl.run.epos.R:128` — `datatype` is computed and never used.
Failure scenario (baseline test 9): `testset.gs` runs and returns a
history built from presence/absence scores.
Proposed change: stop on SilicoDArT input.
**Consequence: SilicoDArT input errors.**

**F8 [MEDIUM, confidence: high] — an unnamed user `sfs` fails (DOC5)**
`R/gl.run.epos.R:252-253` — class numbers are read from names like `d1`.
Failure scenario (baseline test 10): `sfs = c(6, 6, 3, ...)` (a plain
vector, as from another program) stops with "arguments imply differing
number of rows: 0, 30" after EPOS has run.
Proposed change: take class numbers from `d<k>` names when present,
otherwise from `minbinsize` upwards; document the format.

**F9 [LOW, confidence: high] — runs cannot be made reproducible (principle: reproducibility)**
`R/gl.run.epos.R:172,178` — bootSfs and epos seed from the system clock.
Failure scenario (probe): two identical calls with `boot = 5` return
different histories.
Proposed change: add `seed = NULL`; when set, pass `-s seed` to epos and
bootSfs.
**Consequence: new argument; default behaviour unchanged.**

**F10 [LOW, confidence: high] — messaging ignores verbosity; SFS plot ignores `plot.theme` (FS3, VRB1, VRB3, PLT1)**
`R/gl.run.epos.R:122-124,247,261` — `utils.flag.start` is commented out,
"Output written to" prints at every verbosity, EPOS warnings are shown
unconditionally, `sfs_plot` uses `theme_bw()`.
Failure scenario (baseline test 4): `verbose = 0` prints the output path.
Proposed change: `utils.flag.start`; gate messages at `verbose >= 2`;
`sfs_plot` uses `plot.theme`.

**F11 [LOW, confidence: high] — documentation gaps and mismatches (DOC1, DOC5, DOC7 (proposed rule))**
`R/gl.run.epos.R:1-86` — no `@name`/`@title`/`@family`; `@author` has no
Author(s) part; `outfile` default documented as 'genepop.gen' (it is
'epos.out'); `@return` lists history as "generation, median, low and high"
(the order is generation, low, median, high); the EPOS description
paragraph is repeated; `-l` and the zero class are not explained;
typos (programm, paramter, "provide via the sfs parameter (see below)").
Failure scenario: a user looks for `genepop.gen` or reads the columns in
the documented order.
Proposed change: fix the items above and add the zero-class and
`other.options` notes.

**F12 [INFO, confidence: high] — EPOS doubles the middle class of a folded SFS on input**
EPOS reports 180 polymorphic sites for 176 loci: the class n/2 (4 loci)
is counted twice. EPOS's own example `testF.dat` shows the same (24,148
sites in the file, 24,384 reported). This is EPOS's convention, not a
wrapper defect. No change.

## Proposed changes

1. Send EPOS the SFS with its true class numbers: class 0 when kept,
   `-x` for classes below `minbinsize`, `-U` and classes 1..(2n - 1) for
   `folded = FALSE` (F1, F2, F3).
   **Consequence: estimates change for `folded = FALSE` and for
   `minbinsize >= 2`; `minbinsize = 0` returns results instead of an
   error. Default calls are unchanged.**
2. Pass `upper`/`lower` to epos2plot (F4).
   **Consequence: `low`/`high` change for non-default `upper`/`lower`.**
3. Command and failure handling: omit `-u` when NULL, require `L` without
   a zero class, space before `other.options`, clear errors for missing
   binaries and EPOS failures, `Sys.chmod(0755)` (F5, F6).
4. Reject SilicoDArT input (F7).
   **Consequence: SilicoDArT input errors.**
5. Accept an unnamed user `sfs` (F8).
6. New `seed` argument passed to epos and bootSfs (F9).
   **Consequence: new argument; default behaviour unchanged.**
7. Verbosity and `sfs_plot` theme (F10).
8. Documentation fixes (F11).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run
- Spec: behaviour vs roxygen with the real binaries (defaults,
  `minbinsize` 0 and 2, `folded = FALSE`, `L`/`u` missing,
  `other.options`, `boot = 5`, `upper`/`lower`, SilicoDArT, unnamed
  `sfs`, wrong path, unknown method, `verbose = 0`, cleanup) — run
- Correct-input comparison: EPOS run directly on correctly labelled SFS
  files (`-x 1`; `-U` with classes 1..59; class 0) — run
- Callers (API3): dartR Shiny (`dartr2shiny/shiny_fun/Fun_gl.run.epos.R`)
  passes `minbinsize`, `folded`, `upper`, `lower`, `method`, `depth`,
  `outpath`, `outfile`, `plot.theme`; no sibling dartR.* package calls it
  (dartRverse only documents it in `gl.download.binary`) — run
- Windows and Linux: SKIPPED — not run; Windows paths (`.exe`, dlls) read
  only
- `method = "exhaustive"`: not run (command construction read only)
- DAT1-DAT4, FS8: not applicable — no genlight returned
- DAT6 (FBM): SKIPPED — no FBM fixture; `gl.sfs` densifies, to be covered
  in the `gl.sfs` review
- Missing data in the SFS: same issue as gl.run.stairway2 F12, left to the
  `gl.sfs` review
- Google Group / GitHub issues: not searched in this session

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | |
| 2 | approved | Luis | |
| 3 | approved | Luis | |
| 4 | approved | Luis | |
| 5 | approved | Luis | |
| 6 | approved | Luis | |
| 7 | approved | Luis | |
| 8 | approved | Luis | |

## Outcome

Branch `review-epos` (from origin/dev 1e90a1e). Evidence:
`tests/testthat/test-gl.run.epos.R`, 11 tests, 37 expectations, all passing
with EPOS, epos2plot and bootSfs (dartRverse `epos_mac.zip`, macOS);
without EPOS_DIR 3 tests run and 8 skip. Every changed assertion is tagged
with its approved change.

- 1 applied: the SFS is written with true classes (class 0 when
  `minbinsize = 0`; all classes from 1 with `-x` below `minbinsize`; `-U`
  and classes 1..2n-1 when unfolded). Each result equals EPOS run directly
  on a correctly labelled file (seed 1): `minbinsize = 2` Ne 31,100 /
  83,100 generations (was 29,000 / 77,300); `folded = FALSE` Ne 9,440 /
  37,100 (was 15,000); `minbinsize = 0` Ne 1.5e7 with 24 monomorphic
  sites (was an error). Side effect: the returned `sfs` no longer carries
  `d1..dn` row names (values identical).
- 2 applied: `-l lower -u upper` passed to epos2plot; with `boot = 5,
  seed = 3`, 0.25/0.75 limits lie inside the 0.025/0.975 limits (test 11).
- 3 applied: `-u` omitted when NULL; `L` required unless `minbinsize = 0`;
  programs called by full path with `system2()`, arguments as a vector (so
  `other.options` needs no leading space); a non-zero exit status stops
  with the program's message ("epos failed (exit status 139)" for
  `-L abc`); missing binaries named with a `gl.download.binary` pointer;
  `Sys.chmod(0755)`.
- 4 applied: SilicoDArT errors (test 9).
- 5 applied: an unnamed `sfs` gives the same history as the computed one
  (test 10).
- 6 applied: `seed` passed as `-s` to epos and bootSfs; two seeded
  bootstrap runs identical (test 11).
- 7 applied: `utils.flag.start`; messages at `verbose >= 2`, EPOS stderr
  at `verbose >= 3`; `sfs_plot` uses `plot.theme`. `verbose = 0` prints
  nothing (test 4).
- 8 applied: roxygen rewritten; `devtools::document()` also rewrote seven
  unrelated Rd files (stale "population structure" family links already
  on dev); those were reverted, only `gl.run.epos.Rd` and
  `gl.run.stairway2.Rd` (new shared family link) are committed.
- Default calls unchanged: with seed 1 the history and diagnostics are
  identical to the pre-change run.
- `devtools::check()`: 0 errors; install WARNING from packages built under
  R 4.4.3; 3 NOTEs (hidden `.git`, time check, NEWS parse) present before
  this change.
- NEWS entry added. Callers: dartR Shiny passes `minbinsize`, `folded`,
  `upper`, `lower` (now effective); it does not read `sfs` row names.
- Windows and Linux not run.
- PR: #105

## Addendum (2026-11-11, follow-up review)

- Origin: existing; follow-up after reading the function against the EPOS
  1.7 and epos2plot 1.2 manuals, Lynch et al. (2020) and the "Running Ne
  methods" workshop slides
- Reviewer: Claude (claude-fable-5-2), dartr-function-review v3.0.0
- Package commit: 8262192 (origin/dev); `R/gl.run.epos.R` unchanged since
  PR #105
- Datasets: platypus.gl TENTERFIELD (41 individuals; 484 polymorphic loci
  with call rate >= 0.95, 104 of them with missing calls);
  possums.gl[1:30, 1:200] after gl.filter.monomorphs (176 loci);
  possums.gl[c(1:15, 31:45), ] (two populations); the EPOS example
  `data/testF.dat` (folded SFS, 30 sequences); simulated folded SFS for 50,
  100 and 200 individuals (search path and timing only)
- External software: epos, epos2plot and bootSfs in `~/programs` (the
  binaries of the first review), macOS arm64
- Baseline: five tests appended to tests/testthat/test-gl.run.epos.R, tagged
  [follow-up baseline]; 16 tests, 54 expectations, all passing on 8262192
  before any change

### Verdict

**Standards: Ready** — structure, documentation tags, guards and verbosity
gating conform since PR #105; the one gap, an exclusion that changes the
result with no message, is F13.
**Spec: Needs work** — with missing calls the estimates are too low by the
share of loci dropped and nothing says so; Ne = 1 placeholders read as
estimates; `method = "exhaustive"` runs the greedy search unless `depth` is
raised.

### Findings

**F13 [HIGH, confidence: high] — loci with missing calls are dropped without a message and `L` still counts them (DOC5; VRB4 (proposed rule))**
`R/gl.run.epos.R:196,237,45-48` — since PR #106 `gl.sfs` builds the
spectrum only from loci scored in every individual. `gl.run.epos` calls it
with `verbose = 0`, so the exclusion message never prints, and `-l L` is
passed unchanged. Ne and time both scale with 1/(u L) (probe: `L` x 10
divides both by 10), so both are too low by the share of loci dropped. The
help text still describes the earlier behaviour ("counted in a lower
frequency class").
Failure scenario (follow-up baseline test 12): platypus.gl TENTERFIELD,
0.7% missing genotypes: 104 of 484 loci are dropped with no message at
`verbose = 2`; with `L = nLoc(x) * 69` the history is Ne 8,600 then 76,700.
On the 380 complete loci with `L = 380 * 69` it is 11,000 then 97,700; the
first run is 22% lower (380/484).
Proposed change: (a) report the number of loci excluded at `verbose >= 1`;
(b) when the SFS is computed from `x` and `L` is used, pass EPOS
`L` x (loci used / loci in `x`), rounded, and state the value in the same
message; (c) rewrite the help text: `L` refers to the loci in `x`. A
user-supplied `sfs` leaves `L` as given.
Alternative to (b): leave `L` as given and print the factor for the user to
apply.
**Consequence of (b): estimates rise for every call on data with missing
calls.**

**F14 [MEDIUM, confidence: high] — Ne = 1 placeholders are returned and plotted as estimates (principle: silent API trap; VRB4 (proposed rule))**
`R/gl.run.epos.R:258-269` — EPOS sets a negative size estimate to 1
(manual p. 2: "negative values are set to the smallest possible population
size, N = 1"). The wrapper passes these on to `history` and to the plot
with no notice.
Failure scenario (follow-up baseline test 13; also baseline test 3): the
default run on 30 possums returns Ne = 1 from 0 to 1.27 generations ago,
then 29,700; with `boot = 200` the median there is still 1. The plot shows a
collapse to one individual at present. The Ne = 1 interval does not move
when `L` changes (1.27 generations for `L` 1e5 and 1e6), unlike every
estimate. Each of the three searches (greedy, `-E 3`, `-E 4`) on a
simulated spectrum of 50 individuals returned a level with N = 1.
Proposed change: after parsing, if `median` (or, for bootstrap runs, `low`)
equals 1 anywhere, warn at `verbose >= 1` with the range of generations,
stating that 1 is the value EPOS writes for a size it could not estimate;
add the same to `@details`. `history` and the plot are not changed.

**F15 [MEDIUM, confidence: high] — `method = "exhaustive"` with the default `depth = 2` is the greedy search (DOC5)**
`R/gl.run.epos.R:75-76,128,241` — `depth` is passed as `-E`. EPOS help:
"-E > 2 differs from greedy"; the first two rounds are the same under both
strategies (manual p. 6).
Failure scenario (follow-up baseline test 14): on the EPOS example
`testF.dat`, `method = "exhaustive"` returns the greedy history
(log-likelihood -106.137, levels {2,11}); `depth = 3` gives -99.505 and a
different history. dartR Shiny offers the method as a button with depth 2
preset.
Proposed change: document that `depth` must be 3 or more, that the rounds
after `depth` continue greedily, and the cost; warn at `verbose >= 1` when
`method = "exhaustive"` runs with `depth < 3`. The default stays 2.
Not proposed: a default of 3. On a simulated spectrum of 50 individuals
`-E 3` ends at log-likelihood -425.8 where greedy reaches -364.0, because
the greedy rounds that follow start from another three-level
configuration ({2,4,18} instead of {2,3,28}). `-E 4` reaches -362.8 but
takes 33 s for 50 individuals and 469 s for 100 (greedy 0.4 s and 1 s).

**F16 [LOW, confidence: high] — populations are pooled without a message; the example pools ten (DOC5)**
`R/gl.run.epos.R:196,58-61,107` — the SFS is always computed with
`singlepop = TRUE`. The method assumes "no substantial population
structure" (Lynch et al. 2020, p. 211). The help text says only that a
multidimensional SFS is not supported, and `@examples` runs all of
possums.gl (ten populations).
Failure scenario (follow-up baseline test 16): two possum populations
pooled at `verbose = 2`: no message.
Proposed change: message at `verbose >= 2` when `nPop(x) > 1`; state the
assumption in `@details`; run the example on one population.

**F17 [LOW, confidence: high] — an unfolded `sfs` from `gl.sfs` keeps the class of fixed loci (DOC5)**
`R/gl.run.epos.R:204-205` — the computed path drops class 2n (loci fixed
for the alternative allele); a user-supplied named `sfs` does not, so EPOS
reads 2n + 1 sequences.
Failure scenario (follow-up baseline test 15):
`sfs = gl.sfs(x, folded = FALSE, singlepop = TRUE)` with `folded = FALSE`
on 30 possums sends classes 1..60 and returns Ne 9,400 from 37,000
generations; the computed path sends 1..59 and returns 9,440 from 37,100.
Proposed change: with `folded = FALSE`, drop class `2 * nInd(x)` from a
user-supplied `sfs` when present, and say so in the `sfs` parameter text.
**Consequence: estimates change slightly for calls that pass an unfolded
`gl.sfs` result.**

**F18 [LOW, confidence: high] — the reference is dated 2019 (DOC5)**
`R/gl.run.epos.R:5,8,112` — Lynch et al. is G3 10(1): 211-223, January
2020; the text cites "Lynch et al. 2019".
Failure scenario: a user cites the wrong year.
Proposed change: 2020 in the three places.

### Proposed changes

9. Report the loci excluded for missing calls at `verbose >= 1`, scale `L`
   to the loci used when the SFS is computed from `x`, and rewrite the help
   text (F13).
   **Consequence: estimates rise for data with missing calls, by loci in
   `x` / loci used; calls on complete data are unchanged.**
10. Warn at `verbose >= 1` when the history contains Ne = 1, and document
    it (F14).
11. Document `depth >= 3` and warn at `verbose >= 1` when
    `method = "exhaustive"` runs with `depth < 3`; default unchanged (F15).
12. Message for pooled populations at `verbose >= 2`, the single-population
    assumption in `@details`, example on one population (F16).
13. Drop the class of fixed loci from a user-supplied unfolded `sfs` (F17).
    **Consequence: estimates change slightly for those calls.**
14. Correct the year of the reference to 2020 (F18).

### Coverage

- Standards walk: FS, DOC, VRB, DEP, PLT, STY, API on the file as of
  8262192 — run; no finding beyond F13 (VRB4)
- Spec: behaviour against the roxygen, the EPOS 1.7 manual, the epos2plot
  1.2 manual and Lynch et al. (2020) — run
- `L` scaling (`L` x 10): Ne and generations divide by 10; the Ne = 1
  interval stays — run
- `method = "exhaustive"` with `depth` 2, 3, 5 on `testF.dat`; `depth` 3
  and 4 timed on simulated spectra (200 individuals: `depth` 3 only) — run
- `L` below the number of SNPs: EPOS's error is passed on, no finding — run
- Callers (API3): dartR Shiny
  (`dartr2shiny/shiny_fun/Fun_gl.run.epos.R`) passes `L` from a numeric
  field, `method` and `depth` (preset 2), and reads `$history`; no sibling
  dartR.* package calls the function — run
- GitHub issues of dartR.popgen for "epos": none found
- Google Group: SKIPPED — not searched in this session
- Windows and Linux: SKIPPED — not run
- DAT6 (FBM): SKIPPED — the densification is in `gl.sfs`
- Workshop slides: text read; figures not viewed

### Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 9 | approved | Luis | "report and scale L" chosen over "report only"; the rise in estimates for data with missing calls approved in those words |
| 10 | approved | Luis | |
| 11 | approved | Luis | default `depth` stays 2 |
| 12 | approved | Luis | |
| 13 | approved | Luis | the slight change in estimates for those calls approved in those words |
| 14 | approved | Luis | |

### Outcome

Branch `review-epos-followup` (from origin/dev 8262192). Evidence:
`tests/testthat/test-gl.run.epos.R`, 18 tests, 78 expectations, all passing
with epos, epos2plot and bootSfs (macOS, `EPOS_DIR`); without the binaries 3
tests run and 15 skip. The 11 tests of the first review pass unchanged.
Every changed assertion is tagged with its approved change.

- 9 applied: the number of loci used is the sum of the spectrum that
  `gl.sfs` returns; when loci were excluded, a warning prints at
  `verbose >= 1` and EPOS is given `round(L * loci used / nLoc(x))`.
  platypus.gl TENTERFIELD: "104 of 484 loci have missing calls and were
  excluded from the sfs ... L scaled to the loci used: 33396 x 380 / 484 =
  26220"; EPOS reports 26,220 sites; the history is Ne 11,000 then 97,700
  (was 8,600 then 76,700) and is identical to the run on the 380 complete
  loci with `L = 380 * 69`. Complete data: no message, `L` as given. With
  `minbinsize = 0` the message prints and nothing is scaled. A user-supplied
  `sfs` leaves `L` as given. `verbose = 0` prints nothing.
- 10 applied: "the history has Ne = 1 for generations 0 to 1.27 before
  present. EPOS writes 1 where it could not estimate a population size ..."
  at `verbose >= 1`; a second wording when only the lower bootstrap limit is
  1. `history` and the plot are unchanged (same values as the baseline); no
  warning for a history without Ne = 1 (`testF.dat`).
- 11 applied: warning at `verbose >= 1` when `method = "exhaustive"` runs
  with `depth < 3`; none for `depth = 3` or for greedy. Results unchanged:
  on `testF.dat` greedy and exhaustive with the default depth are identical
  (-106.137), `depth = 3` gives -99.505. `method` and `depth` rewritten in
  the help.
- 12 applied: "x has 2 populations, which are pooled into one sfs ..." at
  `verbose >= 2`; the assumption stated in `@details`; the example runs
  possums.gl[1:30, ].
- 13 applied: with `folded = FALSE`, a last class equal to `2 * nInd(x)` is
  dropped from a user-supplied `sfs`, named or unnamed. The unfolded
  `gl.sfs` result now gives the history of the computed path (Ne 9,440; was
  9,400). An unfolded `sfs` of another sample size is sent as given.
- 14 applied: 2020 in the description (twice) and in `@references`.
- `devtools::document()` changed only `man/gl.run.epos.Rd`. `checkRd`: four
  "Lost braces" notes on the `@return` items, present before this change.
- End to end at `verbose = 3` on the platypus and possum data, with
  `plot.display = TRUE` and `plot.file`: runs, plot saved.
- `devtools::check(args = "--no-manual")`: 0 errors; the install WARNING
  (packages built under R 4.4.3) and 2 NOTEs (hidden `.git` of the
  worktree, time check) do not come from this change.
- NEWS entry added. Callers: dartR Shiny passes `L`, `method` and `depth`;
  its estimates rise for data with missing calls as for any caller. Its
  depth field is preset to 2, so "exhaustive" there still runs the greedy
  search; the new warning goes to the console only. That is for dartR Shiny
  to settle and is not changed here.
- While applying: `minbinsize = 0` on a genlight without monomorphic loci
  stops with EPOS's own message ("Please include either the zero-clas in
  the SFS or enter the sequence length"). The message is passed on as
  designed in PR #105; no change.
- Windows and Linux not run.
- PR: #135

## Machine block

```json
{
  "function": "gl.run.epos",
  "package": "dartR.popgen",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "1e90a1e",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 3},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 3},
    {"id": "F7", "severity": "MEDIUM", "confidence": "high", "rule": "FS4", "status": "approved", "change": 4},
    {"id": "F8", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 5},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "principle: reproducibility", "status": "approved", "change": 6},
    {"id": "F10", "severity": "LOW", "confidence": "high", "rule": "VRB1", "status": "approved", "change": 7},
    {"id": "F11", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 8},
    {"id": "F12", "severity": "INFO", "confidence": "high", "rule": "none (EPOS convention)", "status": "no-change", "change": null},
    {"id": "F13", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 9},
    {"id": "F14", "severity": "MEDIUM", "confidence": "high", "rule": "principle: silent API trap", "status": "approved", "change": 10},
    {"id": "F15", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 11},
    {"id": "F16", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 12},
    {"id": "F17", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 13},
    {"id": "F18", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 14}
  ],
  "followup": {"date": "2026-11-11", "skill_version": "3.0.0", "commit": "8262192", "verdict_standards": "ready", "verdict_spec": "needs_work", "previous_pr": 105},
  "coverage_skipped": ["Windows/Linux not run", "DAT6: no FBM fixture", "Google Group not searched"],
  "status": "pr-open",
  "pr": 135
}
```
