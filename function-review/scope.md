# What the function review covers

Written 22 September 2026, when the pending functions were added to
`manifest.csv`. Counts come from `origin/dev` of each package on that date.

## A row in the manifest

One row per function in scope. In scope means either of:

- the function is exported in the package's `NAMESPACE`, or
- its name begins with `utils.`, whether it is exported or not.

That is the scope the campaign had already been using: of the 166 functions
reviewed before this date, all were exported functions or `utils.` helpers.

## Out of scope

Top-level definitions that are neither exported nor named `utils.*`. These are
vendored third-party code (the OutFLANK routines in dartR.popgen), simulation
internals in dartR.sim, package startup helpers in dartRverse, and unexported
helpers elsewhere. **88** definitions across the suite, listed in the last
column below. They can be brought in later; nothing about the scope rule
stops it.

Functions defined inside another function are not rows. Two such names look
like candidates and are deliberately absent: `gl.fbm2gen` (a fallback stub
inside `gl.ibd`, dartR.spatial) and `gl.report.nall.pop` (a helper inside
`gl.report.nall`, dartR.sim).

## Coverage on 22 September 2026

| Package | In scope | Manifest rows | Pending | Out of scope |
|---|---|---|---|---|
| dartR.base | 173 | 182 | 25 | 18 |
| dartR.popgen | 37 | 37 | 30 | 33 |
| dartR.captive | 27 | 27 | 27 | 19 |
| dartR.sim | 13 | 13 | 11 | 9 |
| dartR.spatial | 10 | 10 | 9 | 0 |
| dartR.sexlinked | 5 | 5 | 5 | 3 |
| dartRverse | 2 | 2 | 2 | 6 |
| dartR.data | 0 | 0 | 0 | 0 |
| **total** | **267** | **276** | **109** | **88** |

dartR.base carries nine rows more than it has functions in scope. Seven name
functions that no longer exist: `gl.mahal.assign`, removed in the Mahalanobis
consolidation, and six `utils.*` helpers folded into their callers after being
reviewed. Those rows are kept, since they record work that was done. The other
two are `gl.He` and `gl.Ho`, whose rows sit under dartR.base and whose
`moved-to-dartR.sim` status records where they went.

`utils.read.ped` has a row in each of dartR.base and dartR.popgen. These are
two different implementations of the same name: dartR.base holds a vendored
copy of `snpStats::read.pedfile` in `R/utils.read.ped.r` (252 lines, reviewed),
and dartR.popgen defines a shorter one inside `R/utils.ld.r` (183 lines, not
reviewed). Worth deciding whether the second should call the first, since
dartR.popgen depends on dartR.base.

`dartR.data` carries data only and defines no functions.

## Statuses

`pending` not yet claimed, `in-review` claimed, `awaiting-approval` report
written, `pr-open` with a number, `done` once merged. `moved-to-<package>`
means the function was relocated rather than fixed in place, and `applied`
means changes went in without a recorded pull request.

The `family` column on a pending row is a guess from the function name. Set it
from the skill's family table when the review is claimed.
