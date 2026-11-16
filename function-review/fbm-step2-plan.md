# FBM step 2: file-backed functions read from disk

Branch: `dev_fbm` (from `dev` after the step 1 commit 3391f29).
Source: Luis Mijangos, *File-backed storage (FBM) in dartRverse: summary*,
14 Nov 2026, step 2.

## Goal

A file-backed object can be converted, imputed, analysed by PCA and printed
without holding the decoded genotype matrix in memory. Results stay
**identical** to the gen-backed path, including random draws with the same
seed.

## Status at branch start

| Item | Status |
|---|---|
| Printing (`show`) | Done in step 1: missing data counted in blocks of loci |
| `gl.pcoa` | Core done (dev_luis, 6ed0ea8): cross-product in blocks, locus-mean imputation in a temporary copy. Verification only |
| `gl.gen2fbm` | To do |
| `gl.fbm2gen` | To do |
| `gl.impute` | To do (largest item) |
| Temporary files | To do |

## Principles

1. **One memory knob.** Size every block in bytes with bigstatsr's own
   option `bigstatsr.block.sizeGB` (default 1 GB) through
   `bigstatsr::block_size()`, not with fixed row or column counts.
2. **Read by loci (columns).** The FBM is stored column by column. Reading a
   block of loci is one contiguous read; reading a block of individuals
   touches the whole file.
3. **Missing values are written as code 3.** Assigning `NA` to an
   `FBM.code256` stores 0 (see the note in `gl.impute`, line 207). Every new
   write path must replace `NA` by 3 first.
4. **Never write into the caller's backing file.** FBMs are references;
   write into a copy (as `gl.impute` and `gl.sim.mutate` do now).
5. **Same results in both storage types**, tested for every changed function.

## Work packages

### WP0 Shared helpers (0.5 day)

- `.fbm_tempfile(prefix)`: `tempfile(prefix, tmpdir = getOption("FBM.dir",
  tempdir()))`. bigstatsr's own default already honours `FBM.dir`; dartR's
  explicit `tempfile("geno_")` calls bypass it. Call sites to replace:
  `gl.gen2fbm`, `gl.impute`, `gl.pcoa`, `gl.randomize.snps`, `gl.sample`,
  `[` (subset), `cbind.dartR`, `rbind.dartR` (dartR.base); `gl.sim.mutate`,
  `utils.sims.r` (dartR.sim); and the readers called with `fbm = TRUE`.
- `.fbm_new(n, p, ...)` creating the FBM through `.fbm_tempfile()` and
  registering it for cleanup (WP4).

### WP1 `gl.gen2fbm` (1 day)

Now: each block of loci calls `x[, ind.col]`, which subsets every SNPbin and
so decodes all loci of all individuals for every block (678 Mb for
2,000 x 20,000; a smaller `chunk` is slower and saves nothing).

Plan: work through blocks of **individuals**, which is how SNPbins are stored.
For each block, decode only those SNPbins (`as.integer`), replace `NA` by 3,
and write the rows `G[rows, ]`. The block size is computed from the byte
budget (`rows = budget / (p * 8)`). Each SNPbin is decoded once.

Alternative (Luis's test version): split by loci by slicing the packed raw
bytes of each SNPbin. This is faster but depends on adegenet's internal
SNPbin layout (bits per allele copy, `NA.posi`). Benchmark both and choose.

Deprecate `chunk` in favour of the byte budget, but keep it for
compatibility.

### WP2 `gl.fbm2gen` (0.5 day)

Now: blocks of 256 individuals, decoded as doubles. At 1 million loci that is
2 GB per block, and each block of rows reads across the whole file.

Plan: size the row block from the byte budget. The output (a list of SNPbins)
is in memory by design, so `gl.fbm2gen` is only for data whose gen-backed
form fits in memory. Say so in the help page.

### WP3 `gl.impute` (4 to 5 days)

Now: `as.matrix()` of the whole object or of each population, `seppop()`
plus `rbind()` (each a file copy), and `obj@fbm[] <- matrix` (a full
write). The verbose >= 3 report decodes the matrix twice more.

Plan, per method:

- **frequency / HW.** Allele frequencies per population from counts
  (`big_counts(G, ind.row = pop rows)`), with no decode. Then, for each
  population and each block of loci, read the population's rows, fill the
  missing values and write them into **one** output copy. No `seppop`,
  `rbind` or reorder copies.
- **random / fill.residual.** One pass over blocks of loci.
- **neighbour.** The distances are accumulated over blocks of loci. For the
  Euclidean distance with pairwise-complete loci, used by
  `utils.dist.ind.snp`:
  `sumsq = Z^2 M' + M (Z^2)' - 2 Z Z'`, where `Z` is the block with `NA` set
  to 0 and `M` marks the scored genotypes. Each individual's neighbour order
  depends only on the distances, so filling is then done block by block with
  the precomputed order.
- **beagle.** It already goes through VCF; check that `gl2vcf` does not
  decode the whole FBM.
- Diagnostics: the counts before and after come from `big_counts`.

### WP4 Temporary files (1.5 days)

- Honour `FBM.dir` (WP0).
- Remove backing files that no object uses: register a finalizer
  (`reg.finalizer(fbm, ..., onexit = TRUE)`) on each FBM that dartR creates
  in the temporary folder. It deletes the `.bk` file when the last reference
  is garbage collected. Keep a registry of files that could not be deleted
  and retry later (see caveat 3). Optionally add an exported
  `gl.fbm.cleanup()` for an explicit sweep.
- Never delete a backing file that the user named (`backingfile =`
  argument), only those created through `.fbm_tempfile()`.

### WP5 `gl.pcoa` verification (0.5 day)

- Benchmark against Luis's figures (7 s and 112 Mb for 1,000 x 20,000).
- Check that no step before the FBM block decodes the matrix (the filters
  now read from file; the comment at lines 792-795 is out of date).
- Note: the mean imputation stores the mean rounded to 2 decimals
  (`snp_fastImputeSimple(method = "mean2")`). The difference from glPca is
  tiny but not zero, so tests need a tolerance.

### WP6 Tests and benchmarks (ongoing, 2 days)

- For every changed function, compare gen and FBM with `identical()`. For
  HW, random and fill.residual, also compare with the same `set.seed()`.
- Benchmark script (in `function-review/`, not run on CRAN): simulated
  2,000 x 20,000 and 2,000 x 200,000 data. Record time and peak R memory
  (`gc(reset = TRUE)` then `gc()` "max used").
- Run all examples with `dartR_fbm=true` (Luis's script) before merging.

**Total: about 2 weeks.**
Order: WP0, WP4, WP1, WP2, WP3, WP5, with WP6 alongside.

## Caveats

1. **Identical random draws (gl.impute).** The gen-backed path draws in
   column-major order: per population for HW, over the whole matrix for
   random and fill.residual. A block version is identical only if it keeps
   that order: populations in the outer loop, blocks of loci in the inner
   loop, column-major within a block. For HW that means one pass over the
   file per population. Swapping the loops (one pass over the file) changes
   the random stream: the results are statistically equivalent but not
   identical. **Decision needed.** Splitting `sample(..., replace = TRUE)`
   into pieces consumes the RNG exactly like one call (checked: identical),
   so blocking by itself does not change the draws; only the loop order does.
2. **Row access is slow.** Anything per individual (`gl.fbm2gen`, `gl.gen2fbm`
   by individuals, neighbour imputation) reads across the whole file. Blocks
   must be large, and the timing measured at 200k loci, not only 20k.
3. **Windows cannot delete a memory-mapped file.** A finalizer may run
   before bigstatsr's external pointer releases the mapping, and `unlink()`
   then fails without an error (the same problem `gl.pcoa` works around with
   `gc()`). Hence the retry registry. Test the cleanup on Windows
   specifically.
4. **Hidden decodes.** Helpers called at the end of functions may decode the
   matrix again: `gl.compliance.check`, `utils.reset.flags`, metric
   recalculation, `seppop` (copies), `rbind` (copies). Audit each changed
   function end to end with the benchmark, not only its main loop.
5. **Neighbour imputation keeps an n x n distance matrix in memory** (8 bytes
   per pair: 32 Mb at 2,000 individuals, 3.2 Gb at 20,000). This is
   acceptable for typical DArT data but should be stated in the help page.
6. **Disk use.** Imputation writes a full copy, and `gl.pcoa` makes a
   temporary imputed copy. Without WP4, these accumulate in `tempdir()` in
   long sessions.
7. **Exactness of block sums.** The distance and count sums are sums of small
   integers held in doubles, so block order does not change the result
   (exact below 2^53). This keeps neighbour ties identical. Do not introduce
   scaling or centring before the sums.
8. **Do not use parallel bigstatsr (`ncores > 1`)** in paths with random
   draws, to stay deterministic.
9. **Branch hygiene.** `dev` is being released as 1.2.7. Merge `dev` into
   `dev_fbm` after the release, and merge `dev_fbm` back only for 1.3.0.
   The class change (filtering without copies) is step 3 and stays out of
   this branch.

## Measurements (16 Nov 2026, Ryzen 7 5800H, 32 GB, Windows)

`gl.gen2fbm` variants, 100 Mb block budget, 10 % missing; all three give
byte-identical backing files:

| Size | A current | B blocks of individuals | C raw bytes by loci |
|---|---|---|---|
| 250 x 5,000 | 0.3 s, 69 Mb | 0.1 s, 48 Mb | 0.1 s, 58 Mb |
| 2,000 x 20,000 | 9.8 s, 655 Mb | 2.4 s, 299 Mb | 2.4 s, 348 Mb |
| 2,000 x 200,000 | 71.8 s, 6,372 Mb | 17.7 s, 2,992 Mb | 21.4 s, 1,063 Mb |

Neither B nor C keeps to the 100 Mb budget yet; the peak needs tuning in
the implementation.

Automatic cleanup (finalizer on the FBM environment, Windows): the first
`unlink()` fails because the file is still mapped, and the retry after the
next `gc()` succeeds. A file stays while any copy of the object exists. 8
chained subsets left 9 files without cleanup and 1 with it. Registering
2,000 finalizers takes 0.01 s.

## Decisions (agreed 16 Nov 2026)

1. HW imputation keeps the identical random stream (one pass per
   population).
2. `gl.gen2fbm` uses raw-byte slicing by loci (variant C), with a test
   against decoded genotypes to catch changes in the SNPbin layout.
3. Cleanup is automatic (finalizers with retry registry) plus an exported
   `gl.fbm.cleanup()`.
