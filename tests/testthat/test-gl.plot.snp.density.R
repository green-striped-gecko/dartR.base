# Characterization tests for gl.plot.snp.density
# Baseline snapshotted before review (dev_luis at 562befa, origin/dev
# fadab59 merged). Assertions tagged [approved diff, change n] were flipped
# in Phase C to reflect the approved behaviour changes (report:
# function-review/reports/dartR.base/gl.plot.snp.density.md).

t1 <- platypus.gl
t1$chromosome <- t1$other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1
t1$position   <- t1$other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1

quiet_plot <- function(...) {
  pdf(NULL)
  on.exit(dev.off())
  out <- capture.output(res <- gl.plot.snp.density(..., verbose = 0))
  list(res = res, out = out)
}

# independent binning of the same data
indep_bins <- function(x, bin.size, min.snps, min.length) {
  df <- data.frame(chr = as.character(x@chromosome), pos = x@position)
  df <- df[!is.na(df$chr) & df$chr != "" & !is.na(df$pos) & df$pos > 0, ]
  n <- tapply(df$pos, df$chr, length)
  mx <- tapply(df$pos, df$chr, max)
  keep <- names(n)[n >= min.snps & mx >= min.length]
  df <- df[df$chr %in% keep, ]
  df$bin <- floor(df$pos / bin.size) * bin.size + bin.size / 2
  list(chr = keep, n = n[keep], mx = mx[keep],
       counts = as.data.frame(table(chr = df$chr, bin = df$bin)))
}

test_that("platypus: returns a ggplot invisibly, binned counts match an independent computation", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  v <- withVisible(gl.plot.snp.density(t1, bin.size = 5e6, min.snps = 10,
                                       min.length = 2e6, verbose = 0))
  expect_false(v$visible)
  p <- v$value
  expect_s3_class(p, "ggplot")
  d <- p$data
  expect_equal(colnames(d), c("chr_label", "bin_center", "n_snps"))
  ib <- indep_bins(t1, 5e6, 10, 2e6)
  expect_equal(nlevels(d$chr_label), 25)
  expect_length(ib$chr, 25)
  expect_equal(sum(d$n_snps), 888)
  expect_equal(sum(ib$n), 888)
  # [approved diff, change 1] baseline: 331 rows, only occupied bins
  expect_equal(nrow(d), sum(floor(ib$mx / 5e6) + 1))   # [approved diff, change 1]
  expect_equal(range(d$n_snps), c(0, 13))              # [approved diff, change 1]
  # every occupied bin agrees with the independent count
  ibc <- ib$counts[ib$counts$Freq > 0, ]
  dd <- d[d$n_snps > 0, ]
  key_d <- paste(sub(" \\(.*$", "", as.character(dd$chr_label)), dd$bin_center)
  key_i <- paste(as.character(ibc$chr), as.numeric(as.character(ibc$bin)))
  expect_setequal(key_d, key_i)
  expect_equal(dd$n_snps[match(key_i, key_d)], ibc$Freq)
})

test_that("chromosome labels carry the SNP count and the last SNP position in Mb", {
  r <- quiet_plot(t1, bin.size = 5e6, min.snps = 10, min.length = 2e6)
  lv <- levels(r$res$data$chr_label)
  expect_true("NC_041728.1_chromosome_1 (83 SNPs, 186.2 Mb)" %in% lv)
  expect_true("NC_041753.1_chromosome_X5 (23 SNPs, 69.3 Mb)" %in% lv)
  r2 <- quiet_plot(t1, bin.size = 5e6, min.snps = 10, min.length = 2e6,
                   chr.info = FALSE)
  expect_true("NC_041728.1_chromosome_1" %in% levels(r2$res$data$chr_label))
})

test_that("chromosomes are ordered alphabetically, numeric-aware, first name at the top", {
  # [approved diff, change 2] baseline: sort(unique(label), decreasing =
  # TRUE), a plain string sort (chr10 between chr1 and chr2); the size
  # ordering promised by @details was dead code.
  r <- quiet_plot(t1, bin.size = 5e6, min.snps = 10, min.length = 2e6,
                  chr.info = FALSE)
  lv <- levels(r$res$data$chr_label)
  # platypus NCBI names: natural and plain order coincide
  expect_equal(lv, sort(lv, decreasing = TRUE))
  expect_equal(lv[1], "NC_041753.1_chromosome_X5")
  expect_equal(lv[25], "NC_041728.1_chromosome_1")
  # unpadded numbers: chr2 before chr10 (reading top to bottom)
  t2 <- t1
  t2$chromosome <- factor(paste0("chr", rep(1:12, length.out = nLoc(t2))))
  r2 <- quiet_plot(t2, bin.size = 5e6, min.snps = 10, min.length = 2e6,
                   chr.info = FALSE)
  lv2 <- levels(r2$res$data$chr_label)
  expect_equal(lv2, paste0("chr", 12:1))                 # [approved diff, change 2]
})

test_that("bins with no SNPs are present with a count of zero", {
  # [approved diff, change 1] baseline: only occupied bins were in the data
  # (minimum n_snps 1, 331 rows), so empty bins were drawn as background.
  r <- quiet_plot(t1, bin.size = 5e6, min.snps = 10, min.length = 2e6)
  d <- r$res$data
  expect_equal(min(d$n_snps), 0)                        # [approved diff, change 1]
  expect_gt(nrow(d), 331)                               # [approved diff, change 1]
  # each chromosome runs from the first bin to the bin of its last SNP
  ib <- indep_bins(t1, 5e6, 10, 2e6)
  for (ch in ib$chr) {
    bc <- d$bin_center[sub(" \\(.*$", "", as.character(d$chr_label)) == ch]
    expect_equal(sort(bc), seq(2.5e6, floor(ib$mx[[ch]] / 5e6) * 5e6 + 2.5e6, by = 5e6))
  }
})

test_that("verbose = 0 is silent; verbose 2 reports filters; verbose 3 prints the table", {
  r <- quiet_plot(t1, bin.size = 5e6, min.snps = 10, min.length = 2e6)
  expect_length(r$out, 0)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  o2 <- capture.output(invisible(gl.plot.snp.density(t1, bin.size = 5e6,
    min.snps = 10, min.length = 2e6, verbose = 2)))
  expect_true(any(grepl("Retained 921 SNPs after initial filtering", o2)))
  # [approved diff, change 4] baseline: no filter summary, verbose 3 == 2
  expect_true(any(grepl("Chromosomes with positioned SNPs: 46", o2)))  # [approved diff, change 4]
  expect_true(any(grepl("dropped 21 with fewer than 10 SNPs and 0 with the last SNP below 2,000,000 bp; retained 25", o2)))
  o3 <- capture.output(invisible(gl.plot.snp.density(t1, bin.size = 5e6,
    min.snps = 10, min.length = 2e6, verbose = 3)))
  expect_gt(length(o3), length(o2))                     # [approved diff, change 4]
  expect_true(any(grepl("NC_041728.1_chromosome_1 +83 +186.2 +38", o3)))
  # the table reads top to bottom: chromosome_1 first, X5 last
  i1 <- grep("NC_041728.1_chromosome_1 ", o3)
  i5 <- grep("NC_041753.1_chromosome_X5 ", o3)
  expect_lt(i1, i5)
})

test_that("input object is not modified", {
  t0 <- t1
  invisible(quiet_plot(t1, bin.size = 5e6, min.snps = 10, min.length = 2e6))
  expect_identical(t1, t0)
})

test_that("an object without chromosome information stops with a clear message", {
  # [approved diff, change 3] baseline: "arguments imply differing number
  # of rows: 0, 255" from data.frame().
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_error(capture.output(gl.plot.snp.density(testset.gl, verbose = 0)),
               "chromosome name and a position for every locus")  # [approved diff, change 3]
})

test_that("SilicoDArT is rejected by the datatype check", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_error(capture.output(gl.plot.snp.density(testset.gs, verbose = 0)),
               "SilicoDArT")
})

test_that("argument checks: bin.size, min.snps, min.length", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_error(capture.output(gl.plot.snp.density(t1, bin.size = 0, verbose = 0)),
               "bin.size must be > 0")
  # [approved diff, change 5] baseline: the messages said "> 1" for a < 1 check
  expect_error(capture.output(gl.plot.snp.density(t1, min.snps = 0, verbose = 0)),
               "min.snps must be >= 1")                  # [approved diff, change 5]
  expect_error(capture.output(gl.plot.snp.density(t1, min.length = 0, verbose = 0)),
               "min.length must be >= 1")                # [approved diff, change 5]
  expect_error(capture.output(gl.plot.snp.density(t1, min.snps = 1e6, verbose = 0)),
               "No chromosomes meet")
})

test_that("FBM-backed object plots identically", {
  tf <- gl.gen2fbm(t1, verbose = 0)
  r <- quiet_plot(t1, bin.size = 5e6, min.snps = 10, min.length = 2e6)
  rf <- quiet_plot(tf, bin.size = 5e6, min.snps = 10, min.length = 2e6)
  expect_equal(rf$res$data, r$res$data)
})

test_that("plot.display controls printing; plot.file saves an RDS; save2tmp is gone", {
  # [approved diff, change 6] baseline: print() unconditional, save2tmp
  # wrote to tempdir() for gl.print.reports().
  f <- names(formals(gl.plot.snp.density))
  expect_true(all(c("plot.display", "plot.file", "plot.dir") %in% f))  # [approved diff, change 6]
  expect_false("save2tmp" %in% f)                                      # [approved diff, change 6]
  td <- tempfile("snpdens"); dir.create(td)
  # whether anything was drawn, read from the device display list; counting
  # png files does not work on Windows, where png() writes a file on
  # dev.off() even when no page was drawn
  drawn <- function(display) {
    pdf(NULL)
    dev.control("enable")
    on.exit(dev.off(), add = TRUE)
    capture.output(invisible(gl.plot.snp.density(t1, bin.size = 5e6,
      min.snps = 10, min.length = 2e6, plot.display = display, verbose = 0)))
    length(recordPlot()[[1]]) > 0
  }
  expect_false(drawn(FALSE))
  expect_true(drawn(TRUE))
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  capture.output(invisible(gl.plot.snp.density(t1, bin.size = 5e6, min.snps = 10,
    min.length = 2e6, plot.display = FALSE, plot.file = "dens", plot.dir = td,
    verbose = 0)))
  expect_true(file.exists(file.path(td, "dens.RDS")))
  expect_s3_class(readRDS(file.path(td, "dens.RDS")), "ggplot")
})
