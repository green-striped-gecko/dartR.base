# Characterization tests for utils.plink.run
# Baseline snapshotted before review (infrastructure wave, dev at
# ddaed27). Assertions marked [approved diff] were flipped in Phase C.
# The composed command string is the contract; no PLINK binary needed.

test_that("the composed command is well formed", {
  td <- tempdir()
  # [approved diff, PLINK-call fix] a missing executable now stops with an
  # error that shows the command; baseline returned the command silently
  err <- tryCatch(
    utils.plink.run(dir.in = td, plink.cmd = "nonexistent_exe_999",
                    syntax = "--file hapmap1", verbose = 0),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "PLINK exited with status")
  # [approved diff I11] baseline: the default plink.path produced a
  # literal "path/" prefix, and a missing space glued syntax onto --out:
  # "path/nonexistent_exe_999 --file hapmap1--out hapmap1".
  expect_false(grepl("hapmap1--out", err))  # [approved diff I11]
  expect_false(grepl("path/nonexistent", err))  # [approved diff I11]
})

# Stub PLINK (shell script), so these run without the real binary
make_stub <- function(dir, body) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(dir, "plink")
  writeLines(c("#!/bin/sh", body), f)
  Sys.chmod(f, "755")
  dir
}

test_that("a failing PLINK run stops with PLINK's message", {
  skip_on_os("windows")
  td <- tempfile(); dir.create(td)
  stub <- make_stub(file.path(td, "bin"),
                    c("echo 'Error: simulated failure' >&2", "exit 2"))
  expect_error(
    utils.plink.run(dir.in = td, plink.path = stub, syntax = "--file x",
                    verbose = 0),
    "simulated failure"
  )
})

test_that("executable and output paths with spaces stay single arguments", {
  skip_on_os("windows")
  td <- file.path(tempfile(), "my data"); dir.create(td, recursive = TRUE)
  args.file <- file.path(td, "args.txt")
  stub <- make_stub(file.path(td, "my bin"),
                    c(paste0("for a in \"$@\"; do echo \"$a\"; done > '",
                             args.file, "'"), "echo 'stub log'", "exit 0"))
  out <- capture.output(
    cmd <- utils.plink.run(dir.in = td, plink.path = stub, out = "my out",
                    syntax = paste("--file", shQuote("my file")),
                    verbose = 0)
  )
  expect_length(out, 0L)
  args <- readLines(args.file)
  expect_equal(args[which(args == "--out") + 1], "my out")
  expect_equal(args[which(args == "--file") + 1], "my file")
  expect_true(any(grepl("stub log", capture.output(
    utils.plink.run(dir.in = td, plink.path = stub, syntax = "--file x",
                    verbose = 3)
  ))))
})
