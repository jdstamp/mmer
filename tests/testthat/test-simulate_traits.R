test_that("simulate_traits works", {
  # given
  plink_file <-
    gsub("\\.bed", "", system.file("testdata", "test.bed", package = "famer"))
  temp_dir <- tempdir(check=TRUE)
  out_file <- paste0(temp_dir, "/test.phen")
  heritability <- 0.5
  rho <- 0.8
  n_additive_snps <- 50
  gxg_indices <- 1:6
  target_mean <- 0
  target_var <- 1
  column_names <- "FID IID TRAIT"
  # when
  sim <- simulate_traits(plink_file,
                         out_file,
                         heritability,
                         rho,
                         n_additive_snps,
                         gxg_indices)
  lines <- readLines(out_file, n = 1)
  # then
  expect_equal(lines[1], column_names)
  on.exit(file.remove(out_file))
})
