test_that("get_residuals is standardized", {
  # given
  pheno_file <- system.file("testdata", "test_h2_0.5.pheno", package="mmer")
  covariate_file <- system.file("testdata", "test.cov", package="mmer")

  # when
  residuals <- get_residuals(pheno_file, covariate_file)
  
  # then
  expect_equal(mean(residuals), 0, tolerance = 1e-10)
  expect_equal(sd(residuals), 1, tolerance = 1e-10)
})

test_that("adjusting phenotypes writes to file", {
  # given
  pheno_file <- system.file("testdata", "test_h2_0.5.pheno", package="mmer")
  covariate_file <- system.file("testdata", "test.cov", package="mmer")
  
  # when
  pheno_data <- adjust_phenotype(pheno_file, covariate_file)
  
  # then
  expect_equal(nrow(pheno_data), 200)
  expect_equal(ncol(pheno_data), 3)
})
