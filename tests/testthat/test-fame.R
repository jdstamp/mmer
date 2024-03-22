test_that("fame end-to-end", {
  # given
  covariate_file = ""
  plink_file <- gsub("\\.bed", "", system.file("testdata", "test.bed", package="famer"))
  annotation_file <- system.file("testdata", "test_single_annot.txt", package="famer")
  pheno_file <- system.file("testdata", "test_h2_0.5.pheno", package="famer")
  gxg_bin = 0
  n_randvecs = 10
  n_blocks = 10

  snp_indices <- c(3, 8, 9, 13, 16, 19, 29, 34, 93, 97)
  expected <- matrix(
    c(
      0.414218, 0.232245, 0.363624,
      0.425354, -0.0394186, 0.594756,
      0.42462, -0.0568018, 0.607989,
      0.42503, -0.0384796, 0.587712,
      0.424906, 0.0268949, 0.521808,
      0.415192, 0.0974395, 0.43878,
      0.424643, -0.0158059, 0.565519,
      0.42001, 0.146484, 0.390797,
      0.424789, -0.0199609, 0.570214,
      0.423682, 0.0074807, 0.5438
    ), ncol = 3, byrow = TRUE
  )
  # when
  observed <- NULL
  for (i in snp_indices) {
    result <- fame(plink_file,
                   pheno_file,
                   annotation_file,
                   covariate_file,
                   gxg_bin,
                   n_randvecs,
                   i,
                   n_blocks)
    observed <- rbind(observed, result$Est)
  }

  # then
  expect_equal(observed, expected, tolerance = 1e-5)
})

test_that("fame end-to-end with covariate file", {
  # given
  plink_file <- gsub("\\.bed", "", system.file("testdata", "test.bed", package="famer"))
  annotation_file <- system.file("testdata", "test_single_annot.txt", package="famer")
  pheno_file <- system.file("testdata", "test_h2_0.5.pheno", package="famer")
  covariate_file <- system.file("testdata", "test.cov", package="famer")
  gxg_bin = 0
  n_randvecs = 10
  snp_index = 19
  n_blocks = 10

  snp_indices <- c(13, 18, 29, 34, 37, 74, 82, 88, 94, 96)
  expected <- matrix(
    c(
      0.420443, -0.0484235, 0.622903,
      0.419517, -0.0105727, 0.584645,
      0.421285, -0.0511097, 0.624927,
      0.414692, 0.163413, 0.397894,
      0.415631, 0.0506359, 0.529439,
      0.419678, -0.0951365, 0.679929,
      0.422461, 0.13161, 0.432249,
      0.417551, 0.0682423, 0.503778,
      0.400282, 0.166734, 0.413728,
      0.416637, 0.0423788, 0.536325
    ), ncol = 3, byrow = TRUE
  )
  # when
  observed <- NULL
  for (i in snp_indices) {
    result <- fame(plink_file,
                   pheno_file,
                   annotation_file,
                   covariate_file,
                   gxg_bin,
                   n_randvecs,
                   i,
                   n_blocks)
    observed <- rbind(observed, result$Est)
  }

  # then
  expect_equal(observed, expected, tolerance = 1e-5)
})
