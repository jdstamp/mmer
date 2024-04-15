test_that("fame end-to-end no covariates no mask", {
  # given
  plink_file <- gsub("\\.bed", "", system.file("testdata", "test.bed", package="famer"))
  pheno_file <- system.file("testdata", "test_h2_0.5.pheno", package="famer")
  covariate_file = ""
  mask_file <- ""
  n_randvecs = 10
  n_blocks = 10
  rand_seed = 123

  snp_indices <- c(3, 8, 9, 13, 16, 19, 29, 34, 93, 97)

  expected_est <- matrix(
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
  expected_se <- matrix(
    c(
      0.1207738, 0.16184765, 0.15113731,
      0.1211554, 0.02243126, 0.08744855,
      0.1209582, 0.04043788, 0.09893784,
      0.1214335, 0.04355752, 0.09343718,
      0.1207895, 0.05948262, 0.10581516,
      0.1201660, 0.08556147, 0.13006928,
      0.1214145, 0.04715497, 0.09027384,
      0.1219379, 0.16492895, 0.17881043,
      0.1210386, 0.03876772, 0.09262924,
      0.1209970, 0.05537885, 0.09436732
    ), ncol = 3, byrow = TRUE
  )
  id <- sprintf("gxg_%d", snp_indices)
  vc_names <- c("id", "grm", "gxg", "error")
  vc_df <- cbind(id, as.data.frame(expected_est))
  colnames(vc_df) <- vc_names
  vc_df <- pivot_output(vc_df,
                        "component",
                        "vc_estimate",
                        vc_names[2:4])
  se_df <- cbind(id, as.data.frame(expected_se))
  colnames(se_df) <- vc_names
  se_df <- pivot_output(se_df,
                        "component",
                        "vc_se",
                        vc_names[2:4])
  # when
  result <- fame(plink_file,
                 pheno_file,
                 covariate_file,
                 mask_file,
                 n_randvecs,
                 snp_indices,
                 n_blocks,
                 rand_seed)
  observed_est <- result$vc_estimate
  observed_se <- result$vc_se

  # then
  expect_equal(observed_est, vc_df, tolerance = 1e-5)
  expect_equal(observed_se, se_df, tolerance = 1e-5)
})

test_that("fame end-to-end with covariate file no mask", {
  # given
  plink_file <- gsub("\\.bed", "", system.file("testdata", "test.bed", package="famer"))
  pheno_file <- system.file("testdata", "test_h2_0.5.pheno", package="famer")
  covariate_file <- system.file("testdata", "test.cov", package="famer")
  mask_file <- system.file("testdata", "test.h5", package="famer")
  n_randvecs = 10
  n_blocks = 10
  rand_seed = 123

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
  id <- sprintf("gxg_%d", snp_indices)
  vc_names <- c("id", "grm", "gxg", "error")
  vc_df <- cbind(id, as.data.frame(expected))
  colnames(vc_df) <- vc_names
  vc_df <- pivot_output(vc_df,
                        "component",
                        "vc_estimate",
                        vc_names[2:4])

  # when
  result <- fame(plink_file,
                 pheno_file,
                 covariate_file,
                 mask_file,
                 n_randvecs,
                 snp_indices,
                 n_blocks,
                 rand_seed)
  observed <- result$vc_estimate

  # then
  expect_equal(observed, vc_df, tolerance = 1e-5)
})

test_that("fame end-to-end no covariates but with mask", {
  # given
  plink_file <- gsub("\\.bed", "", system.file("testdata", "test.bed", package="famer"))
  pheno_file <- system.file("testdata", "test_h2_0.5.pheno", package="famer")
  covariate_file = ""
  mask_file <- system.file("testdata", "test.h5", package="famer")
  n_randvecs = 10
  n_blocks = 10
  rand_seed = 123

  snp_indices <- c(3, 8, 9, 13, 16, 19, 29, 34, 93, 97)

  expected_est <- matrix(
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
  expected_se <- matrix(
    c(
      0.1207738, 0.16184765, 0.15113731,
      0.1211554, 0.02243126, 0.08744855,
      0.1209582, 0.04043788, 0.09893784,
      0.1214335, 0.04355752, 0.09343718,
      0.1207895, 0.05948262, 0.10581516,
      0.1201660, 0.08556147, 0.13006928,
      0.1214145, 0.04715497, 0.09027384,
      0.1219379, 0.16492895, 0.17881043,
      0.1210386, 0.03876772, 0.09262924,
      0.1209970, 0.05537885, 0.09436732
    ), ncol = 3, byrow = TRUE
  )
  id <- sprintf("gxg_%d", snp_indices)
  vc_names <- c("id", "grm", "gxg", "error")
  vc_df <- cbind(id, as.data.frame(expected_est))
  colnames(vc_df) <- vc_names
  vc_df <- pivot_output(vc_df,
                        "component",
                        "vc_estimate",
                        vc_names[2:4])
  se_df <- cbind(id, as.data.frame(expected_se))
  colnames(se_df) <- vc_names
  se_df <- pivot_output(se_df,
                        "component",
                        "vc_se",
                        vc_names[2:4])
  # when
  result <- fame(plink_file,
                 pheno_file,
                 covariate_file,
                 mask_file,
                 n_randvecs,
                 snp_indices,
                 n_blocks,
                 rand_seed)
  observed_est <- result$vc_estimate
  observed_se <- result$vc_se

  # then
  expect_equal(observed_est, vc_df, tolerance = 2e-2)
  expect_equal(observed_se, se_df, tolerance = 2e-2)
})
