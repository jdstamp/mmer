test_that("mme end-to-end no mask", {
  # given
  plink_file <- gsub("\\.bed", "", system.file("testdata", "test.bed", package="mmer"))
  pheno_file <- system.file("testdata", "test_h2_0.5.pheno", package="mmer")
  mask_file <- ""
  gxg_h5_group <- "gxg"
  ld_h5_group <- "ld"
  chunksize <- 3
  n_randvecs <- 10
  n_blocks <- 10
  rand_seed <- 123
  n_threads <- 3
  log_level <- "DEBUG"

  snp_indices <- c(3, 8, 9, 13, 16, 19, 29, 34, 93, 97)

  expected_est <- matrix(
    c(
      0.422807472, 0.030770829, 0.525349947,
      0.424659762, -0.034297866, 0.586712315,
      0.423531782, -0.057827236, 0.610711311,
      0.424752022, -0.040029493, 0.589124766,
      0.422525234, -0.038636283, 0.590873964,
      0.423146054, 0.022630877, 0.525568865,
      0.424032481, -0.007207079, 0.557439240,
      0.423783848, 0.097367174, 0.448923311,
      0.424160268, -0.014771179, 0.565349071,
      0.423546323, 0.017600500, 0.535181598
    ), ncol = 3, byrow = TRUE
  )
  expected_se <- matrix(
    c(
      0.12076409, 0.07483502, 0.10296534,
      0.12092711, 0.02880959, 0.08871578,
      0.12084017, 0.04131487, 0.10085628,
      0.12135710, 0.04443874, 0.09379178,
      0.12059002, 0.04714957, 0.09737975,
      0.12084237, 0.04549436, 0.09319999,
      0.12096347, 0.05326928, 0.09342799,
      0.12158889, 0.10457918, 0.12324774,
      0.12089678, 0.04467541, 0.09617960,
      0.12099050, 0.07292024, 0.10320580
    ), ncol = 3, byrow = TRUE
  )
  id <- sprintf("rs%d", snp_indices)
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
  result <- mme(plink_file,
                 pheno_file,
                 mask_file,
                 grm_file = "",
                 snp_indices,
                 chunksize,
                 n_randvecs,
                 n_blocks,
                 n_threads,
                 gxg_h5_group,
                 ld_h5_group,
                 rand_seed,
                 log_level)
  observed_est <- result$vc_estimate
  observed_se <- result$vc_se

  # then
  expect_equal(observed_est, vc_df, tolerance = 1e-2)
  expect_equal(observed_se, se_df, tolerance = 1e-1)
})

test_that("mme end-to-end with mask", {
  # given
  plink_file <- gsub("\\.bed", "", system.file("testdata", "test.bed", package="mmer"))
  pheno_file <- system.file("testdata", "test_h2_0.5.pheno", package="mmer")
  mask_file <- system.file("testdata", "test.h5", package="mmer")
  gxg_h5_group <- "gxg"
  ld_h5_group <- "ld"
  chunksize <- 3
  n_randvecs <- 10
  n_blocks <- 10
  rand_seed <- 123
  n_threads <- 1

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
      0.12076409, 0.07483502, 0.10296534,
      0.12092711, 0.02880959, 0.08871578,
      0.12084017, 0.04131487, 0.10085628,
      0.12135710, 0.04443874, 0.09379178,
      0.12059002, 0.04714957, 0.09737975,
      0.12084237, 0.04549436, 0.09319999,
      0.12096347, 0.05326928, 0.09342799,
      0.12158889, 0.10457918, 0.12324774,
      0.12089678, 0.04467541, 0.09617960,
      0.12099050, 0.07292024, 0.10320580
    ), ncol = 3, byrow = TRUE
  )

  id <- sprintf("rs%d", snp_indices)
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
  result <- mme(plink_file,
                 pheno_file,
                 mask_file,
                 grm_file = "",
                 snp_indices,
                 chunksize,
                 n_randvecs,
                 n_blocks,
                 n_threads,
                 gxg_h5_group,
                 ld_h5_group,
                 rand_seed)
  observed_est <- result$vc_estimate
  observed_se <- result$vc_se

  # then
  expect_equal(observed_est, vc_df, tolerance = 1e-1)
  expect_equal(observed_se, se_df, tolerance = 1e-1)
})

test_that("mme end-to-end no mask only one gxg idx", {
  # given
  plink_file <- gsub("\\.bed", "", system.file("testdata", "test.bed", package="mmer"))
  pheno_file <- system.file("testdata", "test_h2_0.5.pheno", package="mmer")
  gxg_h5_group <- "gxg"
  ld_h5_group <- "ld"
  mask_file <- ""
  chunksize <- 3
  n_randvecs <- 10
  n_blocks <- 10
  rand_seed <- 123
  n_threads <- 3
  log_level <- "DEBUG"

  snp_indices <- c(3)

  expected_est <- matrix(
    c(
      0.422807472, 0.030770829, 0.525349947
    ), ncol = 3, byrow = TRUE
  )
  expected_se <- matrix(
    c(
      0.12076409, 0.07483502, 0.10296534
    ), ncol = 3, byrow = TRUE
  )
  id <- sprintf("rs%d", snp_indices)
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
  result <- mme(plink_file,
                pheno_file,
                mask_file,
                grm_file = "",
                snp_indices,
                chunksize,
                n_randvecs,
                n_blocks,
                n_threads,
                gxg_h5_group,
                ld_h5_group,
                rand_seed,
                log_level)
  observed_est <- result$vc_estimate
  observed_se <- result$vc_se

  # then
  expect_equal(observed_est, vc_df, tolerance = 1e-2)
  expect_equal(observed_se, se_df, tolerance = 1e-1)
})

test_that("mme end-to-end no mask - chunksize 1", {
  # given
  plink_file <- gsub("\\.bed", "", system.file("testdata", "test.bed", package="mmer"))
  pheno_file <- system.file("testdata", "test_h2_0.5.pheno", package="mmer")
  mask_file <- ""
  gxg_h5_group <- "gxg"
  ld_h5_group <- "ld"
  chunksize <- 1
  n_randvecs <- 10
  n_blocks <- 10
  rand_seed <- 123
  n_threads <- 3
  log_level <- "DEBUG"

  snp_indices <- c(3, 8, 9, 13, 16, 19, 29, 34, 93, 97)

  expected_est <- matrix(
    c(
      0.422807472, 0.030770829, 0.525349947,
      0.424659762, -0.034297866, 0.586712315,
      0.423531782, -0.057827236, 0.610711311,
      0.424752022, -0.040029493, 0.589124766,
      0.422525234, -0.038636283, 0.590873964,
      0.423146054, 0.022630877, 0.525568865,
      0.424032481, -0.007207079, 0.557439240,
      0.423783848, 0.097367174, 0.448923311,
      0.424160268, -0.014771179, 0.565349071,
      0.423546323, 0.017600500, 0.535181598
    ), ncol = 3, byrow = TRUE
  )
  expected_se <- matrix(
    c(
      0.12076409, 0.07483502, 0.10296534,
      0.12092711, 0.02880959, 0.08871578,
      0.12084017, 0.04131487, 0.10085628,
      0.12135710, 0.04443874, 0.09379178,
      0.12059002, 0.04714957, 0.09737975,
      0.12084237, 0.04549436, 0.09319999,
      0.12096347, 0.05326928, 0.09342799,
      0.12158889, 0.10457918, 0.12324774,
      0.12089678, 0.04467541, 0.09617960,
      0.12099050, 0.07292024, 0.10320580
    ), ncol = 3, byrow = TRUE
  )
  id <- sprintf("rs%d", snp_indices)
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
  result <- mme(plink_file,
                pheno_file,
                mask_file,
                grm_file = "",
                snp_indices,
                chunksize,
                n_randvecs,
                n_blocks,
                n_threads,
                gxg_h5_group,
                ld_h5_group,
                rand_seed,
                log_level)
  observed_est <- result$vc_estimate
  observed_se <- result$vc_se

  # then
  expect_equal(observed_est, vc_df, tolerance = 1e-2)
  expect_equal(observed_se, se_df, tolerance = 1e-1)
})

test_that("mme end-to-end but with mask - chunksize 1", {
  # given
  plink_file <- gsub("\\.bed", "", system.file("testdata", "test.bed", package="mmer"))
  pheno_file <- system.file("testdata", "test_h2_0.5.pheno", package="mmer")
  mask_file <- system.file("testdata", "test.h5", package="mmer")
  gxg_h5_group <- "gxg"
  ld_h5_group <- "ld"
  chunksize <- 1
  n_randvecs <- 10
  n_blocks <- 10
  rand_seed <- 123
  n_threads <- 1

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
      0.12076409, 0.07483502, 0.10296534,
      0.12092711, 0.02880959, 0.08871578,
      0.12084017, 0.04131487, 0.10085628,
      0.12135710, 0.04443874, 0.09379178,
      0.12059002, 0.04714957, 0.09737975,
      0.12084237, 0.04549436, 0.09319999,
      0.12096347, 0.05326928, 0.09342799,
      0.12158889, 0.10457918, 0.12324774,
      0.12089678, 0.04467541, 0.09617960,
      0.12099050, 0.07292024, 0.10320580
    ), ncol = 3, byrow = TRUE
  )

  id <- sprintf("rs%d", snp_indices)
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
  result <- mme(plink_file,
                pheno_file,
                mask_file,
                grm_file = "",
                snp_indices,
                chunksize,
                n_randvecs,
                n_blocks,
                n_threads,
                gxg_h5_group,
                ld_h5_group,
                rand_seed)
  observed_est <- result$vc_estimate
  observed_se <- result$vc_se

  # then
  expect_equal(observed_est, vc_df, tolerance = 1e-1)
  expect_equal(observed_se, se_df, tolerance = 1e-1)
})

test_that("mme end-to-end no mask but grm_file given", {
  # given
  plink_file <- gsub("\\.bed", "", system.file("testdata", "test.bed", package="mmer"))
  pheno_file <- system.file("testdata", "test_h2_0.5.pheno", package="mmer")
  mask_file <- ""
  grm_file <- gsub("\\.bed", "", system.file("testdata", "test_2.bed", package="mmer"))
  gxg_h5_group <- "gxg"
  ld_h5_group <- "ld"
  chunksize <- 3
  n_randvecs <- 10
  n_blocks <- 10
  rand_seed <- 123
  n_threads <- 1
  log_level <- "DEBUG"

  snp_indices <- c(3, 8, 9, 13, 16, 19, 29, 34, 93, 97)

  expected_est <- matrix(
    c( 0.312522861,  0.046777132,  0.622023339,
    0.309669317, -0.026756092, 0.693556216,
    0.308059043, -0.052450517,  0.720612821,
    0.310735042,  -0.032323038,  0.695527045,
    0.307407862, -0.035542460,  0.702862358,
     0.312659400,  0.028673626,  0.629217578,
     0.311631355,  0.003075571, 0.659642252,
     0.317221962,  0.099933139,  0.552807241,
     0.310994395, -0.005147438,  0.668621299,
     0.311537460,  0.027011965, 0.638850872
    ), ncol = 3, byrow = TRUE
  )
  expected_se <- matrix(
    c(
        0.13025910, 0.07776531, 0.13485656, 0.12899953, 0.02963136, 0.12461910,
        0.12921256, 0.04123380, 0.13290286, 0.12947492, 0.04419796, 0.12994167,
        0.12844011, 0.04655299, 0.12859591, 0.12903855, 0.04598168, 0.12000908,
        0.12903300, 0.05619646, 0.12648850, 0.13037978, 0.10439670, 0.14663854,
        0.12972055, 0.04413555, 0.13109683, 0.12965935, 0.07485820, 0.13425578
    ), ncol = 3, byrow = TRUE
  )
  id <- sprintf("rs%d", snp_indices)
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
  result <- mme(plink_file,
                pheno_file,
                mask_file,
                grm_file = grm_file,
                snp_indices,
                chunksize,
                n_randvecs,
                n_blocks,
                n_threads,
                gxg_h5_group,
                ld_h5_group,
                rand_seed,
                log_level)
  observed_est <- result$vc_estimate
  observed_se <- result$vc_se

  # then
  expect_equal(observed_est, vc_df, tolerance = 1e-2)
  expect_equal(observed_se, se_df, tolerance = 1e-1)
})
