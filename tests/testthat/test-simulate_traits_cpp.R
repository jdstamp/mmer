test_that("trait simulation works", {
  # given
  plink_file <-
    gsub("\\.bed",
         "",
         system.file("testdata", "test.bed", package = "famer"))
  additive_heritability <- 0.3
  gxg_heritability <- 0.1
  n_additive_snps <- 50
  group_size <- 10
  gxg_group_1 <- 1:3
  gxg_group_2 <- 7:11

  target_mean <- 0
  target_var <- 1
  target_error_variance <- 1 - additive_heritability - gxg_heritability
  # when
  sim <- simulate_traits_cpp(plink_file,
                             additive_heritability,
                             gxg_heritability,
                             n_additive_snps,
                             gxg_group_1,
                             gxg_group_2)
  observed_mean <- mean(sim$trait)
  observed_var <- var(as.numeric(sim$trait))

  # then
  expect_equal(observed_var, target_var, tolerance = 2e-1)
  expect_equal(sim$epistatic_variance, gxg_heritability, tolerance = 1e-2)
  expect_equal(sim$additive_variance, additive_heritability, tolerance = 1e-2)
  expect_equal(sim$error_variance, target_error_variance, tolerance = 1e-2)
})
