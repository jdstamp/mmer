test_that("simulate_a_traits works", {
  # given
  plink_file <-
    gsub("\\.bed", "", system.file("testdata", "test.bed", package = "famer"))
  out_file <- tempfile()
  additive_heritability <- 0.3
  gxg_heritability <- 0.1
  additive_snps <- sort(sample(1:100, 50, replace = F))
  gxg_group_1 <- sort(sample(additive_snps, 10, replace = F))
  gxg_group_2 <-sort(sample(setdiff(additive_snps, gxg_group_1), 10, replace = F))
  target_mean <- 0
  target_var <- 1
  column_names <- c("FID", "IID", "TRAIT")
  n_samples <- 200
  # when
  simulate_a_trait(plink_file,
                         out_file,
                         additive_heritability,
                         gxg_heritability,
                         additive_snps,
                         gxg_group_1,
                         gxg_group_2)
  observed <- read.table(out_file, header = T)
  # then
  expect_equal(colnames(observed), column_names)
  expect_equal(nrow(observed), n_samples)
})

test_that("simulate_a_traits works for zero gxg heritability", {
  # given
  plink_file <-
    gsub("\\.bed", "", system.file("testdata", "test.bed", package = "famer"))
  out_file <- tempfile()
  additive_heritability <- 0.3
  gxg_heritability <- 0.0
  additive_snps <- sort(sample(1:100, 50, replace = F))
  gxg_group_1 <- sort(sample(additive_snps, 10, replace = F))
  gxg_group_2 <-sort(sample(setdiff(additive_snps, gxg_group_1), 10, replace = F))
  target_mean <- 0
  target_var <- 1
  column_names <- c("FID", "IID", "TRAIT")
  n_samples <- 200
  # when
  simulate_a_trait(plink_file,
                   out_file,
                   additive_heritability,
                   gxg_heritability,
                   additive_snps,
                   gxg_group_1,
                   gxg_group_2)
  observed <- read.table(out_file, header = T)
  # then
  expect_equal(mean(observed$TRAIT), target_mean)
  expect_equal(var(observed$TRAIT), target_var, tolerance = 1e-1)
})


test_that("simulate_a_traits works on ukbb data", {
  # given
  plink_file <- "/Users/jds/data/ukbb/ukbb_c12"
  out_file <- tempfile()
  additive_heritability <- 0.4
  gxg_heritability <- 0.0
  additive_snps <- sort(sample(1:100, 50, replace = F))
  gxg_group_1 <- sort(sample(additive_snps, 10, replace = F))
  gxg_group_2 <-sort(sample(setdiff(additive_snps, gxg_group_1), 10, replace = F))
  target_mean <- 0
  target_var <- 1
  column_names <- c("FID", "IID", "TRAIT")
  n_samples <- 200
  # when
  simulate_a_trait(plink_file,
                   out_file,
                   additive_heritability,
                   gxg_heritability,
                   additive_snps,
                   gxg_group_1,
                   gxg_group_2)
  observed <- read.table(out_file, header = T)
  print("yesss")
  print(nrow(observed))
  # then
  expect_equal(mean(observed$TRAIT), target_mean)
  expect_equal(var(observed$TRAIT), target_var, tolerance = 1e-1)

  result <- fame(plink_file,
                 out_file,
                 "",
                 "",
                 c(1),
                 10,
                 370)
})
