# test_that("Compare to MAPIT inference", {
#   # given
#   plink_file <- gsub("\\.bed", "", system.file("testdata", "test.bed", package="famer"))
#   phen_file <- system.file("testdata", "test_h2_0.5.pheno", package="famer")
#   covariate_file = ""
#   mask_file <- ""
#   n_randvecs <- 1000
#   n_blocks <- 10
#   rand_seed <- 123
#   focal_snp <- sample(1:100, 1)
#
#   test_phen <- read.table(phen_file, header = T)
#   test_data <- read_plink(plink_file, verbose = F)
#   y <- test_phen$pheno
#   y <- y - mean(y)
#   X <- t(test_data$X)
#   X_mean <- apply(X, 2, mean)
#   X_sd <- 1 / sqrt(X_mean * (1 - 0.5 * X_mean))
#   # X <- scale(X)
#   X <- t((t(X) - X_mean) * X_sd)
#   apply(X, 2, mean)
#   apply(X, 2, sd)
#
#   n_snps <- ncol(X)
#   n_samples <- nrow(X)
#
#   K <- X %*% t(X) / n_snps
#
#   ERROR <- NULL
#
#   for (i in 1:1) {
#     focal_snp <- sample(1:100, 1)
#     xk <- X[, focal_snp]
#     G <- (n_snps * K - xk %*% t(xk)) / (n_snps - 1)
#     G <- t(t(G * xk) * xk)
#
#
#     S11 <- sum(diag(K %*% K))
#     S21 <- sum(diag(K %*% G))
#     S31 <- sum(diag(K))
#     S22 <- sum(diag(G %*% G))
#     S32 <- sum(diag(G))
#     S33 <- n_samples
#     S_expected <- matrix(c(S11, S21, S31,
#                            S21, S22, S32,
#                            S31, S32, S33), nrow = 3)
#
#     q1 <- t(y) %*% K %*% y
#     q2 <- t(y) %*% G %*% y
#     q3 <- t(y) %*% y
#     q_expected <- c(q1, q2, q3)
#     sigma_expected <- solve(S_expected, q_expected)
#     invS <- solve(S_expected)
#     I <- diag(1, n_samples, n_samples)
#     var_sig_expected <- c()
#     V <- sigma_expected[1] * K + sigma_expected[2] * G + sigma_expected[3] * I
#     for (i in 1:3) {
#       H <- invS[1, i] * K + invS[2, i] * G + invS[3, i] * I
#       var_sig_expected <- c(var_sig_expected, 2 * t(y) %*% H %*% V %*% H %*% y)
#     }
#
#     # when
#     observed <- fame(plink_file,
#                    phen_file,
#                    covariate_file,
#                    mask_file,
#                    focal_snp,
#                    n_randvecs,
#                    n_blocks#,
#                    #rand_seed
#     )
#
#     S_observed <- observed$S
#     q_observed <- observed$q
#     sigma_observed <- observed$vc_estimate$vc_estimate
#     var_sig_observed <- observed$vc_se$vc_se
#
#     error_S <- S_observed - S_expected
#     error_q <- q_observed - q_expected
#     error_sig <- sigma_observed - sigma_expected
#     error_var <- var_sig_observed - var_sig_expected
#
#     error_S_KK <- abs(error_S[1,1])
#     error_S_K <- abs(error_S[3,1])
#     error_S_KG <- abs(error_S[2,1])
#     error_S_G <- abs(error_S[3,2])
#     error_S_GG <- abs(error_S[2,2])
#     error_q_K <- abs(error_q[1])
#     error_q_G <- abs(error_q[2])
#     error_df <- data.frame(
#       element = c("trKK", "trK", "trKG", "trG", "trGG", "yKy", "yGy"),
#       error = c(error_S_KK, error_S_K, error_S_KG, error_S_G, error_S_GG, error_q_K, error_q_G),
#       delta = c(error_S[1,1], error_S[3,1], error_S[2,1], error_S[3,2], error_S[2,2], error_q[1], error_q[2])
#     )
#     ERROR <- rbind(ERROR, error_df)
#   }
#   # then
#   ER_SUMMARY <- ERROR %>%
#     group_by(element) %>%
#     summarize(mean_err = mean(error),
#               sd_err = sd(error),
#               mean_delta = mean(delta),
#               sd_delta = sd(delta)) %>%
#     arrange(desc(mean_err))
#
# })
