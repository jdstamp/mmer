test_that("fame works", {
  plink_file = "/Users/jds/git/FAME/example/small"
  pheno_file = "/Users/jds/git/FAME/example/h2_0.5.pheno"
  annotation_file = "/Users/jds/git/FAME/example/single-bin/single.annot.small.txt"
  covariate_file = ""

  plink_file <- gsub("\\.bed", "", system.file("testdata", "test.bed", package="famer"))
  annotation_file <- system.file("testdata", "test_single_annot.txt", package="famer")
  pheno_file <- system.file("testdata", "test_h2_0.5.pheno", package="famer")
  print(pheno_file)
  gxg_bin = 0
  num_evec = 10
  snp_index = 10
  jack_number = 10
  fame(plink_file,
       pheno_file,
       annotation_file,
       covariate_file,
       gxg_bin,
       num_evec,
       snp_index,
       jack_number)
})
