/*
 * This file uses the Catch unit testing library, alongside
 * testthat's simple bindings, to test a C++ function.
 *
 * For your own packages, ensure that your test files are
 * placed within the `src/` folder, and that you include
 * `LinkingTo: testthat` within your DESCRIPTION file.
 */

// All test files should include the <testthat.h>
// header file.
#include "read_genotypes.h"
#include <testthat.h>

// Normally this would be a function from your package's
// compiled library -- you might instead just include a header
// file providing the definition, and let R CMD INSTALL
// handle building and linking.
int mapValues(int value) {
  switch (value) {
  case 0:
    return 0;
  case 1:
    return -1;
  case 2:
    return 1;
  case 3:
    return 2;
  default:
    // Handle invalid input
    return -1; // or any other default value you prefer
  }
}

int oldMapping(int val, float p_j, bool allow_missing) {
  if (val == 1 && !allow_missing) {
    val = impute_genotype(p_j);
    val++;
    val = (val == 1) ? 0 : val;
  }
  val--;
  val = (val < 0) ? 0 : val;
  return (val);
}

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.
context("C++ test plink parsing") {

  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  test_that("allow missing true y[l] = 0") {
    // given
    bool allow_missing = true;
    float p_j = 0;
    int yl = 0;
    int expected = oldMapping(yl, p_j, allow_missing);
    // when
    int observed = mapValues(yl);
    // then
    expect_true(expected == observed);
  }

  test_that("allow missing true y[l] = 1") {
    // given
    bool allow_missing = false;
    float p_j = 0;
    int yl = 1;
    int expected = oldMapping(yl, p_j, allow_missing);
    // when
    int val = mapValues(yl);
    std::cout << "@@@@ " << val << std::endl;
    int observed = (val == -1) ? impute_genotype(p_j) : val;
    std::cout << "observed " << observed << std::endl;
    // then
    expect_true(expected == observed);
  }

  test_that("allow missing true y[l] = 2") {
    // given
    bool allow_missing = true;
    float p_j = 0;
    int yl = 2;
    int expected = oldMapping(yl, p_j, allow_missing);
    // when
    int observed = mapValues(yl);
    // then
    expect_true(expected == observed);
  }

  test_that("allow missing true y[l] = 3") {
    // given
    bool allow_missing = true;
    float p_j = 0;
    int yl = 3;
    int expected = oldMapping(yl, p_j, allow_missing);
    // when
    int observed = mapValues(yl);
    // then
    expect_true(expected == observed);
  }
}

context("C++ get_observed_allelefreq") {

  test_that("allow missing true y[l] = 0") {
    // given
    bool allow_missing = true;
    float p_j = 0;
    int yl = 0;
    int expected = oldMapping(yl, p_j, allow_missing);
    // when
    int observed = mapValues(yl);
    // then
    expect_true(expected == observed);
  }

  test_that("allow missing true y[l] = 1") {
    // given
    bool allow_missing = false;
    float p_j = 0;
    int yl = 1;
    int expected = oldMapping(yl, p_j, allow_missing);
    // when
    int val = mapValues(yl);
    std::cout << "@@@@ " << val << std::endl;
    int observed = (val == -1) ? impute_genotype(p_j) : val;
    std::cout << "observed " << observed << std::endl;
    // then
    expect_true(expected == observed);
  }

  test_that("allow missing true y[l] = 2") {
    // given
    bool allow_missing = true;
    float p_j = 0;
    int yl = 2;
    int expected = oldMapping(yl, p_j, allow_missing);
    // when
    int observed = mapValues(yl);
    // then
    expect_true(expected == observed);
  }

  test_that("allow missing true y[l] = 3") {
    // given
    bool allow_missing = true;
    float p_j = 0;
    int yl = 3;
    int expected = oldMapping(yl, p_j, allow_missing);
    // when
    int observed = mapValues(yl);
    // then
    expect_true(expected == observed);
  }
}
