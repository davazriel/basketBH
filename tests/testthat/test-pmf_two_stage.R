test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("PMF sums to one", {

  pmf <- pmf_two_stage(
    n1 = 11,
    r1 = 2,
    n2 = 19,
    p = 0.15
  )

  expect_equal(
    sum(pmf),
    1,
    tolerance = 1e-12
  )

})
