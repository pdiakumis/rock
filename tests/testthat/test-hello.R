context("hello")

test_that("atomics are greeted", {
  expect_equal(length(hello(c("foo", "bar"))), 2)
  expect_equal(length(hello(1:10)), 10)
})

test_that("non-atomics are not greeted", {
  expect_error(hello(mtcars))
  expect_error(hello(list(1, 2, 3)))
})

test_that("can run on internal data", {
  expect_equal(length(hello(cat_names)), 3)
})
