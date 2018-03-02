context("hello")

test_that("hello atomics", {
  expect_equal(length(hello(c("foo", "bar"))), 2)
  expect_equal(length(hello(1:10)), 10)
})

test_that("hello non-atomics", {
  expect_error(hello(mtcars))
  expect_error(hello(1:10))
})
