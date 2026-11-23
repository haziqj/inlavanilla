test_that("the version is a date, as in R-INLA", {
  ver <- utils::packageDescription("inlavanilla", fields = "Version")
  # YY.MM.DD, with -1, -2 for a second or third version on one day
  expect_match(ver, "^[0-9]{2}\\.[0-9]{2}\\.[0-9]{2}(-[0-9]+)?$")
  day <- as.Date(sub("-[0-9]+$", "", ver), format = "%y.%m.%d")
  expect_false(is.na(day))
  # A day of slack for time zones
  expect_lte(as.numeric(day), as.numeric(Sys.Date()) + 1)
})
