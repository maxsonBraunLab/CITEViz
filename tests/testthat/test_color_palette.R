#light testing of get_palette function

test_that("returns correct color palette given a number", {
  expect_equal(get_palette(4), 
               c("#00C5CD", "#B68E55", "#5E6E5C", "#CD00CD"))
})
