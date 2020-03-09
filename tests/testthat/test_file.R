test_that("loading sample data", {

  data("long_test")
  data("long_train")
  data("surv_train")
  data("surv_test")

})

test_that("fit model", {

  data("long_train")
  data("surv_train")
  w = funest_fit(long_train, surv_train, tv_names = list("Y1", "Y2", "Y3"), noftree = 10,
               fv_names = list("W"), t_star = 5.5, t_pred = 11)

})

