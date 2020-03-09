#' Predicting survival probability with time-varing covariates
#' @description The function funest_pred takes the functional ensemble survival tree object
#' from funest_fit() to produce predicted survival probability at user specified t_star
#'  and t_pred along with prediction accuracy measures.
#'  Must run "predictSurvProb.ranger = predictor_loader()" before calling this function.
#'
#' @param funest.fit returned object from funest_fit() function
#' @param long_test long form of survival data from the testing set
#' @param surv_test short form of survival data from the testing set
#' @param tv_names a list of names of time-varying covariates
#' @param fv_names a list of names of fixed covariates
#' @param t_star time for the last observed biomarker measurement
#' @param t_pred time at prediction
#'
#' @return A list of three items. The first is a matrix of individual ID and
#' their corresponding predicted survival probability. The second is
#' the estimated Brier score. The third is the estimated area under the ROC curve.
#' \itemize{
#' \item pred_pb - predicted survival probability at t_pred for each individual
#' conditional on being alive at t_star
#' \item bs - Brier score
#' \item AUC - area under the receiver operating characteristic (ROC) curve
#' }
#'
#' @references
#' \insertRef{auc}{funest}
#'
#' \insertRef{bs}{funest}
#'
#' @importFrom survival Surv
#' @importFrom prodlim Hist
#' @importFrom pec pec
#' @export
#' @examples
#' library(funest)
#' data("long_train")
#' data("surv_train")
#' data("long_test")
#' data("surv_test")
#' # must run the following line before calling funest_pred()
#' predictSurvProb.ranger = predictor_loader()
#' w = funest_fit(long_train, surv_train, tv_names = list("Y1", "Y2", "Y3"),noftree = 10,
#'  fv_names = list("W"), t_star = 5.5, t_pred = 11)
#' pred = funest_pred(w, long_test, surv_test, tv_names = list("Y1", "Y2", "Y3"),
#'  fv_names = list("W"), t_star = 5.5, t_pred = 11)
#' pred$bs
#' pred$AUC


funest_pred = function(funest.fit, long_test, surv_test,
                     tv_names, fv_names,
                     t_star, t_pred){
  # combine data together
  max_time = max(surv_test$time)
  if(t_star >= t_pred){
    stop("t_star must be smaller than t_pred!")
  }
  if(t_pred > max_time){
    stop("t_pred must be smaller than maximum observation time in the dataset!")
  }
  model = funest.fit[[1]]
  rg = funest.fit[[2]]

  long_train = model[[1]]
  surv_train = model[[2]]
  fmla = model[[3]]
  score_names = model[[4]]
  nofp = model[[5]]
  train_data.sub = model[[6]]

  train_patID = surv_train$id
  train_npat = length(train_patID)

  n_of_tv = length(tv_names)
  long_combine = rbind(long_train, long_test)
  combine_patID = unique(long_combine$id)
  combine_npat = length(combine_patID)
  combine_tv_list = gen_tv_list(n_of_tv)
  combine_time_list = list()
  for(i in 1:combine_npat){
    tid = combine_patID[i]
    tdat = long_combine[long_combine$id == tid,]
    combine_time_list = c(combine_time_list, list(tdat$obstime))
    temp_tv = tdat[, unlist(tv_names)]
    for(j in 1:n_of_tv){
      ttv = as.vector(temp_tv[j])[,1]
      combine_tv_list[[j]] = c(combine_tv_list[[j]], list(ttv))
    }
  }

  combine_tfobj = gen_tv_list(n_of_tv)
  for(k in 1:n_of_tv){
    combine_tfobj[[k]] = funData::as.funData(funData::irregFunData(
      argvals = combine_time_list, X = combine_tv_list[[k]]))
  }

  combine_mfobj = funData::multiFunData(combine_tfobj)
  arg_list = gen_arg_list(n_of_tv, nofp)
  combine_mfpca = MFPCA::MFPCA(combine_mfobj, M = nofp,
                        uniExpansions = arg_list)

  combine_scores = combine_mfpca$scores

  test_scores = combine_scores[-c(1:train_npat),]

  test_covs = cbind(test_scores, surv_test[, unlist(fv_names)])
  colnames(test_covs) = score_names

  test_data.sub = cbind(surv_test, test_covs)

  pred.mod = predict(rg, data = test_data.sub, importance = "none")

  predictSurvProb.ranger <- function (object, newdata, times, ...) {
    #ptemp <- ranger:::predict.ranger(object, data = newdata, importance = "none")$survival
    ptemp <- predict(object, data = newdata, importance = "none")$survival
    pos <- prodlim::sindex(jump.times = object$unique.death.times,
                           eval.times = times)
    p <- cbind(1, ptemp)[, pos + 1, drop = FALSE]
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
      stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",
                 NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ",
                 NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
    p
  }

  dummy = predictSurvProb.ranger(rg, test_data.sub, test_data.sub$time)

  SB = pec(object = rg, formula = fmla,
           traindata = train_data.sub,
           data=test_data.sub)

  upper = which(SB$time > t_pred)[1]
  lower = which(SB$time <= t_pred)[length(which(SB$time <= t_pred))]

  btw = c(lower, upper)

  #this is the closest SB to t_pred
  take_this = which(abs(c(SB$time[lower], SB$time[upper]) - t_pred) == min(abs(c(SB$time[lower], SB$time[upper]) - t_pred)))

  surv_time = pred.mod$unique.death.times
  sb_time = SB$time[btw[take_this]]
  surv_ind = max(which(surv_time <= sb_time))

  pred_pb = pred.mod$survival[,surv_ind]

  surv_pb = data.frame(ID = surv_test$id,
                       pred_pb = pred_pb)

  rf.sb = SB$AppErr$ranger[btw[take_this]]

  timeEvent = test_data.sub[, c("time", "event")]
  #pred.mod = predict(rg, data = test_data.sub, importance = "none")


  roc = tdROC::tdROC(
    X = 1- cond.prob(pred.mod = pred.mod, newdata = test_data.sub, Tstart = t_star, Tpred = t_pred),
    Y= test_data.sub$time,
    delta = test_data.sub$event,
    tau = t_pred, span = 0.05,
    nboot = 0, alpha = 0.05,
    n.grid = 1000, cut.off = 0.5)

  auc.m = roc$AUC$value

  return(list(pred_pb = surv_pb, bs = rf.sb, AUC = auc.m))
}
