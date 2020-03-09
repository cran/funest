#' Fitting functional ensemble survival tree model
#' @description The function funest_fit takes a long and a short form of the survival data,
#' among other arguments for a random survival forest, to fit an functional ensemble survival tree
#' model for predicting survival probability.
#'
#' @param long_train long form of survival data from the training set
#' @param surv_train short form of survival data from the training set
#' @param noftree number of trees in the random survival forest
#' @param nofcov number of covariates selected in each survival tree
#' @param split_rule binary splitting rule for random survival forest, default is "maxstat"
#' @param tv_names a list of names of time-varying covariates
#' @param fv_names a list of names of fixed covariates
#' @param nofp number of multivariate principal components
#' @param t_star time for the last observed biomarker measurement
#' @param t_pred time at prediction
#' @param ... extra arguments that can be passed to ranger()
#'
#' @return A list compose two items. The first item is a list
#' of necessary information for prediction used in funest_pred()
#' function. The second item is the ranger object of the fitted
#' random survival forest.
#' \itemize{
#' \item misc - a list composed of 1) long_train: long form of survival data from the training set,
#' 2) surv_train: short form of survival data from the training set,
#' 3) fmla: covariates passed into the ensemble survival tree
#' 4) score_names: intermediate names for the covariates
#' 5) nofp: number of multivariate principal components
#' 6) train_data.sub: data frame of all covariates after MFPCA been performed
#' \item rg - functional ensemble survival tree model
#' }
#'
#' @references
#' \insertRef{nestpaper}{funest}
#'
#' \insertRef{ranger}{funest}
#'
#' @importFrom survival Surv
#' @export
#' @examples
#' library(funest)
#' data("long_train")
#' data("surv_train")
#' w = funest_fit(long_train, surv_train, tv_names = list("Y1", "Y2", "Y3"), fv_names = list("W"),
#'  noftree = 10, t_star = 5.5, t_pred = 11)

funest_fit = function(long_train, surv_train,
                    noftree = 500, nofcov = 2, split_rule = "maxstat",
                    tv_names, fv_names,
                    nofp = 3, t_star, t_pred, ...){
  # find number of unique individuals in training sample and test sample
  max_time = max(surv_train$time)
  if(t_star >= t_pred){
    stop("t_star must be smaller than t_pred!")
  }
  if(t_pred > max_time){
    stop("t_pred must be smaller than maximum observation time in the dataset!")
  }
  train_patID = surv_train$id
  train_npat = length(train_patID)

  #
  n_of_tv = length(tv_names)
  train_tv_list = gen_tv_list(n_of_tv)
  train_time_list = list()
  for(i in 1:train_npat){
    tid = train_patID[i]
    tdat = long_train[long_train$id == tid,]
    train_time_list = c(train_time_list, list(tdat$obstime))
    temp_tv = tdat[, unlist(tv_names)]
    for(j in 1:n_of_tv){
      ttv = as.vector(temp_tv[j])[,1]
      train_tv_list[[j]] = c(train_tv_list[[j]], list(ttv))
    }
  }

  train_tfobj = gen_tv_list(n_of_tv)
  for(k in 1:n_of_tv){
    train_tfobj[[k]] = funData::as.funData(funData::irregFunData(
      argvals = train_time_list, X = train_tv_list[[k]]))
  }

  train_mfobj = funData::multiFunData(train_tfobj)

  arg_list = gen_arg_list(n_of_tv, nofp)

  train_mfpca = MFPCA::MFPCA(train_mfobj, M = nofp,
                      uniExpansions = arg_list)

  train_scores = train_mfpca$scores
  train_funs = train_mfpca$vectors

  train_covs = cbind(train_scores, surv_train[,unlist(fv_names)])
  d_delta = ncol(train_covs)

  score_names = c()
  for(q in 1:(d_delta)){
    tname = paste("score", as.character(q), sep = "")
    score_names = c(score_names, tname)
  }
  colnames(train_covs) = score_names

  train_data.sub = cbind(surv_train, train_covs)

  fmla = stats::as.formula(paste("Surv(time,event) ~ ",
                          paste(score_names, collapse= "+")))


  rg = ranger::ranger(fmla, data = train_data.sub, num.trees = noftree,
              mtry = nofcov, splitrule = split_rule, ...)

  misc = list(
    long_train = long_train,
    surv_train = surv_train,
    fmla = fmla,
    score_names = score_names,
    nofp = nofp,
    train_data.sub = train_data.sub
  )

  return(list(misc = misc,
              rg = rg
              ))
}
