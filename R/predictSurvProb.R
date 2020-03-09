#' predictor_loader
#' @description An intermediate function for loading the necessary function into .GlobalEnv
#' @return None
#' @examples
#' # must run the following code before calling funest_pred()
#' predictSurvProb.ranger = predictor_loader()
#' @export
predictor_loader <- function () {
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
  return(predictSurvProb.ranger)
}


predict.cr_survreg <- function(object, newdata = NULL, times, type = 'survival', start = NULL, ...) {
  if (is.null(newdata))
    newdata <- attr(object,'cr_survreg')$data
  outl <- purrr::map(object, ~do.call(predict,
                                      args = c(
                                        list(.x, newdata = newdata, type = type, times = times, start = start),
                                        list(...))
  ))
  as.matrix(as.data.frame(outl))
}

print.cr_survreg <- function(x, ...) {
  purrr::walk2(.x = x, .y = names(x), .f = function(ob,name) {
    cat("\n", name, "\n======\n", sep="")
    print(ob)
  })
}

summary.cr_survreg <- function(object, ...) {
  print(object)
}


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
