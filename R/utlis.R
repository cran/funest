gen_tv_list = function(n){
  tv_list = list()
  for(i in 1:n){
    tname = paste("Y", as.character(i), sep = '')
    tv_list[[tname]] = list()
  }
  return(tv_list)
}

cond.prob = function(pred.mod, newdata, Tstart, Tpred){
  #risk.Tstart = as.numeric(summary(survfit(model, newdata = newdata, se.fit = F, conf.int = F), times = Tstart)$surv)
  #risk.Tpred = as.numeric(summary(survfit(model, newdata = newdata, se.fit = F, conf.int = F), times = Tpred)$surv)

  T1 = which(pred.mod$unique.death.times <= Tstart)
  T2 = which(pred.mod$unique.death.times <= Tpred)

  T1 = max(T1)
  T2 = max(T2)

  risk.Tstart = pred.mod$survival[,T1]
  risk.Tpred = pred.mod$survival[,T2]

  return(risk.Tpred/risk.Tstart)
}

gen_arg_list = function(n_of_tv, nofp){
  arg_list = list()
  for(i in 1:n_of_tv){
    arg_list[i] = list(list(type = "uFPCA",nbasis = nofp))
  }
  return(arg_list)
}
