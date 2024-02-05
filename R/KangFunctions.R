# Functions provided by Kang et al. (2016):
# L. Kang, A. Liu, and L. Tian. Linear combination methods to improve diagnostic/prognostic accuracy on future observations. 
# Statistical Methods in Medical Research, 25(4):1359–1380, 2016. doi: 10.1177/0962280213481053.
# URL: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4180010/#APP1


## Mann-Whitney U stat for AUC: continuous data

.nonp.auc <- function(u, v) {
  n1 = length(u)
  n2 = length(v)
  return(sum(sapply(u,function(x) sum(x<v)))/n1/n2)
}



## Su and Liu method

.suliu <- function(new.1, new.2) {
  a = stats::var(new.1) + stats::var(new.2)
  b = colMeans(new.2) - colMeans(new.1)
  est.coef = as.numeric(solve(a) %*% b)
  check.sign = .nonp.auc(new.1 %*% est.coef, new.2 %*% est.coef)
  if(check.sign >= 0.5) return(list(coef = est.coef, auc.combined = check.sign)) else return(list(coef = -est.coef, auc.combined = 1 - check.sign))
}



## Logistic regression approach

.logistic <- function(new.1, new.2) {
  n1 = nrow(new.1)
  n2 = nrow(new.2)
  dat.lr <- data.frame(cbind(response = rep(c(0,1), times = c(n1,n2)), rbind(new.1,new.2)))
  obj.lr <- try(glm(response ~., data = dat.lr, family = binomial(link="logit")), silent = T)
  if (any(class(obj.lr) == "try-error")) return(list(coef = rep(0,ncol(new.1)), auc.combined = 0.5)) else est.coef = as.numeric(obj.lr$coef[-1])
  check.sign = .nonp.auc(new.1 %*% est.coef, new.2 %*% est.coef)
  if(check.sign>=0.5) return(list(coef = est.coef, auc.combined = check.sign)) else return(list(coef = -est.coef, auc.combined = 1 - check.sign))
}



####### data.1, data.2 must be of two-column (Modified)

.nonpar.combine2.auc <- function(alpha, rate, data.1, data.2) {
  n1 = nrow(data.1)
  n2 = nrow(data.2)
  new.1 = data.1 %*% c(alpha, rate)
  new.2 = data.2 %*% c(alpha, rate)
  .nonp.auc(new.1, new.2)
}

.nonpar.combine2.coef <- function(new.1, new.2, evalnum = 201) {
  rate = seq(-1, 1, length = evalnum)
  alpha = rev(rate)[-1]
  auc.rate_x = sapply(rate, .nonpar.combine2.auc, alpha = 1, data.1 = new.1, data.2 = new.2)
  auc.alpha_x = sapply(alpha, .nonpar.combine2.auc, rate = 1, data.1 = new.1, data.2 = new.2)
  auc.rateinv_x = sapply(rate, .nonpar.combine2.auc, alpha = -1, data.1 = new.1, data.2 = new.2)
  auc.alphainv_x = sapply(alpha, .nonpar.combine2.auc, rate = -1, data.1 = new.1, data.2 = new.2)
  auc.0 = c(auc.rate_x, auc.alpha_x, auc.rateinv_x, auc.alphainv_x)
  amax.idx = which.max(auc.0)
  if((amax.idx/evalnum) <= 1){
    return(c(alpha = 1, rate = rate[amax.idx], auc.max = auc.0[amax.idx]))
  }else if((amax.idx/evalnum) <= 2){
    return(c(alpha = alpha[amax.idx - evalnum], rate = 1, auc.max = auc.0[amax.idx]))
  }else if((amax.idx/evalnum) <= 3){
    return(c(alpha = -1,rate = rate[amax.idx - 2*evalnum], auc.max = auc.0[amax.idx]))
  }else if((amax.idx/evalnum) <= 4){
    return(c(alpha = alpha[amax.idx - 3*evalnum], rate = -1, auc.max = auc.0[amax.idx]))
  }
}
 
.nonp.auc.check <- function(health, middle) {
  auc.i = numeric(ncol(health))
  for (i in 1:ncol(health)) {
    new.1 = health[,i]
    new.2 = middle[,i]
    auc.i[i] = .nonp.auc(new.1, new.2)}
  auc.i
}



### Step-wise method

.step.coef <- function(new.1, new.2, design = 'step-down') {
  n1 = nrow(new.1)
  n2 = nrow(new.2)
  VARnum = ncol(new.1)
  combcoef = matrix(0, nrow = VARnum - 1, ncol = 2)
  if (design == 'step-down') {
    auc.order = sort(.nonp.auc.check(health = new.1, middle = new.2), index.return = T, decreasing = T)$ix} else {auc.order = sort(.nonp.auc.check(health = new.1, middle = new.2), index.return = T, decreasing = F)$ix}
  combmarker.1 = new.1[,auc.order[1]]
  combmarker.2 = new.2[,auc.order[1]]
  nal.coef = 1
  for (i in 2:VARnum) {
    combmarker.1 = cbind(combmarker.1, new.1[,auc.order[i]])
    combmarker.2 = cbind(combmarker.2, new.2[,auc.order[i]])
    temp.info = .nonpar.combine2.coef(combmarker.1, combmarker.2)
    combcoef[i-1,] = temp.info[1:2]
    nal.coef = c(nal.coef * combcoef[i-1,1], combcoef[i-1,2])
    combmarker.1 = combmarker.1 %*% temp.info[1:2]
    combmarker.2 = combmarker.2 %*% temp.info[1:2]
  }
  nal.coef = nal.coef[sort(auc.order, index.return = T)$ix]
  check.sign = .nonp.auc(new.1 %*% nal.coef, new.2 %*% nal.coef)
  if(check.sign <= 0.5) nal.coef=-nal.coef
  return(list(coef = as.numeric(nal.coef),
              auc.combined = as.numeric(temp.info[3]),
              check = (max(check.sign, 1-check.sign) == temp.info[3])
  ))
}



######### Min-Max method

.liu.coef <- function(data.1, data.2) {
  max_min.1 = cbind(apply(data.1, 1, max), apply(data.1, 1, min))
  max_min.2 = cbind(apply(data.2, 1, max), apply(data.2, 1, min))
  est.coef = .nonpar.combine2.coef(max_min.1, max_min.2)[1:2]
  check.sign = .nonp.auc(max_min.1 %*% est.coef, max_min.2 %*% est.coef)
  if(check.sign >= 0.5) return(list(coef = est.coef, auc.combined = check.sign)) else return(list(coef = -est.coef, auc.combined = 1 - check.sign))
}





### Pepe&Thompson method for multivariate markers (Added)

.pepethompson.multi <- function(X, D){
  init.var <- 1:ncol(X)
  final.var <- NULL; final.beta <- NULL
  comb.table <- utils::combn(init.var, 2)
  AUCs.pepethompson <- sapply(1:ncol(comb.table), function(col){
    output <- .nonpar.combine2.coef(X[D==0,comb.table[,col]], X[D==1,comb.table[,col]])
    biroc.pepethompson <- biROC(X[,comb.table[,col]], D, method = "fixedLinear", coefLinear = output[1:2])
    c(output[1:2], biroc.pepethompson$auc)})
  final.var <- c(final.var, comb.table[,which.max(AUCs.pepethompson[3,])])
  final.beta <- cbind(final.beta, AUCs.pepethompson[1:2,which.max(AUCs.pepethompson[3,])])
  S <- X[,final.var] %*% final.beta
  for(k in 2:(ncol(X)-1)){
    temp.var <- setdiff(init.var, final.var)
    X.temp <- matrix(X[,temp.var], ncol = length(temp.var))
    AUCs.pepethompson <- sapply(1:length(temp.var), function(col){
      output <- .nonpar.combine2.coef(cbind(S[D==0,], X.temp[D==0,col]), cbind(S[D==1,], X.temp[D==1,col]))
      biroc.pepethompson <- biROC(cbind(S, X.temp[,col]), D, method = "fixedLinear", coefLinear = output[1:2])
      c(output[1:2], biroc.pepethompson$auc)})
    final.var <- c(final.var, temp.var[which.max(AUCs.pepethompson[3,])])
    final.beta <- c(AUCs.pepethompson[1,which.max(AUCs.pepethompson[3,])] * final.beta, AUCs.pepethompson[2,which.max(AUCs.pepethompson[3,])])
    S <- X[,final.var] %*% final.beta
  }
  list(final.var = final.var, final.beta = final.beta, S = S)
}