gROC <- function(X, D, ...) {
  UseMethod("gROC")
}

gROC.default <- function(X, D, side = c("right", "left", "both", "both2"),
                         restric = FALSE, optim = TRUE, t0 = NULL, t0max = FALSE, verbose = FALSE, ...){

  levels.names <- levels(as.factor(D))
  controls <- split(X,D)[[levels.names[1]]]; cases <- split(X,D)[[levels.names[2]]]
  D <- ifelse(as.factor(D)==levels.names[1], 0, 1); levels <- levels(as.factor(D))

  m <- length(controls)

  t <- seq(0,1,1/m)
  N <- length(t)

  XX <- sort(c(controls,cases))
  e <- ifelse(length(unique(XX))>1, min(unique(XX)[-1]-unique(XX)[-length(unique(XX))])/2, sqrt(.Machine$double.eps))

  main <- function(side){
    if(side=='right'){
      c <- as.numeric(quantile(controls,1-t,type=3))
      roc <- 1-ecdf(cases)(c)
      results <- list(roc=roc, c=c)
    }
    if(side=='left'){
      c <- as.numeric(quantile(controls,t,type=3)) - e*as.numeric(ecdf(controls)(quantile(controls,t,type=3)) > t)
      roc <- ecdf(cases)(c)
      results <- list(roc=roc, c=c)
    }
    if(side=='both'){
      A <- sapply(1:N, function(i){
        if(i == N){
          roc <- 1; xl <- xu <- max(controls)
        }else{
          gamma <- seq(1,i,1)
          index.gamma.t <- which.max(ecdf(cases)(sort(controls)[gamma]-e) + 1 - ecdf(cases)(sort(controls)[m-i+gamma]))
          gamma.t <- gamma[index.gamma.t]
          xl <-  sort(controls)[gamma.t]; xu <- sort(controls)[m-i+gamma.t]
          roc <- ecdf(cases)(xl-e) + 1 - ecdf(cases)(xu)
        }
        c(roc, xl, xu)
      })
      results <- list(roc=A[1,], xl=A[2,], xu=A[3,])
    }
    if(side=='both2'){
      A <- sapply(1:N, function(i){
        if(i >= m){
          if(i == m){
            xl.opt <- c(min(min(controls)-e, min(cases)), min(controls))
            xu.opt <- c(max(controls), max(max(controls)+e, max(cases)))
            index.opt.t <- which.max(ecdf(cases)(xu.opt-e) - ecdf(cases)(xl.opt))
            xl <- xl.opt[index.opt.t]; xu <- xu.opt[index.opt.t]
          }
          if(i == N){
            xl <- min(min(controls)-e, min(cases)); xu <- max(max(controls)+e, max(cases))
          }
        }else{
          gamma <- seq(1,m-i,1)
          index.gamma.t <- which.max(ecdf(cases)(sort(controls)[gamma+i]-e) - ecdf(cases)(sort(controls)[gamma]))
          gamma.t <- gamma[index.gamma.t]
          xl <-  sort(controls)[gamma.t]; xu <- sort(controls)[gamma.t+i]
        }
        roc <- max(ecdf(cases)(xu-e) - ecdf(cases)(xl), 0)
        c(roc, xl, xu)
      })
      results <- list(roc=A[1,], xl=A[2,], xu=A[3,])
    }
    return(results)
  }

  side <- match.arg(side)
  mainres <- main(side)

  roc <- mainres$roc
  if(side=='right' || side=='left') c <- mainres$c
  if(side=='both' || side=='both2') xl <- mainres$xl; xu <- mainres$xu

  # auc <- mean(roc[-1] + roc[-N])/2
  auc <- sum(roc[-N]*(t[-1] - t[-N]))

  if((side!='both' & side!='both2') & restric==TRUE){
    warning('The non-parametric estimation with restriction is just computed for the generalization (side="both").')
    restric <- FALSE
  }

  if(restric){

    aucfree <- auc

    if(side=='both'){

      gamma.t <- function(i){
        gamma <- seq(1,i,1)
        index.gamma.t <- which.max(ecdf(cases)(sort(controls)[gamma]-e) + 1 - ecdf(cases)(sort(controls)[m-i+gamma]))
        gamma[index.gamma.t]
      }

      auc.rest.i0 <- function(i0){
        TT <- rep(0,N); TT[i0] <- t[i0]
        ROC <- rep(0,N); ROC[i0] <- roc[i0]
        XL <- rep(0,N); XL[i0] <- xl[i0]
        XU <- rep(0,N); XU[i0] <- xu[i0]
        GAMMA <- rep(0,N); GAMMA[i0] <- gamma.t(min(i0,m))

        for(i in seq(max(i0-1,1),1,-1)){
          TT[i] <- (i-1)/m
          if(ROC[i+1]==0 || i==1){
            XL[i] <- min(controls); XU[i] <- max(controls)
          }else{
            gamma <- c(max(GAMMA[i+1]-1,1), GAMMA[i+1])
            index.gamma.t <- which.max(ecdf(cases)(sort(controls)[gamma]-e) + 1 - ecdf(cases)(sort(controls)[m-i+gamma]))
            GAMMA[i] <- gamma[index.gamma.t]
            XL[i] <-  sort(controls)[GAMMA[i]]; XU[i] <- sort(controls)[m-i+GAMMA[i]]
          }
          ROC[i] <- ecdf(cases)(XL[i]-e) + 1 - ecdf(cases)(XU[i])
        }

        for(i in seq(min(i0+1,m+1),m+1,1)){
          TT[i] <- (i-1)/m
          if(ROC[i-1]==1 || i==m+1){
            XL[i] <- XL[i-1]; XU[i] <- XU[i-1]
          }else{
            gamma <- c(GAMMA[i-1], GAMMA[i-1]+1)
            index.gamma.t <- which.max(ecdf(cases)(sort(controls)[gamma]-e) + 1 - ecdf(cases)(sort(controls)[m-i+gamma]))
            GAMMA[i] <- gamma[index.gamma.t]
            XL[i] <-  sort(controls)[GAMMA[i]]; XU[i] <- sort(controls)[m-i+GAMMA[i]]
          }
          ROC[i] <- ecdf(cases)(XL[i]-e) + 1 - ecdf(cases)(XU[i])
        }

        results <- list(roc=ROC, t=TT, xl=XL, xu=XU, auc=sum(ROC[-N]*(TT[-1] - TT[-N])))
        return(results)
      }

    }

    if(side=='both2'){

      gamma.t <- function(i){
        gamma <- seq(1,m-i,1)
        index.gamma.t <- which.max(ecdf(cases)(sort(controls)[gamma+i]-e) - ecdf(cases)(sort(controls)[gamma]))
        gamma.t <- gamma[index.gamma.t]
      }

      auc.rest.i0 <- function(i0){
        TT <- rep(0,N); ROC <- rep(0,N); XL <- rep(0,N); XU <- rep(0,N); GAMMA <- rep(0,N);
        if(i0 >= m) i0 <- m-1
        TT[i0] <- t[i0]; ROC[i0] <- roc[i0]; XL[i0] <- xl[i0]; XU[i0] <- xu[i0]; GAMMA[i0] <- gamma.t(i0)

        for(i in seq(min(i0+1,m-1),m-1,1)){
          TT[i] <- (i-1)/m
          if(ROC[i-1]==1 || i==m-1){
            XL[i] <- xl[m-1]; XU[i] <- xl[m-1]
          }else{
            gamma <- c(max(GAMMA[i-1]-1,1), GAMMA[i-1])
            index.gamma.t <- which.max(ecdf(cases)(sort(controls)[gamma+i]-e) - ecdf(cases)(sort(controls)[gamma]))
            GAMMA[i] <- gamma[index.gamma.t]
            XL[i] <-  sort(controls)[GAMMA[i]]; XU[i] <- sort(controls)[GAMMA[i]+i]
          }
          ROC[i] <- ecdf(cases)(XU[i]-e) - ecdf(cases)(XL[i])
        }

        for(i in seq(max(i0-1,1),1,-1)){
          TT[i] <- (i-1)/m
          gamma <- c(GAMMA[i+1], GAMMA[i+1]+1)
          index.gamma.t <- which.max(ecdf(cases)(sort(controls)[gamma+i]-e) - ecdf(cases)(sort(controls)[gamma]))
          GAMMA[i] <- gamma[index.gamma.t]
          XL[i] <-  sort(controls)[GAMMA[i]]; XU[i] <- sort(controls)[GAMMA[i]+i]
          ROC[i] <- max(ecdf(cases)(XU[i]-e) - ecdf(cases)(XL[i]), 0)
        }

        for(i in (m-1):N){
          TT[i] <- t[i]; ROC[i] <- roc[i]
          XL[i] <- xl[i]; XU[i] <- xu[i]
        }

        results <- list(roc=ROC, t=TT, xl=XL, xu=XU, auc=sum(ROC[-N]*(TT[-1] - TT[-N])))
        return(results)
      }

    }

    if(optim){

      if(side == 'both'){

        X0 <- sort(controls)
        X1 <- sort(cases)
        m <- length(X0); n <- length(X1)

        index.pairs <- combinations(m,2)
        pair.points <- matrix(X0[index.pairs], ncol = 2)
        N <- nrow(pair.points)
        A <- sapply(1:N, function(i){
          xl <- pair.points[i,1]; xu <- pair.points[i,2]
          c(sum(X0 %[]% c(xl,xu)), sum(X1 %[]% c(xl,xu)))
        })
        info <- cbind(1:N, pair.points, A[1,]/m, (n-A[2,])/n, index.pairs[,2] - index.pairs[,1])
        info

        final.point <- which(info[,4] == 1)

        Cmat  <- matrix(Inf, N+1, N+1)
        bar <- txtProgressBar(min = 1, max = N, style = 3)
        for(i in 1:N){
          for(j in 1:N){
            base <- info[j,6]-info[i,6]
            if(sum(pair.points[i,] %[]% pair.points[j,])==2 & base <= 1){
              Cmat[i,j] <- 1/abs(info[j,5]*(info[j,4]-info[i,4]))
            }
          }
          setTxtProgressBar(bar,i)
        }
        close(bar)

        MAX <- max(Cmat[Cmat!=Inf])+1
        Cmat[N+1, 1:N] <- ifelse(info[,6]==1, MAX, Inf)

        diag(Cmat) <- NA
        Cmat

        obj <- allShortestPaths(Cmat)
        output.complete <- list(walk.nodes = extractPath(obj, N+1, final.point))
        output.complete$walk.arcs <- t(sapply(1:(length(output.complete$walk.nodes)-1), function(i){
          node.s <- output.complete$walk.nodes[i];
          node.f <- output.complete$walk.nodes[i+1];
          c(node.s, node.f, Cmat[node.s, node.f])
        }))
        output.complete$walk.arcs

        output.temp <- info[output.complete$walk.nodes[-1],]

        output <- rbind(
          c(0, rep((output.temp[1,2] + output.temp[1,3])/2, 2), 0, 1, 0),
          output.temp,
          c(0, range(cases), 1, 0, 0)
        )
        lengthroc <- nrow(output)

        results <- list(t = rev(1 - output[,4]), roc =  rev(output[,5]), xl = rev(output[,2]), xu = rev(output[,3]), auc = sum(rev(output[,5])[-lengthroc]*(rev(1 - output[,4])[-1] - rev(1 - output[,4])[-lengthroc])))
        # results <- list(roc=ROC, t=TT, xl=XL, xu=XU, auc=sum(ROC[-N]*(TT[-1] - TT[-N])))

      }else{

        warning("Optimal ROC curve under restriction (C) is only computed for side = 'both'. Check your input parameters.")

      }

    }else{

      if(t0max){
        if(verbose){
          cat("Progress bar: Estimation of the optimal initial point of Sp\n"); flush.console()
          bar <- txtProgressBar(min = 0, max = N, style = 3)
        }
        aucsi0 <- sapply(1:N, function(i0){
          if(verbose) setTxtProgressBar(bar, i0)
          auc.rest.i0(i0)$auc
        })
        if(verbose) close(bar)
        i0max <- which.max(aucsi0)
        results <- auc.rest.i0(i0max)
      }else{
        if(is.null(t0)){
          i0 <- which.max(roc - t)
        }else{
          if(t0<1 | t0>N | t0%%1!=0){
            stop("t0 should be an integer number between 1 and control size + 1")
          }else{
            i0 <- t0
          }
        }
        results <- auc.rest.i0(i0)
      }

    }



    roc <- results$roc; t <- results$t
    xl <- results$xl; xu <- results$xu
    auc <- results$auc

  }

  if(side=='right' || side=='left'){
    results <- list(levels = levels.names, controls = controls, cases = cases, side = side, t = t, roc = roc, auc = auc, c = c, param = FALSE)
  }
  if(side=='both' || side=='both2'){
    if(restric){
      if(t0max){
        results <- list(levels = levels.names, controls = controls, cases = cases, side = side, t = t, roc = roc, auc = auc, aucfree = aucfree, aucs = aucsi0, xl = xl, xu = xu, param = FALSE)
      }else{
        results <- list(levels = levels.names, controls = controls, cases = cases, side = side, t = t, roc = roc, auc = auc, aucfree = aucfree, xl = xl, xu = xu, param = FALSE)
      }
    }else{
      results <- list(levels = levels.names, controls = controls, cases = cases, side = side, t = t, roc = roc, auc = auc, xl = xl, xu = xu, param = FALSE)
    }
  }

  attr(results, 'class') <- 'groc'

  return(results)

}


gROC_param <- function(X, D, side = c("right", "left", "both", "both2"), N = NULL, ...){

  D <- as.factor(D); levels <- levels(D)
  controls <- split(X,D)[[levels[1]]]; cases <- split(X,D)[[levels[2]]]
  a <- (mean(cases) - mean(controls))/sd(cases)
  b <- sd(controls)/sd(cases)
  constant <- abs(2*a*b/(1-b^2))
  p0 <- mean(controls) + sd(controls)*a*b/(b^2-1)

  if(is.null(N)){
    t <- seq(0,1,1/1000)
  }else{
    t <- seq(0,1,1/N)
  }
  N <- length(t)

  inf <- mean(controls) + sd(controls)*qnorm(.Machine$double.eps); sup <- mean(controls) + sd(controls)*qnorm(1-.Machine$double.eps)

  main <- function(side){
    if(side=='right'){
      roc <- pnorm(a + b*qnorm(t))
      c <- mean(controls) + sd(controls)*qnorm(1-t)
      c[1] <- sup; c[N] <- inf
      results <- list(roc=roc, c=c)
    }
    if(side=='left'){
      roc <- pnorm(-a + b*qnorm(t))
      c <- mean(controls) + sd(controls)*qnorm(t)
      c[1] <- inf; c[N] <- sup
      results <- list(roc=roc, c=c)
    }
    if(side=='both'){
      gamma <- seq(0,1,0.001)
      roc <- apply(pnorm(a + b*qnorm(gamma%*%t(t))) + 1 - pnorm(a + b*qnorm(1 - (1-gamma)%*%t(t))), 2, max)
      gammat <- apply(pnorm(a + b*qnorm(gamma%*%t(t))) + 1 - pnorm(a + b*qnorm(1 - (1-gamma)%*%t(t))), 2, which.max)
      xl <-  mean(controls) + sd(controls)*qnorm((1-gamma[gammat])*t); xu <- mean(controls) + sd(controls)*qnorm(1-gamma[gammat]*t)
      xl[xl==-Inf] <- inf; xu[xu==Inf] <- sup
      results <- list(roc=roc, xl=xl, xu=xu)
    }
    if(side=='both2'){
      gamma <- seq(0,1,0.001)
      roc <- apply(pnorm(a + b*qnorm(1 - (1-gamma)%*%t(1-t))) - pnorm(a + b*qnorm(gamma%*%(1-t(t)))), 2, max)
      gammat <- apply(pnorm(a + b*qnorm(1 - (1-gamma)%*%t(1-t))) - pnorm(a + b*qnorm(gamma%*%(1-t(t)))), 2, which.max)
      xl <- mean(controls) + sd(controls)*qnorm((1-gamma[gammat])*(1-t)); xu <- mean(controls) + sd(controls)*qnorm(1-gamma[gammat]*(1-t))
      xl[xl==-Inf] <- inf; xu[xu==Inf] <- sup
      results <- list(roc=roc, xl=xl, xu=xu)
    }
    return(results)
  }

  side <- match.arg(side)
  mainres <- main(side)

  auc <- mean(mainres$roc[-1] + mainres$roc[-N])/2
  # auc <- sum(mainres$roc[-N]*(t[-1] - t[-N]))

  if(side=='right' || side=='left'){
    results <- list(levels = levels, controls = controls, cases = cases, side = side, t = t, roc = mainres$roc, auc = auc, c = mainres$c, a = a, b = b, constant = constant, p0 = p0, param = TRUE)
  }
  if(side=='both' || side=='both2'){
    results <- list(levels = levels, controls = controls, cases = cases, side = side, t = t, roc = mainres$roc, auc = auc, xl = mainres$xl, xu = mainres$xu, a = a, b = b, constant = constant, p0 = p0, param = TRUE)
  }

  attr(results, 'class') <- 'groc'

  return(results)

}
