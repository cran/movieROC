print.groc <- function(x, ...){

  obj <- x
  cat("Data was encoded with", obj$levels[1], "(controls) and", obj$levels[2], "(cases).\n")

  printside <- function(side){
    switch(side,
           right = cat("It is assumed that larger values of the marker indicate larger confidence that a given subject is a case.\n"),
           left = cat("It is assumed that lower values of the marker indicate larger confidence that a given subject is a case.\n"),
           both = cat("It is assumed that both lower and larges values of the marker indicate larger confidence that a given subject is a case.\n"),
           both2 = cat("It is assumed that both lower and larges values of the marker indicate larger confidence that a given subject is a control.\n")
    )
  }
  printside(obj$side)

  cat("There are", length(obj$controls),"controls and", length(obj$cases),"cases.\n")

  index.youden <- which.max(obj$roc - obj$t)
  subset <- format(predict(obj, FPR=obj$t[index.youden])$ClassSubsets, digits = 3, trim = TRUE)
  sentence.temp <- paste("The specificity and sensitivity reported by the Youden index are ", format(1-obj$t[index.youden], digits = 3, nsmall = 2), " and ", format(obj$roc[index.youden], digits = 3, nsmall = 2), ", respectively, corresponding to the following classification subset: ", sep = "")
  if(obj$side!="both"){
    cat(sentence.temp, "(", subset[1], ", ", subset[2], ").\n", sep="")
  }else if(subset[1,1] == subset[1,2] | subset[2,1] == subset[2,2]){
    if(subset[1,1] == subset[1,2]) subset <- subset[-1,] else subset <- subset[-2,]
    cat(sentence.temp, "(", subset[1], ", ", subset[2], ").\n", sep="")
  }else{
    cat(sentence.temp, "(", subset[1,1], ", ", subset[1,2], ") U (", subset[2,1], ", ", subset[2,2], ").\n", sep="")
  }

  printauc <- function(side){
    switch(side,
           right = cat("The area under the right-sided ROC curve (AUC) is ", format(obj$auc, digits = 3, nsmall = 3),".\n", sep=""),
           left = cat("The area under the left-sided ROC curve (AUC) is ", format(obj$auc, digits = 3, nsmall = 3),".\n", sep=""),
           both = cat("The area under the gROC curve (gAUC) is ", format(obj$auc, digits = 3, nsmall = 3),".\n", sep=""),
           both2 = cat("The area under the opposite gROC curve (gAUC) is ", format(obj$auc, digits = 3, nsmall = 3),".\n", sep="")
    )
  }
  printauc(obj$side)

}


print.hroc <- function(x, ...){

  obj <- x
  cat("Data was encoded with", obj$levels[1], "(controls) and", obj$levels[2], "(cases).\n")

  cat("There are", sum(obj$D==0),"controls and", sum(obj$D==1),"cases.\n")

  printtype <- function(type){
    switch(type,
           lrm = cat("A logistic regression model of the form", obj$formula, "has been performed.\n"),
           overfitting = cat("The overfitted ROC curve is reported.\n"),
           standard = cat("The standard ROC curve is reported.\n")
    )
  }
  printtype(obj$type)

  if(obj$type == 'lrm'){
    cat("The estimated parameters of the model are the following:\n")
    print(format(obj$model, digits = 4, trim = TRUE))
  }

  index.youden <- which.max(obj$Sp + obj$Se)
  subset <- predict(obj, FPR=1-obj$Sp[index.youden])$ClassSubsets
  if(length(subset) == 0 & index.youden == 1) subset <- predict(obj, FPR = 1e-16)$ClassSubsets
  printsubset <- paste("(", format(subset[1,1], digits = 3),", ",format(subset[1,2], digits = 3),")", sep="")
  if(nrow(subset) > 1){
    for(i in 2:nrow(subset)){
      printsubset <- paste(printsubset, paste("U (",format(subset[i,1], digits = 3),", ",format(subset[i,2], digits = 3),")", sep=""))
    }
  }
  cat("The specificity and sensitivity reported by the Youden index are ", format(obj$Sp[index.youden], digits = 3, nsmall = 2), " and ", format(obj$Se[index.youden], digits = 3, nsmall = 2), ", respectively, corresponding to the following classification subset: ", printsubset, sep="")

  cat(".\nThe area under the ROC curve (AUC) is ", format(obj$auc, digits = 3, nsmall = 3),".\n", sep="")

}



print.biroc <- function(x, ...){

  obj <- x
  cat("Data was encoded with", obj$levels[1], "(controls) and", obj$levels[2], "(cases).\n")

  cat("There are", sum(obj$D==0),"controls and", sum(obj$D==1),"cases.\n")

  index.youden <- which.max(obj$roc - obj$t)

  if(obj$method == "lrm"){

    if(obj$stepModel){
      cat("A stepwise logistic regression model from the initial formula", obj$formula, "has been performed.\n")
      cat("The estimated parameters of the resulting model are the following:\n")
      print(obj$coefModel)
    }else{
      cat("A logistic regression model of the form", obj$formula, "has been performed.\n")
      cat("The estimated parameters of the model are the following:\n")
      print(obj$coefModel)
    }

    sent <- " for the transformation h(X.1,X.2) = ";
    for(i in 1:length(obj$coefModel)){
      if(i == 1){
        sent <- paste0(sent, format(obj$coefModel[i], digits = 3))
      }else{
        sent <- paste0(sent, ifelse(obj$coefModel[i]>0, " + ", " - "), format(abs(obj$coefModel[i]), digits = 3), "*", names(obj$coefModel)[i])
      }
    }

  }else if(obj$method == "fixedLinear"){

    if(obj$methodLinear == "coefLinear"){
      cat("A linear combination with fixed parameters in `coefLinear` has been considered.\n")
    }else{
      cat("A linear combination with fixed parameters estimated by", obj$methodLinear,"approach has been considered.\n")
    }

    sent <- " for the transformation h(X.1,X.2) = ";
    part <- c("X.1", "X.2")
    for(i in 1:ncol(obj$X)){
      if(i == 1){
        sent <- paste0(sent, ifelse(obj$coefLinear[i]>0, " ", " - "), format(abs(obj$coefLinear[i]), digits = 3), "*", part[i])
      }else{
        sent <- paste0(sent, ifelse(obj$coefLinear[i]>0, " + ", " - "), format(abs(obj$coefLinear[i]), digits = 3), "*", part[i])
      }
    }

    }else if(obj$method == "fixedQuadratic"){

    sent <- " for the transformation h(X.1,X.2) = ";
    part <- c("X.1", "X.2", "X.1*X.2", "X.1^2", "X.2^2")
    for(i in 1:length(obj$coefQuadratic)){
      if(i == 1){
        sent <- paste0(sent, ifelse(obj$coefQuadratic[i]>0, " ", " - "), format(abs(obj$coefQuadratic[i]), digits = 3), "*", part[i])
      }else{
        sent <- paste0(sent, ifelse(obj$coefQuadratic[i]>0, " + ", " - "), format(abs(obj$coefQuadratic[i]), digits = 3), "*", part[i])
      }
    }

    }else if(obj$method == "dynamicMeisner" | obj$method == "dynamicEmpirical"){

      cat("A linear combination with dynamic parameters has been considered.\n")

      coefs <- obj$coefLinear[,index.youden]
      sent <- " for the transformation h(X.1,X.2) = ";
      part <- c("X.1", "X.2")
      for(i in 1:ncol(obj$X)){
        if(i == 1){
          sent <- paste0(sent, ifelse(coefs[i]>0, " ", " - "), format(abs(coefs[i]), digits = 3), "*", part[i])
        }else{
          sent <- paste0(sent, ifelse(coefs[i]>0, " + ", " - "), format(abs(coefs[i]), digits = 3), "*", part[i])
        }
      }

    }else{ # obj$method == "kernelOptimal"
      
      cat("Optimal transformation based on bivariate kernel density estimation has been considered.\n")
      
      cat(obj$H.method,"method has been used to estimate the bandwidth matrix in each group.\n")
      
      sent <- " for the transformation computed in optimalT(X)"
      
    }
  
  cat("The specificity and sensitivity reported by the Youden index are ", format(1-obj$t[index.youden], digits = 3, nsmall = 2), " and ", format(obj$roc[index.youden], digits = 3, nsmall = 2), ", respectively, corresponding to the cut-off point ", format(obj$c[index.youden], digits=3), sent, ".\n", sep="")

  cat("The area under the ROC curve (AUC) is ", format(obj$auc, digits = 3, nsmall = 3),".\n", sep="")

}



print.multiroc <- function(x, ...){

  obj <- x
  cat("Data was encoded with", obj$levels[1], "(controls) and", obj$levels[2], "(cases).\n")

  cat("There are", sum(obj$D==0),"controls and", sum(obj$D==1),"cases.\n")

  cat("A total of", ncol(obj$X), "variables have been considered.\n")


  index.youden <- which.max(obj$roc - obj$t)

  if(obj$method == "lrm"){

    if(obj$stepModel){
      cat("A stepwise logistic regression model from the initial formula", obj$formula, "has been performed.\n")
      cat("The estimated parameters of the resulting model are the following:\n")
      print(obj$coefModel)
    }else{
      cat("A logistic regression model of the form", obj$formula, "has been performed.\n")
      cat("The estimated parameters of the model are the following:\n")
      print(obj$coefModel)
    }
    
    sent <- " for the transformation h(X) in the formula above"

  }else if(obj$method == "fixedLinear"){

    if(obj$methodLinear == "coefLinear"){
      cat("A linear combination with fixed parameters in `coefLinear` has been considered.\n")
    }else{
      cat("A linear combination with fixed parameters estimated by", obj$methodLinear,"approach has been considered.\n")
    }
    
    sent <- " for the transformation h(X) = ";
    if(length(colnames(obj$X))!=ncol(obj$X)) part <- paste0("X.", 1:ncol(obj$X)) else part <- colnames(obj$X)
    for(i in 1:ncol(obj$X)){
      if(i == 1){
        sent <- paste0(sent, ifelse(obj$coefLinear[i]>0, " ", " - "), format(abs(obj$coefLinear[i]), digits = 3), "*", part[i])
      }else{
        sent <- paste0(sent, ifelse(obj$coefLinear[i]>0, " + ", " - "), format(abs(obj$coefLinear[i]), digits = 3), "*", part[i])
      }
    }

  }else if(obj$method == "dynamicMeisner"){

    cat("A linear combination with dynamic parameters has been considered.\n")

    coefs <- obj$coefLinear[,index.youden]
    sent <- " for the transformation h(X) = ";
    if(length(colnames(obj$X))!=ncol(obj$X)) part <- paste0("X.", 1:ncol(obj$X)) else part <- colnames(obj$X)
    for(i in 1:ncol(obj$X)){
      if(i == 1){
        sent <- paste0(sent, ifelse(coefs[i]>0, " ", " - "), format(abs(coefs[i]), digits = 3), "*", part[i])
      }else{
        sent <- paste0(sent, ifelse(coefs[i]>0, " + ", " - "), format(abs(coefs[i]), digits = 3), "*", part[i])
      }
    }

  }else{ # obj$method == "kernelOptimal"
    
    cat("Optimal transformation based on bivariate kernel density estimation has been considered.\n")
    
    cat(obj$H.method,"method has been used to estimate the bandwidth matrix in each group.\n")
    
    sent <- " for the transformation computed in optimalT(X)"
    
  }
  
  cat("The specificity and sensitivity reported by the Youden index are ", format(1-obj$t[index.youden], digits = 3, nsmall = 2), " and ", format(obj$roc[index.youden], digits = 3, nsmall = 2), ", respectively, corresponding to the cut-off point ", format(obj$c[index.youden], digits=3), sent, ".\n", sep="")
  
  cat("The area under the ROC curve (AUC) is ", format(obj$auc, digits = 3, nsmall = 3),".\n", sep="")

}
