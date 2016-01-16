#' Summarizing Penalized Multinomial Logistic Model Fits
#' @description This function is for class \code{pmlr} object.
#'
#' @param object an object of class \code{"pmlr"}, a result of a call to \code{\link{pmlr}}.
#' @param ... further arguments passed to or from other methods
#' @return \code{summary.pmlr} returns an object of class \code{"summary.pmlr"}, a list with components
#' \item{call}{the matched call of the \code{object}}
#' \item{method}{which method was used for hypothesis testing and computation of confidence intervals}
#' \item{coef}{an array containing the coefficient estimates, standard errors, and
#' test statistics and their p-values asssociated with the chosen method for the p parameters for the J categories}
#' \item{joint.test}{an array contatining the test statistics and p-values from constrained hypothesis tests under all betas = 0,
#' all betas are equal, and betas are proportional.}
#' \item{test.all0.vs.constraint}{Returned only if joint hypothesis testing was done:
#' An array containing likelihood ratio test statistics and p-values testing all \eqn{H_0}: betas=0 vs. other constraints (\eqn{H_C}),
#' which are 'all betas are equal' and 'betas are proportion'.}
#' @export
summary.pmlr <- function(object, ...) {
  ret=list()
  ret$call <- object$call
  ret$method <- object$method

  #print(object$call)

  p <- dim(object$coefficients)[2]
  J <- dim(object$coefficients)[3]

  warnSep <- FALSE;

  if(object$method == "wald") {
    cat("\n Wald Confidence Intervals and P-values \n\n")
    output1 <- array(data = NA, dim = c(p,9,J))
    dimnames(output1) <- list(dimnames(object$coefficients)[[2]],
                              c("Estimate","Wald Std. Err.","Lower CI","Upper CI","OR","Lower CI (OR)","Upper CI (OR)","ChiSq","Pr(>ChiSq)"),
                              dimnames(object$coefficients)[[3]])
    for (i in 1:J) {
      output1[,1,i] <- ifelse(is.infinite(object$separation[,,i]), object$separation[,,i], object$coefficients[,,i])
      output1[,2,i] <- ifelse(is.infinite(object$separation[,,i]), NA, sqrt(diag(object$var[,,i])))
      output1[,5,i] <- ifelse(is.infinite(object$separation[,,i]), exp(object$separation[,,i]), exp(object$coefficients[,,i]))
      output1[,8,i] <- ifelse(is.infinite(object$separation[,,i]), NA, object$statistic[,,i])
      output1[,9,i] <- ifelse(is.infinite(object$separation[,,i]), NA, object$pvalue[,,i])

      if(object$CI.calculate){
        output1[,3,i] <- ifelse(is.infinite(object$separation[,,i]), NA, object$CI.lower[,,i])
        output1[,4,i] <- ifelse(is.infinite(object$separation[,,i]), NA, object$CI.upper[,,i])
        output1[,6,i] <- ifelse(is.infinite(object$separation[,,i]), NA, exp(object$CI.lower[,,i]))
        output1[,7,i] <- ifelse(is.infinite(object$separation[,,i]), NA, exp(object$CI.upper[,,i]))
      }

      if (is.infinite(as.vector(object$separation)[i]))
      {
        warnSep <- TRUE;
      }
    }
  }

  if (object$method == "likelihood")
  {
    #cat("\n Likelihood Confidence Intervals and P-values \n\n ")
    output1 <- array(data = NA, dim = c(p,9,J))
    dimnames(output1) <- list(dimnames(object$coefficients)[[2]],
                              c("Estimate","Wald Std. Err.","Lower CI","Upper CI","OR Estimate","Lower CI (OR)","Upper CI (OR)","ChiSq","Pr(>ChiSq)"),
                              dimnames(object$coefficients)[[3]])
    for (i in 1:J) {
      output1[,1,i] <- ifelse(is.infinite(object$separation[,,i]), object$separation[,,i], object$coefficients[,,i])
      output1[,2,i] <- ifelse(is.infinite(object$separation[,,i]), NA, sqrt(diag(object$var[,,i])))
      output1[,5,i] <- ifelse(is.infinite(object$separation[,,i]), exp(object$separation[,,i]), exp(object$coefficients[,,i]))
      output1[,8,i] <- object$statistic[,,i]
      output1[,9,i] <- object$pvalue[,,i]

      if(object$CI.calculate){
        output1[,3,i] <- object$CI.lower[,,i]
        output1[,4,i] <- object$CI.upper[,,i]
        output1[,6,i] <- exp(object$CI.lower[,,i])
        output1[,7,i] <- exp(object$CI.upper[,,i])
      }

      if (is.infinite(as.vector(object$separation)[i]))
      {
        warnSep <- TRUE;
      }
    }
  }

  if (object$method == "score") ## NOT WORKING AT THE MOMENT
  {
    #cat("\n Score Hypothesis Tests and P-values \n\n ")
    output1 <- array(data = NA, dim = c(p,4,J))
    dimnames(output1) <- list(dimnames(object$coefficients)[[2]], c("Estimate","Wald Std. Err.", "ChiSq","Pr(>ChiSq)"),
                              dimnames(object$coefficients)[[3]])
    for (i in 1:J) {
      output1[,1,i] <- ifelse(is.infinite(object$separation[,,i]), object$separation[,,i],object$coefficients[,,i])
      output1[,2,i] <- ifelse(is.infinite(object$separation[,,i]), NA, sqrt(diag(object$var[,,i])))
      output1[,3,i] <- object$statistic[,,i]
      output1[,4,i] <- object$pvalue[,,i]
      if (is.infinite(as.vector(object$separation)[i]))
      {
        warnSep <- TRUE;
      }
    }
  }

  if (object$method == "none")
  {
    output1 <- array(data = NA, dim = c(p,2,J))
    cat("\n Parameter Estimates and Wald Standard Errors \n \n");
    dimnames(output1) <- list(dimnames(object$coefficients)[[2]], c("Estimate","Wald Std. Error"), dimnames(object$coefficients)[[3]])
    for (i in 1:J) {
      output1[,1,i] <- ifelse(is.infinite(object$separation[,,i]), object$separation[,,i], object$coefficients[,,i])
      output1[,2,i] <- ifelse(is.infinite(object$separation[,,i]), NA, sqrt(diag(object$var[,,i])))
      if (is.infinite(as.vector(object$separation)[i]))
      {
        warnSep <- TRUE;
      }
    }
  }
  ## joint test results
  output2 <- NULL
  if (object$joint) {
    #cat("\n\n Joint Statistical Tests \n\n")
    output2.dimnames <- c(paste("Estimate (b_{i,1})",sep=""), "ChiSq", "Pr(>ChiSq)")
    output3.dimnames <-  output2.dimnames[-1]
    output2 <- array(data = NA, dim = c(p,length(output2.dimnames),3))
    output3 <- array(data = NA, dim = c(p,length(output3.dimnames),2))
    dimnames(output2) <-
      list(dimnames(object$coefficients)[[2]],output2.dimnames,
           c("H_0: b_{i,1} = b_{i,2} = ... = b_{i,J} = 0 (all zero)",
             "H_0: b_{i,1} = b_{i,2} = ... = b_{i,J} (all equal)",
             "H_0: b_{i,J} = J*b_{i,1}, b_{i,J-1} = (J-1)*b_{i,1}, ... , b_{i,2} = 2*b_{i,1} (ith covariate, proporionality)"))
    dimnames(output3) <-
      list(dimnames(object$coefficients)[[2]],output3.dimnames,
           c("H_0: b_{i,1} = b_{i,2} = ... = b_{i,J} (all equal)",
             "H_0: b_{i,J} = J*b_{i,1}, b_{i,J-1} = (J-1)*b_{i,1}, ... , b_{i,2} = 2*b_{i,1} (ith covariate, proporionality)"))

    output2[,1,1] <- NA

    testres.all0 <- object$joint.test.all0$test.h0
    print(object$joint.test.all0)
    testres.allequal <- object$joint.test.allequal$test.h0
    testres.proportion <- object$joint.test.proportion$test.h0

    output2[,"ChiSq",1] <- testres.all0[,"ChiSq"]
    output2[,"Pr(>ChiSq)",1] <- testres.all0[,"Pr(>ChiSq)"]

    for(i in 1:p){
      ## taking i-th row (corresponding to the parameter) of the first column from the i-th array
      output2[i,1,2] <- object$beta0allequal[i,1,i]
    }
    #output2[,2,2] <- NA
    output2[,"ChiSq",2] <- testres.allequal[,"ChiSq"]
    output2[,"Pr(>ChiSq)",2] <- testres.allequal[,"Pr(>ChiSq)"]
    output3[,"ChiSq",1] <- object$joint.test.allequal$test.all0.vs.constraint[,"ChiSq"]
    output3[,"Pr(>ChiSq)",1] <- object$joint.test.allequal$test.all0.vs.constraint[,"Pr(>ChiSq)"]

    for(i in 1:p){
      ## taking i-th row (corresponding to the parameter) of the first column from the i-th array
      output2[i,1,3] <- object$beta0proportion[i,1,i]
    }
    #output2[,2,3] <- NA ## Std.Errors not sure what to do
    output2[,"ChiSq",3] <- testres.proportion[,"ChiSq"]
    output2[,"Pr(>ChiSq)",3] <- testres.proportion[,"Pr(>ChiSq)"]
    output3[,"ChiSq",2] <- object$joint.test.proportion$test.all0.vs.constraint[,"ChiSq"]
    output3[,"Pr(>ChiSq)",2] <- object$joint.test.proportion$test.all0.vs.constraint[,"Pr(>ChiSq)"]

    #Proportionality is only valid with >= 2 categories...Set to NA if J < 2 by pmlr.R
    if ( is.na(output2[,1,3]) || is.na(output2[,2,3])  )
    {
      cat( "\n The Propotionality test: H_0: b_{i,J} = J*b_{i,1}, b_{i,J-1} = (J-1)*b_{i,1}, ... , b_{i,2} = 2*b_{i,1} requires at least 2 categories. \n \n")
    }
    else {
      #cat("\nH_0: b_{i,J} = J*b_{i,1}, b_{i,J-1} = (J-1)*b_{i,1}, ... , b_{i,2} = 2*b_{i,1} (ith covariate, proporionality)\n\n")
      #print(output2[,,3])
    }

  }

  if (warnSep)
  {
    cat("\n Warning: Separation has likely occurred for at least one coefficient. \n \n")
  }

  ret$coef = output1

  if(!is.null(output2)){
    ret$joint.test = output2
    ret$test.all0.vs.constraint = output3
  }

  class(ret) <- "summary.pmlr"
  return(ret)

}

## the coefficient estimates for the joint tests are ... taken from each test beta_ij
