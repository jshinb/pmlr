#' Penalized maximum likelihood estimation for multinomial logistic regression using the Jeffreys prior
#' @aliases pmlr
#' @aliases print.pmlr
#'
#' @description Extends the approach proposed by Firth (1993) for bias reduction of MLEs
#' in exponential family models to the multinomial logistic regression model with general
#' covariate types.  Modification of the logistic regression score function to remove
#' first-order bias is equivalent to penalizing the likelihood by the Jeffreys prior,
#' and yields penalized maximum likelihood estimates (PLEs) that always exist.
#' Constrained hypothesis testing is conducted via likelihood ratio statistics.
#' Profile or Wald confidence intervals (CI) are constructed for the PLEs.
#'
#' @param formula an object of class \code{\link{formula}} (or one that can be coerced to that class):
#' a symbolic description of the model to be fitted.
#' Typically, \code{formula} has the form \code{response ~ terms} where \code{response} is either a factor
#' with \eqn{J+1} levels (the first is used as the baseline category) or a \eqn{J}-column indicator matrix,
#' and \code{terms} is a series of terms which specifies a linear predictor for the response.
#' @param data a data frame containing the variables in the model
#' @param weights an optional vector of weights to be used in the fitting process.
#' Should be a numeric vector of counts in the case of frequency data.
#' @param penalized a logical variable indicating whether penalized maximum
#' likelihood should be used (default)
#' @param method a character string specifying whether p-values and confidence intervals should be based
#' on the profile likelihood, score or the Wald method. Must be one of \dQuote{likelihood} (default),
#' \dQuote{wald}, \dQuote{score} or \dQuote{none}. Note that the \dQuote{none} option provides parameter
#' estimates and Wald standard errors only, without confidence intervals or hypothesis tests.
#' This option may be useful for diagnostic purposes. In addition, the \dQuote{score} option provides
#' hypothesis tests only due to the rarity of score based confidence limits in practice.
#' @param joint a logical variable indicating whether joint hypothesis tests should be performed
#' in addition to individual parameter tests.
#' If \code{TRUE}, \deqn{H_{0}: \beta_{1i} = \cdots = \beta_{Ji} = 0} and \eqn{H_0: \beta_{1i} = \cdots = \beta_{Ji}} are also tested for covariate \eqn{i, i = 1, \ldots, P}.
#' In addition, the proportional model: \eqn{H_0: b_{i,J} = J*b_{i,1}, b_{i,J-1} = (J-1)*b_{i,1}, ... , b_{i,2} = 2*b_{i,1} }, for each \eqn{i} is tested.
#' Note that the consideration of joint hypothesis tests is only meaningful for \eqn{J > 2} categories.
#' As such, an error will occur if this is set to TRUE and binomial data is entered into the function.
#' @param CI.calculate a logical variable whether the confidence limits were computed. Default is \code{FALSE}.
#' @param CI.alpha the significance level for the profile or Wald confidence intervals (default is 0.05).
#' Relevant only if \code{CI.calcuate=TRUE}.
#' @param tol The tolerance for convergence can be specified for diagnostic purposes.
#' The default is 1e-6.  Note that this does not correspond (exactly) to the precision in each estimate,
#' as the tolerance is defined as the sum of the absolute value of the elements of the `step'
#' \eqn{=A^{-1}_{(t)}U(B^*_{(t)}} in the modified Newton-Raphson algorithm.
#' Nonetheless, higher tolerance values may reduce computation times and facilitate convergence in some cases.
#' @param verbose Logical: TRUE displays parameter values after each iteration and states when completion of each routine has occurred.
#' Default is \code{FALSE}.
#'
#' @details
#' Logistic regression is one of the most widely used regression models in practice,
#' but alternatives to conventional maximum likelihood estimation methods may be more
#' appropriate for small or sparse samples.
#' It is well known that the usual maximum likelihood estimates (MLEs) of the log odds ratios
#' are biased in finite samples, and there is a non-zero probability that an MLE is infinite
#' (i.e., does not exist).
#' This corresponds to the problem of separation (Lesaffre and Albert, 1989).
#'
#' In this package, we extend the approach proposed by Firth (1993)
#' for bias reduction of MLEs in exponential family models to the multinomial logistic regression
#' model with general covariate types.  Modification of the logistic regression score function to
#' remove first order bias is equivalent to penalizing the likelihood by the Jeffreys prior,
#' and yields penalized likelihood estimates (PLEs) that always exist, even in samples
#' in which MLEs are infinite.  PLEs are an attractive alternative in small to moderate-sized
#' samples, and are preferred to exact conditional MLEs when there are continuous covariates.
#'
#' We consider a multicategory outcome \eqn{y} that is a multinomial variable with \eqn{J + 1}
#' categories.  For each category \eqn{j (j = 1, \ldots, J)} there is a regression function
#' in which the log odds of response in category \eqn{j}, relative to category 0, is a linear
#' function of regression parameters and a vector \eqn{\bold{x}} of \eqn{P} covariates
#' (including a constant):
#' \eqn{\log\{\mathrm{prob}(y = j| \bold{x})/\mathrm{prob}(y = 0 | \bold{x})\} = \bold{\beta}_{j}^T \bold{x}}.
#' Let \eqn{\bold{y}_i} be a \eqn{J \times 1} vector of indicators for the observed response
#' category for observation \eqn{i}, with the corresponding \eqn{J \times 1}
#' vector of probabilities \eqn{\bold{\Theta}_{i} = (\Theta_{i1}, \ldots, \Theta_{iJ})^T}.
#' The vector of MLEs, \eqn{\hat{B} = vec[(\bold{\hat{\beta}}_1, \ldots, \bold{\hat{\beta}}_J)^T]},
#' is estimated from observations \eqn{(\bold{y}_i,\bold{x}_i), i = 1, \ldots, n},
#' by solving the score equations of the log-likelihood \eqn{l(B)}.
#' We denote the score function by \eqn{U(B)} and the \eqn{PJ \times PJ}
#' Fisher information matrix by \eqn{A = A(B)}.
#'
#' The order \eqn{n^{-1}} bias of estimates based on the usual likelihood \eqn{L(B)}
#' is removed by applying the penalty \eqn{|A|^{1/2}}, and basing estimation on the
#' penalized likelihood \eqn{L^*(B) = L(B) |A|^{1/2}}.
#' The vector \eqn{\hat{B}^*} of penalized estimates (PLEs) is the solution to
#' the score equations of the penalized log-likelihood \eqn{l^*(B) = l(B) + \frac{1}{2} \log{|A|}}.
#' The introduction of bias into the score function through the penalty removes the leading
#' term in the asymptotis bias of the MLEs.
#' The modified score function proposed by Firth for the binomial logistic model extends directly
#' to the multinomial model as \eqn{U^*(B) = U(B) - A \, b(B)}.
#' The bias term \eqn{b(B)} is the leading term in the asymptotic bias of the multinomial MLEs,
#' obtained from the Taylor series expansion of the log-likelihood (Cox and Snell, 1968)
#' and is a function of the matrix of third derivatives of \eqn{l(B)} with respect to \eqn{B}
#' (Bull et al., 2002).
#'
#' The PLEs are obtained by a modified Fisher scoring algorithm.  Using \eqn{t} to
#' denote the iteration number, the modified iterative updating equations are:
#'  \deqn{B^*_{(t+1)} = B^*_{(t)} + A^{-1}_{(t)}U^*(B^*_{(t)}) = B^*_{(t)} + b(B^*_{(t)}) + A^{-1}_{(t)}U(B^*_{(t)}).}
#'  Thus, in comparison to the usual updating equation used to obtain the MLEs,
#'  the score function modification operates by applying the asymptotic bias correction at each step in the iterative process.
#'  This prevents estimates from going off to infinity and failing to converge when there is separation in the data.
#'
#'  Symmetric Wald-type CIs for \eqn{\beta_{jp}} can be constructed using \eqn{Var^{* 1/2}(\hat{\beta}^*_{jp})},
#'  obtained from the inverse of \eqn{A^*}, however, performance is expected to be poor in situations,
#'  where separation is likely to occur.
#'  Asymmetric CIs for the PLEs can be constructed from the profile log-likelihood for
#'  \eqn{\beta_{jp}}, which is the function \eqn{l^*_0(B(s))}, where \eqn{B(s)}
#'  is the argument that maximizes \eqn{l^*} under the single-parameter constraint
#'  \eqn{H_0: \beta_{jp} = s}.  The \eqn{100(1 - \alpha)\%} CI for \eqn{\beta_{jp}}
#'  is given by all parameter values that are compatible with the data
#'  (i.e., all \eqn{s} such that the likelihood ratio statistic \eqn{LR_P(s) \le q},
#'  where \eqn{q} is the \eqn{1 - \alpha} percentile of the \eqn{\chi^2} distribution).
#'  This is equivalent to \eqn{l^*_0(B(s)) \ge l^*(\hat{B}^*) - \frac{1}{2} q}.
#'  The endpoints of the interval are then found by numerically solving the equality for
#'  values of \eqn{s}.  Based on the algorithm employed in SAS PROC LOGISTIC for MLEs,
#'  our method for finding these roots does not require computing \eqn{l^*_0(B(s))},
#'  which in itself would involve maximizing \eqn{l^*(B)},
#'  but proceeds directly by solving the constrained optimization problem:
#'  maximize \eqn{l^*(B)} such that \eqn{l^*_(B) = l^*(\hat{B}^*) - \frac{1}{2} q} and \eqn{\beta_{jp} = s}.
#'  We, however, modify the starting values in the iterative scheme used by SAS to obtain a
#'  new algorithm that is slower, but simple and more robust (see Bull et al. 2007 for details).
#'
#' @return
#' \code{pmlr} returns an object of class \dQuote{\code{pmlr}}, a list with the following components:
#'   \item{call}{the matched call.}
#'   \item{method}{a character string indicating the method used to carry out hypothesis testing and/or confidence limit calcuation.}
#'   \item{separation}{an array indicating coefficients for which separation occured;
#'   \code{NA} returned in the absense of separation.}
#'   \item{converged}{a logical value indicating whehter the IWLS algorithm judged to have converged.}
#'   \item{coefficients}{an array containing the coefficients of the \eqn{p} parameters for the \eqn{J} categories.}
#'   \item{var}{an array containing the variance-covariance matrices for the \eqn{J} categories. If \code{penalized=TRUE}, \code{var} is obtained based on \eqn{A^*}, the information matrix for the PLEs.}
#'   \item{CI.calculate}{a logical value whether the confidence limits were computed.}
#'   \item{CI.lower}{returned only if \code{CI.calculate=TRUE}: an array containing the lower confidence limits from the individual parameter tests.}
#'   \item{CI.upper}{returned only if \code{CI.calculate=TRUE}: an array containing the upper confidence limits from the individual parameter tests.}
#'   \item{statistic}{an array containing the test statistics from the individual parameter tests.}
#'   \item{pvalue}{an array containing the p-values from the individual parameter tests.}
#'   \item{logLik}{the value of the log-likelihood function for the fitted model (under no constraints).}
#'   \item{df}{the degrees of freedom, i.e., the number of estimated parameters in the model.}
#'   \item{converged.H0}{an array containing the logical values whether the fitting algorithm for each null model is jusdged to have converged.}
#'   \item{logLik.H0}{an array containing the value of the log-likelihood function for the fitted model under each null hypothesis.}
#'   \item{joint}{a logical value indicating whether the joint hypothesis tests were performed.}
#'   \item{beta0all0}{When a joint likelihood test of all betas = 0 is called, the estimated betas are provided.  This is for information only and not displayed in the output.}
#'   \item{var0all0}{When a joint likelihood test of all betas = 0 is called, the estimated variances are also provided.  This is for information only and not displayed in the output.}
#'   \item{beta0allequal}{same as \code{'beta0all0'} except for the null hypothesis that all betas are equal.}
#'   \item{var0allequal}{same as var0all0 except for the null hypothesis that all betas are equal.}
#'   \item{beta0proportion}{same as beta0all0 except for the null hypothesis that all betas are proportional.}
#'   \item{var0proportion}{same as var0all0 except for the null hypothesis that all betas are proportional.}
#'   \item{logLik.C}{array contatining the values of log likelihood function under constraints: betas = 0; all betas are equal; and betas are proportional.}
#'   \item{joint.test.all0}{a list containing the following three components evaluated for the constrained hypothesis that all betas = 0.}
#'   \item{joint.test.all0$h0}{a character string describing the constrained hypothesis}
#'   \item{joint.test.all0$converged}{an array contatining whether the fitting algorithm achieved convergence for each constrained hypothesis.}
#'   \item{joint.test.all0$test.h0}{an array contatining the test statistics and p-values from constrained hypothesis tests all betas = 0. }
#'
#'   \item{joint.test.allequal}{same as \code{'joint.test.all0'} except for the constrained hypothesis that that all betas are equal components evaluated for
#'   the constrained hypothesis that all betas are equal, and for the additional component \code{'test.all0.vs.constraint'} (see below).}
#'   \item{joint.test.allequal$test.all0.vs.constraint}{a data frame with test statistics and p-values for comparing the likelihoods for
#'   all betas are equal vs. all betas are zero.}
#'   \item{joint.test.proportion}{same as \code{'joint.test.allequal'} except for the constrained hypothesis
#'   that betas are proportional.}
#'
#' @note This implementation is not a true scoring or Newton-type algorithm because
#' it updates with the inverse of \eqn{A}, the Fisher information matrix for the MLEs,
#' rather than the information for the PLEs, \eqn{A^*}, which includes an additional term
#' corresponding to the second derivatives of the penalty: \eqn{\frac{1}{2} \log |A|}.
#' As a result, convergence using the modified scoring algorithm for the PLEs is slower
#' than a scoring algorithm based on \eqn{A^*}, which converges at a quadratic rate.
#' For well behaved and larger datasets this usually means no more than 2 or 3 steps
#' beyond that required for the MLEs.  Starting values of \eqn{\beta{jp} = 0} are used
#' and are generally satisfactory.  For smaller datasets (i.e., less than 50 observations)
#' and especially datasets in which there are infinite MLEs, convergence is slower
#' and could take up to 35 or 40 iterations.  In datasets with separations, starting values
#' other than 0 can lead to divergence problems.
#'
#' The first argument to the \code{pmlr} function is a formula of the form \code{response ~ terms}
#' where \code{response} can be either a \eqn{J}-column indicator matrix or a factor with
#' \eqn{J + 1} levels.  In the case of frequency data (e.g., the \code{hepatitis} dataset)
#' the baseline category is determined by the \eqn{J}-column indicator matrix with baseline
#' category 0.  In the case of data with individual records (e.g., the \code{enzymes} dataset)
#' the baseline category is determined by the outcome, which is coded as a factor.
#' The first level (i.e., \code{levels(x)[1]} is taken to be the baseline.
#' Note that the default levels attribute of a factor is the unique set of values taken
#' \code{as.character(x)}, sorted into increasing order of \code{x}.
#' For example, the enzymes dataset has categories 1, 2, and 3 with baseline category 1.
#' The baseline category can be changed to category 2 via
#' \code{enzymes$Group <- factor(enzymes$Group, levels-c("2","1","3"))}.
#'
#' @references Bull, S. B., Greenwood, C. M. T., and Hauck, W. W. (1997)
#' Jacknife bias reduction for polychotomous logistic regression. \emph{Statistics in Medicine},
#' \bold{16}, 545--560; Bull, S. B. (1997) Correction.
#' \emph{Statistics in Medicine}, \bold{16}, 2928.
#'
#' Bull, S. B., Mak, C. and Greenwood, C. M. T. (2002) A modified score function estimator
#' for multinomial logistic regression in small samples.
#' \emph{Computational Statistics & Data Analysis}, \bold{39}, 57--74.
#'
#' Bull, S. B., Lewinger, J. P., Lee, S. S. F. (2005) Penalized maximum likelihood estimation for multinomial logistic regression using the Jeffreys prior.  \emph{Technical Report No. 0505}, Department of Statistics, University of Toronto.
#'
#' Bull, S. B., Lewinger, J. P. and Lee, S. S. F. (2007)
#' Confidence intervals for multinomial logistic regression in sparse data.
#' \emph{Statistics in Medicine}, \bold{26}, 903--918.
#'
#' Cox, D. R. and Snell, E. J. (1968)
#' A general definition of residuals.
#' \emph{Journal of the Royal Statistical Society, Series B}, \bold{30}, 248--275.
#'
#' Firth, D. (1993)
#' Bias reduction of maximum likelihood estimates.
#' \emph{Biometrika}, \bold{80}, 27--38.
#'
#' Lesaffre, E. and Albert, A. (1989)
#' Partial separation in logistic discrimination.
#' \emph{Journal of the Royal Statistical Society, Series B}, \bold{51}, 109--116.
#'
#' SAS Institute Inc. (1999)
#' \emph{SAS OnlineDoc, Version 8, The LOGISTIC Procedure},
#' Confidence Intervals for Parameters, Chapter 39, Section 26, Cary, NC.
#'
#' @examples
#' data(hepatitis)
#'
#' # As reported in Bull et al. (2007)
#' fit <- pmlr(cbind(HCV, nonABC) ~ group + time + group:time, data = hepatitis,
#' weights = counts, method="score")
#' summary(fit)
#'
#' data(enzymes)
#' # Exclude patients in Group 4 (post-necrotic cirrhosis)
#' enzymes <- enzymes[enzymes$Group != 4,]
#'
#' # Center and scale covariates
#' AST <- scale(log(enzymes$AST))
#' ALT <- scale(log(enzymes$ALT))
#' GLDH <- scale(log(enzymes$GLDH))
#' OCT <- scale(log(enzymes$OCT))
#' enzymes <- data.frame(Patient = enzymes$Patient,
#'                       Group = enzymes$Group, AST, ALT, GLDH, OCT)
#'
#' # Remove 10 observations to create separation
#' enzymes <- enzymes[-c(9, 18, 33, 58, 61, 77, 94, 97, 99, 100),]
#'
#' # Multinomial: Acute viral hepatitis and aggressive chronic hepatits
#' # vs. persistent chronic hepatitis
#' # Assign Group 2 (persistent chronic hepatitis) as baseline category
#' enzymes$Group <- factor(enzymes$Group, levels=c("2","1","3"))
#' fit <- pmlr(Group ~ AST + GLDH, data = enzymes, method="wald")
#' summary(fit)
#'
#' # Binomial: Acute viral hepatitis vs. persistent chronic hepatitis
#' # Exclude patients in Group 3 (agressive chronic hepatitis)
#' enzymes.1vs2 <- enzymes[enzymes$Group != 3,]
#' # Assign Group 2 (persistent chronic hepatitis) as baseline category
#' enzymes.1vs2$Group <- factor(enzymes.1vs2$Group, levels=c("2","1"))
#' fit <- pmlr(Group ~ AST + GLDH, data = enzymes.1vs2, method="none")
#' summary(fit)
#'
#' # Binomial: Aggressive chronic hepatitis vs. persistent chronic hepatitis
#' # Exclude patients in Group 1 (acute viral hepatitis)
#' enzymes.3vs2 <- enzymes[enzymes$Group != 1,]
#' # Assign Group 2 (persistent chronic hepatitis) as baseline category
#' enzymes.3vs2$Group <- factor(enzymes.3vs2$Group, levels=c("2","3"))
#' fit <- pmlr(Group ~ AST + GLDH, data = enzymes.3vs2, method="none")
#' summary(fit)
#'
#' @seealso \code{\link{summary.pmlr}}
#'
#' @keywords models
#' @keywords regression
#' @keywords htest
#' @export
#-------------------------------------------------

pmlr <- function(formula, data, weights = NULL, penalized = TRUE,
                 method = c("likelihood", "wald", "score", "none")[1], joint = FALSE,
                 CI.calculate = FALSE, CI.alpha = 0.05,
                 tol = 1e-6, verbose=FALSE) {

  if(!method %in% c("likelihood","score","wald","none")) stop("invalid method: please select a valid one")

  ret <- list()
  # change this!
  if (penalized) useAstar.wald <- TRUE
  else useAstar.wald <- FALSE

  useAstar.LRCI <- FALSE

  # check the data and update the formula


  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weights"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  y <- model.response(mf); if (is.factor(y)) y <- indicator.matrix(y)
  x <- model.matrix(mt, mf)
  J0 <- ncol(y)

  # check if there are no observations for
  # response matrix:
  exclude.cols.y <- exclude.cols.x <- NULL
  for(y.i in 1:ncol(y)){
    exclude.cols.y <- c(exclude.cols.y, var(y[,y.i],na.rm=T)==0)
  }
  if( any(exclude.cols.y) ) y <- y[,!exclude.cols.y,drop=FALSE]

  for(x.i in 1:ncol(x)){
    exclude.cols.x <- c(exclude.cols.x, var(x[,x.i],na.rm=T)==0)
  }
  if( any(exclude.cols.x) ) x <- x[,c(1,which(!exclude.cols.x)),drop=FALSE]


  wt <- as.vector(model.weights(mf)); if (is.null(wt)) wt <- rep(1, times = dim(y)[1])
  p <- ncol(x) # p = number of covariates
  J <- ncol(y) # J =  (response-categories - 1)

  if(J0>J & J<=1 ) {
    warning('The response variable is not multinomial; we will let joint=FALSE')
    joint = FALSE
  }

  if ( (J <=1) && (joint) )
    stop("Sorry, it is not possible to perform joint hypothesis tests for binomial data.  Please set joint=FALSE to proceed.")

  if (penalized) {
    fit <- getPLEs(x, y, wt, tol=tol, verbose=verbose); B <- fit$B; B.inf <- NULL; Ainv <- fit$Ainv; Astarinv <- fit$Astarinv; l.max <- fit$lstar.max; fit.conv <- fit$conv
    if (verbose) { cat("GetPLEs Complete. \n") }
  } else {
    fit <- getMLEs(x, y, wt, tol=tol, verbose=verbose); B <- fit$B; B.inf <- fit$B.inf; Ainv <- fit$Ainv; l.max <- fit$l.max; fit.conv <- fit$conv
    if (verbose) { cat("GetMLEs Complete. \n") }
  }

  separation <- array(data = NA, dim = c(1,p,J))
  dimnames(separation)<-list("", colnames(x), colnames(y))
  if (!penalized) {
    for (i in 1:J) {
      separation[1,,i] <- t(is.infinite(B.inf))[i,]
      separation[1,,i] <- ifelse(separation[1,,i], t(B.inf)[i,], NA)
    }
  }

  # parameter estimation
  coef <- array(data = NA, dim = c(1,p,J))
  dimnames(coef)<-list("", colnames(x), colnames(y))
  for (i in 1:J) coef[1,,i] <- t(B)[i,]
  # variance estimation
  var <- array(data = NA, dim = c(p,p,J))
  dimnames(var)<-list(colnames(x), colnames(x), colnames(y))
  for (i in 1:J) var[,,i] <- Ainv[seq(from = i, by = J, length = p), seq(from = i, by = J, length = p)]
  # profile-likleihood CI estimation
  if(CI.calculate) {
    CI.lower <- CI.upper <- array(data = NA, dim = c(1,p,J))
    dimnames(CI.lower) <- dimnames(CI.upper) <-list("", colnames(x), colnames(y))
  }
  # test
  stat <- pval <- array(data = NA, dim = c(1,p,J))
  dimnames(stat) <- dimnames(pval) <- list("", colnames(x), colnames(y))

  if (joint){
    ## JS Added 2015-12-15
    ## to create the following objects only when performing the joint test
    beta0all0 <- var0all0 <- beta0allequal <- var0allequal <- beta0proportion <- var0proportion <- array(data = NA, dim = c(1,p,J))
    dimnames(beta0all0) <- dimnames(var0all0) <- dimnames(beta0allequal) <- dimnames(var0allequal) <- dimnames(beta0proportion) <- dimnames(var0proportion) <-
      list("", colnames(x), colnames(y))
    list("", colnames(x), c("H_0: b_{i,1} = b_{i,2} = ... = b_{i,J} = 0",
                            "H_0: b_{i,1} = b_{i,2} = ... = b_{i,J}",
                            "H_0: b_{i,J} = J*b_{i,1}, b_{i,J-1} = (J-1)*b_{i,1}, ... , b_{i,2} = 2*b_{i,1}"))
    ## stat.all0.vs.constraint and pval.all0.vs.constraint need not be evalueated
    ## under all0-constrianed hypothesis - "H_0: b_{i,1} = b_{i,2} = ... = b_{i,J} = 0"
    stat.all0.vs.constraint <- pval.all0.vs.constraint <- stat.joint <- pval.joint <- l0.joint <- array(data = NA, dim = c(1,p,3))
    dimnames(stat.all0.vs.constraint) <- dimnames(pval.all0.vs.constraint) <- dimnames(stat.joint) <- dimnames(pval.joint) <- dimnames(l0.joint) <-
      list("", colnames(x),
           c("H_0: b_{i,1} = b_{i,2} = ... = b_{i,J} = 0",
             "H_0: b_{i,1} = b_{i,2} = ... = b_{i,J}",
             "H_0: b_{i,J} = J*b_{i,1}, b_{i,J-1} = (J-1)*b_{i,1}, ... , b_{i,2} = 2*b_{i,1}"))

    stat.all0.vs.constraint <- pval.all0.vs.constraint <- pval.all0.vs.constraint[,,-1,drop=FALSE]
    test.h0.all0 <- test.h0.allequal <- test.h0.proportion <- list()
  }

  converged <- NA
  if (method == "likelihood") {
    ## Likelihood ratio test for H_0: b_{i,j} = 0 (ith covariate, jth category)
    testRun <- test.LR(x, y, wt, mt, B, h0 = 1, penalized, tol = tol, verbose=verbose);
    if( !any(testRun$conv) )
      warning("test.LR: algorithm did not converge\n")
    if (verbose) { cat("test.LR Complete. \n") }
    for (i in 1:J) {
      stat[1,,i] <- t(testRun$statistic)[i,]
      pval[1,,i] <- t(testRun$pvalue)[i,]
    }
    test.unconstrained <- list()
    test.unconstrained$statistic <- stat
    test.unconstrained$pvalue <- pval
    test.unconstrained$logLik.H0 <- NA
    test.unconstrained$converged <- testRun$conv
    if(J==1){
      #if J > 2 - individual l0 is rather meaningless
      #cat("J==1 \n")
      l0 <- array(data = NA, dim = c(1,p,J))
      dimnames(l0) <- list("", colnames(x), colnames(y))
      l0[,,1] <- testRun$l0
      test.unconstrained$logLik.H0 <- l0
    }

    if(CI.calculate){
      # Profile confidence intervals - only when method = "likelihood" ??
      profileCI.lower <- profileCIs(x, y, wt, B, B.inf, side = -1, CI.alpha, l.max, step = 0.05,
                                    useAstar = useAstar.LRCI, penalized, tol = tol, verbose=verbose)
      if (verbose) { cat("profileCI.lower Complete. \n") }
      profileCI.upper <- profileCIs(x, y, wt, B, B.inf, side = 1, CI.alpha, l.max, step = 0.05,
                                    useAstar = useAstar.LRCI, penalized, tol = tol, verbose=verbose)
      if (verbose) { cat("profileCI.upper Complete. \n") }
      for (i in 1:J) {
        CI.lower[1,,i] <- t(profileCI.lower)[i,]
        CI.upper[1,,i] <- t(profileCI.upper)[i,]
      }
    }#if(profileCI.calculate) ends
    if (joint) {
      # h0=2: Likelihood ratio test for H_0: b_{i,1} = b_{i,2} = ... = b_{i,J} = 0 (ith covariate)
      testRun <- test.LR(x, y, wt, mt, B, h0 = 2, penalized, tol = tol, verbose=verbose);

      test.h0.colnames = c("ChiSq","Pr(>ChiSq)")
      test.h0.all0$h0 <- "H_0: b_{i,1} = b_{i,2} = ... = b_{i,J} = 0 (ith covariate)"
      test.h0.all0$converged <- testRun$conv
      test.h0.all0$test.h0 <- matrix(NA,ncol=length(test.h0.colnames),nrow=p)
      dimnames(test.h0.all0$test.h0) <- list( colnames(x), test.h0.colnames )
      #test.h0.all0$test.h0[,"logLik"] <-
      l0.joint[,,1] <- t(testRun$l0)
      test.h0.all0$test.h0[,"ChiSq"] <- t(testRun$statistic)
      test.h0.all0$test.h0[,"Pr(>ChiSq)"] <- t(testRun$pvalue)
      beta0all0 <- testRun$beta0.array
      var0all0 <- testRun$var0.array

      # h0=3: Likelihood ratio test for H_0: b_{i,1} = b_{i,2} = ... = b_{i,J} (ith covariate)
      testRun <- test.LR(x, y, wt, mt, B, h0 = 3, penalized, tol = tol, verbose=verbose);
      beta0allequal <- testRun$beta0.array
      var0allequal <- testRun$var0.array
      test.h0.allequal$h0 <- "H_0: b_{i,1} = b_{i,2} = ... = b_{i,J} (ith covariate)"
      test.h0.allequal$converged <- testRun$conv
      test.h0.colnames <- c("ChiSq","Pr(>ChiSq)")
      test.h0.allequal$test.all0.vs.constraint <-
        test.h0.allequal$test.h0 <- matrix(NA,ncol=length(test.h0.colnames),nrow=p)
      dimnames(test.h0.allequal$test.all0.vs.constraint) <- dimnames(test.h0.allequal$test.h0) <- list(colnames(x), test.h0.colnames)
      #test.h0.allequal$test.h0[,"logLik"] <-
      l0.joint[,,2] <- t(testRun$l0)
      test.h0.allequal$test.h0[,"ChiSq"] <- t(testRun$statistic)
      test.h0.allequal$test.h0[,"Pr(>ChiSq)"] <- t(testRun$pvalue)
      chisq <- -2*(l0.joint[,,1] - t(testRun$l0))
      test.h0.allequal$test.all0.vs.constraint[,"ChiSq"] <- chisq
      test.h0.allequal$test.all0.vs.constraint[,"Pr(>ChiSq)"] <- pchisq(chisq, df=1, lower.tail=F)


      # h0=4: Likelihood ratio test for proportionality: H_0: b_{i,J} = J*b_{i,1}, b_{i,J-1} = (J-1)*b_{i,1}, ... , b_{i,2} = 2*b_{i,1} (ith covariate)
      if (J >= 2) {
        testRun <- test.LR(x, y, wt, mt, B, h0 = 4, penalized, tol = tol, verbose=verbose);
        l0.joint[,,3] <- t(testRun$l0)
        test.h0.proportion$h0 <- "H_0: b_{i,J} = J*b_{i,1}, b_{i,J-1} = (J-1)*b_{i,1}, ... , b_{i,2} = 2*b_{i,1} (ith covariate)"
        test.h0.proportion$converged <- testRun$conv
        test.h0.colnames <- c("ChiSq","Pr(>ChiSq)")
        test.h0.proportion$test.all0.vs.constraint <-
          test.h0.proportion$test.h0 <- matrix(NA,ncol=length(test.h0.colnames),nrow=p)
        dimnames(test.h0.proportion$test.all0.vs.constraint) <-
          dimnames(test.h0.proportion$test.h0) <-
          list(colnames(x), test.h0.colnames)

        test.h0.proportion$test.h0[,"ChiSq"] <- t(testRun$statistic)
        test.h0.proportion$test.h0[,"Pr(>ChiSq)"] <- t(testRun$pvalue)
        chisq <- -2*(l0.joint[,,1] - l0.joint[,,3])
        test.h0.proportion$test.all0.vs.constraint[,"ChiSq"] <- chisq
        test.h0.proportion$test.all0.vs.constraint[,"Pr(>ChiSq)"] <- pchisq(chisq, df=1, lower.tail=F)
        beta0proportion <- testRun$beta0.array
        var0proportion <- testRun$var0.array
      }
      else {#if (J < 2) ?? no test ??
        #test.h0.proportion$test.h0[,"logLik"] <- NA
        test.h0.proportion$test.h0[,"ChiSq"] <- NA
        test.h0.proportion$test.h0[,"Pr(>ChiSq)"] <- NA
        test.h0.proportion$test.all0.vs.constraint[,"ChiSq"] <- NA
        test.h0.proportion$test.all0.vs.constraint[,"Pr(>ChiSq)"] <- NA
      }
      if (verbose) { cat("test.LR (joint) Complete. \n") }
    }
  }

  #*** NEED TO THINK ABOUT THIS

  if (useAstar.wald) Ainv <- Astarinv# ** when is it being used?
  # Wald test
  if (method == "wald") {
    ## Wald test for H_0: b_{i,j} = 0 (ith covariate, jth category)
    for (i in 1:J) {
      stat[,,i] <- (coef[,,i]/sqrt(diag(var[,,i])))^2
      pval[,,i] <- pchisq(stat[,,i], df = 1, lower.tail = FALSE)
    }
    if(CI.calculate){
      CI.lower <-  CI.upper <- array(data = NA, dim = c(1,p,J))
      dimnames(CI.lower) <- dimnames(CI.upper) <-list("", colnames(x), colnames(y))
      for (i in 1:J) {
        CI.lower[,,i] <- coef[,,i] - qnorm(p = 1 - CI.alpha/2) * sqrt(diag(var[,,i]))
        CI.upper[,,i] <- coef[,,i] + qnorm(p = 1 - CI.alpha/2) * sqrt(diag(var[,,i]))
      }
    }
    test.unconstrained <- list()
    test.unconstrained$statistic <- stat
    test.unconstrained$pvalue <- pval

    if (joint) {
      # Wald test for H_0: b_{i,1} = b_{i,2} = ... = b_{i,J} = 0 (ith covariate)
      for (i in 1:p) {
        C <- matrix(data = 0, nrow = J, ncol = p * J); C[1:J,(((i - 1) * J) + 1):(i * J)] <- diag(J)
        stat.joint[1,i,1] <- test.wald(vec(t(B)), Ainv, C)
        pval.joint[1,i,1] <- pchisq(stat.joint[1,i,1], df = J, lower.tail = FALSE)
      }
      test.h0.all0$h0 <- "H_0: b_{i,1} = b_{i,2} = ... = b_{i,J} = 0 (ith covariate)"
      test.h0.colnames <- c("ChiSq","Pr(>ChiSq)")
      test.h0.all0$test.h0 <- matrix(NA,ncol=length(test.h0.colnames),nrow=p)
      dimnames(test.h0.all0$test.h0) <- list(colnames(x), test.h0.colnames)
      test.h0.all0$test.h0[,"ChiSq"] <- stat.joint[,,1]
      test.h0.all0$test.h0[,"Pr(>ChiSq)"] <- pval.joint[,,1]

      # Wald test for H_0: b_{i,1} = b_{i,2} = ... = b_{i,J} for the ith covariate
      for (i in 1:p) {
        C <- matrix(data = 0, nrow = J-1, ncol = p *J); C[1:(J - 1), ((i - 1) * J + 1):((i * J) - 1)] <- diag(J - 1); for (j in (1:J - 1)) C[j,((i - 1) * J) + 1 + j] <- (-1)
        stat.joint[1,i,2] <- test.wald(vec(t(B)), Ainv, C)
        pval.joint[1,i,2] <- pchisq(stat.joint[1,i,2], df = J - 1, lower.tail = FALSE)
      }
      test.h0.allequal$h0 <- "H_0: b_{i,1} = b_{i,2} = ... = b_{i,J} (ith covariate)"
      test.h0.allequal$test.h0 <- matrix(NA,ncol=length(test.h0.colnames),nrow=p)
      dimnames(test.h0.allequal$test.h0) <- list(colnames(x), test.h0.colnames)
      test.h0.allequal$test.h0[,"ChiSq"] <- stat.joint[,,2]
      test.h0.allequal$test.h0[,"Pr(>ChiSq)"] <- pval.joint[,,2]

      # Wald test for proportionality, H_0: b_{i,J} = J*b_{i,1}, b_{i,J-1} = (J-1)*b_{i,1}, ... , b_{i,2} = 2*b_{i,1} (ith covariate)
      if (J >= 2) {
        for (i in 1:p) {
          C <- matrix(data = 0, nrow = J-1, ncol = p *J)
          C[1:(J - 1), ((i - 1) * J + 1):((i * J) - 1)] <- diag(J - 1);
          for (j in (1:J - 1)) C[j,((i - 1) * J) + 1 + j] <- (-j/(j+1))
          stat.joint[1,i,3] <- test.wald(vec(t(B)), Ainv, C)
          pval.joint[1,i,3] <- pchisq(stat.joint[1,i,3], df = J - 1, lower.tail = FALSE)
        } }
      else {
        stat.joint[1,,3] <- NA
        pval.joint[1,,3] <- NA
      }
      test.h0.proportion$h0 <- "H_0: b_{i,J} = J*b_{i,1}, b_{i,J-1} = (J-1)*b_{i,1}, ... , b_{i,2} = 2*b_{i,1} (ith covariate)"
      test.h0.proportion$test.h0 <- matrix(NA,ncol=length(test.h0.colnames),nrow=p)
      dimnames(test.h0.proportion$test.h0) <- list(colnames(x), test.h0.colnames)
      test.h0.proportion$test.h0[,"ChiSq"] <- stat.joint[,,3]
      test.h0.proportion$test.h0[,"Pr(>ChiSq)"] <- pval.joint[,,3]

      if (verbose) { cat("Wald Tests (joint) Complete. \n") }
    }#if(joint) ends
  }

  ####------------------------------ NEED TO FIX!!! ------------------------------####
  if (method == "score") {
    ## Score test for H_0: b_{i,j} = 0 (ith covariate, jth category)
    ## carry out estimation and score test undr the unconstrained hypothesis when the user intends to
    testRun <- test.score(x, y, wt, mt, B, h0 = 1, penalized, verbose=verbose);
    for (i in 1:J) {
      stat[1,,i] <- t(testRun$statistic)[i,];
      pval[1,,i] <- t(testRun$pvalue)[i,];
    }
    test.unconstrained <- list()
    test.unconstrained$statistic <- stat
    test.unconstrained$pvalue <- pval
    test.unconstrained$converged <- testRun$conv

    if(J==1){
      #if J > 2 - individual l0 is rather meaningless
      #cat("J==1 \n")
      l0 <- array(data = NA, dim = c(1,p,J))
      dimnames(l0) <- list("", colnames(x), colnames(y))
      l0[,,1] <- testRun$l0
      test.unconstrained$logLik.H0 <- l0
    }

    ##***NEW - added by JS (Oct.28, 2015)
    if( !any(testRun$conv) )
      warning("test.score: algorithm did not converge\n")
    if (verbose) { cat("Score Tests Complete. \n") }

    if (joint) {
      ## Score test for H_0: b_{i,1} = b_{i,2} = ... = b_{i,J} = 0 (ith covariate)
      testRun <- test.score(x, y, wt, mt, B, h0 = 2, penalized, verbose = verbose)
      test.h0.all0$h0 <- "H_0: b_{i,1} = b_{i,2} = ... = b_{i,J} = 0 (ith covariate)"
      test.h0.all0$converged <- testRun$conv
      test.h0.colnames <- c("ChiSq","Pr(>ChiSq)")
      test.h0.all0$test.h0 <- matrix(NA,ncol=length(test.h0.colnames),nrow=p)
      dimnames(test.h0.all0$test.h0) <- list(colnames(x), test.h0.colnames)
      #test.h0.all0$test.h0[,"logLik"] <-
      l0.joint[,,1] <- t(testRun$l0)
      test.h0.all0$test.h0[,"ChiSq"] <- t(testRun$statistic)
      test.h0.all0$test.h0[,"Pr(>ChiSq)"] <- t(testRun$pvalue)
      beta0all0 <- testRun$beta0.array
      var0all0 <- testRun$var0.array
      if(!any(testRun$conv))
        warning("test.score: algorithm did not converge - joint test for all beta=0 \n")

      ## Score test for H_0: b_{i,1} = b_{i,2} = ... = b_{i,J} (ith covariate)
      testRun <- test.score(x, y, wt, mt, B, h0 = 3, penalized, verbose = verbose);
      test.h0.allequal$h0 <- "H_0: b_{i,1} = b_{i,2} = ... = b_{i,J} (ith covariate)"
      test.h0.allequal$converged <- testRun$conv
      test.h0.colnames <- c("ChiSq","Pr(>ChiSq)")
      test.h0.allequal$test.all0.vs.constraint <-
        test.h0.allequal$test.h0 <- matrix(NA,ncol=length(test.h0.colnames),nrow=p)
      dimnames(test.h0.allequal$test.all0.vs.constraint) <-
        dimnames(test.h0.allequal$test.h0) <- list(colnames(x), test.h0.colnames)
      #test.h0.allequal$test.h0[,"logLik"] <-
      l0.joint[,,2] <- t(testRun$l0)
      test.h0.allequal$test.h0[,"ChiSq"] <- t(testRun$statistic)
      test.h0.allequal$test.h0[,"Pr(>ChiSq)"] <- t(testRun$pvalue)
      chisq <- -2*(l0.joint[,,1] - l0.joint[,,2])
      test.h0.allequal$test.all0.vs.constraint[,"ChiSq"] <- chisq
      test.h0.allequal$test.all0.vs.constraint[,"Pr(>ChiSq)"] <- pchisq(chisq, df=1, lower.tail=F)
      beta0allequal <- testRun$beta0.array
      var0allequal <- testRun$var0.array

      if(!any(testRun$conv))
        warning("test.score: algorithm did not converge - joint test for all beta equal \n")

      ## Score test for proportionality, H_0: b_{i,J} = J*b_{i,1}, b_{i,J-1} = (J-1)*b_{i,1}, ... , b_{i,2} = 2*b_{i,1} (ith covariate)
      if (J >= 2) {
        testRun <- test.score(x, y, wt, mt, B, h0 = 4, penalized, verbose = verbose);
        test.h0.proportion$h0 <- "H_0: b_{i,J} = J*b_{i,1}, b_{i,J-1} = (J-1)*b_{i,1}, ... , b_{i,2} = 2*b_{i,1} (ith covariate)"
        test.h0.proportion$converged <- testRun$conv
        test.h0.colnames <- c("ChiSq","Pr(>ChiSq)")
        test.h0.proportion$test.all0.vs.constraint <- test.h0.proportion$test.h0 <- matrix(NA,ncol=length(test.h0.colnames),nrow=p)
        dimnames(test.h0.proportion$test.all0.vs.constraint) <-
          dimnames(test.h0.proportion$test.h0) <- list(colnames(x), test.h0.colnames)
        #test.h0.proportion$test.h0[,"logLik"] <-
        l0.joint[,,3] <- t(testRun$l0)
        test.h0.proportion$test.h0[,"ChiSq"] <- t(testRun$statistic)
        test.h0.proportion$test.h0[,"Pr(>ChiSq)"] <- t(testRun$pvalue)
        chisq <- -2*(l0.joint[,,1]- l0.joint[,,3])
        test.h0.proportion$test.all0.vs.constraint[,"ChiSq"] <- chisq
        test.h0.proportion$test.all0.vs.constraint[,"Pr(>ChiSq)"] <- pchisq(chisq, df=1, lower.tail=F)
        beta0proportion <- testRun$beta0.array
        var0proportion <- testRun$var0.array
        if( !any(testRun$conv) )
          warning("test.score: algorithm did not converge - joint test for proportion beta \n")

      }
      else {
        #test.h0.proportion$test.h0[,"logLik"] <- NA
        test.h0.proportion$test.h0[,"ChiSq"] <- NA
        test.h0.proportion$test.h0[,"Pr(>ChiSq)"] <- NA
        test.h0.proportion$test.all0.vs.constraint[,"ChiSq"] <- NA
        test.h0.proportion$test.all0.vs.constraint[,"Pr(>ChiSq)"] <- NA
      }

      if (verbose) { cat("Score Tests (joint) Complete. \n") }
    }
  }#if(method == "score") ends


  if (method == "none")  { #No CIs or Hypothesis tests are performed.
    ret <- list(coefficients = coef, var = var, separation = separation, call = call, method = method, joint = joint)
  }
  else{ #(method != "none")
    ret$call = call
    ret$method = method
    ret$separation <- separation
    ret$converged = fit.conv
    ret$coefficients = coef
    ret$var = var
    ret$CI.calculate = CI.calculate
    if(CI.calculate){
      ret$CI.lower <- CI.lower
      ret$CI.upper <- CI.upper
    }
    ##test
    ret$statistic = test.unconstrained$statistic
    ret$pvalue = test.unconstrained$pvalue
    ret$logLik = l.max
    ret$df=length(colnames(x))*J ## log-likelihood under no constraints ** NEW
    ret$converged.H0 = test.unconstrained$converged
    ret$logLik.H0 = test.unconstrained$logLik.H0

    ret$joint = joint
    if(joint){
      ret$beta0all0 = beta0all0
      ret$var0all0 = var0all0

      ret$beta0allequal = beta0allequal
      ret$var0allequal = var0allequal

      ret$beta0proportion = beta0proportion
      ret$var0proportion = var0proportion

      ## Joint test results
      ret$logLik.joint = l0.joint
      ## joint test more info
      ret$joint.test.all0 = test.h0.all0
      ret$joint.test.allequal = test.h0.allequal
      ret$joint.test.proportion = test.h0.proportion
    } #if(joint) ends

  }

  #attr(ret, "class") <- c("pmlr")
  class(ret) <- "pmlr"
  return(ret)
}
