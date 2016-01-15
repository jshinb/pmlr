getPLEs <- function(x, y, wt, tol=1e-6, verbose=FALSE,...) {
  
  p <- ncol(x) # p covariates
  J <- ncol(y) # J categories
  
  B <- matrix(data = 0, nrow = p, ncol = J)
  vecB <- vec(t(B))
  
  ys <- 1 - apply(y, 1, sum)  
  
  iter1 <- 1; stop <- FALSE; 
  
  while (!stop) {
    f1 <- texp(x %*% B); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))
    score <- vec(t((t(x * wt)) %*% (y - f2)))
    temp <- getW(x, wt, B); A <- temp$A; W <- temp$W
    lstar <- apply(t(tln(cbind(f2, fref0)) * cbind(y, ys)) %*% wt, 2, sum) + 0.5  * tln(det(A))
    
    # CAN THIS BE INCORPORATED INTO getW?
    # WHY ARE WE ADDING 1e-06 vs. 1e-05?
    if (det(A) < 1e-10) diag(A) <- diag(A) + 1e-11
    #if (det(A) < 1e-05) diag(A) <- diag(A) + 1e-06
    Ainv <- solve(A)
    # Newton-Raphson method to get the next PLEs
    
    modScore <- score + 0.5 * (W %*% vec(Ainv))
    step <- Ainv %*% modScore
    iter2 <- 1; increase <- FALSE
    
    
    ## WHY 55??? - variance blows up when the number of iteration is too high (JS - based on conversation with SBB
    ## during a meeting on Sep 17, 2015)
    while (!increase & iter2 <= 55) {
      vecB.temp <- vecB + step
      B.temp <- matrix(data = vecB.temp, nrow = p, ncol = J, byrow = TRUE)
      f1 <- texp(x %*% B.temp); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))
      temp <- getW(x, wt, B.temp); A <- temp$A; W <- temp$W 
      lstar.temp <- apply(t(tln(cbind(f2, fref0)) * cbind(y, ys)) %*% wt, 2, sum) + 0.5 * tln(det(A))
      increase <- lstar.temp > lstar
      if (!increase) { step <- step/2; iter2 <- iter2 + 1 }
    }
    
    vecB <- vecB + step
    if (verbose) { cat(vecB, ": iter2 -", iter2,": step-",step, "\n", sep=" ") }
    B <- matrix(data = vecB, nrow = p, ncol = J, byrow=TRUE)
    stop <- sum(abs(step)) <= tol
    if (!stop) iter1 <- iter1 + 1
  }
  
  f1 <- texp(x %*% B); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))
  temp <- getV(x, wt, B); A <- temp$A; W <- temp$W; V <- temp$V
  lstar.max <- apply(t(tln(cbind(f2, fref0)) * cbind(y, ys)) %*% wt, 2, sum) + 0.5 * tln(det(A))
  #if (det(A) < 1e-05) diag(A) <- diag(A) + 1e-05
  if (det(A) < 1e-10) diag(A) <- diag(A) + 1e-11 #JS
  Ainv <- solve(A)
  Astar <- getAstar(p, J, A, W, V)
  if (det(Astar) < 1e-05) diag(Astar) <- diag(Astar) + 1e-06
  #if (det(Astar) < 1e-10) diag(Astar) <- diag(Astar) + 1e-11 #JS
  Astarinv <- solve(Astar)
  
  ## NOTE: when J>0, the columns/rows of the matrices A and Astar correspond to
  ## Intercept.1, Intercetp.2, ... Intercept.3, cov_11, cov_12,...cov_1J, ..., cov_p1, ..., cov_pJ
  fit <- list(B = B, Ainv = Ainv, Astarinv = Astarinv, lstar.max = lstar.max, iter1 = iter1, conv=stop)
  
  fit
  
}
