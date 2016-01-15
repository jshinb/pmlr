test.score <- function(x, y, wt, mt, B, h0, penalized, tol = 1e-6, verbose=FALSE,...) {
  p <- nrow(B) # p covariates
  J <- ncol(B) # J categories
  
  ys <- 1 - apply(y,1,sum)
  
  if (h0 == 1) 
  {
    B0.array <- var0.array <- NULL
    
    ntests <- p*J; size <- J;  L <- 1;
    l0 <- pvalue <- matrix(data=0, nrow=ntests, ncol=1)
    scoreT <- rep(NA, ntests)  
    AinvT <- matrix(0, nrow=ntests, ncol=ntests)
    
    for (i in 1:ntests) {
      iter1 <- 0
      e <- matrix(data=0, nrow=p * J, ncol=L); e[i,1] <- 1
      
      f1 <- texp(x %*% B); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))
      temp <- getW(x, wt, B); A <- temp$A; W <- temp$W
      la <- apply(t(tln(cbind(f2, fref0)) * cbind(y, ys)) %*% wt, 2, sum)
      if (penalized) la <- la + 0.5 * tln(det(A))
      
      # calculate the likelihhood under H_0
      B0 <- matrix(data = 0, nrow = p, ncol = J, byrow = TRUE)
      converge <- FALSE
      conv.vec <- NULL

      while (!converge & (iter1 <= 20)) {
        #cat("(!converge & (iter1 <= 20)):\n");print((!converge & (iter1 <= 20)))
        f1 <- texp(x %*% B0); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))
        temp <- getW(x, wt, B0); A <- temp$A; W <- temp$W
        l0[i,1] <- apply(t(tln(cbind(f2, fref0)) * cbind(y, ys)) %*% wt, 2, sum)
        if (penalized) l0[i,1] <- l0[i,1] + 0.5 * tln(det(A))
        
        if (det(A) < 0.01) diag(A) <- diag(A) + 0.01
        Ainv <- solve(A)
        
        score <- vec(t((t(x * wt)) %*% (y - f2)))
        if (penalized) score <- score + 0.5 * (W %*% vec(Ainv))
        #cat("# iter1 - ",iter1,"; score: \n", sep="");print(score)
        lambda <- -solve(t(e) %*% Ainv %*% e) %*% (t(e) %*% Ainv %*% score)
        delta <- Ainv %*% (score + e %*% lambda); delta <- matrix(data = delta, nrow = p, ncol = J, byrow = TRUE)
        
        iter2 <- 1; increase <- FALSE;
        while (!increase & iter2 < 5) {
          B0.temp <- B0 + delta
          f1 <- texp(x %*% B0.temp); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))
          temp <- getW(x, wt, B0.temp); A <- temp$A; W <- temp$W
          l.temp <- apply(t(tln(cbind(f2, fref0)) * cbind(y, ys)) %*% wt, 2, sum)
          if (penalized) l.temp <- l.temp + 0.5 * tln(det(A))
          increase <- (l.temp > l0[i,1])
          if (!increase) delta <- delta/2
          iter2 <- iter2 + 1
        }
        
        B0 <- B0 + delta
        
        f1 <- texp(x %*% B0); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))
        temp <- getW(x, wt, B0); A <- temp$A; W <- temp$W
        l0[i,1] <- apply(t(tln(cbind(f2, fref0)) * cbind(y, ys)) %*% wt, 2, sum)
        if (penalized) l0[i,1] <- l0[i,1] + 0.5 * tln(det(A))
        converge <- max(abs(delta)) < tol
        if (verbose) { cat(l0, "\n") }
        iter1 <- iter1 + 1
        #cat("(!converge & (iter1 <= 20)):\n"); print((!converge & (iter1 <= 20)))
      }
      scoreT[i] <- score[i];
      AinvT[i,i] <- Ainv[i,i]; # ORIGINAL
      conv.vec <- c(conv.vec,converge)
    }
    ## diag(AinvT) <- diag(Ainv); # REVISED BY JS: THIS IS WRONG!
    # faster computation
    #cat("\n scoreT: \n",sep=""); print(scoreT)
    #cat("\n AinvT: \n",sep=""); print(AinvT)
    statistic <- matrix(data = scoreRes <- AinvT %*% scoreT^2 , nrow = p, ncol = size, byrow = TRUE)
    pvalue <- pchisq(statistic, L, lower.tail = FALSE)
  }
  
  else if ((h0 ==2) || (h0 == 3) || (h0 == 4) ){
    B0.array <- var0.array <- array( NA,dim = c(p,J,p) )
    ntests <- p
    scoreRes <- rep(0, ntests);
    
    ###############
    l0 <- pvalue <- matrix(data=0, nrow=ntests, ncol=1)
    for (i in 1:ntests) {
      
      if (h0 == 2) {
        size <- 1;  L <- J;
        C1 <- matrix(data=0, nrow=L, ncol=p * J); C1[1:L, (((i - 1) * J) + 1):(i * J)] <- diag(L); e <- t(C1);
      }
      
      if (h0 == 3) {
        size <- 1;  L <- J - 1;
        C1 <- matrix(data=0, nrow=L, ncol=p * J)
        C1[1:L, (((i - 1) * J) + 1):((i * J) - 1)] <- diag(L)
        for (m in 1:(J - 1)) {C1[m, ((i - 1) * J) + 1 + m] <- (-1)}
        e <- t(C1)
      }
      
      if (h0 == 4) {
        size <- 1;  L <- J - 1
        C1 <- matrix(data=0, nrow=L, ncol=p * J)
        C1[1:L, (((i - 1) * J) + 1):((i * J) - 1)] <- diag(L)
        for (m in 1:(J - 1)) {C1[m, ((i - 1) * J) + 1 + m] <- (-m/(m + 1))}
        e <- t(C1)
      }
      
      iter1 <- 0
      
      f1 <- texp(x %*% B); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))
      temp <- getW(x, wt, B); A <- temp$A; W <- temp$W
      la <- apply(t(tln(cbind(f2, fref0)) * cbind(y, ys)) %*% wt, 2, sum)
      if (penalized) la <- la + 0.5 * tln(det(A))
      
      ## calculate the likelihhod under H_0
      B0 <- matrix(data = 0, nrow = p, ncol = J, byrow = TRUE)
      converge <- FALSE
      conv.vec <- NULL
      
      while (!converge & iter1 <= 20) {
        f1 <- texp(x %*% B0); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))
        temp <- getW(x, wt, B0); A <- temp$A; W <- temp$W
        l0[i,1] <- apply(t(tln(cbind(f2, fref0)) * cbind(y, ys)) %*% wt, 2, sum) ## check
        if (penalized) l0[i,1] <- l0[i,1] + 0.5 * tln(det(A))
        
        if (det(A) < 0.01) diag(A) <- diag(A) + 0.01
        Ainv <- solve(A)
        
        score <- vec(t((t(x * wt)) %*% (y - f2)))
        
        if (penalized) score <- score + 0.5 * (W %*% vec(Ainv))
        lambda <- -solve(t(e) %*% Ainv %*% e) %*% (t(e) %*% Ainv %*% score)
        
        delta <- Ainv %*% (score + e %*% lambda); delta <- matrix(data = delta, nrow = p, ncol = J, byrow = TRUE)
        
        iter2 <- 1; increase <- FALSE;
        while (!increase & iter2 < 5) {
          B0.temp <- B0 + delta
          f1 <- texp(x %*% B0.temp); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))
          temp <- getW(x, wt, B0.temp); A <- temp$A; W <- temp$W
          l.temp <- apply(t(tln(cbind(f2, fref0)) * cbind(y, ys)) %*% wt, 2, sum)
          if (penalized) l.temp <- l.temp + 0.5 * tln(det(A))
          increase <- (l.temp > l0[i,1])
          if (!increase) delta <- delta/2
          iter2 <- iter2 + 1
        }
        
        B0 <- B0 + delta
        f1 <- texp(x %*% B0); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))
        temp <- getW(x, wt, B0); A <- temp$A; W <- temp$W
        l0[i,1] <- apply(t(tln(cbind(f2, fref0)) * cbind(y, ys)) %*% wt, 2, sum)
        if (penalized) l0[i,1] <- l0[i,1] + 0.5 * tln(det(A))
        converge <- max(abs(delta)) < tol
        iter1 <- iter1 + 1
      }
      
      scoreRes[i] <- (t(score) %*% e) %*% (t(e) %*% Ainv %*% e) %*% (t(e) %*% score)
      B0.array[,,i] <- B0
      var0.array[,,i] <- diag(Ainv)
      # cat("h0 - ", h0, "\n", sep="")
      # scoreRes[i] <- t(score) %*% (Ainv) %*% score; # only correct when 
      # print( t(score) %*% (Ainv) %*% score )
      # print(scoreRes[i])
      conv.vec <- c(conv.vec,converge)
    }
    
    if(h0 > 1) dimnames(B0.array) <- dimnames(var0.array)<- list(colnames(x), paste(paste0((attr(mt,"variables")[[2]]),collapse=" "),"=",colnames(y),sep=""), colnames(x))
    
    
    ##<- t(as.matrix(scoreT[,,1])) %*% as.matrix(AinvT[,,1]) %*% as.matrix(scoreT[,,1])
    statistic <- matrix(data = scoreRes , nrow = p, ncol = size, byrow = TRUE)
    pvalue <- pchisq(statistic, L, lower.tail = FALSE)
    
  }
  else{stop("not available h0")}
  
  ## Added by JS (Oct. 28, 2015): conv=converge
  ##list(statistic = statistic, pvalue = pvalue, Ainv=Ainv, score=score)#Original
  ##list(statistic = statistic, pvalue = pvalue, Ainv=Ainv, score=scoreT, l0=l0, W=W, conv=converge)## revised
  conv.mat <- matrix(conv.vec, nrow = p, ncol = (ntests/p), byrow = TRUE)
  rownames(conv.mat) <- colnames(x)
  
  list(statistic = statistic, pvalue = pvalue, Ainv = Ainv, beta0.array = B0.array, var0.array = var0.array,
       score = scoreRes, l0 = l0, la = la, conv = conv.mat)## revised:W = W, 
}
