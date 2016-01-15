test.LR <- function(x, y, wt, mt, B, h0, penalized, tol = 1e-6, verbose = FALSE, ...) {
   p <- nrow(B) # p covariates
   J <- ncol(B) # J categories

  
   ys <- 1 - apply(y,1,sum)

# under H_a
# is this the likelihood from fit$lstar.max?

   f1 <- texp(x %*% B); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))
   temp <- getW(x, wt, B); A <- temp$A; W <- temp$W
   la <- apply(t(tln(cbind(f2, fref0)) * cbind(y, ys)) %*% wt, 2, sum)
   if (penalized) la <- la + 0.5 * tln(det(A))

# calculate the likelihhod under H_0

   if (h0 == 1) {ntests <- p * J;   size <- J;  L <- 1}
   if (h0 == 2) {ntests <- p;       size <- 1;  L <- J}
   if (h0 == 3) {ntests <- p;       size <- 1;  L <- J - 1}
   if (h0 == 4) {ntests <- p;       size <- 1;  L <- J - 1}

   l0 <- pvalue <- matrix(data=0, nrow=ntests, ncol=1)
   conv.vec <- NULL
  ## var0.array - need to think about!
   if(h0 == 1) {
       B0.array <- var0.array <- NULL
   }

   ## ** MODIFIED

   if(h0 > 1) {
       B0.array <- var0.array <-  array( NA,dim = c(p,J,p) )
       #conv.array <- rep(NA,p)
   }
  
for (i in 1:ntests) {
    # How do we know and do we have to know which term is being tested?
      iter1 <- 0
      if (h0 == 1) {
          e <- matrix(data=0, nrow=p * J, ncol=L);
          e[i,1] <- 1
          ##cat("e matrix when h0 - ",h0,":\n",sep="")
          ##print(e)
      }

      if (h0 == 2) {
          C1 <- matrix(data=0, nrow=L, ncol=p * J);
          C1[1:L, (((i - 1) * J) + 1):(i * J)] <- diag(L);
          e <- t(C1)
          #cat("e matrix when h0 - ",h0,":\n",sep="")
          #print(e)
      }

      if (h0 == 3) {
         C1 <- matrix(data=0, nrow=L, ncol=p * J)
         C1[1:L, (((i - 1) * J) + 1):((i * J) - 1)] <- diag(L)
         for (m in 1:(J - 1)) {C1[m, ((i - 1) * J) + 1 + m] <- (-1)}
         e <- t(C1)
         #cat("e matrix when h0 - ",h0,":\n",sep="")
         #print(e)
     }
      
      if (h0 == 4) {
         C1 <- matrix(data=0, nrow=L, ncol=p * J)
         C1[1:L, (((i - 1) * J) + 1):((i * J) - 1)] <- diag(L)
         for (m in 1:(J - 1)) {C1[m, ((i - 1) * J) + 1 + m] <- (-m/(m + 1))}
         e <- t(C1)
         #cat("e matrix when h0 - ",h0,":\n",sep="")
         #print(e)
      }

      if(verbose){
          cat('e matrix for h0=',h0,'\n',sep="")
          print(e)
      }
      
      B0 <- matrix(data = 0, nrow = p, ncol = J, byrow = TRUE)
      converge <- FALSE

      niter1 <- 25 #original value was 20

      while (!converge & iter1 <= niter1) {
         f1 <- texp(x %*% B0); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))
         temp <- getW(x, wt, B0); A <- temp$A; W <- temp$W
         
         l0[i,1] <- apply(t(tln(cbind(f2, fref0)) * cbind(y, ys)) %*% wt, 2, sum)
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
         converge <- max(abs(delta)) < tol #convergence of beta
         if (verbose) {
             cat('l0-value:', l0, "\n",sep=" ")
             cat('max(abs(delta)), compared to \'tol\': ', max(abs(delta)), "\n",sep="")
         }
         iter1 <- iter1 + 1
     }
      if(h0 > 1) {
          B0.array[,,i] <- B0
          var0.array[,,i] <- diag(Ainv)
      }
      conv.vec <- c(conv.vec,converge)
      #cat("test - ",i," has been done: the inverse of information matrix is:\n",sep="")
      #print(Ainv)
  }

   conv.mat <- matrix(conv.vec, nrow = p, ncol = (ntests/p), byrow = TRUE)
   rownames(conv.mat) <- colnames(x)
   

   if(h0>1){
       dimnames(B0)[[1]] <- colnames(x)
   }

   if(h0 > 1) dimnames(B0.array) <- dimnames(var0.array)<- list(colnames(x), paste(paste0((attr(mt,"variables")[[2]]),collapse=" "),"=",colnames(y),sep=""), colnames(x))
   
   statistic <- matrix(data = LR <- 2 * (la - l0), nrow = p, ncol = size, byrow = TRUE)
   pvalue <- pchisq(statistic, L, lower.tail = FALSE)

   #conv=converge: added by JS, Sep 17, 2015
   list(l0 = l0, la = la, statistic = statistic, pvalue = pvalue, beta0.array = B0.array, var0.array = var0.array,
        Ainv0 = Ainv, conv = conv.mat)

}

