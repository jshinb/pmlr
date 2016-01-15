#' @export
print.pmlr <- function(x,digits = max(3, getOption("digits") - 3), ...) {

    cat("Call: ")
    print(x$call)
    cat("\n")

    cat("Method: ",x$method,"\n",sep="")

    p <- dim(x$coefficients)[2]
    J <- dim(x$coefficients)[3]
    summary.table.colnames <- c("Estimate","Wald Std.Error","Chisq","P(>ChiSq)")
    output1 <- array(data = NA, dim = c(p,length(summary.table.colnames),J))
    dimnames(output1) <- list(dimnames(x$coefficients)[[2]],
                              summary.table.colnames,
                              dimnames(x$coefficients)[[3]])
    for (i in 1:J) {
        output1[,1,i] <- ifelse(is.infinite(x$separation[,,i]), x$separation[,,i],
                                format(round(x$coefficients[,,i],digits),justify="right"))
        output1[,2,i] <- ifelse(is.infinite(x$separation[,,i]), NA, format(round(sqrt(diag(x$var[,,i])),digits),justify="right"))
        ch = x$statistic[,,i]
        output1[,3,i] <- ifelse(is.infinite(x$separation[,,i]), NA,
                                format(round(ch, max(0, digits - log10(ch))),justify="right"))
        output1[,4,i] <- ifelse(is.infinite(x$separation[,,i]), NA,format.pval( x$pvalue[,,i], digits, eps = 0,justify="right"))
    }
    print(output1,quote=F,right=T)

    warnSep <- FALSE;
    for (i in 1:J) {
        if (is.infinite(as.vector(x$separation)[i]))
            {
                warnSep <- TRUE;
            }
    }

    if (warnSep)
        {
            cat("\n Warning: Separation has likely occurred for at least one coefficient. \n \n")
        }
}

