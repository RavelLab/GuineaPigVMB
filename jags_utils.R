##
## JAGS utilities
##
## Author: Pawel Gajer
## July 28, 2012
##

## return the mean and 95% CI
## modification of coda summary.mcmc()
## from output.R in coda/R
ci.mcmc <- function (object,short=TRUE,qs=NA, with.se=FALSE, n=NA)
{
    x <- mcmc.list(object)

    if ( length(qs)==1 && is.na(qs) ){
        if ( short && !with.se ){
            qs <- c(0.5, 0.025, 0.975)
            statnames <- c("Mean", "Median", "2.5%", "97.5%")
        } else if (short && with.se) {
            qs <- c(0.5, 0.025, 0.975)
            statnames <- c("Mean", "S.E.", "Median", "2.5%", "97.5%")
        } else if (!short && !with.se) {
            qs <- c(0.5, 0.025, 0.975, 0.25, 0.75)
            statnames <- c("Mean", "Median", "2.5%", "97.5%", "25%", "75%")
        } else {
            qs <- c(0.5, 0.025, 0.975, 0.25, 0.75)
            statnames <- c("Mean", "S.E.", "Median", "2.5%", "97.5%", "25%", "75%")
        }
    } else {
        statnames <- c("Mean", paste(100*qs[1:length(qs)],"%",sep=""))
    }

    if ( !is.numeric(qs) )
        stop(paste("qs has to be a numeric vector; qs:",qs))

    ##xtsvar <- matrix(nrow = nchain(x), ncol = nvar(x))
    if (is.matrix(x[[1]])) {
        xlong <- do.call("rbind", x)
    }
    else {
        xlong <- as.matrix(x)
    }

    ## xmean <- apply(xlong, 2, mean)
    ## varquant <- t(apply(xlong, 2, quantile, qs))

    nc <- ncol(xlong)

    xmean <- numeric(nc)
    for ( i in 1:nc )
        xmean[i] <- mean(xlong[,i],na.rm=TRUE)

    varquant <- matrix(nrow=nc,ncol=length(qs))
    for ( i in 1:nc )
        varquant[i,] <- quantile(xlong[,i],probs=qs,na.rm=TRUE)

    varstats <- matrix(nrow = nc, ncol = length(statnames),
                       dimnames = list(varnames(x), statnames))
    varstats[, 1] <- xmean

    if ( with.se )
    {
        xsd <- numeric(nc)
        for ( i in 1:nc )
             xsd[i] <- sd(xlong[,i])/sqrt(n)

        varstats[, 2] <- xsd
        varstats[, 3:(length(qs)+2)] <- varquant

    } else {
        varstats[, 2:(length(qs)+1)] <- varquant
    }

    varstats <- drop(varstats)
    return(varstats)
}
