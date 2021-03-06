---
  "Modeling relative abundances of phylotypes using zero-inflated negative binomial models"

---
  
  
  setwd("../../data")

library(rjags)
source("../R/jags_utils.R") # load ci.mcmc()

raTbl <- read.csv("gp_relAb_tbl.csv",header=TRUE)
dim(raTbl) # 177 198

## meta-data and Chlamydia columns
colnames(raTbl)[1:9]


## Extracting absolute abundances of all phylotypes excluding Chlamydia
## (getting rid of meta-data columns and Chlamydia)
ct <- raTbl[,10:ncol(raTbl)]


## Identifying phylotypes present in at least 25% of all samples
bt <- apply(ct, 2, function(x) as.numeric(x>0))
tbt <- colSums(bt)
phProp <- 100 * tbt / nrow(bt)
selPh.25 <- names(phProp)[phProp>25]
length(selPh.25) # 67

## Creating day variable in raTbl that combines Cycle and Day
raTbl$day <- paste(raTbl$Cycle,".",raTbl$Day,sep="")
raTbl$day.f <- factor(raTbl$day, levels=c("1.2", "1.5", "1.8", "1.11",
                                          "1.14", "1.16", "2.2", "2.5", "2.8", "2.11", "2.14", "2.16"))
length(levels(raTbl$day.f)) # 12


## for each selected phylotype fit the following zero-inflated negative binomial model
##
## y  ~ treatment * day*cycle + (1|animal)
##

## directory to store plots of standardized residuals
zinbin.stRes.relAb.v4.dir <- "../pics/zinbin_relAb_v4_dir/"

## inverse of logit()
expit <- function(x) 1/(1 + exp(-x))

## function to run a model corresponding to the i-th element of selPh.25
jagsFn.v4 <- function(i)
{
  ph <- selPh.25[i]
  j.dat <- list(y=as.vector(ct[,ph]),
                log.total=log.total,
                treatment=as.integer(factor(raTbl$Experimental.group,
                                            levels=c("Infected","Non-infected","Mock"))),
                day=as.integer(raTbl$day.f),
                animal=as.integer(factor(raTbl$Animal)),
                nAnimals=length(unique(raTbl$Animal)))
  str(j.dat)
  j.dat
  
  j.m <- jags.model("../R/zinbin_withOffset_v4_ri.bug", j.dat, n.chains=2)
  nIter <- 10000
  update(j.m, nIter)
  j.pars <- c ("a","mur","mut","stRes")
  j.out <- coda.samples(j.m, j.pars, thin=10, n.iter=nIter)
  gelman.diag(j.out)
  j.ci <- ci.mcmc(j.out)
  
  ## plot standardized residual
  iRes <- grep("stRes",rownames(j.ci))
  resBar <- j.ci[iRes,1]
  res95 <- j.ci[iRes,3:4]
  n <- length(resBar)
  ylim <- range(as.vector(res95))
  if ( ylim[2] > 50 )
    ylim[2] <- 50
  
  pdf(paste("../",zinbin.stRes.relAb.v4.dir,ph,"_stRes_v4_plot.pdf",sep=""))
  plot(1:n, resBar, ylim=c(-2, ylim[2]), pch=19, xlab="Sample Index", ylab="Standardized Residual",cex=0.5, main="")
  for ( j in 1:n )
    segments(j,res95[j,1],j, res95[j,2], lwd=0.5)
  abline(h=2,col='red')
  abline(h=-2,col='red')
  dev.off()
  
  list(j.out=j.out, j.ci=j.ci)
}

## run jagsFn.v4() in parallel
## on my MacBook Pro it took about 1.5-2h to run the code on 8 cores
library(snowfall)
nProc <- 8        # number of cores/CPUs to use
ptm <- proc.time()
sfInit(parallel=TRUE, cpus=nProc, type="SOCK")
sfLibrary(rjags)
sfExport("selPh.25","ct","raTbl","log.total","ci.mcmc","zinbin.stRes.relAb.v4.dir")
results <- sfLapply(1:length(selPh.25), jagsFn.v4)
sfStop()
etm <- (proc.time() - ptm)
print(etm)

save(results, file="gp_results_relAb_v4.RData")
