
------------------------------------------------------------------------------------------------
  "Comparison of evenness between Guinea pig vaginal, human vaginal and Guinea pig stool samples"
------------------------------------------------------------------------------------------------
  #Set the working directory to a folder containg all data and files needed to run the code 
setwd("../data/")
  #Install Packages needed
  
library("knitr",lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
library("memoise",lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
library("plyr",lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
library("scales",lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
library("markdown",lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")

eTbl <- read.csv("evenness.csv", header=T, row.names=1)

str(eTbl)

unique(eTbl$small_cluster)

(uqLoc <- unique(eTbl$location))

## Visual inspection of evenness quantile-quantile plots to see if they are
## normally distributed

e.1 <- eTbl$even[eTbl$small_cluster=="1"]
e.3A <- eTbl$even[eTbl$small_cluster=="3A"]
e.3B <- eTbl$even[eTbl$small_cluster=="3B"]
e.4A <- eTbl$even[eTbl$small_cluster=="4A"]
e.4B <- eTbl$even[eTbl$small_cluster=="4B"]
e.H <- eTbl$even[eTbl$small_cluster=="H"]
e.MK <- eTbl$even[eTbl$small_cluster=="MK"]
e.P <- eTbl$even[eTbl$small_cluster=="P"]
e.S <- eTbl$even[eTbl$small_cluster=="S"]



boxplot( eTbl$even ~ eTbl$small_cluster, col=2)

qqnorm(e.1)
qqline(e.1, col=2)

qqnorm(e.3A)
qqline(e.3A, col=2)

qqnorm(e.3B)
qqline(e.3B, col=2)

qqnorm(e.4A)
qqline(e.4A, col=2)

qqnorm(e.4B)
qqline(e.4B, col=2)

qqnorm(e.H)
qqline(e.H, col=2)

qqnorm(e.MK)
qqline(e.MK, col=2)

qqnorm(e.P)
qqline(e.P, col=2)

qqnorm(e.S)
qqline(e.S, col=2)




----------------------------------------------------------
  ## Setting up the normal model
  ----------------------------------------------------------
  
  
  
  ## Small_clusters
  
  m.dat <- list(y=eTbl$logit,
                gr=as.integer(factor(eTbl$small_cluster)),
                nGr=length(unique(eTbl$small_cluster)),
                subjID=as.integer(factor(eTbl$animal)),
                nSubj=length(unique(eTbl$animal)))


table(m.dat$gr, m.dat$subjID) 

## Only guinea pig has multiple samples per subject


## Creating guinea pig indicttor variable

gr3Idx <- ifelse(m.dat$gr==3, 1, 0)

table(m.dat$gr)



----------------------------------------------------------------------
  ## Small cluster 
  ----------------------------------------------------------------------
  
  table(eTbl$small_cluster)


m.dat <- list(y=eTbl$logit,
              gr=as.integer(factor(eTbl$small_cluster)),
              nGr=length(unique(eTbl$small_cluster)),
              subjID=as.integer(factor(eTbl$animal)),
              nSubj=length(unique(eTbl$animal)),
              subjIdx=gr3Idx)

str(m.dat)

```
---------------------------------------------------------------------------------------------------------
  ## Logit(evenness) vs categoris with animal is a random effect included in the analysis of the samples. 
  ## Only included when there are multiple samples per subject.
  ---------------------------------------------------------------------------------------------------------
  
  ## Running JAGS model
  
  m <- jags.model("norm_evenness_ri.bug", m.dat, n.chains=2)
nIter <- 10000

update(m, nIter)

m.pars <- c ("a","tau","stRes","y.rep","tau.b")
m.out <- coda.samples(m, m.pars, thin=10, n.iter=nIter)

## Checking convergence
g <- gelman.diag(m.out, multivariate=FALSE)$psrf
cbind(rownames(g)[g[,1]>1.1], g[g[,1]>1.1,1])


## Parameters mean/median values and their 95% CIs
m.ci <- ci.mcmc(m.out,short=T)


## --------------------
## y.rep
## --------------------

i.rep <- grep("y.rep", rownames(m.ci))
m.ci.yrep <- m.ci[i.rep,]

## Percentage of samples with response values outside of the 95% of predicted values

outSide <- 0
for ( i in 1:length(m.dat$y) )
{
  if ( m.dat$y[i] < m.ci.yrep[i,3] || m.dat$y[i] > m.ci.yrep[i,4] )
  {
    outSide <- outSide + 1
  }
}

perc.out <- 100*outSide/length(m.dat$y) # 5.367793

o <- order(m.dat$gr, m.dat$y)


op <- par(mar=c(3.8,3.5,0.1,0.1),mgp=c(2.75,0.6,0),tcl=-0.3)
plot(m.dat$y[o], las=1, ylab="logit(evenness)", type='n', axes=F, xlab="Small Cluster")
axis(2, las=1)
grPos <- c()

for ( i in 1:m.dat$nGr )
{
  grPos[i] <- mean(which(m.dat$gr[o]==i))
}
axis(1, at=grPos, labels=1:m.dat$nGr, tcl=0)
box()

for ( i in 1:nrow(m.ci.yrep) )
  segments(i,m.ci.yrep[o[i],3], i, m.ci.yrep[o[i],4], col=2)
points(m.dat$y[o], las=1, ylab="",pch=1)
text(275, 5, labels=paste("Percentage of response values outside of predicted 95% intervals: ",
                          formatC(perc.out, digits=2)), cex = 0.6)
par(op)

## In the figure 1 = CST1, 2 = CST 3A, 3 = CST 3B, 4 = CST 4A , 
## 5 = CST 4B , 6 = Healthy Guinea pig vaginal samples, 7 = Mock-infected 
## Guinea pig vaginal samples, 8 = Infected Guinea pig vaginal samples and 
## Guinea pig stool samples. 


------------------------------------------------------------
  ## Small cluster logit(evenness) values and their 95% CI
  ------------------------------------------------------------
  
  i.a <- grep("a\\[", rownames(m.ci))
m.ci.a <- m.ci[i.a,]


op <- par(mar=c(2,3.5,0.1,0.1),mgp=c(2.75,0.6,0),tcl=-0.3)
plot(c(0.5, m.dat$nGr+0.5), range(c(as.vector(m.ci.a), m.dat$y)), type='n',
     las=1, axes=FALSE, xlab="", ylab="logit(evenness)")
axis(2, las=1)
axis(1, at=1:m.dat$nGr, labels=levels(factor(eTbl$small_cluster)), tcl=0)
box()

## Plotting values
for ( i in 1:m.dat$nGr )
{
  points(jitter(rep(i, length(which(m.dat$gr==i))), amount=0.2), 
         m.dat$y[m.dat$gr==i])
  points(i, m.ci.a[i,1], col=2, cex=1, pch=19)
  segments(i, m.ci.a[i,3], i, m.ci.a[i,4], col=2, lwd=3)
}
par(op)


## Evenness scale

expit <- function(x) 1/(1 + exp(-x)) # inverse of logit()

y <- expit(m.dat$y) # evenness values

## Transforming m.ci.a to evenness scale
m.ci.e <- m.ci.a
for ( i in 1:nrow(m.ci.a) )
  m.ci.e[i,] <- expit(m.ci.a[i,])

op <- par(mar=c(2,3.5,0.1,0.1),mgp=c(2.75,0.6,0),tcl=-0.3)
plot(c(0.5, m.dat$nGr+0.5), range(c(as.vector(m.ci.e), y)), type='n',
     las=1, axes=FALSE, xlab="", ylab="Evenness")
axis(2, las=1)
axis(1, at=1:m.dat$nGr, labels=levels(factor(eTbl$small_cluster)), tcl=0)
box()

##Ploting values
for ( i in 1:m.dat$nGr )
{
  points(jitter(rep(i, length(which(m.dat$gr==i))), amount=0.2), y[m.dat$gr==i])
  points(i, m.ci.e[i,1], col=2, cex=1, pch=19)
  segments(i, m.ci.e[i,3], i, m.ci.e[i,4], col=2, lwd=3)
}
par(op)

```
-------------------------------------------------------------------------------
  ## Comparison of evenness between Guinea pig vaginal and human vaginal samples.
  -------------------------------------------------------------------------------
  
  
  eTbl <- read.csv("evenness_loc.csv", header=T, row.names=1)

str(eTbl)

(uqLoc <- unique(eTbl$location))

## Visual inspection of evenness quantile-quantile plots to see if they are
## normally distributed

e.vgp <- eTbl$even[eTbl$location=="gp"]
e.vh <- eTbl$even[eTbl$location=="h"]

boxplot( eTbl$even ~ eTbl$location, col=2)


qqnorm(e.vgp)
qqline(e.vgp, col=2)

qqnorm(e.vh)
qqline(e.vh, col=2)



table(m.dat$gr, m.dat$subjID) 

# Only guinea pig has multiple samples per subject

## Creating guinea pig indicttor variable

gr3Idx <- ifelse(m.dat$gr==3, 1, 0)

table(m.dat$gr)

--------------------------------
  ## Location guinea pig or human 
  --------------------------------
  
  
  table(eTbl$location)


m.dat <- list(y=eTbl$logit,
              gr=as.integer(factor(eTbl$location)),
              nGr=length(unique(eTbl$location)),
              subjID=as.integer(factor(eTbl$animal)),
              nSubj=length(unique(eTbl$animal)),
              subjIdx=gr3Idx)

str(m.dat)

------------------------------------------------------
  ## Running JAGS model
  ------------------------------------------------------
  
  
  m <- jags.model("norm_evenness_ri.bug", m.dat, n.chains=2)
nIter <- 10000
update(m, nIter)
m.pars <- c ("a","tau","stRes","y.rep","tau.b")
m.out <- coda.samples(m, m.pars, thin=10, n.iter=nIter)

## Checking convergence
g <- gelman.diag(m.out, multivariate=FALSE)$psrf
cbind(rownames(g)[g[,1]>1.1], g[g[,1]>1.1,1])

## Parameters mean/median values and their 95% CIs
m.ci <- ci.mcmc(m.out,short=T)


## --------------------
## y.rep
## --------------------

i.rep <- grep("y.rep", rownames(m.ci))
m.ci.yrep <- m.ci[i.rep,]

## Percentage of samples with response values outside of the 
## 95% of predicted values

outSide <- 0
for ( i in 1:length(m.dat$y) )
{
  if ( m.dat$y[i] < m.ci.yrep[i,3] || m.dat$y[i] > m.ci.yrep[i,4] )
  {
    outSide <- outSide + 1
  }
}

perc.out <- 100*outSide/length(m.dat$y) # 5.367793


o <- order(m.dat$gr, m.dat$y)


op <- par(mar=c(3.8,3.5,0.1,0.1),mgp=c(2.75,0.6,0),tcl=-0.3)
plot(m.dat$y[o], las=1, ylab="logit(evenness)", type='n', axes=F, 
     xlab="Location")
axis(2, las=1)
grPos <- c()

##Ploting values
for ( i in 1:m.dat$nGr )
{
  grPos[i] <- mean(which(m.dat$gr[o]==i))
}
axis(1, at=grPos, labels=1:m.dat$nGr, tcl=0)
box()

for ( i in 1:nrow(m.ci.yrep) )
  segments(i,m.ci.yrep[o[i],3], i, m.ci.yrep[o[i],4], col=2)
points(m.dat$y[o], las=1, ylab="",pch=1)

text(275, 5, labels=paste("Percentage of response values outside of predicted 95% intervals: ",
                          formatC(perc.out, digits=2)), cex = 0.6)
par(op)
## In the figure  1 = Guinea pig vaginal samples and 2 = Human vaginal samples.

```
--------------------------------------------------------------
  ## Location logit(evenness) values and their 95% CI
  --------------------------------------------------------------
  
  i.a <- grep("a\\[", rownames(m.ci))
m.ci.a <- m.ci[i.a,]


op <- par(mar=c(2,3.5,0.1,0.1),mgp=c(2.75,0.6,0),tcl=-0.3)

plot(c(0.5, m.dat$nGr+0.5), range(c(-1, 3)), type='n', las=1, 
     axes=FALSE, xlab="", ylab="logit(evenness)")
axis(2, las=1)
axis(1, at=1:m.dat$nGr, labels=levels(factor(eTbl$location)), tcl=0)
box()

##Ploting values
for ( i in 1:m.dat$nGr )
{
  points(jitter(rep(i, length(which(m.dat$gr==i))), amount=0.2), m.dat$y[m.dat$gr==i])
  points(i, m.ci.a[i,1], col=2, cex=1, pch=19)
  segments(i, m.ci.a[i,3], i, m.ci.a[i,4], col=2, lwd=3)
}
par(op)

## ----------------------
## Evenness scale
## ----------------------
expit <- function(x) 1/(1 + exp(-x)) # inverse of logit()

y <- expit(m.dat$y) # evenness values

## Transforming m.ci.a to evenness scale
m.ci.e <- m.ci.a
for ( i in 1:nrow(m.ci.a) )
  m.ci.e[i,] <- expit(m.ci.a[i,])

op <- par(mar=c(2,3.5,0.1,0.1),mgp=c(2.75,0.6,0),tcl=-0.3)
plot(c(0.5, m.dat$nGr+0.5), range(c(as.vector(m.ci.e), y)), type='n', 
     las=1, axes=FALSE, xlab="", ylab="Evenness")

axis(2, las=1)
axis(1, at=1:m.dat$nGr, labels=levels(factor(eTbl$location)), tcl=0)
box()

for ( i in 1:m.dat$nGr )
{
  points(jitter(rep(i, length(which(m.dat$gr==i))), amount=0.2), y[m.dat$gr==i])
  points(i, m.ci.e[i,1], col=2, cex=1, pch=19)
  segments(i, m.ci.e[i,3], i, m.ci.e[i,4], col=2, lwd=3)
}
par(op)
```
-----------------------------------------------------------------------------
  ## Comparison of evenness between Guinea pig vaginal samples during the 
  ## hight of infection and after in health and infected exprimental groups.
  -----------------------------------------------------------------------------
  
  eTbl <- read.csv("evenness_gp4.csv", header=T, row.names=1)

str(eTbl)


unique(eTbl$grouping)


## Visual inspection of evenness quantile-quantile 
## plots to see if they are normally distributed

e.P1 <- eTbl$even[eTbl$grouping=="P_1"]
e.P2 <- eTbl$even[eTbl$grouping=="P_2"]
e.H1 <- eTbl$even[eTbl$grouping=="H_1"]
e.H2 <- eTbl$even[eTbl$grouping=="H_2"]

## P1 = Infected animal samples from cycle 1 days 2 , 5 and 8
## P2 = Infected animal samples from cycle 2 days 2 , 5 and 8
## P1 = Non-infected animal samples from cycle 1 days 2 , 5 and 8
## P2 = Non-infected animal samples from cycle 2 days 2 , 5 and 8

boxplot( eTbl$even ~ eTbl$grouping, col=2)

qqnorm(e.P1)
qqline(e.P1, col=2)

qqnorm(e.P2)
qqline(e.P2, col=2)



qqnorm(e.H1)
qqline(e.H1, col=2)

qqnorm(e.H2)
qqline(e.H2, col=2)

----------------------------------------------------
  ## Setting up normal model
  ----------------------------------------------------
  
  
  ## ------------------
## Guinea pig groups
## ------------------

table(eTbl$grouping)


m.dat <- list(y=eTbl$logit,
              gr=as.integer(factor(eTbl$grouping)),
              nGr=length(unique(eTbl$grouping)),
              subjID=as.integer(factor(eTbl$animal)),
              nSubj=length(unique(eTbl$animal)),
              subjIdx=gr3Idx)

str(m.dat)

## Creating guinea pig indicttor variable

gr3Idx <- ifelse(m.dat$gr==3, 1, 0)

table(m.dat$gr)


-------------------------------------------------
  ## Running JAGS model
  -------------------------------------------------
  
  
  m <- jags.model("norm_evenness_ri.bug", m.dat, n.chains=2)
nIter <- 10000
update(m, nIter)
m.pars <- c ("a","tau","stRes","y.rep","tau.b")
m.out <- coda.samples(m, m.pars, thin=10, n.iter=nIter)

## Checking convergence
g <- gelman.diag(m.out, multivariate=FALSE)$psrf
cbind(rownames(g)[g[,1]>1.1], g[g[,1]>1.1,1])


## Parameters mean/median values and their 95% CIs
m.ci <- ci.mcmc(m.out,short=T)


## --------------------
## y.rep
## --------------------
i.rep <- grep("y.rep", rownames(m.ci))
m.ci.yrep <- m.ci[i.rep,]

## Percentage of samples with response values outside of the 95% of predicted values
outSide <- 0
for ( i in 1:length(m.dat$y) )
{
  if ( m.dat$y[i] < m.ci.yrep[i,3] || m.dat$y[i] > m.ci.yrep[i,4] )
  {
    outSide <- outSide + 1
  }
}

perc.out <- 100*outSide/length(m.dat$y) # 5.367793


o <- order(m.dat$gr, m.dat$y)
op <- par(mar=c(3.8,3.5,0.1,0.1),mgp=c(2.75,0.6,0),tcl=-0.3)
plot(m.dat$y[o], las=0, ylab="logit(evenness)", type='n', axes=F, 
     xlab="Guinea pig group")
axis(2, las=2)
grPos <- c()

for ( i in 1:m.dat$nGr )
{
  grPos[i] <- mean(which(m.dat$gr[o]==i))
}

axis(1, at=grPos, labels=1:m.dat$nGr, tcl=0)
box()

for ( i in 1:nrow(m.ci.yrep) )
  segments(i,m.ci.yrep[o[i],3], i, m.ci.yrep[o[i],4], col=2)
points(m.dat$y[o], las=1, ylab="",pch=1)
text(200, 0, labels=paste("Percentage of response values outside of predicted 95% intervals: ",
                          formatC(perc.out, digits=2)), cex = 0.8)
par(op)
## In the figure 1 = Healthy cycle 1 days 2 ,5 and 8, 
## 2 = Healthy infected cycle 2 days 2 ,5 and 8, 3 = Infected cycle 1 days 2 ,5 and 8 and 4 =  Infected cycle 2 days 2 ,5 and 8.



----------------------------------------------------------
  ## Grouping logit(evenness) values and their 95% CI
  ----------------------------------------------------------
  
  i.a <- grep("a\\[", rownames(m.ci))
m.ci.a <- m.ci[i.a,]


op <- par(mar=c(2,3.5,0.1,0.1),mgp=c(2.75,0.6,0),tcl=-0.3)
plot(c(0.5, m.dat$nGr+0.5), range(c(-1, 3)), type='n', las=1, axes=FALSE, 
     xlab="", ylab="logit(evenness)")
axis(2, las=1)
axis(1, at=1:m.dat$nGr, labels=levels(factor(eTbl$grouping)), tcl=0)
box()

##Ploting values
for ( i in 1:m.dat$nGr )
{
  points(jitter(rep(i, length(which(m.dat$gr==i))), amount=0.2),
         m.dat$y[m.dat$gr==i])
  points(i, m.ci.a[i,1], col=2, cex=1, pch=19)
  segments(i, m.ci.a[i,3], i, m.ci.a[i,4], col=2, lwd=3)
}
par(op)


## Evenness scale

expit <- function(x) 1/(1 + exp(-x)) # inverse of logit()

y <- expit(m.dat$y) # evenness values

## Transforming m.ci.a to evenness scale
m.ci.e <- m.ci.a
for ( i in 1:nrow(m.ci.a) )
  m.ci.e[i,] <- expit(m.ci.a[i,])

op <- par(mar=c(2,3.5,0.1,0.1),mgp=c(2.75,0.6,0),tcl=-0.3)
plot(c(0.5, m.dat$nGr+0.5), range(c(as.vector(m.ci.e), y)), 
     type='n', las=1, axes=FALSE, xlab="", ylab="Evenness")
axis(2, las=1)
axis(1, at=1:m.dat$nGr, labels=levels(factor(eTbl$grouping)), tcl=0)
box()

##Ploting values
for ( i in 1:m.dat$nGr )
{
  points(jitter(rep(i, length(which(m.dat$gr==i))), amount=0.2), 
         y[m.dat$gr==i])
  points(i, m.ci.e[i,1], col=2, cex=1, pch=19)
  segments(i, m.ci.e[i,3], i, m.ci.e[i,4], col=2, lwd=3)
}
par(op)
