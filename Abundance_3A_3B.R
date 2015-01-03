
##
## Code for Figures 3A and 3B
##

library(rjags)
source("../R/jags_utils.R")   # load ci.mcmc(), plot.ci()

load("logTotal16S_v2.RData")  # logTotal16S - log total 16S rRNA abundance
load("raTbl.RData")           # raTbl - absolute abundnace table with metadata and row names
load("aTbl.RData")            # at, aTbl - - absolute abundnace table
load("aTbl2.RData")           # ct2, aTbl2, qt4
load("logOmpA_Nov4b.RData")   # logOmpA
load("aTbl.ompA_Nov4b.RData") # aTbl.ompA, ompA.qt3 - absolute abundnace table restricted to samples for which we have ompA data
load("phColorTbl.RData")      # phColorTbl - phylotype color table


at <- aTbl2[,9:(ncol(aTbl)-4)]

cs <- colSums(at)
o <- order(cs, decreasing=TRUE)
at <- at[,o]

pt2 <- t(apply(ct2+0.1, 1, function(x) x/sum(x) ))

total16S.n <- as.numeric( exp(logTotal16S) )

all(rownames(ct2)==names(logTotal16S)) # TRUE

at2 <- pt2
for ( i in 1:nrow(at2) )
    at2[i,] <- pt2[i,]*total16S.n[i]

log.at15 <- log10(at2[,1:15])

r1 <- c(3,8)

## --------------------------
## Code for Figure 3A
## --------------------------

ompA <- as.integer( exp(logOmpA) )

idx <- aTbl.ompA$Experimental.group=="Infected"

ompA.dat <- list(y=as.vector(ompA[idx]),
                  day=as.integer(aTbl.ompA$day.f[idx]),
                  animal=as.integer(factor(aTbl.ompA$animal[idx])),
                  nAnimals=length(unique(aTbl.ompA$animal[idx])))
str(ompA.dat)

ompA.m <- jags.model("nbin_ri_day.bug", ompA.dat, n.chains=3) ## no offset
nIter <- 20000
update(ompA.m, nIter)
ompA.pars <- c ("a","mur","stRes")
ompA.out <- coda.samples(ompA.m, ompA.pars, thin=20, n.iter=nIter)
## gelman.diag(ompA.out)
ompA.ci <- ci.mcmc(ompA.out)

## standardized residual
iRes <- grep("stRes",rownames(ompA.ci))
resBar <- ompA.ci[iRes,1]
res95 <- ompA.ci[iRes,3:4]
n <- length(resBar)
ylim <- range(as.vector(res95))
if ( ylim[2] > 50 )
    ylim[2] <- 50

plot(1:n, resBar, ylim=ylim, pch=19, xlab="Sample Index", ylab="Standardized Residual",cex=0.5, main="Total16S")
for ( j in 1:n )
    segments(j,res95[j,1],j, res95[j,2], lwd=0.5)
abline(h=2,col='red')
abline(h=-2,col='red')


i.a <- grep("a\\[",rownames(ompA.ci))
aBar <- ompA.ci[i.a,1]
aCI <- ompA.ci[i.a,3:4]

t.dat <- list(y=as.integer(exp(logTotal16S[rownames(raTbl)]-1)),
               treatment=as.integer(factor(aTbl$Experimental.group,
               levels=c("Infected","Non-infected","Mock"))),
               day=as.integer(aTbl$day.f),
               animal=as.integer(factor(raTbl$animal)),
               nAnimals=length(unique(raTbl$animal)))
str(t.dat)

t.m <- jags.model("nbin_ri.bug", t.dat, n.chains=2) ## no offset
nIter <- 10000
update(t.m, nIter)
t.pars <- c ("a")
t.out <- coda.samples(t.m, t.pars, thin=10, n.iter=nIter)
##gelman.diag(t.out)
t.ci <- ci.mcmc(t.out)

t.ci <- t.ci + 1 
## correction accounting for -1 in y=as.integer(exp(logTotal16S.v1-1)),

i.a1 <- grep("a\\[1",rownames(t.ci))
a1Bar <- t.ci[i.a1,1]
a1CI <- t.ci[i.a1,3:4]

i.a2 <- grep("a\\[2",rownames(t.ci))
a2Bar <- t.ci[i.a2,1]
a2CI <- t.ci[i.a2,3:4]

i.a3 <- grep("a\\[3",rownames(t.ci))
a3Bar <- t.ci[i.a3,1]
a3CI <- t.ci[i.a3,3:4]


r <- range(c(as.vector(aCI), as.vector(a1CI),as.vector(a2CI),as.vector(a3CI)))/log(10)
dy <- diff(r)/30
dx <- 0.3

##pdf("../pics/log10_ompA_total16S_v3.pdf", width=11, height=8.5)
op <- par(mar=c(4,4,0.5,0.5), mgp=c(2.75,0.6,0),tcl = -0.3)
plot(c(1,12), r, type='n', xlab="Days", ylab="",main="", axes=F)
axis(1, at=1:12, labels=c(2,5,8,11,14,16,2,5,8,11,14,16))
segments(1,2.3,6,2.3, xpd=NA, col='gray90')
mtext("1st estrous cycle", side=1, at=3.5, line=2.3)
segments(7,2.3,12,2.3, xpd=NA, col='gray90')
mtext("2nd estrous cycle", side=1, at=9.5, line=2.3)
aty <- axTicks(2)
labels <- sapply(aty,function(i) as.expression(bquote(10^ .(i))) )
axis(2,at=aty,labels=labels, las=1)
box()
## ompA
polygon(c(1:6, 6:1), c(aCI[1:6,1]/log(10), aCI[6:1,2]/log(10)), col=adjustcolor('red', alpha.f=0.05) , border=NA)
lines(1:6, aBar[1:6]/log(10), type="l", col='red', lwd=2)
polygon(c(7:12, 12:7), c(aCI[7:12,1]/log(10), aCI[12:7,2]/log(10)), col=adjustcolor('red', alpha.f=0.05) , border=NA)
lines(7:12, aBar[7:12]/log(10), type="l", col='red', lwd=2)
## infected
alpha <- 0.1
polygon(c(1:6, 6:1), c(a1CI[1:6,1]/log(10), a1CI[6:1,2]/log(10)), col=adjustcolor('orange', alpha.f=alpha) , border=NA)
lines(1:6, a1Bar[1:6]/log(10), type="l", col='orange', lwd=2)
polygon(c(7:12, 12:7), c(a1CI[7:12,1]/log(10), a1CI[12:7,2]/log(10)), col=adjustcolor('orange', alpha.f=alpha) , border=NA)
lines(7:12, a1Bar[7:12]/log(10), type="l", col='orange', lwd=2)
## mock-infected
polygon(c(1:6, 6:1), c(a3CI[1:6,1]/log(10), a3CI[6:1,2]/log(10)), col=adjustcolor('blue', alpha.f=alpha) , border=NA)
lines(1:6, a3Bar[1:6]/log(10), type="l", col='blue', lwd=2)
polygon(c(7:12, 12:7), c(a3CI[7:12,1]/log(10), a3CI[12:7,2]/log(10)), col=adjustcolor('blue', alpha.f=alpha) , border=NA)
lines(7:12, a3Bar[7:12]/log(10), type="l", col='blue', lwd=2)
## non-infected
alpha.g <- 0.15
polygon(c(1:6, 6:1), c(a2CI[1:6,1]/log(10), a2CI[6:1,2]/log(10)), col=adjustcolor('green', alpha.f=alpha.g) , border=NA)
lines(1:6, a2Bar[1:6]/log(10), type="l", col='green', lwd=2)
polygon(c(7:12, 12:7), c(a2CI[7:12,1]/log(10), a2CI[12:7,2]/log(10)), col=adjustcolor('green', alpha.f=alpha.g) , border=NA)
lines(7:12, a2Bar[7:12]/log(10), type="l", col='green', lwd=2)
legend("bottomleft",legend=c("C. caviae ompA","Infected guinea pigs","Non-infected guinea pigs","Mock-infected guinea pigs"), lwd=2, col=c("red","orange","green","blue"), inset=0.05, cex=1)
par(op)
##dev.off()


## --------------------------
## Code for Figure 3B
## --------------------------


##
## The top 15 or so phylotypes absolute abundance over each estrous
## cycle in 1. Infected, 2. Non-infected animals and 3. Mock-infected
##

##pdf("../pics/3B_3panels_v2.pdf", width=8.5, height=11)
op <- par(mar=c(0,4,0.5,0.5), mgp=c(2.75,0.6,0), tcl=-0.3)
lhei <- c(5,5,5,2)
layout(matrix(1:4, 4, 1, byrow = TRUE), heights = lhei,respect =F)
##layout.show(4)

plot(c(1,12), r1, type='n', xlab="Days", ylab="",main="", axes=F)
aty <- axTicks(2)
labels <- sapply(aty,function(i) as.expression(bquote(10^ .(i))) )
axis(2,at=aty,labels=labels, las=1)
box()
text(1,7.8, "A", cex=1.8)
## ompA
lines(1:6, aBar[1:6]/log(10), type="l", col='red', lwd=3)
lines(7:12, aBar[7:12]/log(10), type="l", col='red', lwd=3)
for ( i in 1:15 )
{
    x <- log.at15[,i]
    idx <- aTbl2$Experimental.group=="Infected"
    x <- x[idx]
    d <- aTbl2$day[idx]
    xx <- aggregate( x ~ d, FUN=mean)
    y <- xx$x
    names(y) <- xx$d
    y <- y[c("1.2", "1.5", "1.8", "1.11", "1.14", "1.16", "2.2", "2.5", "2.8", "2.11", "2.14", "2.16")]
    lines(1:6, y[1:6], col=phColorTbl[i])
    lines(7:12, y[7:12], col=phColorTbl[i])
}

## non-infected
op2 <- par(mar=c(0,4,0.5,0.5), mgp=c(2.75,0.6,0), tcl=-0.3)
plot(c(1,12), r1, type='n', xlab="", ylab="",main="", axes=F)
aty <- axTicks(2)
labels <- sapply(aty,function(i) as.expression(bquote(10^ .(i))) )
axis(2,at=aty,labels=labels, las=1)
box()
text(1,7.8, "B", cex=1.8)
for ( i in 1:15 )
{
    x <- log.at15[,i]
    idx <- aTbl2$Experimental.group=="Non-infected" ## [1] "Mock"         "Infected"     "Non-infected"
    x <- x[idx]
    d <- aTbl2$day[idx]
    xx <- aggregate( x ~ d, FUN=mean)
    y <- xx$x
    names(y) <- xx$d
    y <- y[c("1.2", "1.5", "1.8", "1.11", "1.14", "1.16", "2.2", "2.5", "2.8", "2.11", "2.14", "2.16")]
    lines(1:6, y[1:6], col=phColorTbl[i])
    lines(7:12, y[7:12], col=phColorTbl[i])
}
par(op2)

## mock-infected
op2 <- par(mar=c(4,4,0.5,0.5), mgp=c(2.75,0.6,0), tcl=-0.3, cex.lab=1.5)
plot(c(1,12), r1, type='n', xlab="Days", ylab="",main="", axes=F)
axis(1, at=1:12, labels=c(2,5,8,11,14,16,2,5,8,11,14,16))
mtext("1st estrous cycle", side=1, at=3.5, line=2.3)
mtext("2nd estrous cycle", side=1, at=9.5, line=2.3)
aty <- axTicks(2)
labels <- sapply(aty,function(i) as.expression(bquote(10^ .(i))) )
axis(2,at=aty,labels=labels, las=1)
box()
text(1,7.8, "C", cex=1.8)
for ( i in 1:15 )
{
    x <- log.at15[,i]
    idx <- aTbl2$Experimental.group=="Mock" ## [1] "Mock"         "Infected"     "Non-infected"
    x <- x[idx]
    d <- aTbl2$day[idx]
    xx <- aggregate( x ~ d, FUN=mean)
    y <- xx$x
    names(y) <- xx$d
    y <- y[c("1.2", "1.5", "1.8", "1.11", "1.14", "1.16", "2.2", "2.5", "2.8", "2.11", "2.14", "2.16")]
    lines(1:6, y[1:6], col=phColorTbl[i])
    lines(7:12, y[7:12], col=phColorTbl[i])
}
par(op2)
op <- par(mar=c(0,4,0.5,0.5), mgp=c(2.75,0.6,0), tcl=-0.3)
plot(c(1,12), r1, type='n', xlab="", ylab="",main="", axes=F)
legend("top", legend=c("C. caviae ompA", colnames(log.at15)), lty=1, col=c('red',phColorTbl[1:15]), inset=0.1, ncol=3, lwd=3)
par(op)
##dev.off()
