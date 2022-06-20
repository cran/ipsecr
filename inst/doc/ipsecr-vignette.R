## ----setup, message = FALSE, results = 'hide'---------------------------------
library(ipsecr)
setNumThreads(2)   # adjust to number of available cores

## ----options, echo = FALSE----------------------------------------------------
runall <- FALSE

## ----retrieve, eval = !runall, echo = FALSE-----------------------------------
# previously saved models ... see chunk 'saveall' at end
load(system.file("example", "fittedmodels.RData", package = "ipsecr"))

## ----ip.single.output---------------------------------------------------------
predict(ip.single)

## ----exampleproxy-------------------------------------------------------------
proxy.ms(captdata)

## ----compareproxy1------------------------------------------------------------
# secr function 'collate' works for both secr and ipsecr fits 
collate(ip.single, ip.single.0)[1,,,]

## ----nontargetdata------------------------------------------------------------
set.seed(123)
ch <- captdata
attr(ch, 'nontarget') <- (1-t(apply(ch,2:3,sum))) * (runif(500)>0.5)
summary(ch)$nontarget

## ----proxyfn2demo-------------------------------------------------------------
proxy.ms(ch)

## ----nontargetdemoresults-----------------------------------------------------
predict(ip.single.nontarget)

## ----simch, eval = TRUE-------------------------------------------------------
tr <- traps(captdata)
mask <- make.mask(tr)
covariates(mask) <- data.frame(D = (mask$x-265)/20)  # for sim.pop
set.seed(1237)
pop <- sim.popn(D = 'D', core = mask, model2D = 'IHP', buffer = 100)
ch <- sim.capthist(tr, popn = pop, detectfn = 'HHN', noccasions = 5, 
  detectpar = list(lambda0 = 0.2, sigma = 25))
# show east-west trend
table(tr[trap(ch),'x'])

## ----Dsurfaceoutput, eval = TRUE----------------------------------------------
coef(ipx)
predict(ipx)
plot(predictDsurface(ipx))
plot(tr, add = TRUE)
plotMaskEdge(ipx$mask, add=T)

## ----Dsurfacetrend, fig.width = 6, fig.height = 4-----------------------------
oldpar <- par(mar = c(4,6,4,4))
# model refers to scaled values of x, so repeat here
m <- mean(traps(ch)$x); s <- sd(traps(ch)$x)
xval <- seq(270,730,10)
xvals <- scale(xval, m, s)
pred <- predict(ipx, newdata = data.frame(x = xvals))
plot(0,0, type='n', xlim= c(270,730), ylim = c(0,40), xlab = 'x', ylab = 'Density')
lines(xval, sapply(pred, '[', 'D','estimate'), col = 'red', lwd=2)
abline(-265/20,0.05) # true linear trend
rug(unique(tr$x))    # trap locations
par(oldpar)

## ----fractional1--------------------------------------------------------------
collate(ip.single, ip.Fr)[1,,,]
ip.single$proctime
ip.Fr$proctime

## ----fractional2, eval = TRUE, message = FALSE--------------------------------
library(FrF2, quietly = TRUE)
NP <- 3
boxsize <- rep(0.2,3)
design <- FrF2(2^(NP-1),NP, factor.names = c('D','lambda0','sigma'), ncenter = 2)
data.frame(design)
# recast factors as numeric
design <- sapply(design, function(x) as.numeric(as.character(x)))
design <- sweep(design, MAR=2, STATS = boxsize, FUN='*')
design
# apply to beta
beta <- log(c(5,0.2,25))
designbeta <- sweep(design, MAR=2, STATS=beta, FUN='+')
designbeta 
