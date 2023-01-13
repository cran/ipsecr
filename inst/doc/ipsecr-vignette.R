## ----schematic, fig.width=6.5, fig.height=6.5, message = FALSE, echo = FALSE----
library(ipsecr)
if (requireNamespace('plot3D')) {
    par(mfrow=c(2,2), mar=c(1,1,1,1), oma=c(1,1,2,1))
    oldplot <- plot3D.IP(ipsecrdemo)
    plot3D.IP(ipsecrdemo, box=2, oldplot)
    mtext(outer = TRUE, side = 3, c('Parameter space','Proxy space'), 
        adj = c(0.21,0.77))
} else warning ('install package plot3D to generate Fig. 2')

## ----setup, message = FALSE, results = 'hide'---------------------------------
library(ipsecr)
if (!require("spatstat")) warning ("install spatstat to run vignette code")
setNumThreads(2)   # adjust to number of available cores

## ----options, echo = FALSE----------------------------------------------------
runall <- FALSE
secr459 <- packageVersion('secr') >= '4.5.9'

## ----retrieve, eval = !runall, echo = FALSE-----------------------------------
# previously saved models ... see chunk 'saveall' at end
load(system.file("example", "fittedmodels.RData", package = "ipsecr"))

## ----ip.single.output---------------------------------------------------------
predict(ip.single)

## ----exampleproxy-------------------------------------------------------------
proxy.ms(captdata)

## ----compareproxy1------------------------------------------------------------
# secr function 'collate' works for both secr and ipsecr fits 
collate(ip.single, ip.single.1)[1,,,]

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

## ----nontargetdata------------------------------------------------------------
set.seed(123)
ch <- captdata
attr(ch, 'nontarget') <- (1-t(apply(ch,2:3,sum))) * (runif(500)>0.5)
summary(ch)$nontarget

## ----proxyfn2demo-------------------------------------------------------------
proxy.ms(ch)

## ----nontargetdemoresults-----------------------------------------------------
predict(ip.single.nontarget)

## ----fractional0, eval = runall-----------------------------------------------
#  ip.Fr <- ipsecr.fit(captdata, detectfn = 'HHN', details = list(factorial = 'fractional'))

## ----retrieveipFr, eval = !runall, echo = FALSE-------------------------------
# previously saved model ... see chunk 'saveall' at end
suppressMessages(load(system.file("example", "ip.Fr.RData", package = "ipsecr")))

## ----fractional1--------------------------------------------------------------
collate(ip.single, ip.Fr)[1,,,]
ip.single$proctime
ip.Fr$proctime

## ----fractional2, eval = FALSE, message = FALSE, warning = FALSE--------------
#  if (require('FrF2')) {
#    NP <- 3
#    boxsize <- rep(0.2,3)
#    design <- FrF2(2^(NP-1),NP, factor.names = c('D','lambda0','sigma'), ncenter = 2)
#    # recast factors as numeric
#    design <- sapply(design, function(x) as.numeric(as.character(x)))
#    design <- sweep(design, MAR=2, STATS = boxsize, FUN='*')
#    # apply to beta
#    beta <- log(c(5,0.2,25))
#    designbeta <- sweep(design, MAR=2, STATS=beta, FUN='+')
#    round(designbeta,3)
#  }

## ----extraparamdata, eval = secr459-------------------------------------------
grid <- make.grid(nx = 10, ny = 10, spacing = 20, detector = 'proximity')
msk <- make.mask(grid, buffer = 100)
set.seed(123)
pop <- sim.popn(D = 20, core = grid, buffer = 100, model2D = 'cluster', 
    details = list(mu = 5, hsigma = 1))
ch <- sim.capthist(grid, pop, detectfn = 14, detectpar = 
        list(lambda0 = 0.2, sigma = 20), noccasions = 5)
plot(ch, border = 20)

## ----extraparamfn-------------------------------------------------------------
# user function to simulate Thomas (Neyman-Scott) distribution of activity centres
# expect parameters mu and hsigma in list 'details$extraparam'
simclusteredpop <- function (mask, D, N, details) {
    secr::sim.popn(
        D = D[1], 
        core = mask, 
        buffer = 0, 
        Ndist = 'poisson',    # necessary for N-S cluster process
        model2D = 'cluster', 
        details = details$extraparam)
}

## ----extraparamdemo, eval = requireNamespace("spatstat") && runall------------
#  # extend the built-in proxy with clumping argument mu
#  # spatstat fits Thomas process parameters kappa and scale = hsigma^2
#  # mu is a model parameter derived from mu = D / kappa
#  clusterproxyT <- function (capthist, ...) {
#      pr <- ipsecr::proxy.ms(capthist)
#      pp <- spatstat.geom::as.ppp(secr::centroids(capthist),
#          W = as.numeric(apply(secr::traps(capthist),2,range)))
#      tfit <- spatstat.core::thomas.estK(pp)
#      c(pr, logmu = log(tfit$modelpar['mu']))
#  }
#  clusterfitT <- ipsecr.fit(ch, proxyfn = clusterproxyT, mask = msk,
#      detectfn = 'HHN', details = list(popmethod = simclusteredpop,
#      extraparam = list(mu = 5, hsigma = NA)), fixed = list(hsigma = 1))

## ----clusterresults-----------------------------------------------------------
predict(clusterfitT)

## ----clusterfitML, eval = runall----------------------------------------------
#  clusterfitML <- secr.fit(ch, mask = msk, detectfn = 'HHN', trace = FALSE)

## ----clustercompareML---------------------------------------------------------
predict(clusterfitML)

## ----adjustvard---------------------------------------------------------------
secr::adjustVarD(clusterfitML)

## ----saveall, echo = FALSE, eval = runall-------------------------------------
#  save(ip.single, ip.single.1, ip.single.nontarget, ipx, clusterfitT, clusterfitML,
#        file = '../inst/example/fittedmodels.RData')
#  save(ip.Fr, file = '../inst/example/ip.Fr.RData')
#  tools::resaveRdaFiles(paste0('../inst/example'),'xz')

