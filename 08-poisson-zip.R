### Poisson Zero Inflated example ###

library(glmmTMB)
library(TMB)
precompile()
library(tmbstan)
library(aghq)
library(parallel)
options(mc.cores = parallel::detectCores())

# globalpath <- "/storage/phd/projects/aghq-softwarepaper/paper/"
# savepath <- paste0(globalpath,"data/")
# plotpath <- paste0(globalpath,"figures/poisson-zip/")

globalpath <- tempdir()
plotpath <- file.path(globalpath,"poisson-zip")
if (!dir.exists(plotpath)) dir.create(plotpath)

savestamp <- "-2021-05-06"

# MCMC code included at the end, but not part of paper
domcmc <- FALSE

## ZIP example from glmmTMB ##
data("Salamanders",package = "glmmTMB")
# First, fit the model for real
tm <- Sys.time()
zipmod <- glmmTMB(
  count ~ mined + (1|site),
  zi = ~mined,
  disp = ~DOY,
  data = Salamanders,
  family = nbinom2,
  doFit = TRUE
)
glmmTMBtime <- difftime(Sys.time(),tm,units = 'secs')
# glmmTMB creates all the necessary information
zipmodinfo <- glmmTMB(
  count ~ mined + (1|site),
  zi = ~mined,
  disp = ~DOY,
  data = Salamanders,
  family = nbinom2,
  doFit = FALSE
)
# Use this to create the TMB template
ff <- with(zipmodinfo,{
  MakeADFun(
    data = data.tmb,
    parameters = parameters,
    random = names(parameters)[grep('theta',names(parameters),invert = TRUE)],
    DLL = "glmmTMB",
    silent = TRUE
  )
})

# Fit using AGHQ
tm <- Sys.time()
zipquad <- marginal_laplace_tmb(ff,3,0)
zipquadsamps <- sample_marginal(zipquad,1e03)
zipsigmapdf <- compute_pdf_and_cdf(
  zipquad,list(totheta = log,fromtheta = exp),
  finegrid = seq(-10,0.2,length.out = 1000)
)
aghqtime <- difftime(Sys.time(),tm,units = 'secs')

# plot
pdf(file.path(plotpath,paste0("sigma-plot-aghq",savestamp,".pdf")),width = 7,height = 7)
with(zipsigmapdf[[1]],plot(transparam,pdf_transparam,type='l',main='',ylab='',xlab=expression(sigma),xaxt='n',xlim = c(0,1.2),ylim = c(0,3),cex.lab=1.5,cex.axis = 1.5))
title(ylab="Density",cex.lab=1.5,line=2.3)
axis(1,cex.axis=1.5,at=seq(0,1.2,by=.2))
abline(v = sqrt(as.numeric(summary(zipmod)$varcor$cond$site)),lty='dashed')
dev.off()

# confidence/credible intervals for random effects
rr <- as.data.frame(ranef(zipmod,condVar = TRUE))
Usamps <- zipquadsamps$samps[rownames(zipquadsamps$samps) == 'b', ]
Uest <- apply(Usamps,1,mean)
Ulower <- apply(Usamps,1,quantile,probs = .025)
Uupper <- apply(Usamps,1,quantile,probs = .975)

xx <- 1:nrow(rr)
offset <- .2
pdf(file.path(plotpath,paste0("interval-plot-aghq",savestamp,".pdf")),width = 7,height = 7)
with(rr,plot(xx-offset,condval,type = 'p',pch=20,xaxt='n',xlab='Site',ylab='',ylim = c(-1.7,1.7)))
points(xx+offset,Uest,type='p',pch=20)
with(rr,axis(1,at = xx,labels = as.character(grp)))
with(rr,arrows(x0=xx-offset,y0=condval-qnorm(.975)*condsd,y1=condval+qnorm(.975)*condsd,length=.05,angle=90,code=3))
arrows(x0=xx+offset,y0=Ulower,y1=Uupper,length=.05,angle=90,code=3)
dev.off()

if (domcmc) {
  # Not part of the paper. For interested readers.
  stanmod <- tmbstan(
    ff,
    chains = 4,
    cores = 4,
    iter = 1e04,
    warmup = 1e03,
    seed = 567984,
    algorithm = 'NUTS'
  )
  save(stanmod,file = file.path(savepath,paste0("stanmod",savestamp,".RData")))
  
  # plot
  standat <- as.data.frame(stanmod)
  
  pdf(file.path(plotpath,paste0("sigma-plot-mcmc",savestamp,".pdf")),width = 7,height = 7)
  hist(exp(standat$theta),freq=FALSE,breaks=100,main='',ylab='',xlab=expression(sigma),xaxt='n',cex.lab=1.5,cex.axis = 1.5)
  title(ylab="Density",cex.lab=1.5,line=2.3)
  axis(1,cex.axis=1.5,at=seq(0,1.6,by=.2))
  
  with(zipsigmapdf[[1]],lines(transparam,pdf_transparam))
  
  dev.off()
}


