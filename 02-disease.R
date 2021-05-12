### Infectious disease example for AGHQ software paper ###

## Setup ----

library(aghq)
library(TMB)
precompile()
library(tmbstan)
library(parallel)
options(mc.cores = parallel::detectCores())

set.seed(8097968)

# change these to wherever you put the code and want the plots and data saved
# globalpath <- "/storage/phd/projects/aghq-softwarepaper/paper/"
# plotpath <- paste0(globalpath,"figures/disease/")
# datapath <- paste0(globalpath,"data/")
# use temp dirs
globalpath <- tempdir()
plotpath <- file.path(globalpath,"disease")
if (!dir.exists(plotpath)) dir.create(plotpath)
datapath <- plotpath

# the TMB template is part of the package. move it to a temp dir
# for compiling since this generates a bunch of new files
file.copy(system.file('extsrc/02_disease.cpp',package='aghq'),globalpath)

# Get the Tomato disease data
data("tswv", package = "EpiILMCT")

# Compile TMB template-- only need to do once
compile(file.path(globalpath,"02_disease.cpp"))
dyn.load(dynlib(file.path(globalpath,"02_disease")))


# Create the functions
dat <- tswv$tswvsir
dat$epidat <- dat$epidat[order(dat$epidat[ ,4]), ]

I <- dat$epidat[ ,4]
R <- dat$epidat[ ,2]
infected <- !is.infinite(I)

datlist <- list(
  D = as.matrix(dist(dat$location[dat$epidat[ ,1], ])),
  I = I,
  R = R,
  infected = as.numeric(infected[infected])
)

ff <- MakeADFun(data = datlist,
                parameters = list(theta1 = 0,theta2 = 0),
                DLL = "02_disease",
                ADreport = FALSE,
                silent = TRUE)

## Inference ----

# AGHQ
tm <- Sys.time()
quadmod <- aghq(ff,9,c(0,0),control = default_control(negate = TRUE))
aghqtime <- difftime(Sys.time(),tm,units='secs')

# STAN
stanmod <- tmbstan(
  ff,
  chains = 4,
  cores = 4,
  iter = 1e04,
  warmup = 1e03,
  init = quadmod$optresults$mode,
  seed = 124698,
  algorithm = "NUTS"
)
save(stanmod,file = file.path(datapath,"disease-stanmod-20210503.RData"))
# load(file.path(globalpath,"data/disease-stanmod-20210405.RData"))
## Summarize ----
pdf(file = file.path(plotpath,"stanmod-trace.pdf"),width=7,height=7)
traceplot(stanmod,window = c(9000,10000))
dev.off()

# Run time
max(get_elapsed_time(stanmod)[,2])
# Number of iterations
as.numeric(aghqtime) * stanmod@sim$iter / max(get_elapsed_time(stanmod)[,2])

stansamps <- as.data.frame(stanmod)
stansamps$alpha <- exp(stansamps$`par[1]`)
stansamps$beta <- exp(stansamps$`par[2]`)

posttrans <- list(totheta = log,fromtheta = exp)
quaddens <- compute_pdf_and_cdf(quadmod,posttrans)

quaddensalpha <- quaddens[[1]]
quaddensbeta <- quaddens[[2]]

# alpha
pdf(file.path(plotpath,"alpha-postplot.pdf"),width=7,height=7)
hist(stansamps$alpha,freq=FALSE,breaks=50,main = "",xlab = "",cex.lab=1.5,cex.axis = 1.5)
with(quaddensalpha,lines(transparam,pdf_transparam))
dev.off()

# beta
pdf(file.path(plotpath,"beta-postplot.pdf"),width=7,height=7)
hist(stansamps$beta,freq=FALSE,breaks=50,main = "",xlab = "",cex.lab=1.5,cex.axis = 1.5)
with(quaddensbeta,lines(transparam,pdf_transparam))
dev.off()

# summary stats

moms <- compute_moment(quadmod,function(x) exp(x))
getks <- function(x,y) {
  suppressWarnings(capture.output(ks <- ks.test(x,y)))
  unname(ks$statistic)
}
M <- nrow(stansamps)

quadsamps <- sample_marginal(quadmod,M)

summstats <- data.frame(
  stat = c('mean','sd','q2.5','q97.5','KS'),
  alphamcmc = c(
    mean(stansamps$alpha),
    sd(stansamps$alpha),
    quantile(stansamps$alpha,c(.025,.975)),
    NA
  ),
  alphaaghq = c(
    moms[1],
    sqrt( compute_moment(quadmod$normalized_posterior,function(x) ( (exp(x)[1] - moms[1])^2 )) ),
    exp(compute_quantiles(quadmod$marginals[[1]])),
    getks(stansamps$`par[1]`,quadsamps[[1]])
  ),
  
  betamcmc = c(
    mean(stansamps$beta),
    sd(stansamps$beta),
    quantile(stansamps$beta,c(.025,.975)),
    NA
  ),
  betaaghq = c(
    moms[2],
    sqrt( compute_moment(quadmod$normalized_posterior,function(x) ( (exp(x)[2] - moms[2])^2 )) ),
    exp(compute_quantiles(quadmod$marginals[[2]])),
    getks(stansamps$`par[2]`,quadsamps[[2]])
  )
)

readr::write_csv(summstats,file.path(plotpath,"summstattable.csv"))


# Joint moment
compute_moment(
  quadmod,
  function(x) exp(x)[1] * 2^(-exp(x)[2])
)
mean(stansamps$alpha * 2^(-stansamps$beta))


