### IMPLEMENTING APPROXIMATE BAYESIAN INFERENCE USINGN ADAPTIVE QUADRATURE: THE AGHQ PACKAGE ###
# This script reproduces the results in that paper
# See https://github.com/awstringer1/aghq-software-paper-code for README
# Alex Stringer
# 2021/05/11
#
# IMPORTANT NOTES:
#   - Use CRAN version 2.1 of aghq.
#   - The script depends on external resources only through the following two lines (1024, 1025):
#         cameroonBorderLL = raster::getData("GADM", country=c('CMR'), level=2)
#         nigeriaBorderLL = raster::getData("GADM", country=c('NGA'), level=2)
#   I have observed a small number of times that the server where this pulls data from has been down.
#   If this happens, it usually comes back up in a day or so.
#   - The compilation of the TMB templates creates a ton of output from the compiler. Each
#   template can take a few minutes at least to compile.
#
# TIMING:
#   - The most recent testing on a private Ubuntu server yielded run times as follows:
#      - Example 2: <10s
#      - Example 4.1: 6m24s
#      - Example 4.2: 5m43s
#      - Example 5.1: 47m50s
#      - Example 5.2: 10m44s
#      - Example 6.1: 11s
#
#   - The only operation that takes a long time is the MCMC within Example 5.1. 
#     To skip Example 5.1, set dofast = TRUE
dofast = FALSE
#
# VARIABLES TO SET:
# Example 5.2:
# Set the resolution for the spatial interpolations.
# The results shown in the paper use:
# reslist <- list(nrow = 200,ncol = 400)
# but this takes a couple hours. Here I set:
reslist <- list(nrow = 50,ncol = 100)
# which is faster but produces a blockier map.
# So some plots show up properly in the spinning:
par(mar = c(5,5,0,0))

## Check for missing packages ----
# The astro example requires ipoptr, and IPOPT (https://coin-or.github.io/Ipopt/INSTALL.html).
# This is a laborious installation.
# If you want to not run the astro example set doastro = FALSE
doastro <- TRUE
if (!('ipoptr' %in% installed.packages()[ ,'Package'])) {
  warning("No installation of ipoptr found. Skipping the astro example.\n")
  doastro <- FALSE
}

# Obtain the disease example data
# If can't obtain, then skip that example
dodisease <- TRUE
if (dodisease) {
  # Get the Tomato disease data
  if (require('EpiILMCT')) {
    data("tswv", package = "EpiILMCT")
  } else {
    download.file("https://github.com/waleedalmutiry/EpiILMCT/blob/master/data/tswv.RData?raw=true",destfile = file.path(globalpath,"tomato.RData"))
    load(file.path(globalpath,'tomato.RData'))
  }
  if (!exists('tswv')) {
    warning('You asked to do the disease example but the data cannot be obtained. This is either because you do not have the EpiILMCT package installed or because something is preventing downloading the data from github on https using download.file(). In any event, skipping the disease example.\n')
    dodisease <- FALSE
  }
}

# See if INLA is installed.
# If INLA is not installed then cannot compare to it in the loaloa example so skip it.
doloaloa <- TRUE
if (!('INLA' %in% installed.packages()[ ,'Package'])) {
  warning("No installation of INLA found. Skipping the loaloa example.\n")
  doloaloa <- FALSE
}
# But, don't do loaloa if the dofast flag is TRUE
if (dofast) doloaloa <- FALSE


## Install and load packages ----

# Check if installed
neededpackages <- c(
  'aghq',
  'TMB',
  'tmbstan',
  'parallel',
  'glmmTMB',
  'geostatsp',
  'PrevMap',
  'geoR',
  'trustOptim',
  'numDeriv'
)

missingpackages <- setdiff(neededpackages,installed.packages()[ ,'Package'])
if (length(missingpackages)) stop(paste0("Missing the following packages, please install: ",missingpackages,"\n"))

# Code for users to run interactively, if they like
if (FALSE) {
  install.packages(missingpackages)
}

library(aghq)
library(TMB)
precompile()
library(tmbstan)
library(parallel)
options(mc.cores = parallel::detectCores())
library(Matrix)
library(glmmTMB)
library(geostatsp)
library(PrevMap)
library(geoR)

## Set up the directory structure ----
# Each example gets its own directory within the tempdir()
globalpath <- tempdir()

## Example 2: Basic Use ----

plotpath <- file.path(globalpath,"basic-use")
if (!dir.exists(plotpath)) dir.create(plotpath)

set.seed(84343124)
y <- rpois(10,5) # True lambda = 5, n = 10


# Define the log posterior, log(pi(theta,y)) here
logpithetay <- function(theta,y) {
  sum(y) * theta - (length(y) + 1) * exp(theta) - sum(lgamma(y+1)) + theta
}
objfunc <- function(x) logpithetay(x,y)
objfuncgrad <- function(x) numDeriv::grad(objfunc,x)
objfunchess <- function(x) numDeriv::hessian(objfunc,x)
# Now create the list to pass to aghq()
funlist <- list(
  fn = objfunc,
  gr = objfuncgrad,
  he = objfunchess
)

# AGHQ with k = 3
# Use theta = 0 as a starting value
thequadrature <- aghq::aghq(ff = funlist,k = 3,startingvalue = 0)

summary(thequadrature)
plot(thequadrature)

# The posterior
thequadrature$normalized_posterior$nodesandweights
# The log normalization constant:
thequadrature$normalized_posterior$lognormconst
# Compare to the truth: 
lgamma(1 + sum(y)) - (1 + sum(y)) * log(length(y) + 1) - sum(lgamma(y+1))
# Quite accurate with only n = 10 and k = 3; this example is very simple.
# The mode found by the optimization:
thequadrature$optresults$mode
# The true mode:
log((sum(y) + 1)/(length(y) + 1))

# Compute the pdf for theta
transformation <- list(totheta = log,fromtheta = exp)
pdfwithlambda <- compute_pdf_and_cdf(
  thequadrature,
  transformation = transformation
)[[1]]
head(pdfwithlambda,n = 2)
lambdapostsamps <- sample_marginal(thequadrature,1e04,transformation = transformation)[[1]]
# Plot along with the true posterior
# pdf(file = file.path(plotpath,'lambda-post-plot.pdf'))
with(pdfwithlambda,{
  hist(lambdapostsamps,breaks = 50,freq = FALSE,main = "",xlab = expression(lambda))
  lines(transparam,pdf_transparam)
  lines(transparam,dgamma(transparam,1+sum(y),1+length(y)),lty='dashed')
})
# dev.off()

# Check if the posterior integrates to 1, by computing the "moment" of "1":
compute_moment(thequadrature$normalized_posterior,
               ff = function(x) 1)
# Posterior mean for theta:
compute_moment(thequadrature$normalized_posterior,
               ff = function(x) x)
# Posterior mean for lambda = exp(theta)
compute_moment(thequadrature$normalized_posterior,
               ff = function(x) exp(x))
# Compare to the truth:
(sum(y) + 1)/(length(y) + 1)

# Quantiles
compute_quantiles(
  thequadrature,
  q = c(.01,.25,.50,.75,.99),
  transformation = transformation
)[[1]]
# The truth:
qgamma(c(.01,.25,.50,.75,.99),1+sum(y),1+length(y))

#### END EXAMPLE 2 ####

## Example 4.1: Infectious Disease Modelling ----

if (dodisease) {
  
  set.seed(8097968)
  
  # use temp dirs
  plotpath <- file.path(globalpath,"disease")
  if (!dir.exists(plotpath)) dir.create(plotpath)
  datapath <- plotpath
  
  # the TMB template is part of the package. move it to a temp dir
  # for compiling since this generates a bunch of new files
  file.copy(system.file('extsrc/02_disease.cpp',package='aghq'),globalpath)
  
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
  # save(stanmod,file = file.path(datapath,"disease-stanmod-20210503.RData"))
  # load(file.path(globalpath,"data/disease-stanmod-20210405.RData"))
  ## Summarize ----
  # pdf(file = file.path(plotpath,"stanmod-trace.pdf"),width=7,height=7)
  # traceplot(stanmod,window = c(9000,10000))
  # dev.off()
  
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
  # pdf(file.path(plotpath,"alpha-postplot.pdf"),width=7,height=7)
  hist(stansamps$alpha,freq=FALSE,breaks=50,main = "",xlab = "",cex.lab=1.5,cex.axis = 1.5)
  with(quaddensalpha,lines(transparam,pdf_transparam))
  # dev.off()
  
  # beta
  # pdf(file.path(plotpath,"beta-postplot.pdf"),width=7,height=7)
  hist(stansamps$beta,freq=FALSE,breaks=50,main = "",xlab = "",cex.lab=1.5,cex.axis = 1.5)
  with(quaddensbeta,lines(transparam,pdf_transparam))
  # dev.off()
  
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
  
  # readr::write_csv(summstats,file.path(plotpath,"summstattable.csv"))
  knitr::kable(summstats,digits = 3)
  
  # Joint moment
  compute_moment(
    quadmod,
    function(x) exp(x)[1] * 2^(-exp(x)[2])
  )
  mean(stansamps$alpha * 2^(-stansamps$beta))
  
  #### END EXAMPLE 4.1 ####
}

if (doastro) {
## Example 4.2: Galactic Mass Estimation ----
  
  set.seed(563478)
  plotpath <- file.path(globalpath,"astro")
  if (!dir.exists(plotpath)) dir.create(plotpath)
  # the TMB template is part of the package. move it to a temp dir
  # for compiling since this generates a bunch of new files
  file.copy(system.file('extsrc/01_astro.cpp',package='aghq'),globalpath)
  
  # Compile TMB template-- only need to do once
  compile(file.path(globalpath,"01_astro.cpp"))
  
  data("gcdatalist",package = 'aghq')
  dyn.load(dynlib(file.path(globalpath,"01_astro")))
  
  
  # Function and its derivatives
  ff <- MakeADFun(data = gcdatalist,
                  parameters = list(theta1 = 0,
                                    theta2 = 0,
                                    theta3 = 0,
                                    theta4 = 0
                  ),
                  DLL = "01_astro",
                  ADreport = FALSE,
                  silent = TRUE)
  # Nonlinear constraints and their jacobian
  Es <- MakeADFun(data = gcdatalist,
                  parameters = list(theta1 = 0,
                                    theta2 = 0,
                                    theta3 = 0,
                                    theta4 = 0
                  ),
                  DLL = "01_astro",
                  ADreport = TRUE,
                  silent = TRUE)
  ## Parameter transformations ##
  parambounds <- list(
    Psi0 = c(1,200),
    gamma = c(.3,.7),
    alpha = c(3.0,3.7),
    beta = c(-.5,1)
  )
  
  get_psi0 <- function(theta) {
    # theta = -log( (Psi0 - 1) / (200 - 1) )
    (parambounds$Psi0[2] - parambounds$Psi0[1]) * 
      exp(-exp(theta)) + parambounds$Psi0[1]
  }
  get_theta1 <- function(Psi0) log(
    -log( 
      (Psi0 - parambounds$Psi0[1]) / (parambounds$Psi0[2] - parambounds$Psi0[1]) 
    )
  )
  
  get_gamma <- function(theta) {
    # theta = -log( (gamma - .3) / (.7 - .3) )
    (parambounds$gamma[2] - parambounds$gamma[1]) * 
      exp(-exp(theta)) + parambounds$gamma[1]
  }
  get_theta2 <- function(gamma) log(
    -log( 
      (gamma - parambounds$gamma[1]) / (parambounds$gamma[2] - parambounds$gamma[1]) 
    )
  )
  
  get_alpha <- function(theta) {
    # theta = log(alpha - 3)
    exp(theta) + parambounds$alpha[1]
  }
  get_theta3 <- function(alpha) log(alpha - parambounds$alpha[1])
  
  get_beta <- function(theta) {
    # theta = -log( (beta - (-.5)) / (1 - (-.5)) )
    (parambounds$beta[2] - parambounds$beta[1]) * 
      exp(-exp(theta)) + parambounds$beta[1]
  }
  get_theta4 <- function(beta) log(
    -log( 
      (beta - parambounds$beta[1]) / (parambounds$beta[2] - parambounds$beta[1]) 
    )
  )
  
  ## Optimization using IPOPT ##
  # The template returns the NEGATIVE log posterior
  # So leave these as negatives. IPOPT will minimize.
  ipopt_objective <- function(theta) ff$fn(theta)
  ipopt_objective_gradient <- function(theta) ff$gr(theta)
  ipopt_objective_hessian <- function(theta) {
    H <- forceSymmetric(ff$he(theta))
    H <- as(H,"dsTMatrix")
    H
  }
  ipopt_objective_hessian_x <- function(theta,obj_factor,hessian_lambda) 
    obj_factor * ipopt_objective_hessian(theta)@x
  ipopt_objective_hessian_structure <- function(theta) {
    H <- ipopt_objective_hessian(theta)
    H <- as(forceSymmetric(H),'dsTMatrix')
    forStruct = cbind(H@i+1, H@j+1)
    tapply(forStruct[,1], forStruct[,2], c)
  }
  
  
  # Box constraints, to improve stability of optimization
  lowerbounds <- c(
    get_theta1(parambounds$Psi0[2] - .001),
    get_theta2(parambounds$gamma[2] - .001),
    get_theta3(parambounds$alpha[1] + .001),
    get_theta4(parambounds$beta[2] - .001)
  )
  
  upperbounds <- c(
    get_theta1(parambounds$Psi0[1] + 1),
    get_theta2(parambounds$gamma[1] + .01),
    get_theta3(parambounds$alpha[2] - .01),
    get_theta4(parambounds$beta[1] + .01)
  )
  
  thetastart <- (lowerbounds + upperbounds)/2 # Start in the middle
  
  # Nonlinear constraints, specified as a function
  ipopt_nonlinear_constraints <- function(theta) Es$fn(theta)
  
  ipopt_nonlinear_constraints_jacobian <- function(theta) {
    J <- Es$gr(theta)
    as(J,"dgTMatrix")
  }
  ipopt_nonlinear_constraints_jacobian_x <- function(theta) 
    ipopt_nonlinear_constraints_jacobian(theta)@x
  ipopt_nonlinear_constraints_jacobian_structure <- function(theta) {
    J <- ipopt_nonlinear_constraints_jacobian(theta)
    J <- as(J,'dgTMatrix')
    forStruct = cbind(J@i+1, J@j+1)
    tapply(forStruct[,2], forStruct[,1], c)
  }
  
  nonlinear_lowerbound <- rep(0,nrow(gcdatalist$y)+2)
  nonlinear_upperbound <- rep(Inf,nrow(gcdatalist$y)+2)
  
  stopifnot(all(ipopt_nonlinear_constraints(thetastart) > 0))
  
  tm <- Sys.time()
  ipopt_result <- ipoptr::ipoptr(
    x0 = thetastart,
    eval_f = ipopt_objective,
    eval_grad_f = ipopt_objective_gradient,
    eval_h = ipopt_objective_hessian_x,
    eval_h_structure = ipopt_objective_hessian_structure(thetastart),
    eval_g = ipopt_nonlinear_constraints,
    eval_jac_g = ipopt_nonlinear_constraints_jacobian_x,
    eval_jac_g_structure = ipopt_nonlinear_constraints_jacobian_structure(thetastart),
    lb = lowerbounds,
    ub = upperbounds,
    constraint_lb = nonlinear_lowerbound,
    constraint_ub = nonlinear_upperbound,
    opts = list(obj_scaling_factor = 1,
                tol = 1e-03)
  )
  optruntime <- difftime(Sys.time(),tm,units = 'secs')
  cat('Run time for mass model optimization:',optruntime,'seconds.\n')
  
  ## AGHQ ----
  # Create the optimization template
  useropt <- list(
    ff = list(
      fn = function(theta) -1*ff$fn(theta),
      gr = function(theta) -1*ff$gr(theta),
      he = function(theta) -1*ff$he(theta)
    ),
    mode = ipopt_result$solution,
    hessian = ff$he(ipopt_result$solution)
  )
  # Do the quadrature
  tm <- Sys.time()
  astroquad <- aghq::aghq(ff,5,thetastart,optresults = useropt,control = default_control(negate=TRUE))
  quadruntime <- difftime(Sys.time(),tm,units = 'secs')
  cat("Run time for mass model quadrature:",quadruntime,"seconds.\n")
  
  # Total runtime
  aghqtime <- optruntime + quadruntime
  
  ## MCMC ----
  tm <- Sys.time()
  stanmod <- tmbstan(
    ff,
    chains = 4,
    cores = 4,
    iter = 1e04,
    warmup = 1e03,
    init = thetastart,
    seed = 48645,
    algorithm = "NUTS"
  )
  stantime <- difftime(Sys.time(),tm,units = 'secs')
  # Save the traceplot
  # pdf(file = file.path(plotpath,"stan-trace.pdf"),width = 7,height = 7)
  # traceplot(stanmod)
  # dev.off()
  
  ## TMB ----
  tm <- Sys.time()
  tmbsd <- TMB::sdreport(ff)
  tmbtime <- difftime(Sys.time(),tm,units = "secs")
  tmbsddat <- data.frame(var = paste0('theta',1:4),est = tmbsd$par.fixed,sd = sqrt(diag(tmbsd$cov.fixed)))
  rownames(tmbsddat) <- NULL
  
  # Times
  # AGHQ
  as.numeric(aghqtime) * stanmod@sim$iter / as.numeric(stantime)
  # TMB
  as.numeric(optruntime) * stanmod@sim$iter / as.numeric(stantime)
  
  
  # Redefine parameters functions for plotting
  get_psi0 <- function(theta)
    (parambounds$Psi0[2] - parambounds$Psi0[1]) * 
    exp(-exp(theta)) + parambounds$Psi0[1]
  get_theta1 <- function(Psi0) 
    log(-log( (Psi0 - parambounds$Psi0[1]) / 
                (parambounds$Psi0[2] - parambounds$Psi0[1]) ))
  
  get_gamma <- function(theta)  
    (parambounds$gamma[2] - parambounds$gamma[1]) * 
    exp(-exp(theta)) + parambounds$gamma[1]
  # Add a little buffer, for stability
  get_theta2 <- function(gamma) 
    log(-log( (gamma - parambounds$gamma[1] + 1e-03) / 
                (parambounds$gamma[2] - parambounds$gamma[1] + 1e-03) ))
  
  get_alpha <- function(theta)
    exp(theta) + parambounds$alpha[1]
  # Add a little buffer, for stability
  get_theta3 <- function(alpha) 
    log(alpha - parambounds$alpha[1] + 1e-03)
  
  
  get_beta <- function(theta)
    (parambounds$beta[2] - parambounds$beta[1]) * 
    exp(-exp(theta)) + parambounds$beta[1]
  get_theta4 <- function(beta) 
    log(-log( (beta - parambounds$beta[1]) / 
                (parambounds$beta[2] - parambounds$beta[1]) ))
  ## Compute the transformed pdfs ##
  translist1 <- list(totheta = get_theta1,fromtheta = get_psi0)
  translist2 <- list(totheta = get_theta2,fromtheta = get_gamma)
  translist3 <- list(totheta = get_theta3,fromtheta = get_alpha)
  translist4 <- list(totheta = get_theta4,fromtheta = get_beta)
  
  psi0pdf <- compute_pdf_and_cdf(astroquad$marginals[[1]],translist1)
  gammapdf <- compute_pdf_and_cdf(astroquad$marginals[[2]],translist2)
  alphapdf <- compute_pdf_and_cdf(astroquad$marginals[[3]],translist3)
  betapdf <- compute_pdf_and_cdf(astroquad$marginals[[4]],translist4)
  
  Psi0prior <- function(Psi0) dunif(Psi0,parambounds$Psi0[1],parambounds$Psi0[2],log = FALSE)
  gammaprior <- function(gamma) dunif(gamma,parambounds$gamma[1],parambounds$gamma[2],log = FALSE)
  alphaprior <- function(alpha) dgamma(alpha - parambounds$alpha[1],shape = 1,rate = 4.6,log = FALSE)
  betaprior <- function(beta) dunif(beta,parambounds$beta[1],parambounds$beta[2],log = FALSE)
  
  # STAN
  standata <- as.data.frame(stanmod)
  standata$psi0 <- get_psi0(standata[ ,1])
  standata$gamma <- get_gamma(standata[ ,2])
  standata$alpha <- get_alpha(standata[ ,3])
  standata$beta <- get_beta(standata[ ,4])
  
  # TMB
  tmbpsi0 <- data.frame(psi0 = psi0pdf$transparam,
                        pdf = dnorm(psi0pdf$theta,tmbsddat[1,2],tmbsddat[1,3]) * abs(numDeriv::grad(get_theta1,psi0pdf$transparam)))
  tmbgamma <- data.frame(gamma = gammapdf$transparam,
                         pdf = dnorm(gammapdf$theta,tmbsddat[2,2],tmbsddat[2,3]) * abs(numDeriv::grad(get_theta2,gammapdf$transparam)))
  tmbalpha <- data.frame(alpha = alphapdf$transparam,
                         pdf = dnorm(alphapdf$theta,tmbsddat[3,2],tmbsddat[3,3]) * abs(numDeriv::grad(get_theta3,alphapdf$transparam)))
  tmbbeta <- data.frame(beta = betapdf$transparam,
                        pdf = dnorm(betapdf$theta,tmbsddat[4,2],tmbsddat[4,3]) * abs(numDeriv::grad(get_theta4,betapdf$transparam)))
  
  
  # pdf(file.path(plotpath,"psi0-plot.pdf"),width = 7,height = 7)
  hist(standata$psi0,freq=FALSE,breaks=50,main = "",xlab = "",cex.lab=1.5,cex.axis = 1.5)
  with(psi0pdf,lines(transparam,pdf_transparam,lwd=2))
  with(psi0pdf,lines(transparam,Psi0prior(transparam),lty = 'dashed',lwd=2))
  with(tmbpsi0,lines(psi0,pdf,lty='dotdash',lwd=2))
  # dev.off()
  
  # pdf(file.path(plotpath,"gamma-plot.pdf"),width = 7,height = 7)
  hist(standata$gamma,freq=FALSE,breaks=50,main = "",xlab = "",cex.lab=1.5,cex.axis = 1.5)
  with(gammapdf,lines(transparam,pdf_transparam,lwd=2))
  with(gammapdf,lines(transparam,gammaprior(transparam),lty = 'dashed',lwd=2))
  with(tmbgamma,lines(gamma,pdf,lty='dotdash',lwd=2))
  # dev.off()
  
  # pdf(file.path(plotpath,"alpha-plot.pdf"),width = 7,height = 7)
  hist(standata$alpha,freq=FALSE,breaks=50,main = "",xlab = "",cex.lab=1.5,cex.axis = 1.5)
  with(alphapdf,lines(transparam,pdf_transparam,lwd=2))
  with(alphapdf,lines(transparam,alphaprior(transparam),lty = 'dashed',lwd=2))
  with(tmbalpha,lines(alpha,pdf,lty='dotdash',lwd=2))
  # dev.off()
  
  # pdf(file.path(plotpath,"beta-plot.pdf"),width = 7,height = 7)
  hist(standata$beta,freq=FALSE,breaks=50,main = "",xlab = "",cex.lab=1.5,cex.axis = 1.5)
  with(betapdf,lines(transparam,pdf_transparam,lwd=2))
  with(betapdf,lines(transparam,betaprior(transparam),lty = 'dashed',lwd=2))
  with(tmbbeta,lines(beta,pdf,lty='dotdash',lwd=2))
  # dev.off()
  
  # KS distances
  getks <- function(x,y) {
    suppressWarnings(capture.output(ks <- ks.test(x,y)))
    unname(ks$statistic)
  }
  M <- nrow(standata)
  
  aghqsamps <- sample_marginal(astroquad,M)
  
  psi0samps <- data.frame(
    mcmc = standata[ ,1],
    aghq = aghqsamps[[1]],
    tmb = rnorm(M,tmbsddat[1,2],tmbsddat[1,3])
  )
  gammasamps <- data.frame(
    mcmc = standata[ ,2],
    aghq = aghqsamps[[2]],
    tmb = rnorm(M,tmbsddat[2,2],tmbsddat[2,3])
  )
  alphasamps <- data.frame(
    mcmc = standata[ ,3],
    aghq = aghqsamps[[3]],
    tmb = rnorm(M,tmbsddat[3,2],tmbsddat[3,3])
  )
  betasamps <- data.frame(
    mcmc = standata[ ,4],
    aghq = aghqsamps[[4]],
    tmb = rnorm(M,tmbsddat[4,2],tmbsddat[4,3])
  )
  
  kstable <- data.frame(
    var = c('psi0','gamma','alpha','beta'),
    ks_aghq = c(
      getks(psi0samps$mcmc,psi0samps$aghq),
      getks(gammasamps$mcmc,gammasamps$aghq),
      getks(alphasamps$mcmc,alphasamps$aghq),
      getks(betasamps$mcmc,betasamps$aghq)
    ),
    ks_tmb = c(
      getks(psi0samps$mcmc,psi0samps$tmb),
      getks(gammasamps$mcmc,gammasamps$tmb),
      getks(alphasamps$mcmc,alphasamps$tmb),
      getks(betasamps$mcmc,betasamps$tmb)
    )
  )
  
  
  # readr::write_csv(kstable,file.path(plotpath,"astro-ks-table.csv"))
  knitr::kable(kstable,digits = 3)
  
  
  # Inference for the mass profile
  Mr <- function(r,theta) {
    p = get_psi0(theta[1])
    g = get_gamma(theta[2])
    
    # Manual unit conversion into "mass of one trillion suns" (so awesome)
    g*p*r^(1-g) * 2.325e09 * 1e-12
  }
  
  rtodo <- 1:150
  Mrout <- numeric(length(rtodo))
  Mrsdout <- numeric(length(rtodo))
  for (rr in 1:length(rtodo)) {
    r <- rtodo[rr]
    
    Mrout[rr] <- compute_moment(
      astroquad,
      function(x) Mr(r,x)
    )
    Mrsdout[rr] <- sqrt(compute_moment(
      astroquad,
      function(x) (Mr(r,x) - Mrout[rr])^2
    ))
  }
  # pdf(file.path(plotpath,"massplot-aghq.pdf"),width=7,height=7)
  plot(rtodo,Mrout,type='l',lwd=2,xlab="",ylab="",xaxt='n',cex.axis=1.5)
  title(ylab=bquote('M(r) ('~10^12~M[sun]~')'),cex.lab=1.5,line=2.3)
  axis(1,at=seq(0,150,by=25),cex.axis=1.5)
  lines(rtodo,Mrout - 2*Mrsdout,lty='dashed',lwd=2)
  lines(rtodo,Mrout + 2*Mrsdout,lty='dashed',lwd=2)
  # dev.off()
  
  # With MCMC
  Mrlist <- list()
  for (i in 1:length(rtodo)) Mrlist[[i]] <- apply(standata[ ,c(1,2)],1,function(x) Mr(rtodo[i],x))
  meanvals <- Reduce(c,Map(mean,Mrlist))
  lowervals <- Reduce(c,Map(quantile,Mrlist,probs = .025))
  uppervals <- Reduce(c,Map(quantile,Mrlist,probs = .975))
  
  
  # pdf(file.path(plotpath,"massplot-mcmc.pdf"),width=7,height=7)
  plot(rtodo,meanvals,type='l',lwd=2,xlab="",ylab="",xaxt='n',cex.axis=1.5)
  title(ylab=bquote('M(r) ('~10^12~M[sun]~')'),cex.lab=1.5,line=2.3)
  axis(1,at=seq(0,150,by=25),cex.axis=1.5)
  lines(rtodo,lowervals,lty='dashed',lwd=2)
  lines(rtodo,uppervals,lty='dashed',lwd=2)
  # dev.off()
  
  # Empirical RMSE
  sqrt(mean( (Mrout - meanvals)^2 )) # 0.0004119291
  sqrt(mean( ((Mrout - 2*Mrsdout) - lowervals)^2 )) # 0.007427573
  sqrt(mean( ((Mrout + 2*Mrsdout) - uppervals)^2 )) # 0.006063218
  
  #### END EXAMPLE 4.2 ####
}

## Example 5.1: Loaloa, without zero-inflation ----
if (doloaloa) {

  ilogit <- function(x) 1 / (1 + exp(-x))
  
  set.seed(78968)
  
  # Set flags
  savestamp <- "20210504-v1"
  plotpath <- file.path(globalpath,"loaloa")
  if (!dir.exists(plotpath)) dir.create(plotpath)
  savepath <- plotpath
  
  # Initialize time variables
  aghqtime <- 0
  inlatime <- 0
  mcmltime <- 0
  mcmctime <- 0
  
  aghqsimtime <- 0
  mcmlsimtime <- 0
  mcmcsimtime <- 0
  
  ## Load the data ##
  # loaloa data appears in two packages
  data(loaloa,package = "PrevMap")
  loaloa2 <- loaloa
  rm(loaloa)
  data(loaloa,package = "geostatsp")
  
  ## Prepare the "inner" model ##
  
  # Design matrices
  Amat <- Diagonal(nrow(loaloa))
  
  design <- cbind(Amat,rep(1,nrow(Amat))) # No covariates, for comparison
  # Response
  y <- loaloa@data$y
  N <- loaloa@data$N
  
  ## Dimensions
  n <- nrow(design)
  p <- 1
  m <- ncol(Amat)
  Wd <- ncol(design)
  
  ## Prior distributions ##
  sigma_u <- 1
  sigma_alpha <- .025
  rho_u <- 2e05
  rho_alpha <- .975
  # rho is "range", WITH the sqrt(8*nu) factor. So PC prior is exponential on 1/rho
  
  # PC Prior for kappa,tau
  maternconstants <- list()
  maternconstants$d <- 2 # Dimension of spatial field, fixed
  maternconstants$nu <- 1 # Shape param, fixed
  get_kappa <- function(sigma,rho)
    sqrt(8*maternconstants$nu)/rho
  get_tau <- function(sigma,rho)
    sigma * get_kappa(sigma,rho)^(maternconstants$nu) *
    sqrt(gamma(maternconstants$nu +
                 maternconstants$d/2) *
           (4*pi)^(maternconstants$d/2) /
           gamma(maternconstants$nu))
  get_sigma <- function(kappa,tau)
    tau /
    (kappa^(maternconstants$nu) *
       sqrt(gamma(maternconstants$nu +
                    maternconstants$d/2) *
              (4*pi)^(maternconstants$d/2) /
              gamma(maternconstants$nu)))
  get_rho <- function(kappa,tau)
    sqrt(8*maternconstants$nu) / kappa
  
  log_prior_theta <- function(theta) {
    # theta = (log(kappa),log(tau))
    kappa <- exp(theta[1])
    tau <- exp(theta[2])
    
    lambda1 <- -(rho_u / sqrt(8*maternconstants$nu)) ^
      (maternconstants$d/2) * log(rho_alpha)
    lambda2 <- -kappa^(-maternconstants$nu) *
      sqrt( gamma(maternconstants$nu) /
              ( gamma(maternconstants$nu +
                        maternconstants$d/2) *
                  (4*pi)^(maternconstants$d/2) ) ) *
      log(sigma_alpha) / sigma_u
    
    log(maternconstants$d) -
      log(2) +
      log(lambda1) +
      log(lambda2) +
      (maternconstants$d/2 - 1) *
      log(kappa) -
      lambda1 *
      kappa^(maternconstants$d/2) -
      lambda2 * tau +
      sum(theta)
  }
  
  beta_prec <- .001
  
  Q_matrix <- function(theta) {
    # theta = log(kappa), log(tau)
    
    kappa <- as.numeric(unname(exp(theta[1])))
    tau <- as.numeric(unname(exp(theta[2])))
    
    sig <- get_sigma(kappa,tau)
    rho <- get_rho(kappa,tau)
    
    # Matern
    mm <- geostatsp::matern(
      loaloa,
      param = c("variance" = sig^2,"range" = rho,"shape" = maternconstants$nu),
      type = "precision"
    )
    
    bb <- beta_prec * Diagonal(p)
    
    rbind(
      cbind(mm,Matrix(0,nrow = nrow(mm),ncol = p,sparse = FALSE)),
      cbind(Matrix(0,nrow = p,ncol = ncol(mm),sparse = FALSE),bb)
    )
  }
  
  log_prior_W <- function(W,theta,Q = NULL) {
    if (is.null(Q)) Q <- Q_matrix(theta)
    -(1/2) * as.numeric(crossprod(W,crossprod(Q,W))) + 
        .5 * determinant(Q,logarithm=TRUE)$modulus
  }
  
  ## Likelihood ----
  
  make_eta <- function(W) as.numeric(design %*% W)
  
  logit <- function(x) log((x) / (1-x))
  ilogit <- function(x) 1 / (1 + exp(-x))
  
  log_likelihood <- function(W) {
    eta <- make_eta(W)
    p <- ilogit(eta)
    sum(y * log(p) + (N - y) * log(1 - p))
  }
  
  grad_log_likelihood <- function(W) {
    eta <- make_eta(W)
    p <- ilogit(eta)
    y - N * p
  }
  
  # NEGATIVE hessian
  hessian_log_likelihood <- function(W) {
    eta <- make_eta(W)
    p <- ilogit(eta)
    Diagonal(length(N),N * p * (1-p))
  }
  
  ## Posterior
  
  log_posterior_W <- function(W,theta,Q = NULL) {
    if (is.null(Q)) Q <- Q_matrix(theta)
    log_prior_W(W,theta,Q) + log_likelihood(W)
  }
  
  grad_log_posterior_W <- function(W,theta,Q = NULL) {
    if (is.null(Q)) Q <- Q_matrix(theta)
    as.numeric(-Q %*% cbind(W) + t(design) %*% grad_log_likelihood(W))
  }
  
  H_matrix <- function(W,theta,Q = NULL) {
    # minus the hessian of the log posterior
    CE <- hessian_log_likelihood(W)
    if (is.null(Q)) Q <- Q_matrix(theta)
    
    CW <- crossprod(design,crossprod(CE,design))
    
    as(Q + CW,"dgCMatrix")
  }
  
  # Simulate the spatial fields
  simulate_spatial_fields <- function(U,
                                      theta,
                                      pointsdata,
                                      resolution = list(nrow = 100,ncol = 100)) {
    # U: matrix of samples, each column is a sample
    # theta: data.frame of theta values
    # Draw from U*|U
    
    # Compute matrix of var, range, shape
    modpar <- cbind(
      var = get_sigma(exp(theta$theta1),exp(theta$theta2))^2,
      range = get_rho(exp(theta$theta1),exp(theta$theta2)),
      shape = maternconstants$nu
    )
    
    fielddat <- pointsdata
    fielddat@data <- as.data.frame(U)
    
    geostatsp::RFsimulate(
      model = modpar,
      data = fielddat,
      x = raster(fielddat,nrows = resolution$nrow,ncols = resolution$ncol)
    )
    
  }
  
  ## AGHQ ----
  
  log_posterior_joint <- function(W,theta) {
    Q <- Q_matrix(theta) # Prior precision matrix
    log_prior_theta(theta) + log_prior_W(W,theta,Q) + log_likelihood(W)
  }
  
  ff <- list(
    fn = log_posterior_joint,
    gr = grad_log_posterior_W,
    he = function(W,theta) -1 * H_matrix(W,theta)
  )
  
  # Starting values taken from Brown (2011)
  startingsig <- .988
  startingrho <- 4.22*1e04
  startingtheta <- c(
    log(get_kappa(startingsig,startingrho)),
    log(get_tau(startingsig,startingrho))
  )
  Wdim <- dim(Q_matrix(c(0,0)))[1]
  
  # AGHQ
  tm <- Sys.time()
  cat("Doing AGHQ, time = ",format(tm),"\n")
  loaloaquad <- aghq::marginal_laplace(
    ff,
    5,
    list(W = rep(0,Wdim),theta = startingtheta)
  )
  aghqtime <- difftime(Sys.time(),tm,units = 'secs')
  # save(loaloaquad,file = file.path(savepath,paste0("loaloaquad",savestamp,".RData")))
  cat("AGHQ took: ",format(aghqtime),"\n")
  
  # Spatial interpolation for AGHQ
  tm <- Sys.time()
  cat("Doing AGHQ simulation, time = ",format(tm),"\n")
  samps <- sample_marginal(loaloaquad,1e02) # very fast
  
  # predictions (prevalence)
  UU <- samps$samps[1:190, ]
  beta <- samps$samps[191, ]
  
  simbrick <- simulate_spatial_fields(
    U = UU,
    theta = samps$theta,
    pointsdata = loaloa,
    resolution = list(nrow = 25,ncol = 50)
  )
  simbrick <- simbrick + beta
  aghqsimtime <- difftime(Sys.time(),tm,units = 'secs')
  # raster::writeRaster(simbrick,file = file.path(savepath,paste0("loaloa-simbrick-aghq",savestamp,".grd")),overwrite = TRUE)
  cat("AGHQ simulation took: ",format(aghqsimtime),"\n")
  
  
  ## INLA ----
  # Using geostatsp::glgm()
  
  tm <- Sys.time()
  cat("Doing INLA, time = ",format(tm),"\n")
  loaFit = glgm(formula = y ~ 1,
                data = loaloa, grid = 50,
                family = "binomial",
                Ntrials = loaloa$N, shape = maternconstants$nu, buffer = 1e05,
                prior = list(sd = c(u = sigma_u, alpha = sigma_alpha), range = c(u = rho_u, alpha = rho_alpha))) 
  inlatime <- difftime(Sys.time(),tm,units = 'secs')
  # save(loaFit,file = file.path(savepath,paste0("loaloainla",savestamp,".RData")))
  cat("INLA took: ",format(inlatime),"\n")
  
  
  # prediction
  # Border for Cameroon and Nigeria
  cameroonBorderLL = raster::getData("GADM", country=c('CMR'), level=2)
  nigeriaBorderLL = raster::getData("GADM", country=c('NGA'), level=2)
  cameroonBorder = spTransform(cameroonBorderLL, projection(loaloa))
  nigeriaBorder = spTransform(nigeriaBorderLL, projection(loaloa))
  cameroonBorderouter <- rgeos::gUnaryUnion(cameroonBorder)
  nigeriaBorderouter <- rgeos::gUnaryUnion(nigeriaBorder)
  
  
  fullborder <- raster::bind(cameroonBorder,nigeriaBorder)
  fullborderouter <- raster::bind(cameroonBorderouter,nigeriaBorderouter)
  
  fullborder <- crop(fullborder,loaloa)
  fullborderouter <- crop(fullborderouter,loaloa)
  
  plot_loaloa <- function(plotraster,breaks) {
    predcols <- mapmisc::colourScale(
      plotraster,
      breaks = breaks,
      style = "fixed",
      col = "Spectral",
      rev = TRUE
    )
    
    plotraster <- mask(plotraster,fullborderouter)
    
    mapmisc::map.new(loaloa)
    plot(plotraster,
         col = predcols$col,
         breaks = predcols$breaks,
         legend=FALSE, add=TRUE)
    points(loaloa,pch = 4)
    plot(fullborder,add = TRUE,border=mapmisc::col2html("black", 0.5), lwd=0.5)
    plot(fullborderouter,add = TRUE)
    mapmisc::legendBreaks('right', predcols, cex=1, bty='o',inset = .05)
  }
  
  
  
  ## MCML ----
  # prevmap
  tm <- Sys.time()
  cat("Doing MCML, time = ",format(tm),"\n")
  initvaluemodel <- glm(cbind(y,N-y)~1,data = loaloa@data,family=binomial)
  initbeta0 <- coef(initvaluemodel)
  loaloa2$logit <- log((loaloa2$NO_INF + 0.5)/(loaloa2$NO_EXAM - loaloa2$NO_INF + 0.5))
  
  vari <- variog(coords = as.matrix(loaloa2[, c("LONGITUDE", "LATITUDE")]),
                 data = loaloa2$logit,
                 uvec = c(0, 0.1, 0.15, 0.2, 0.4, 0.8, 1.4, 1.8, 2, 2.5, 3))
  vari.fit <- variofit(vari, 
                       ini.cov.pars = c(2, 0.2),
                       cov.model = "matern",
                       fix.nugget = FALSE, 
                       nugget = 0 ,
                       fix.kappa = TRUE, 
                       kappa = 0.5)
  
  # change "phi" to be in meters, 10^4 starting value
  par0 <- c(initbeta0, vari.fit$cov.pars, vari.fit$nugget)
  c.mcmc <- control.mcmc.MCML(
    n.sim = 10000, burnin = 2000,
    thin = 8, h = (1.65)/(nrow(loaloa2) ^ (1/6)))
  
  tmp <- capture.output({
    fit.MCML1 <- binomial.logistic.MCML(formula = y ~ 1,units.m = ~ N,
                                        par0 = par0[1:3],
                                        coords = ~ long + lat, data = as.data.frame(loaloa),
                                        control.mcmc = c.mcmc,
                                        kappa = maternconstants$nu, 
                                        fixed.rel.nugget = 0,
                                        start.cov.pars = c(1e04))
  })
  mcmltime <- difftime(Sys.time(),tm,units = 'secs')
  cat("MCML took: ",format(mcmltime),"\n")
  
  tm <- Sys.time()
  cat("Doing MCML Sims, time = ",format(tm),"\n")
  tmp <- capture.output({
    preds_mcml <- spatial.pred.binomial.MCML(
      fit.MCML1,
      as.data.frame(raster(loaFit$raster$random.mean),xy=TRUE),
      control.mcmc = c.mcmc,
      type = 'marginal',
      scale.predictions = 'prevalence',
      standard.errors = TRUE,
      thresholds = .2,
      scale.thresholds = 'prevalence'
    )
  })
  
  # rasterize the grid
  mcmcl_predraster <- raster(loaFit$raster$random.mean)
  values(mcmcl_predraster) <- preds_mcml$prevalence$predictions
  mcmlsimtime <- difftime(Sys.time(),tm,units = 'secs')
  # save(fit.MCML1,preds_mcml,file = file.path(savepath,paste0("loaloamcml",savestamp,".RData")))
  # raster::writeRaster(mcmcl_predraster,file.path(savepath,paste0("loaloamcml-preds",savestamp,".grd")),overwrite = TRUE)
  cat("MCML Sims took: ",format(mcmlsimtime),"\n")
  
  ## MCMC ----
  tm <- Sys.time()
  cat("Doing MCMC, time = ",format(tm),"\n")
  # priors
  lam <- -log(rho_alpha) * rho_u
  thepriors <- control.prior(
    beta.mean = 0,
    beta.covar = 1/beta_prec,
    log.prior.sigma2 = function(x) stats::dexp(sqrt(x), rate=0.922219863528484,log = TRUE) - (1/2)*log(x),
    log.prior.phi = function(x) -2*log(x) + dexp(1/x,sqrt(8)*lam,log = TRUE),
    log.normal.nugget = c(-15,1e-06) # Basically no nugget
  )
  
  mcmc.Bayes <- control.mcmc.Bayes(
    n.sim = 6000, burnin = 1000, thin = 1,h.theta1 = 1, h.theta2 = 0.7, h.theta3 = 0.05,
    L.S.lim = c(5,50), epsilon.S.lim = c(0.03,0.06),
    start.beta = -2.3, 
    start.sigma2 = 2.6,
    start.phi = 0.8, 
    start.nugget = exp(-15),
    start.S = predict(initvaluemodel)
  )
  
  tmp <- capture.output({
    fit.Bayes <- binomial.logistic.Bayes(
      formula = y ~ 1,
      units.m = ~ N,
      coords = ~ long + lat,
      data = as.data.frame(loaloa), control.prior = thepriors,
      control.mcmc = mcmc.Bayes, kappa = maternconstants$nu)
    summary(fit.Bayes, hpd.coverage = 0.95)
  })
  
  mcmctime <- difftime(Sys.time(),tm,units = 'secs')
  cat("MCMC took: ",format(mcmctime),"\n")
  
  tm <- Sys.time()
  cat("Doing MCMC Sims, time = ",format(tm),"\n")
  
  tmp <- capture.output({
    pred.Bayes <- spatial.pred.binomial.Bayes(
      fit.Bayes, 
      as.data.frame(raster(loaFit$raster$random.mean),xy=TRUE),
      type = "marginal",
      scale.predictions = "prevalence", quantiles = NULL,
      standard.errors = TRUE, thresholds = 0.2,
      scale.thresholds = "prevalence")
  })
  
  bayesian_predraster <- raster(loaFit$raster$random.mean)
  values(bayesian_predraster) <- pred.Bayes$prevalence$predictions
  
  mcmcsimtime <- difftime(Sys.time(),tm,units = 'secs')
  # save(fit.Bayes,pred.Bayes,file = file.path(savepath,paste0("loaloamcmc",savestamp,".RData")))
  # raster::writeRaster(bayesian_predraster,file.path(savepath,paste0("loaloamcmc-preds",savestamp,".grd")),overwrite = TRUE)
  cat("MCMC Sims took: ",format(mcmcsimtime),"\n")
  
  
  # Plot them!
  
  
  br <- c(0,.05,.1,.15,.2,.25,.3,.4,.5,.6,1)
  
  aghqpostmean <- calc(ilogit(simbrick),mean)
  
  # png(file.path(plotpath,paste0("aghq-mean-map",savestamp,".png")))
  plot_loaloa(aghqpostmean,br)
  # dev.off()
  
  # png(file.path(plotpath,paste0("mcml-mean-map",savestamp,".png")))
  plot_loaloa(mcmcl_predraster,br)
  # dev.off()
  
  # png(file.path(plotpath,paste0("inla-mean-map",savestamp,".png")))
  plot_loaloa(loaFit$raster$predict.invlogit,br)
  # dev.off()
  
  # png(file.path(plotpath,paste0("mcmc-mean-map",savestamp,".png")))
  plot_loaloa(bayesian_predraster,br)
  # dev.off()
  # difference
  difftable <- data.frame(
    comparison = c('mean','max'),
    mcml_mcmc = c(
      mean(abs(values(mcmcl_predraster) - values(bayesian_predraster))),
      max(abs(values(mcmcl_predraster) - values(bayesian_predraster)))
    ),
    inla_mcmc = c(
      mean(abs(values(loaFit$raster$predict.invlogit) - values(bayesian_predraster))),
      max(abs(values(loaFit$raster$predict.invlogit) - values(bayesian_predraster)))
    ),
    inla_mcml = c(
      mean(abs(values(loaFit$raster$predict.invlogit) - values(mcmcl_predraster))),
      max(abs(values(loaFit$raster$predict.invlogit) - values(mcmcl_predraster)))
    ),
    aghq_mcml = c(
      mean(abs(values(aghqpostmean) - values(mcmcl_predraster))),
      max(abs(values(aghqpostmean) - values(mcmcl_predraster)))
    ),
    aghq_mcmc = c(
      mean(abs(values(aghqpostmean) - values(bayesian_predraster))),
      max(abs(values(aghqpostmean) - values(bayesian_predraster)))
    ),
    aghq_inla = c(
      mean(abs(values(aghqpostmean) - values(loaFit$raster$predict.invlogit))),
      max(abs(values(aghqpostmean) - values(loaFit$raster$predict.invlogit)))
    )
  )
  # readr::write_csv(difftable,file.path(savepath,paste0("diff-table-loaloa",savestamp,".csv")))
  knitr::kable(difftable,digits = 3)
  
  # covariance params
  aghqsigmamean <- compute_moment(loaloaquad,function(x) get_sigma(exp(x[1]),exp(x[2])))
  aghqrhomean <- compute_moment(loaloaquad,function(x) get_rho(exp(x[1]),exp(x[2])))
  aghqlogsigmamean <- compute_moment(loaloaquad,function(x) log(get_sigma(exp(x[1]),exp(x[2]))))
  aghqlogrhomean <- compute_moment(loaloaquad,function(x) log(get_rho(exp(x[1]),exp(x[2]))))
  
  aghqlogsigmasd <- sqrt(compute_moment(loaloaquad,function(x) (log(get_sigma(exp(x[1]),exp(x[2]))) - aghqlogsigmamean)^2))
  aghqlogrhosd <- sqrt(compute_moment(loaloaquad,function(x) (log(get_rho(exp(x[1]),exp(x[2]))) - aghqlogrhomean)^2))
  
  covparamtable <- data.frame(
    method = c('AGHQ','INLA','MCML','MCMC'),
    sigmamean = c(
      aghqsigmamean,
      loaFit$parameters$summary$mean[3],
      exp(summary(fit.MCML1)$cov.par[1,1] / 2),
      mean(sqrt(fit.Bayes$estimate[ ,'sigma^2']))
    ),
    sigmalower = c(
      exp(aghqlogsigmamean - 2*aghqlogsigmasd),
      loaFit$parameters$summary$`0.025quant`[3],
      exp((summary(fit.MCML1)$cov.par[1,1] - 2*summary(fit.MCML1)$cov.par[1,2])/2),
      quantile(sqrt(fit.Bayes$estimate[ ,'sigma^2']),probs = .025)
    ),
    sigmaupper = c(
      exp(aghqlogsigmamean + 2*aghqlogsigmasd),
      loaFit$parameters$summary$`0.975quant`[3],
      exp((summary(fit.MCML1)$cov.par[1,1] + 2*summary(fit.MCML1)$cov.par[1,2])/2),
      quantile(sqrt(fit.Bayes$estimate[ ,'sigma^2']),probs = .975)
    ),
    
    rhomean = c(
      aghqrhomean,
      loaFit$parameters$summary$mean[2] * 1000,
      exp(summary(fit.MCML1)$cov.par[2,1]) * sqrt(8),
      mean(fit.Bayes$estimate[ ,'phi']) * sqrt(8)
    ),
    rholower = c(
      exp(aghqlogrhomean - 2*aghqlogrhosd),
      loaFit$parameters$summary$`0.025quant`[2] * 1000,
      exp((summary(fit.MCML1)$cov.par[2,1] - 2*summary(fit.MCML1)$cov.par[2,2]))*sqrt(8),
      quantile(fit.Bayes$estimate[ ,'phi'] * sqrt(8),probs = .025)
    ),
    rhoupper = c(
      exp(aghqlogrhomean + 2*aghqlogrhosd),
      loaFit$parameters$summary$`0.975quant`[2] * 1000,
      exp((summary(fit.MCML1)$cov.par[2,1] + 2*summary(fit.MCML1)$cov.par[2,2]))*sqrt(8),
      quantile(fit.Bayes$estimate[ ,'phi'] * sqrt(8),probs = .975)
    )
  )
  # readr::write_csv(covparamtable,file.path(savepath,paste0("covparam-table-loaloa",savestamp,".csv")))
  knitr::kable(covparamtable,digits = 3)
  
  # Write the timing table
  timingtable <- data.frame(
    task = c("AGHQ","INLA","MCML","MCMC","AGHQSim","MCMLSim","MCMCSim"),
    time = c(aghqtime,inlatime,mcmltime,mcmctime,aghqsimtime,mcmlsimtime,mcmcsimtime)
  )
  # readr::write_csv(timingtable,file.path(savepath,paste0("timing-table-loaloa",savestamp,".csv")))
  knitr::kable(timingtable,digits = 3)
  
  # Total and num iter
  totaltable <- c(
    mcmc = as.numeric(timingtable[timingtable$task == 'MCMC','time']) + as.numeric(timingtable[timingtable$task == 'MCMCSim','time']),
    aghq = as.numeric(timingtable[timingtable$task == 'AGHQ','time']) + as.numeric(timingtable[timingtable$task == 'AGHQSim','time']),
    inla = as.numeric(timingtable[timingtable$task == 'INLA','time']),
    mcml = as.numeric(timingtable[timingtable$task == 'MCML','time']) + as.numeric(timingtable[timingtable$task == 'MCMLSim','time'])
  )
  effiter <- totaltable[2:4] * 6000 / totaltable[1]
  effiter
  
  #### END OF EXAMPLE 5.1 ####
}

## Example 5.2: Loaloa with zero-inflation

savestamp <- "20210505-v1"
plotpath <- file.path(globalpath,"loaloazip")
if (!dir.exists(plotpath)) dir.create(plotpath)
savepath <- plotpath

file.copy(system.file('extsrc/05_loaloazip.cpp',package='aghq'),globalpath)

# Compile TMB template-- only need to do once
compile(file.path(globalpath,"05_loaloazip.cpp"))
dyn.load(dynlib(file.path(globalpath,"05_loaloazip")))

# Initialize time variables
aghqtime <- 0
aghqsimtime <- 0
mcmctime <- 0
mcmcsimtime <- 0

## Prepare the "inner" model ##

# Design matrices
Amat <- Diagonal(nrow(loaloa))

Xmat <- cbind(rep(1,nrow(Amat)))
# Design matrix: zip model and risk model are the same
design <- bdiag(
  # ZIP
  cbind(
    Amat,
    Xmat
  ),
  # Risk
  cbind(
    Amat,
    Xmat
  )
)

# Response
y <- loaloa@data$y
N <- loaloa@data$N

## Dimensions
n <- nrow(Xmat) # Number of obs
p <- ncol(Xmat) * 2 # Number of betas
m <- ncol(Amat) * 2 # Number of spatial points
Wd <- ncol(design) # Number of total params
# Check
stopifnot(Wd == m + p)

## Prior distributions ##
# Use the same prior for both sets of Matern params
sigma_u <- 1
sigma_alpha <- .025
rho_u <- 2e05
rho_alpha <- .975

# PC Prior for kappa,tau
maternconstants <- list()
maternconstants$d <- 2 # Dimension of spatial field, fixed
maternconstants$nu <- 1 # Shape param, fixed
get_kappa <- function(sigma,rho) sqrt(8*maternconstants$nu)/rho
get_tau <- function(sigma,rho) sigma * get_kappa(sigma,rho)^(maternconstants$nu) * sqrt(gamma(maternconstants$nu + maternconstants$d/2) * (4*pi)^(maternconstants$d/2) / gamma(maternconstants$nu))
get_sigma <- function(kappa,tau) tau / (kappa^(maternconstants$nu) * sqrt(gamma(maternconstants$nu + maternconstants$d/2) * (4*pi)^(maternconstants$d/2) / gamma(maternconstants$nu)))
get_rho <- function(kappa,tau) sqrt(8*maternconstants$nu) / kappa

# Precision for betas

beta_prec <- .001

## Log Posterior ----

startingsig <- 1
startingrho <- 4.22*1e04

datlist <- list(
  y = y,
  N = N,
  design = design,
  nu = maternconstants$nu,
  rho_u = rho_u,
  rho_alpha = rho_alpha,
  sigma_u = sigma_u,
  sigma_alpha = sigma_alpha,
  D = raster::pointDistance(loaloa,lonlat = FALSE),
  betaprec = beta_prec
)
# NOTE: for some initial values of W, TMB's inner optimization seems to fail
# This was tried over a bunch of random intializations and most worked, and all
# gave the same optimum. But this is why we set the seed here and use a random start.
set.seed(4564)
paraminit <- list(
  W = rnorm(ncol(design)),
  logkappa = log(get_kappa(startingsig,startingrho)),
  logtau = log(get_tau(startingsig,startingrho))
)

ff <- MakeADFun(data = datlist,
                parameters = paraminit,
                random = "W",
                DLL = "05_loaloazip",
                ADreport = FALSE,
                silent = TRUE)

tm <- Sys.time()
cat("Doing AGHQ, time = ",format(tm),"\n")
loaloazipquad <- aghq::marginal_laplace_tmb(
  ff,
  3,
  startingvalue = c(paraminit$logkappa,paraminit$logtau)
)
aghqtime <- difftime(Sys.time(),tm,units = 'secs')
# save(loaloazipquad,file = file.path(savepath,paste0("loaloazipquad",savestamp,".RData")))
cat("AGHQ took: ",format(aghqtime),"\n")

simulate_spatial_fields <- function(U,
                                    theta,
                                    pointsdata,
                                    resolution = list(nrow = 100,ncol = 100)) {
  # U: matrix of samples, each column is a sample
  # theta: data.frame of theta values
  # Draw from U*|U
  
  # Compute matrix of var, range, shape
  modpar <- cbind(
    var = get_sigma(exp(theta$theta1),exp(theta$theta2))^2,
    range = get_rho(exp(theta$theta1),exp(theta$theta2)),
    shape = maternconstants$nu
  )
  
  fielddat <- pointsdata
  fielddat@data <- as.data.frame(U)
  
  geostatsp::RFsimulate(
    model = modpar,
    data = fielddat,
    x = raster(fielddat,nrows = resolution$nrow,ncols = resolution$ncol)
  )
}

## Post samples ----

cat("Doing AGHQ field simulations, time = ",format(tm),"\n")
loazippostsamples <- sample_marginal(loaloazipquad,100)
# Extract out the U and V
postU <- loazippostsamples$samps[c(1:190), ]
postV <- loazippostsamples$samps[c(192:381), ]

postBeta <- loazippostsamples$samps[c(191,382), ]

tm <- Sys.time()
fieldbrickzip <- simulate_spatial_fields(
  postU,
  loazippostsamples$theta,
  loaloa,
  resolution = reslist
)
fieldbrickrisk <- simulate_spatial_fields(
  postV,
  loazippostsamples$theta,
  loaloa,
  resolution = reslist
)

fieldbrickzip_withint <- fieldbrickzip + postBeta[1, ]
fieldbrickrisk_withint <- fieldbrickrisk + postBeta[2, ]

simfieldsmeanzip <- mean(1 / (1 + exp(-fieldbrickzip_withint)))

simfieldsmeanrisk <- mean(
  (1 / (1 + exp(-fieldbrickzip_withint))) *
    (1 / (1 + exp(-fieldbrickrisk_withint)))
)
aghqsimtime <- difftime(Sys.time(),tm,units = 'secs')
# raster::writeRaster(fieldbrickzip,file.path(savepath,paste0("fieldbrickzipaghq",savestamp,".grd")),overwrite = TRUE)
# raster::writeRaster(fieldbrickrisk,file.path(savepath,paste0("fieldbrickzipaghq",savestamp,".grd")),overwrite = TRUE)

# cameroonBorderLL = raster::getData("GADM", country=c('CMR'), level=2)
# nigeriaBorderLL = raster::getData("GADM", country=c('NGA'), level=2)
# cameroonBorder = spTransform(cameroonBorderLL, projection(loaloa))
# nigeriaBorder = spTransform(nigeriaBorderLL, projection(loaloa))
# cameroonBorderouter <- rgeos::gUnaryUnion(cameroonBorder)
# nigeriaBorderouter <- rgeos::gUnaryUnion(nigeriaBorder)
# 
# 
# fullborder <- raster::bind(cameroonBorder,nigeriaBorder)
# fullborderouter <- raster::bind(cameroonBorderouter,nigeriaBorderouter)
# 
# fullborder <- crop(fullborder,loaloa)
# fullborderouter <- crop(fullborderouter,loaloa)

plot_loaloa <- function(plotraster,breaks) {
  predcols <- mapmisc::colourScale(
    plotraster,
    breaks = breaks,
    style = "fixed",
    col = "Spectral",
    rev = TRUE
  )
  
  plotraster <- mask(plotraster,fullborderouter)
  
  mapmisc::map.new(loaloa)
  plot(plotraster,
       col = predcols$col,
       breaks = predcols$breaks,
       legend=FALSE, add=TRUE)
  points(loaloa,pch = 4)
  plot(fullborder,add = TRUE,border=mapmisc::col2html("black", 0.5), lwd=0.5)
  plot(fullborderouter,add = TRUE)
  mapmisc::legendBreaks('right', predcols, cex=1, bty='o',inset = .05)
}

br <- c(0,.05,.1,.15,.2,.25,.3,.4,.5,.6,1)
brzero <- c(.35,.5,.6,.7,.8,.9,.91,.92,.93,.94,.95,1)

# pdf(file.path(plotpath,paste0("loaloa-zip-postmean.pdf")),width=7,height=7)
plot_loaloa(simfieldsmeanzip,brzero)
# dev.off()
# pdf(file.path(plotpath,paste0("loaloa-risk-postmean.pdf")),width=7,height=7)
plot_loaloa(simfieldsmeanrisk,br)
# dev.off()
cat("AGHQ simulations took: ",format(aghqsimtime),"\n")

# Write the timing table
timingtable <- data.frame(
  task = c("AGHQ","MCMC","AGHQSim","MCMCSim"),
  time = c(aghqtime,mcmctime,aghqsimtime,mcmcsimtime)
)
# readr::write_csv(timingtable,file.path(savepath,paste0("timing-table",savestamp,".csv")))
knitr::kable(timingtable,digits = 3)


#### END EXAMPLE 5.2 ####

## Example 6.1: the zero-inflated overdispersed Poisson from glmmTMB ----

plotpath <- file.path(globalpath,"poisson-zip")
if (!dir.exists(plotpath)) dir.create(plotpath)

savestamp <- "-2021-05-06"

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
# pdf(file.path(plotpath,paste0("sigma-plot-aghq",savestamp,".pdf")),width = 7,height = 7)
with(zipsigmapdf[[1]],plot(transparam,pdf_transparam,type='l',main='',ylab='',xlab=expression(sigma),xaxt='n',xlim = c(0,1.2),ylim = c(0,3),cex.lab=1.5,cex.axis = 1.5))
title(ylab="Density",cex.lab=1.5,line=2.3)
axis(1,cex.axis=1.5,at=seq(0,1.2,by=.2))
abline(v = sqrt(as.numeric(summary(zipmod)$varcor$cond$site)),lty='dashed')
# dev.off()

# confidence/credible intervals for random effects
rr <- as.data.frame(ranef(zipmod,condVar = TRUE))
Usamps <- zipquadsamps$samps[rownames(zipquadsamps$samps) == 'b', ]
Uest <- apply(Usamps,1,mean)
Ulower <- apply(Usamps,1,quantile,probs = .025)
Uupper <- apply(Usamps,1,quantile,probs = .975)

xx <- 1:nrow(rr)
offset <- .2
# pdf(file.path(plotpath,paste0("interval-plot-aghq",savestamp,".pdf")),width = 7,height = 7)
with(rr,plot(xx-offset,condval,type = 'p',pch=20,xaxt='n',xlab='Site',ylab='',ylim = c(-1.7,1.7)))
points(xx+offset,Uest,type='p',pch=20)
with(rr,axis(1,at = xx,labels = as.character(grp)))
with(rr,arrows(x0=xx-offset,y0=condval-qnorm(.975)*condsd,y1=condval+qnorm(.975)*condsd,length=.05,angle=90,code=3))
arrows(x0=xx+offset,y0=Ulower,y1=Uupper,length=.05,angle=90,code=3)
# dev.off()

#### END EXAMPLE 6.1 ####

sessionInfo()