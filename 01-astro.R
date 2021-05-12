### Astro example for AGHQ software paper ###

## Setup ----

library(aghq)
library(TMB)
precompile()
library(tmbstan)
library(parallel)
options(mc.cores = parallel::detectCores())
library(Matrix)

# set the seed
set.seed(563478)

# change these to where the code is and where you want the plots saved
# globalpath <- "/storage/phd/projects/aghq-softwarepaper/paper/"
# plotpath <- paste0(globalpath,"figures/astro/")
globalpath <- tempdir()
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
pdf(file = file.path(plotpath,"stan-trace.pdf"),width = 7,height = 7)
traceplot(stanmod)
dev.off()

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


pdf(file.path(plotpath,"psi0-plot.pdf"),width = 7,height = 7)
hist(standata$psi0,freq=FALSE,breaks=50,main = "",xlab = "",cex.lab=1.5,cex.axis = 1.5)
with(psi0pdf,lines(transparam,pdf_transparam,lwd=2))
with(psi0pdf,lines(transparam,Psi0prior(transparam),lty = 'dashed',lwd=2))
with(tmbpsi0,lines(psi0,pdf,lty='dotdash',lwd=2))
dev.off()

pdf(file.path(plotpath,"gamma-plot.pdf"),width = 7,height = 7)
hist(standata$gamma,freq=FALSE,breaks=50,main = "",xlab = "",cex.lab=1.5,cex.axis = 1.5)
with(gammapdf,lines(transparam,pdf_transparam,lwd=2))
with(gammapdf,lines(transparam,gammaprior(transparam),lty = 'dashed',lwd=2))
with(tmbgamma,lines(gamma,pdf,lty='dotdash',lwd=2))
dev.off()

pdf(file.path(plotpath,"alpha-plot.pdf"),width = 7,height = 7)
hist(standata$alpha,freq=FALSE,breaks=50,main = "",xlab = "",cex.lab=1.5,cex.axis = 1.5)
with(alphapdf,lines(transparam,pdf_transparam,lwd=2))
with(alphapdf,lines(transparam,alphaprior(transparam),lty = 'dashed',lwd=2))
with(tmbalpha,lines(alpha,pdf,lty='dotdash',lwd=2))
dev.off()

pdf(file.path(plotpath,"beta-plot.pdf"),width = 7,height = 7)
hist(standata$beta,freq=FALSE,breaks=50,main = "",xlab = "",cex.lab=1.5,cex.axis = 1.5)
with(betapdf,lines(transparam,pdf_transparam,lwd=2))
with(betapdf,lines(transparam,betaprior(transparam),lty = 'dashed',lwd=2))
with(tmbbeta,lines(beta,pdf,lty='dotdash',lwd=2))
dev.off()

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


readr::write_csv(kstable,file.path(plotpath,"astro-ks-table.csv"))

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
pdf(file.path(plotpath,"massplot-aghq.pdf"),width=7,height=7)
plot(rtodo,Mrout,type='l',lwd=2,xlab="",ylab="",xaxt='n',cex.axis=1.5)
title(ylab=bquote('M(r) ('~10^12~M[sun]~')'),cex.lab=1.5,line=2.3)
axis(1,at=seq(0,150,by=25),cex.axis=1.5)
lines(rtodo,Mrout - 2*Mrsdout,lty='dashed',lwd=2)
lines(rtodo,Mrout + 2*Mrsdout,lty='dashed',lwd=2)
dev.off()

# With MCMC
Mrlist <- list()
for (i in 1:length(rtodo)) Mrlist[[i]] <- apply(standata[ ,c(1,2)],1,function(x) Mr(rtodo[i],x))
meanvals <- Reduce(c,Map(mean,Mrlist))
lowervals <- Reduce(c,Map(quantile,Mrlist,probs = .025))
uppervals <- Reduce(c,Map(quantile,Mrlist,probs = .975))


pdf(file.path(plotpath,"massplot-mcmc.pdf"),width=7,height=7)
plot(rtodo,meanvals,type='l',lwd=2,xlab="",ylab="",xaxt='n',cex.axis=1.5)
title(ylab=bquote('M(r) ('~10^12~M[sun]~')'),cex.lab=1.5,line=2.3)
axis(1,at=seq(0,150,by=25),cex.axis=1.5)
lines(rtodo,lowervals,lty='dashed',lwd=2)
lines(rtodo,uppervals,lty='dashed',lwd=2)
dev.off()

# Empirical RMSE
sqrt(mean( (Mrout - meanvals)^2 )) # 0.0004119291
sqrt(mean( ((Mrout - 2*Mrsdout) - lowervals)^2 )) # 0.007427573
sqrt(mean( ((Mrout + 2*Mrsdout) - uppervals)^2 )) # 0.006063218

