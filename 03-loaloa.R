### Loa-loa example: not zero-inflated ###
# Comparison with INLA, no MCMC

library(geostatsp)
library(aghq)
library(PrevMap)
library(geoR)

ilogit <- function(x) 1 / (1 + exp(-x))

set.seed(78968)

# Set flags
# globalpath <- "/storage/phd/projects/aghq-softwarepaper/paper/"
# savepath <- paste0(globalpath,"data/")
savestamp <- "20210504-v1"
# plotpath <- paste0(globalpath,"figures/loaloa/")
globalpath <- tempdir()
plotpath <- file.path(globalpath,"loaloa")
if (!dir.exists(plotpath)) dir.create(plotpath)
savepath <- plotpath



# Flags, which analysis to do?
doaghq <- TRUE
doinla <- TRUE
domcml <- TRUE
domcmc <- TRUE
dopostsamplingaghq <- TRUE
doplotting <- TRUE


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

if (doaghq) {
  tm <- Sys.time()
  cat("Doing AGHQ, time = ",format(tm),"\n")
  loaloaquad <- aghq::marginal_laplace(
    ff,
    5,
    list(W = rep(0,Wdim),theta = startingtheta)
  )
  aghqtime <- difftime(Sys.time(),tm,units = 'secs')
  save(loaloaquad,file = file.path(savepath,paste0("loaloaquad",savestamp,".RData")))
  cat("AGHQ took: ",format(aghqtime),"\n")
}


if (dopostsamplingaghq){
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
  raster::writeRaster(simbrick,file = file.path(savepath,paste0("loaloa-simbrick-aghq",savestamp,".grd")),overwrite = TRUE)
  cat("AGHQ simulation took: ",format(aghqsimtime),"\n")
}

## INLA ----
# Using geostatsp::glgm()

if (doinla) {
  tm <- Sys.time()
  cat("Doing INLA, time = ",format(tm),"\n")
  loaFit = glgm(formula = y ~ 1,
                data = loaloa, grid = 50,
                family = "binomial",
                Ntrials = loaloa$N, shape = maternconstants$nu, buffer = 1e05,
                prior = list(sd = c(u = sigma_u, alpha = sigma_alpha), range = c(u = rho_u, alpha = rho_alpha))) 
  inlatime <- difftime(Sys.time(),tm,units = 'secs')
  save(loaFit,file = file.path(savepath,paste0("loaloainla",savestamp,".RData")))
  cat("INLA took: ",format(inlatime),"\n")
}

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
if (domcml) {
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
  
  fit.MCML1 <- binomial.logistic.MCML(formula = y ~ 1,units.m = ~ N,
                                      par0 = par0[1:3],
                                      coords = ~ long + lat, data = as.data.frame(loaloa),
                                      control.mcmc = c.mcmc,
                                      kappa = maternconstants$nu, 
                                      fixed.rel.nugget = 0,
                                      start.cov.pars = c(1e04))
  mcmltime <- difftime(Sys.time(),tm,units = 'secs')
  cat("MCML took: ",format(mcmltime),"\n")
  
  tm <- Sys.time()
  cat("Doing MCML Sims, time = ",format(tm),"\n")
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
  
  # rasterize the grid
  mcmcl_predraster <- raster(loaFit$raster$random.mean)
  values(mcmcl_predraster) <- preds_mcml$prevalence$predictions
  mcmlsimtime <- difftime(Sys.time(),tm,units = 'secs')
  save(fit.MCML1,preds_mcml,file = file.path(savepath,paste0("loaloamcml",savestamp,".RData")))
  raster::writeRaster(mcmcl_predraster,file.path(savepath,paste0("loaloamcml-preds",savestamp,".grd")),overwrite = TRUE)
  cat("MCML Sims took: ",format(mcmlsimtime),"\n")
}
## MCMC ----
if (domcmc) {
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
  
  fit.Bayes <- binomial.logistic.Bayes(
    formula = y ~ 1,
    units.m = ~ N,
    coords = ~ long + lat,
    data = as.data.frame(loaloa), control.prior = thepriors,
    control.mcmc = mcmc.Bayes, kappa = maternconstants$nu)
  summary(fit.Bayes, hpd.coverage = 0.95)
  
  mcmctime <- difftime(Sys.time(),tm,units = 'secs')
  cat("MCMC took: ",format(mcmctime),"\n")
  
  tm <- Sys.time()
  cat("Doing MCMC Sims, time = ",format(tm),"\n")
  
  pred.Bayes <- spatial.pred.binomial.Bayes(
    fit.Bayes, 
    as.data.frame(raster(loaFit$raster$random.mean),xy=TRUE),
    type = "marginal",
    scale.predictions = "prevalence", quantiles = NULL,
    standard.errors = TRUE, thresholds = 0.2,
    scale.thresholds = "prevalence")
  
  bayesian_predraster <- raster(loaFit$raster$random.mean)
  values(bayesian_predraster) <- pred.Bayes$prevalence$predictions
  
  mcmcsimtime <- difftime(Sys.time(),tm,units = 'secs')
  save(fit.Bayes,pred.Bayes,file = file.path(savepath,paste0("loaloamcmc",savestamp,".RData")))
  raster::writeRaster(bayesian_predraster,file.path(savepath,paste0("loaloamcmc-preds",savestamp,".grd")),overwrite = TRUE)
  cat("MCMC Sims took: ",format(mcmcsimtime),"\n")
}


# Plot them!

if (FALSE) {
  # Option to load the result from disk
  
  # AGHQ
  load(file.path(savepath,paste0("loaloaquad",savestamp,".RData")))
  simbrick <- raster::brick(file.path(savepath,paste0("loaloa-simbrick-aghq",savestamp,".grd")),overwrite = TRUE)
  # INLA
  load(file.path(savepath,paste0("loaloainla",savestamp,".RData")))
  # MCML
  load(file.path(savepath,paste0("loaloamcml",savestamp,".RData")))
  mcmcl_predraster <- raster::raster(file.path(savepath,paste0("loaloamcml-preds",savestamp,".grd")),overwrite = TRUE)
  # MCMC
  load(file.path(savepath,paste0("loaloamcmc",savestamp,".RData")))
  bayesian_predraster <- raster::raster(file.path(savepath,paste0("loaloamcmc-preds",savestamp,".grd")),overwrite = TRUE)
}

if (doplotting) {
  br <- c(0,.05,.1,.15,.2,.25,.3,.4,.5,.6,1)
  
  aghqpostmean <- calc(ilogit(simbrick),mean)
  
  png(file.path(plotpath,paste0("aghq-mean-map",savestamp,".png")))
  plot_loaloa(aghqpostmean,br)
  dev.off()
  
  png(file.path(plotpath,paste0("mcml-mean-map",savestamp,".png")))
  plot_loaloa(mcmcl_predraster,br)
  dev.off()
  
  png(file.path(plotpath,paste0("inla-mean-map",savestamp,".png")))
  plot_loaloa(loaFit$raster$predict.invlogit,br)
  dev.off()
  
  png(file.path(plotpath,paste0("mcmc-mean-map",savestamp,".png")))
  plot_loaloa(bayesian_predraster,br)
  dev.off()
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
  readr::write_csv(difftable,file.path(savepath,paste0("diff-table-loaloa",savestamp,".csv")))
  
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
  readr::write_csv(covparamtable,file.path(savepath,paste0("covparam-table-loaloa",savestamp,".csv")))
}

# Write the timing table
timingtable <- data.frame(
  task = c("AGHQ","INLA","MCML","MCMC","AGHQSim","MCMLSim","MCMCSim"),
  time = c(aghqtime,inlatime,mcmltime,mcmctime,aghqsimtime,mcmlsimtime,mcmcsimtime)
)
readr::write_csv(timingtable,file.path(savepath,paste0("timing-table-loaloa",savestamp,".csv")))
# Total and num iter
totaltable <- c(
  mcmc = as.numeric(timingtable[timingtable$task == 'MCMC','time']) + as.numeric(timingtable[timingtable$task == 'MCMCSim','time']),
  aghq = as.numeric(timingtable[timingtable$task == 'AGHQ','time']) + as.numeric(timingtable[timingtable$task == 'AGHQSim','time']),
  inla = as.numeric(timingtable[timingtable$task == 'INLA','time']),
  mcml = as.numeric(timingtable[timingtable$task == 'MCML','time']) + as.numeric(timingtable[timingtable$task == 'MCMLSim','time'])
)
effiter <- totaltable[2:4] * 6000 / totaltable[1]


cat("All done!\n")