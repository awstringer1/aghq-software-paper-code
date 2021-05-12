### Basic Use example ###

# globalpath <- "/storage/phd/projects/aghq-softwarepaper/paper/"
# plotpath <- paste0(globalpath,"/figures/basic-use/")
globalpath <- tempdir()
plotpath <- file.path(globalpath,"basic-use")
if (!dir.exists(plotpath)) dir.create(plotpath)

library(aghq)

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
pdf(file = file.path(plotpath,'lambda-post-plot.pdf'))
with(pdfwithlambda,{
  hist(lambdapostsamps,breaks = 50,freq = FALSE,main = "",xlab = expression(lambda))
  lines(transparam,pdf_transparam)
  lines(transparam,dgamma(transparam,1+sum(y),1+length(y)),lty='dashed')
})
dev.off()

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