#' Reparameterized Weighted Lindley Distribution
#'
#'@description Density function, cumulative distribution function, quantile function, random number generation and hazard function of the reparameterized weighted Lindley distribution.
#'
#'@usage
#'dRWL(x,mu,sigma,log.d=FALSE)
#'pRWL(q,mu,sigma,lower.tail=TRUE, log.p=FALSE)
#'qRWL(p, mu, sigma, lower.tail=TRUE, log.q=FALSE)
#'rRWL(n, mu, sigma)
#'hRWL(x, mu, sigma, log.h = FALSE)
#'
#'@param x,q vector of quantiles.
#'@param mu,sigma positive parameters.
#'@param log.d,log.p,log.q,log.h logical; If TRUE are given as log().
#'@param lower.tail logical; If TRUE, (default), P(X ≤ x) are returned, otherwise P(X > x).
#'@param p vector of probabilities.
#'@param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#'
#'
#'
#'@import LindleyR
#'
#'@export

dRWL <- function(x,mu,sigma,log.d=FALSE)
{
  theta <- (sigma*(1-mu) + sqrt(sigma*sigma*(mu-1)*(mu-1) + 4*mu*sigma*(sigma+1)))/(2*mu)
  dens <- dwlindley(x=x, theta=theta,alpha=sigma,log = log.d)

  dens
}


#' Reparameterized Weighted Lindley Distribution
#'
#'@description Density function, cumulative distribution function, quantile function, random number generation and hazard function of the reparameterized weighted Lindley distribution.
#'
#'@usage
#'dRWL(x,mu,sigma,log.d=FALSE)
#'pRWL(q,mu,sigma,lower.tail=TRUE, log.p=FALSE)
#'qRWL(p, mu, sigma, lower.tail=TRUE, log.q=FALSE)
#'rRWL(n, mu, sigma)
#'hRWL(x, mu, sigma, log.h = FALSE)
#'
#'@param x,q vector of quantiles.
#'@param mu,sigma positive parameters.
#'@param log.d,log.p,log.q,log.h logical; If TRUE are given as log().
#'@param lower.tail logical; If TRUE, (default), P(X ≤ x) are returned, otherwise P(X > x).
#'@param p vector of probabilities.
#'@param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#'
#'
#'
#'@import LindleyR
#'
#'@export

pRWL <- function(q,mu,sigma,lower.tail=TRUE, log.p=FALSE)
{
  theta <- (sigma*(1-mu) + sqrt(sigma*sigma*(mu-1)*(mu-1) + 4*mu*sigma*(sigma+1)))/(2*mu)
  cum <- pwlindley(q,theta=theta,alpha=sigma,lower.tail = lower.tail,log.p = log.p)

  cum
}

#' Reparameterized Weighted Lindley Distribution
#'
#'@description Density function, cumulative distribution function, quantile function, random number generation and hazard function of the reparameterized weighted Lindley distribution.
#'
#'@usage
#'dRWL(x,mu,sigma,log.d=FALSE)
#'pRWL(q,mu,sigma,lower.tail=TRUE, log.p=FALSE)
#'qRWL(p, mu, sigma, lower.tail=TRUE, log.q=FALSE)
#'rRWL(n, mu, sigma)
#'hRWL(x, mu, sigma, log.h = FALSE)
#'
#'@param x,q vector of quantiles.
#'@param mu,sigma positive parameters.
#'@param log.d,log.p,log.q,log.h logical; If TRUE are given as log().
#'@param lower.tail logical; If TRUE, (default), P(X ≤ x) are returned, otherwise P(X > x).
#'@param p vector of probabilities.
#'@param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#'
#'
#'
#'@import LindleyR
#'
#'@export

qRWL <- function(p, mu, sigma, lower.tail=TRUE, log.q=FALSE)
{

theta <- (sigma*(1-mu) + sqrt(sigma*sigma*(mu-1)*(mu-1) + 4*mu*sigma*(sigma+1)))/(2*mu)
quant <- qwlindley(q,theta=theta,alpha=phi,lower.tail = lower.tail, log.p = log.q)

quant

}

#' Reparameterized Weighted Lindley Distribution
#'
#'@description Density function, cumulative distribution function, quantile function, random number generation and hazard function of the reparameterized weighted Lindley distribution.
#'
#'@usage
#'dRWL(x,mu,sigma,log.d=FALSE)
#'pRWL(q,mu,sigma,lower.tail=TRUE, log.p=FALSE)
#'qRWL(p, mu, sigma, lower.tail=TRUE, log.q=FALSE)
#'rRWL(n, mu, sigma)
#'hRWL(x, mu, sigma, log.h = FALSE)
#'
#'@param x,q vector of quantiles.
#'@param mu,sigma positive parameters.
#'@param log.d,log.p,log.q,log.h logical; If TRUE are given as log().
#'@param lower.tail logical; If TRUE, (default), P(X ≤ x) are returned, otherwise P(X > x).
#'@param p vector of probabilities.
#'@param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#'
#'
#'
#'@import LindleyR
#'
#'@export

rRWL <- function(n, mu, sigma)
{
  theta <- (sigma*(1-mu) + sqrt(sigma*sigma*(mu-1)*(mu-1) + 4*mu*sigma*(sigma+1)))/(2*mu)
  rand <- rwlindley(n,theta=theta,alpha=sigma)

  rand
}

#' Reparameterized Weighted Lindley Distribution
#'
#'@description Density function, cumulative distribution function, quantile function, random number generation and hazard function of the reparameterized weighted Lindley distribution.
#'
#'@usage
#'dRWL(x,mu,phi,log.d=FALSE)
#'pRWL(q,mu,phi,lower.tail=TRUE, log.p=FALSE)
#'qRWL(p, mu, phi, lower.tail=TRUE, log.q=FALSE)
#'rRWL(n, mu, phi)
#'hRWL(x, mu, phi, log.h = FALSE)
#'
#'@param x,q vector of quantiles.
#'@param mu,phi positive parameters.
#'@param log.d,log.p,log.q,log.h logical; If TRUE are given as log().
#'@param lower.tail logical; If TRUE, (default), P(X ≤q x) are returned, otherwise P(X > x).
#'@param p vector of probabilities.
#'@param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#'
#'
#'
#'@import LindleyR
#'
#'@export

hRWL <- function(x, mu, sigma, log.h = FALSE)
{

  theta <- (sigma*(1-mu) + sqrt(sigma*sigma*(mu-1)*(mu-1) + 4*mu*sigma*(sigma+1)))/(2*mu)
  hf <- hwlindley(n,theta=theta,alpha=sigma,log=log.h)

  hf
}

#'Reparameterized Weighted Lindley (RWL) distribution for fitting a GAMLSS
#'
#'@description The fuction \code{RWL()} defines the BS distribution, a two paramenter
#'distribution, for a gamlss.family object to be used in GAMLSS fitting using using the
#'function \code{gamlss()}, with mean equal to the parameter \code{mu} and \code{sigma}
#'equal the precision parameter. The functions \code{dRWL}, \code{pWL}, \code{qRWL} and
#'\code{rRWL} define the density, distribution function, quantile function and random
#'genetation for the \code{RWL} parameterization of the RWL distribution.
#'
#'@usage RWL(mu.link = "identity", sigma.link = "identity")
#'dRWL(x, mu = 1, sigma = 1, log = FALSE)
#'pRWL(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'qRWL(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'rRWL(n, mu = 1, sigma = 1)
#'plotRWL(mu = .5, sigma = 1, from = 0, to = 0.999, n = 101, ...)
#'meanRWL(obj)
#'
#' @param mu.link object for which the extraction of model residuals is meaningful.
#' @param sigma.link type of residual to be used.
#' @param x,q vector of quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from
#' @param to up to where to plot the distribution
#' @param obj a fitted RBS object
#' @param ... other graphical parameters for plotting
#'
#'
#'@return returns a \code{gamlss.family} object which can be used to fit a normal distribution in the \code{gamlss()} function.
#'
#'@note For the function RWL(), mu is the mean and sigma is the precision parameter of the reparameterized weighted Lindley distribution.
#'
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}
#'
#'@export
#'
RWL <- function (mu.link = "log" , sigma.link="log")
{
  mstats = checklink("mu.link", "Reparameterized Weighted Lindley", substitute(mu.link),
                     c("sqrt","log","identity"))
  dstats = checklink("sigma.link", "Reparameterized Weighted Lindley", substitute(sigma.link)
                     ,c("sqrt", "log", "identity"))
  structure(list(family = c("RWL","Reparameterized Weighted Lindley"),
                 parameters = list(mu=TRUE,sigma=TRUE),
                 nopar = 2,
                 type = "Continuous",
                 mu.link = as.character(substitute(mu.link)),
                 sigma.link = as.character(substitute(sigma.link)),
                 mu.linkfun = mstats$linkfun,
                 sigma.linkfun = dstats$linkfun,
                 mu.linkinv = mstats$linkinv,
                 sigma.linkinv = dstats$linkinv,
                 mu.dr = mstats$mu.eta,
                 sigma.dr = dstats$mu.eta,
                 #the first derivative of the likelihood with respect to the location parameter mu
                 dldm = function(y,mu,sigma) #first derivate of log-density respect to mu
                 {

                  A_i <- -sigma + ((sigma^2)*(mu-1) + 2*sigma*(sigma+1))/sqrt((sigma^2)*((mu-1)^2) + 4*mu*sigma*(sigma+1) )
                  b_mu_sigma <- (sigma*(1-mu) + sqrt(sigma*sigma*(mu-1)*(mu-1) + 4*mu*sigma*(sigma+1)))
                   mustart = 1/(2*mu)
                   ystart =  ((sigma+1)*y)/(4*mu*mu) - (sigma^2)/(4*(sigma+1)*y) + sigma/((sigma*y) + y + (sigma*mu))

                   dldm = ystart-mustart
                   dldm
                 },
                 #the expected second derivative of the likelihood with respect to the location parameter mu
                 d2ldm2 = function(mu,sigma) {        #expected of second derivate of log-density respect to mu
                   d2ldm2 =  - sigma/(2*mu*mu) - ((sigma/(sigma+1))^2)*Ims(mu,sigma)
                   d2ldm2
                 },

                 #the first derivative of the likelihood with respect to the scale parameter sigma
                 dldd = function(y,mu,sigma) {      #first derivate log-density respect to sigma
                   sigmastart  = -(sigma)/(2*(sigma+1))
                   y2start   = (y+mu)/((sigma*y) + y + (sigma*mu)) - y/(4*mu) - (mu*sigma*(sigma+2))/(4*y*((sigma+1)^2))
                   dldd  = y2start-sigmastart
                   dldd
                 },
                 #the expected second derivative of the likelihood with respect to the scale parameter sigma ok
                 d2ldd2 = function(mu,sigma) {      #expected of second derivate log-density respect to sigma
                   lss =  ((sigma^2) + (3*sigma) + 1)/(2*sigma*sigma*((sigma+1)^2))
                   d2ldd2 = -lss - ((mu^2)/((sigma+1)^4))*Ims(mu,sigma)
                   d2ldd2
                 },
                 #the expected cross derivative of the likelihood with respect to both the location mu and scale parameter sigma
                 d2ldmdd = function(mu,sigma) {   #expected of partial derivate of log-density respect to mu and sigma
                   lms = 1/(2*mu*(sigma+1))
                   d2ldmdd = - lms - ((mu*sigma)/((sigma+1)^3))*Ims(mu,sigma)
                   d2ldmdd
                 },

                 G.dev.incr = function(y,mu,sigma,...) -2*dRBS(y,mu,sigma,log=TRUE),
                 rqres = expression(rqres(pfun = "pRBS", type = "Continuous", y = y, mu = mu, sigma = sigma)),
                 mu.initial = expression({mu <- y + mean(y)/2 }),
                 sigma.initial = expression({sigma <- rep(1,length(y)) }),
                 mu.valid = function(mu) all(mu>0) ,
                 sigma.valid = function(sigma) all(sigma > 0),
                 y.valid = function(y) all(y > 0),
                 mean = function(mu,sigma) mu,
                 variance = function(mu, sigma) (mu*mu)*((2*sigma+5)/((sigma+1)^2)) ),
            class = c("gamlss.family","family"))
}



