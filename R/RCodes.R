#' Reparameterized Weighted Lindley Distribution
#' 
#'@description Density function, cumulative distribution function, quantile function, random number generation and hazard function of the reparameterized weighted Lindley distribution.
#'
#'@usage 
#'drwl(x,mu,phi,log.d=FALSE)
#'prwl(q,mu,phi,lower.tail=TRUE, log.p=FALSE) 
#'qrwl(p, mu, phi, lower.tail=TRUE, log.q=FALSE)
#'rrWL(n, mu, phi)
#'hrwl(x, mu, phi, log.h = FALSE)
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


library(LindleyR)

drwl <- function(x,mu,phi,log.d=FALSE)
{
theta <- (phi*(1-mu) + sqrt(phi*phi(mu-1)*(mu-1) + 4*mu*phi*(phi+1)))/(2*mu)
dens <- dwlindley(x=x, theta=theta,alpha=phi,log = log.d)    
  
dens
}

prwl <- function(q,mu,phi,lower.tail=TRUE, log.p=FALSE)
{
  theta <- (phi*(1-mu) + sqrt(phi*phi(mu-1)*(mu-1) + 4*mu*phi*(phi+1)))/(2*mu)
  cum <- pwlindley(q,theta=theta,alpha=phi,lower.tail = lower.tail,log.p = log.p)  
  
  cum
}

qrwl <- function(p, mu, phi, lower.tail=TRUE, log.q=FALSE)
{
  theta <- (phi*(1-mu) + sqrt(phi*phi(mu-1)*(mu-1) + 4*mu*phi*(phi+1)))/(2*mu)
  quant <- qwlindley(q,theta=theta,alpha=phi,lower.tail = lower.tail, log.p = log.q)  
  
  quant
  
}


rrWL <- function(n, mu, phi)
{
theta <- (phi*(1-mu) + sqrt(phi*phi(mu-1)*(mu-1) + 4*mu*phi*(phi+1)))/(2*mu)
rand <- rwlindley(n,theta=theta,alpha=phi)  
  
rand  
}  


hrwl <- function(x, mu, phi, log.h = FALSE)
{

  theta <- (phi*(1-mu) + sqrt(phi*phi(mu-1)*(mu-1) + 4*mu*phi*(phi+1)))/(2*mu)
  hf <- hwlindley(n,theta=theta,alpha=phi,log=log.h)  
  
  hf   
  
  
}  
