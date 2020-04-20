#'@name RWL
#'
#'@aliases  RWL
#'@aliases  dRWL
#'@aliases  pRWL
#'@aliases  qRWL
#'@aliases  rRWL
#'@aliases  hRWL
#'@aliases  plotRWL
#'@aliases  meanRWL
#'
#'@title Reparameterized Weighted Lindley (RWL) distribution for fitting a GAMLSS
#'
#'@description The fuction \code{RWL()} defines the RWL distribution, a two paramenter
#'distribution, for a gamlss.family object to be used in GAMLSS fitting using using the
#'function \code{gamlss()}, with mean equal to the parameter \code{mu} and \code{sigma}
#'equal the precision parameter. The functions \code{dRWL}, \code{pWL}, \code{qRWL},
#'\code{rRWL} and \code{hRWL} define the density, distribution function, quantile function, random
#'genetation and hazard function for the \code{RWL} parameterization of the RWL distribution.
#'
#'@usage RWL(mu.link = "log", sigma.link = "log")
#'
#' @param mu.link object for which the extraction of model residuals is meaningful.
#' @param sigma.link type of residual to be used.
#' @param x,q vector of quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log.d, log.p logical; if TRUE, density d are given as log(d).
#' @param log.p, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param log.q, log.p logical; if TRUE, quantiles q are given as log(q).
#' @param log.h, log.p logical; if TRUE, hazard h are given as log(h).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from.
#' @param to up to where to plot the distribution.
#' @param obj a fitted RWL object.
#' @param ... other graphical parameters for plotting.
#' @param title title of the plot.
#'
#'
#'@return returns a \code{gamlss.family} object which can be used to fit a reparameterized weighted Lindley distribution in the \code{gamlss()} function.
#'
#'@note For the function RWL(), mu is the mean and sigma is the precision parameter of the reparameterized weighted Lindley distribution.
#'
#'@author
#'Manoel Santos-Neto \email{santosnetoce at protemail.com}
#'
#'@importFrom gamlss.dist checklink
#'
#'@export

RWL <- function (mu.link = "log" , sigma.link="log")
{
  mstats <- checklink("mu.link", "RWL", substitute(mu.link),
                      c("sqrt","log","identity"))
  dstats <- checklink("sigma.link", "RWL", substitute(sigma.link)
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

                   part_1 <- (sigma+1)/b_mu_sigma
                   part_2 <- 1/(b_mu_sigma + 2*mu*sigma)
                   part_3 <- y/(2*mu)
                   part_4 <- (2*sigma*mu) - (y*b_mu_sigma)
                   part_4_p <- 1/(2*mu*mu)
                   part_5 <- (2*sigma)/(b_mu_sigma + (2*mu*sigma))

                   dldm <- A_i*(part_1 - part_2 - part_3) - part_4_p*part_4 - part_5

                   dldm
                 },
                 #the expected second derivative of the likelihood with respect to the location parameter mu
                 d2ldm2 = function(mu,sigma) {        #expected of second derivate of log-density respect to mu

                   A_i <- -sigma + ((sigma^2)*(mu-1) + 2*sigma*(sigma+1))/sqrt((sigma^2)*((mu-1)^2) + 4*mu*sigma*(sigma+1) )
                   b_mu_sigma <- (sigma*(1-mu) + sqrt(sigma*sigma*(mu-1)*(mu-1) + 4*mu*sigma*(sigma+1)))
                   dA_i <- (sigma^2)/sqrt((sigma^2)*((mu-1)^2) + 4*mu*sigma*(sigma+1) ) - ((sigma*sigma*(mu-1) + 2*sigma*(sigma+1))^2)/((sqrt((sigma^2)*((mu-1)^2) + 4*mu*sigma*(sigma+1) ))^3)

                   part_1 <- (sigma+1)*(dA_i*b_mu_sigma - A_i*A_i)
                   part_1_d <- b_mu_sigma*b_mu_sigma
                   part_2 <- sigma
                   part_2_d <- mu*mu
                   part_3 <- dA_i*(b_mu_sigma + 2*mu*sigma) - (A_i + 2*sigma)^2
                   part_3_d <- (b_mu_sigma + 2*mu*sigma)^2
                   part_4 <- dA_i*mu - A_i
                   part_4_d <- 2*mu
                   part_5 <- A_i*mu - 2*b_mu_sigma
                   part_5_d <- 2*mu*mu

                   d2ldm2 <-  part_1/part_1_d + part_2/part_2_d - part_3/part_3_d - part_4/part_4_d + part_5/part_5_d

                   d2ldm2
                 },

                 #the first derivative of the likelihood with respect to the scale parameter sigma
                 dldd = function(y,mu,sigma) {      #first derivate log-density respect to sigma

                   b_mu_sigma <- (sigma*(1-mu) + sqrt(sigma*sigma*(mu-1)*(mu-1) + 4*mu*sigma*(sigma+1)))
                   B_i <- 1 - mu + (sigma*((mu-1)^2) + (2*mu)*((2*sigma)+1))/sqrt(sigma*sigma*(mu-1)*(mu-1) + 4*mu*sigma*(sigma+1))

                   part_1 <- (sigma+1)/b_mu_sigma
                   part_2 <- 1/(b_mu_sigma + 2*mu*sigma)
                   part_3 <- y/(2*mu)
                   part_4 <- (2*mu)/(b_mu_sigma + (2*mu*sigma))

                   dldd  <- log(b_mu_sigma) + B_i*(part_1 - part_2 - part_3) - log(2*mu) - digamma(sigma) + log(y) - part_4


                   dldd
                 },

                 #the expected second derivative of the likelihood with respect to the scale parameter sigma ok
                 d2ldd2 = function(mu,sigma) {      #expected of second derivate log-density respect to sigma

                   b_mu_sigma <- (sigma*(1-mu) + sqrt(sigma*sigma*(mu-1)*(mu-1) + 4*mu*sigma*(sigma+1)))
                   B_i <- 1 - mu + (sigma*((mu-1)^2) + (2*mu)*((2*sigma)+1))/sqrt(sigma*sigma*(mu-1)*(mu-1) + 4*mu*sigma*(sigma+1))
                   dB_i <-  ((mu-1)^2 + 4*mu)/sqrt((sigma^2)*((mu-1)^2) + 4*mu*sigma*(sigma+1) ) - ((sigma*(mu-1)*(mu-1) + 2*mu*(2*sigma+1))^2)/((sqrt((sigma^2)*((mu-1)^2) + 4*mu*sigma*(sigma+1) ))^3)

                   part_1 <- (B_i + (sigma+1)*dB_i)*b_mu_sigma - (sigma+1)*B_i*B_i
                   part_1_d <- b_mu_sigma*b_mu_sigma
                   part_2 <- dB_i*(b_mu_sigma + 2*mu*sigma) - (B_i + 2*mu)^2
                   part_2_d <- (b_mu_sigma + 2*mu*sigma)^2


                   d2ldd2 <-  B_i/b_mu_sigma + part_1/part_1_d - trigamma(sigma) - dB_i/2 - part_2/part_2_d

                   d2ldd2
                 },

                 #the expected cross derivative of the likelihood with respect to both the location mu and scale parameter sigma
                 d2ldmdd = function(mu,sigma) {   #expected of partial derivate of log-density respect to mu and sigma

                   A_i <- -sigma + ((sigma^2)*(mu-1) + 2*sigma*(sigma+1))/sqrt((sigma^2)*((mu-1)^2) + 4*mu*sigma*(sigma+1) )
                   B_i <- 1 - mu + (sigma*((mu-1)^2) + (2*mu)*((2*sigma)+1))/sqrt(sigma*sigma*(mu-1)*(mu-1) + 4*mu*sigma*(sigma+1))
                   C_i <- -1 + (2*sigma*(mu-1) + 2*(2*sigma+1))/sqrt((sigma^2)*((mu-1)^2) + 4*mu*sigma*(sigma+1) ) - ( (sigma*(mu-1)*(mu-1) + 2*mu*(2*sigma+1))*(sigma*sigma*(mu-1) + 2*sigma*(sigma+1)))/((sqrt((sigma^2)*((mu-1)^2) + 4*mu*sigma*(sigma+1) ))^3)
                   b_mu_sigma <- (sigma*(1-mu) + sqrt(sigma*sigma*(mu-1)*(mu-1) + 4*mu*sigma*(sigma+1)))

                   part_1 <- (sigma+1)*(C_i*b_mu_sigma - B_i*A_i)
                   part_1_d <- b_mu_sigma*b_mu_sigma
                   part_2 <- C_i*mu - B_i
                   part_2_d <- 2*mu
                   part_3 <- (C_i +2)*(b_mu_sigma + 2*mu*sigma) - (A_i + 2*sigma)*(B_i + 2*mu)
                   part_3_d <- (b_mu_sigma + 2*mu*sigma)^2

                   d2ldmdd <- A_i/b_mu_sigma + part_1/part_1_d - 1/mu - part_2/part_2_d - part_3/part_3_d

                   d2ldmdd
                 },

                 G.dev.incr = function(y,mu,sigma,...) -2*dRWL(y,mu,sigma,log.d = TRUE),
                 rqres = expression(rqres(pfun = "pRWL", type = "Continuous", y = y, mu = mu, sigma = sigma)),
                 mu.initial = expression({mu <- (y+mean(y))  }),
                 sigma.initial = expression({sigma <- rep(1,length(y)) }),
                 mu.valid = function(mu) all(mu>0) ,
                 sigma.valid = function(sigma) all(sigma > 0),
                 y.valid = function(y) all(y > 0),
                 mean = function(mu,sigma) mu,
                 variance = function(mu, sigma){
                   b_mu_sigma <- (sigma*(1-mu) + sqrt(sigma*sigma*(mu-1)*(mu-1) + 4*mu*sigma*(sigma+1)))

                   part_1 <- (2*mu)^2
                   part_1_d <- b_mu_sigma^2
                   part_2 <- (b_mu_sigma + 2*mu*(sigma+2))*(sigma*(sigma+1))
                   part_2_d <- b_mu_sigma + 2*mu*sigma

                   v_y  <- (part_1/part_1_d)*(part_2/part_2_d) - mu^2

                  v_y
                 }
                ),
            class = c("gamlss.family","family"))
}

#'@rdname RWL
#'@importFrom LindleyR dwlindley
#'
#'@export


dRWL <- function(x,mu=1,sigma=1,log.d=FALSE)
{
  theta <- (sigma*(1-mu) + sqrt(sigma*sigma*(mu-1)*(mu-1) + 4*mu*sigma*(sigma+1)))/(2*mu)

  dens <- dwlindley(x=x, theta=theta,alpha=sigma,log = log.d)

  dens
}


#'@rdname RWL
#'@importFrom LindleyR pwlindley
#'
#'@export


pRWL <- function(q,mu=1,sigma=1,lower.tail=TRUE, log.p=FALSE)
{
  theta <- (sigma*(1-mu) + sqrt(sigma*sigma*(mu-1)*(mu-1) + 4*mu*sigma*(sigma+1)))/(2*mu)
  cum <- pwlindley(q,theta=theta,alpha=sigma,lower.tail = lower.tail,log.p = log.p)

  cum
}

#'@rdname RWL
#'@importFrom LindleyR qwlindley
#'
#'@export


qRWL <- function(p, mu=1, sigma=1, lower.tail=TRUE, log.q=FALSE)
{

theta <- (sigma*(1-mu) + sqrt(sigma*sigma*(mu-1)*(mu-1) + 4*mu*sigma*(sigma+1)))/(2*mu)
quant <- qwlindley(p,theta=theta,alpha=sigma,lower.tail = lower.tail, log.p = log.q)

quant

}

#'@rdname RWL
#'@importFrom LindleyR rwlindley
#'
#'@export

rRWL <- function(n, mu=1, sigma=1)
{
  theta <- (sigma*(1-mu) + sqrt(sigma*sigma*(mu-1)*(mu-1) + 4*mu*sigma*(sigma+1)))/(2*mu)
  rand <- rwlindley(n,theta=theta,alpha=sigma)

  rand
}

#'@rdname RWL
#'@importFrom LindleyR hwlindley
#'
#'@export

hRWL <- function(x, mu=1, sigma=1, log.h = FALSE)
{

  theta <- (sigma*(1-mu) + sqrt(sigma*sigma*(mu-1)*(mu-1) + 4*mu*sigma*(sigma+1)))/(2*mu)
  hf <- hwlindley(x=x,theta=theta,alpha=sigma,log=log.h)

  hf
}


#'@rdname RWL
#'@importFrom graphics plot
#'
#'@export
plotRWL = function(mu = .5, sigma = 1, from = 0, to = 0.999, n = 101, title="title", ...)
{
  y = seq(from = 0.001, to = to, length.out = n)
  pdf = dRWL(y, mu = mu, sigma = sigma)
  plot(pdf ~ y, main=title, ylim = c(0, max(pdf)), type = "l",lwd=3)

}

#'@rdname RWL
#'@importFrom stats fitted
#'
#'@export
meanRWL = function (obj)
{
  if (obj$family[1] != "RWL")
    stop("the object do not have a RWL distribution")
  meanofY = fitted(obj, "mu")
  meanofY
}




