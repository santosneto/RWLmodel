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

                   A_i <- -sigma + ((sigma^2)*(mu - 1) + 2*sigma*(sigma+1))/sqrt((sigma^2)*((mu-1)^2) + 4*mu*sigma*(sigma+1) )
                   b_mu_sigma <- (sigma*(1 - mu) + sqrt(sigma*sigma*(mu-1)*(mu-1) + 4*mu*sigma*(sigma+1)))
                   dA_i <- (sigma^2)/sqrt((sigma^2)*((mu-1)^2) + 4*mu*sigma*(sigma+1) ) - ((sigma*sigma*(mu-1) + 2*sigma*(sigma+1))^2)/((sqrt((sigma^2)*((mu-1)^2) + 4*mu*sigma*(sigma+1) ))^3)

                   part_1 <- (sigma + 1)*(dA_i*b_mu_sigma - A_i*A_i)
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
                 mu.initial = expression({mu <- (y + mean(y)/2)  }),
                 sigma.initial = expression({sigma <- rep(sd(y)/((mean(y))^1.5), length(y)) }),
                 mu.valid = function(mu) all(mu > 0) ,
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


dRWL <- function(x, mu = 1, sigma = 1, log.d = FALSE)
{
  theta <- (sigma*(1 - mu) + sqrt(sigma*sigma*(mu - 1)*(mu - 1) + 4*mu*sigma*(sigma + 1)))/(2*mu)

  dens <- dwlindley(x = x, theta = theta, alpha = sigma,log = log.d)

  dens
}


#'@rdname RWL
#'@importFrom LindleyR pwlindley
#'
#'@export


pRWL <- function(q, mu=1, sigma = 1,lower.tail = TRUE, log.p = FALSE)
{
  theta <- (sigma*(1 - mu) + sqrt(sigma*sigma*(mu - 1)*(mu - 1) + 4*mu*sigma*(sigma + 1)))/(2*mu)
  cum <- pwlindley(q,theta = theta,alpha = sigma,lower.tail = lower.tail,log.p = log.p)

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

#'@name score_test_rwl
#'
#'@aliases score_test_rwl
#'
#'@title Precision test
#'
#'@description Tests the null hypothesis of precision fixed in RBS models against the alternative of precision variable.
#'
#'@usage score_test_rwl(modelh0,modelh1)
#'
#' @param modelh0 model under null hypothesis.
#' @param modelh1 model under alternative hypothesis.
#'
#' @return A list with class "htest" containing the following components:
#' @return \code{statistic}	the value of the test statistic.
#' @return \code{parameter}	the degrees of freedom for the test statistic.
#' @return \code{p.value}	the p-value for the test.
#' @return \code{method}	a character string indicating what type of likelihood ratio test was performed.
#' @return \code{data.name} a character string giving the name(s) of the data
#'
#'@author
#'Manoel Santos-Neto \email{manoelferreira@uaest.ufcg.edu.br}
#'
#'
#'@importFrom stats pchisq
#'@importFrom stats pchisq
#'@importFrom Deriv Deriv
#'@export

score_test_rwl <- function(modelh0, modelh1)
{

    UalphaH0.gam <- function(modelh0, modelh1){
    phi_linkstr <- modelh1$sigma.link
    phi_linkobj <- make.link(phi_linkstr)
    phi_mu.eta <- phi_linkobj$mu.eta
    tau <- phi_linkobj$linkfun
    tauH0 <- tau(modelh0$sigma.fv)
    muH0 <- modelh0$mu.fv
    deltaH0 <- modelh0$sigma.fv
    z <- modelh1$sigma.x
    vt <- modelh0$y

    LL <- function(theta,y){
      L <- dRWL(x = y,mu = theta[1],sigma = theta[2],log.d = TRUE)
      return(L)
    }


    U <- matrix(NA,nrow = length(vt),2)
    for (i in 1:length(vt))
    {
      U[i,] <- pracma::grad(LL,x0 = cbind(muH0,deltaH0)[i,],y = vt[i])
    }


    dH0 <- U[,2]
    rval <- dH0 * phi_mu.eta(tauH0) * z
    colSums(rval)
  }

  vcovH0 <- function(modelh0, modelh1){
    mu <- modelh0$mu.fv
    sigma <-  modelh0$sigma.fv
    x <- modelh1$mu.x
    z <- modelh1$sigma.x
    linkstr <- modelh0$mu.link
    linkobj <- make.link(linkstr)
    mu.eta <- linkobj$mu.eta
    eta <- linkobj$linkfun
    etaH0 <- eta(modelh0$mu.fv)
    phi_linkstr <- modelh1$sigma.link
    phi_linkobj <- make.link(phi_linkstr)
    phi_mu.eta <- phi_linkobj$mu.eta
    tau <- phi_linkobj$linkfun
    tauH0 <- tau(modelh0$sigma.fv)
    y <- modelh0$y

    dai <- function(link)
    {
      switch(link, log = {mu.eta.2 <- function(eta) rep.int(1, length(eta))},
             identity = {mu.eta.2 <- function(eta) rep.int(0, length(eta))},
             sqrt = {mu.eta.2 <- function(eta) 1/eta} )
    }

    mu.eta.2    <- dai(linkstr)
    sigma.eta.2 <- dai(phi_linkstr)
    d_ai        <- mu.eta.2(etaH0)
    d_bi        <- sigma.eta.2(tauH0)
    ai          <- mu.eta(etaH0)
    bi          <- phi_mu.eta(tauH0)

    ll       <- function(y,mu,sigma){
      sigma2 <- sigma^2
      a_us   <- sigma*(1 - mu) + sqrt(sigma2*((mu - 1)^2) + 4*mu*sigma*(sigma + 1))
      fy     <- (sigma + 1)*log(a_us) + (sigma - 1)*log(y) + log(1 + y) - (a_us/(2*mu))*y - sigma*log(2*mu) - log(a_us + 2*mu*sigma) - lgamma(sigma)

      fy
    }


    dm   <- Deriv::Deriv(ll,'mu')
    ds   <- Deriv::Deriv(ll,'sigma')
    dmm  <- Deriv::Deriv(Deriv(ll,'mu'),'mu')
    dss  <- Deriv::Deriv(Deriv(ll,'sigma'),'sigma')
    dms  <- Deriv::Deriv(Deriv(ll,'mu'),'sigma')

    d_mu    <- dm(y,mu,sigma)
    d_sigma <- ds(y,mu,sigma)
    v       <- dmm(y,mu,sigma)
    u       <- dss(y,mu,sigma)
    s       <- dms(y,mu,sigma)

    ci   <- v*(ai^2) + d_mu*ai*d_ai
    mi   <- s*ai*bi
    wi   <- u*(bi^2) + d_sigma*bi*d_bi

    kbb  <- crossprod(ci * x, x)
    kaa  <- crossprod(wi * z, z)
    kba  <- crossprod(mi * x, z)
    hess <- cbind(rbind(kbb, t(kba)), rbind(kba, kaa))
    vcov <- solve(-hess)
    return(vcov)
  }

  METHOD <- "Rao score test"
  DNAME <- deparse(substitute(modelh0) )
  DNAME <- paste(DNAME, "vs", deparse(substitute(modelh1)))
  p <- modelh0$mu.df
  p1 <- p + 1
  varalpha <- vcovH0(modelh0, modelh1)[-(1:p1), -(1:p1)]
  Ua. <- UalphaH0.gam(modelh0, modelh1)[-1]
  SC <- t(Ua.) %*% varalpha %*% Ua.
  q <- modelh1$sigma.df
  gl <- q - 1
  PVAL <- pchisq(SC, q - 1, lower.tail = F)
  names(gl) <- "df"
  names(SC) <- "SC"
  RVAL <- list(statistic = SC, parameter = gl, p.value = PVAL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}


#'@name envelope.RWL
#'
#'@aliases envelope.RWL
#'
#'@title Envelopes
#'
#'@description Computes simulation envelopes.
#'
#' @param model object of class \code{gamlss}.
#' @param k number of replications for envelope construction. Default is 19.
#' @param xlabel a label for the x axis.
#' @param ylabel a label for the y axis.
#' @param color a specification for the default envelope color.
#' @param font the name of a font family for x and y axis.
#'
#'
#'@return A simulated envelope of the class BP, GA, IG, RBS and WEI3.
#'
#'@author
#'Manoel Santos-Neto \email{manoelferreira@uaest.ufcg.edu.br}
#'
#'@references
#'Atkinson, A. C. (1985) Plots, transformations and regression : an introduction to graphical methods of diagnostic regression analysis. Oxford Science Publications, Oxford.
#'Bourguignon, M., Santos-Neto, M. and Castro, M. (2018). A new regression model for positive data. arXiv.
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A, Barros, M. (2014)  Birnbaum-Saunders statistical modeling: a new approach. \emph{Statistical Modelling}, v. 14, p. 21-48, 2014.
#'
#'Santos-Neto, M., Cysneiros, F.J.A., Leiva, V., Barros, M. (2016) Reparameterized Birnbaum-Saunders
#'regression models with varying precision. \emph{Electronic Journal of Statistics}, 10, 2825--2855. doi: \email{10.1214/16-EJS1187}.
#'
#'@importFrom graphics par points polygon
#'@importFrom stats qqnorm rnorm
#'@importFrom gamlss gamlss gamlss.control
#'@import ggplot2
#'@import dplyr
#'@export

envelope.RWL <- function(model, k = 100, color = "grey50", xlabel = "Theorical Quantile", ylabel = "Empirical Quantile", font = "serif")
{

  n   <- model$N
  td  <-  model$residuals
  re  <- matrix(0,n,k)

  for (i in 1:k)
  {
    y1 <- rnorm(n)
    re[, i] <- sort(y1)
  }

  e10 <- numeric(n)
  e20 <- numeric(n)
  e11 <- numeric(n)
  e21 <- numeric(n)
  e12 <- numeric(n)
  e22 <- numeric(n)

  for (l in 1:n)
  {
        eo <- sort(re[l,])
    e10[l] <- eo[ceiling(k*0.01)]
    e20[l] <- eo[ceiling(k*(1 - 0.01))]
    e11[l] <- eo[ceiling(k*0.05)]
    e21[l] <- eo[ceiling(k*(1 - 0.05))]
    e12[l] <- eo[ceiling(k*0.1)]
    e22[l] <- eo[ceiling(k*(1 - 0.1))]
  }

  a   <- qqnorm(e10, plot.it = FALSE)$x
  r   <- qqnorm(td, plot.it = FALSE)$x
  xb  <- apply(re, 1, mean)
  rxb <- qqnorm(xb, plot.it = FALSE)$x

  df  <- data.frame(r = r,xab = a,emin = cbind(e10,e11,e12),emax = cbind(e20,e21,e22),xb = xb,td = td,rxb = rxb)
  ggplot(df,aes(r,td)) + geom_ribbon(aes(x = xab, ymin = emin.e10, ymax = emax.e20),fill = color,alpha = 0.5)  + geom_ribbon(aes(x = xab, ymin = emin.e11, ymax = emax.e21),fill = color,alpha = 0.5) + geom_ribbon(aes(x = xab, ymin = emin.e12, ymax = emax.e22),fill = color,alpha = 0.5) + scale_fill_gradient(low = "grey25", high = "grey75") + geom_point() + geom_line(aes(rxb,xb),lty = 2) + xlab(xlabel) + ylab(ylabel) + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())  + theme(text = element_text(size = 10,family = font))
}



#'@name diag.RWL
#'
#'@aliases diag.RWL
#'
#'@title Diagnostic Analysis - Local Influence
#'
#'@description Diagnostics for the BP model
#'
#'@param model Object of class \code{gamlss} holding the fitted model.
#'@param scheme Default is "case.weight". But, can be "response".
#'@param mu.link  Defines the mu.link, with "identity" link as the default for the mu parameter.
#'@param sigma.link Defines the sigma.link, with "identity" link as the default for the sigma parameter.

#'
#'@return Local influence measures.
#'
#'@author
#'Manoel Santos-Neto \email{manoelferreira@uaest.ufcg.edu.br}
#'
#'@references
#'Bourguignon, M., Santos-Neto, M. and Castro, M. (2018). A new regression model for positive data. arXiv.
#'
#'@importFrom pracma hessian ones
#'@importFrom stats make.link sd
#'@export

diag.RWL <- function(model,mu.link = "log", sigma.link = "log", scheme="case.weight")
{

  x <- model$mu.x
  z <- model$sigma.x
  y <- model$y
  p <- ncol(x)
  q <- ncol(z)

  linkstr <- mu.link
  linkobj <- make.link(linkstr)
  linkfun <- linkobj$linkfun
  linkinv <- linkobj$linkinv
  mu.eta  <- linkobj$mu.eta

  sigma_linkstr <- sigma.link
  sigma_linkobj <- make.link(sigma_linkstr)
  sigma_linkfun <- sigma_linkobj$linkfun
  sigma_linkinv <- sigma_linkobj$linkinv
  sigma_mu.eta  <- sigma_linkobj$mu.eta


  B <- function(Delta,I,M)
  {
    B <- (t(Delta) %*% (I - M) %*% Delta)
    return(B)
  }

  loglik <- function(vP)
  {
    betab <- vP[1:p]
    alpha <- vP[-(1:p)]
    eta   <- as.vector(x %*% betab)
    tau   <- as.vector(z %*% alpha)
    mu    <- linkinv(eta)
    sigma <- sigma_linkinv(tau)

    a <- mu*(1 + sigma)
    b <- 2 + sigma

    sigma2 <- sigma^2
    a_us   <- sigma*(1 - mu) + sqrt(sigma2*((mu - 1)^2) + 4*mu*sigma*(sigma + 1))
    fy     <- (sigma + 1)*log(a_us) + (sigma - 1)*log(y) + log(1 + y) - (a_us/(2*mu))*y - sigma*log(2*mu) - log(a_us + 2*mu*sigma) - lgamma(sigma)

    return(sum(fy))
  }

  muest <- model$mu.coefficients
  sigmaest <- model$sigma.coefficients
  x0 <- c(muest,sigmaest)
  h0 <- hessian(loglik,x0)

  Ldelta <- h0[(p + 1):(p + q),(p + 1):(p + q)]
  Lbeta <- h0[1:p,1:p]
  b11 <- cbind(matrix(0, p, p), matrix(0, p, q))
  b12 <- cbind(matrix(0, q, p), solve(Ldelta))
  B1 <- rbind(b11, b12)  #parameter beta
  b211 <- cbind(solve(Lbeta), matrix(0, p, q))
  b212 <- cbind(matrix(0, q, p), matrix(0, q, q))
  B2 <- rbind(b211,b212)  # parameter delta

  b311 <- cbind(matrix(0, p, p), matrix(0, p, q))
  b312 <- cbind(matrix(0, q, p), matrix(0, q, q))
  B3 <- rbind(b311,b312)  # parameter theta

  if (scheme == "case.weight")
  {
    ############################Case Weight####################################

    mu <- model$mu.fv
    sigma <- model$sigma.fv
    eta <- linkfun(mu)
    ai <- mu.eta(eta)

    ll       <- function(y,mu,sigma){
      sigma2 <- sigma^2
      a_us   <- sigma*(1 - mu) + sqrt(sigma2*((mu - 1)^2) + 4*mu*sigma*(sigma + 1))
      fy     <- (sigma + 1)*log(a_us) + (sigma - 1)*log(y) + log(1 + y) - (a_us/(2*mu))*y - sigma*log(2*mu) - log(a_us + 2*mu*sigma) - lgamma(sigma)

      fy
    }

    dldm       <- Deriv::Deriv(ll,'mu')
    dldd       <- Deriv::Deriv(ll,'sigma')
    Deltamu    <- crossprod(x,diag(ai*dldm(y,mu,sigma) ))
    Deltasigma <- crossprod(z,diag(bi*dldd(y,mu,sigma)))
    Delta      <- rbind(Deltamu,Deltasigma)

    ##################theta#########################
    BT              <- B(Delta,solve(h0),B3)
    autovmaxthetaPC <- eigen(BT,symmetric = TRUE)$val[1]
    vetorpcthetaPC  <- eigen(BT,symmetric = TRUE)$vec[,1]
    dmaxG.theta     <- abs(vetorpcthetaPC)
    vCithetaPC      <- 2*abs(diag(BT))
    Cb0             <- vCithetaPC
    Cb.theta        <- Cb0/sum(Cb0)
    ######################betas########################
    BM              <- B(Delta,solve(h0),B1)
    autovmaxbetaPC  <- eigen(BM,symmetric = TRUE)$val[1]
    vetorpcbetaPC   <- eigen(BM,symmetric = TRUE)$vec[,1]
    dmaxG.beta      <- abs(vetorpcbetaPC)
    vCibetaPC       <- 2*abs(diag(BM))
    Cb1             <- vCibetaPC
    Cb.beta         <- Cb1/sum(Cb1)
    ####################alphas#########################
    BD              <- B(Delta,solve(h0),B2)
    autovmaxdeltaPC <- eigen(BD,symmetric = TRUE)$val[1]
    vetordeltaPC    <- eigen(BD,symmetric = TRUE)$vec[,1]
    dmaxG.alpha     <- abs(vetordeltaPC)
    vCideltaPC      <- 2*abs(diag(BD))
    Cb2             <- vCideltaPC
    Cb.alpha        <- Cb2/sum(Cb2)

    result <- list(dmax.beta = dmaxG.beta,
                   dmax.alpha = dmaxG.alpha,
                   dmax.theta = dmaxG.theta,
                   Ci.beta = Cb.beta,
                   Ci.alpha = Cb.alpha,
                   Ci.theta = Cb.theta)
    return(result)
  }

  if (scheme == "response")
  {
    ############################Response####################################
    mu    <- model$mu.fv
    sigma <- model$sigma.fv
    eta   <- linkfun(mu)

    ai    <- mu.eta(eta)
    tau   <- sigma_linkfun(sigma)
    bi    <- sigma_mu.eta(tau)
    sigma2 <- sigma^2
    a_us   <- sigma*(1 - mu) + sqrt(sigma2*((mu - 1)^2) + 4*mu*sigma*(sigma + 1))
    sy    <- (((2*mu)/a_us)^2)*(((sigma + 1)*((a_us + 2*sigma)^2) - a_us^2)/((a_us + 2*mu*sigma)^2))


    ll       <- function(y,mu,sigma){
      sigma2 <- sigma^2
      a_us   <- sigma*(1 - mu) + sqrt(sigma2*((mu - 1)^2) + 4*mu*sigma*(sigma + 1))
      fy     <- (sigma + 1)*log(a_us) + (sigma - 1)*log(y) + log(1 + y) - (a_us/(2*mu))*y - sigma*log(2*mu) - log(a_us + 2*mu*sigma) - lgamma(sigma)

      fy
    }

    dymu     <- Deriv::Deriv(Deriv(ll,'mu'),'y')
    dysigma  <- Deriv::Deriv(Deriv(ll,'sigma'),'y')
    Deltamu  <- crossprod(x,diag(ai*dymu(y,mu,sigma)*sy))
    p        <- ncol(x)
    q        <- ncol(z)
    Deltasigma <- crossprod(z,diag(bi*dysigma(y,mu,sigma)*sy))
    Delta <- rbind(Deltamu,Deltasigma)

    ###############thetas###########################
    BT              <- B(Delta,solve(h0),B3)
    autovmaxthetaPC <- eigen(BT,symmetric = TRUE)$val[1]
    vetorthetaRP    <- eigen(BT,symmetric = TRUE)$vec[,1]
    dmaxG.theta     <- abs(vetorthetaRP)
    vCithetaRP      <- 2*abs(diag(BT))
    Cb0             <- vCithetaRP
    Cb.theta        <- Cb0/sum(Cb0)

    #################betas##########################
    BM              <- B(Delta,solve(h0),B1)
    autovmaxbetaRP  <- eigen(BM,symmetric = TRUE)$val[1]
    vetorbetaRP     <- eigen(BM,symmetric = TRUE)$vec[,1]
    dmaxG.beta      <- abs(vetorbetaRP)
    vCibetaRP       <- 2*abs(diag(BM))
    Cb1             <- vCibetaRP
    Cb.beta         <- Cb1/sum(Cb1)
    ####################alpha#######################
    BD              <- B(Delta,solve(h0),B2)
    autovmaxdeltaRP <- eigen(BD,symmetric = TRUE)$val[1]
    vetordeltaRP    <- eigen(BD,symmetric = TRUE)$vec[,1]
    dmaxG.alpha     <- abs(vetordeltaRP)
    vCideltaRP      <- 2*abs(diag(BD))
    Cb2             <- vCideltaRP
    Cb.alpha        <- Cb2/sum(Cb2)


    result <- list(dmax.beta = dmaxG.beta,
                   dmax.alpha = dmaxG.alpha,
                   dmax.theta = dmaxG.theta,
                   Ci.beta = Cb.beta,
                   Ci.alpha = Cb.alpha,
                   Ci.theta = Cb.theta)
    return(result)
  }


}
