% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RCodes.R
\name{RWL}
\alias{RWL}
\alias{dRWL}
\alias{pRWL}
\alias{qRWL}
\alias{rRWL}
\alias{hRWL}
\alias{plotRWL}
\alias{meanRWL}
\title{Reparameterized Weighted Lindley (RWL) distribution for fitting a GAMLSS}
\usage{
RWL(mu.link = "log", sigma.link = "log")

dRWL(x, mu = 1, sigma = 1, log.d = FALSE)

pRWL(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)

qRWL(p, mu = 1, sigma = 1, lower.tail = TRUE, log.q = FALSE)

rRWL(n, mu = 1, sigma = 1)

hRWL(x, mu = 1, sigma = 1, log.h = FALSE)

plotRWL(
  mu = 0.5,
  sigma = 1,
  from = 0,
  to = 0.999,
  n = 101,
  title = "title",
  ...
)

meanRWL(obj)
}
\arguments{
\item{mu.link}{object for which the extraction of model residuals is meaningful.}

\item{sigma.link}{type of residual to be used.}

\item{x, q}{vector of quantiles}

\item{mu}{vector of scale parameter values}

\item{sigma}{vector of shape parameter values}

\item{log.d, }{log.p logical; if TRUE, density d are given as log(d).}

\item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]}

\item{log.p, }{log.p logical; if TRUE, probabilities p are given as log(p).}

\item{p}{vector of probabilities.}

\item{log.q, }{log.p logical; if TRUE, quantiles q are given as log(q).}

\item{n}{number of observations. If \code{length(n) > 1}, the length is taken to be the number required.}

\item{log.h, }{log.p logical; if TRUE, hazard h are given as log(h).}

\item{from}{where to start plotting the distribution from.}

\item{to}{up to where to plot the distribution.}

\item{title}{title of the plot.}

\item{...}{other graphical parameters for plotting.}

\item{obj}{a fitted RWL object.}
}
\value{
returns a \code{gamlss.family} object which can be used to fit a reparameterized weighted Lindley distribution in the \code{gamlss()} function.
}
\description{
The fuction \code{RWL()} defines the RWL distribution, a two paramenter
distribution, for a gamlss.family object to be used in GAMLSS fitting using using the
function \code{gamlss()}, with mean equal to the parameter \code{mu} and \code{sigma}
equal the precision parameter. The functions \code{dRWL}, \code{pWL}, \code{qRWL},
\code{rRWL} and \code{hRWL} define the density, distribution function, quantile function, random
genetation and hazard function for the \code{RWL} parameterization of the RWL distribution.
}
\note{
For the function RWL(), mu is the mean and sigma is the precision parameter of the reparameterized weighted Lindley distribution.
}
\author{
Manoel Santos-Neto \email{santosnetoce at protemail.com}
}
