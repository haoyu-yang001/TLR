#' Estimate empirical-null location and scale from Z-values
#'
#' Estimates the mean and standard deviation of an empirical null
#' distribution from a vector of Z-statistics using an empirical
#' characteristic-function approach.
#'
#' @param x A numeric vector of Z-values.
#' @param gamma A positive tuning parameter controlling the frequency
#'   level used in empirical-null estimation. The default is
#'   \code{0.1}.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{\code{mu}}{The estimated empirical-null mean.}
#'   \item{\code{s}}{The estimated empirical-null standard deviation.}
#' }
#'
#' @details
#' This function provides an alternative implementation of empirical-null
#' estimation based on the empirical characteristic function. In this
#' package, it is mainly used by [JC_prop_estmate()] when estimating
#' component proportions from paired Z-statistics.
#'
#' It is not part of the main TLR procedure proposed in the paper, but is
#' useful for comparison procedures and empirical-null based auxiliary
#' calculations.
#'
#' @references
#' Jin J, Cai TT (2007). \emph{Estimating the null and the proportion of
#' nonnull effects in large-scale multiple comparisons}. Journal of the
#' American Statistical Association, 102(478), 495--506.
#'
#' @export

EstNull.func<-function (x,gamma=0.1){
  # x is a vector of z-values
  # gamma is a parameter, default is 0.1
  # output the estimated mean and standard deviation

  n = length(x)
  t = c(1:1000)/200

  gan    = n^(-gamma)
  that   = 0
  shat   = 0
  uhat   = 0
  epshat = 0

  phiplus   = rep(1,1000)
  phiminus  = rep(1,1000)
  dphiplus  = rep(1,1000)
  dphiminus = rep(1,1000)
  phi       = rep(1,1000)
  dphi      = rep(1,1000)

  for (i in 1:1000) {
    s = t[i]
    phiplus[i]   = mean(cos(s*x))
    phiminus[i]  = mean(sin(s*x))
    dphiplus[i]  = -mean(x*sin(s*x))
    dphiminus[i] = mean(x*cos(s*x))
    phi[i]       = sqrt(phiplus[i]^2 + phiminus[i]^2)
  }

  ind = min(c(1:1000)[(phi - gan) <= 0])
  tt = t[ind]
  a  = phiplus[ind]
  b  = phiminus[ind]
  da = dphiplus[ind]
  db = dphiminus[ind]
  c  = phi[ind]

  that   = tt
  shat   = -(a*da + b*db)/(tt*c*c)
  shat   = sqrt(shat)
  uhat   = -(da*b - db*a)/(c*c)
  epshat = 1 - c*exp((tt*shat)^2/2)

  return(musigma=list(mu=uhat,s=shat))
}

#' Estimate the non-null proportion from Z-values
#'
#' Estimates the proportion of non-null effects in a large collection of
#' Z-statistics using an empirical characteristic-function method.
#'
#' @param x A numeric vector of Z-values or other approximately
#'   standardized test statistics.
#' @param u A numeric value specifying the null mean.
#' @param sigma A numeric value specifying the null standard deviation.
#'
#' @return
#' A numeric scalar giving the estimated proportion of non-null signals.
#'
#' @details
#' This function is an alternative implementation of the non-null
#' proportion estimator used in Jin--Cai style empirical-null methods. In
#' this package, it is mainly used as a helper for
#' [JC_prop_estmate()] when estimating mixture proportions from paired
#' Z-statistics.
#'
#' @references
#' Jin J, Cai TT (2007). \emph{Estimating the null and the proportion of
#' nonnull effects in large-scale multiple comparisons}. Journal of the
#' American Statistical Association, 102(478), 495--506.
#'
#' @export

epsest.func <- function(x,u,sigma){
  # x is a vector
  # u is the mean
  # sigma is the standard deviation

  z  = (x - u)/sigma
  xi = c(0:100)/100
  tmax=sqrt(log(length(x)))
  tt=seq(0,tmax,0.1)

  epsest=NULL

  for (j in 1:length(tt)) {

    t=tt[j]
    f  = t*xi
    f  = exp(f^2/2)
    w  = (1 - abs(xi))
    co  = 0*xi

    for (i in 1:101) {
      co[i] = mean(cos(t*xi[i]*z));
    }
    epshat = 1 - sum(w*f*co)/sum(w)
    epsest=c(epsest,epshat)
  }
  return(epsest=max(epsest))
}

#' Estimate joint mixture proportions from paired Z-statistics
#'
#' Estimates the proportions of the four joint components induced by two
#' marginal testing problems using Jin--Cai style empirical-null and
#' non-null proportion estimation.
#'
#' @param zvalues A numeric matrix with two columns of Z-statistics. The
#'   first column corresponds to the first marginal test and the second
#'   column corresponds to the second marginal test.
#'
#' @return
#' A numeric vector of length four giving the estimated proportions of the
#' joint components. The entries correspond to:
#' \describe{
#'   \item{1}{non-null in the first margin and null in the second margin;}
#'   \item{2}{null in the first margin and non-null in the second margin;}
#'   \item{3}{null in both margins;}
#'   \item{4}{non-null in both margins.}
#' }
#'
#' @details
#' This function first estimates an empirical null distribution for each
#' margin using [EstNull.func()], then estimates the marginal non-null
#' proportions using [epsest.func()], and finally combines the marginal
#' estimates into four joint mixture proportions under a working
#' independence construction.
#'
#' In the context of large-scale mediation testing, these component
#' proportions can be useful for auxiliary estimation and comparison
#' procedures, although they are not part of the main TLR ranking and
#' thresholding steps proposed in the paper.
#'
#' @references
#' Jin J, Cai TT (2007). \emph{Estimating the null and the proportion of
#' nonnull effects in large-scale multiple comparisons}. Journal of the
#' American Statistical Association, 102(478), 495--506.
#'
#'
#' @export

JC_prop_estmate <- function(zvalues){
  musigma1=EstNull.func(zvalues[,1])
  musigma2=EstNull.func(zvalues[,2])
  pi1_first <- epsest.func(zvalues[,1],musigma1$mu,musigma1$s)
  pi2_first <- epsest.func(zvalues[,2],musigma2$mu,musigma2$s)
  return(c(pi1_first*(1-pi2_first),pi2_first*(1-pi1_first),(1-pi1_first)*(1-pi2_first),pi1_first*pi2_first))
}

# Function to reset temporary directory for each replication
