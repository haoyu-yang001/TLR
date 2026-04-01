#' Jin--Cai empirical-null correction for p-values
#'
#' Applies a Jin--Cai style empirical-null correction to a vector of
#' p-values by first transforming them to Z-scores, estimating the
#' empirical null location and scale, and then recalculating tail
#' probabilities under the estimated null distribution.
#'
#' @param pval A numeric vector of p-values.
#'
#' @return
#' A numeric vector of empirically corrected p-values.
#'
#' @details
#' In large-scale testing problems, the theoretical null distribution may
#' be distorted by dependence, model misspecification, or other sources of
#' deviation from ideal assumptions. This function uses the empirical-null
#' estimation procedure implemented in [nullParaEst()] to adjust p-values
#' accordingly.
#'
#' This correction is not part of the main TLR procedure proposed in the
#' paper, but can be used as an auxiliary correction step for comparison
#' methods such as DACT-style procedures.
#'
#' @references
#' Jin J, Cai TT (2007). \emph{Estimating the null and the proportion of
#' nonnull effects in large-scale multiple comparisons}. Journal of the
#' American Statistical Association, 102(478), 495--506.
#'
#'
#' @export
JCCorrect = function(pval){
  z = stats::qnorm(pval,lower.tail = F)
  res= nullParaEst(z)
  pval.JC = stats::pnorm(z,mean = res$mu,sd = res$s,lower.tail = F)
  return(pval.JC)
}


#' Efron empirical-null correction for p-values
#'
#' Applies Efron's empirical-null adjustment to a vector of p-values using
#' the \pkg{locfdr} framework. The function estimates the empirical-null
#' mean and standard deviation from transformed Z-scores and recomputes
#' the p-values under the estimated null distribution.
#'
#' @param pval A numeric vector of p-values.
#'
#' @return
#' A numeric vector of empirically corrected p-values.
#'
#' @details
#' This function is useful when the theoretical standard normal null may
#' not provide an adequate approximation, for example in large-scale
#' testing settings with dependence or inflation. It is mainly included as
#' an auxiliary correction tool for comparison procedures rather than as a
#' component of the main TLR method.
#'
#' @references
#' Efron B (2004). \emph{Large-scale simultaneous hypothesis testing: the
#' choice of a null hypothesis}. Journal of the American Statistical
#' Association, 99(465), 96--104.
#'
#'
#' @export
EfronCorrect = function(pval){
  z = stats::qnorm(1-pval)
  res <- locfdr::locfdr(z,nulltype = 1)
  mean.emp = res$fp0["mlest","delta"]
  sd.emp = res$fp0["mlest","sigma"]
  pval.emp = stats::pnorm(z,mean = mean.emp,sd = sd.emp,lower.tail = F)
  return(pval.emp)
}

#' Estimate the proportion of non-null signals
#'
#' Estimates the proportion of non-null effects in a large collection of
#' test statistics using a Fourier-based method. This function can be used
#' to estimate marginal null proportions when constructing comparison
#' procedures for large-scale mediation testing.
#'
#' @param x A numeric vector of test statistics, typically Z-scores.
#' @param u A numeric value specifying the null mean.
#' @param sigma A numeric value specifying the null standard deviation.
#'
#' @return
#' A numeric scalar giving the estimated proportion of non-null signals.
#'
#' @details
#' The estimator is based on the empirical characteristic function and is
#' used here to estimate the marginal non-null proportion from a large set
#' of test statistics. In this package, it is primarily used as a
#' computational component in DACT-style combination rules and related
#' empirical-null adjustments.
#'
#' @references
#' Jin J, Cai TT (2007). \emph{Estimating the null and the proportion of
#' nonnull effects in large-scale multiple comparisons}. Journal of the
#' American Statistical Association, 102(478), 495--506.
#'
#' @export
nonnullPropEst <- function(x,u,sigma){
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

#' Estimate empirical-null mean and scale from Z-scores
#'
#' Estimates the location and scale parameters of an empirical null
#' distribution from a large collection of Z-scores using the empirical
#' characteristic function.
#'
#' @param x A numeric vector of Z-scores.
#' @param gamma A tuning parameter controlling the frequency level used in
#'   the empirical characteristic-function estimator. The default is
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
#' This function implements an empirical-null estimation step that can be
#' used when the theoretical standard normal null is not adequate. In this
#' package, it mainly serves as a helper for empirical-null correction
#' procedures such as [JCCorrect()].
#'
#' @references
#' Jin J, Cai TT (2007). \emph{Estimating the null and the proportion of
#' nonnull effects in large-scale multiple comparisons}. Journal of the
#' American Statistical Association, 102(478), 495--506.
#'
#' @export
nullParaEst<-function (x,gamma=0.1){
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


#' DACT-style p-value combination for large-scale mediation testing
#'
#' Combines two marginal p-value vectors into a single p-value for testing
#' the composite null hypothesis of no mediation effect using a
#' divide-aggregate composite-null testing (DACT) style construction.
#'
#' @param p_a A numeric vector of p-values for testing the first marginal
#'   null component, typically \eqn{H_{0\cdot}: \gamma_j = 0}.
#' @param p_b A numeric vector of p-values for testing the second marginal
#'   null component, typically \eqn{H_{\cdot 0}: \beta_j = 0}.
#' @param Z_a An optional numeric vector of Z-scores corresponding to
#'   \code{p_a}. If set to \code{"NULL"}, the Z-scores are computed
#'   internally from \code{p_a}.
#' @param Z_b An optional numeric vector of Z-scores corresponding to
#'   \code{p_b}. If set to \code{"NULL"}, the Z-scores are computed
#'   internally from \code{p_b}.
#' @param correction A character string specifying whether to apply an
#'   empirical-null correction to the combined p-values. Supported values
#'   are \code{"NULL"}, \code{"Efron"}, and \code{"JC"}.
#'
#' @return
#' A numeric vector of DACT-style combined p-values. If the requested
#' empirical-null correction fails, the function returns \code{NULL}.
#'
#' @details
#' The DACT framework combines information across the three composite-null
#' cases by estimating marginal null proportions and then forming a
#' weighted combination of \code{p_a}, \code{p_b}, and
#' \eqn{\max(p_a, p_b)^2}. This function provides a convenient
#' implementation for method comparison in large-scale mediation testing.
#'
#' In the paper, DACT is treated as an existing method for comparison with
#' the proposed TLR procedure. The paper notes that DACT can perform well
#' under weak and sparse signals but may exhibit substantial error-rate
#' inflation in dense and strong signal settings.
#'
#' @references
#' Liu Z, Shen J, Barfield R, Schwartz J, Baccarelli AA, Lin X (2022).
#' \emph{Large-scale hypothesis testing for causal mediation effects with
#' applications in genome-wide epigenetic studies}. Journal of the
#' American Statistical Association.
#'
#'
#' @export
DACT_liu <- function (p_a, p_b,Z_a= "NULL", Z_b= "NULL", correction = "NULL") {

  if(all(Z_a== "NULL" & Z_b == "NULL")){
    Z_a = stats::qnorm(p_a, lower.tail = F)
    Z_b = stats::qnorm(p_b, lower.tail = F)
  }

  pi0a = 1 - nonnullPropEst(Z_a, 0, 1)
  pi0b = 1 - nonnullPropEst(Z_b, 0, 1)
  if (pi0a > 1) {
    pi0a = 1
  }
  if (pi0b > 1) {
    pi0b = 1
  }
  p.mat = cbind(p_a, p_b)
  p3 = (apply(p.mat, 1, max))^2
  wg1 = pi0a * (1 - pi0b)
  wg2 = (1 - pi0a) * pi0b
  wg3 = pi0a * pi0b
  wg.sum = wg1 + wg2 + wg3
  wg.std = c(wg1, wg2, wg3)/wg.sum
  p_dact = wg.std[1] * p_a + wg.std[2] * p_b + wg.std[3] * p3
  if (correction == "Efron") {
    x <- tryCatch(EfronCorrect(p_dact), error = function(e) e)
    if(inherits(x, "error")){
      p_dact = NULL
    }else{
      p_dact = EfronCorrect(p_dact)
    }
  }
  if (correction == "JC") {
    x <- tryCatch(JCCorrect(p_dact), error = function(e) e)
    if(inherits(x, "error")){
      p_dact = NULL
    }else{
      p_dact = JCCorrect(p_dact)
    }
  }
  return(p_dact)
}
