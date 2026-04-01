#' Empirical correction under one-sided partial-null configurations
#'
#' Constructs empirical CDF-based corrections for paired p-values under
#' the one-sided partial-null configurations arising in large-scale
#' mediation testing.
#'
#' @param input_pvalues A two-column numeric matrix of paired p-values.
#'   The first column contains p-values for testing
#'   \eqn{H_{0\cdot}: \gamma_j = 0}, and the second column contains
#'   p-values for testing \eqn{H_{\cdot 0}: \beta_j = 0}.
#' @param alpha1 A numeric value giving the estimated marginal null
#'   proportion for \eqn{\gamma_j = 0}.
#' @param alpha2 A numeric value giving the estimated marginal null
#'   proportion for \eqn{\beta_j = 0}.
#' @param verbose Logical; if \code{TRUE}, prints progress information
#'   during the piecewise interpolation step.
#'
#' @return
#' A two-column numeric matrix of corrected p-values. The first column
#' gives corrected values under the \eqn{H_{\gamma 0}} configuration, and
#' the second column gives corrected values under the \eqn{H_{0 \beta}}
#' configuration.
#'
#' @details
#' This function estimates empirical CDF corrections for the partial-null
#' components of the composite null hypothesis by using shape-constrained
#' estimators of the marginal p-value distributions. In this package, it
#' is used as an internal helper for the empirical version of the HDMT
#' procedure.
#'
#' @keywords internal

p_value_underH1 = function(input_pvalues, alpha1, alpha2, verbose = FALSE) {

  # input_pvalues contains two columns.
  # The first column is p_1j (p for alpha=0). The second column is p_2j (p for beta=0).
  # alpha1: proportion of null alpha=0
  # alpha2: proportion of null beta=0

  pmax <- apply(input_pvalues,1,max)
  nmed <- length(pmax)
  efdr1 <- rep(0,nmed)

  nmed  <- nrow(input_pvalues)
  cdf12 <- input_pvalues

  xx1 <- c(0,input_pvalues[order(input_pvalues[,1]),1])
  yy1 <- c(0,seq(1,nmed,by=1)/nmed)
  unique_xx1 <- unique(xx1)
  unique_yy1 <- yy1[which(!duplicated(xx1))]
  gfit1<- fdrtool::gcmlcm(unique_xx1,unique_yy1,type="lcm")
  xknots1 <- gfit1$x.knots[-1]
  Fknots1 <- cumsum(diff(gfit1$x.knots)*gfit1$slope.knots)

  xx2 <- c(0,input_pvalues[order(input_pvalues[,2]),2])
  yy2 <- c(0,seq(1,nmed,by=1)/nmed)
  unique_xx2 <- unique(xx2)
  unique_yy2 <- yy2[which(!duplicated(xx2))]
  gfit2<- fdrtool::gcmlcm(unique_xx2,unique_yy2,type="lcm")
  xknots2 <- gfit2$x.knots[-1]
  Fknots2 <- cumsum(diff(gfit2$x.knots)*gfit2$slope.knots)

  if (alpha1!=1) Fknots1 <- (Fknots1 - alpha1*xknots1)/(1-alpha1) else Fknots1 <- rep(0,length(xknots1))
  if (alpha2!=1) Fknots2 <- (Fknots2 - alpha2*xknots2)/(1-alpha2) else Fknots2 <- rep(0,length(xknots2))


  orderq1 <- pmax
  orderq2 <- pmax

  gcdf1 <- pmax
  gcdf2 <- pmax
  for (i in 1:length(xknots1)) {
    if (i==1) {
      gcdf1[orderq1<=xknots1[i]] <- (Fknots1[i]/xknots1[i])*orderq1[orderq1<=xknots1[i]]
    } else {
      if (sum(orderq1>xknots1[i-1] & orderq1<=xknots1[i])>0){
        if(verbose) {print(i)}
        temp <- orderq1[orderq1>xknots1[i-1] & orderq1<=xknots1[i]]
        gcdf1[orderq1>xknots1[i-1] & orderq1<=xknots1[i]] <- Fknots1[i-1] + (Fknots1[i]-Fknots1[i-1])/(xknots1[i]-xknots1[i-1])*(temp-xknots1[i-1])
      }
    }
  }

  for (i in 1:length(xknots2)) {
    if (i==1) {
      gcdf2[orderq2<=xknots2[i]] <- (Fknots2[i]/xknots2[i])*orderq2[orderq2<=xknots2[i]]
    } else {
      if (sum(orderq2>xknots2[i-1] & orderq2<=xknots2[i])>0){
        temp <- orderq2[orderq2>xknots2[i-1] & orderq2<=xknots2[i]]
        gcdf2[orderq2>xknots2[i-1] & orderq2<=xknots2[i]] <- Fknots2[i-1] + (Fknots2[i]-Fknots2[i-1])/(xknots2[i]-xknots2[i-1])*(temp-xknots2[i-1])
      }
    }
  }


  gcdf1 <- ifelse(gcdf1>1,1,gcdf1)
  gcdf2 <- ifelse(gcdf2>1,1,gcdf2)

  cdf12[,1] <- pmax(gcdf1,2) # p_1j under H10
  cdf12[,2] <- pmax(gcdf2,2) # p_2j under H01

  return(cdf12 = cdf12)

}

#' Empirical HDMT p-values for composite-null mediation testing
#'
#' Computes p-values from the empirical version of the
#' high-dimensional mediation testing (HDMT) procedure for testing the
#' composite null hypothesis of no mediation effect.
#'
#' @param p.alpha A numeric vector of p-values for testing
#'   \eqn{H_{0\cdot}: \gamma_j = 0}.
#' @param p.beta A numeric vector of p-values for testing
#'   \eqn{H_{\cdot 0}: \beta_j = 0}.
#' @param estws A list of estimated null proportions. It should contain
#'   \code{alpha00}, \code{alpha10}, \code{alpha01}, \code{alpha1}, and
#'   \code{alpha2}.
#'
#' @return
#' A numeric vector of empirical HDMT p-values.
#'
#' @details
#' HDMT is an existing method for large-scale mediation testing that
#' models the composite-null distribution of the MaxP statistic rather
#' than assuming a uniform null distribution. This function implements the
#' empirically corrected version, which uses estimated partial-null
#' distributions to improve calibration relative to the purely asymptotic
#' approximation.
#'
#' In the paper, empirical HDMT is treated as a comparison method to the
#' proposed TLR procedure.
#'
#' @references
#' Dai JY, Stanford JL, LeBlanc M (2022). \emph{A multiple-testing
#' procedure for high-dimensional mediation hypotheses}. Journal of the
#' American Statistical Association.
#'
#'
#' @export
HDMT_emp = function(p.alpha,p.beta,estws){
  input_pvalues = cbind(p.alpha,p.beta)
  pmax <- apply(input_pvalues,1,max)
  c = estws$alpha10+estws$alpha01+estws$alpha00
  pi.01.est = estws$alpha01/c
  pi.10.est = estws$alpha10/c
  pi.00.est = estws$alpha00/c

  x <- tryCatch(p_value_underH1(input_pvalues = input_pvalues, alpha1 = estws$alpha1, alpha2 = estws$alpha2), error = function(e) e)
  if(!inherits(x, "error")){
    correction = p_value_underH1(input_pvalues = input_pvalues, alpha1 = estws$alpha1, alpha2 = estws$alpha2)
    p_1j_H10 = correction[,1]
    p_2j_H01 = correction[,2]
    p.HDMT = pi.10.est*pmax*p_1j_H10+pi.01.est*pmax*p_2j_H01+pi.00.est*pmax^2
  }else if(inherits(x, "error")){
    p.HDMT = pi.10.est*pmax+pi.01.est*pmax+pi.00.est*pmax^2
  }

  p.HDMT[which(p.HDMT<=0)] = min(p.HDMT[which(p.HDMT>0)])
  p.HDMT[which(p.HDMT>=1)] = 1-1e-8

  return(p.HDMT)
}

#' Asymptotic HDMT p-values for composite-null mediation testing
#'
#' Computes p-values from the asymptotic version of the
#' high-dimensional mediation testing (HDMT) procedure for testing the
#' composite null hypothesis of no mediation effect.
#'
#' @param p.alpha A numeric vector of p-values for testing
#'   \eqn{H_{0\cdot}: \gamma_j = 0}.
#' @param p.beta A numeric vector of p-values for testing
#'   \eqn{H_{\cdot 0}: \beta_j = 0}.
#' @param estws A list of estimated null proportions. It should contain
#'   \code{alpha00}, \code{alpha10}, and \code{alpha01}.
#'
#' @return
#' A numeric vector of asymptotic HDMT p-values.
#'
#' @details
#' This function implements the asymptotic HDMT approximation based on the
#' estimated proportions of the three null components under the composite
#' null. Compared with [HDMT_emp()], it does not use empirical correction
#' for the partial-null cases.
#'
#' In the paper, the asymptotic HDMT procedure is included as an existing
#' comparison method and is noted to be more conservative than its
#' empirical counterpart.
#'
#' @references
#' Dai JY, Stanford JL, LeBlanc M (2022). \emph{A multiple-testing
#' procedure for high-dimensional mediation hypotheses}. Journal of the
#' American Statistical Association.
#'
#'
#' @export
HDMT_asy = function(p.alpha,p.beta,estws){
  input_pvalues = cbind(p.alpha,p.beta)
  pmax <- apply(input_pvalues,1,max)
  nullprop <- estws
  c = nullprop$alpha10+nullprop$alpha01+nullprop$alpha00
  pi.01.est = nullprop$alpha01/c
  pi.10.est = nullprop$alpha10/c
  pi.00.est = nullprop$alpha00/c

  p.HDMT = pi.10.est*pmax+pi.01.est*pmax+pi.00.est*pmax^2

  p.HDMT[which(p.HDMT<=0)] = min(p.HDMT[which(p.HDMT>0)])
  p.HDMT[which(p.HDMT>=1)] = 1-1e-8

  return(p.HDMT)
}

#' Multiple-testing rule based on HDMT
#'
#' Applies the HDMT procedure to paired mediation p-values and returns the
#' selected hypotheses under either size control or false discovery rate
#' (FDR) control.
#'
#' @param p.M A numeric vector of p-values for testing
#'   \eqn{H_{0\cdot}: \gamma_j = 0}.
#' @param p.Y A numeric vector of p-values for testing
#'   \eqn{H_{\cdot 0}: \beta_j = 0}.
#' @param estws A list of estimated null proportions.
#' @param coefs Included for interface compatibility with other methods.
#'   It is not used directly in the HDMT calculation.
#' @param exact Numeric or logical indicator specifying which HDMT version
#'   to use: \code{1} for the empirical version and \code{0} for the
#'   asymptotic version.
#' @param significance_upper A numeric value specifying the target size or
#'   FDR level.
#' @param control.method A character string specifying the error-control
#'   target. Supported values are \code{"size"} and \code{"FDR"}.
#'
#' @return
#' An integer vector containing the indices of selected hypotheses.
#'
#' @details
#' This function computes either empirical HDMT p-values or asymptotic
#' HDMT p-values and then applies a multiple-testing rule based on the
#' requested error criterion. It is primarily included for comparison with
#' the proposed TLR procedure in simulation studies and empirical data
#' analyses.
#'
#' @references
#' Dai JY, Stanford JL, LeBlanc M (2022). \emph{A multiple-testing
#' procedure for high-dimensional mediation hypotheses}. Journal of the
#' American Statistical Association.
#'
#'
#' @export

HDMT_control = function(p.M,p.Y,estws,coefs, exact = 1,significance_upper, control.method){
  sim.num <- length(p.M)

  F_emp <- function(x){x}

  if(exact==1){
    Ts <- HDMT_emp(p.M,p.Y,estws)
  }
  if(exact==0){
    Ts <- HDMT_asy(p.M,p.Y,estws)
  }

  ss <- search_sorted_index(F_emp,Ts,control.method,significance_upper)

  return(ss)
}

#' Binary search for a multiple-testing cutoff
#'
#' Uses a binary search over sorted statistics to identify the rejection
#' cutoff under either size control or false discovery rate (FDR) control.
#'
#' @param F_emp A function giving the null cumulative distribution
#'   function evaluated at a statistic value.
#' @param Ts A numeric vector of statistics or p-values to be thresholded.
#' @param control.method A character string specifying the error-control
#'   target. Supported values are \code{"size"} and \code{"FDR"}.
#' @param significance_upper A numeric value specifying the target level.
#'
#' @return
#' An integer vector containing the indices of selected hypotheses.
#'
#' @details
#' This internal helper sorts the input statistics and searches for the
#' largest rejection threshold satisfying the requested criterion. Under
#' size control, it uses the null CDF directly; under FDR control, it uses
#' the ratio of the estimated null fraction to the empirical rejection
#' fraction.
#'
#' @keywords internal

search_sorted_index <- function(F_emp, Ts, control.method,significance_upper) {

  sorted_index <- order(Ts) # Get the indices that sort 'Ts'
  lower_bound <- 1
  upper_bound <- length(Ts)

  if(control.method == "size"){
    while (upper_bound - lower_bound > 1) {
      mid_index <- floor((lower_bound + upper_bound) / 2)
      mid_value <- Ts[sorted_index[mid_index]]

      if (F_emp(mid_value) < significance_upper) {
        lower_bound <- mid_index
      } else {
        upper_bound <- mid_index
      }
    }

    t_s <- Ts[sorted_index[lower_bound]]

    if(F_emp(t_s) - significance_upper < 0){
      return(which(Ts <= t_s))
    }else{
      return(which(Ts < t_s))
    }
  }

  if(control.method == "FDR"){
    while (upper_bound - lower_bound > 1) {
      mid_index <- floor((lower_bound + upper_bound) / 2)
      mid_value <- Ts[sorted_index[mid_index]]

      if (F_emp(mid_value)*length(Ts)/mid_index < significance_upper ) {
        lower_bound <- mid_index
      } else {
        upper_bound <- mid_index
      }
    }

    t_s <- Ts[sorted_index[lower_bound]]

    if(F_emp(t_s)*length(Ts)/mid_index - significance_upper < 0){
      return(which(Ts <= t_s))
    }else{
      return(which(Ts < t_s))
    }

  }
}
