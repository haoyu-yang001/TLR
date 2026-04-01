#' Thresholding step for the MDACT comparison procedure
#'
#' Computes the rejection set for an MDACT-style multiple-testing
#' procedure under the composite null hypothesis of no mediation effect.
#'
#' @param p.M A numeric vector of p-values for testing
#'   \eqn{H_{0\cdot}: \gamma_j = 0}.
#' @param p.Y A numeric vector of p-values for testing
#'   \eqn{H_{\cdot 0}: \beta_j = 0}.
#' @param pi.10.est A numeric value giving the estimated proportion of the
#'   partial-null configuration \eqn{H_{\gamma 0}}.
#' @param pi.01.est A numeric value giving the estimated proportion of the
#'   partial-null configuration \eqn{H_{0 \beta}}.
#' @param pi.00.est A numeric value giving the estimated proportion of the
#'   complete null configuration \eqn{H_{00}}.
#' @param significance_upper A numeric value specifying the target size or
#'   false discovery rate (FDR) level.
#' @param control.method A character string specifying the error-control
#'   target. Supported values are \code{"size"} and \code{"FDR"}.
#' @param ifcond A numeric indicator controlling the calibration rule used
#'   under FDR control.
#'
#' @return
#' An integer vector containing the indices of selected hypotheses.
#'
#' @details
#' This function evaluates an MDACT-style composite test statistic based
#' on weighted combinations of the two marginal p-values and
#' \eqn{\max(p_M, p_Y)^2}, then determines a rejection cutoff using an
#' estimated null distribution under the three null components of the
#' composite null hypothesis.
#'
#' This function is included as part of a comparison procedure rather than
#' as a component of the main TLR method proposed in the paper.
#'
#' @keywords internal

DACT_thr <- function(p.M, p.Y,pi.10.est,pi.01.est,pi.00.est,significance_upper,control.method,ifcond){
  sim.num <- length(p.M)
  c <- pi.01.est + pi.10.est + pi.00.est
  pi.11.est <- max(1 - c, 0)

  Ts <- pi.01.est*p.M + pi.10.est*p.Y + pi.00.est*pmax(p.M,p.Y)^2

  F00 <- function(t){
    F_0int <- function(p.y){
      c0 <- pi.00.est * p.y^2 + pi.01.est * p.y
      c1 <- t- pi.10.est * p.y
      c2 <- (t - pi.00.est * p.y^2 - pi.10.est * p.y)/pi.01.est
      m2 <- (-pi.01.est+sqrt(pmax(0,pi.01.est^2 + 4*pi.00.est*(t-pi.10.est*p.y))))/(2*pi.00.est)
      s1 <- pmin(pmax(c2,0), 1)
      s2 <- pmin(pmax(m2,0), 1)
      s <- rep(0, length(p.y))
      s[c0>=c1] <- s1[c0>=c1]
      s[c0<c1] <- s2[c0<c1]
      s
    }

    x <- tryCatch(integrate(F_0int, 0, 1, rel.tol = .Machine$double.eps), error = function(e) e)
    if(!inherits(x, "error")){
      res <- integrate(F_0int, 0, 1, rel.tol = .Machine$double.eps)
    }else{
      res <- integrate(F_0int, 0, 1)
    }
    res$value
  }

  F01 <- function(t){
    c0 <- pi.00.est * p.Y^2 + pi.01.est * p.Y
    c1 <- t- pi.10.est * p.Y
    c2 <- (t - pi.00.est * p.Y^2 - pi.10.est * p.Y)/pi.01.est
    m2 <- (-pi.01.est+sqrt(pmax(0,pi.01.est^2 + 4*pi.00.est*(t-pi.10.est*p.Y))))/(2*pi.00.est)
    s1 <- pmin(pmax(c2,0), 1)
    s2 <- pmin(pmax(m2,0), 1)

    (sum(s1[c0>=c1]) + sum(s2[c0<c1])) / length(p.Y)
  }

  F10 <- function(t){
    c0 <- pi.00.est * p.M^2 + pi.10.est * p.M
    c1 <- t- pi.01.est * p.M
    c2 <- (t - pi.00.est * p.M^2 - pi.01.est * p.M)/pi.10.est
    m2 <- (-pi.10.est+sqrt(pmax(0,pi.10.est^2 + 4*pi.00.est*(t-pi.01.est*p.M))))/(2*pi.00.est)
    s1 <- pmin(pmax(c2,0), 1)
    s2 <- pmin(pmax(m2,0), 1)

    (sum(s1[c0>=c1]) + sum(s2[c0<c1])) / length(p.M)
  }

  if(control.method == "FDR"){
    if(ifcond == 0){
      thr <- function(t) {
        c * (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*F10(t) +
                                pi.01.est/(pi.01.est+pi.11.est)*F01(t) +
                                (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) -
                                   pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*F00(t))) / max(mean(Ts < t), 1 / sim.num) - significance_upper
      }
    }
    if(ifcond == 1){
      if(length(which(p.adjust(p.M,method = "fdr") < 0.1)) < 50 | length(which(p.adjust(p.Y,method = "fdr") < 0.1)) < 50){
        est10 <- 1
        est01 <- 1
        thr <- function(t) {
          (F10(t) + F01(t) - F00(t)) / max(mean(Ts < t), 1 / sim.num) - significance_upper
        }
      }else{
        est10 <- mean(p.adjust(p.M,method = "fdr") < 0.1 & p.Y > 0.9)/stats::punif(1-0.9)/mean(p.adjust(p.M,method = "fdr") < 0.1)
        est01 <- mean(p.adjust(p.Y,method = "fdr") < 0.1 & p.M > 0.9)/stats::punif(1-0.9)/mean(p.adjust(p.Y,method = "fdr") < 0.1)
      }
      thr <- function(t) {
        (est10*F10(t) + est01*F01(t) +
           (pi.00.est - est10*(pi.00.est + pi.01.est) - est01*(pi.00.est + pi.10.est))*F00(t)) / max(mean(Ts < t), 1 / sim.num) - significance_upper
      }
    }
    #t_s <- uniroot(thr, c(min(c(p.M, p.Y)), 1))$root
    sorted_index <- order(Ts) # Get the indices that sort 'Ts'
    lower_bound <- 1
    upper_bound <- length(Ts)
    thr1 <- 50
    while (upper_bound - lower_bound > thr1) {
      mid_index <- floor((lower_bound + upper_bound) / 2)
      mid_value <- Ts[sorted_index[mid_index]]

      x <- tryCatch(thr(mid_value), error = function(e) e)
      k <- 10
      while(inherits(x, "error")){
        x <- tryCatch(thr(round(mid_value,k)), error = function(e) e)
        k <- k-1
      }
      if (x < 0) {
        lower_bound <- mid_index
      } else {
        upper_bound <- mid_index
      }
    }
    t_s <- Ts[sorted_index[lower_bound]]

    if(thr(t_s) < 0){
      return(which(Ts <= t_s))
    }else{
      return(which(Ts < t_s))
    }

  }
  if(control.method == "size"){
    thr <- function(t) {
      (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*F10(t) +
                          pi.01.est/(pi.01.est+pi.11.est)*F01(t) +
                          (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) -
                             pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*F00(t))) - significance_upper
    }
    x <- tryCatch(uniroot(thr, c(min(c(p.M, p.Y)), 1))$root, error = function(e) e)
    F_emp <- function(t) {
      F_emp <- (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*F10(t) +
                                   pi.01.est/(pi.01.est+pi.11.est)*F01(t) +
                                   (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) -
                                      pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*F00(t)))
      return(F_emp)
    }
    if(!inherits(x, "error")){
      t_s <- x
    }else{
      sorted_index <- order(Ts) # Get the indices that sort 'Ts'
      lower_bound <- 1
      upper_bound <- length(Ts)

      while (upper_bound - lower_bound > 1) {
        mid_index <- floor((lower_bound + upper_bound) / 2)
        mid_value <- Ts[sorted_index[mid_index]]

        x <- tryCatch(F_emp(round(mid_value,6)), error = function(e) e)

        if(!inherits(x, "error")){
          fn <- F_emp(round(mid_value,6))
        }else{
          fn <- 0
        }

        if (fn < significance_upper) {
          lower_bound <- mid_index
        } else {
          upper_bound <- mid_index
        }
      }
      t_s <- Ts[sorted_index[lower_bound]]
    }

    x <- tryCatch(F_emp(t_s), error = function(e) e)
    if(!inherits(x, "error")){
      if(F_emp(t_s) < significance_upper){
        return(which(Ts <= t_s))
      }else{
        return(which(Ts < t_s))
      }
    }else{
      return(which(Ts <= t_s))
    }
  }
}

#' MDACT p-values for composite-null mediation testing
#'
#' Computes p-values associated with an MDACT-style test statistic for the
#' composite null hypothesis of no mediation effect.
#'
#' @param p.M A numeric vector of p-values used to estimate the partial-null
#'   distributions associated with the first marginal test.
#' @param p.Y A numeric vector of p-values used to estimate the partial-null
#'   distributions associated with the second marginal test.
#' @param p.M_one A numeric vector of p-values entering the MDACT
#'   statistic for the first marginal test.
#' @param p.Y_one A numeric vector of p-values entering the MDACT
#'   statistic for the second marginal test.
#' @param estws A list of estimated null proportions. It should contain
#'   \code{alpha01}, \code{alpha10}, and \code{alpha00}.
#'
#' @return
#' A numeric vector of MDACT p-values.
#'
#' @details
#' The function evaluates an MDACT-style statistic of the form
#' \eqn{\pi_{0\beta} p_M + \pi_{\gamma 0} p_Y + \pi_{00} \max(p_M, p_Y)^2}
#' and computes p-values using an estimated null distribution obtained by
#' combining the three null components of the composite null hypothesis.
#'
#' In the paper, procedures of this type are considered as comparison
#' methods for large-scale mediation testing rather than the primary
#' proposed approach.
#'
#' @export

MDACT_pvalues = function(p.M,p.Y,p.M_one,p.Y_one,estws){
  sim.num <- length(p.M)

  pi.01.est = max(1e-3, estws$alpha01)
  pi.10.est = max(1e-3, estws$alpha10)
  pi.00.est = max(1e-3, estws$alpha00)
  c <- pi.01.est + pi.10.est + pi.00.est
  pi.11.est <- 1-max(c,1-1e-6)

  Ts <- pi.01.est*p.M_one + pi.10.est*p.Y_one + pi.00.est*pmax(p.M_one,p.Y_one)^2

  F00 <- function(t){
    F_0int <- function(p.y){
      c0 <- pi.00.est * p.y^2 + pi.01.est * p.y
      c1 <- t- pi.10.est * p.y
      c2 <- (t - pi.00.est * p.y^2 - pi.10.est * p.y)/pi.01.est
      m2 <- (-pi.01.est+sqrt(pmax(0,pi.01.est^2 + 4*pi.00.est*(t-pi.10.est*p.y))))/(2*pi.00.est)
      s1 <- pmin(pmax(c2,0), 1)
      s2 <- pmin(pmax(m2,0), 1)
      s <- rep(0, length(p.y))
      s[c0>=c1] <- s1[c0>=c1]
      s[c0<c1] <- s2[c0<c1]
      s
    }

    x <- tryCatch(integrate(F_0int, 0, 1, rel.tol = .Machine$double.eps), error = function(e) e)
    if(!inherits(x, "error")){
      res <- integrate(F_0int, 0, 1, rel.tol = .Machine$double.eps)
    }else{
      res <- integrate(F_0int, 0, 1)
    }
    res$value
  }

  F01 <- function(t){
    c0 <- pi.00.est * p.Y^2 + pi.01.est * p.Y
    c1 <- t- pi.10.est * p.Y
    c2 <- (t - pi.00.est * p.Y^2 - pi.10.est * p.Y)/pi.01.est
    m2 <- (-pi.01.est+sqrt(pmax(0,pi.01.est^2 + 4*pi.00.est*(t-pi.10.est*p.Y))))/(2*pi.00.est)
    s1 <- pmin(pmax(c2,0), 1)
    s2 <- pmin(pmax(m2,0), 1)

    (sum(s1[c0>=c1]) + sum(s2[c0<c1])) / length(p.Y)
  }

  F10 <- function(t){
    c0 <- pi.00.est * p.M^2 + pi.10.est * p.M
    c1 <- t- pi.01.est * p.M
    c2 <- (t - pi.00.est * p.M^2 - pi.01.est * p.M)/pi.10.est
    m2 <- (-pi.10.est+sqrt(pmax(0,pi.10.est^2 + 4*pi.00.est*(t-pi.01.est*p.M))))/(2*pi.00.est)
    s1 <- pmin(pmax(c2,0), 1)
    s2 <- pmin(pmax(m2,0), 1)

    (sum(s1[c0>=c1]) + sum(s2[c0<c1])) / length(p.M)
  }

  Fnull <- function(t) {
    (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*F10(t) +
                        pi.01.est/(pi.01.est+pi.11.est)*F01(t) +
                        (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) -
                           pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*F00(t)))
  }

  pv <- vapply(Ts, Fnull, numeric(1))

  return(pv)
}

#' Adaptive Benjamini--Hochberg procedure with estimated null proportion
#'
#' Applies an adaptive Benjamini--Hochberg (BH) procedure using a simple
#' plug-in estimator of the null proportion.
#'
#' @param p_values A numeric vector of p-values.
#' @param q A numeric value specifying the target false discovery rate
#'   (FDR) level.
#' @param delta A positive numeric value specifying the step size used for
#'   the lambda grid in null-proportion estimation.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{\code{rejection_set}}{The indices of rejected hypotheses.}
#'   \item{\code{lambda_opt}}{The selected tuning parameter value.}
#'   \item{\code{pi_hat_opt}}{The estimated null proportion at the chosen
#'   tuning parameter.}
#' }
#'
#' @details
#' This function implements a simple adaptive BH rule by first estimating
#' the null proportion over a grid of tuning parameters and then applying
#' the BH procedure with the adjusted rejection threshold. It is included
#' as a general multiple-testing helper and is not specific to the TLR
#' methodology.
#'
#' @references
#' Benjamini Y, Hochberg Y (1995). \emph{Controlling the false discovery
#' rate: a practical and powerful approach to multiple testing}. Journal
#' of the Royal Statistical Society: Series B, 57(1), 289--300.
#'
#' @export

adaptive_bh_as <- function(p_values, q, delta = 0.01) {
  n <- length(p_values)

  # Step 1: Find the optimal lambda using the stopping rule
  lambda_grid <- seq(q, 1, by = delta)
  pi_hat <- numeric(length(lambda_grid))

  for (i in seq_along(lambda_grid)) {
    lambda <- lambda_grid[i]
    pi_hat[i] <- (1 + sum(p_values >= lambda)) / (n * (1 - lambda))
  }

  ambda_opt_index <- 1
  for (i in 2:length(pi_hat)) {
    if (pi_hat[i] > pi_hat[i - 1]) {
      lambda_opt_index <- i - 1
      break
    }
  }
  # Find the smallest lambda where pi_hat stops decreasing
  lambda_opt <- lambda_grid[lambda_opt_index]
  pi_hat_opt <- pi_hat[lambda_opt_index]

  # Step 2: Apply the BH procedure with the adjusted threshold
  p_sorted <- sort(p_values)
  k <- max(which(p_sorted <= (1:n) * q / (n * pi_hat_opt)))
  threshold <- p_sorted[k]

  # Return the results
  rejection_set <- which(p_values <= threshold)
  list(rejection_set = rejection_set, lambda_opt = lambda_opt, pi_hat_opt = pi_hat_opt)
}

#' MDACT multiple-testing procedure for mediation analysis
#'
#' Applies an MDACT-style multiple-testing rule to paired mediation
#' p-values and returns the selected hypotheses.
#'
#' @param p.M A numeric vector of p-values for testing
#'   \eqn{H_{0\cdot}: \gamma_j = 0}.
#' @param p.Y A numeric vector of p-values for testing
#'   \eqn{H_{\cdot 0}: \beta_j = 0}.
#' @param estws A list of estimated null proportions. It should contain
#'   \code{alpha01}, \code{alpha10}, and \code{alpha00}.
#' @param significance_upper A numeric value specifying the target size or
#'   false discovery rate (FDR) level.
#' @param control.method A character string specifying the error-control
#'   target. Supported values are \code{"size"} and \code{"FDR"}.
#' @param ifcond A numeric indicator controlling the calibration rule used
#'   in the thresholding step under FDR control.
#'
#' @return
#' An integer vector containing the indices of selected hypotheses.
#'
#' @details
#' This function computes an MDACT-style composite statistic based on the
#' estimated proportions of the three null components and then applies the
#' thresholding rule implemented in [DACT_thr()] to determine the
#' rejection set.
#'
#' This procedure is included for comparison with the proposed tail
#' likelihood ratio (TLR) method in large-scale mediation testing.
#'
#' @export

MDACT = function(p.M,p.Y,estws,significance_upper,control.method,ifcond){
  sim.num <- length(p.M)

  pi.01.est = max(1e-3, estws$alpha01)
  pi.10.est = max(1e-3, estws$alpha10)
  pi.00.est = max(1e-3, estws$alpha00)
  c <- pi.01.est + pi.10.est + pi.00.est

  pi.01.est.sd = pi.01.est/c
  pi.10.est.sd = pi.10.est/c
  pi.00.est.sd = pi.00.est/c

  Ts <- pi.01.est*p.M + pi.10.est*p.Y + pi.00.est*pmax(p.M,p.Y)^2

  ss <- DACT_thr(p.M, p.Y,pi.10.est,pi.01.est,pi.00.est,significance_upper,control.method,ifcond)

  return(ss)
}
