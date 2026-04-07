#' Estimate tail approximation parameters for the TLR statistic
#'
#' Estimates the left-tail approximation parameters used in the tail
#' likelihood ratio (TLR) procedure for large-scale mediation testing.
#' Under the paper's tail model, the p-value distribution under the
#' alternative is approximated near zero by
#' \eqn{F_1(t) \approx C t^\alpha}, which implies the density approximation
#' \eqn{f_1(t) \approx C \alpha t^{\alpha - 1}}. This function estimates
#' the tail exponents \eqn{\alpha_1, \alpha_2} and the corresponding tail
#' constants \eqn{C_1, C_2} from the paired p-values.
#'
#' @param X A two-column numeric matrix of paired p-values. The first
#'   column contains \eqn{p_{\gamma,j}} and the second column contains
#'   \eqn{p_{\beta,j}}.
#' @param estws A list of estimated null proportions returned by the
#'   null-proportion estimation step. It should contain at least
#'   \code{alpha1} and \code{alpha2}, the estimated marginal null
#'   proportions for \eqn{\gamma_j = 0} and \eqn{\beta_j = 0}.
#' @param m0 An integer specifying the tail cutoff. The \code{m0}-th
#'   order statistic in each p-value column is used to define the region
#'   over which the tail parameters are estimated.
#'
#' @return A numeric vector of length 4 containing
#'   \code{c(alpha1_hat, alpha2_hat, C1, C2)}.
#'
#' @details
#' The TLR method approximates the alternative p-value distributions in
#' the extreme left tail rather than modeling the full distribution over
#' \eqn{(0,1)}. This function implements the plug-in estimation step for
#' the tail parameters used in the TLR ranking statistic.
#'
#'
#' @export

alphas_estimate_tail <- function(X,estws,m0){
  X  <- as.matrix(X)
  input_pvalues = X
  nullprop <- estws

  m <- nrow(X)
  m0 <- ceiling(m0)
  r1 <- sort(X[, 1])[m0]
  r2 <- sort(X[, 2])[m0]

  solve.alpha <- function(x) {
    sum(log(X[X[, 1] <= r1, 1])) / m - min(0.99, nullprop$alpha1) * (r1*log(r1)-r1) - ((m0 / m  - r1 * nullprop$alpha1) / r1^{x} / (1 - nullprop$alpha1))*(1 - min(0.99, nullprop$alpha1)) * (r1^{x} * log(r1) - x^{-1} * r1^{x})
  }

  alpha1_hat <- try(uniroot(solve.alpha, c(1e-3, 1))$root, silent = TRUE)
  if(inherits(alpha1_hat, "try-error")){
    x_vals <- seq(0.001, 1, by = 0.01)
    y_vals <- sapply(x_vals, solve.alpha)
    closest_index <- which.min(abs(y_vals))
    alpha1_hat <- x_vals[closest_index]
  }

  solve.alpha <- function(x) {
    sum(log(X[X[, 2] <= r2, 2])) / m - min(0.99, nullprop$alpha2) * (r2*log(r2)-r2) - ((m0 / m  - r2 * nullprop$alpha2) / r2^{x} / (1 - nullprop$alpha2))*(1 - min(0.99, nullprop$alpha2)) * (r2^{x} * log(r2) - x^{-1} * r2^{x})
  }

  alpha2_hat <- try(uniroot(solve.alpha, c(1e-3, 1))$root, silent = TRUE)
  if(inherits(alpha2_hat, "try-error")){
    x_vals <- seq(0.001, 1, by = 0.01)
    y_vals <- sapply(x_vals, solve.alpha)
    closest_index <- which.min(abs(y_vals))
    alpha2_hat <- x_vals[closest_index]
  }

  C1 <- max(0.1,min((m0 / m  - r1 * nullprop$alpha1) / r1^{alpha1_hat} / (1 - nullprop$alpha1),3))
  C2 <- max(0.1,min((m0 / m  - r2 * nullprop$alpha2) / r2^{alpha2_hat} / (1 - nullprop$alpha2),3))

  return(c(alpha1_hat,alpha2_hat, C1, C2))
}

#' Compute an LFDR-style ranking statistic for composite-null mediation testing
#'
#' Computes a local-false-discovery-rate-style ranking statistic for the
#' composite null hypothesis of no mediation effect. The statistic is
#' formed from the likelihood-ratio representation of the composite-null
#' posterior probability and uses beta densities as working models for
#' the alternative p-value distributions.
#'
#' @param p.M A numeric vector of p-values for testing
#'   \eqn{H_{0\cdot}: \gamma_j = 0}.
#' @param p.Y A numeric vector of p-values for testing
#'   \eqn{H_{\cdot 0}: \beta_j = 0}.
#' @param estws A list of estimated null proportions. It should contain
#'   \code{alpha00}, \code{alpha10}, and \code{alpha01}, corresponding to
#'   the estimated proportions of the three null cases.
#' @param alpha_true A positive numeric value specifying the first shape
#'   parameter of the working beta alternative density.
#' @param beta_true A positive numeric value specifying the second shape
#'   parameter of the working beta alternative density.
#'
#' @return A numeric vector of LFDR-style ranking statistics, one for
#'   each p-value pair.
#'
#' @details
#' For the composite null
#' \eqn{H_{0,j}: \gamma_j \beta_j = 0}, the paper expresses the local FDR
#' as a monotone function of a likelihood ratio involving the three null
#' components and the alternative component. This function evaluates a
#' working version of that ranking quantity using beta models for the
#' marginal alternative p-value densities.
#'
#' Smaller values of the returned statistic correspond to stronger
#' evidence against the composite null.
#'
#'
#' @export

lfdr_stat <- function(p.M, p.Y, estws, alpha_true, beta_true){
  B_val <- beta(alpha_true, beta_true)
  f1gamma <- 1/B_val * p.M^{alpha_true - 1} * (1-p.M)^{beta_true - 1}
  f1beta <- 1/B_val * p.Y^{alpha_true - 1} * (1-p.Y)^{beta_true - 1}
  lfdr_statistic <- estws$alpha00/(f1gamma * f1beta) + estws$alpha10/f1beta + estws$alpha01/f1gamma
  return(lfdr_statistic)
}

#' Compute an LFDR-based threshold for composite-null mediation testing
#'
#' Computes a threshold on an LFDR-style ranking statistic for paired
#' mediation p-values under the composite null hypothesis of no mediation
#' effect.
#'
#' @param p.M A numeric vector of p-values for testing
#'   \eqn{H_{0\cdot}: \gamma_j = 0}.
#' @param p.Y A numeric vector of p-values for testing
#'   \eqn{H_{\cdot 0}: \beta_j = 0}.
#' @param estws A list of estimated null proportions. It should contain
#'   \code{alpha00}, \code{alpha10}, and \code{alpha01}.
#' @param alpha_true A positive numeric value specifying the first shape
#'   parameter of the working beta alternative density.
#' @param beta_true A positive numeric value specifying the second shape
#'   parameter of the working beta alternative density.
#' @param sets A list of index sets used to evaluate the empirical null
#'   proportion among selected hypotheses. The first three elements are
#'   assumed to correspond to the three null components.
#' @param sig_level A numeric value specifying the target upper bound for
#'   the empirical null proportion among selected hypotheses.
#'
#' @return A numeric threshold value on the LFDR-style ranking statistic.
#'
#' @details
#' This function computes the LFDR-style ranking statistic from the
#' paired p-values and then uses a bisection search over the sorted
#' statistic values to identify a cutoff satisfying the target empirical
#' error criterion determined by \code{sets} and \code{sig_level}.
#'
#' This function is primarily useful for simulation studies or
#' comparisons with LFDR-based rules; the main method proposed in the
#' paper is the TLR procedure based on tail approximation and the
#' composite-null distribution of the TLR statistic.
#'
#'
#' @export

lfdr_thre <- function(p.M, p.Y, estws, alpha_true, beta_true,sets,sig_level){
  B_val <- beta(alpha_true, beta_true)
  f1gamma <- 1/B_val * p.M^{alpha_true - 1} * (1-p.M)^{beta_true - 1}
  f1beta <- 1/B_val * p.Y^{alpha_true - 1} * (1-p.Y)^{beta_true - 1}
  lfdr_statistic <- estws$alpha00/(f1gamma * f1beta) + estws$alpha10/f1beta + estws$alpha01/f1gamma

  sorted_index <- order(lfdr_statistic)
  lower_bound <- 1
  upper_bound <- length(lfdr_statistic)

  thr <- function(lfdr_statistic,sets,sig_level,t) {
    selected <- which(lfdr_statistic < t)
    if (length(selected) == 0) return(Inf)
    prop <- length(intersect(selected, unlist(sets[1:3]))) / length(selected)
    prop - sig_level
  }

  while (upper_bound - lower_bound > 1) {
    mid_index <- floor((lower_bound + upper_bound) / 2)
    mid_value <- lfdr_statistic[sorted_index[mid_index]]
    if (thr(lfdr_statistic,sets,sig_level,mid_value) < 0) {
      lower_bound <- mid_index
    } else {
      upper_bound <- mid_index
    }
  }
  t_s <- lfdr_statistic[sorted_index[lower_bound]]
  return(t_s)
}

#' Tail likelihood ratio multiple-testing procedure for mediation analysis
#'
#' Applies the tail likelihood ratio (TLR) procedure for large-scale
#' causal mediation testing and returns the selected hypotheses.
#'
#' @param p.M A numeric vector of p-values for testing
#'   \eqn{H_{0\cdot}: \gamma_j = 0}.
#' @param p.Y A numeric vector of p-values for testing
#'   \eqn{H_{\cdot 0}: \beta_j = 0}.
#' @param estws A list of estimated proportions of the three composite-null
#'   components. It should contain \code{alpha01}, \code{alpha10}, and
#'   \code{alpha00}.
#' @param coefs A numeric vector containing the estimated tail exponents,
#'   typically \eqn{(\alpha_1, \alpha_2)}.
#' @param coefs2 A numeric vector containing the estimated tail constants,
#'   typically \eqn{(C_1, C_2)}.
#' @param significance_upper A numeric value specifying the target error
#'   control level, such as the desired FDR level or size level.
#' @param control.method A character string specifying the type of error
#'   control. Supported values are \code{"FDR"} and \code{"size"}.
#' @param ifcond An indicator controlling which calibration rule is used
#'   in the thresholding step under FDR control.
#'
#' @return An integer vector of selected hypothesis indices.
#'
#' @details
#' The TLR method ranks hypotheses using a tail likelihood ratio statistic
#' derived from the composite-null likelihood ratio, where the marginal
#' alternative p-value densities are approximated in the left tail by
#' \eqn{f_1(t) \approx C \alpha t^{\alpha - 1}}. The rejection threshold
#' is then determined using the derived null distribution of the TLR
#' statistic under the composite null, which improves robustness across
#' sparse and dense signal settings. See Sections 3 and 4 of the paper.
#'
#' @export

TLR <- function(p.M,p.Y,estws,coefs,coefs2,significance_upper,control.method,ifcond){
  sim.num <- length(p.M)

  pi.01.est = max(1e-3, estws$alpha01)
  pi.10.est = max(1e-3, estws$alpha10)
  pi.00.est = max(1e-3, estws$alpha00)
  c <- pi.01.est + pi.10.est + pi.00.est
  pi.11.est = 1-c
  pi.01.est.sd = pi.01.est/c
  pi.10.est.sd = pi.10.est/c
  pi.00.est.sd = pi.00.est/c

  Ts <- pi.10.est * coefs2[2]^(-1) * coefs[2]^(-1) * (p.Y^(1-coefs[2])) +
    pi.01.est * coefs2[1]^(-1) * coefs[1]^(-1) * (p.M^(1-coefs[1])) +
    pi.00.est * coefs2[1]^(-1) * coefs2[2]^(-1) * coefs[1]^(-1) * coefs[2]^(-1) * (p.M^(1-coefs[1]) * p.Y^(1-coefs[2]))

  ss <- TLR_thre(p.M, p.Y,coefs,coefs2,pi.10.est,pi.01.est,pi.00.est,significance_upper,control.method,ifcond)

  return(ss)
}

#' Thresholding step for TLR
#'
#' Internal worker used by [TLR()] to compute the rejection set.
#'
#' @param p.M,p.Y Numeric vectors of p-values.
#' @param coefs,coefs2 Tail-approximation parameters.
#' @param pi.10.est,pi.01.est,pi.00.est Estimated null proportions.
#' @param significance_upper Target size or FDR level.
#' @param control.method Either `"size"` or `"FDR"`.
#' @param ifcond Numeric indicator controlling conditional calibration.
#'
#' @return Integer indices of selected hypotheses.
#' @export

TLR_thre <- function(p.M, p.Y,coefs,coefs2,pi.10.est,pi.01.est,pi.00.est,significance_upper,control.method,ifcond){
  sim.num <- length(p.M)
  c <- pi.01.est + pi.10.est + pi.00.est
  pi.11.est <- max(1 - c, 0)

  Ts <- pi.10.est * coefs2[2]^(-1) * coefs[2]^(-1) * (p.Y^(1-coefs[2])) +
    pi.01.est * coefs2[1]^(-1) * coefs[1]^(-1) * (p.M^(1-coefs[1])) +
    pi.00.est * coefs2[1]^(-1) * coefs2[2]^(-1) * coefs[1]^(-1) * coefs[2]^(-1) * (p.M^(1-coefs[1]) * p.Y^(1-coefs[2]))

  ER00 <- function(t) {
    s <- max((t-pi.01.est*coefs2[1]^(-1)/coefs[1])/(pi.00.est*coefs2[1]^(-1)*coefs2[2]^(-1)/coefs[1]+pi.10.est*coefs2[2]^(-1)),0)
    f <- function(x, pi.00.est, pi.01.est, pi.10.est, coefs,coefs2, t) {
      c <- coefs[1]^(1/(1-coefs[1]))*coefs[2]^(1/(1-coefs[2]))*(1/(1-coefs[2]))
      integrand <- ((t-pi.10.est*coefs2[2]^(-1)*x)/(pi.00.est*coefs2[1]^(-1)*coefs2[2]^(-1)*x+pi.01.est*coefs2[1]^(-1)))^(1/(1-coefs[1])) * x ^(coefs[2]/(1-coefs[2]))
      return(integrand*c)
    }
    up <- min((t/(pi.10.est*coefs2[2]^(-1))),1/coefs[2])
    re2 <- (s*coefs[2])^(1/(1-coefs[2]))
    if(s<=up){
      re1 <- integrate(f, lower = s, upper = up, pi.00.est = pi.00.est, pi.01.est = pi.01.est, pi.10.est = pi.10.est, coefs=coefs, coefs2=coefs2, t = t)
      return(re1$value + re2)
    }else{
      return(1)
    }
  }

  ER01 <- function(t) {
    s <- pmin(pmax((t - pi.10.est * coefs2[2]^(-1) * coefs[2]^{-1} * p.Y^{1 - coefs[2]}) /
                     (pi.01.est * coefs2[1]^(-1) * coefs[1]^(-1) +
                        pi.00.est * coefs2[1]^(-1) * coefs2[2]^(-1) * coefs[1]^(-1) * coefs[2]^(-1) * p.Y^{1 - coefs[2]}), 0), 1)
    mean(s^{1 / (1 - coefs[1])})
  }

  ER10 <- function(t) {
    s <- pmin(pmax((t - pi.01.est * coefs2[1]^(-1) * coefs[1]^{-1} * p.M^{1 - coefs[1]}) /
                     (pi.10.est * coefs2[2]^(-1) * coefs[2]^(-1) +
                        pi.00.est * coefs2[1]^(-1) * coefs2[2]^(-1) * coefs[1]^(-1) * coefs[2]^(-1) * p.M^{1 - coefs[1]}), 0), 1)
    mean(s^{1 / (1 - coefs[2])})
  }

  if(control.method == "FDR"){
    if(ifcond == 0){
      thr <- function(t) {
        c * (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*ER10(t) +
                                pi.01.est/(pi.01.est+pi.11.est)*ER01(t) +
                                (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) -
                                   pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*ER00(t))) / max(mean(Ts < t), 1 / sim.num) - significance_upper
      }
    }
    if(ifcond == 1){
      if(length(which(p.adjust(p.M,method = "fdr") < 0.1)) < 50 | length(which(p.adjust(p.Y,method = "fdr") < 0.1)) < 50){
        est10 <- 1
        est01 <- 1
        thr <- function(t) {
          (ER10(t) + ER01(t) - ER00(t)) / max(mean(Ts < t), 1 / sim.num) - significance_upper
        }
      }else{
        est10 <- mean(p.adjust(p.M,method = "fdr") < 0.1 & p.Y > 0.9)/stats::punif(1-0.9)/mean(p.adjust(p.M,method = "fdr") < 0.1)
        est01 <- mean(p.adjust(p.Y,method = "fdr") < 0.1 & p.M > 0.9)/stats::punif(1-0.9)/mean(p.adjust(p.Y,method = "fdr") < 0.1)
      }
      thr <- function(t) {
        (est10*ER10(t) + est01*ER01(t) +
           (pi.00.est - est10*(pi.00.est + pi.01.est) - est01*(pi.00.est + pi.10.est))*ER00(t)) / max(mean(Ts < t), 1 / sim.num) - significance_upper
      }
    }
    #t_s <- uniroot(thr, c(min(c(p.M, p.Y)), 1))$root
    sorted_index <- order(Ts) # Get the indices that sort 'Ts'
    lower_bound <- 1
    upper_bound <- length(Ts)

    while (upper_bound - lower_bound > 1) {
      mid_index <- floor((lower_bound + upper_bound) / 2)
      mid_value <- Ts[sorted_index[mid_index]]

      if (thr(mid_value) < 0) {
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
      (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*ER10(t) +
                          pi.01.est/(pi.01.est+pi.11.est)*ER01(t) +
                          (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) -
                             pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*ER00(t))) - significance_upper
    }
    x <- tryCatch(uniroot(thr, c(min(c(p.M, p.Y)), 1))$root, error = function(e) e)

    F_emp <- function(t) {
      F_emp <- (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*ER10(t) +
                                   pi.01.est/(pi.01.est+pi.11.est)*ER01(t) +
                                   (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) -
                                      pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*ER00(t)))
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

        if (F_emp(mid_value) < significance_upper) {
          lower_bound <- mid_index
        } else {
          upper_bound <- mid_index
        }
      }
      t_s <- Ts[sorted_index[lower_bound]]
    }

    if(F_emp(round(t_s,6)) < significance_upper){
      return(which(Ts <= t_s))
    }else{
      return(which(Ts < t_s))
    }
  }
}

#' TLR p-values
#'
#' Computes p-values associated with the TLR statistic.
#'
#' @param p.M,p.Y Marginal p-values used to estimate partial-null components.
#' @param p.M_one,p.Y_one P-values used in the TLR statistic.
#' @param estws Estimated null proportions.
#' @param coefs,coefs2 Tail-approximation parameters.
#'
#' @return A numeric vector of p-values.
#' @export

TLR_pvalues <- function(p.M,p.Y,p.M_one,p.Y_one,estws,coefs,coefs2){
  sim.num <- length(p.M)

  pi.01.est = max(1e-3, estws$alpha01)
  pi.10.est = max(1e-3, estws$alpha10)
  pi.00.est = max(1e-3, estws$alpha00)
  c <- pi.01.est + pi.10.est + pi.00.est
  pi.11.est <- 1-min(c,0.5)

  Ts <- pi.10.est * coefs2[2]^(-1) * coefs[2]^(-1) * (p.Y_one^(1-coefs[2])) +
    pi.01.est * coefs2[1]^(-1) * coefs[1]^(-1) * (p.M_one^(1-coefs[1])) +
    pi.00.est * coefs2[1]^(-1) * coefs2[2]^(-1) * coefs[1]^(-1) * coefs[2]^(-1) * (p.M_one^(1-coefs[1]) * p.Y_one^(1-coefs[2]))

  ER00 <- function(t) {
    s <- max((t-pi.01.est*coefs2[1]^(-1)/coefs[1])/(pi.00.est*coefs2[1]^(-1)*coefs2[2]^(-1)/coefs[1]+pi.10.est*coefs2[2]^(-1)),0)
    f <- function(x, pi.00.est, pi.01.est, pi.10.est, coefs,coefs2, t) {
      c <- coefs[1]^(1/(1-coefs[1]))*coefs[2]^(1/(1-coefs[2]))*(1/(1-coefs[2]))
      integrand <- ((t-pi.10.est*coefs2[2]^(-1)*x)/(pi.00.est*coefs2[1]^(-1)*coefs2[2]^(-1)*x+pi.01.est*coefs2[1]^(-1)))^(1/(1-coefs[1])) * x ^(coefs[2]/(1-coefs[2]))
      return(integrand*c)
    }
    up <- min((t/(pi.10.est*coefs2[2]^(-1))),1/coefs[2])
    re2 <- (s*coefs[2])^(1/(1-coefs[2]))
    if(s<=up){
      re1 <- integrate(f, lower = s, upper = up, pi.00.est = pi.00.est, pi.01.est = pi.01.est, pi.10.est = pi.10.est, coefs=coefs, coefs2=coefs2, t = t)
      return(re1$value + re2)
    }else{
      return(1)
    }
  }

  ER01 <- function(t) {
    s <- pmin(pmax((t - pi.10.est * coefs2[2]^(-1) * coefs[2]^{-1} * p.Y^{1 - coefs[2]}) /
                     (pi.01.est * coefs2[1]^(-1) * coefs[1]^(-1) +
                        pi.00.est * coefs2[1]^(-1) * coefs2[2]^(-1) * coefs[1]^(-1) * coefs[2]^(-1) * p.Y^{1 - coefs[2]}), 0), 1)
    mean(s^{1 / (1 - coefs[1])})
  }

  ER10 <- function(t) {
    s <- pmin(pmax((t - pi.01.est * coefs2[1]^(-1) * coefs[1]^{-1} * p.M^{1 - coefs[1]}) /
                     (pi.10.est * coefs2[2]^(-1) * coefs[2]^(-1) +
                        pi.00.est * coefs2[1]^(-1) * coefs2[2]^(-1) * coefs[1]^(-1) * coefs[2]^(-1) * p.M^{1 - coefs[1]}), 0), 1)
    mean(s^{1 / (1 - coefs[2])})
  }

  Fnull <- function(t) {
    (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*ER10(t) +
                        pi.01.est/(pi.01.est+pi.11.est)*ER01(t) +
                        (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) -
                           pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*ER00(t)))
  }

  pv <- vapply(Ts, Fnull, numeric(1))

  return(pv)
}

#' TLR p-values with trimmed empirical moments
#'
#' Computes TLR p-values using trimmed means in the partial-null components.
#'
#' @param p.M,p.Y Marginal p-values used to estimate component expectations.
#' @param p.M_one,p.Y_one P-values used in the TLR statistic.
#' @param estws Estimated null proportions.
#' @param coefs,coefs2 Tail-approximation parameters.
#'
#' @return A numeric vector of p-values.
#' @export

TLR_pvalues_trimmed <- function(p.M,p.Y,p.M_one,p.Y_one,estws,coefs,coefs2){
  sim.num <- length(p.M)

  pi.01.est = max(1e-3, estws$alpha01)
  pi.10.est = max(1e-3, estws$alpha10)
  pi.00.est = max(1e-3, estws$alpha00)
  c <- pi.01.est + pi.10.est + pi.00.est
  pi.11.est <- 1-min(c,0.5)

  Ts <- pi.10.est * coefs2[2]^(-1) * coefs[2]^(-1) * (p.Y_one^(1-coefs[2])) +
    pi.01.est * coefs2[1]^(-1) * coefs[1]^(-1) * (p.M_one^(1-coefs[1])) +
    pi.00.est * coefs2[1]^(-1) * coefs2[2]^(-1) * coefs[1]^(-1) * coefs[2]^(-1) * (p.M_one^(1-coefs[1]) * p.Y_one^(1-coefs[2]))

  ER00 <- function(t) {
    s <- max((t-pi.01.est*coefs2[1]^(-1)/coefs[1])/(pi.00.est*coefs2[1]^(-1)*coefs2[2]^(-1)/coefs[1]+pi.10.est*coefs2[2]^(-1)),0)
    f <- function(x, pi.00.est, pi.01.est, pi.10.est, coefs,coefs2, t) {
      c <- coefs[1]^(1/(1-coefs[1]))*coefs[2]^(1/(1-coefs[2]))*(1/(1-coefs[2]))
      integrand <- ((t-pi.10.est*coefs2[2]^(-1)*x)/(pi.00.est*coefs2[1]^(-1)*coefs2[2]^(-1)*x+pi.01.est*coefs2[1]^(-1)))^(1/(1-coefs[1])) * x ^(coefs[2]/(1-coefs[2]))
      return(integrand*c)
    }
    up <- min((t/(pi.10.est*coefs2[2]^(-1))),1/coefs[2])
    re2 <- (s*coefs[2])^(1/(1-coefs[2]))
    if(s<=up){
      re1 <- integrate(f, lower = s, upper = up, pi.00.est = pi.00.est, pi.01.est = pi.01.est, pi.10.est = pi.10.est, coefs=coefs, coefs2=coefs2, t = t)
      return(re1$value + re2)
    }else{
      return(1)
    }
  }

  ER01 <- function(t) {
    s <- pmin(pmax((t - pi.10.est * coefs2[2]^(-1) * coefs[2]^{-1} * p.Y^{1 - coefs[2]}) /
                     (pi.01.est * coefs2[1]^(-1) * coefs[1]^(-1) +
                        pi.00.est * coefs2[1]^(-1) * coefs2[2]^(-1) * coefs[1]^(-1) * coefs[2]^(-1) * p.Y^{1 - coefs[2]}), 0), 1)
    mean(s^{1 / (1 - coefs[1])}, trim = 0.1)
  }

  ER10 <- function(t) {
    s <- pmin(pmax((t - pi.01.est * coefs2[1]^(-1) * coefs[1]^{-1} * p.M^{1 - coefs[1]}) /
                     (pi.10.est * coefs2[2]^(-1) * coefs[2]^(-1) +
                        pi.00.est * coefs2[1]^(-1) * coefs2[2]^(-1) * coefs[1]^(-1) * coefs[2]^(-1) * p.M^{1 - coefs[1]}), 0), 1)
    mean(s^{1 / (1 - coefs[2])}, trim = 0.1)
  }

  Fnull <- function(t) {
    (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*ER10(t) +
                        pi.01.est/(pi.01.est+pi.11.est)*ER01(t) +
                        (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) -
                           pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*ER00(t)))
  }

  pv <- vapply(Ts, Fnull, numeric(1))

  return(pv)
}
