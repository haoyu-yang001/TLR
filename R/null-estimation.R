#' Estimate composite-null mixture proportions from paired p-values
#'
#' Estimates the proportions of the three null components in the
#' composite null hypothesis of no mediation effect, together with the
#' two marginal null proportions, from paired p-values.
#'
#' @param input_pvalues A numeric matrix or data frame with exactly two
#'   columns of p-values. The first column should contain p-values for
#'   testing \eqn{H_{0\cdot}: \gamma_j = 0}, and the second column should
#'   contain p-values for testing \eqn{H_{\cdot 0}: \beta_j = 0}.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{\code{alpha10}}{Estimated proportion of the partial-null
#'   configuration \eqn{H_{\gamma 0}}.}
#'   \item{\code{alpha01}}{Estimated proportion of the partial-null
#'   configuration \eqn{H_{0 \beta}}.}
#'   \item{\code{alpha00}}{Estimated proportion of the complete null
#'   configuration \eqn{H_{00}}.}
#'   \item{\code{alpha1}}{Estimated marginal null proportion for
#'   \eqn{\gamma_j = 0}.}
#'   \item{\code{alpha2}}{Estimated marginal null proportion for
#'   \eqn{\beta_j = 0}.}
#' }
#'
#' @details
#' In large-scale causal mediation testing, each hypothesis is associated
#' with a pair of marginal p-values corresponding to the
#' exposure--mediator and mediator--outcome associations. The composite
#' null hypothesis of no mediation effect consists of three null cases:
#' \eqn{H_{00}}, \eqn{H_{\gamma 0}}, and \eqn{H_{0 \beta}}. This function
#' estimates the proportions of these three components, along with the
#' marginal null proportions, using tail-area summaries of the paired
#' p-value distribution and a QQ-slope-based selection rule.
#'
#' These estimated proportions serve as plug-in quantities for downstream
#' procedures such as HDMT, MDACT, and the proposed TLR method.
#'
#'
#' @export

null_estimation <- function (input_pvalues){
  if (is.null(ncol(input_pvalues)))
    stop("input_pvalues should be a matrix or data frame")
  if (ncol(input_pvalues) != 2)
    stop("inpute_pvalues should have 2 column")
  input_pvalues <- matrix(as.numeric(input_pvalues), nrow = nrow(input_pvalues))
  if (sum(complete.cases(input_pvalues)) < nrow(input_pvalues))
    warning("input_pvalues contains NAs to be removed from analysis")
  input_pvalues <- input_pvalues[complete.cases(input_pvalues),
  ]
  if (!is.null(nrow(input_pvalues)) & nrow(input_pvalues) <
      1)
    stop("input_pvalues doesn't have valid p-values")
  pcut <- seq(0.1, 0.8, 0.1)
  frac1 <- rep(0, 8)
  frac2 <- rep(0, 8)
  frac12 <- rep(0, 8)
  for (i in 1:8) {
    frac1[i] <- mean(input_pvalues[, 1] >= pcut[i])/(1 -
                                                       pcut[i])
    frac2[i] <- mean(input_pvalues[, 2] >= pcut[i])/(1 -
                                                       pcut[i])
    frac12[i] <- mean(input_pvalues[, 2] >= pcut[i] & input_pvalues[,
                                                                    1] >= pcut[i])/(1 - pcut[i])^2
  }
  alphaout <- matrix(0, 4, 5)
  ll <- 1
  qqslope <- rep(0, 4)
  for (lambda in c(0.5, 0.6, 0.7, 0.8)) {

    alpha00 <- min(frac12[pcut >= lambda][1], 1)
    if (stats::ks.test(input_pvalues[, 1], "punif", 0, 1, alternative = "greater")$p >
        0.05){
      alpha1 <- 1
    }else{
      alpha1 <- min(frac1[pcut >= lambda][1], 1)}

    if (stats::ks.test(input_pvalues[, 2], "punif", 0, 1, alternative = "greater")$p >
        0.05){
      alpha2 <- 1}else{
        alpha2 <- min(frac2[pcut >= lambda][1], 1)}


    if (alpha00 == 1) {
      alpha01 <- 0
      alpha10 <- 0
      alpha11 <- 0
    }else {
      if (alpha1 == 1 & alpha2 == 1) {
        alpha01 <- 0
        alpha10 <- 0
        alpha11 <- 0
        alpha00 <- 1
      }
      if (alpha1 == 1 & alpha2 != 1) {
        alpha10 <- 0
        alpha11 <- 0
        alpha01 <- alpha1 - alpha00
        alpha01 <- max(0, alpha01)
        alpha00 <- 1 - alpha01
      }
      if (alpha1 != 1 & alpha2 == 1) {
        alpha01 <- 0
        alpha11 <- 0
        alpha10 <- alpha2 - alpha00
        alpha10 <- max(0, alpha10)
        alpha00 <- 1 - alpha10
      }
      if (alpha1 != 1 & alpha2 != 1) {
        alpha10 <- alpha2 - alpha00
        alpha10 <- max(0, alpha10)
        alpha01 <- alpha1 - alpha00
        alpha01 <- max(0, alpha01)
        if ((1 - alpha00 - alpha01 - alpha10) < 0) {
          alpha11 <- 0
          alpha10 <- 1 - alpha1
          alpha01 <- 1 - alpha2
          alpha00 <- 1 - alpha10 - alpha01
        }
        else {
          alpha11 <- 1 - alpha00 - alpha01 - alpha10
        }
      }
    }

    pmax <- apply(input_pvalues, 1, max)
    pmax <- pmax[order(pmax)]
    nnulls <- sum(pmax > 0.8)
    nmed <- nrow(input_pvalues)
    pexp <- rep(0, nnulls)

    b <- alpha01 + alpha10
    a <- 1 - b
    if(a == 0){
      alpha.null <- list(alpha10 = alpha10,alpha01 = alpha01,
                         alpha00 = alpha00, alpha1 = alpha1, alpha2 = alpha2)
      return(alpha.null)
    }

    for (i in 1:nmed) {
      c <- (-i/nmed)
      pexp[i] <- (-b + sqrt(b^2 - 4 * a * c))/(2 * a)
    }

    xx <- -log(pexp[(nmed - nnulls + 1):nmed], base = 10)
    yy <- -log(pmax[(nmed - nnulls + 1):nmed], base = 10)
    fit1 <- stats::lm(yy ~ xx - 1)
    qqslope[ll] <- fit1$coef[1]

    alphaout[ll, 1] <- alpha10
    alphaout[ll, 2] <- alpha01
    alphaout[ll, 3] <- alpha00
    alphaout[ll, 4] <- alpha1
    alphaout[ll, 5] <- alpha2
    ll <- ll + 1
  }
  bestslope <- which.min(qqslope)
  alpha.null <- list(alpha10 = alphaout[bestslope, 1], alpha01 = alphaout[bestslope,
                                                                          2], alpha00 = alphaout[bestslope, 3], alpha1 = alphaout[bestslope,
                                                                                                                                  4], alpha2 = alphaout[bestslope, 5])
  return(alpha.null)
}
