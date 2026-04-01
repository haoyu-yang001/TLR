#' Summarize FDR and power from a selected rejection set
#'
#' Computes false discovery rate (FDR) and power by comparing a selected
#' set of hypotheses with the true signal set in simulation studies.
#'
#' @param sat_index An integer vector of indices selected by a testing
#'   procedure.
#' @param true_set An integer vector of indices corresponding to the true
#'   non-null hypotheses.
#'
#' @return
#' If \code{true_set} is non-empty, returns a numeric vector of length two
#' containing:
#' \describe{
#'   \item{FDR}{The empirical false discovery rate.}
#'   \item{power}{The empirical power.}
#' }
#' If \code{true_set} is empty, returns a numeric scalar indicating
#' whether at least one false positive is made.
#'
#' @details
#' This function is designed for simulation studies in large-scale
#' mediation testing, where one compares the rejection set produced by a
#' method with the known truth. When there are no true signals, the
#' function returns a false-positive indicator instead of FDR and power.
#'
#' @export
index_fpn_fdr_power_func <- function(sat_index,true_set){
  if(length(true_set)==0){
    fpn <- ifelse(length(sat_index)>0,1,0)
    return(fpn)
  }else{
    fdr <-length(which(sat_index %!in% true_set))/max(length(sat_index),1)
    power <- length(which(sat_index %in% true_set))/length(true_set)
    return(c(fdr,power))
  }
}

#' Summarize size and power from a selected rejection set
#'
#' Computes empirical size and power by comparing a selected set of
#' hypotheses with the true signal set in simulation studies.
#'
#' @param sat_index An integer vector of indices selected by a testing
#'   procedure.
#' @param true_set An integer vector of indices corresponding to the true
#'   non-null hypotheses.
#' @param sim.num A positive integer giving the total number of simulated
#'   hypotheses.
#'
#' @return
#' If \code{true_set} is non-empty, returns a numeric vector of length two
#' containing:
#' \describe{
#'   \item{size}{The empirical size, computed as the proportion of false
#'   rejections among all hypotheses.}
#'   \item{power}{The empirical power.}
#' }
#' If \code{true_set} is empty, returns the total number of false
#' positives.
#'
#' @details
#' This function is useful in simulation studies where one wants to report
#' size instead of false discovery rate. In the context of large-scale
#' mediation testing, it can be used to summarize the error-rate and power
#' performance of competing methods under the known data-generating
#' mechanism.
#'
#' @export
size_func <- function(sat_index,true_set,sim.num){
  if(length(true_set)==0){
    fpn <- length(sat_index)
    return(fpn)
  }else{
    size <- length(which(sat_index %!in% true_set))/sim.num
    power <- length(which(sat_index %in% true_set))/length(true_set)
    return(c(size,power))
  }
}


#' Simulate paired p-values under a four-component mixture model
#'
#' Generates paired Z-statistics and corresponding two-sided p-values
#' under a four-component mixture model for large-scale mediation testing.
#'
#' @param sim.num A positive integer giving the total number of simulated
#'   hypotheses.
#' @param ws A numeric vector of length four giving the mixture
#'   proportions of the four components. These correspond to the two
#'   marginal null/non-null configurations:
#'   null/non-null in the first margin, null/non-null in the second
#'   margin, complete null, and simultaneous non-null.
#' @param m A numeric value controlling the overall signal magnitude.
#' @param delta A numeric value controlling the relative signal scaling
#'   between the two margins.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{\code{pvalues}}{A two-column matrix of paired p-values.}
#'   \item{\code{zvalues}}{A two-column matrix of paired Z-statistics.}
#' }
#'
#' @details
#' This function simulates paired marginal test statistics under a
#' four-component mixture structure, which is useful for evaluating
#' large-scale mediation testing procedures under controlled settings.
#' Signal strengths in the two margins are governed jointly by
#' \code{m} and \code{delta}, and two-sided p-values are computed from the
#' simulated Z-statistics.
#'
#' The simulated output can be used to compare the operating
#' characteristics of TLR, HDMT, DACT, Sobel, and related methods in terms
#' of error-rate control and power.
#'
#' @export
generate_pvalues_indp_from_norm <- function(sim.num,ws,m,delta){
  pras <- c(m/sqrt(1+delta^2),m/sqrt(1+delta^2)*delta)
  combinations <- as.matrix(expand.grid(rep(list(c(0, 1)),2)))[c(2,3,1,4),]
  # Calculate the number of units to be assigned to each set based on probabilities
  sim.num_per_set <- round(ws * sim.num)
  sets <- vector("list", length = length(ws))
  current_index <- 1
  for (i in 1:length(ws)) {
    if (sim.num_per_set[i] > 0) {
      sets[[i]] <- current_index:(current_index + sim.num_per_set[i] - 1)
      current_index <- current_index + sim.num_per_set[i]
    } else {
      sets[[i]] <- integer(0)  # Empty set
    }
  }
  mu_M <- sample(c(1,-1),length(c(unlist(sets[1]),unlist(sets[4]))),replace = T)*rnorm(length(c(unlist(sets[1]),unlist(sets[4]))),pras[1],1)
  mu_Y <- sample(c(1,-1),length(c(unlist(sets[2]),unlist(sets[4]))),replace = T)*rnorm(length(c(unlist(sets[2]),unlist(sets[4]))),pras[2],1)
  Z.M <- rnorm(sim.num,0,1); Z.M[c(unlist(sets[1]),unlist(sets[4]))] <- rnorm(length(c(unlist(sets[1]),unlist(sets[4]))),mu_M,1);
  Z.Y <- rnorm(sim.num,0,1); Z.Y[c(unlist(sets[2]),unlist(sets[4]))] <- rnorm(length(c(unlist(sets[2]),unlist(sets[4]))),mu_Y,1);
  p.M <- 2* (1 - stats::pnorm(abs(Z.M))); p.Y <- 2* (1 - stats::pnorm(abs(Z.Y)))
  p.M[p.M==0] <- 1e-17 ; p.Y[p.Y==0] <- 1e-17
  return(list(pvalues = cbind(p.M,p.Y),zvalues=cbind(Z.M,Z.Y)))
}
