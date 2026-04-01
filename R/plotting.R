#' QQ plot against the uniform null distribution
#'
#' Draws a quantile-quantile (QQ) plot comparing observed p-values with
#' their expected values under the null uniform distribution. The plot can
#' be used to assess calibration of p-values, visualize departure from the
#' null, and examine enrichment of small p-values in the tail.
#'
#' @param pvalues A numeric vector of p-values, or a named list of numeric
#'   vectors of p-values. If a list is supplied, each element is plotted as
#'   a separate group.
#' @param should.thin Logical; if \code{TRUE}, nearby points are thinned by
#'   rounding the observed and expected coordinates before plotting. This
#'   is useful for large-scale testing problems with many p-values.
#' @param thin.obs.places Integer; number of decimal places used when
#'   rounding observed \eqn{-\log_{10}(p)} values for thinning.
#' @param thin.exp.places Integer; number of decimal places used when
#'   rounding expected \eqn{-\log_{10}(p)} values for thinning.
#' @param xlab Label for the x-axis. By default this is the expected
#'   \eqn{-\log_{10}(p)} under the uniform null distribution.
#' @param ylab Label for the y-axis. By default this is the observed
#'   \eqn{-\log_{10}(p)}.
#' @param draw.conf Logical; if \code{TRUE}, adds a pointwise confidence
#'   envelope under the null uniform distribution.
#' @param conf.points Integer; number of points used to draw the
#'   confidence envelope.
#' @param conf.col Color used to fill the confidence envelope.
#' @param conf.alpha Significance level used for the confidence envelope.
#' @param already.transformed Logical; if \code{FALSE}, the input is
#'   interpreted as raw p-values and transformed internally to
#'   \eqn{-\log_{10}(p)}. If \code{TRUE}, the input is assumed to already
#'   be on the \eqn{-\log_{10}(p)} scale.
#' @param pch Plotting character passed to [lattice::xyplot()].
#' @param aspect Aspect ratio passed to [lattice::xyplot()].
#' @param prepanel A prepanel function passed to [lattice::xyplot()]. The
#'   default sets matched x- and y-axis ranges.
#' @param par.settings A list of graphical parameters passed to
#'   [lattice::xyplot()].
#' @param ... Additional arguments passed to [lattice::xyplot()].
#'
#' @return
#' A \code{"trellis"} object produced by [lattice::xyplot()].
#'
#' @details
#' For a single vector of p-values, the function plots the observed
#' ordered \eqn{-\log_{10}(p)} values against their expected order
#' statistics under the null uniform distribution. If a list of p-value
#' vectors is provided, each list element is displayed as a separate
#' group, allowing visual comparison across methods or simulation
#' settings.
#'
#' In large-scale mediation testing and epigenome-wide studies, this plot
#' is useful for checking whether p-values are well calibrated under the
#' null and for visualizing excess small p-values that may indicate true
#' signals or model misspecification.
#'
#' @examples
#' set.seed(1)
#' p0 <- stats::runif(1000)
#' qqunif.plot(p0)
#'
#' p.list <- list(
#'   Method1 = stats::runif(1000),
#'   Method2 = c(stats::runif(950), stats::runif(50, 0, 1e-3))
#' )
#' qqunif.plot(p.list)
#'
#' @export

qqunif.plot<-function(pvalues,
                      should.thin=T, thin.obs.places=2, thin.exp.places=2,
                      xlab=expression(paste("Expected (",-log[10],~p[beta[j]],")")),
                      ylab=expression(paste("Observed (",-log[10],~p[beta[j]],")")),
                      draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
                      already.transformed=FALSE, pch=20, aspect="iso", prepanel=prepanel.qqunif,
                      par.settings=list(superpose.symbol=list(pch=pch)), ...) {


  #error checking
  if (length(pvalues) == 0) stop("pvalue vector is empty, can't draw plot")
  if (!(is.numeric(pvalues) ||
        (is.list(pvalues) && all(vapply(pvalues, is.numeric, logical(1)))))) {
    stop("pvalue vector is not numeric, can't draw plot")
  }
  if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
  if (!already.transformed) {
    if (any(unlist(pvalues) == 0)) stop("pvalue vector contains zeros, can't draw plot")
  } else {
    if (any(unlist(pvalues) < 0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
  }


  grp<-NULL
  n<-1
  exp.x<-c()
  if(is.list(pvalues)) {
    nn<-sapply(pvalues, length)
    rs<-cumsum(nn)
    re<-rs-nn+1
    n<-min(nn)
    if (!is.null(names(pvalues))) {
      grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
      names(pvalues)<-NULL
    } else {
      grp=factor(rep(1:length(pvalues), nn))
    }
    pvo<-pvalues
    pvalues<-numeric(sum(nn))
    exp.x<-numeric(sum(nn))
    for(i in 1:length(pvo)) {
      if (!already.transformed) {
        pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
        exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
      } else {
        pvalues[rs[i]:re[i]] <- pvo[[i]]
        exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
      }
    }
  } else {
    n <- length(pvalues)+1
    if (!already.transformed) {
      exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
      pvalues <- -log10(pvalues)
    } else {
      exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
    }
  }


  #this is a helper function to draw the confidence interval
  panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
    conf.points = min(conf.points, n-1);
    mpts<-matrix(nrow=conf.points*2, ncol=2)
    for(i in seq(from=1, to=conf.points)) {
      mpts[i,1]<- -log10((i-.5)/n)
      mpts[i,2]<- -log10(stats::qbeta(1-conf.alpha/2, i, n-i))
      mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
      mpts[conf.points*2+1-i,2]<- -log10(stats::qbeta(conf.alpha/2, i, n-i))
    }
    grid::grid.polygon(x=mpts[,1],y=mpts[,2], gp=grid::gpar(fill=conf.col, lty=0), default.units="native")
  }

  #reduce number of points to plot
  if (should.thin) {
    if (!is.null(grp)) {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places),
                                grp=grp))
      grp = thin$grp
    } else {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places)))
    }
    pvalues <- thin$pvalues
    exp.x <- thin$exp.x
  }
  gc()

  prepanel.qqunif= function(x,y,...) {
    A = list()
    A$xlim = range(x, y)*1.02
    A$xlim[1]=0
    A$ylim = A$xlim
    return(A)
  }

  #draw the plot
  lattice::xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
         prepanel=prepanel, scales=list(axs="i"), pch=pch,
         panel = function(x, y, ...) {
           if (draw.conf) {
             panel.qqconf(n, conf.points=conf.points,
                          conf.col=conf.col, conf.alpha=conf.alpha)
           };
           lattice::panel.xyplot(x,y, ...);
           lattice::panel.abline(0,1);
         }, par.settings=par.settings, ...
  )
}

