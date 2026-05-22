
#' Plot exploitable-biomass to spawning-biomass correction curves
#'
#' Plots the relationship between spawning biomass depletion,
#' \eqn{SB/SB_0}, and the exploitable-biomass to spawning-biomass ratio,
#' \eqn{EB/SB}, produced by \code{\link{make_ebsb}}. Points show the raw
#' equilibrium calculations and, optionally, lines show the fitted asymptotic
#' approximation used for the JABBA-Select selectivity correction.
#'
#' @param x An object returned by \code{\link{make_ebsb}}.
#' @param show_fit Logical. If \code{TRUE}, fitted \eqn{EB/SB} approximation
#'   curves are added when available. Default is \code{TRUE}.
#' @param add_maturity_line Logical. If \code{TRUE}, a horizontal reference line
#'   is drawn at \eqn{EB/SB = 1}. Default is \code{TRUE}.
#' @param ylim Optional numeric vector of length two giving y-axis limits. If
#'   \code{NULL}, limits are chosen to include the data and fitted curves, with
#'   at least the range \code{c(0.7, 1.3)} shown.
#' @param lwd Line width for fitted curves and the legend. Default is \code{2}.
#' @param pch Plotting character for raw equilibrium points. Default is
#'   \code{21}.
#' @param lty Line type for fitted curves and legend. Default is \code{1}.
#' @param xlab Label for the x-axis. Default is
#'   \code{expression(SB / SB[0])}.
#' @param ylab Label for the y-axis. Default is \code{"EB/SB"}.
#' @param legend_pos Position of the legend passed to \code{\link[graphics]{legend}}.
#'   Default is \code{"topright"}. Set to \code{NULL} to suppress the legend.
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}.
#'
#' @return Invisibly returns \code{x}.
#'
#' @details
#' Values above one indicate that exploitable biomass is proportionally larger
#' than spawning biomass at a given depletion level. This can occur when the
#' selectivity curve selects fish before or around maturity. Values below one
#' occur when the selected component is smaller than the spawning component,
#' for example when selectivity is shifted toward larger or older fish.
#'
#' @seealso \code{\link{make_ebsb}}, \code{\link{jselex}}
#'
#' @examples
#' g <- jgrow(
#'   list(
#'     Linf = c(80, 65),
#'     k = c(0.2, 0.25),
#'     t0 = c(-0.5, -0.7)
#'   ),
#'   age = 0:20,
#'   nsexes = 2
#' )
#'
#' stk <- jbio(
#'   grow = g,
#'   aW = 0.000001,
#'   bW = 3.04,
#'   mat = c(40, 45, 0),
#'   M = c(0.3, 0.3),
#'   h = 0.7
#' )
#'
#' sel1 <- jselex(stk$grow, selpars = c(35, 45, 999, 0.15, 1))
#' sel2 <- jselex(stk$grow, selpars = c(50, 60, 999, 0.15, 1))
#'
#' ebsb <- make_ebsb(
#'   stk,
#'   selex = list(Sel1 = sel1, Sel2 = sel2)
#' )
#'
#' plot_ebsb(ebsb)
#'
#' @export
plot_ebsb <- function(x,
                      show_fit = TRUE,
                      add_maturity_line = TRUE,
                      ylim = NULL,
                      lwd = 2,
                      pch=21,
                      lty=1,
                      xlab = expression(SB / SB[0]),
                      ylab = "EB/SB",
                      legend_pos = "topright",
                      ...) {
  
  if (is.null(x$curve)) {
    stop("x must be an object returned by make_ebsb()")
  }
  
  dat <- x$curve
  sels <- unique(dat$selectivity)
  
  if (is.null(ylim)) {
    yy <- dat$EB_SB
    
    if (show_fit && "EB_SB_fit" %in% names(dat)) {
      yy <- c(yy, dat$EB_SB_fit)
    }
    
    ylim <- range(c(0.7, 1.3, yy), finite = TRUE)
  }
  
  plot(
    NA, NA,
    xlim = range(dat$P, finite = TRUE),
    ylim = ylim,
    xlab = xlab,
    ylab = ylab,
    ...
  )
  
  if (add_maturity_line) {
    abline(h = 1, lty = 3)
  }
  
  for (i in seq_along(sels)) {
    dd <- dat[dat$selectivity == sels[i], ]
    dd <- dd[order(dd$P), ]
    
    points(dd$P, dd$EB_SB, pch=pch,col = i)
    
    if (show_fit && "EB_SB_fit" %in% names(dd) && any(is.finite(dd$EB_SB_fit))) {
      lines(dd$P, dd$EB_SB_fit, lwd = lwd, lty = lty, col = i)
    }
  }
  
  if (!is.null(legend_pos)) {
    legend(
      legend_pos,
      legend = sels,
      col = seq_along(sels),
      lty = lty,
      lwd = lwd,
      bty = "n"
    )
  }
  
  invisible(x)
}
