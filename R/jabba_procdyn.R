#' fit_ebsb()
#'
#' Fits an externally computed equilibrium ratio of exploitable to spawning biomass, EB/SB, as a saturating function of
#' spawning-biomass depletion P = SB/SB0, for use in \code{build_jabba(index.ebsb = )} (Selex module):
#' EB/SB(P) = u1 + (u2 - u1) * (1 - exp(-u3 * (P - P1))) / (1 - exp(-u3 * (P2 - P1)))
#'
#' Same curve as \code{make_ebsb()} (which computes it from \code{jbio}/\code{jselex} life history), but for an
#' equilibrium curve computed elsewhere, e.g. from FLBRP or from length-based selectivity estimated with FLicc. With FLBRP:
#' \preformatted{
#' fbar(b) <- FLQuant(seq(0, 4*Fmsy, length.out = 120))
#' sel <- catch.sel(b) %/% apply(catch.sel(b), 2:6, max)
#' sb <- c(ssb(b)); eb <- c(quantSums(stock.n(b) %*% stock.wt(b) %*% sel))
#' fit_ebsb(P = sb/sb[1], ebsb = eb/sb)
#' }
#' @param P spawning biomass depletion SB/SB0 along the equilibrium curve
#' @param ebsb ratio EB/SB at each P
#' @param pmin smallest depletion used in the fit (default 0.005)
#' @param plot logical, plot the fit
#' @return list(pars = c(u1, u2, u3, P1, P2), fit = nls object); pass \code{$pars} to \code{index.ebsb}
#' @export
fit_ebsb <- function(P, ebsb, pmin = 0.005, plot = FALSE){
  d = data.frame(P = as.numeric(P), EB_SB = as.numeric(ebsb))
  d = d[is.finite(d$P) & is.finite(d$EB_SB) & d$P > pmin, ]
  if(nrow(d) < 5) stop("fit_ebsb needs at least 5 points along the equilibrium curve")
  P1 = min(d$P); P2 = max(d$P)
  fit = stats::nls(EB_SB ~ u1 + (u2 - u1) * (1 - exp(-u3 * (P - P1))) / (1 - exp(-u3 * (P2 - P1))), data = d,
                   start = list(u1 = d$EB_SB[which.min(d$P)], u2 = d$EB_SB[which.max(d$P)], u3 = 3),
                   control = stats::nls.control(maxiter = 200, warnOnly = TRUE))
  pars = c(stats::coef(fit), P1 = P1, P2 = P2)
  if(plot){
    graphics::plot(d$P, d$EB_SB, xlab = "SB/SB0", ylab = "EB/SB", pch = 16, col = "grey")
    graphics::lines(d$P, ebsb_ratio(d$P, pars), lwd = 2, col = 4)
  }
  list(pars = pars, fit = fit)
}

#' ebsb_ratio()
#'
#' Evaluates the fitted EB/SB curve from \code{make_ebsb()} or \code{fit_ebsb()} at spawning biomass depletion P
#' @param P SB/SB0
#' @param pars c(u1, u2, u3, P1, P2)
#' @return EB/SB at P
#' @export
ebsb_ratio <- function(P, pars){
  pars = as.numeric(pars)
  Pc = pmin(pmax(P, pars[4]), pars[5])
  pars[1] + (pars[2] - pars[1]) * (1 - exp(-pars[3] * (Pc - pars[4]))) / (1 - exp(-pars[3] * (pars[5] - pars[4])))
}

#' ebsb_convert()
#'
#' Post-hoc conversion of exploitable-biomass status from a standard JABBA fit to spawning-biomass status:
#' SB/SBmsy = (EB/EBmsy) * r(Pmsy) / r(P), with r = EB/SB and P = SB/SB0, solved by fixed-point iteration.
#' Recommended when selectivity is logistic and fish are selected before maturity; for selection at or after
#' maturity or dome-shaped selectivity, fit the Selex module (\code{index.ebsb}) instead.
#' @param BBmsy vector of EB/EBmsy (e.g. \code{fit$timeseries[,"mu","BBmsy"]})
#' @param ebsb c(u1, u2, u3, P1, P2) from \code{make_ebsb()} or \code{fit_ebsb()}
#' @param SBmsyK SBmsy/SB0 of the equilibrium curve
#' @param iter number of fixed-point iterations
#' @return SB/SBmsy
#' @export
ebsb_convert <- function(BBmsy, ebsb, SBmsyK, iter = 50){
  rmsy = ebsb_ratio(SBmsyK, ebsb)
  x = BBmsy
  for(i in seq_len(iter)) x = BBmsy * rmsy / ebsb_ratio(SBmsyK * x, ebsb)
  x
}

#' jabba_tier2()
#'
#' Arguments for Tier 2 of the BioProcDyn protocol: a two-step tuned, depletion-dependent process error
#' with a production lag. Step 1 is a JABBA fit with constant, freely estimated process error (current practice
#' or the Selex module without correction). Step 2 fixes the level so that the dynamic form
#' sd(x) = sig.msy * x^-k reproduces the step-1 variance over the step-1 trajectory,
#' sig.msy = sigma1 / sqrt(mean(x1^(-2*k0))), and estimates k with a prior centred on a life-history value
#' k0 = max(0, 0.064 - 0.155 * log(r)), derived from age-structured simulations of 150 FishLife species.
#' This is analogous to the iterative tuning of recruitment variability in Stock Synthesis.
#'
#' Usage: \code{jb2 = do.call(build_jabba, c(args_step1, jabba_tier2(fit1, r = 0.4, amat = 3)))}
#' @param step1 fitted JABBA model (output of \code{fit_jabba}) with constant process error
#' @param r intrinsic rate of increase (e.g. the mean of the r prior); sets k0 if k0 is NULL
#' @param amat age at 50\% maturity (years) for the production lag (Tier 1); NULL = no lag
#' @param k0 prior mean of k; NULL = life-history value from r; 0 for short-lived species
#' @param k.sd prior sd of k (default 0.15)
#' @param k.bounds bounds of k (default c(-0.5, 2))
#' @param xlim B/Bmsy range over which the dynamic form is evaluated (default c(0.2, 3))
#' @param biascor logical, mean-unbiased lognormal process error (\code{proc.biascor})
#' @return list of build_jabba() arguments: proc.floor, proc.k.prior, proc.xlim, sigma.proc, proc.lag, proc.biascor
#' @export
jabba_tier2 <- function(step1, r = NULL, amat = NULL, k0 = NULL, k.sd = 0.15, k.bounds = c(-0.5, 2),
                        xlim = c(0.2, 3), biascor = FALSE){
  if(is.null(k0)){
    if(is.null(r)) stop("provide r or k0")
    k0 = max(0, 0.064 - 0.155 * log(r))
  }
  sig1 = sqrt(step1$pars["sigma2", "Median"])
  x1 = pmin(pmax(step1$timeseries[, "mu", "BBmsy"], xlim[1]), xlim[2])
  sig.msy = sig1 / sqrt(mean(x1^(-2 * k0)))
  out = list(proc.floor = c(sig.msy, k0), proc.k.prior = c(k0, k.sd, k.bounds), proc.xlim = xlim,
             sigma.proc = 0.001, proc.biascor = biascor)
  if(!is.null(amat)) out$proc.lag = max(0, round(amat))
  attr(out, "tuning") = c(sigma.step1 = sig1, sig.msy = sig.msy, k0 = k0)
  out
}

# internal: posterior draws of spawning biomass and depletion from kbtrj
.selex_draws <- function(jabba){
  kb = jabba$kbtrj
  if(is.null(kb)) stop("requires fit_jabba(..., save.trj = TRUE)")
  kb = kb[kb$type=="fit" & as.character(kb$run)==as.character(jabba$scenario),]
  kb$PK = kb$BB0/kb$stock   # SBmsy/K per draw
  kb
}
.qband <- function(x, years, col, lwd=2, CI=TRUE, lty=1){
  q = suppressWarnings(apply(x, 2, quantile, c(0.5,0.025,0.975), na.rm=TRUE))
  q[, colSums(!is.na(x))==0] = NA
  # draw each run of years with data separately (gaps = years without catch for that fleet)
  ok = !is.na(q[1,]); r = rle(ok); end = cumsum(r$lengths); start = end - r$lengths + 1
  for(k in which(r$values)){ i = start[k]:end[k]
    if(CI) polygon(c(years[i],rev(years[i])), c(q[2,i],rev(q[3,i])), col=grDevices::adjustcolor(col,0.2), border=NA)
    lines(years[i], q[1,i], col=col, lwd=lwd, lty=lty) }
  invisible(q)
}

#' jbplot_fleetF()
#'
#' Harvest rate by fleet. Each fleet's harvest rate acts on its own exploitable biomass,
#' h_f = C_f / EB_f with EB_f = SB * (EB/SB)_f(SB/SB0) (Selex) or EB_f = B (standard JABBA, single EB state).
#' Panel 1: h_f by fleet (median, 95\% CI). Panel 2: total F/Fmsy with the contribution of each fleet,
#' apportioned by its share of the removals from the state (C_f/(EB/SB)_f).
#' @param jabba output of fit_jabba (save.trj = TRUE)
#' @param ebsb.fleets optional named list of EB/SB curves c(u1,u2,u3,P1,P2) by fleet (names = catch columns);
#'        default: the curves used in the fit, by catch column (Selex), or none (standard JABBA)
#' @param catch optional data.frame(year, fleet1, fleet2, ...); default jabba$inputseries$catch
#' @param cols fleet colours
#' @param plotCIs show 95\% credibility intervals
#' @param add if TRUE, do not reset par
#' @param fmax upper y limit for F/Fmsy
#' @return invisible list with harvest-rate quantiles by fleet and the F/Fmsy shares
#' @export
jbplot_fleetF <- function(jabba, ebsb.fleets=NULL, catch=NULL, cols=NULL, plotCIs=TRUE, add=FALSE, fmax=NULL){
  if(is.null(catch)) catch = jabba$inputseries$catch
  catch = as.data.frame(catch); fleets = names(catch)[-1]; nf = length(fleets)
  kb = .selex_draws(jabba); years = sort(unique(kb$year))
  catch = catch[match(years, catch[,1]), , drop=FALSE]
  catch[!is.na(catch) & catch <= 0] = NA          # zero or missing catch: fleet not fishing that year
  selex = !is.null(jabba$selex)
  if(is.null(ebsb.fleets) & selex){   # curves used in the fit, by catch column
    cm = ebsb_matrix(jabba$selex$ebsb); map = jabba$selex$setup$catch
    ebsb.fleets = setNames(lapply(fleets, function(f) cm[ifelse(is.null(map) || is.na(map[f]), 1, map[f]),]), fleets)
  }
  if(!is.null(ebsb.fleets) && is.null(names(ebsb.fleets))) names(ebsb.fleets) = fleets
  if(is.null(cols)) cols = grDevices::hcl.colors(max(nf,2), "Dark 3")[1:nf]
  B = t(sapply(split(kb$B, kb$iter), identity)); P = t(sapply(split(kb$BB0, kb$iter), identity))
  Hf = Rf = list()
  for(f in fleets){
    rho = if(!is.null(ebsb.fleets)) matrix(ebsb_ratio(P, ebsb.fleets[[f]]), nrow(P)) else 1
    Cf = matrix(catch[[f]], nrow(B), length(years), byrow=TRUE)
    Hf[[f]] = Cf/(B*rho); Rf[[f]] = ifelse(is.na(Cf), 0, Cf/rho)
  }
  Rtot = Reduce(`+`, Rf); share = lapply(Rf, function(x) apply(x/Rtot, 2, median, na.rm=TRUE))
  FF = t(sapply(split(kb$harvest, kb$iter), identity)); FFq = apply(FF, 2, quantile, c(0.5,0.025,0.975))
  if(!add) par(mfrow=c(1,2), mar=c(3.5,3.5,1,0.5), mgp=c(2,0.5,0), tck=-0.02, cex=0.8)
  ymax = max(sapply(Hf, function(x) quantile(x, ifelse(plotCIs,0.975,0.5), na.rm=TRUE)), na.rm=TRUE)
  plot(years, years, type="n", ylim=c(0, ymax*1.05), xlab="Year", ylab=ifelse(selex,"Harvest rate on EB by fleet","Harvest rate by fleet"))
  Hq = Map(function(x, cl) .qband(x, years, cl, CI=plotCIs), Hf, cols)
  legend("topleft", fleets, col=cols, lwd=2, bty="n", cex=0.8)
  ymax2 = if(is.null(fmax)) max(FFq[ifelse(plotCIs,3,1),], 1.2, na.rm=TRUE)*1.05 else fmax
  plot(years, years, type="n", ylim=c(0, ymax2), xlab="Year", ylab=expression(F/F[MSY]))
  lo = rep(0, length(years))
  for(j in seq_along(fleets)){ hi = lo + FFq[1,]*share[[j]]
    polygon(c(years,rev(years)), c(lo,rev(hi)), col=grDevices::adjustcolor(cols[j],0.6), border=NA); lo = hi }
  if(plotCIs){ lines(years, FFq[2,], lty=2); lines(years, FFq[3,], lty=2) }
  lines(years, FFq[1,], lwd=2); abline(h=1, lty=5)
  invisible(list(H=Hq, FFmsy=FFq, share=share))
}

#' jbplot_ebsb_status()
#'
#' Spawning-biomass and exploitable-biomass status from a Selex fit (index.ebsb): SB/SBmsy (state) and
#' EB_f/EB_f,msy = (SB/SBmsy) * (EB/SB)_f(P) / (EB/SB)_f(Pmsy) for the model curve and, optionally, by fleet.
#' Divergence between the lines shows where status read from exploitable biomass misrepresents spawning-stock status.
#' @param jabba output of fit_jabba with index.ebsb (save.trj = TRUE)
#' @param ebsb.fleets optional named list of EB/SB curves; default: the curves used in the fit
#' @param cols colours for fleets
#' @param plotCIs show 95\% credibility intervals
#' @param add if TRUE, do not reset par
#' @return invisible list of quantiles
#' @export
jbplot_ebsb_status <- function(jabba, ebsb.fleets=NULL, cols=NULL, plotCIs=TRUE, add=FALSE){
  if(is.null(jabba$selex)) stop("jbplot_ebsb_status requires a Selex fit (build_jabba(index.ebsb = ...))")
  kb = .selex_draws(jabba); years = sort(unique(kb$year))
  SBB = t(sapply(split(kb$stock, kb$iter), identity)); P = t(sapply(split(kb$BB0, kb$iter), identity))
  Pm = t(sapply(split(kb$PK, kb$iter), identity))
  if(is.null(ebsb.fleets)){ cm = ebsb_matrix(jabba$selex$ebsb); ebsb.fleets = setNames(lapply(seq_len(nrow(cm)), function(i) cm[i,]), rownames(cm)) }
  curves = ebsb.fleets
  if(is.null(cols)) cols = grDevices::hcl.colors(max(length(curves),2), "Dark 3")[seq_along(curves)]
  EBB = lapply(curves, function(u) SBB*matrix(ebsb_ratio(P,u),nrow(P))/matrix(ebsb_ratio(Pm,u),nrow(P)))
  if(!add) par(mfrow=c(1,1), mar=c(3.5,3.5,1,0.5), mgp=c(2,0.5,0), tck=-0.02, cex=0.8)
  ymax = max(apply(SBB,2,quantile,ifelse(plotCIs,0.975,0.5)), sapply(EBB, function(x) max(apply(x,2,median))), 1.2)*1.05
  plot(years, years, type="n", ylim=c(0,ymax), xlab="Year", ylab="Status relative to MSY level")
  out = list(SBB = .qband(SBB, years, "black", CI=plotCIs))
  for(j in seq_along(curves)) out[[names(curves)[j]]] = .qband(EBB[[j]], years, cols[j], CI=FALSE)
  abline(h=1, lty=5)
  legend("topright", c(expression(SB/SB[MSY]), paste0("EB/EBmsy: ", names(curves))), col=c("black",cols), lwd=2,
         bty="n", cex=0.8)
  invisible(out)
}

#' ebsb_matrix()
#'
#' Standardises EB/SB curve input for \code{build_jabba(index.ebsb = )}: a single vector c(u1,u2,u3,P1,P2),
#' a named list of such vectors, a matrix, or the \code{$pars} data.frame of \code{make_ebsb()}
#' @param x curve(s)
#' @return numeric matrix (curves x 5) with columns u1, u2, u3, P1, P2 and curve names as row names
#' @export
ebsb_matrix <- function(x){
  nm = c("u1","u2","u3","P1","P2")
  if(is.data.frame(x)){
    rn = if("selectivity" %in% names(x)) as.character(x$selectivity) else rownames(x)
    m = as.matrix(x[, nm, drop=FALSE]); rownames(m) = rn
  } else if(is.list(x)){
    m = do.call(rbind, lapply(x, function(z) as.numeric(unlist(z))[1:5]))
    rownames(m) = if(is.null(names(x))) paste0("sel", seq_along(x)) else names(x)
  } else if(is.matrix(x)){
    m = x[, 1:5, drop=FALSE]; if(is.null(rownames(m))) rownames(m) = paste0("sel", seq_len(nrow(m)))
  } else {
    m = matrix(as.numeric(x)[1:5], nrow=1, dimnames=list("sel1", NULL))
  }
  storage.mode(m) = "double"; colnames(m) = nm
  if(any(!is.finite(m))) stop("index.ebsb must contain finite c(u1,u2,u3,P1,P2) for each curve")
  m
}

#' ebsb_eff()
#'
#' Effective EB/SB of combined removals: 1 / sum_c w_c / (EB/SB)_c(P), with w the catch allocation among curves
#' @param P SB/SB0 (vector or matrix)
#' @param curves curve matrix from \code{ebsb_matrix()}
#' @param w allocation among curves (sums to 1)
#' @return effective EB/SB with the dimensions of P
#' @export
ebsb_eff <- function(P, curves, w){
  curves = ebsb_matrix(curves); w = rep(w, length.out=nrow(curves)); w = w/sum(w)
  out = 0*P
  for(c in seq_len(nrow(curves))) out = out + w[c]/ebsb_ratio(P, curves[c,])
  1/out
}
