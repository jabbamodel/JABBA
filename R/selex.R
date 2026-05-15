#' Compute length-at-age from a growth model
#'
#' Computes expected length-at-age for one or more sexes using a specified
#' growth model. The current implementation supports the von Bertalanffy
#' growth model, with optional support for specifying length at age zero
#' through \code{L0} instead of \code{t0}.
#'
#' @param param A list of growth parameters. For \code{model = "vb"}, this should
#'   contain:
#'   \itemize{
#'     \item \code{Linf}: asymptotic length.
#'     \item \code{k}: von Bertalanffy growth coefficient.
#'     \item \code{t0}: theoretical age at zero length.
#'     \item \code{L0}: optional length at age zero, used as an alternative to
#'       \code{t0} when provided.
#'   }
#'   Parameters may be scalar or sex-specific vectors of length \code{nsexes}.
#' @param age Numeric vector of ages. Default is \code{0:20}.
#' @param nsexes Integer. Number of sexes. Default is \code{1}.
#' @param model Character. Growth model to use. Currently \code{"vb"} is
#'   implemented; \code{"schnute"} is reserved for future use.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{age}}{Vector of ages.}
#'   \item{\code{nsexes}}{Number of sexes.}
#'   \item{\code{La}}{Matrix of expected length-at-age with dimensions
#'   age by sex.}
#'   \item{\code{model}}{Growth model used.}
#'   \item{\code{param}}{Input growth parameters.}
#' }
#'
#' @examples
#' pars <- list(
#'   Linf = c(80, 65),
#'   k    = c(0.2, 0.25),
#'   t0   = c(-0.5, -0.7)
#' )
#'
#' g <- jgrow(
#'   param = pars,
#'   age = 0:20,
#'   nsexes = 2,
#'   model = "vb"
#' )
#'
#' head(g$La)
#'
#' @export
jgrow <- function(param=list(Linf=80,k=0.2,t0=-0.5),age=0:12,nsexes=NULL, model = "vb"){
  
  model <- match.arg(model, c("vb", "schnute"))
  
  if (!is.list(param)) stop("'param' must be a named list")
  
  nsexes <- if(is.null(nsexes)) length(param[[1]]) else nsexes
  
  if (is.null(age)) stop("age must be provided")
  
  La <- matrix(NA_real_, nrow = length(age), ncol = 2)
  
  if (model == "vb") {
    
    Linf <- param$Linf
    k    <- param$k
    t0   <- param$t0
    L0   <- param$L0
    
    if (is.null(Linf) || is.null(k)) {
      stop("For model = 'vb', param$Linf and param$K must be provided")
    }
    
    vb_fun <- function(age, Linf, k, t0 = NULL, L0 = NULL) {
      if (is.null(t0) && !is.null(L0)) {
        t0 <- (1 / k) * log((Linf - L0) / Linf)
      } else if (is.null(L0) && !is.null(t0)) {
        L0 <- Linf * (1 - exp(-k * (0 - t0)))
      } else if (is.null(t0) && is.null(L0)) {
        stop("Either 't0' or 'L0' must be supplied for VB growth")
      }
      
      Linf * (1 - exp(-k * (age - t0)))
    }
    
    La[, 1] <- vb_fun(
      age  = age,
      Linf = Linf[1],
      k    = k[1],
      t0   = if (is.null(t0)) NULL else t0[1],
      L0   = if (is.null(L0)) NULL else L0[1]
    )
    
    if (nsexes > 1) {
      La[, 2] <- vb_fun(
        age  = age,
        Linf = Linf[length(Linf)],
        k    = k[length(k)],
        t0   = if (is.null(t0)) NULL else t0[length(t0)],
        L0   = if (is.null(L0)) NULL else L0[length(L0)]
      )
    } else {
      La[, 2] <- La[, 1]
    }
  }
  
  if (model == "schnute") {
    stop("model = 'schnute' is not yet implemented")
  }
  
  La[La < 0.01] <- 0.01
  
  out <- list(
    model  = model,
    age    = age,
    nsexes = nsexes,
    La     = La,
    param  = param
  )
  
  return(out)
}


#' jbio()
#'
#' Creates a stock object for JABBA-Select life history
#' @param amin minimum age
#' @param amax maximum age
#' @param nsexes number of sexes
#' @param grow optional output list from \code{jgrow()}
#' @param model growth model c("vb","schnute") if \code{grow=NULL}
#' @param Linf asymptotic length
#' @param k growth coefficient
#' @param t0 theoretical age at zero length
#' @param L0 optional alternative to \code{t0}
#' @param aW length-weight parameter a
#' @param bW length-weight parameter b
#' @param mat maturity parameters c(m50,m95,type), where type = 0 for length and 1 for age
#' @param M natural mortality
#' @param h steepness
#' @param fec optional fecundity at age
#' @param empty logical option to create empty object
#' @return list with age, La, Wa, Mat, M, h and optional fec
#' @export
#' @examples
#' stk <- jbio(
#'   amin = 0, amax = 20, nsexes = 2,
#'   Linf = c(80,65), k = c(0.2,0.25), t0 = c(-0.5,-0.7),
#'   aW = 0.01, bW = 3.04,
#'   mat = c(40,45,0), M = c(0.3,0.4), h = 0.7
#' )
#' head(stk$La)
#'
#' g <- jgrow(list(age = 0:20, nsexes = 2,
#'                 Linf = c(80,65), k = c(0.2,0.25), t0 = c(-0.5,-0.7)))
#' stk2 <- jbio(grow = g, aW = 0.01, bW = 3.04,
#'              mat = c(40,45,0), M = c(0.3,0.4), h = 0.7)
jbio <- function(amin = 0, amax = 10, nsexes = 2,
                 grow = NULL, model = "vb",
                 Linf = c(80,65), k = c(0.2,0.25), t0 = c(-0.5,-0.7), L0 = NULL,
                 aW = 0.01, bW = 3.04, mat = c(40,45,0),
                 M = c(0.3,0.4), h = 0.7, fec = NA, empty = FALSE){
  
  if(empty){
    stk <- list()
    stk$nsexes <- nsexes
    stk$age <- amin:amax
    stk$h <- h
    return(stk)
  }
  
  if(is.null(grow)){
    age <- amin:amax
    grow <- jgrow(
      param = list(age = age, nsexes = nsexes, Linf = Linf, k = k, t0 = t0, L0 = L0),
      model = model
    )
  }
  
  age <- grow$age
  nages <- length(age)
  
  stk <- list()
  stk$nsexes <- nsexes
  stk$nseas  <- 1
  stk$age    <- age
  stk$La     <- grow$La
  stk$Wa     <- matrix(NA, nrow = nages, ncol = 2)
  stk$Mat    <- matrix(NA, nrow = nages, ncol = 2)
  stk$M      <- matrix(NA, nrow = 1, ncol = 2)
  stk$h      <- h
  stk$grow   <- grow
  
  if(length(aW) == 1) aW <- rep(aW, 2)
  if(length(bW) == 1) bW <- rep(bW, 2)
  if(length(M)  == 1) M  <- rep(M, 2)
  
  stk$Wa[,1] <- aW[1] * stk$La[,1]^bW[1]
  stk$Wa[,2] <- aW[length(aW)] * stk$La[,2]^bW[length(bW)]
  
  m50  <- mat[1]
  m95  <- mat[2]
  mtyp <- mat[3]
  
  if(mtyp == 0){
    stk$Mat[,1] <- 1 / (1 + exp(-log(19) * (stk$La[,1] - m50) / (m95 - m50)))
    stk$Mat[,2] <- 1 / (1 + exp(-log(19) * (stk$La[,2] - m50) / (m95 - m50)))
  }
  
  if(mtyp == 1){
    stk$Mat[,1] <- 1 / (1 + exp(-log(19) * (age - m50) / (m95 - m50)))
    stk$Mat[,2] <- stk$Mat[,1]
  }
  
  stk$M[1,] <- M
  
  if(!all(is.na(fec))){
    if(length(fec) == nages){
      stk$fec <- cbind(fec, fec)
    } else if(is.matrix(fec)) {
      stk$fec <- fec
    } else {
      warning("fec ignored because dimensions do not match age structure")
    }
  }
  
  return(stk)
}


#' jselex()
#'
#' Computes selectivity-at-age from growth
#' @param grow growth object from \code{jgrow()}
#' @param selpars vector with selectivity parameters:
#' \itemize{
#'   \item \code{SL50}: length at 50\% selection
#'   \item \code{SL95}: length at 95\% selection
#'   \item \code{SLdesc}: mean of descending limb
#'   \item \code{CVdesc}: CV of descending limb
#'   \item \code{Smin}: minimum selectivity on descending limb
#' }
#' @return list with Sel and selpars
#' @export
#' @examples
#' g <- jgrow()
#'   
#' sel <- jselex(g, selpars = c(35,45,70,0.15,0.05))
#' head(sel$Sel)
jselex <- function(grow, selpars){
  
  if(is.null(grow$La))
    stop("grow$La is missing")
  
  if(missing(selpars) || is.null(selpars))
    stop("selpars must be provided")
  
  if(length(selpars) < 5)
    stop("selpars must be c(SL50, SL95, SLdesc, CVdesc, Smin)")
  
  La     <- grow$La
  nages  <- nrow(La)
  nsexes <- if(is.null(grow$nsexes)) 1 else grow$nsexes
  
  Sel <- matrix(NA, nrow = nages, ncol = 2)
  
  SL50   <- selpars[1]
  SL95   <- selpars[2]
  SLdesc <- selpars[3]
  CVdesc <- selpars[4]
  Smin   <- selpars[5]
  
  get_sel <- function(Li, SL50, SL95, SLdesc, CVdesc, Smin){
    
    sel_up <- 1 / (1 + exp(-log(19) * (Li - SL50) / (SL95 - SL50)))
    
    sel_dn <- dnorm(Li, mean = SLdesc, sd = SLdesc * CVdesc)
    sel_dn <- sel_dn / max(sel_dn, na.rm = TRUE)
    
    sel_dn <- 1 + (Smin - 1) * sel_dn
    
    sel <- ifelse(Li <= SLdesc, sel_up, sel_dn)
    sel[sel < 0] <- 0
    sel[sel > 1] <- 1
    as.numeric(sel)
  }
  
  Sel[,1] <- get_sel(La[,1], SL50, SL95, SLdesc, CVdesc, Smin)
  
  if(nsexes > 1){
    Sel[,2] <- get_sel(La[,2], SL50, SL95, SLdesc, CVdesc, Smin)
  } else {
    Sel[,2] <- Sel[,1]
  }
  
  return(list(Sel = Sel, selpars = selpars))
}


#' Calculate exploitable-biomass to spawning-biomass ratios
#'
#' Computes the equilibrium relationship between exploitable biomass and
#' spawning biomass, \eqn{EB/SB}, as a function of spawning biomass depletion
#' \eqn{P = SB/SB_0} for one or more selectivity curves. This helper is intended
#' for the JABBA-Select \code{select_mode = "selex"} option, where the
#' observation model is corrected for the fact that abundance indices generally
#' observe exploitable biomass rather than spawning biomass.
#'
#' The calculation uses the life-history information in a stock object created by
#' \code{\link{jbio}} and one or more selectivity objects created by
#' \code{\link{jselex}}. For each fishing mortality value in \code{Fgrid}, the
#' function computes survivorship-at-age, spawning biomass per recruit,
#' exploitable biomass per recruit, equilibrium recruitment relative to unfished
#' recruitment under a Beverton-Holt stock-recruitment relationship, and the
#' resulting depletion \eqn{SB/SB_0}. The resulting \eqn{EB/SB} curve can
#' optionally be approximated with a three-parameter asymptotic curve for use as
#' fixed data in the JAGS observation model.
#'
#' @param stk A life-history stock object returned by \code{\link{jbio}}.
#'   Required elements are \code{Wa}, \code{Mat}, \code{M}, \code{h},
#'   \code{nsexes}, and either \code{age} or \code{grow$age}.
#' @param selex A single selectivity object returned by \code{\link{jselex}},
#'   or a named list of such objects. If a single object is supplied it is
#'   internally converted to a one-element list. If a list is unnamed, default
#'   names \code{sel1}, \code{sel2}, ... are assigned.
#' @param Fgrid Numeric vector of fishing mortality values used to generate the
#'   equilibrium curve. The default is \code{seq(0, 5, length.out = 100)}.
#' @param sex_ratio Optional numeric vector giving the relative contribution of
#'   each sex. If \code{NULL}, equal sex ratios are assumed.
#' @param plus_group Logical. Should the last age be treated as a plus group?
#'   Default is \code{TRUE}.
#' @param fit_curve Logical. Should the equilibrium \eqn{EB/SB} curve be fitted
#'   with the three-parameter asymptotic approximation used for the JABBA-Select
#'   observation correction? Default is \code{TRUE}.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{by_selectivity}}{A named list with one element per selectivity
#'   object. Each element contains the raw curve, the fitted curve object, and
#'   fitted parameters.}
#'   \item{\code{curve}}{A combined data frame with columns \code{selectivity},
#'   \code{F}, \code{SBPR}, \code{EBPR}, \code{Rrel}, \code{P}, \code{EB_SB},
#'   and, when \code{fit_curve = TRUE}, \code{EB_SB_fit}.}
#'   \item{\code{pars}}{A data frame of fitted approximation parameters:
#'   \code{selectivity}, \code{u1}, \code{u2}, \code{u3}, \code{P1}, and
#'   \code{P2}. These can be passed as fixed data to the JABBA observation
#'   model.}
#'   \item{\code{selex}}{The selectivity object or list of selectivity objects
#'   used in the calculation.}
#'   \item{\code{settings}}{A list of settings used for the calculation.}
#' }
#'
#' @details
#' The fitted approximation has the form
#' \deqn{
#' EB/SB(P) =
#' u_1 + (u_2 - u_1)
#' \frac{1 - \exp[-u_3(P - P_1)]}
#'      {1 - \exp[-u_3(P_2 - P_1)]}
#' }
#' where \eqn{P = SB/SB_0}, \eqn{P_1} and \eqn{P_2} are the minimum and maximum
#' depletion values in the generated equilibrium curve, and \eqn{u_1},
#' \eqn{u_2}, and \eqn{u_3} are estimated by nonlinear least squares.
#'
#' @seealso \code{\link{jgrow}}, \code{\link{jbio}}, \code{\link{jselex}},
#'   \code{\link{plot_ebsb}}
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
#' sel_early <- jselex(stk$grow, selpars = c(35, 45, 999, 0.15, 1))
#' sel_late  <- jselex(stk$grow, selpars = c(50, 60, 999, 0.15, 1))
#'
#' ebsb <- make_ebsb(
#'   stk,
#'   selex = list(early = sel_early, late = sel_late)
#' )
#'
#' head(ebsb$curve)
#' ebsb$pars
#'
#' @export
make_ebsb <- function(stk, selex,
                      Fgrid = seq(0, 5, length.out = 100),
                      sex_ratio = NULL,
                      plus_group = TRUE,
                      fit_curve = TRUE) {
  
  # ------------------------------------------------------------
  # Allow either a single jselex object or a named list of jselex objects
  # ------------------------------------------------------------
  
  is_single_selex <- is.list(selex) && !is.null(selex$Sel)
  
  if (is_single_selex) {
    selex <- list(selex = selex)
  }
  
  if (is.null(names(selex)) || any(names(selex) == "")) {
    names(selex) <- paste0("sel", seq_along(selex))
  }
  
  # ------------------------------------------------------------
  # Internal single-selectivity calculator
  # ------------------------------------------------------------
  
  make_ebsb_one <- function(stk, sx, sx_name) {
    
    if (is.null(stk$Wa))  stop("stk$Wa is missing")
    if (is.null(stk$Mat)) stop("stk$Mat is missing")
    if (is.null(stk$M))   stop("stk$M is missing")
    if (is.null(stk$h))   stop("stk$h is missing")
    if (is.null(sx$Sel))  stop("selex object is missing $Sel")
    
    age <- stk$age
    if (is.null(age)) age <- stk$grow$age
    
    nages <- length(age)
    nsexes <- stk$nsexes
    
    if (is.null(sex_ratio)) {
      sr <- rep(1 / nsexes, nsexes)
    } else {
      sr <- sex_ratio
    }
    
    if (length(sr) != nsexes) {
      stop("sex_ratio must have length equal to stk$nsexes")
    }
    
    Wa  <- stk$Wa[, seq_len(nsexes), drop = FALSE]
    Mat <- stk$Mat[, seq_len(nsexes), drop = FALSE]
    Sel <- sx$Sel[, seq_len(nsexes), drop = FALSE]
    
    M <- stk$M
    if (is.matrix(M) || is.data.frame(M)) {
      M <- as.numeric(M[1, seq_len(nsexes)])
    } else {
      M <- rep(as.numeric(M), length.out = nsexes)
    }
    
    h <- stk$h
    
    calc_N <- function(F, M_s, Sel_s) {
      Z <- M_s + F * Sel_s
      
      N <- numeric(nages)
      N[1] <- 1
      
      if (nages > 1) {
        for (a in 2:nages) {
          N[a] <- N[a - 1] * exp(-Z[a - 1])
        }
      }
      
      if (plus_group) {
        N[nages] <- N[nages] / max(1e-12, 1 - exp(-Z[nages]))
      }
      
      N
    }
    
    calc_pr <- function(F) {
      SBPR <- 0
      EBPR <- 0
      
      for (s in seq_len(nsexes)) {
        N_s <- calc_N(F = F, M_s = M[s], Sel_s = Sel[, s])
        
        SBPR <- SBPR + sr[s] * sum(N_s * Wa[, s] * Mat[, s])
        EBPR <- EBPR + sr[s] * sum(N_s * Wa[, s] * Sel[, s])
      }
      
      c(SBPR = SBPR, EBPR = EBPR)
    }
    
    pr0 <- calc_pr(F = 0)
    SBPR0 <- pr0["SBPR"]
    
    out <- data.frame(
      selectivity = sx_name,
      F = Fgrid,
      SBPR = NA_real_,
      EBPR = NA_real_,
      Rrel = NA_real_,
      P = NA_real_,
      EB_SB = NA_real_
    )
    
    for (i in seq_along(Fgrid)) {
      pr <- calc_pr(Fgrid[i])
      
      SBPR <- pr["SBPR"]
      EBPR <- pr["EBPR"]
      
      # Beverton-Holt equilibrium recruitment relative to R0
      Rrel <- (4 * h * SBPR - (1 - h) * SBPR0) /
        ((5 * h - 1) * SBPR)
      
      Rrel <- max(0, Rrel)
      
      P <- (SBPR * Rrel) / SBPR0
      
      out$SBPR[i] <- SBPR
      out$EBPR[i] <- EBPR
      out$Rrel[i] <- Rrel
      out$P[i] <- P
      out$EB_SB[i] <- EBPR / SBPR
    }
    
    out <- out[is.finite(out$P) & is.finite(out$EB_SB) & out$P > 0, ]
    out <- out[order(out$P), ]
    out <- out[!duplicated(round(out$P, 6)), ]
    
    fit <- NULL
    pars <- NULL
    
    if (fit_curve && nrow(out) > 5) {
      
      P1 <- min(out$P)
      P2 <- max(out$P)
      
      start <- list(
        u1 = out$EB_SB[which.min(out$P)],
        u2 = out$EB_SB[which.max(out$P)],
        u3 = 3
      )
      
      nls_fit <- try(
        nls(
          EB_SB ~ u1 + (u2 - u1) *
            (1 - exp(-u3 * (P - P1))) /
            (1 - exp(-u3 * (P2 - P1))),
          data = out,
          start = start,
          control = nls.control(maxiter = 200, warnOnly = TRUE)
        ),
        silent = TRUE
      )
      
      if (!inherits(nls_fit, "try-error")) {
        pars <- coef(nls_fit)
        
        out$EB_SB_fit <-
          pars["u1"] + (pars["u2"] - pars["u1"]) *
          (1 - exp(-pars["u3"] * (out$P - P1))) /
          (1 - exp(-pars["u3"] * (P2 - P1)))
        
        fit <- list(
          pars = pars,
          P1 = P1,
          P2 = P2,
          nls = nls_fit
        )
      } else {
        warning("EB/SB curve fit failed for selectivity: ", sx_name)
        out$EB_SB_fit <- NA_real_
      }
    }
    
    list(
      curve = out,
      fit = fit,
      pars = if (!is.null(fit)) {
        data.frame(
          selectivity = sx_name,
          u1 = unname(fit$pars["u1"]),
          u2 = unname(fit$pars["u2"]),
          u3 = unname(fit$pars["u3"]),
          P1 = fit$P1,
          P2 = fit$P2,
          row.names = NULL
        )
      } else {
        data.frame(
          selectivity = sx_name,
          u1 = NA_real_,
          u2 = NA_real_,
          u3 = NA_real_,
          P1 = NA_real_,
          P2 = NA_real_,
          row.names = NULL
        )
      }
    )
  }
  
  # ------------------------------------------------------------
  # Apply across selectivity objects
  # ------------------------------------------------------------
  
  res <- lapply(seq_along(selex), function(i) {
    make_ebsb_one(stk, selex[[i]], names(selex)[i])
  })
  
  names(res) <- names(selex)
  
  curve <- do.call(rbind, lapply(res, `[[`, "curve"))
  pars  <- do.call(rbind, lapply(res, `[[`, "pars"))
  
  list(
    by_selectivity = res,
    curve = curve,
    pars = pars,
    selex = selex,
    settings = list(
      Fgrid = Fgrid,
      plus_group = plus_group,
      sex_ratio = sex_ratio,
      fit_curve = fit_curve
    )
  )
}



