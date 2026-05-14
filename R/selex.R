
#' jgrow()
#'
#' Computes length-at-age from a growth model
#' @param param list with growth parameters:
#' \itemize{
#'   \item \code{Linf}: asymptotic length
#'   \item \code{k}: growth coefficient
#'   \item \code{t0}: theoretical age at zero length
#'   \item \code{L0}: length at age 0, optional alternative to \code{t0}
#' }
#' @param age vector of ages
#' @param nsexes number of sexes
#' @param model growth model c("vb","schnute")
#' @return list with age, nsexes, La, model and param
#' @export
#' @examples
#' pars <- list(
#'  Linf   = c(80,65),
#'   k     = c(0.2,0.25),
#'   t0     = c(-0.5,-0.7)
#' )
#' g <- jgrow(pars, model = "vb")
#' head(g$La)
#' 
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
#'                 Linf = c(80,65), K = c(0.2,0.25), t0 = c(-0.5,-0.7)))
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

