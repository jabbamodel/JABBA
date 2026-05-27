
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
    age <- seq(amin, amax, by = 1)
    nages <- length(age)
    ncol_stk <- max(2, nsexes)
    return(list(
      nsexes = nsexes,
      nseas  = 1,
      age    = age,
      La     = matrix(NA_real_, nrow = nages, ncol = ncol_stk),
      Wa     = matrix(NA_real_, nrow = nages, ncol = ncol_stk),
      Mat    = matrix(NA_real_, nrow = nages, ncol = ncol_stk),
      M      = matrix(NA_real_, nrow = nages, ncol = ncol_stk),
      h      = h,
      grow   = NULL,
      fec    = NULL
    ))
  }
  
  if(is.null(grow)){
    age <- seq(amin, amax, by = 1)
    grow <- jgrow(
      param  = list(Linf = Linf, k = k, t0 = t0, L0 = L0),
      age    = age,
      nsexes = nsexes,
      model  = model
    )
  } else {
    if(is.null(grow$age) || is.null(grow$La))
      stop("grow must be an object returned by jgrow(), with $age and $La")
    if(!is.null(grow$nsexes)) nsexes <- grow$nsexes
  }
  
  age <- grow$age
  nages <- length(age)
  ncol_stk <- max(2, nsexes)
  
  rep_to_cols <- function(x, n = ncol_stk) {
    if(length(x) == 1) rep(x, n) else rep(x, length.out = n)
  }
  
  stk <- list()
  stk$nsexes <- nsexes
  stk$nseas  <- 1
  stk$age    <- age
  stk$La     <- grow$La[, seq_len(ncol_stk), drop = FALSE]
  stk$Wa     <- matrix(NA_real_, nrow = nages, ncol = ncol_stk)
  stk$Mat    <- matrix(NA_real_, nrow = nages, ncol = ncol_stk)
  stk$M      <- matrix(NA_real_, nrow = nages, ncol = ncol_stk)
  stk$h      <- h
  stk$grow   <- grow
  
  aW <- rep_to_cols(aW)
  bW <- rep_to_cols(bW)
  M  <- rep_to_cols(M)
  
  for(s in seq_len(ncol_stk)){
    stk$Wa[,s] <- aW[s] * stk$La[,s]^bW[s]
    stk$M[,s]  <- M[s]
  }
  
  m50  <- mat[1]
  m95  <- mat[2]
  mtyp <- mat[3]
  
  if(mtyp == 0){
    for(s in seq_len(ncol_stk)){
      stk$Mat[,s] <- 1 / (1 + exp(-log(19) * (stk$La[,s] - m50) / (m95 - m50)))
    }
  }
  
  if(mtyp == 1){
    mat_age <- 1 / (1 + exp(-log(19) * (age - m50) / (m95 - m50)))
    for(s in seq_len(ncol_stk)) stk$Mat[,s] <- mat_age
  }
  
  stk$fec <- NULL
  if(!all(is.na(fec))){
    if(is.matrix(fec)){
      if(nrow(fec) != nages) stop("fec matrix must have one row per age")
      stk$fec <- fec[, seq_len(ncol_stk), drop = FALSE]
    } else if(length(fec) == 1){
      # scalar fecundity = offspring per mature female; apply maturity ogive
      stk$fec <- stk$Mat * fec
    } else if(length(fec) == nages){
      stk$fec <- matrix(rep(fec, ncol_stk), nrow = nages, ncol = ncol_stk)
    } else {
      warning("fec ignored because dimensions do not match age structure")
    }
  }
  
  return(stk)
}

#' jbspr0()
#'
#' Robust spawning output per recruit at unfished conditions.
#'
#' Computes unfished spawning output per recruit, commonly denoted
#' \eqn{SPR_0}, from survivorship under natural mortality, maturity-at-age
#' and weight-at-age:
#'
#' \deqn{
#'   SPR_0 = \sum_x l_x Mat_x W_x
#' }
#'
#' where \eqn{l_x} is survivorship from recruitment to age \eqn{x},
#' \eqn{Mat_x} is maturity-at-age, and \eqn{W_x} is weight-at-age.
#'
#' This is a robust version of the original \code{jbspr0()} that accepts both
#' older \code{jbio()} objects with \code{M} and \code{Wa} stored as age-by-sex
#' matrices and newer objects where these quantities may be vectors or
#' one-row matrices expanded internally.
#'
#' If \code{plus_group = TRUE}, the final age is treated as a plus group using:
#'
#' \deqn{
#'   l_A^+ = \frac{l_A}{1 - \exp(-M_A)}
#' }
#'
#' This matches the original implementation:
#'
#' \preformatted{
#' survivors[A] <- survivors[A] * (-1 / (exp(-M[A]) - 1))
#' }
#'
#' @param bio list returned by \code{jbio()}. Must contain \code{age},
#'   \code{M}, \code{Mat}, and \code{Wa}.
#' @param byage logical. If \code{TRUE}, return age-specific contributions
#'   \eqn{l_x Mat_x W_x}. If \code{FALSE}, return their sum. Default is
#'   \code{FALSE}.
#' @param female_col integer. Column to use when \code{M}, \code{Mat}, or
#'   \code{Wa} are age-by-sex matrices. Default is \code{1}.
#' @param plus_group logical. If \code{TRUE}, apply the terminal plus-group
#'   correction. Default is \code{TRUE}.
#'
#' @return If \code{byage = FALSE}, a single numeric value giving total
#'   unfished spawning output per recruit. If \code{byage = TRUE}, a numeric
#'   vector of age-specific contributions.
#'
#' @examples
#' # One-sex example
#' bio <- jbio(
#'   amin = 0,
#'   amax = 20,
#'   nsexes = 1,
#'   Linf = 80,
#'   k = 0.2,
#'   t0 = -0.5,
#'   aW = 0.00001,
#'   bW = 3.0,
#'   mat = c(40, 45, 0),
#'   M = 0.3,
#'   h = 0.7
#' )
#'
#' # Total SPR0
#' jbspr0(bio)
#'
#' # SPR0 contribution by age
#' jbspr0(bio, byage = TRUE)
#'
#' # Without terminal plus group
#' jbspr0(bio, plus_group = FALSE)
#'
#' # Sex-specific biology, using the female column
#' bio_2sex <- jbio(
#'   amin = 0,
#'   amax = 20,
#'   nsexes = 2,
#'   Linf = c(90, 75),
#'   k = c(0.18, 0.22),
#'   t0 = c(-0.5, -0.5),
#'   aW = c(0.00001, 0.00001),
#'   bW = c(3.0, 3.0),
#'   mat = c(45, 55, 0),
#'   M = c(0.25, 0.32),
#'   h = 0.7
#' )
#'
#' jbspr0(bio_2sex, female_col = 1)
#'
#' @export
jbspr0 <- function(
    bio = jbio(),
    byage = FALSE,
    female_col = 1,
    plus_group = TRUE
){
  
  get_age_col <- function(x, nages, col = 1, name = "object"){
    
    if(is.null(x)){
      stop(name, " is missing.")
    }
    
    if(is.matrix(x) || is.data.frame(x)){
      
      x <- as.matrix(x)
      
      if(ncol(x) < col){
        stop(name, " has fewer columns than female_col.")
      }
      
      if(nrow(x) == 1){
        return(rep(as.numeric(x[1, col]), nages))
      }
      
      if(nrow(x) != nages){
        stop(name, " has ", nrow(x), " rows but expected ", nages, ".")
      }
      
      return(as.numeric(x[, col]))
    }
    
    x <- as.numeric(x)
    
    if(length(x) == 1){
      return(rep(x, nages))
    }
    
    if(length(x) != nages){
      stop(name, " has length ", length(x), " but expected ", nages, ".")
    }
    
    x
  }
  
  if(is.null(bio$age)){
    stop("bio$age is missing.")
  }
  
  age <- as.numeric(bio$age)
  nages <- length(age)
  
  M   <- get_age_col(bio$M,   nages, female_col, "bio$M")
  Mat <- get_age_col(bio$Mat, nages, female_col, "bio$Mat")
  Wa  <- get_age_col(bio$Wa,  nages, female_col, "bio$Wa")
  
  if(any(!is.finite(M)) || any(M < 0)){
    stop("bio$M must contain finite non-negative natural mortality values.")
  }
  
  if(any(!is.finite(Mat)) || any(Mat < 0)){
    stop("bio$Mat must contain finite non-negative maturity values.")
  }
  
  if(any(!is.finite(Wa)) || any(Wa < 0)){
    stop("bio$Wa must contain finite non-negative weight-at-age values.")
  }
  
  survivors <- numeric(nages)
  survivors[1] <- 1
  
  if(nages > 1){
    for(a in 2:nages){
      survivors[a] <- survivors[a - 1] * exp(-M[a - 1])
    }
  }
  
  if(plus_group){
    expZ <- exp(-M[nages])
    survivors[nages] <- survivors[nages] / max(1e-12, 1 - expZ)
  }
  
  out <- survivors * Mat * Wa
  
  if(!byage){
    out <- sum(out, na.rm = TRUE)
  }
  
  return(out)
}


#' jbleslie()
#'
#' Estimate intrinsic population growth rate from a Leslie matrix.
#'
#' This function uses demographic life-history information returned by
#' \code{jbio()} to approximate the intrinsic rate of population increase
#' \code{r}. It follows the demographic logic described by McAllister et al.
#' (2001), where fecundity-at-age enters the first row of the Leslie matrix
#' and survivorship is represented separately by the age-transition elements.
#'
#' If \code{bio$fec} is absent, the function switches to a steepness-based
#' mode. In this case, fecundity-at-age is not supplied directly. Instead,
#' an implied replacement scaling is derived from Beverton-Holt steepness
#' \code{h} and unfished spawning output per recruit:
#'
#' \deqn{
#'   R_s = \frac{4h}{SPR_0(1-h)}
#' }
#'
#' and the implied fecundity-at-age schedule is:
#'
#' \deqn{
#'   m_x = R_s W_x Mat_x
#' }
#'
#' where \eqn{W_x} is weight-at-age and \eqn{Mat_x} is maturity-at-age.
#' This mode is mainly intended for teleost-style life-history priors where
#' fecundity is represented indirectly through spawning biomass and steepness,
#' rather than explicit pups or eggs per female.
#'
#' In contrast, if \code{bio$fec} is supplied, \code{jbleslie()} uses the
#' fecundity-based demographic approach directly and ignores \code{bio$h}.
#'
#' The key quantities are:
#' \itemize{
#'   \item \code{lambda}: dominant eigenvalue of the Leslie matrix. This is the
#'   finite annual population growth rate.
#'   \item \code{r}: intrinsic rate of increase, calculated as \code{log(lambda)}.
#'   \item \code{gt}: generation time, calculated as the mean age of reproduction,
#'   weighted by survivorship and female offspring production.
#'   \item \code{L}: the Leslie matrix, optionally returned when
#'   \code{return_L = TRUE}.
#' }
#'
#' If \code{bio$fec} is present, it is interpreted as fecundity-at-age before
#' applying the female offspring ratio. In the usual \code{jbio()} workflow,
#' scalar fecundity is expanded as maturity-at-age multiplied by fecundity of
#' mature females. If \code{bio$fec} is absent, fecundity-at-age is derived from
#' Beverton-Holt steepness \code{bio$h} and unfished spawning output per recruit.
#'
#' @param bio list returned by \code{jbio()}.
#' @param female_col integer. Column to use for female biology when matrices
#'   are sex-specific. Default is \code{1}.
#' @param offspring_female_ratio numeric. Proportion of offspring that are
#'   female. Use \code{0.5} if fecundity is total offspring per mature female.
#'   Use \code{1} if fecundity is already expressed as female offspring per
#'   mature female. Default is \code{0.5}.
#' @param plus_group logical. If \code{TRUE}, the terminal age is treated as a
#'   plus group with survivors remaining in the final age class. If \code{FALSE},
#'   the terminal age is treated as a true maximum age with senescence. Default
#'   is \code{TRUE}.
#' @param return_L logical. If \code{TRUE}, return the Leslie matrix. Default is
#'   \code{FALSE}.
#'
#' @return A list with:
#' \tabular{ll}{
#'   \code{r} \tab intrinsic rate of increase, \code{log(lambda)} \cr
#'   \code{gt} \tab generation time \cr
#'   \code{lambda} \tab dominant eigenvalue of the Leslie matrix \cr
#'   \code{L} \tab Leslie matrix, only returned if \code{return_L = TRUE}
#' }
#'
#' @examples
#' # Fecundity-based example: total pups or offspring per mature female
#' bio_fec <- jbio(
#'   amin = 0,
#'   amax = 20,
#'   nsexes = 1,
#'   Linf = 154.7,
#'   k = 0.129,
#'   t0 = -1,
#'   aW = 0.01,
#'   bW = 3.04,
#'   mat = c(5.5, 6, 1),
#'   M = 0.26,
#'   fec = 32
#' )
#'
#' jbleslie(bio_fec)
#'
#' # Return the Leslie matrix
#' out <- jbleslie(bio_fec, return_L = TRUE)
#' out$L
#' 
#' # Steepness-based teleost-style example
#' # Here fecundity is not supplied directly. Instead, jbleslie()
#' # derives an implied fecundity-at-age schedule from h, maturity,
#' # weight-at-age and unfished spawning output per recruit.
#' 
#' bio_h <- jbio(
#'   amin = 0,
#'   amax = 12,
#'   nsexes = 1,
#'   Linf = 80,
#'   k = 0.2,
#'   t0 = -0.5,
#'   aW = 0.00001,
#'   bW = 3.0,
#'   mat = c(35, 38, 0),
#'   M = 0.3,
#'   h = 0.76,
#'   fec = NULL
#' )
#'
#' jbleslie(bio_h)
#'
#' # Inspect the implied Leslie matrix
#' out_h <- jbleslie(bio_h, return_L = TRUE)
#' out_h$r
#' out_h$gt
#' out_h$lambda
#' out_h$L
#' 
#' # Sex-specific biology; use female column
#' bio_2sex <- jbio(
#'   amin = 0,
#'   amax = 12,
#'   nsexes = 2,
#'   Linf = c(90, 75),
#'   k = c(0.17, 0.22),
#'   t0 = c(-0.5, -0.5),
#'   aW = c(0.000012, 0.000015),
#'   bW = c(3.07, 3.01),
#'   mat = c(45, 55, 0),
#'   M = c(0.25, 0.32),
#'   h = 0.76
#' )
#'
#' jbleslie(bio_2sex, female_col = 1)
#'
#' @export
jbleslie <- function(
    bio,
    female_col = 1,
    offspring_female_ratio = 0.5,
    plus_group = TRUE,
    return_L = FALSE
){
  
  get_age_col <- function(x, nages, col = 1, name = "object"){
    if(is.null(x)) return(NULL)
    
    if(is.matrix(x) || is.data.frame(x)){
      x <- as.matrix(x)
      
      if(ncol(x) < col){
        stop(name, " has fewer columns than female_col.")
      }
      
      if(nrow(x) == 1){
        return(rep(as.numeric(x[1, col]), nages))
      }
      
      if(nrow(x) != nages){
        stop(name, " has ", nrow(x), " rows but expected ", nages, ".")
      }
      
      return(as.numeric(x[, col]))
    }
    
    x <- as.numeric(x)
    
    if(length(x) == 1){
      return(rep(x, nages))
    }
    
    if(length(x) != nages){
      stop(name, " has length ", length(x), " but expected ", nages, ".")
    }
    
    x
  }
  
  if(is.null(bio$age)){
    stop("bio$age is missing.")
  }
  
  age <- as.numeric(bio$age)
  nages <- length(age)
  
  M   <- get_age_col(bio$M,   nages, female_col, "bio$M")
  mat <- get_age_col(bio$Mat, nages, female_col, "bio$Mat")
  wa  <- get_age_col(bio$Wa,  nages, female_col, "bio$Wa")
  
  if(any(!is.finite(M)) || any(M < 0)){
    stop("bio$M must contain finite non-negative natural mortality values.")
  }
  
  # Survivorship lx, with l0 = 1.
  lx <- numeric(nages)
  lx[1] <- 1
  
  if(nages > 1){
    for(a in 2:nages){
      lx[a] <- lx[a - 1] * exp(-M[a - 1])
    }
  }
  
  # Optional plus-group correction for life-table calculations.
  lx_gt <- lx
  if(plus_group){
    lx_gt[nages] <- lx_gt[nages] / max(1e-12, 1 - exp(-M[nages]))
  }
  
  # Fecundity-at-age entering the first row of the Leslie matrix.
  # This is m_x in McAllister et al. (2001): m_x = s * g_x * f_x.
  if(is.null(bio$fec)){
    
    if(is.null(bio$h)){
      stop("bio$fec is NULL and bio$h is missing. Need either fecundity or steepness.")
    }
    
    spr0   <- jbspr0(bio, byage = FALSE)
    spr0_a <- jbspr0(bio, byage = TRUE)
    
    if(!is.finite(spr0) || spr0 <= 0){
      stop("Unfished spawning output per recruit, spr0, must be positive.")
    }
    
    # Beverton-Holt replacement scaling used for steepness-based priors.
    # This creates an implied fecundity-at-age schedule.
    rs <- 4 * bio$h / (spr0 * (1 - bio$h))
    
    fec_at_age <- rs * wa * mat
    
    gt <- sum(age * spr0_a, na.rm = TRUE) / spr0
    
  } else {
    
    fec <- get_age_col(bio$fec, nages, female_col, "bio$fec")
    
    # Convert total offspring to female offspring unless fecundity is already
    # expressed as female offspring per female.
    fec_at_age <- offspring_female_ratio * fec
    
    R0 <- sum(lx_gt * fec_at_age, na.rm = TRUE)
    
    if(!is.finite(R0) || R0 <= 0){
      stop("Lifetime female reproductive output R0 must be positive.")
    }
    
    gt <- sum(age * lx_gt * fec_at_age, na.rm = TRUE) / R0
  }
  
  if(any(!is.finite(fec_at_age)) || any(fec_at_age < 0)){
    stop("fecundity-at-age must contain finite non-negative values.")
  }
  
  # Leslie matrix.
  # First row = fecundity-at-age m_x.
  # Subdiagonal = annual survival S_x = exp(-M_x).
  L <- matrix(0, nrow = nages, ncol = nages)
  L[1, ] <- fec_at_age
  
  if(nages > 1){
    for(a in 2:nages){
      L[a, a - 1] <- exp(-M[a - 1])
    }
  }
  
  # Optional terminal plus group.
  # If FALSE, terminal age is treated as senescence / true maximum age.
  if(plus_group){
    L[nages, nages] <- exp(-M[nages])
  }
  
  lambda <- max(Re(eigen(L, only.values = TRUE)$values), na.rm = TRUE)
  r <- log(lambda)
  
  out <- list(
    r = r,
    gt = gt,
    lambda = lambda
  )
  
  if(return_L){
    out$L <- L
  }
  
  out
}

