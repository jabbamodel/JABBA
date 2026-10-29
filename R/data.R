# Documentation of the JABBA example data. All datasets share one structure: a list (or a list of stocks) with
# catch = data.frame(year, ...) by fleet, cpue = data.frame(year, ...) by index, se = data.frame(year, ...) (log-scale
# SE/CV) or NULL, and info = list(species, latin, stock, code, years, indices, units, ...); simulated stocks add truth.
# Build scripts: data-raw/jabba_data.R and data-raw/kobsim.R

#' ICCAT example stocks
#'
#' Catch and abundance-index data for four ICCAT stocks assessed with JABBA, in the common JABBA data format.
#' \itemize{
#'   \item \code{bet}: Atlantic bigeye tuna (\emph{Thunnus obesus}), 1950-2017, 2 indices
#'   \item \code{whm}: Atlantic white marlin (\emph{Kajikia albida}), 1959-2017, 13 indices
#'   \item \code{swos}: South Atlantic swordfish (\emph{Xiphias gladius}), 1950-2015, 7 indices
#'   \item \code{bum}: Atlantic blue marlin (\emph{Makaira nigricans}), 1959-2016, 10 indices
#' }
#' @format A named list of stocks; each stock is a list with
#' \describe{
#'   \item{catch}{data.frame(year, catch): total catch (t)}
#'   \item{cpue}{data.frame(year, index1, ...): relative abundance indices (NA = missing)}
#'   \item{se}{data.frame(year, index1, ...): log-scale standard errors of the indices}
#'   \item{info}{list(species, latin, stock, rfmo, code, years, indices, catch.unit, cpue.unit, se.type)}
#' }
#' @source ICCAT stock assessments (\url{https://www.iccat.int})
#' @examples
#' data(iccat)
#' names(iccat)
#' bet <- iccat$bet
#' jb <- build_jabba(catch = bet$catch, cpue = bet$cpue, se = bet$se, model.type = "Fox")
"iccat"

#' Simulated pollack example
#'
#' Data simulated with an FLR age-structured operating model (pollack-like life history), with catch, a vulnerable
#' biomass index (lognormal error, sd 0.25), catch-at-age for Z-based auxiliary indices (\code{zage()}), and truth.
#' @format A list with
#' \describe{
#'   \item{catch}{data.frame(year, Catch), years 70-115}
#'   \item{cpue}{data.frame(year, Index)}
#'   \item{se}{NULL}
#'   \item{catch.n}{catch-at-age in long FLR format (age, year, unit, season, area, iter, data), years 85-115}
#'   \item{truth}{data.frame(year, SB_SBmsy, F_Fmsy, catch) of the operating model}
#'   \item{refpts}{c(Fmsy, Bmsy, MSY)}
#'   \item{info}{list(species, latin, stock, code, years, indices, catch.unit, cpue.unit, note)}
#' }
#' @examples
#' data(pol.sim)
#' z <- zage(pol.sim$catch.n, ages = 4:12)
"pol.sim"

#' Simulated silver kob example (JABBA-Select)
#'
#' Simulated data for silver kob (\emph{Argyrosomus inodorus}) from the JABBA-Select example
#' (Winker et al. 2020, Fisheries Research 222: 105355). An age-structured operating model (40 years,
#' recruitment sd 0.6, Beverton-Holt steepness 0.8, M = 0.18) is fished by two fleets with logistic length
#' selectivity: S1 (SL50 = 300 mm) in years 1-25 and S2 (SL50 = 500 mm) from year 26 after a size-limit increase.
#' Two abundance indices of exploitable biomass are split by selectivity regime (A1.S1, A1.S2, A2.S1, A2.S2).
#'
#' @format A list with
#' \describe{
#'   \item{catch}{data.frame(year, S1, S2): catch by fleet (t)}
#'   \item{cpue, se}{indices by selectivity regime and their log-scale SEs}
#'   \item{fleets}{fleet table of JABBA-Select (selectivity and q assignment)}
#'   \item{selpars}{list of selectivity parameters c(SL50, SL95, SLdesc, CVdesc, Smin) for \code{jselex()}}
#'   \item{lh, stk}{life-history parameters and the \code{jbio()} object}
#'   \item{ebsb}{output of \code{make_ebsb()} for S1 and S2; \code{$pars} gives the curves for \code{index.ebsb}}
#'   \item{truth}{data.frame(year, SB, EB, F, catch, sel): operating-model truth}
#'   \item{SB0}{unfished spawning biomass}
#'   \item{info}{list(species, latin, stock, code, years, indices, units, note)}
#' }
#' @source \url{https://github.com/jabbamodel/JABBA-Select}
#' @examples
#' data(kobsim)
#' ebsb_matrix(kobsim$ebsb$pars)
"kobsim"
