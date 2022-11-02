#' catch
#'
#' Poly.sim catch data
#'
#'@docType data
#'
#'@format data.frame of 46 obs. of 8 variables
#'\describe{
#'\item{age}{age}
#'\item{year}{year}
#'\item{unit}{unit}
#'\item{season}{season}
#'\item{area}{area}
#'\item{iter}{iter}
#'\item{data}{data}
#'\item{qname}{qname}
#'}
#'
"catch"

#' catch.n
#' 
#' Ploy.sim catch by n
#' 
#' @docType data
#' 
#' @format data.frame with 1271 obs. of 7 variables
#'\describe{
#'\item{age}{age}
#'\item{year}{year}
#'\item{unit}{unit}
#'\item{season}{season}
#'\item{area}{area}
#'\item{iter}{iter}
#'\item{data}{data}
#' }
"catch.n"


#' index
#'
#' Poly.sim index
#'
#'@docType data
#'
#'@format data.frame of 46 obs. of 8 variables
#'\describe{
#'\item{age}{age}
#'\item{year}{year}
#'\item{unit}{unit}
#'\item{season}{season}
#'\item{area}{area}
#'\item{iter}{iter}
#'\item{data}{data}
#'\item{qname}{qname}
#'}
#'
"index"


#' om
#' 
#' poly.sim om
#' 
#' @docType data
#' 
#' @format List of 6
#' \describe{
#' \item{stock}{data.frame of 130 rows and 7 cols}
#' \item{harvest}{data.frame of 130 rows and 7 cols}
#' \item{catch}{data.frame of 130 rows and 7 cols}
#' \item{fmsy}{fmsy}
#' \item{bmsy}{bmsy}
#' \item{msy}{msy}
#' }
"om"


#' iccat
#'
#' iccat list
#' 
#' @docType data
#' 
#' @format List of 4 
#' \describe{
#' \item{bet}{Cpue, se, and catch data.frames from JR2 }
#' \item{whm}{Cpue, se, and catch data.frames from JPNLL priors, JPNLL, CTPLL, BRARR, BRALL, USARR, USALL, VENLL, and VENGL }
#' \item{swos}{Cpue, se, and catch data.frames from BRA_LL, SPA_LL, JPN_LL, URU_LL, SA_LL}
#' \item{bum}{Cpue, se, and catch data.frames from JAP.LL.hist, TAI, USA.LL, VEN.LL, VEN.GIL, VEN.Rec, BRA.LL, BRA.Rec, GHA.GIL}
#' }
"iccat"