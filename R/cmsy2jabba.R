#' cmsy.bkprior()
#'
#' approximates initial b/k prior using CMSY++15c.R
#' @param catch  catch time series, requires data.frame(year, catch)
#' @param bw kernel bandwith for catch smoothing (default=3)
#' @return list with initial start bk range and c(mu,cv) 
#' @export
cmsy.bkprior = function(catch,bw=3,prior.r=c(0.05,0.5)){

yr <- catch[,1]
ct.raw <- catch[,2]
nyr          <- length(yr) # number of years in the time series
# apply kernel smoothing with a bandwidth of bw
ct           <- ksmooth(x=yr,y=ct.raw,kernel="normal",n.points=length(yr),bandwidth=bw)$y
# get index of years with lowest and highest catch
min.yr.i     <- which.min(ct)
max.yr.i     <- which.max(ct)
min.ct       <- ct[min.yr.i]
max.ct       <- ct[max.yr.i]
mean.ct      <- mean(ct)
sd.ct        <- sd(ct)

# use mean of first 3 years to determine startbio
ct.3  <- mean(ct[1:3])
start.yr <- min(yr) # initialize / reset start.yr.new with NA
# use initial biomass range from input file if stated

# in recent years >= 1980 and for medium or high resilience species
# if catch < 0.1 max catch, assume very low biomass
if(ct.3 < 0.1*max.ct && yr[1]>=1980 && prior.r[2]>0.3) { startbio <- c(0.01,0.2)
# if catch < 0.25 max catch, assume low biomass
} else if(ct.3 < 0.25*max.ct && yr[1]>=1980 && prior.r[2]>0.3) { startbio <- c(0.01,0.3)
# if catch < 0.33 max catch, assume low biomass
} else if(ct.3 < 0.33*max.ct && yr[1]>=1980 && prior.r[2]>0.3) { startbio <- c(0.01,0.4)
# if time series started before 1980
# if catch < 0.1 max catch, assume nearly unexploited biomass
} else if(ct.3 < 0.1*max.ct) { startbio <- c(0.9,1)
# if catch < 0.25 max catch, assume high biomass
} else if(ct.3 < 0.25*max.ct) { startbio <- c(0.8,1)
# if catch < 0.33 max catch, assume high biomass
} else if(ct.3 < 0.33*max.ct) { startbio <- c(0.6,1)
# if catch < 0.66 max catch, assume medium to high biomass
} else if(ct.3 < 0.66*max.ct | start.yr <=1960) { startbio <- c(0.4,0.8)
# otherwise assume low to medium biomass
} else startbio <- c(0.2,0.6) 

   sd.bk = (startbio[2]-startbio[1])/(4*0.98)
   mu.bk = mean(startbio)
   cv.bk = sd.bk/mu.bk

  bkprior = list()
  bkprior$range = startbio
  bkprior$mu.cv = c(mu.bk,cv.bk)
  return(bkprior)
} # }}}

#' range2prior()
#'
#' approximates mu, cv, logsd from a range limits
#' @param lo lower limit
#' @param hi upper limit
#' @return  mu,cv and logsd of range limits  
#' @export
range2prior<-function(lo,hi){
  sdev = mean(hi-lo)/(4*0.98)
  mu = mean(c(lo,hi))
  cv = sdev/mu
  logsd = plot_lnorm(mu,CV=cv,Plot=F)[2]
  return(c(mu,cv,logsd))
}

#' cmsy.prior()
#'
#' Computes r range from FishBase resiliance category
#' @param res Resiliance: Very low, Low, Medium, High
#' @return  r prior range
#' @export
cmsy.rprior <- function(res){
  # initial range of r based on resilience
  if(res == "High") {
    prior.r <- c(0.6,1.5)} else if(res == "Medium") {
      prior.r <- c(0.2,0.8)}    else if(res == "Low") {
        prior.r <- c(0.05,0.5)}  else { # i.e. res== "Very low"
          prior.r <- c(0.015,0.1)} 
  return(prior.r) 
  }





