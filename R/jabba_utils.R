#' gets beta prior paramenters
#'
#' gets beta prior paramenters and plots distribution
#' @param mu mean of beta distribution
#' @param CV cv of beta distribution
#' @param Min min of x-axis range (default = 0)
#' @param Prior X-axis label
#' @param Plot c(TRUE,FALSE)
#' @return a and b parameter (shape and scale)
#' @export
get_beta <- function(mu,CV,Min=0,Prior="x",Plot=FALSE){
  a = seq(0.0001,1000,0.001)
  b= (a-mu*a)/mu
  s2 = a*b/((a+b)^2*(a+b+1))
  sdev = sqrt(s2)
  # find beta )parameter a
  CV.check = (sdev/mu-CV)^2
  a = a[CV.check==min(CV.check)]
  #find beta parameter b
  b = (a-mu*a)/mu
  x = seq(Min,1,0.001)
  pdf = dbeta(x,a,b)
  if(Plot==TRUE){
    plot(x,pdf,type="l",xlim=range(x[pdf>0.01]),xlab=paste(Prior),ylab="",yaxt="n")
    polygon(c(x,rev(x)),c(rep(0,length(x)),rev(ifelse(pdf==Inf,100000,pdf))),col="grey")
  }
  return(c(a,b))
}

#' gets gamma prior paramenters
#'
#' gets gamma prior paramenters and plots distribution
#' @param mu mean of gamma distribution
#' @param CV cv of gamma distribution
#' @param Prior X-axis label
#' @param Plot c(TRUE,FALSE)
#' @return a and b parameter (shape and scale)
#' @export
get_gamma <- function(mu,CV,Prior="x", Plot=FALSE){
  a = seq(0.00001,10000,0.0001)
  b = a/mu
  s2 = (a/b^2)
  sdev = sqrt(s2)
  # find beta )parameter a
  CV.check = (sdev/mu-CV)^2
  a = a[CV.check==min(CV.check)]
  #find beta parameter b
  b = a/mu
  x = sort(rgamma(1000,a,b))
  pdf = dgamma(x,a,b)
  if(Plot==TRUE){
    plot(x,pdf,type="l",xlim=range(x[pdf>0.01]),xlab=paste(Prior),ylab="",yaxt="n")
    polygon(c(x,rev(x)),c(rep(0,length(x)),rev(ifelse(pdf==Inf,100000,pdf))),col="grey")
  }
  return(c(a,b))
}


#' gets lognormal prior paramenters
#'
#' gets lognormal prior paramenters mean and log.sd with bias correction and plots distribution
#' @param mu mean of lognormal distribution
#' @param CV cv of lognoram distribution
#' @param Prior X-axis label
#' @param Plot c(TRUE,FALSE)
#' @return mean and lod.sd
#' @export
plot_lnorm <- function(mu,CV,Prior="x",Plot=FALSE){
  sdev= sqrt(log(CV^2+1))
  rand.pr = rlnorm(1000,log(mu)-0.5*sdev^2,sdev)
  x = seq(min(rand.pr),quantile(rand.pr,0.995),max(rand.pr/500))
  pdf = dlnorm(x,log(mu),sdev)
  if(Plot==TRUE){
    plot(x,pdf,type="l",xlim=range(x),xlab=paste(Prior),ylab="",yaxt="n")
    polygon(c(x,rev(x)),c(rep(0,length(x)),rev(ifelse(pdf==Inf,100000,pdf))),col="grey")
  }
  return(c(mu,sdev))
}

#' Gets lognormal prior paramenters
#'
#' Gets lognormal prior paramenters mean and log.sd with bias correction and plots distribution
#' @param mu mean of lognormal distribution
#' @param CV cv of lognoram distribution
#' @param Prior X-axis label
#' @param Plot c(TRUE,FALSE)
#' @return mean and lod.sd
#' @export
plot_lnorm <- function(mu,CV,Prior="x",Plot=FALSE){
  sdev= sqrt(log(CV^2+1))
  rand.pr = rlnorm(1000,log(mu)-0.5*sdev^2,sdev)
  x = seq(min(rand.pr),quantile(rand.pr,0.995),max(rand.pr/500))
  pdf = dlnorm(x,log(mu),sdev)
  if(Plot==TRUE){
    plot(x,pdf,type="l",xlim=range(x),xlab=paste(Prior),ylab="",yaxt="n")
    polygon(c(x,rev(x)),c(rep(0,length(x)),rev(ifelse(pdf==Inf,100000,pdf))),col="grey")
  }
  return(c(mu,sdev))
}

#' Function kobeJabba for FLR
#'
#' Function to convert kobe posteriors into KOBE FLR input object
#' @param x posterior array dims(iter,year,stock,harvest)
#' @export
kobeJabba<-function(x){

  out=cbind(reshape::melt(x[,,2]),c(x[,,3]))
  names(out)=c("iter","year","stock","harvest")
  out$year=out$year
  out}

#' Function to convert JABBA projections into  FLR object
#'
#' Function to convert kobe projection matrix posteriors into Kobe FLR input object
#' @param x posterior array dims(iter,year,tac,stock,harvest)
#' @export
kobeJabbaProj<-function(x){

  out=cbind(reshape::melt(x[,,,2]),reshape::melt(x[,,,3])[,4],reshape::melt(x[,,,1])[,4])
  names(out)=c("iter","year","tac","stock","harvest","bk")
  out$year=out$year

  out}

#' Function to do runs.test and 3 x sigma limits
#'
#' runs test is conducted with library(snpar)
#' @param x residuals from CPUE fits
#' @param type only c("resid","observations")
#' @return runs p value and 3 x sigma limits
#' @export
runs_sig3 <- function(x,type=NULL) {
  if(is.null(type)) type="resid"
  if(type=="resid"){mu = 0}else{mu = mean(x, na.rm = TRUE)}
  # Average moving range
  mr  <- abs(diff(x - mu))
  amr <- mean(mr, na.rm = TRUE)
  # Upper limit for moving ranges
  ulmr <- 3.267 * amr
  # Remove moving ranges greater than ulmr and recalculate amr, Nelson 1982
  mr  <- mr[mr < ulmr]
  amr <- mean(mr, na.rm = TRUE)
  # Calculate standard deviation, Montgomery, 6.33
  stdev <- amr / 1.128
  # Calculate control limits
  lcl <- mu - 3 * stdev
  ucl <- mu + 3 * stdev
  if(nlevels(factor(sign(x)))>1){
    runstest = snpar::runs.test(x)
    pvalue = round(runstest$p.value,3)} else {
      pvalue = 0.001
    }

  return(list(sig3lim=c(lcl,ucl),p.runs= pvalue))
}

#' normIndex()
#'
#' Function function to normalize CPUE indices
#' @param cpue first column must be year
#' @return normalized cpue
#' @export
normIndex <- function(cpue){
  cpue[,-1] = t(t(as.matrix(cpue[,-1]))/apply(as.matrix(cpue[,-1]),2,mean,na.rm=TRUE))
  return(cpue)
}


#' rc4 r4ss color palette 
#' @param n number of colors
#' @param alpha transluscency 
#' @return vector of color codes
#' @export
rc4 <- function(n,alpha=1){
  # a subset of rich.colors by Arni Magnusson from the gregmisc package
  # a.k.a. rich.colors.short, but put directly in this function
  # to try to diagnose problem with transparency on one computer
  x <- seq(0, 1, length = n)
  r <- 1/(1 + exp(20 - 35 * x))
  g <- pmin(pmax(0, -0.8 + 6 * x - 5 * x^2), 1)
  b <- dnorm(x, 0.25, 0.15)/max(dnorm(x, 0.25, 0.15))
  rgb.m <- matrix(c(r, g, b), ncol = 3)
  rich.vector <- apply(rgb.m, 1, function(v) rgb(v[1], v[2], v[3], alpha=alpha))
}



#' ss3col r4ss color generator 
#' @param n number of colors
#' @param alpha transluscency 
#' @return vector of color codes
#' @export
ss3col <- function(n,alpha=1){
  if(n>3) col <- rc4(n+1)[-1]
  if(n<3)  col <- rc4(n)
  if(n==3) col <- c("blue","red","green3")
  if(alpha<1){
    # new approach thanks to Trevor Branch
    cols <- adjustcolor(col, alpha.f=alpha)
  } else {
    cols=col
  }
  return(cols)
}


#' jbmase()
#'
#' Computes Mean Absolute Scaled Errors as a measure of prediction skill
#' @param hc object list of hindcasts from hindcast_jabba() or jbhcxval()
#' @param naive.min minimum MASE denominator (naive predictions) for MASE.adj (default = 0.1)
#' @param index option to compute for specific indices (numeric & in order)
#' @param verbose if FALSE run silent
#' @return hc containing estimates of key joint results from all hindcast run 
#' @export

jbmase <- function(hc,naive.min=0.1,index=NULL,verbose=TRUE){
  MASE = NULL
  d. = do.call(rbind,lapply(hc,function(x){
    x$diags}))
  
  xmin= min(d.$year)
  all.indices = unique(d.$name)  
  if(is.null(index)) index = 1:length(all.indices)
  # subset 
  d. = d.[d.$name%in%all.indices[index],]
  peels = unique(d.$retro.peels)
  styr = max(hc[[1]]$yr)-max(peels)
  years = min(d.$year):max(d.$year)
  yr = unique(d.$year)
  endyrvec = rev(sort(years[length(years)-peels]))
  
  if(verbose)cat("\n","><> Only including indices that have years overlapping hind-cast horizan","\n")
  # check in index
  indices = unique(d.$name)
  valid = NULL
  for(i in 1:length(indices)){
    if(nrow(d.[d.$name%in%indices[i] & d.$year>styr & d.$retro.peels%in%peels[1],])>1){ # Only run if overlap
      valid=c(valid,paste(indices[i]))}
  }
  if(verbose) cat("\n","><> Including indices:",valid,"\n")
  n.indices = length(valid)  
  for(i in 1:length(indices)){
    if(nrow(d.[d.$name%in%indices[i] & d.$year>styr & d.$retro.peels%in%peels[1],])>1){ # Only run if overlap
      xv = d.[d.$name%in%indices[i],]
      yr = unique(xv$year)
      yr.eval <- sort(endyrvec)
      yr.obs <- yr.eval%in%yr
      pe.eval = which(yr.eval%in%yr)[-1]
      if(length(which(yr.eval%in%yr))-length(pe.eval)<1){
        pe.eval = pe.eval[-1]
      } 
      npe <- length(pe.eval)  # number of prediction errors
      obs.eval <- rep(NA,length(yr.eval))
      obs.eval[yr.eval%in%yr] = xv$obs[xv$retro.peels==min(xv$retro.peels)][yr%in%yr.eval]
      if(is.na(obs.eval[1]))
        nhc = length(endyrvec)-1
      
      py = xv$year[xv$retro.peels==min(xv$retro.peels) & xv$year>styr-xmin]
      obs =xv$obs[xv$retro.peels==min(xv$retro.peels) & xv$year>styr-xmin]
      
      if(verbose) cat(paste("\n","Computing MASE with",ifelse(npe<(length(endyrvec)-1),"only","all"),
                            npe,"of",length(endyrvec)-1," prediction residuals for Index",xv$name[1]),"\n")
      if(verbose & npe<(length(endyrvec)-1))cat(paste("\n","Warning:  Unequal spacing of naive predictions residuals may influence the interpretation of MASE","\n","\n"))
      
      naive.eval=log(obs.eval[is.na(obs.eval)==F][-length(obs.eval[is.na(obs.eval)==F])])-log(obs.eval[is.na(obs.eval)==F][-1])
      nhc = length(endyrvec)-1
      pred.resid = NULL
      for(j in 1:(length(peels)-1)){
        if(endyrvec[j] %in% xv$year){
          
          x <- min(py):max(yr.eval)
          x <- x[1:(length(x)-peels[j])]
          x = x[x%in%unique(xv$year)]
          y <- xv[xv$retro.peels==peels[j+1] & xv$year%in%x,]$hat
          pred.resid = c(pred.resid,log(y[length(x)])-log(obs[length(x)])) # add log() for v1.1
        }}
      
      maepr =  mean(abs(pred.resid))
      if(is.na(obs.eval[1])) obs.eval[1] =  rev(obs[obs%in%obs.eval==F])[1] 
      scaler = mean(abs(naive.eval))
      scaler.adj = mean(pmax(abs(naive.eval),naive.min))
      
      MASE.i = NULL
      MASE.i = data.frame(Index=unique(xv$name)[1], MASE=maepr/scaler,MASE.adj=maepr/scaler.adj,MAE.PR=maepr,MAE.base=scaler,n.eval=npe)
      
    } else{
      xv = d.[d.$name%in%indices[i],]
      if(verbose) cat(paste0("\n","No observations in evaluation years to compute prediction residuals for Index ",xv $name[1]),"\n")
      MASE.i = NULL
      MASE.i = data.frame(Index=unique(xv$name)[1], MASE=NA,MASE.adj=NA,MAE.PR=NA,MAE.base=NA,n.eval=0)  
    }
    MASE = rbind(MASE,MASE.i)
    
  } # end of index loop
  
  
  return(MASE)
}


#' zage()
#'
#' Computes the Z[t] from catch-at-age data
#' @param ca data.frame input with column names year, age, data
#' @param ages for which the z slope is taken
#' @return data.frame with data = z
#' @export
zage =  function(ca,ages = "missing"){
  out=NULL
  CN = catch.n
  years=unique(CN$year)
  age = unique(CN$age)
  
  if(missing(ages)){
    Amin = age[which(aggregate(data~age,CN,mean)$data==max(aggregate(data~age,CN,mean)$data))] 
    Amax = max(CN$age)-1
    ages= Amin:Amax
  }
  for(t in 1:length(years)){
    Ct= CN[CN$year==years[t]&CN$age%in%ages,]
    Z = -coef(lm(log(data)~age,Ct))[[2]]
    zt = Ct[1,]
    zt$data = Z 
    zt$age = "all"
    out = rbind(out,zt)
  }
  return(out)
}

