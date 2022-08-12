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
#' runs test is conducted with library(randtests)
#' @param x residuals from CPUE fits
#' @param type only c("resid","observations")
#' @param mixing c("less","greater","two.sided"). Default less is checking for postive autocorrelation only    
#' @return runs p value and 3 x sigma limits
#' @export
#' @author Henning Winker (JRC-EC) and Laurence Kell (Sea++)
jbruns_sig3 <- function(x,type=NULL,mixing="less") {
  if(is.null(type)) type="resid"
  if(type=="resid"){
    mu = 0}else{mu = mean(x, na.rm = TRUE)}
  alternative=c("two.sided","left.sided","right.sided")[which(c("two.sided", "less","greater")%in%mixing)]
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
    # Make the runs test non-parametric
    runstest = randtests::runs.test(x,threshold = 0,alternative = alternative)
    if(is.na(runstest$p.value)) p.value =0.001
    pvalue = round(runstest$p.value,3)} else {
      pvalue = 0.001
    }
  
  return(list(sig3lim=c(lcl,ucl),p.runs= pvalue))
}



#' JABBA runs test 
#'
#' Residual diagnostics with runs test p-value 
#' @param jabba output list from fit_jabba
#' @param mixing c("less","greater","two.sided"). Default "less" is checking for positive autocorrelation only
#' @param index option to plot specific indices (numeric & in order)
#' @export
#' @examples 
#' data(iccat)
#' bet= iccat$bet
#' jb = build_jabba(catch=bet$catch,cpue=bet$cpue,se=bet$se,assessment="BET",scenario = "Ref",model.type = "Pella",igamma = c(0.001,0.001),verbose=FALSE)
#' fit = fit_jabba(jb,quickmcmc=TRUE,verbose=FALSE)
#' jbrunstest(fit)
#' jbrunstest(fit,index=2)
#' jbplot_runstest(fit,verbose=FALSE)

jbrunstest <- function(jabba,index=NULL,mixing="less"){
    
    all.indices = unique(jabba$diags$name)
    if(is.null(index)) index = 1:length(all.indices)
    indices = unique(jabba$diags$name)[index]
    n.indices = length(indices)
    Resids = jabba$residuals
    years = jabba$yr
    check.yrs = abs(apply(jabba$residuals,2,sum,na.rm=TRUE))
    cpue.yrs = years[check.yrs>0]
    
      runs = NULL
      for(i in 1:n.indices){
        resid = (Resids[index[i],is.na(Resids[index[i],])==F])
        res.yr = years[is.na(Resids[index[i],])==F]
        if(length(resid)>3){
        runstest = jbruns_sig3(x=as.numeric(resid),type="resid",mixing=mixing)
        runs = rbind(runs,c(runstest$p.runs,runstest$sig3lim))
        } else {
          runs = rbind(runs,c(NA,NA,NA))
        }}
        
        runstable = data.frame(Index=indices,runs.p=as.matrix(runs)[,1],Test=ifelse(is.na(as.matrix(runs)[,1]),"Excluded",ifelse(as.matrix(runs)[,1]<0.05,"Failed","Passed")),sigma3.lo=as.matrix(runs)[,2],sigma3.hi=as.matrix(runs)[,3]) 
        colnames(runstable) = c("Index","runs.p","test","sigma3.lo","sigma3.hi")
        
        return(runstable)
      
} # end of runstest function


#' jbretro() computes Mohn's rho and forecast rho
#'
#' Quantities retrospective pattern of B, F, BBmsy, FFmsy, BB0 and SP #'
#' @param hc output list from hindast_jabba()
#' @param type option c("B","F","BBmsy","FFmsy","BB0","SP")
#' @param forecast  includes retrospective forecasting if TRUE
#' @return Mohn's rho statistic for several quantities
#' @export
#' @examples 
#' data(iccat)
#' bet= iccat$bet
#' jb = build_jabba(catch=bet$catch,cpue=bet$cpue,se=bet$se,assessment="BET",scenario = "Ref",model.type = "Pella",igamma = c(0.001,0.001),verbose=FALSE)
#' fit = fit_jabba(jb,quickmcmc=TRUE,verbose=FALSE)
#' hc = hindcast_jabba(jbinput=jb,fit=fit,peels=1:3)
#' jbretro(hc)
#' jbplot_retro(hc)
#' jbplot_retro(hc,forecast=TRUE) # with retro forecasting

jbretro <- function(hc,type=c("B","F","BBmsy","FFmsy","procB","SP"),forecast=TRUE){
  
  hc.ls = hc 
  
  peels = as.numeric(do.call(c,lapply(hc.ls,function(x){x$diags$retro.peels[1]})))
  Ref = hc.ls[[1]]
  hc = list(scenario = Ref$scenario, yr=Ref$yr,catch=Ref$catch,peels=NULL,timeseries = NULL,refpts=NULL,pfunc=NULL,diags=NULL,settings=Ref$settings)
  for(i in 1:length(peels)){
    hc.ls[[i]]$pfunc$level = peels[i] 
    hc.ls[[i]]$refpts$level = peels[i]
    hc$timeseries$mu = rbind(hc$timeseries$mu,data.frame(factor=hc.ls[[i]]$diags[1,1],level=peels[i],hc.ls[[i]]$timeseries[,"mu",])) 
    hc$timeseries$lci = rbind(hc$timeseries$lci,data.frame(factor=hc.ls[[i]]$diags[1,1],level=peels[i],hc.ls[[i]]$timeseries[,"lci",])) 
    hc$timeseries$uci = rbind(hc$timeseries$uci,data.frame(factor=hc.ls[[i]]$diags[1,1],level=peels[i],hc.ls[[i]]$timeseries[,"uci",])) 
    hc$diags = rbind(hc$diags,hc.ls[[i]]$diags)
    hc$refpts= rbind(hc$refpts,hc.ls[[i]]$refpts[1,])
    hc$pfunc= rbind(hc$pfunc ,hc.ls[[i]]$pfunc)
  }
  
  retros = unique(peels)
  runs= hc$timeseries$mu$level
  years= hc$yr
  nyrs = length(years)
  FRP.rho = c("B","F", "Bmsy", "Fmsy", "procB","MSY")  
  rho = data.frame(mat.or.vec(length(retros)-1,length(FRP.rho)))
  colnames(rho) = FRP.rho
  fcrho = rho
  
   for(k in 1:length(type)){
      j = which(c("B","F","BBmsy","FFmsy","BB0","procB","SP")%in%type[k])
      if(type[k]%in%c("B","F","BBmsy","FFmsy","procB")){
        y = hc$timeseries$mu[,j+2]
        ref = hc$timeseries$mu[runs%in%retros[1],j+2]
        ylc = hc$timeseries$lci[runs%in%retros[1],j+2]
        yuc = hc$timeseries$uci[runs%in%retros[1],j+2]
        for(i in 1:length(retros)){
          
          if(i>1){
            rho[i-1,k] =  (y[runs%in%retros[i]][(nyrs-retros[i])]-ref[(nyrs-retros[i])])/ref[(nyrs-retros[i])]
            fcrho[i-1,k] = (y[runs%in%retros[i]][(nyrs+1-retros[i])]-ref[(nyrs+1-retros[i])])/ref[(nyrs+1-retros[i])]
            if(type[k]=="procB"){
              rho[i-1,k] =  (exp(y[runs%in%retros[i]][(nyrs-retros[i])])-exp(ref[(nyrs-retros[i])]))/exp(ref[(nyrs-retros[i])])
              fcrho[i-1,k] =  (exp(y[runs%in%retros[i]][(nyrs+1-retros[i])])-exp(ref[(nyrs+1-retros[i])]))/exp(ref[(nyrs+1-retros[i])])
            }
          }
        }
      } else {
        # Plot SP
        for(i in 1:length(retros)){
          if(i>1){
            rho[i-1,6] =  (hc$refpts$msy[hc$refpts$level==retros[i]]-hc$refpts$msy[hc$refpts$level==retros[1]])/hc$refpts$msy[hc$refpts$level==retros[1]]
            fcrho[i-1,6] = NA
          }
        }}
   } # end k    
  rho = rbind(rho,apply(rho,2,mean))
  rownames(rho) = c(rev(years)[retros[-1]],"rho.mu")
  
  fcrho = rbind(fcrho,apply(fcrho,2,mean))
  rownames(fcrho) = c(rev(years)[retros[-1]],"forecastrho.mu")
  if(forecast){
    out = list()
    out$Mohns.rho = rho
    out$Forecast.rho = fcrho
  } else {
    out = rho
  }
  
  return(out)
} # end of Retrospective Plot



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
    if(n==1) cols = "darkgrey"#adjustcolor("darkgrey", alpha.f=alpha)
  } else {
    cols=col
    if(n==1) cols <- "black"
  }
  
  return(cols)
}


#' jbmase()
#'
#' Computes Mean Absolute Scaled Errors as a measure of prediction skill
#' @param hc object list of hindcasts from hindcast_jabba() or jbhcxval()
#' @param naive.min minimum MASE denominator (naive predictions) for MASE.adj (default = 0.1)
#' @param index option to compute for specific indices (numeric & in order)
#' @param residuals if TRUE, outputs individual prediction and naive residuals   
#' @param verbose if FALSE run silent
#' @return hc containing estimates of key joint results from all hindcast run 
#' @export

jbmase <- function(hc,naive.min=0.1,index=NULL,residuals=FALSE,verbose=TRUE){
  MASE = Residuals = NULL
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
        if(endyrvec[1:length( naive.eval)][j] %in% xv$year){
          
          x <- min(py):max(yr.eval)
          x <- x[1:(length(x)-peels[j])]
          x = x[x%in%unique(xv$year)]
          y <- xv[xv$retro.peels==peels[j+1] & xv$year%in%x,]$hat
          pred.resid = c(pred.resid,log(y[length(x)])-log(obs[length(x)])) # add log() for v1.1
        }}
      
      pred.resid = pred.resid[1:length( naive.eval)]
      maepr =  mean(abs(pred.resid))
      if(is.na(obs.eval[1])) obs.eval[1] =  rev(obs[obs%in%obs.eval==F])[1] 
      scaler = mean(abs(naive.eval))
      scaler.adj = mean(pmax(abs(naive.eval),naive.min))
      res.i = NULL
      res.i = data.frame(Index=unique(xv$name)[1],Year=yr.eval[pe.eval],Pred.Res=pred.resid,Naive.Res=naive.eval,n.eval=npe) 
      MASE.i = NULL
      MASE.i = data.frame(Index=unique(xv$name)[1], MASE=maepr/scaler,MASE.adj=maepr/scaler.adj,MAE.PR=maepr,MAE.base=scaler,n.eval=npe)
      Residuals = rbind(Residuals,res.i)  
    } else{
      xv = d.[d.$name%in%indices[i],]
      if(verbose) cat(paste0("\n","No observations in evaluation years to compute prediction residuals for Index ",xv $name[1]),"\n")
      MASE.i = NULL
      MASE.i = data.frame(Index=unique(xv$name)[1], MASE=NA,MASE.adj=NA,MAE.PR=NA,MAE.base=NA,n.eval=0)  
    }
    MASE = rbind(MASE,MASE.i)
    
  } # end of index loop
  jstats = apply(abs(Residuals[c("Pred.Res","Naive.Res")]),2,mean)
  joint = data.frame(Index="joint",
                     MASE=jstats[1]/jstats[2],MAE.PR=jstats[1],MAE.base=jstats[2],
                     MASE.adj=jstats[1]/pmax(jstats[2],naive.min),n.eval=nrow(Residuals))  
  MASE = rbind(MASE,joint)
  if(!residuals) return(MASE)
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

#' addBfrac()
#'
#' add biomass reference to kb ouput as fraction Bmsy or B0, e.g. for Blim or MSST
#' @param jabba bfrac fraction of Bmsy or B0
#' @param base defines biomass base "bmsy" or "b0"
#' @param quantiles default is 95CIs as c(0.025,0.975)
#' @return 
#' @export
addBfrac <- function(kb, bfrac=0.5, bref = c("bmsy","b0"),quantiles = c(0.025,0.975)){
  if(!is.null(kb$settings)){ 
    if(is.null(jabba$kbtrj)) stop("rerun with fit_jabba(...,save.trj = TRUE)")
    kb = kb$kbtrj
  } 
  if(bref[1]=="bmsy"){
    kb$BBfrac = kb$stock/bfrac
    kb$Bref =  kb$B/(kb$stock/bfrac)
  }
  
  if(bref[1]=="b0"){
    kb$BBfrac = kb$BB0/bfrac
    kb$Bref =  kb$B/(kb$BB0/bfrac)
  }
  Bref = data.frame(t(quantile(kb$Bref,c(0.5,quantiles[1],quantiles[2]))))
  colnames(Bref) = c("mu","lci","uci")
  return(list(kb=kb,bref=Bref))
}


