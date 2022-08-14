
#' jbpar()
#'
#' Set the par() to options suitable for jabba multi plots   
#' @param mfrow determines plot frame set up
#' @param plot.cex cex graphic option
#' @export
jbpar <- function(mfrow=c(1,1),mex=0.8,plot.cex=0.8,mai=c(0.5,0.5,0.15,.15),omi = c(0.2,0.2,0.2,0),mar=c(2.5,2.5, 0.7, 0.7),labs=TRUE){
  if(!labs)  mai=c(0.35,0.15,0,.15)
  par(list(mfrow=mfrow,mai = mai, mgp =c(1.5,0.5,0),omi = omi + 0.1,mar=mar, tck = -0.02,cex=plot.cex))
}

#' jbplot_indices
#'
#' Plots indices with assumed mimumum observation errors  
#' @param jabbainput list from build_jabba() 
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param width plot width
#' @param height plot hight
#' @param add if TRUE par is not called to enable manual multiplots
#' @param add.legend option to add legend
#' @param legend.loc place lengend
#' @param ylab option to change y-axis label
#' @param xlab option to change x-axis label
#' @param plot.cex cex graphic option
#' @param cols option to choose own colour palette
#' @export
jbplot_indices <- function(jabbainput, output.dir=getwd(),as.png=FALSE,width=5,height=4.5,plot.cex=0.9,xlab="Year",ylab="Abundance index",add=FALSE,add.legend=TRUE,legend.loc="topright",cols=NULL){
  
  years =   jabbainput$data$yr  
  abundance =jabbainput$settings$model.type
  dat = normIndex(jabbainput$data$cpue)
  se = as.matrix(sqrt(jabbainput$jagsdata$SE2)[1:length(years),])
  y = as.matrix(log(jabbainput$jagsdata$I)[1:length(years),])
  nI = ncol(y) 
  if(is.null(cols)) cols = jabbainput$settings$cols
  Par = list(mfrow=c(1,1),mar = c(4, 4, 1, 1), mgp =c(2.5,1,0),mai = c(0.6, 0.6, 0.1, 0.1),mex=0.8, tck = -0.02,cex=plot.cex)
  if(as.png==TRUE){png(file = paste0(output.dir,"/Index_",jabbainput$settings$assessment,"_",jabbainput$settings$scenario,".png"), width = width, height = height,
                       res = 200, units = "in")}
  if(add==FALSE) par(Par)
  
  Ylim = c(0, max(exp(log(jabbainput$jagsdata$I)+sqrt(jabbainput$jagsdata$SE2)*1.05)  ,na.rm =T))
  plot(years,years,type="n",xlim=c(min(years-1),max(years+2)),ylab=ylab,xlab=xlab,ylim=Ylim, frame = TRUE,xaxs="i",yaxs="i",xaxt="n")
  iv = c(-0.25,0.25)
  for(j in 1:nI){
    ds =runif(1,-0.2,0.2)
    for(t in 1:length(years)){
      lines(rep(years[t]+ds,2),c(exp(y[t,j]-1.96*se[t,j]),exp(y[t,j]+1.96*se[t,j])))
      lines(years[t]+ds+iv,rep(exp(y[t,j]-1.96*se[t,j]),2))  
      lines(years[t]+ds+iv,rep(exp(y[t,j]+1.96*se[t,j]),2))
    }  
    lines(years+ds,dat[,j+1],type="b",pch=21,cex=1.2,col=grey(0.4,0.7),bg=cols[j])
  }
  
  axis(1,at=seq(min(dat[,1]),max(dat[,1]),ceiling(length(dat[,1])/8)),tick=seq(min(dat[,1]),max(dat[,1]),ceiling(length(dat[,1])/8)),cex.axis=0.9)
  
  if(add.legend) legend(legend.loc,paste(names(dat)[2:(nI+1)]), lty = c(1, rep(nI)), lwd = c(rep(-1,nI)),pch=c(rep(21,nI)), pt.bg = c(cols[1:nI]), bty = "n", cex = 0.9,y.intersp = 0.8)
  if(as.png==TRUE) dev.off()
} # End of index plot




#' Plots Total Catch
#'
#' jbplot_catch()
#'
#' @param jabba output list from fit_jabba
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param add if true don't call par() to allow construction of multiplots
#' @param width plot width
#' @param height plot hight
#' @param verbose silent option
#' @export
jbplot_catch <- function(jabba,output.dir=getwd(),as.png = FALSE,add=FALSE, width = 5, height = 3.5,verbose=TRUE){
  if(verbose) cat(paste0("\n","><> jbplot_catch()  <><","\n"))
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  if(as.png==TRUE){png(file = paste0(output.dir,"/Landings_",jabba$assessment,"_",jabba$scenario,".png"), width = width, height = height,
      res = 200, units = "in")}
  if(add==FALSE) {par(Par)}

  cord.x <- c(jabba$yr,rev(jabba$yr))
  y<-rep(0,length(jabba$yr))
  plot(jabba$yr,(jabba$catch),type="l",ylim=c(0,max(jabba$catch,na.rm=T)), xaxs="i", yaxs="i",lty=1,lwd=1.3,xlab="Year",ylab=paste0("Catch ",jabba$settings$catch.metric),main="")
  polygon(cord.x,c(jabba$catch,rev(y)),col="gray",border=1,lty=1)
  if(as.png==TRUE){ dev.off()}
}


#' Plots estimated catch + CIs
#'
#' jbplot_catcherror() only works if add.catch.cv = TRUE in build_jabba
#'
#' @param jabba output list from fit_jabba
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param add if true don't call par() to allow construction of multiplots
#' @param width plot width
#' @param height plot height
#' @param verbose silent option
#' @export

jbplot_catcherror <- function(jabba,output.dir=getwd(),as.png = FALSE,add=FALSE, width = 5, height = 3.5,verbose=TRUE){
    if(jabba$settings$add.catch.CV==TRUE){
    if(verbose) cat(paste0("\n","><> jbplot_catcherror()  <><","\n"))
    Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.7,0), tck = -0.02,cex=0.8)
    if(as.png) { png(file = paste0(output.dir,"/Catch.fit_",jabba$assessment,"_",jabba$scenario,".png"), width = width, height = height,
        res = 200, units = "in")}
    if(add==FALSE){par(Par)}
    # estimated Catch
    predC = jabba$est.catch
    years = jabba$yr
    cord.x <- c(jabba$yr,rev(jabba$yr))
    cord.y<-c(predC[,3],rev(predC[,4]))
    plot(years,(jabba$catch),type="n",ylim=c(0,max(predC,na.rm=T)),lty=1,lwd=1.3,xlab="Year",ylab=paste0("Catch ",jabba$settings$catch.metric),main="")
    polygon(cord.x,cord.y,col="gray",border=0,lty=1)
    lines(years,predC[,2],lwd=2,col=4)
    points(years,(jabba$catch),pch=21,bg=0,cex=1.5)
    legend("topright",c("Observed","Predicted"),pch=c(21,-1),bg=0,lwd=c(-1,2),col=c(1,4),bty="n")
    if(as.png==TRUE & add==FALSE){dev.off()}
  } else {
   if(verbose) cat(paste0("\n","><> jbplot_catcherror() only available if add.catch.CV=TRUE <><","\n"))
  }
}

#' Plot of prior and posterior distributions
#'
#' Prior and posterior distributions of esimated parameters: K, r, psi (depletion) and variances
#'
#' @param jabba output list from fit_jabba
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param add if true don't call par() to allow construction of multiplots
#' @param width plot width
#' @param height plot hight
#' @param mfrow set up plot frame
#' @param addPP show PPMR and PPVR
#' @param verbose silent option  
#' @export

jbplot_ppdist <- function(jabba, output.dir=getwd(),as.png = FALSE,mfrow=c(round((ncol(jabba$pars_posterior))/3+0.33,0),3),width  = 8, height = 2.5*round(ncol(jabba$pars_posterior)/3,0),cex=0.8,verbose=TRUE,addPP=TRUE){
  if(verbose) cat(paste0("\n","><> jbplot_ppist() - prior and posterior distributions  <><","\n"))
  out =   jabba$pars_posterior
  node_id = names(out)
  #informative priors
  Prs = as.matrix(cbind(jabba$settings$K.pr,jabba$settings$r.pr,c(0,0),jabba$settings$psi.pr))

  #Posteriors
  Par = list(mfrow=mfrow,mai=c(0.4,0.1,0,.1),omi = c(0.3,0.5,0.1,0) + 0.1,mgp=c(1,0.1,0), tck = -0.02,cex=0.8)
  if(as.png){png(file = paste0(output.dir,"/Posteriors_",jabba$assessment,"_",jabba$scenario,".png"),width  = width, height = height,
      res = 200, units = "in")}
  par(Par)

  for(i in 1:length(node_id))
  {

    post.par = as.numeric(unlist(out[paste(node_id[i])]))

    if(i==1){

      rpr =  rlnorm(10000,log(Prs[1,i]),Prs[2,i])
      pdf = stats::density(post.par,adjust=2)
      prior = dlnorm(sort(rpr),log(Prs[1,i]),Prs[2,i])
      plot(pdf,type="l",ylim=range(prior,pdf$y),xlim=range(c(pdf$x,quantile(rpr,c(0.0001,0.95)))),yaxt="n",xlab="K",ylab="",xaxs="i",yaxs="i",main="")

      polygon(c(sort(rpr),rev(sort(rpr))),c(prior,rep(0,length(sort(rpr)))),col=gray(0.4,1))
      polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
      legend('right',c("Prior","Posterior"),pch=22,cex=cex+0.1,pt.cex=cex+0.1,pt.bg = c(grey(0.4,1),grey(0.8,0.6)),bty="n")
      PPVR = round((sd(post.par)/mean(post.par))^2/(sd(rpr)/mean(rpr))^2,3)
      PPVM = round(mean(post.par)/mean(rpr),3)
      legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=cex,bty="n", x.intersp=0.5,y.intersp=0.85)

    }
    if(i==2){
      rpr = rlnorm(10000,log(Prs[1,i]),Prs[2,i])
      pdf = stats::density(post.par,adjust=2)
      prior = dlnorm(sort(rpr),log(Prs[1,i]),Prs[2,i])
      plot(pdf$x,pdf$y,type="l",ylim=range(prior,pdf$y),xlim=range(c(post.par,quantile(rpr,c(0.0001,0.95)))),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i")
      polygon(c(sort(rpr),rev(sort(rpr))),c(prior,rep(0,length(sort(rpr)))),col=gray(0.4,1))
      polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
      PPVR = round((sd(post.par)/mean(post.par))^2/(sd(rpr)/mean(rpr))^2,3)
      PPVM = round(mean(post.par)/mean(rpr),3)
     if(addPP) legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=cex,bty="n", x.intersp=0.5,y.intersp=0.85)


    }

    if(i==3){
      if(jabba$settings$model.id<4){
        plot(1,1,type="n",xlim=range(0.5,2.5),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i")
        abline(v=jabba$pars["m",1],lwd=2)}
      if(jabba$settings$model.id==4){
        mpr = rlnorm(10000,log(jabba$settings$mu.m),jabba$settings$m.CV)
        pdf = stats::density(post.par,adjust=2)
        prior = dlnorm(sort(mpr),log(jabba$settings$mu.m),jabba$settings$m.CV)
        plot(pdf$x,pdf$y,type="l",ylim=range(prior,pdf$y),xlim=range(c(post.par,quantile(rpr,c(0.0001,0.95)))),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i")

        polygon(c(sort(mpr),rev(sort(mpr))),c(prior,rep(0,length(sort(rpr)))),col=gray(0.4,1))
        polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
        PPVR = round((sd(post.par)/mean(post.par))^2/(sd(mpr)/mean(mpr))^2,3)
        PPVM = round(mean(post.par)/mean(mpr),3)
       if(addPP) legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=cex,bty="n", x.intersp=0.5,y.intersp=0.85)

      }
    }


    if(i==4){
      if(jabba$settings$psi.dist=="beta"){
        rpr = rbeta(10000,(Prs[1,4]),Prs[2,4])
        pdf = stats::density(post.par,adjust=2)
        prior = dbeta(sort(rpr),jabba$settings$psi.pr[1],jabba$settings$psi.pr[2])
        PPVR = round((sd(post.par)/mean(post.par))^2/(sd(rpr)/mean(rpr))^2,3)
        PPVM = round(mean(post.par)/mean(rpr),3)
        if(addPP) legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=cex,bty="n", x.intersp=0.5,y.intersp=0.85)

      } else {
        rpr = rlnorm(10000,log(Prs[1,4]),Prs[2,4])
        pdf = stats::density(post.par,adjust=2)
        prior = dlnorm(sort(rpr),log(Prs[1,4]),Prs[2,4])}
      plot(pdf,type="l",ylim=range(quantile(c(prior,pdf$y,c(0,0.95)))),xlim=range(c(0.5,post.par,pdf$x,quantile(rpr,c(0.001,0.999)))),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i",main="")
      polygon(c(sort(rpr),rev(sort(rpr))),c(prior,rep(0,length(sort(rpr)))),col=gray(0.4,1))
      polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
      PPVR = round((sd(post.par)/mean(post.par))^2/(sd(rpr)/mean(rpr))^2,3)
      PPVM = round(mean(post.par)/mean(rpr),3)
      if(addPP) legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=cex,bty="n", x.intersp=0.5,y.intersp=0.85)

      #legend('topright',c("Prior","Posterior"),pch=22,pt.cex=1.5,pt.bg = c(grey(0.4,1),grey(0.8,0.6)),bty="n")
    }

    if(i>4){
      if(jabba$settings$sigma.proc!=TRUE & i==length(node_id)) {
        plot(1,1,type="n",xlim=range(0,0.15^2),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i")
        abline(v=jabba$settings$sigma.proc^2,lwd=2)} else {

          pdf = stats::density(post.par,adjust=2)
          plot(pdf,type="l",xlim=range(0,post.par),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i",main="")
          if(i==length(node_id)& jabba$settings$igamma[1]>0.9){
            rpr = 1/rgamma(10000,jabba$settings$igamma[1],jabba$settings$igamma[2])
            prior = stats::density(rpr,adjust=2)
            polygon(c(prior$x,rev(prior$x)),c(prior$y,rep(0,length(prior$y))),col=gray(0.4,1))
            PPVR = round((sd(post.par)/mean(post.par))^2/(sd(rpr)/mean(rpr))^2,3)
            PPVM = round(mean(post.par)/mean(rpr),3)
           if(addPP) legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=cex,bty="n", x.intersp=0.5,y.intersp=0.85)

          }

          polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
          #legend('topright',c("Posterior"),pch=22,pt.cex=1.5,pt.bg = c(grey(0.8,0.6)),bty="n")
        } }

  }
  mtext(paste("Density"), side=2, outer=TRUE, at=0.5,line=1,cex=0.9)
  if(as.png){dev.off()}
} # End of ppdist plot


#' Plot mcmc chains
#'
#' MCMC chains of esimated parameters: K, r, m (shape), psi (depletion) and variances
#'
#' @param jabba output list from fit_jabba
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param add if true don't call par() to allow construction of multiplots
#' @param width plot width
#' @param height plot hight
#' @param mfrow set up plot frame  
#' @param verbose silent option
#' @export
jbplot_mcmc <- function(jabba, output.dir=getwd(),as.png = FALSE,mfrow=c(round((ncol(jabba$pars_posterior))/3+0.33,0),3),width  = 8, height = 2.5*round(jabba$pars_posterior/3,0),verbose=TRUE){

 if(verbose) cat(paste0("\n","><> jbplot_mcmc() - mcmc chains  <><","\n"))
  out =   jabba$pars_posterior
  node_id = names(out)

  Par = list(mfrow=c(round(length(node_id)/3+0.33,0),3),mai=c(0.4,0.1,0,.1),omi = c(0.3,0.5,0.1,0) + 0.1,mgp=c(1,0.1,0), tck = -0.02,cex=0.8)
  if(as.png==TRUE){png(file = paste0(output.dir,"/MCMC_",jabba$assessment,"_",jabba$scenario,".png"), width = 8, height = 2.5*round(length(node_id)/3,0),
      res = 200, units = "in")}
  par(Par)
  for(i in 1:length(node_id)){

    post.par = as.numeric(unlist(out[paste(node_id[i])]))
    plot(out[,i],xlab=paste(node_id[i]),ylab="",type="l",col=4)
    lines(rep(mean(out[,i]),length(out[,i])),col=2,lwd=2)
  }
  if(as.png==TRUE){dev.off()}
}

#' Plot of  fitted CPUE indices
#'
#' Plots observed and fitted cpue indices with expexted CIs (dark grey) and posterior predictive distribution (light grey)
#'
#' @param jabba output list from fit_jabba
#' @param index option to plot specific indices (numeric & in order)
#' @param output.dir directory to save plots
#' @param add if TRUE is surpresses par() within the plotting function
#' @param as.png save as png file of TRUE
#' @param single.plots if TRUE plot invidual fits else make multiplot
#' @param plotCIs plots error bars on CPUE observations
#' @param width plot width
#' @param height plot height
#' @param verbose silent option
#' @export
jbplot_cpuefits <- function(jabba,index=NULL, output.dir=getwd(),add=FALSE,as.png=FALSE,single.plots=add,width=NULL,height=NULL,plotCIs=TRUE,verbose=TRUE){
  if(as.png==TRUE) add=FALSE
  if(jabba$settings$CatchOnly==FALSE | jabba$settings$Auxiliary==TRUE){
    if(verbose) cat(paste0("\n","><> jbplot_cpue() - fits to CPUE <><","\n"))
    
    
    if(jabba$settings$Auxiliary){
      if(jabba$settings$CatchOnly){
        jabba$settings$nI = jabba$settings$nA
        jabba$settings$I = as.matrix(jabba$settings$A)
        jabba$settings$SE2 = as.matrix(jabba$settings$A.SE2)
        di = dim(jabba$cpue.ppd)[3]
        jabba$cpue.ppd[,,1:(di-1)] =  jabba$cpue.ppd[,,2:(di)] 
        jabba$cpue.hat[,,1:(di-1)] =  jabba$cpue.hat[,,2:(di)] 
        
      } else {
        jabba$settings$nI = jabba$settings$nI+jabba$settings$nA
        jabba$settings$I = cbind(jabba$settings$I,as.matrix(jabba$settings$A))
        jabba$settings$SE2 = cbind(jabba$settings$SE2,as.matrix(jabba$settings$A.SE2))
      }
    }
      if(is.null(index)) index = 1:jabba$settings$nI
    
    
    N = jabba$settings$N
    years= jabba$yr
    indices = unique(jabba$diags$name)[index]
    CPUE = jabba$settings$I
  
    n.indices = length(indices)
    series = 1:n.indices
    #check.yrs = abs(apply(jabba$residuals,2,sum,na.rm=TRUE))
    #cpue.yrs = years[check.yrs>0]
    
    if(single.plots==TRUE){
    if(is.null(width)) width = 5
    if(is.null(height)) height = 3.5
    for(i in 1:n.indices){
      Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
      if(as.png==TRUE){png(file = paste0(output.dir,"/Fits",jabba$assessment,"_",jabba$scenario,"_",indices[i],".png"), width = width, height = height,
                           res = 200, units = "in")}
      if(!add){
      if(as.png==TRUE | i==1) par(Par)
      }
      # set observed vs predicted CPUE
      Yr = jabba$yr
      Yr = min(Yr):max(Yr)
      yr = Yr-min(years)+1

      fit =  t(jabba$cpue.ppd[,c(2,1,3),index[i]])
      fit.hat = t(jabba$cpue.hat[,c(2,1,3),index[i]])
      mufit = mean(fit[2,])
      fit = fit#/mufit
      fit.hat = fit.hat#/mufit

      cpue.i = CPUE[is.na(CPUE[,index[i]])==F,index[i]]
      yr.i = Yr[is.na(CPUE[,index[i]])==F]
      se.i = sqrt(jabba$settings$SE2[is.na(CPUE[,index[i]])==F,index[i]])

      ylim = c(min(fit*0.9,exp(log(cpue.i)-1.96*se.i)), max(fit*1.05,exp(log(cpue.i)+1.96*se.i)))

      cord.x <- c(Yr,rev(Yr))
      cord.y <- c(fit[1,yr],rev(fit[3,yr]))
      cord.yhat <- c(fit.hat[1,yr],rev(fit.hat[3,yr]))
      # Plot Observed vs predicted CPUE
      plot(years,CPUE[,index[i]],ylab="Index",xlab="Year",ylim=ylim,xlim=range(jabba$yr),type='n',xaxt="n",yaxt="n")
      axis(1,labels=TRUE,cex=0.8)
      axis(2,labels=TRUE,cex=0.8)
      polygon(cord.x,cord.y,col=grey(0.5,0.5),border=0,lty=2)
      polygon(cord.x,cord.yhat,col=grey(0.3,0.5),border=grey(0.3,0.5),lty=2)

      lines(Yr,fit[2,yr],lwd=2,col=1)
      if(plotCIs){
        upper <- qlnorm(.975,meanlog=log(cpue.i),sdlog=se.i)
        lower <- qlnorm(.025,meanlog=log(cpue.i),sdlog=se.i)
        arrows(x0=yr.i, y0=lower,
               x1=yr.i, y1=upper,
               length=0.02, angle=90, code=3, col=1)
      } 
      
      
      
      points(yr.i,cpue.i,pch=21,xaxt="n",yaxt="n",bg="white")#}
      legend('top',paste(indices[i]),bty="n",y.intersp = -0.2,cex=0.9)
      #mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
      #mtext(paste("Normalized Index"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
      if(as.png==TRUE) dev.off()
    }} else {

    if(is.null(width)) width = 7
    if(is.null(height)) height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0)
    Par = list(mfrow=c(round(n.indices/2+0.01,0),ifelse(n.indices==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
    if(as.png==TRUE){png(file = paste0(output.dir,"/Fits_",jabba$assessment,"_",jabba$scenario,".png"), width = 7, height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0),
                           res = 200, units = "in")}
    if(!add) par(Par)

    for(i in 1:n.indices){
      # set observed vs predicted CPUE
      Yr = jabba$yr
      Yr = min(Yr):max(Yr)
      yr = Yr-min(years)+1
      
      fit =  t(jabba$cpue.ppd[,c(2,1,3),index[i]])
      fit.hat = t(jabba$cpue.hat[,c(2,1,3),index[i]])
      mufit = mean(fit[2,])
      fit = fit#/mufit
      fit.hat = fit.hat#/mufit
      
      cpue.i = CPUE[is.na(CPUE[,index[i]])==F,index[i]]
      yr.i = Yr[is.na(CPUE[,index[i]])==F]
      se.i = sqrt(jabba$settings$SE2[is.na(CPUE[,index[i]])==F,index[i]])
      
      
      
      ylim = c(min(fit*0.9,exp(log(cpue.i)-1.96*se.i)), max(fit*1.05,exp(log(cpue.i)+1.96*se.i)))
      
      cord.x <- c(Yr,rev(Yr))
      cord.y <- c(fit[1,yr],rev(fit[3,yr]))
      cord.yhat <- c(fit.hat[1,yr],rev(fit.hat[3,yr]))
      # Plot Observed vs predicted CPUE
      plot(years,CPUE[,index[i]],ylab="",xlab="",ylim=ylim,xlim=range(jabba$yr),type='n',xaxt="n",yaxt="n")
      axis(1,labels=TRUE,cex=0.8)
      axis(2,labels=TRUE,cex=0.8)
      polygon(cord.x,cord.y,col=grey(0.5,0.5),border=0,lty=2)
      polygon(cord.x,cord.yhat,col=grey(0.3,0.5),border=grey(0.3,0.5),lty=2)

      lines(Yr,fit[2,yr],lwd=2,col=1)
      #if(jabba$settings$SE.I  ==TRUE | max(jabba$settings$SE2)>0.01){ gplots::plotCI(yr.i,cpue.i/mufit,ui=exp(log(cpue.i)+1.96*se.i)/mufit,li=exp(log(cpue.i)-1.96*se.i)/mufit,add=T,gap=0,pch=21,xaxt="n",yaxt="n",pt.bg = "white")}else{
      if(plotCIs){
        upper <- qlnorm(.975,meanlog=log(cpue.i),sdlog=se.i)
        lower <- qlnorm(.025,meanlog=log(cpue.i),sdlog=se.i)
         arrows(x0=yr.i, y0=lower,
               x1=yr.i, y1=upper,
               length=0.02, angle=90, code=3, col=1)
        } 
      
       points(yr.i,cpue.i,pch=21,xaxt="n",yaxt="n",bg="white")#}

      legend('top',paste(indices[i]),bty="n",y.intersp = -0.2,cex=0.9)
    }
      mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
      mtext(paste("Index"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
      if(as.png==TRUE){dev.off()}
    }
    } else {
      if(verbose) cat(paste0("\n","><> jbplot_cpuefits() not available CatchOnly=TRUE <><","\n"))
    }
} # End of CPUE plot function

#' log(CPUE) fits
#'
#' Plot of fitted CPUE indices on log-scale (r4ss-style)
#'
#' @param jabba output list from fit_jabba
#' @param index option to plot specific indices (numeric & in order)
#' @param output.dir directory to save plots
#' @param add if TRUE is surpresses par() within the plotting function
#' @param as.png save as png file of TRUE
#' @param single.plots if TRUE plot invidual fits else make multiplot
#' @param plotCIs plots error bars on CPUE observations
#' @param width plot width
#' @param height plot height
#' @param verbose silent option
#' @export

jbplot_logfits <- function(jabba,index=NULL, output.dir=getwd(),add=FALSE,as.png=FALSE,single.plots=add,width=NULL,height=NULL,plotCIs=TRUE,verbose=TRUE){
  if(jabba$settings$CatchOnly==FALSE | jabba$settings$Auxiliary==TRUE){
    if(verbose) cat(paste0("\n","><> jbplot_cpue() - fits to CPUE <><","\n"))
    
    if(jabba$settings$Auxiliary){
      if(jabba$settings$CatchOnly){
        jabba$settings$nI = jabba$settings$nA
        jabba$settings$I = as.matrix(jabba$settings$A)
        jabba$settings$SE2 = as.matrix(jabba$settings$A.SE2)
        di = dim(jabba$cpue.ppd)[3]
        jabba$cpue.ppd[,,1:(di-1)] =  jabba$cpue.ppd[,,2:(di)] 
        jabba$cpue.hat[,,1:(di-1)] =  jabba$cpue.hat[,,2:(di)] 
        
        
      } else {
        jabba$settings$nI = jabba$settings$nI+jabba$settings$nA
        jabba$settings$I = cbind(jabba$settings$I,as.matrix(jabba$settings$A))
        jabba$settings$SE2 = cbind(jabba$settings$SE2,as.matrix(jabba$settings$A.SE2))
      }
    }
    if(is.null(index)) index = 1:jabba$settings$nI
    
    N = jabba$settings$N
    years= jabba$yr
    indices = unique(jabba$diags$name)[index]
    CPUE = jabba$settings$I
    n.indices = length(indices)
    series = 1:n.indices
    check.yrs = abs(apply(jabba$residuals,2,sum,na.rm=TRUE))
    cpue.yrs = years[check.yrs>0]

    if(single.plots==TRUE){
      if(is.null(width)) width = 5
      if(is.null(height)) height = 3.5
      for(i in 1:n.indices){
        Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
        if(as.png==TRUE){png(file = paste0(output.dir,"/logFits",jabba$assessment,"_",jabba$scenario,"_",indices[i],".png"), width = width, height = height,
                             res = 200, units = "in")}
        if(!add){
        if(as.png==TRUE | i==1) par(Par)
        }
      Yr = jabba$yr
      Yr = min(Yr):max(Yr)
      yr = Yr-min(years)+1

      fit = t(jabba$cpue.hat[,c(2,1,3),index[i]])
      mufit = mean(fit[2,])
      fit = fit/mufit
      cpue.i = CPUE[is.na(CPUE[,index[i]])==F,index[i]]
      yr.i = Yr[is.na(CPUE[,index[i]])==F]
      se.i = sqrt(jabba$settings$SE2[is.na(CPUE[,index[i]])==F,index[i]])

      ylim = log(c(min(fit[,yr[is.na(CPUE[,index[i]])==F]]*0.8,exp(log(cpue.i)-1.96*se.i)/mufit), max(fit[,yr[is.na(CPUE[,index[i]])==F]]*1.3,exp(log(cpue.i)+1.96*se.i)/mufit)))

      # Plot Observed vs predicted CPUE
      plot(years,CPUE[,index[i]],ylab="log Index",xlab="Year",ylim=ylim,xlim=range(yr.i),type='n',xaxt="n",yaxt="n")
      axis(1,labels=TRUE,cex=0.8)
      axis(2,labels=TRUE,cex=0.8)

      lines(Yr,log(fit[2,yr]),lwd=2,col=4)
      
      if(plotCIs){
        upper <- qnorm(.975,mean=log(cpue.i/mufit),sd=se.i)
        lower <- qnorm(.025,mean=log(cpue.i/mufit),sd=se.i)
        arrows(x0=yr.i, y0=lower,
               x1=yr.i, y1=upper,
               length=0.02, angle=90, code=3, col=1)
      } 
      
      points(yr.i,log(cpue.i/mufit),pch=21,xaxt="n",yaxt="n",bg="white")#}
      legend('topright',paste(indices[i]),bty="n",y.intersp = -0.2,cex=0.8)
      #mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
      #mtext(paste("Normalized Index"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
      if(as.png==TRUE){dev.off()}
      }
      } else { # single.plots = F
        if(is.null(width)) width = 7
        if(is.null(height)) height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0)
        Par = list(mfrow=c(round(n.indices/2+0.01,0),ifelse(n.indices==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
        if(as.png==TRUE){png(file = paste0(output.dir,"/logFits_",jabba$assessment,"_",jabba$scenario,".png"), width = 7, height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0),
                             res = 200, units = "in")}
        if(!add) par(Par)
        for(i in 1:n.indices){
          Yr = jabba$yr
          Yr = min(Yr):max(Yr)
          yr = Yr-min(years)+1
          
          fit = t(jabba$cpue.hat[,c(2,1,3),index[i]])
          mufit = mean(fit[2,])
          fit = fit/mufit
          cpue.i = CPUE[is.na(CPUE[,index[i]])==F,index[i]]
          yr.i = Yr[is.na(CPUE[,index[i]])==F]
          se.i = sqrt(jabba$settings$SE2[is.na(CPUE[,index[i]])==F,index[i]])
          
          ylim = log(c(min(fit[,yr[is.na(CPUE[,index[i]])==F]]*0.8,exp(log(cpue.i)-1.96*se.i)/mufit), max(fit[,yr[is.na(CPUE[,index[i]])==F]]*1.3,exp(log(cpue.i)+1.96*se.i)/mufit)))
          
          # Plot Observed vs predicted CPUE
          plot(years,CPUE[,index[i]],ylab="",xlab="",ylim=ylim,xlim=range(yr.i),type='n',xaxt="n",yaxt="n")
          axis(1,labels=TRUE,cex=0.8)
          axis(2,labels=TRUE,cex=0.8)

          lines(Yr,log(fit[2,yr]),lwd=2,col=4)
          
          if(plotCIs){
            upper <- qnorm(.975,mean=log(cpue.i/mufit),sd=se.i)
            lower <- qnorm(.025,mean=log(cpue.i/mufit),sd=se.i)
            arrows(x0=yr.i, y0=lower,
                   x1=yr.i, y1=upper,
                   length=0.02, angle=90, code=3, col=1)
          } 
          
            points(yr.i,log(cpue.i/mufit),pch=21,xaxt="n",yaxt="n",bg="white")#}
          legend('topright',paste(indices[i]),bty="n",y.intersp = -0.2,cex=0.8)
        }
        mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
        mtext(paste("Log Index"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
        if(as.png==TRUE){dev.off()}

        }
        }  else {
    if(verbose) cat(paste0("\n","><> jbplot_logfit() not available CatchOnly=TRUE <><","\n"))
  }
} # End of logfit


#' JABBA residual plot
#'
#' plots residuals for all indices as boxplot with a loess showing systematic trends
#'
#' @param jabba output list from fit_jabba
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param add if true don't call par() to allow construction of multiplots
#' @param width plot width
#' @param height plot hight
#' @param cols option to add colour palette 
#' @param verbose silent option
#' @export
jbplot_residuals <- function(jabba,output.dir=getwd(),as.png = FALSE,add=FALSE, width = 5, height = 3.5,cols=NULL,verbose=TRUE){

  if(jabba$settings$CatchOnly==FALSE | jabba$settings$Auxiliary==TRUE){
    if(verbose) cat(paste0("\n","><> jbplot_cpue() - fits to CPUE <><","\n"))
    
    if(jabba$settings$Auxiliary){
      if(jabba$settings$CatchOnly){
        jabba$settings$nI = jabba$settings$nA
        jabba$settings$I = as.matrix(jabba$settings$A)
        jabba$settings$SE2 = as.matrix(jabba$settings$A.SE2)
      } else {
        jabba$settings$nI = jabba$settings$nI+jabba$settings$nA
        jabba$settings$I = cbind(jabba$settings$I,as.matrix(jabba$settings$A))
        jabba$settings$SE2 = cbind(jabba$settings$SE2,as.matrix(jabba$settings$A.SE2))
      }
    }
  
    if(is.null(cols)) cols = jabba$settings$cols 
    years = jabba$yr
    check.yrs = abs(apply(jabba$residuals,2,sum,na.rm=TRUE))
    cpue.yrs = years[check.yrs>0]
    Resids = jabba$residuals
    Yr = jabba$yr
    n.years = length(Yr)
    n.indices = jabba$settings$nI
    indices = unique(jabba$diags$name)
    series = 1:jabba$settings$nI

    # JABBA-residual plot
    Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
    if(as.png==TRUE){png(file = paste0(output.dir,"/Residuals_",jabba$assessment,"_",jabba$scenario,".png"), width = width, height = height,
        res = 200, units = "in")}
    if(add==FALSE) par(Par)


    plot(Yr,Yr,type = "n",ylim=ifelse(rep(max(Resids,na.rm = T),2)>0.9,range(1.2*Resids,na.rm = T),range(c(-1.3,1.2))),xlim=range(cpue.yrs),ylab="log residuals",xlab="Year")
    boxplot(Resids,add=TRUE,at=c(Yr),xaxt="n",col=grey(0.8,0.5),notch=FALSE,outline = FALSE)
    abline(h=0,lty=2)
    positions=runif(n.indices,-0.2,0.2)

    for(i in 1:n.indices){
      for(t in 1:n.years){
        lines(rep((Yr+positions[i])[t],2),c(0,Resids[i,t]),col=cols[i])}
      points(Yr+positions[i],Resids[i,],col=1,pch=21,bg=cols[i])}
    mean.res = apply(Resids,2,mean,na.rm =TRUE)[as.numeric(colnames(Resids))%in%cpue.yrs]
    smooth.res = predict(loess(mean.res~cpue.yrs),data.frame(cpue.yrs=cpue.yrs))
    lines(cpue.yrs,smooth.res,lwd=2)
    # get degree of freedom
    Nobs =length(as.numeric(Resids)[is.na(as.numeric(Resids))==FALSE])

    RMSE = round(jabba$stats[5,2],1)

    legend('topright',c(paste0("RMSE = ",RMSE,"%")),bty="n")
    legend('bottomright',c(paste(indices),"Loess"),bty="n",col=1,pt.cex=1.1,cex=0.75,pch=c(rep(21,n.indices),-1),pt.bg=c(jabba$settings$col[series],1),lwd=c(rep(-1,n.indices),2))
    if(as.png==TRUE){dev.off()}
  }  else {
    if(verbose) cat(paste0("\n","><> jbplot_residuals() not available CatchOnly=TRUE <><","\n"))
  }
} # End of functions

#' JABBA standardized residual plot
#'
#' plots standardized residuals for all indices as boxplot with a loess showing systematic trends
#'
#' @param jabba output list from fit_jabba
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param add if true don't call par() to allow construction of multiplots
#' @param width plot width
#' @param height plot hight
#' @param cols option to add colour palette 
#' @param verbose silent option
#' @export
jbplot_stdresiduals <- function(jabba, output.dir=getwd(),as.png=FALSE,add=FALSE,width = 5, height= 3.5,cols=NULL,verbose=TRUE){

  if(jabba$settings$CatchOnly==FALSE){
    if(verbose) cat(paste0("\n","><> jbplot_staresiduals() - standardized residuals  <><","\n"))
    if(is.null(cols)) cols = jabba$settings$cols 
    years = jabba$yr
    check.yrs = abs(apply(jabba$residuals,2,sum,na.rm=TRUE))
    cpue.yrs = years[check.yrs>0]
    Resids = jabba$residuals
    Yr = jabba$yr
    n.years = length(Yr)
    n.indices = jabba$settings$nI
    indices = unique(jabba$diags$name)
    series = 1:jabba$settings$nI
    StResid = jabba$std.residuals

    Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
    if(as.png==TRUE){png(file = paste0(output.dir,"/StandardizedResids_",jabba$assessment,"_",jabba$scenario,".png"), width = width, height = height,
        res = 200, units = "in")}
    if(add==FALSE) par(Par)
    # Standardized Residuals
    plot(Yr,Yr,type = "n",ylim=c(min(-1,-1.2*max(abs(StResid),na.rm = T)),max(1,1.2*max(abs(StResid),na.rm = T))),xlim=range(cpue.yrs),ylab="Standardized residuals",xlab="Year")
    boxplot(StResid,add=TRUE,at=c(Yr),xaxt="n",col=grey(0.8,0.5),notch=FALSE,outline = FALSE)
    abline(h=0,lty=2)
    positions=runif(n.indices,-0.2,0.2)

    for(i in 1:n.indices){
      for(t in 1:n.years){
        lines(rep((Yr+positions[i])[t],2),c(0,StResid[i,t]),col=cols[i])}
      points(Yr+positions[i],StResid[i,],col=1,pch=21,bg=cols[i])}
    mean.res = apply(StResid,2,mean,na.rm =TRUE)[as.numeric(colnames(Resids))%in%cpue.yrs]
    smooth.res = predict(loess(mean.res~cpue.yrs),data.frame(cpue.yrs=cpue.yrs))
    lines(cpue.yrs,smooth.res,lwd=2)
    SDNR = round(sqrt(sum(StResid^2,na.rm =TRUE)/(jabba$stats[1,2]-1)),2)
    Crit.value = (qchisq(.95, df=(jabba$stats[1,2]-1))/(jabba$stats[1,2]-1))^0.5
    legend('topright',c(paste0("SDNR = ",SDNR,"(",round(Crit.value,2),")")),bty="n")
    legend('bottomright',c(paste(indices),"Loess"),bty="n",col=1,cex=0.75,pt.cex=1.1,pch=c(rep(21,n.indices),-1),pt.bg=c(jabba$settings$cols[series],1),lwd=c(rep(-1,n.indices),2))
    if(as.png==TRUE) dev.off()
  }  else {
   if(verbose) cat(paste0("\n","><> jbplot_stdresiduals() not available CatchOnly=TRUE <><","\n"))
  }

} # end of function

#' JABBA runs test plots
#'
#' Residual diagnostics with runs test p-value and 3xsigma limits
#' @param jabba output list from fit_jabba
#' @param mixing c("less","greater","two.sided"). Default "less" is checking for positive autocorrelation only
#' @param index option to plot specific indices (numeric & in order)
#' @param output.dir directory to save plots
#' @param add if true par() is surpressed within the plot function
#' @param as.png save as png file of TRUE
#' @param single.plots if TRUE plot invidual fits else make multiplot
#' @param width plot width
#' @param height plot hight
#' @param verbose silent option
#' @export
#' @examples 
#' data(iccat)
#' bet= iccat$bet
#' jb = build_jabba(catch=bet$catch,cpue=bet$cpue,se=bet$se,assessment="BET",scenario = "Ref",model.type = "Pella",igamma = c(0.001,0.001),verbose=FALSE)
#' fit = fit_jabba(jb,quickmcmc=TRUE,verbose=FALSE)
#' jbrunstest(fit)
#' jbrunstest(fit,index=2)
#' jbplot_runstest(fit,verbose=FALSE)

jbplot_runstest <- function(jabba,index=NULL,mixing="less", output.dir=getwd(),add=FALSE,as.png=FALSE,single.plots=add,width=NULL,height=NULL,verbose=TRUE){
  if(as.png==TRUE) add=FALSE
  if(jabba$settings$CatchOnly==FALSE | jabba$settings$Auxiliary==TRUE){
   if(verbose) cat(paste0("\n","><> jbplot_runstest()   <><","\n"))

    all.indices = unique(jabba$diags$name)
    if(is.null(index)) index = 1:length(all.indices)
    indices = unique(jabba$diags$name)[index]
    n.indices = length(indices)
    Resids = jabba$residuals
    years = jabba$yr
    check.yrs = abs(apply(jabba$residuals,2,sum,na.rm=TRUE))
    cpue.yrs = years[check.yrs>0]
    

    if(single.plots==TRUE){
        if(is.null(width)) width = 5
        if(is.null(height)) height = 3.5
        for(i in 1:n.indices){
          Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
          if(as.png==TRUE){png(file = paste0(output.dir,"/ResRunsTests_",jabba$assessment,"_",jabba$scenario,"_",indices[i],".png"), width = width, height = height,
                               res = 200, units = "in")}

          if(add==FALSE){
          if(as.png==TRUE | i==1) par(Par)
          }


      resid = (Resids[index[i],is.na(Resids[index[i],])==F])
      res.yr = years[is.na(Resids[index[i],])==F]
      runstest = jbruns_sig3(x=as.numeric(resid),type="resid",mixing=mixing)
      # CPUE Residuals with runs test
      plot(res.yr,rep(0,length(res.yr)),type="n",ylim=c(min(-1,runstest$sig3lim[1]*1.25),max(1,runstest$sig3lim[2]*1.25)),lty=1,lwd=1.3,xlab="Year",ylab="Residuals")
      abline(h=0,lty=2)
      lims = runstest$sig3lim
      cols =  c(rgb(1,0,0,0.5),rgb(0,1,0,0.5))[ifelse(runstest$p.runs<0.05,1,2)]
      rect(min(years-1),lims[1],max(years+1),lims[2],col=cols,border=cols) # only show runs if RMSE >= 0.1
      for(j in 1:length(resid)){
        lines(c(res.yr[j],res.yr[j]),c(0,resid[j]))
      }
      points(res.yr,resid,pch=21,bg=ifelse(resid < lims[1] | resid > lims[2],2,"white"),cex=1)
      legend('topright',paste(indices[i]),bty="n",y.intersp = -0.2,cex=0.8)
      #mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
      #mtext(expression(log(cpue[obs])-log(cpue[pred])), side=2, outer=TRUE, at=0.5,line=1,cex=1)
      if(as.png==TRUE){dev.off()}
        } # end of loop
    } else { # single.plot = FALSE
    if(is.null(width)) width = 7
    if(is.null(height)) height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0)
    Par = list(mfrow=c(round(n.indices/2+0.01,0),ifelse(n.indices==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
    if(as.png==TRUE){png(file = paste0(output.dir,"/ResRunsTests_",jabba$assessment,"_",jabba$scenario,".png"), width = 7, height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0),
                         res = 200, units = "in")}
    if(add==FALSE) par(Par)
    for(i in 1:n.indices){
      resid = (Resids[index[i],is.na(Resids[index[i],])==F])
      res.yr = years[is.na(Resids[index[i],])==F]
      runstest = jbruns_sig3(x=as.numeric(resid),type="resid",mixing=mixing)
      # CPUE Residuals with runs test
      plot(res.yr,rep(0,length(res.yr)),type="n",ylim=c(min(-1,runstest$sig3lim[1]*1.25),max(1,runstest$sig3lim[2]*1.25)),lty=1,lwd=1.3,xlab="Year",ylab="Residuals")
      abline(h=0,lty=2)
      lims = runstest$sig3lim
      cols =  c(rgb(1,0,0,0.5),rgb(0,1,0,0.5))[ifelse(runstest$p.runs<0.05,1,2)]
      rect(min(years-1),lims[1],max(years+1),lims[2],col=cols,border=cols) # only show runs if RMSE >= 0.1
      for(j in 1:length(resid)){
        lines(c(res.yr[j],res.yr[j]),c(0,resid[j]))
      }
      points(res.yr,resid,pch=21,bg=ifelse(resid < lims[1] | resid > lims[2],2,"white"),cex=1)
      legend('topright',paste(indices[i]),bty="n",y.intersp = -0.2,cex=0.8)

    }
    mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
    mtext("Residuals", side=2, outer=TRUE, at=0.5,line=1,cex=1)
    if(as.png==TRUE){dev.off()}
  }


  }else {
    if(verbose) cat(paste0("\n","><> jbplot_runstest() not available CatchOnly=TRUE <><","\n"))
  }

} # end of runstest plot function



#' Plot of process error deviation on log(biomass)
#'
#' shows the difference of the expected biomass and its stochastic realization
#'
#' @param jabba output list from fit_jabba
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param add if true don't call par() to allow construction of multiplots
#' @param width plot width
#' @param height plot height
#' @param verbose silent option
#' @export
jbplot_procdev <- function(jabba, output.dir=getwd(),as.png=FALSE,add=FALSE,width = 5, height = 3.5,verbose=TRUE){

  if(verbose) cat(paste0("\n","><> jbplot_procdev() - Process error diviations on log(biomass)  <><","\n"))

  years=jabba$yr
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  if(as.png==TRUE){png(file = paste0(output.dir,"/ProcDev_",jabba$assessment,"_",jabba$scenario,".png"), width = width, height = height,
      res = 200, units = "in")}
  if(add==FALSE) par(Par)
  ylim = c(min(-0.22,jabba$timeseries[,,"procB"]),max(0.22,jabba$timeseries[,,"procB"]))#range(proc.dev)*1.1
  cord.x <- c(years,rev(years))
  cord.y <- c(jabba$timeseries[,2,"procB"],rev(jabba$timeseries[,3,"procB"]))
  # Process Error
  plot(years,jabba$timeseries[,1,"procB"],ylab="Process Error Deviates",xlab="Year",ylim=ylim,type="n")
  polygon(cord.x,cord.y,col='grey',border=0,lty=2)
  lines(years,jabba$timeseries[,1,"procB"],lwd=2)
  lines(years,rep(0,length(years)),lty=5)
  if(as.png){dev.off()}
} # end of plot function


#' Plot of estimated trajectories
#'
#' plots choice of trajectories of Biomass, F, B/Bmsy, F/Fmsy or B/B0
#'
#' @param jabba output list from fit_jabba
#' @param type = c("B","F","BBmsy","FFmsy","BB0")
#' @param ylabs yaxis labels for quants
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param add if true don't call par() to allow construction of multiplots
#' @param mfrow set up plot frame  
#' @param width plot width
#' @param height plot height
#' @param verbose silent option
#' @export
jbplot_trj <-  function(jabba, type = c("B","F","BBmsy","FFmsy","BB0"),ylabs=NULL,output.dir=getwd(),as.png=FALSE,add=FALSE,mfrow=c(1,1),width=5,height=3.5,verbose=TRUE){

  for(i in 1:length(type)){

    Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
    if(as.png==TRUE){png(file = paste0(output.dir,"/",type[i],"_",jabba$assessment,"_",jabba$scenario,".png"), width = width, height = height,
        res = 200, units = "in")}
    if(add==FALSE){par(Par)}
    if(verbose) cat(paste0("\n","><> jbplot_trj() - ", type[i]," trajectory  <><","\n"))

    j = which(c("B","F","BBmsy","FFmsy","BB0")%in%type[i])
    if(is.null(ylabs)) ylabs = c(paste("Biomass",jabba$settings$catch.metric),ifelse(jabba$settings$harvest.label=="Fmsy","Fishing mortality F","Harvest rate H"),expression(B/B[MSY]),ifelse(jabba$settings$harvest.label=="Fmsy",expression(F/F[MSY]),expression(H/h[MSY])),expression(B/B[0]))
    
    trj = jabba$timeseries[,,paste(type[i])]
    years = jabba$yr
    ylim = c(0, max(trj[,3]))
    cord.x <- c(years,rev(years))
    cord.y <- c(trj[,2],rev(trj[,3]))
    plot(years,trj[,1],ylab=ylabs[j],xlab="Year",ylim=ylim,type="n")
    polygon(cord.x,cord.y,col='grey',border=0,lty=2)
    lines(years,trj[,1],lwd=2,col=1)
    if(type[i]=="B") lines(years,rep(jabba$refpts$bmsy[1],length(years)),lty=5)
    if(type[i]=="F") lines(years,rep(jabba$refpts$fmsy[1],length(years)),lty=5)
    if(type[i]%in%c("BBmsy","FFmsy") ) lines(years,rep(1,length(years)),lty=5)
    if(type[i]=="B") text((max(years)-min(years))/30+years[1],jabba$refpts$bmsy[1]*1.11,expression(paste(B[MSY])))
    if(type[i]=="F") text((max(years)-min(years))/30+years[1],jabba$refpts$fmsy[1]*1.11,ifelse(jabba$settings$harvest.label=="Fmsy",expression(F[MSY]),expression(H[MSY])))
    if(type[i]=="BB0"){
      lines(years,rep(jabba$refpts$bmsy[1]/jabba$refpts$k[1],length(years)),lty=5)
      text((max(years)-min(years))/30+years[1],jabba$refpts$bmsy[1]/jabba$refpts$k[1]*1.11,expression(paste(B[MSY])))
    }
    if(as.png==TRUE){
      if(add==FALSE | i==length(type)){dev.off()}}

    }} # end of plot function


#' Projection plot
#'
#' plots choice of projections of B/Bmsy, F/Fmsy or B/B0 over fixed quotas set up in build_jabba()
#'
#' @param jabba output list from fit_jabba
#' @param type Options are c("BBmsy","FFmsy","BB0")
#' @param CIs  Option to show CIs
#' @param flim max ylim value for FFmsy plot
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param add if true don't call par() to allow construction of multiplots
#' @param mfrow set up plot frame  
#' @param width plot width
#' @param height plot hight
#' @param cols option to add colour palette 
#' @export
jbplot_prj <-  function(jabba, type = c("BB0","BBmsy","FFmsy"),CIs=TRUE,flim=6,output.dir=getwd(),as.png=FALSE,add=FALSE,mfrow=c(1,1),width=5,height=3.5,cols=NULL){
  
  for(i in 1:length(type)){
    
    Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
    if(as.png==TRUE){png(file = paste0(output.dir,"/prj",type[i],"_",jabba$assessment,"_",jabba$scenario,".png"), width = width, height = height,
                         res = 200, units = "in")}
    if(add==FALSE){par(Par)}
    cat(paste0("\n","><> jbplot_prj() - ", type[i]," trajectory  <><","\n"))
    nTAC = jabba$settings$nTAC
    TACs = jabba$settings$TAC[nrow(jabba$settings$TAC),]
    imp.yr=jabba$settings$TAC.implementation
    if(is.null(cols)) cols = jabba$settings$cols
    shaded = rev(seq(0.4,0.9,0.5/nTAC))
    j = which(c("BB0","BBmsy","FFmsy")%in%type[i])
    ylabs = c(expression(B/B[0]),expression(B/B[MSY]),ifelse(jabba$settings$harvest.label=="Fmsy",expression(F/F[MSY]),expression(H/h[MSY])))
    prj = jabba$projections[,,,paste(type[i])]
    yrs = as.numeric(dimnames(jabba$projections)[[1]])
    ylim = c(0, max( c(max(jabba$projections[,ifelse(CIs,3,1),,paste(type[i])]),ifelse(type[i]=="BB0",1,2))))
    ylim[2] = min(flim,ylim[2])
    plot(yrs,yrs,ylim=ylim,xlim=c(min(yrs+0.3),max(yrs)),type="n",ylab=ylabs[j],xlab="Projection Years")
    
    if(CIs==TRUE){
    for(j in jabba$settings$nTAC:1)polygon(c(yrs,rev(yrs)),c(prj[,2,j],rev(prj[,3,j])),col=grey(shaded[j],1),border=NA)  
    }
    for(j in 1:jabba$settings$nTAC) lines(yrs,prj[,1,j],col=cols[j],lwd=2)  
    lines(yrs[1:(imp.yr-max(jabba$yr)+ifelse(type[i]=="FFmsy",0,1))],   prj[,1,1][1:(imp.yr-max(jabba$yr)++ifelse(type[i]=="FFmsy",0,1))],col=1,lwd=2.2)
    
    if(type[i]%in%c("BBmsy","FFmsy") ) lines(yrs,rep(1,length(yrs)),lty=5)
    if(type[i]=="BB0"){
      lines(c(yrs[1]-1,yrs[-1]),rep(jabba$refpts$bmsy[1]/jabba$refpts$k[1],length(yrs)),lty=5)
      text((max(yrs)-min(yrs)-1)/30+yrs[1],jabba$refpts$bmsy[1]/jabba$refpts$k[1]*1.11,expression(paste(B[MSY])))
    }
    legend("topleft",paste(TACs,"(t)"),col=(cols[1:nTAC]),lwd=2,cex=0.8)   
  
    if(as.png==TRUE){
      if(add==FALSE | i==length(type)){dev.off()}}
    
  }} # end of plot function


#' JABBA SPdyn Plot
#'
#' plots the production vs biomass (Walters et al 2008) and color-coded kobe phases
#'
#' @param jabba output list from fit_jabba
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param add if true don't call par() to allow construction of multiplots
#' @param width plot width
#' @param height plot height
#' @param verbose silent option
#' @export
jbplot_spdyn <-  function(jabba ,output.dir=getwd(),as.png=FALSE,add=FALSE,width=5,height=4.5,verbose=TRUE){
  if(verbose) cat(paste0("\n","><> jbplot_spdyn() - SPt+1 = Bt+1-Bt+Ct vs B  <><","\n"))

  # extract pars
  m = jabba$pfun$m[1]
  Bit = jabba$pfun$SB_i
  Cmsy = jabba$pfunc$Cmsy
  B = jabba$timeseries[,"mu","B"]
  Hmsy.sp = jabba$pfunc$Hmsy[1]
  SB0.sp =jabba$pfunc$SB0[1]
  SP = jabba$pfunc$SP
  Bmsy.sp = jabba$estimates["SBmsy",1]
  MSY.sp = jabba$estimates["MSY",1]
  N = jabba$settings$N
  years = jabba$yr
  ylim=c(min(c(1.2*jabba$timeseries[-1,1,"SPt"],0)),max(c(max(jabba$timeseries[-1,1,"SPt"],na.rm=T)*1.05,max(MSY.sp*1.1))))
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  if(as.png==TRUE){png(file = paste0(output.dir,"/SPdyn_",jabba$assessment,"_",jabba$scenario,".png"), width = width, height = height,
      res = 200, units = "in")}
  if(add==FALSE) par(Par)

  green.x = c(max(Bit,B*2),max(Bit,B*2),Bmsy.sp,Bmsy.sp,max(Bit,B*2))
  green.y = c(Bmsy.sp,ylim[1],ylim[1],max(SP),max(Cmsy))
  red.x = c(0,0,Bmsy.sp,Bmsy.sp,0)
  red.y = c(SB0.sp,0,max(SP),SB0.sp,SB0.sp)
  plot(Bit,SP,type = "n",ylim,xlim=c(0,max(Bit,B*1.1)),ylab="Surplus Production",xlab="Biomass",xaxs="i",yaxs="i")
  rect(0,ylim[1],SB0.sp*2,SB0.sp*1.1,col="green",border=0)
  rect(0,ylim[1],SB0.sp*2,SB0.sp,col="yellow",border=0)
  rect(0,max(SP),SB0.sp*2,SB0.sp,col="orange",border=0)
  polygon(green.x,green.y,border = 0,col="green")
  polygon(red.x,red.y,border = 0,col="red")
  abline(h=0,lty=2)
  ry.sp = Bit[Bit<=Bmsy.sp]
  for(i in 1:length(ry.sp)){

    lines(rep(Bit[i],2),c(Cmsy[i],SP[i]),col=ifelse(i %% 2== 0,"yellow","red"),lty=3)
    #i = i+1
  }

  lines(Bit,SP,col=4,lwd=2)
  lines(B,jabba$timeseries[,1,"SPt"],lty=1,lwd=1)
  points(B,jabba$timeseries[,1,"SPt"],cex=0.8,pch=16)
  sel.yr = c(1,round(quantile(1:N,0.7),0),N)
  points(B[sel.yr],jabba$timeseries[sel.yr,1,"SPt"],col= 1,pch=c(22,21,24),bg="white",cex=1.7)
  #abline(h=max(SP),col=4,lty=5)
  sel.years =years[sel.yr]
  #lines(rep(Bmsy.sp,2),c(0,max(SP)),lty=2,col=4)

  legend('topright',
         c("Expected","Predicted",paste(sel.years)),
         lty=c(1,1,1,1,1),pch=c(-1,16,22,21,24),pt.bg=c(0,0,rep("white",3)),
         col=c(4,rep(1,4)),lwd=c(2,1,1,1,1),cex=0.8,pt.cex=c(-1,1,rep(1.3,3)),bty="n")

  if(as.png==TRUE) dev.off()
} #end of plotting function


#' JABBA SPphase Plot
#'
#' plots the production function + Catch vs biomass and color-coded kobe phases
#'
#' @param jabba output list from fit_jabba
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param add if true don't call par() to allow construction of multiplots
#' @param width plot width
#' @param height plot hight
#' @param verbose silent option
#' @export

jbplot_spphase <-  function(jabba ,output.dir=getwd(),as.png=FALSE,add=FALSE,width=5,height=4.5,verbose=TRUE){
  if(verbose) cat(paste0("\n","><> jbplot_spphase() - JABBA Surplus Production Phase Plot  <><","\n"))
  
  # extract pars
  m = jabba$pfun$m[1]
  Bit = jabba$pfun$SB_i
  Cmsy = jabba$pfunc$Cmsy
  B = jabba$timeseries[,"mu","B"]
  Hmsy.sp = jabba$pfunc$Hmsy[1]
  SB0.sp =jabba$pfunc$SB0[1]
  SP = jabba$pfunc$SP
  Bmsy.sp = jabba$estimates["SBmsy",1]
  MSY.sp = jabba$estimates["MSY",1]
  N = jabba$settings$N
  years = jabba$yr
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  if(as.png==TRUE){png(file = paste0(output.dir,"/SPphase_",jabba$assessment,"_",jabba$scenario,".png"), width = width, height = height,
                       res = 200, units = "in")}
  if(add==FALSE) par(Par)
  
  green.x = c(max(Bit,B),max(Bit,B),Bmsy.sp,Bmsy.sp,max(Bit))
  green.y = c(Bmsy.sp,0,0,max(SP),max(Cmsy))
  red.x = c(0,0,Bmsy.sp,Bmsy.sp,0)
  red.y = c(SB0.sp,0,max(SP),SB0.sp,SB0.sp)
  plot(Bit,SP,type = "n",ylim=c(0,max(c(max(jabba$catch,na.rm=T)*1.05,max(MSY.sp*1.1)))),xlim=c(0,max(Bit,B)),ylab="Surplus Production",xlab="Biomass",xaxs="i",yaxs="i")
  rect(0,0,SB0.sp*1.1,SB0.sp*1.1,col="green",border=0)
  rect(0,0,SB0.sp,SB0.sp,col="yellow",border=0)
  rect(0,max(SP),SB0.sp,SB0.sp,col="orange",border=0)
  polygon(green.x,green.y,border = 0,col="green")
  polygon(red.x,red.y,border = 0,col="red")
  
  ry.sp = Bit[Bit<=Bmsy.sp]
  for(i in 1:length(ry.sp)){
    
    lines(rep(Bit[i],2),c(Cmsy[i],SP[i]),col=ifelse(i %% 2== 0,"yellow","red"),lty=3)
    #i = i+1
  }
  
  polygon(c(-10000,10^7,10^7,-10000),c(rep(MSY.sp[1],2),rep(MSY.sp[3],2)),border = FALSE,col=rgb(0,0,1,0.4))
  lines(Bit,SP,col=4,lwd=2)
  lines(B,jabba$catch,lty=1,lwd=1)
  points(B,jabba$catch,cex=0.8,pch=16)
  sel.yr = c(1,round(quantile(1:N,0.7),0),N)
  points(B[sel.yr],jabba$catch[sel.yr],col= 1,pch=c(22,21,24),bg="white",cex=1.7)
  abline(h=max(SP),col=4,lty=5)
  sel.years =years[sel.yr]
  lines(rep(Bmsy.sp,2),c(-1000,max(SP)),lty=2,col=4)
  
  legend('topright',
         c(expression(B[MSY]),"MSY","SP","Catch",paste(sel.years)),
         lty=c(2,5,1,1,1,1,1),pch=c(-1,-1,-1,16,22,21,24),pt.bg=c(0,0,0,0,rep("white",3)),
         col=c(4,4,4,rep(1,4)),lwd=c(1,1,2,1,1,1),cex=0.8,pt.cex=c(-1,-1,-1,0.5,rep(1.3,3)),bty="n")
  
  if(as.png==TRUE) dev.off()
} #end of plotting function


#' KOBE phase plot
#'
#' plots the stock status posterior over B/Bmsy and F/Fmsy
#'
#' @param jabba output list from fit_jabba
#' @param output.dir directory to save plots
#' @param ylab yaxis label
#' @param xlab xaxis label
#' @param as.png save as png file of TRUE
#' @param add if true don't call par() to allow construction of multiplots
#' @param width plot width
#' @param height plot height
#' @param verbose silent option
#' @export
jbplot_kobe <-  function(jabba ,ylab=NULL,xlab=NULL, output.dir=getwd(),as.png=FALSE,add=FALSE,width=5,height=4.5,verbose=TRUE){

  if(verbose) cat(paste0("\n","><> jbplot_kobe() - Stock Status Plot  <><","\n"))

  mu.f = jabba$timeseries[,,"FFmsy"]
  mu.b = jabba$timeseries[,,"BBmsy"]
  f = jabba$kobe$harvest
  b = jabba$kobe$stock
  years=jabba$yr
  N = length(years)
  # fit kernel function
  kernelF <- gplots::ci2d(b,f,nbins=151,factor=1.5,ci.levels=c(0.50,0.80,0.75,0.90,0.95),show="none",col=1,xlab= ifelse(jabba$settings$harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))),ylab=expression(paste(B/B[MSY])))

  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  if(as.png==TRUE){ png(file = paste0(output.dir,"/Kobe_",jabba$assessment,"_",jabba$scenario,".png"), width = width, height = height,
      res = 200, units = "in")}
  if(add==FALSE) par(Par)

  #Create plot
  if(is.null(ylab)) ylab=ifelse(jabba$settings$harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY])))
  if(is.null(xlab)) xlab=expression(paste(B/B[MSY]))
  
  plot(1000,1000,type="b", xlim=c(0,max(1/(jabba$refpts$bmsy/jabba$refpts$k)[1],mu.b[,1]) +0.05), ylim=c(0,max(mu.f[,1],quantile(f,0.85),2.)),lty=3,ylab=ylab,xlab=xlab,xaxs="i",yaxs="i")
  c1 <- c(-1,100)
  c2 <- c(1,1)

  # extract interval information from ci2d object
  # and fill areas using the polygon function
  zb2 = c(0,1)
  zf2  = c(1,100)
  zb1 = c(1,100)
  zf1  = c(0,1)
  polygon(c(zb1,rev(zb1)),c(0,0,1,1),col="green",border=0)
  polygon(c(zb2,rev(zb2)),c(0,0,1,1),col="yellow",border=0)
  polygon(c(1,100,100,1),c(1,1,100,100),col="orange",border=0)
  polygon(c(0,1,1,0),c(1,1,100,100),col="red",border=0)

  polygon(kernelF$contours$"0.95",lty=2,border=NA,col="cornsilk4")
  polygon(kernelF$contours$"0.8",border=NA,lty=2,col="grey")
  polygon(kernelF$contours$"0.5",border=NA,lty=2,col="cornsilk2")

  points(mu.b[,1],mu.f[,1],pch=16,cex=1)
  lines(c1,c2,lty=3,lwd=0.7)
  lines(c2,c1,lty=3,lwd=0.7)
  lines(mu.b[,1],mu.f[,1], lty=1,lwd=1.)
  sel.yr = c(1,round(quantile(1:N,0.7),0),N)
  points(mu.b[sel.yr,1],mu.f[sel.yr,1],col=
           1,pch=c(22,21,24),bg="white",cex=1.9)

  # Get Propability
  Pr.green = sum(ifelse(b>1 & f<1,1,0))/length(b)*100
  Pr.red = sum(ifelse(b<1 & f>1,1,0))/length(b)*100
  Pr.yellow = sum(ifelse(b<1 & f<1,1,0))/length(b)*100
  Pr.orange = sum(ifelse(b>1 & f>1,1,0))/length(b)*100



  sel.years = c(years[sel.yr])
  ## Add legend

  legend('topright',
         c(paste(sel.years),"50% C.I.","80% C.I.","95% C.I.",paste0(round(c(Pr.red,Pr.yellow,Pr.orange,Pr.green),1),"%")),
         lty=c(1,1,1,rep(-1,8)),pch=c(22,21,24,rep(22,8)),pt.bg=c(rep("white",3),"cornsilk2","grey","cornsilk4","red","yellow","orange","green"),
         col=1,lwd=1.1,cex=0.9,pt.cex=c(rep(1.3,3),rep(1.7,3),rep(2.2,4)),bty="n")

  if(as.png==TRUE){dev.off()}
} # End of Kobe plot


#---------------------------------------------------------
# Produce 'post-modern' biplot (see Quinn and Collie 2005)
#---------------------------------------------------------


#' Stock Status Biplot
#'
#' plots 'post-modern' biplot (see Quinn and Collie 2005)
#'
#' @param jabba output list from fit_jabba
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param add if true don't call par() to allow construction of multiplots
#' @param width plot width
#' @param height plot hight
#' @param verbose silent option
#' @export
jbplot_biplot <-  function(jabba ,output.dir=getwd(),as.png=FALSE,add=FALSE,width=5,height=4.5,verbose=TRUE){

  if(verbose) cat(paste0("\n","><> jbplot_biplot() - Stock Status Plot  <><","\n"))
  mu.f = jabba$timeseries[,,"FFmsy"]
  mu.b = jabba$timeseries[,,"BBmsy"]
  f = jabba$kobe$harvest
  b = jabba$kobe$stock
  years=jabba$yr
  N = length(years)

  # fit kernel function
  kernelF <- gplots::ci2d(f,b,nbins=201,factor=2,ci.levels=c(0.50,0.80,0.75,0.90,0.95),show="none",col=1,ylab= ifelse(jabba$settings$harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))),xlab=expression(paste(B/B[MSY])))


  Par = list(mfrow=c(1,1),mai=c(0.2,0.15,0,.15),omi = c(0.3,0.25,0.2,0) + 0.1, mgp =c(3,1,0), tck = -0.02,cex=0.8)
  if(as.png==TRUE){ png(file = paste0(output.dir,"/Biplot_",jabba$assessment,"_",jabba$cenario,".png"), width = width, height = height,
      res = 200, units = "in")}
  if(add==FALSE) par(Par)

  #Create plot
  plot(1000,1000,type="b", ylim=c(0,max(1/(jabba$refpts$bmsy/jabba$refpts$k)[1],mu.b[,1]) +0.05), xlim=c(0,max(mu.f[,1],quantile(f,0.85),2.)),lty=3,xlab=ifelse(jabba$settings$harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))),ylab=expression(paste(B/B[MSY])),xaxs="i",yaxs="i")

  # and fill areas using the polygon function
  fint = seq(0.001,100,0.01)
  # read ftarget,bthreshold
  ftarget<-0.8
  bthreshold<-0.2

  #Zone X
  xb=bthreshold+(1.0-bthreshold)/ftarget*fint
  xf =  ifelse(xb>1,0.8,fint)
  polygon(c(0,0,xf),c(max(xb),bthreshold,xb),col="green")
  zb = bthreshold+(1.0-bthreshold)*fint
  zf  = ifelse(zb>1,1,fint)
  polygon(c(zf,rep(max(fint),2),rep(0,2)),c(zb,max(zb),0,0,bthreshold),col="red")

  polygon(c(xf,rev(zf)),c(xb,rev(zb)),col="yellow")

  c1 <- c(-1,100)
  c2 <- c(1,1)

  # extract interval information from ci2d object
  # and fill areas using the polygon function
  polygon(kernelF$contours$"0.95",lty=2,border=NA,col="cornsilk4")
  polygon(kernelF$contours$"0.8",border=NA,lty=2,col="grey")
  polygon(kernelF$contours$"0.5",border=NA,lty=2,col="cornsilk2")
  points(mu.f[2,],mu.b[2,],pch=16,cex=1)

  lines(c1,c2,lty=3,lwd=0.7)
  lines(c2,c1,lty=3,lwd=0.7)
  lines(mu.f[,1],mu.b[,1], lty=1,lwd=1.)
  sel.yr = c(1,round(quantile(1:N,0.7),0),N)
  points(mu.f[sel.yr,1],mu.b[sel.yr,1],col=
           1,pch=c(22,21,24),bg="white",cex=1.9)


  sel.years = years[sel.yr]
  ## Add legend
  legend('topright',
         c(paste(sel.years),"50% C.I.","80% C.I.","95% C.I."),
         lty=c(1,1,1,-1,-1,-1),pch=c(22,21,24,22,22,22),pt.bg=c(rep("white",3),"cornsilk2","grey","cornsilk4"),
         col=1,lwd=1.1,cex=0.9,pt.cex=c(rep(1.3,4),1.7,1.7,1.7),bty="n")

  Zone  = NULL
  Status = NULL
  X  = 0.15
  Y = 0
  Z = -0.15

  for(i  in 1:length(f))
  {
    if(b[i]>1.0){
      if(f[i]<ftarget){
        Zone[i]<-X
      } else if (f[i]>1.0){
        Zone[i]<-Z
      } else {
        Zone[i]<-Y
      }
    } else {
      if(b[i]>bthreshold+(1.0-bthreshold)/ftarget*f[i]){
        Zone[i]<-X
      } else if(b[i]<bthreshold+(1.0-bthreshold)*f[i]){
        Zone[i]<-Z
      } else {
        Zone[i]<-Y
      }


    }}

  perGreen = round(length(Zone[Zone==0.15])/length(Zone)*100,1)
  perYellow = round(length(Zone[Zone==0])/length(Zone)*100,1)
  perRed = round(length(Zone[Zone==-0.15])/length(Zone)*100,1)

  mtext(expression(paste(B/B[MSY])), side=2, outer=TRUE, at=0.5,line=1,cex=0.9)
  mtext(ifelse(jabba$settings$harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))), side=1, outer=TRUE, at=0.5,line=1,cex=0.9)

  text(0.65,2.4,paste0(perGreen,"%"))
  text(0.9,2.4,paste0(perYellow,"%"))
  text(1.2,2.4,paste0(perRed,"%"))

  if(as.png==TRUE) dev.off()

} # End of biplot function

#' Plot of biomass depletion prior vs posterios
#'
#' prior is specified as b/k or b/bmsy in build_jabba()
#'
#' @param jabba output list from fit_jabba
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param add if true don't call par() to allow construction of multiplots
#' @param width plot width
#' @param height plot height
#' @param verbose silent option
#' @export
jbplot_bprior <- function(jabba, output.dir=getwd(),as.png=FALSE,add=FALSE,width = 5, height = 3.5,verbose=TRUE){
  if(jabba$settings$b.pr[3]==0){
   if(verbose) cat(paste0("\n","><> No additional biomass depletion prior specified  <><","\n"))
  } else {
  if(verbose) cat(paste0("\n","><> jbplot_bprior - biomass depletion prior vs posterior   <><","\n"))
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  if(as.png==TRUE){png(file = paste0(output.dir,"/Bprior_",jabba$assessment,"_",jabba$scenario,".png"), width = width, height = height,
                       res = 200, units = "in")}
  if(add==FALSE) par(Par)
  
  xlabs = c(expression(B/K),expression(B/B[MSY]),expression(F/F[MSY])) 
  bpr =  rlnorm(10000,log(jabba$settings$b.pr[1]),jabba$settings$b.pr[2])
  pdf = stats::density(jabba$bppd,adjust=2)
  prior = dlnorm(sort(bpr),log(jabba$settings$b.pr[1]),jabba$settings$b.pr[2])
  plot(pdf,type="l",ylim=c(0,max(prior,pdf$y)*1.1),xlim=range(c(pdf$x,quantile(bpr,c(0.0001,0.95)))),yaxt="n",xlab=xlabs[jabba$settings$b.pr[4]+1],ylab="Density",xaxs="i",yaxs="i",main="")

 polygon(c(sort(bpr),rev(sort(bpr))),c(prior,rep(0,length(sort(bpr)))),col=gray(0.4,1))
 polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
 legend('right',c("Prior","Posterior"),pch=22,pt.cex=1.5,pt.bg = c(grey(0.4,1),grey(0.8,0.6)),bty="n")
 PPVR = round((sd(jabba$bppd)/mean(jabba$bppd))^2/(sd(bpr)/mean(bpr))^2,3)
 PPVM = round(mean(jabba$bppd)/mean(bpr),3)
 legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=1,bty="n")
 legend("top",c(paste(round(jabba$settings$b.pr[3],0))),cex=1.2,bty="n",y.intersp = -0.2)
 if(as.png==TRUE){dev.off()}
  }
} # end of biomass prior plotting function 



#' wrapper jbplot function
#'
#' plots all routine JABBA plot to output.dir if as.png=TRUE (default)
#'
#' @param jabba output list from fit_jabba
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param verbose silent option
#' @export
jabba_plots = function(jabba,output.dir = getwd(),as.png=TRUE,statusplot ="kobe",verbose=TRUE){

  jbplot_catch(jabba,as.png=as.png,output.dir=output.dir,verbose=verbose) # catch.metric
  jbplot_catcherror(jabba,as.png=as.png,output.dir=output.dir,verbose=verbose) 
  jbplot_cpuefits(jabba,as.png=as.png,output.dir=output.dir,verbose=verbose) # check years
  jbplot_logfits(jabba,as.png=as.png,output.dir=output.dir,verbose=verbose) # check n.indices
  jbplot_mcmc(jabba,as.png=as.png,output.dir=output.dir,verbose=verbose)
  jbplot_ppdist(jabba,as.png=as.png,output.dir=output.dir,verbose=verbose) # check m
  jbplot_procdev(jabba,as.png=as.png,output.dir=output.dir,verbose=verbose)
  jbplot_trj(jabba,as.png=as.png,output.dir=output.dir,verbose=verbose)
  jbplot_spphase(jabba,as.png=as.png,output.dir=output.dir,verbose=verbose) # check TC
  jbplot_residuals(jabba,as.png=as.png,output.dir=output.dir,verbose=verbose) # check years
  jbplot_runstest(jabba,as.png=as.png,output.dir=output.dir,verbose=verbose)
  jbplot_summary(jabba,as.png,output.dir=output.dir)
  
  if(statusplot =="kobe"){
    jbplot_kobe(jabba,as.png=as.png,output.dir=output.dir,verbose=verbose)} else {
      jbplot_biplot(jabba,as.png=as.png,output.dir=output.dir,verbose=verbose)}
}



#' jbplot_retro() to plot retrospective pattern
#'
#' Plots retrospective pattern of B, F, BBmsy, FFmsy, BB0 and SP #'
#' @param hc output list from hindast_jabba()
#' @param type  single plot option c("B","F","BBmsy","FFmsy","BB0","SP")
#' @param forecast  includes retrospective forecasting if TRUE
#' @param ylabs yaxis labels for quants
#' @param add  add to multi plot if TRUE
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param single.plots if TRUE plot invidual fits else make multiplot
#' @param width plot width
#' @param height plot hight
#' @param xlim  allows to "zoom-in" requires speficiation Xlim=c(first.yr,last.yr)
#' @param cols option to add colour palette 
#' @param legend.loc location of legend
#' @param verbose if FALSE be silent
#' @return Mohn's rho statistic for several quantaties
#' @export
#' @examples 
#' data(iccat)
#' bet= iccat$bet
#' jb = build_jabba(catch=bet$catch,cpue=bet$cpue,se=bet$se,assessment="BET",scenario = "Ref",model.type = "Pella",igamma = c(0.001,0.001),verbose=FALSE)
#' fit = fit_jabba(jb,quickmcmc=TRUE,verbose=FALSE)
#' hc = hindcast_jabba(jbinput=jb,fit=fit,peels=1:3)
#' jbplot_retro(hc)
#' jbplot_retro(hc,forecast=TRUE) # with retro forecasting

jbplot_retro <- function(hc,type=c("B","F","BBmsy","FFmsy","procB","SP"),forecast=FALSE,ylabs=NULL,
                         add=F,output.dir=getwd(),as.png=FALSE,single.plots=add,width=NULL,height=NULL,xlim=NULL,cols=NULL,legend.loc="topright",verbose=TRUE){
  
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
  
  
  if(verbose)cat(paste0("\n","><> jbplot_retro() - retrospective analysis <><","\n"))
  if(add) single.plots=TRUE
  if(single.plots==F) type=c("B","F","BBmsy","FFmsy","procB","SP")
  
  if(is.null(ylabs)) ylabs = c(paste("Biomass",hc$settings$catch.metric),"Fishing mortality F",expression(B/B[MSY]),expression(F/F[MSY]),expression(B/B[0]),"Process Deviations",paste("Surplus Production",hc$settings$catch.metric))
  retros = unique(peels)
  runs= hc$timeseries$mu$level
  years= hc$yr
  nyrs = length(years)
  if(is.null(cols)) cols = c("black",ss3col(length(peels)-1))
  if(is.null(xlim)){xlim = range(years)}
  FRP.rho = c("B","F", "Bmsy", "Fmsy", "procB","MSY")  
  rho = data.frame(mat.or.vec(length(retros)-1,length(FRP.rho)))
  colnames(rho) = FRP.rho
  fcrho = rho
  
  if(single.plots==TRUE){
    if(is.null(width)) width = 5
    if(is.null(height)) height = 3.5
    for(k in 1:length(type)){
      Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
      if(as.png==TRUE){png(file = paste0(output.dir,"/Retro",hc$scenario ,"_",type[k],".png"), width = width, height = height,
                           res = 200, units = "in")}
      
      
      if(as.png==TRUE | add==FALSE) par(Par)
      
      j = which(c("B","F","BBmsy","FFmsy","BB0","procB","SP")%in%type[k])
      
      
      if(type[k]%in%c("B","F","BBmsy","FFmsy","procB")){
        y = hc$timeseries$mu[,j+2]
        ref = hc$timeseries$mu[runs%in%retros[1],j+2]
        ylc = hc$timeseries$lci[runs%in%retros[1],j+2]
        yuc = hc$timeseries$uci[runs%in%retros[1],j+2]
        if(type[k]=="procB") ylim=c(-max(y[years>=xlim[1] & years<=xlim[2]],yuc[years>=xlim[1] & years<=xlim[2]])
                                    ,max(y[years>=xlim[1] & years<=xlim[2]],yuc[years>=xlim[1] & years<=xlim[2]]))
        if(!type[k]=="procB") ylim=c(0,max(y[years>=xlim[1] & years<=xlim[2]],yuc[years>=xlim[1] & years<=xlim[2]]))
        
        plot(years,years,type="n",ylim=ylim,ylab=ifelse(length(ylabs)>1,ylabs[j],ylabs),xlab="Year",xlim=xlim)
        polygon(c(years,rev(years)),c(ylc,rev(yuc)),col="grey",border="grey")
        for(i in 1:length(retros)){
          lines(years[1:(nyrs-retros[i])],y[runs%in%retros[i]][1:(nyrs-retros[i])],col= cols[i],lwd=ifelse(i==1,2,1.5),lty=1)
          if(forecast){
            lines(years[(nyrs-retros[i]):(nyrs+1-retros[i])],y[runs%in%retros[i]][(nyrs-retros[i]):(nyrs+1-retros[i])],col= cols[i],lwd=1,lty=2)
            points(years[(nyrs+1-retros[i])],y[runs%in%retros[i]][(nyrs+1-retros[i])],pch=16,col= cols[i],cex=0.8)
          }
            
          if(i>1){
            rho[i-1,k] =  (y[runs%in%retros[i]][(nyrs-retros[i])]-ref[(nyrs-retros[i])])/ref[(nyrs-retros[i])]
            fcrho[i-1,k] = (y[runs%in%retros[i]][(nyrs+1-retros[i])]-ref[(nyrs+1-retros[i])])/ref[(nyrs+1-retros[i])]
            if(type[k]=="procB"){
              rho[i-1,k] =  (exp(y[runs%in%retros[i]][(nyrs-retros[i])])-exp(ref[(nyrs-retros[i])]))/exp(ref[(nyrs-retros[i])])
              fcrho[i-1,k] =  (exp(y[runs%in%retros[i]][(nyrs+1-retros[i])])-exp(ref[(nyrs+1-retros[i])]))/exp(ref[(nyrs+1-retros[i])])
            }
          }
        }
        if(type[k]%in%c("BBmsy","FFmsy")) abline(h=1,lty=2)
        if(type[k]%in%c("procB")) abline(h=0,lty=2)
      } 
      else {
        # Plot SP
        plot(years,years,type="n",ylim=c(0,max(hc$pfunc$SP*1.12)),xlim=c(0,max(hc$pfunc$SB_i)),ylab=ifelse(length(ylabs)>1,ylabs[j],ylabs),xlab=ylabs[1])
        for(i in 1:length(retros)){
          lines(hc$pfunc$SB_i[hc$pfunc$level%in%retros[i]],hc$pfunc$SP[hc$pfunc$level%in%retros[i]],col=cols[i],lwd=ifelse(i==1,2,1.5),lty=1)
          points(mean(hc$pfunc$SB_i[hc$pfunc$level%in%retros[i]][hc$pfunc$SP[hc$pfunc$level%in%retros[i]]==max(hc$pfunc$SP[hc$pfunc$level%in%retros[i]])]),max(hc$pfunc$SP[hc$pfunc$level%in%retros[i]]),col=cols[i],pch=16,cex=1.2)
          if(i>1){
            rho[i-1,6] =  (hc$refpts$msy[hc$refpts$level==retros[i]]-hc$refpts$msy[hc$refpts$level==retros[1]])/hc$refpts$msy[hc$refpts$level==retros[1]]
            fcrho[i-1,6] = NA
          }
        }}
      if(single.plots==TRUE | k==1 )  legend(legend.loc,paste(years[nyrs-retros]),col=cols,bty="n",cex=0.7,pt.cex=0.7,lwd=c(2,rep(1.5,length(retros))))
      if(!forecast) legend("top",legend=bquote(rho == .(round(mean(rho[,k]),2))) ,bty="n",x.intersp=-0.2,y.intersp=-0.3,cex=0.8)
      if(forecast)legend("top",legend=bquote(rho == .(round(mean(rho[,k]),2))~"("~.(round(mean(fcrho[,k]),2))~")") ,bty="n",x.intersp=-0.2,y.intersp=-0.3,cex=0.8)
      
      
      if(as.png==TRUE) dev.off()
    } # End type loop
  } else { # Multi plot
    if(is.null(width)) width = 7
    if(is.null(height)) height = 8 
    Par = list(mfrow=c(3,2),mai=c(0.45,0.49,0.1,.15),omi = c(0.15,0.15,0.1,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
    if(as.png==TRUE){png(file = paste0(output.dir,"/Retro_",hc$scenario,".png"), width = width, height = height,
                         res = 200, units = "in")}
    par(Par)
    for(k in 1:length(type)){
      
      j = which(c("B","F","BBmsy","FFmsy","BB0","procB","SP")%in%type[k])
      
      
      if(type[k]%in%c("B","F","BBmsy","FFmsy","procB")){
        y = hc$timeseries$mu[,j+2]
        ref = hc$timeseries$mu[runs%in%retros[1],j+2]
        ylc = hc$timeseries$lci[runs%in%retros[1],j+2]
        yuc = hc$timeseries$uci[runs%in%retros[1],j+2]
        if(type[k]=="procB") ylim=c(-max(y[years>=xlim[1] & years<=xlim[2]],yuc[years>=xlim[1] & years<=xlim[2]])
                                    ,max(y[years>=xlim[1] & years<=xlim[2]],yuc[years>=xlim[1] & years<=xlim[2]]))
        if(!type[k]=="procB") ylim=c(0,max(y[years>=xlim[1] & years<=xlim[2]],yuc[years>=xlim[1] & years<=xlim[2]]))
        
        plot(years,years,type="n",ylim=ylim,ylab=ylabs[j],xlab="Year",xlim=xlim)
        polygon(c(years,rev(years)),c(ylc,rev(yuc)),col="grey",border="grey")
        for(i in 1:length(retros)){
          lines(years[1:(nyrs-retros[i])],y[runs%in%retros[i]][1:(nyrs-retros[i])],col= cols[i],lwd=ifelse(i==1,2,1.5),lty=1)
          if(forecast){
            lines(years[(nyrs-retros[i]):(nyrs+1-retros[i])],y[runs%in%retros[i]][(nyrs-retros[i]):(nyrs+1-retros[i])],col= cols[i],lwd=1,lty=2)
            points(years[(nyrs+1-retros[i])],y[runs%in%retros[i]][(nyrs+1-retros[i])],pch=16,col= cols[i],cex=0.8)
          }
          if(i>1){
            rho[i-1,k] =  (y[runs%in%retros[i]][(nyrs-retros[i])]-ref[(nyrs-retros[i])])/ref[(nyrs-retros[i])]
            fcrho[i-1,k] = (y[runs%in%retros[i]][(nyrs+1-retros[i])]-ref[(nyrs+1-retros[i])])/ref[(nyrs+1-retros[i])]
            if(type[k]=="procB"){
              rho[i-1,k] =  (exp(y[runs%in%retros[i]][(nyrs-retros[i])])-exp(ref[(nyrs-retros[i])]))/exp(ref[(nyrs-retros[i])])
              fcrho[i-1,k] =  (exp(y[runs%in%retros[i]][(nyrs-retros[i])])-exp(ref[(nyrs+1-retros[i])]))/exp(ref[(nyrs+1-retros[i])])
            }
          }
        }
        if(type[k]%in%c("BBmsy","FFmsy")) abline(h=1,lty=2)
        if(type[k]%in%c("procB")) abline(h=0,lty=2)
        if(single.plots==TRUE | k==1 )  legend(legend.loc,paste(years[nyrs-retros]),col=cols,bty="n",cex=0.7,pt.cex=0.7,lwd=c(2,rep(1.5,length(retros))))
        if(!forecast) legend("top",legend=bquote(rho == .(round(mean(rho[,k]),2))) ,bty="n",x.intersp=-0.2,y.intersp=-0.3,cex=0.8)
        if(forecast)legend("top",legend=bquote(rho == .(round(mean(rho[,k]),2))~"("~.(round(mean(fcrho[,k]),2))~")") ,bty="n",x.intersp=-0.2,y.intersp=-0.3,cex=0.8)
        
      }  else {
        # Plot SP
        plot(years,years,type="n",ylim=c(0,max(hc$pfunc$SP*1.15)),xlim=c(0,max(hc$pfunc$SB_i)),ylab=ylabs[j],xlab=ylabs[1])
        for(i in 1:length(retros)){
          lines(hc$pfunc$SB_i[hc$pfunc$level%in%retros[i]],hc$pfunc$SP[hc$pfunc$level%in%retros[i]],col=cols[i],lwd=ifelse(i==1,2,1.5),lty=1)
          points(mean(hc$pfunc$SB_i[hc$pfunc$level%in%retros[i]][hc$pfunc$SP[hc$pfunc$level%in%retros[i]]==max(hc$pfunc$SP[hc$pfunc$level%in%retros[i]])]),max(hc$pfunc$SP[hc$pfunc$level%in%retros[i]]),col=cols[i],pch=16,cex=1.2)
          if(i>1){
            rho[i-1,6] =  (hc$refpts$msy[hc$refpts$level==retros[i]]-hc$refpts$msy[hc$refpts$level==retros[1]])/hc$refpts$msy[hc$refpts$level==retros[1]]
            fcrho[i-1,6] = NA
          }      
        }
        legend("top",legend=bquote(rho == .(round(mean(rho[,k]),2))) ,bty="n",x.intersp=-0.2,y.intersp=-0.3,cex=0.8)
        
      }
      
      
    }
    if(as.png==TRUE) dev.off()
  }
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
    
  if(verbose) return(out)
} # end of Retrospective Plot


#' jbplot_summary() 
#'
#' Compares B, F, BBmsy, FFmsy, BB0 and SP for various model scanarios that have to be saved as rdata 
#' @param jabbas list() of JABBA model 1:n
#' @param type for single plots optional select type=c("B","F","BBmsy","FFmsy","BB0","SP")
#' @param plotCIs Plot Credibilty Interval 
#' @param ylabs yaxis labels for quants
#' @param prefix Plot name specifier
#' @param save.summary option to save a summary of all loaded model runs
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param single.plots if TRUE plot invidual fits else make multiplot
#' @param add adds single plots to mfrow set up
#' @param width plot width
#' @param height plot hight
#' @param xlim allows to "zoom-in" requires speficiation xlim=c(first.yr,last.yr)
#' @param cols option to add colour palette 
#' @param legend.loc location of legend
#' @param legend.cex size of legend
#' @param legend.add show legend
#' @param plot.cex cex setting in par()
#' @param verbose if FALSE be silent
#' @export
jbplot_summary <- function(jabbas,type=c("B","F","BBmsy","FFmsy","BB0","SP"),plotCIs=TRUE,ylabs=NULL,prefix="Summary",save.summary=FALSE,output.dir=getwd(),as.png=FALSE,single.plots=add,add=FALSE,width=NULL,height=NULL,xlim=NULL,cols=NULL,legend.loc="topright",legend.cex=0.8,legend.add=TRUE,plot.cex=0.8,verbose=TRUE){
  
  if(!is.null(jabbas$settings)) jabbas = list(jabbas)
  if(as.png) add=FALSE
  if(is.null(cols)) cols= ss3col(length(jabbas))
  if(verbose) cat(paste0("\n","><> jbplot_summary()"))
  jbs = list(assessment=jabbas[[1]]$assessment,yr= NULL,catch=NULL,timeseries = NULL,refpts=NULL,pfunc=NULL,settings=NULL)
  if(single.plots==F) type=c("B","F","BBmsy","FFmsy","BB0","SP")
  scenarios = NULL
  for(i in 1:length(jabbas)) scenarios = c(scenarios,jabbas[[i]]$scenario)
  if(is.null(names(jabbas))){
    scenarios = scenarios
    if(length(unique(scenarios))<length(scenarios)){
      scenarios = paste0(scenarios,1:length(scenarios))
    } 
  } else {
    scenarios = names(jabbas)
  } 
  
  
  for(i in 1:length(jabbas)){
    jabba = jabbas[[i]]
    if(i==1){
      jbs$yr = jabba$yr
      jbs$catch = jabba$catch
      jbs$settings$cols = jabba$settings$cols
      jbs$settings$harvest = jabba$settings$harvest.label
      jbs$settings$catch.metric = jabba$settings$catch.metric 
      
    }
    jabba$pfunc$level = scenarios[i]
   
    jbs$timeseries$mu = rbind(jbs$timeseries$mu,data.frame(year=jabba$yr,factor=jabba$pfunc[1,1],level=jabba$pfunc[1,2],jabba$timeseries[,"mu",])) 
    jbs$timeseries$lci = rbind(jbs$timeseries$lci,data.frame(year=jabba$yr,factor=jabba$pfunc[1,1],level=jabba$pfunc[1,2],jabba$timeseries[,"lci",])) 
    jbs$timeseries$uci = rbind(jbs$timeseries$uci,data.frame(year=jabba$yr,factor=jabba$pfunc[1,1],level=jabba$pfunc[1,2],jabba$timeseries[,"uci",])) 
    #jbs$diags = rbind(jbs$diags,jabba$diags) #prevents Catch-Only comparison 
    jbs$refpts= rbind(jbs$refpts,jabba$refpts[1,])
    jbs$pfunc= rbind(jbs$pfunc ,jabba$pfunc)
    
   }
  
  if(save.summary){
    save(jbs,file=paste0(output.dir,"/",prefix,"_",assessment,"_summary.rdata"))
  }
  
  
  #type=c("B","F","BBmsy","FFmsy","BB0","SP")
  if(is.null(ylabs)) ylabs = c(paste("Biomass",jabba$settings$catch.metric),ifelse(jabba$settings$harvest=="Fmsy","Fishing mortality F","Harvest rate H"),expression(B/B[MSY]),ifelse(jabba$settings$harvest=="Fmsy",expression(F/F[MSY]),expression(H/H[MSY])),expression(B/B[0]),paste("Surplus Production",jabba$settings$catch.metric))
  runs= jbs$timeseries$mu$level
  years= jbs$yr
  nyrs = length(years)
  if(length(scenarios)>1){cols= c(cols)}else(cols=1)
  
  if(is.null(xlim)){xlim = range(years)}
  
  if(single.plots==TRUE){
    if(is.null(width)) width = 5
    if(is.null(height)) height = 3.5
    for(k in 1:length(type)){
      Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=plot.cex)
      if(as.png==TRUE){png(file = paste0(output.dir,"/",prefix,"_",jbs$assessment,"_",type[k],".png"), width = width, height = height,
                           res = 200, units = "in")}
      
      if(!add) if(as.png==TRUE | k==1) par(Par)
      
      j = which(c("B","F","BBmsy","FFmsy","BB0","SP")%in%type[k])
      
      
      if(type[k]%in%c("B","F","BBmsy","FFmsy","BB0")){
        y = jbs$timeseries$mu[,j+3]
        years = jbs$timeseries$mu$year
        plot(years,years,type="n",ylim=c(0,max(y[years>=xlim[1] & years<=xlim[2]],ifelse(plotCIs==T,max(jbs$timeseries$uci[,j+3][years>=xlim[1] & years<=xlim[2]]),0))),ylab=ifelse(length(ylabs)>1,ylabs[j],ylabs),xlab="Year",xlim=xlim)
        if(plotCIs==TRUE){
          for(i in 1:length(scenarios)){
            yr = jbs$timeseries$uci[runs%in%scenarios[i],"year"]
            ylc = jbs$timeseries$lci[runs%in%scenarios[i],j+3]
            yuc = jbs$timeseries$uci[runs%in%scenarios[i],j+3]
            polygon(c(yr,rev(yr)),c(ylc,rev(yuc)),col=ifelse(length(scenarios)>1,scales::alpha(cols[i],0.2),"grey"),border=ifelse(length(scenarios)>1,scales::alpha(cols[i],0.2),"grey"))
          }}
        for(i in 1:length(scenarios)){
          yr = jbs$timeseries$uci[runs%in%scenarios[i],"year"]
          lines(yr,y[runs%in%scenarios[i]],col= cols[i],lwd=2,lty=1)
        }
        if(type[k]%in%c("BBmsy","FFmsy")) abline(h=1,lty=2)
      }  else {
        # Plot SP
        plot(years,years,type="n",ylim=c(0,max(jbs$pfunc$SP)),xlim=c(0,max(jbs$pfunc$SB_i)),ylab=ifelse(length(ylabs)>1,ylabs[j],ylabs),xlab=ylabs[1])
        for(i in 1:length(scenarios)){
          lines(jbs$pfunc$SB_i[jbs$pfunc$level%in%scenarios[i]],jbs$pfunc$SP[jbs$pfunc$level%in%scenarios[i]],col=cols[i],lwd=2,lty=1)
          points(mean(jbs$pfunc$SB_i[jbs$pfunc$level%in%scenarios[i]][jbs$pfunc$SP[jbs$pfunc$level%in%scenarios[i]]==max(jbs$pfunc$SP[jbs$pfunc$level%in%scenarios[i]])]),max(jbs$pfunc$SP[jbs$pfunc$level%in%scenarios[i]]),col=cols[i],pch=16,cex=1.2)
        }
      }
      if(legend.add==TRUE){
      if(single.plots==TRUE | k==1 )  legend("topright",paste(scenarios),col=cols,bty="n",cex=legend.cex,pt.cex=0.7,lwd=c(2,rep(1,length(scenarios))))
      }
      if(as.png==TRUE) dev.off()
    } # End type loop
  } else { # Multi plot
    if(is.null(width)) width = 7
    if(is.null(height)) height = 8 
    Par = list(mfrow=c(3,2),mai=c(0.45,0.49,0.1,.15),omi = c(0.15,0.15,0.1,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=plot.cex)
    if(as.png==TRUE){png(file = paste0(output.dir,"/",prefix,"_",jbs$assessment,".png"), width = width, height = height,
                         res = 200, units = "in")}
    par(Par)
    for(k in 1:length(type)){
      
      j = which(c("B","F","BBmsy","FFmsy","BB0","SP")%in%type[k])
      
      
      if(type[k]%in%c("B","F","BBmsy","FFmsy","BB0")){
        years= jbs$timeseries$mu$year
        y = jbs$timeseries$mu[,j+3]
        plot(years,years,type="n",ylim=c(0,max(y[years>=xlim[1] & years<=xlim[2]],ifelse(plotCIs==T,max(jbs$timeseries$uci[,j+3][years>=xlim[1] & years<=xlim[2]]),0))),ylab=ylabs[j],xlab="Year",xlim=xlim)
        if(plotCIs==TRUE){
          for(i in 1:length(scenarios)){
            ylc = jbs$timeseries$lci[runs%in%scenarios[i],j+3]
            yuc = jbs$timeseries$uci[runs%in%scenarios[i],j+3]
            yr = jbs$timeseries$uci[runs%in%scenarios[i],"year"]
            polygon(c(yr,rev(yr)),c(ylc,rev(yuc)),col=ifelse(length(scenarios)>1,scales::alpha(cols[i],0.2),"grey"),border=ifelse(length(scenarios)>1,scales::alpha(cols[i],0.2),"grey"))
          }}
        for(i in 1:length(scenarios)){
          yr = jbs$timeseries$uci[runs%in%scenarios[i],"year"]
          lines(yr,y[runs%in%scenarios[i]],col= cols[i],lwd=2,lty=1)
        }
        if(type[k]%in%c("BBmsy","FFmsy")) abline(h=1,lty=2)
        if(single.plots==TRUE | k==1 )  legend(legend.loc,paste(scenarios),col=cols,bty="n",cex=legend.cex,pt.cex=0.7,lwd=c(2,rep(1,length(scenarios))))
        
      }  else {
        # Plot SP
        plot(years,years,type="n",ylim=c(0,max(jbs$pfunc$SP)),xlim=c(0,max(jbs$pfunc$SB_i)),ylab=ylabs[j],xlab=ylabs[1])
        for(i in 1:length(scenarios)){
          lines(jbs$pfunc$SB_i[jbs$pfunc$level%in%scenarios[i]],jbs$pfunc$SP[jbs$pfunc$level%in%scenarios[i]],col=cols[i],lwd=2,lty=1)
          points(mean(jbs$pfunc$SB_i[jbs$pfunc$level%in%scenarios[i]][jbs$pfunc$SP[jbs$pfunc$level%in%scenarios[i]]==max(jbs$pfunc$SP[jbs$pfunc$level%in%scenarios[i]])]),max(jbs$pfunc$SP[jbs$pfunc$level%in%scenarios[i]]),col=cols[i],pch=16,cex=1.2)
        }
      }
      
      
    }
    if(as.png==TRUE) dev.off()
  }
} # end of Summary Plot


#' jbplot_hcxval
#' 
#' Plots and summarizes results from one step head hindcast cross-validation using the output form jabba_hindcast 
#'
#' @param hc output object from jabba_hindcast
#' @param naive.min minimum MASE denominator (naive predictions) for dj (default = 0.1)
#' @param mase.adj if TRUE, show adjusted mase in brackets
#' @param index option to plot specific indices (numeric & in order)
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param single.plots if TRUE plot invidual fits else make multiplot
#' @param add if TRUE plots par is only called for first plot
#' @param width plot width
#' @param height plot hight
#' @param minyr minimum year shown in plot 
#' @param cols option to add colour palette 
#' @param legend.loc location of legend
#' @param legend.cex size of legend
#' @param legend.add show legend
#' @param label.add show index name and MASE
#' @param label.cex label siz
#' @param ymax upper ylim scaler
#' @param verbose if FALSE then silent
#' @return hcxval statistics by index: MASE, MAE.PR predition residuals,MAE.base for random walk, n.eval obs evaluated 
#' @export
jbplot_hcxval <- function(hc,index=NULL,naive.min=0.1,mase.adj=FALSE, output.dir=getwd(),as.png=FALSE,single.plots=add,add=FALSE,width=NULL,height=NULL,minyr=NULL,cols=NULL,legend.loc="topright",legend.cex=0.7,legend.add=TRUE,label.add=TRUE,label.cex=0.9,verbose=TRUE,ymax=1){
  if(as.png==TRUE) add= FALSE
  MASE = jbmase(hc,naive.min=naive.min,verbose=verbose,residuals=FALSE)
  if(is.null(cols)) cols = ss3col(length(hc)-1)
  d. = do.call(rbind,lapply(hc,function(x){
    x$diags}))
  
  
  all.indices = unique(d.$name)  
  if(is.null(index)) index = 1:length(all.indices)
  # subset 
  d. = d.[d.$name%in%all.indices[index],]
  d. = d.[d.$name%in%MASE[is.na(MASE$MASE),]$Index==F,]
  
  peels = unique(d.$retro.peels)
  styr = max(hc[[1]]$yr)-max(peels)
  years = min(d.$year):max(d.$year)
  yr = unique(d.$year)
  endyrvec = rev(sort(years[length(years)-peels]))
  
  if(is.null(minyr)){
    xmin = min(hc[[1]]$yr)} else {
      xmin = min(minyr,min(endyrvec)-3)  
    }
  d.=d.[d.$year>=xmin,]
  
  # check in index
  indices = unique(d.$name)
  
  n.indices = length(indices)  
  if(single.plots==FALSE){  
    Par = list(mfrow=c(round(n.indices/2+0.01,0),ifelse(n.indices==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
    if(as.png==TRUE){png(file = paste0(output.dir,"/hcxaval_",hc[[1]]$assessment,".png"), width = 7, height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0),
                         res = 200, units = "in")}
    if(add==F) par(Par)
  }  
  for(i in 1:length(indices)){
    if(nrow(d.[d.$name%in%indices[i] & d.$year>styr & d.$retro.peels%in%peels[1],])>1){ # Only run if overlap
      if(single.plots==TRUE){  
        if(is.null(width)) width = 5
        if(is.null(height)) height = 3.5
        
        Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
        if(as.png==TRUE){png(file = paste0(output.dir,"/hcxval_",hc$scenario,"_",indices[i],".png"), width = width, height = height,
                             res = 200, units = "in")}
        if(add==F){
          if(as.png==TRUE | indices[i]==valid[1]) par(Par)
        }
      }
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
      
      #if(length(endyrvec[yr%in%endyrvec])>0 & length(which(yr.eval%in%yr))>1){ # ><>
      
      
      py = xv$year[xv$retro.peels==min(xv$retro.peels) & xv$year>styr-xmin]
      obs =xv$obs[xv$retro.peels==min(xv$retro.peels) & xv$year>styr-xmin]
      hat = xv$hat[xv$retro.peels==min(xv$retro.peels) & xv$year>styr-xmin]
      lc = xv$hat.lci[xv$retro.peels==min(xv$retro.peels) & xv$year>styr-xmin]
      uc = xv$hat.uci[xv$retro.peels==min(xv$retro.peels) & xv$year>styr-xmin]
      plot(0, type = "n", xlim = c(max(min(yr),min(endyrvec-xmin+1)),min(c(max(yr),max(endyrvec)))), yaxs = "i", 
           ylim = c(ifelse(min(c(lc,obs))*0.5<0.5,0,min(c(lc,obs))*0.5),max(c(uc,obs)*1.25*ymax)), xlab = "Year", ylab = "Index")
      
      polygon(c(py,rev(py)),c(lc,rev(uc)),col=grey(0.5,0.4),border=grey(0.5,0.4))
      polygon(c(py[py<=min(endyrvec)],rev(py[py<=min(endyrvec)])),c(lc[py<=min(endyrvec)],rev(uc[py<=min(endyrvec)])),col=grey(0.4,0.4),border=grey(0.4,0.4))
      
      #lines(py,obs,pch=21,lty=2,col="white")
      points(py,obs,pch=21,cex=1.6,bg="white")
      
      naive.eval=log(obs.eval[is.na(obs.eval)==F][-length(obs.eval[is.na(obs.eval)==F])])-log(obs.eval[is.na(obs.eval)==F][-1])
      nhc = length(endyrvec)-1
      points(rev(yr.eval[-1])[1:length( naive.eval)][is.na(naive.eval)==F],rev(obs.eval[-1])[1:length( naive.eval)][is.na(naive.eval)==F],pch=21,cex=1.6,bg=((cols[1:(length(peels)-1)]))[is.na(naive.eval)==F])
      
      
      pred.resid = NULL
      for(j in 1:(length(peels)-1)){
        if(endyrvec[1:length( naive.eval)][j] %in% xv$year){
          
          x <- min(py):max(yr.eval)
          x <- x[1:(length(x)-peels[j])]
          x = x[x%in%unique(xv$year)]
          y <- xv[xv$retro.peels==peels[j+1] & xv$year%in%x,]$hat
          pred.resid = c(pred.resid,log(y[length(x)])-log(obs[length(x)])) # add log() for v1.1
          lines(x, y, lwd=2,
                lty=1, col=cols[j], type="l",cex=0.9)
          
          lines(x[(length(x)-1):(length(x))], y[(length(x)-1):(length(x))], lwd=2,
                col=1,lty=2)
          
          points(x[length(x)], y[length(y)],pch=21,
                 bg=cols[j],col=1, type="p",cex=1)
          
        }}
      lines(py,hat,col=1,lwd=2,lty=1,type="l",pch=16)
      
      mase = MASE[MASE$Index%in%indices[i],]$MASE
      maseadj = MASE[MASE$Index%in%indices[i],]$MASE.adj
      lmase = paste0(unique(xv$name)[1], ": MASE = ",round(mase,2))
      if(mase.adj==TRUE) lmase = paste0(lmase,"(",round(maseadj,2),")")
      
      if(label.add) legend("top",lmase,bty="n",y.intersp=-0.2,cex=label.cex)
      
      if(single.plots==TRUE & as.png==TRUE) dev.off()
      if(legend.add==T){
        if(single.plots==TRUE | i==1){ legend(legend.loc,paste(c("Ref",paste0("-",max(hc[[1]]$yr)-peels[-1]+1))),col=c(1,cols),bty="n",cex=legend.cex,lwd=2)}
      }
      
    } else{
      xv = d.[d.$name%in%indices[i],]
    }
    
  } # end of index loop
  if(single.plots==FALSE){
    mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1.)
    mtext(paste("Index"), side=2, outer=TRUE, at=0.5,line=1,cex=1)  
  }
  if(single.plots==FALSE & as.png==TRUE) dev.off()
  
  if(verbose) return(MASE)
}
#}}}

#' jrplot_PPC
#'
#' Plots Posterior Predictive Checks for Chisq Goodness-of-Fit (Omnibus Discrepancy) 
#' @param jabba output list from fit_jabba
#' @param joint.ppc FALSE/TRUE, if TRUE a joint Posterior Predictive Check (indices combined) is plotted
#' @param thin.plot  TRUE/FALSE thinning option for plotting
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param single.plots if TRUE plot invidual fits else make multiplot
#' @param width plot width
#' @param height plot hight
#' @param ylab option to change y-axis label
#' @param xlab option to change x-axis label
#' @param plot.cex graphic option
#' @param indices names of indices to plot (default = "all")
#' @param index.label show index name in plot
#' @param add if TRUE par is not called to enable manual multiplots
#' @param verbose commentary 
#' @return Bayesin p values and n observation per index
#' @export
jbplot_PPC <- function(jabba,joint.ppc=FALSE,thin.plot = TRUE ,output.dir=getwd(),as.png=FALSE,single.plots=add,width=NULL,height=NULL,ylab=expression(paste("Predicted D(",chi^2,")")),xlab=expression(paste("Realized D(",chi^2,")")),plot.cex=0.8,index=NULL,index.label=TRUE,add=FALSE,verbose=TRUE){
  
  if(verbose) cat(paste0("\n","><> jrplot_PPC() - Posterior Predictive Checks <><","\n"))
  if(is.null(jabba$PPC)) stop("To enable PPCs please use option fit_jabba(jabbainput,do.ppc=TRUE)") 
  PPC = jabba$PPC
  Bp = NULL # Bayesian p-value
  name.check = jabba$diags$name
  if(joint.ppc){
    PPC$name ="Combined"
    indices = "Combined"
    n.indices = 1
    
  } else {
    all.indices = unique(jabba$diags$name)
    if(is.null(index)) index = 1:length(all.indices)
    indices = unique(jabba$diags$name)[index]
    n.indices = length(indices)

    
    for(i in 1:n.indices){
      ppc = PPC[PPC$name%in%indices[i],2:3]
      Bp = rbind(Bp,data.frame(Index=indices[i],Bayesian.p=sum(ppc$Dy>ppc$Dx)/length(ppc$Dx),nobs=nrow(jabba$diags[jabba$diags$name==indices[i],]))) 
    }}
  ppc = PPC[PPC$name%in%indices,2:3]
  Bp = rbind(Bp,data.frame(Index="Combined",Bayesian.p=sum(ppc$Dy>ppc$Dx)/length(ppc$Dx),nobs=nrow(jabba$diags[name.check%in%indices,]))) 
  
  
  if(single.plots==TRUE){
    if(is.null(width)) width = 5
    if(is.null(height)) height = 5
    for(i in 1:n.indices){
      Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=plot.cex)
      if(as.png==TRUE){png(file = paste0(output.dir,"/Fits",jabba$assessment,"_",jabba$scenario,"_",indices[i],".png"), width = width, height = height,
                           res = 200, units = "in")}
      if(add==FALSE){
        if(as.png==TRUE | i==1) par(Par)
      }  
      # set observed vs predicted CPUE
      ppc = PPC[PPC$name%in%indices[i],2:3]
      nmc = nrow(ppc)
      if(thin.plot){
        thinning = seq(1,nmc,floor(nmc/1000*0.5))
        ppct = ppc[thinning,]
      } else {
        ppct=ppc
      }
      
      # Plot Observed vs predicted CPUE
      plot(ppct[,1],ppct[,2],ylab="",xlab="",ylim=c(0,max(ppc)*1.05),xlim=c(0,max(ppc)*1.05),type='n',xaxt="n",yaxt="n")
      axis(1,labels=TRUE,cex=0.8)
      axis(2,labels=TRUE,cex=0.8)
      points(ppct,pch=21,col=grey(0.4,0.7),bg=grey(0.7,0.5),cex=0.9)
      abline(0,1,lwd=1.5)
      if(index.label==TRUE)legend('top',paste(indices[i],": p =",round(Bp$Bayesian.p[i],3)),bty="n",y.intersp = -0.2,cex=0.9)
      mtext(xlab, side=1, outer=TRUE, at=0.5,line=0.5,cex=1)
      mtext(ylab, side=2, outer=TRUE, at=0.5,line=0.3,cex=1)
      if(as.png==TRUE) dev.off()
    }
  } else {
    
    if(is.null(width)) width = 7
    if(is.null(height)) height = ifelse(n.indices==1,7,ifelse(n.indices==2,4.,3.5))*round(n.indices/2+0.01,0)
    Par = list(mfrow=c(round(n.indices/2+0.01,0),ifelse(n.indices==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=plot.cex)
    if(as.png==TRUE){png(file = paste0(output.dir,"/Fits_",jabba$assessment,"_",jabba$scenario,".png"), width = 7, height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0),
                         res = 200, units = "in")}
    par(Par)
    
    for(i in 1:n.indices){
      # set observed vs predicted CPUE
      # set observed vs predicted CPUE
      ppc = PPC[PPC$name%in%indices[i],2:3]
      nmc = nrow(ppc)
      if(thin.plot){
        thinning = seq(1,nmc,ceiling(nmc/1000*0.2))
        ppct = ppc[thinning,]
      } else {
        ppct=ppc
      }
      
      # Plot Observed vs predicted CPUE
      plot(ppct[,1],ppct[,2],ylab="",xlab="",ylim=c(0,max(ppc)*1.05),xlim=c(0,max(ppc)*1.05),type='n',xaxt="n",yaxt="n")
      axis(1,labels=TRUE,cex=0.8)
      axis(2,labels=TRUE,cex=0.8)
      points(ppct,pch=21,col=grey(0.4,0.7),bg=grey(0.7,0.5),cex=0.9)
      abline(0,1,lwd=1.5)
      if(index.label==TRUE)legend('top',paste(indices[i],": p =",round(Bp$Bayesian.p[i],3)),bty="n",y.intersp = -0.2,cex=0.9)
      
      
    }
    if(add==FALSE | i==1){
      mtext(xlab, side=1, outer=TRUE, at=0.5,line=0.75,cex=1)
      mtext(ylab, side=2, outer=TRUE, at=0.5,line=1,cex=1)
    }
    if(as.png==TRUE){dev.off()}
  }
  if(verbose) cat("\n","Posterior Predictive Checks with Bayesian p values","\n")
  return(Bp)
} # End of CPUE plot function



