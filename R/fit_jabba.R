#' Fit JABBA model
#'
#' Fits JABBA model in JAGS and produce output object as list()
#' @param jbinput List of input variables as output by build_jabba()
#' MCMC settings
#' @param ni number of iterations
#' @param nt thinning interval of saved iterations
#' @param nb burn-in
#' @param nc number of mcmc chains
#' Initial values
#' @param init.values = FALSE,
#' @param init.K = NULL,
#' @param init.r = NULL,
#' @param init.q = NULL,# vector
#' @param peels = NULL, # retro peel option
#' @param do.ppc conducts and saves posterior predictive checks
#' @param save.trj adds posteriors of stock, harvest and bk trajectories
#' @param save.csvs option to write csv outputs
#' @param save.jabba saves jabba fit as rdata object
#' @param save.all add complete posteriors to fitted object 
#' @param output.dir path to save plot. default is getwd()
#' @param quickmcmc option to run short mcmc
#' @param seed default 123, set random by e.g. sample.int(999,1) 
#' @param verbose option show cat comments and progress
#' @return A result list containing estimates of model input, settings and results
#' @export
#' @examples
#' data(iccat)
#' jbinput <- build_jabba(catch=iccat$bet$catch,cpue=iccat$bet$cpue,se=iccat$bet$se,model.type="Fox")
#' bet1 = fit_jabba(jbinput,quickmcmc=TRUE,verbose=TRUE)
#' jbplot_summary(bet1)

fit_jabba = function(jbinput,
                     # MCMC settings
                     ni = 30000, # Number of iterations
                     nt = 5, # Steps saved
                     nb = 5000, # Burn-in
                     nc = 2, # number of chains
                     # init values
                     init.values = FALSE,
                     init.K = NULL,
                     init.r = NULL,
                     init.q = NULL,# vector
                     peels = NULL, # retro peel option
                     do.ppc = TRUE,
                     save.trj = TRUE,
                     save.all = FALSE,
                     save.jabba = FALSE,
                     save.csvs = FALSE,
                     output.dir = getwd(),
                     quickmcmc = FALSE,
                     seed = 123,
                     verbose=TRUE
){
  
  tmpath <- tempfile()
  dir.create(tmpath)
  if(!verbose) {
    progress.bar="none"
  } else {
    progress.bar="text"
  }
  
  set.seed(seed)
  
  if(jbinput$settings$projection==TRUE) save.trj=TRUE
  
  #write jabba model
  jabba2jags(jbinput, tmpath)
  
  # mcmc saved
  nsaved = (ni-nb)/nt*nc
  # jabba model data
  jbd = jbinput$jagsdata
  
  if(quickmcmc==TRUE){
    ni = 8000
    nb = 2000
    nt = 2
    nc = 2
  }
  
  # a, b pars for beta prior
  ab = get_beta(max(min(jbinput$settings$psi.prior.raw[1],0.95),0.05),CV=0.05/max(min(jbinput$settings$psi.prior.raw[1],0.95),0.05))
  
  qbound = jbinput$jagsdata$q_bounds
  
  # Initial starting values (new Eq)
  if(init.values==FALSE){
    inits = function(){list(K= rlnorm(1,log(jbd$K.pr[1])-0.5*0.3^2,0.3),r = rlnorm(1,log(jbd$r.pr[1]),jbd$r.pr[2]) ,
                            q = pmin(pmax(qbound[1]*1.05,runif(jbd$nq,min(jbd$I,na.rm=T)/max(jbd$TC,na.rm=T),mean(jbd$I,na.rm=T)/max(jbd$TC,na.rm=T))),qbound[2]*0.95)
                            ,psi=rbeta(1,ab[1],ab[2]),isigma2.est=runif(1,20,100), itau2=runif(jbd$nvar,80,200))}
  }else {
    
    if(is.null(init.K))
      stop("\n","\n","><> Provide init.K guess for option init.values=TRUE  <><","\n","\n")
    if(is.null(init.r))
      stop("\n","\n","><> Provide init.r guess for option init.values=TRUE  <><","\n","\n")
    if(is.null(init.q))
      stop("\n","\n","><> Provide init.q vector guess for option init.values=TRUE  <><","\n","\n")
    if(length(init.q)!= jbinput$jagsdata$nq)
      stop("\n","\n","><> init.q vector must match length of estimable q's, length(unique(sets.q))   <><","\n","\n")
    inits = function(){ list(K= init.K,r=init.r,q = init.q,psi=rbeta(1,ab[1],ab[2]), isigma2.est=runif(1,20,100), itau2=runif(jbd$nvar,80,200))}
  }
  
  
  
  out = output.dir
  if(file.exists(out)==FALSE) stop("\n","\n","><> output.dir does not exist <><","\n","\n")
  
  
  # retrospecitive peel
  years = jbinput$data$yr
  if(is.null(peels)) peels = 0
  if(peels > 0){
    jbd$I[(length(years)-peels+1) : length(years),]  = NA
    
    if(jbinput$settings$Auxiliary){
    jbd$A[(length(years)-peels+1) : length(years),]  = NA
    }
    
  }
  jbinput$jagsdata$I = jbd$I 
  jbinput$jagsdata$A = jbd$A 
  # update
  # jabba model building conditions
  params = jbinput$settings$params
  
  
  ptm <- proc.time()
  
  mod <- R2jags::jags(jbd, inits,params,file.path(tmpath,"JABBA.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, quiet=!verbose, progress.bar = progress.bar)  # adapt is burn-in
  
  proc.time() - ptm
  save.time = proc.time() - ptm
  
  # unpack
  settings= c(jbinput$data,jbinput$jagsdata,jbinput$settings)
  catch = settings$catch
  cpue= settings$cpue
  auxiliary = settings$auxiliary
  se = settings$se
  n.years = settings$N
  years = settings$yr
  assessment = settings$assessment
  scenario = settings$scenario
  CPUE = settings$I
  if(settings$Auxiliary) AUXI = settings$A
  n.indices = settings$nI
  
  
  # if run with library(rjags)
  posteriors = mod$BUGSoutput$sims.list
  if(verbose)
    message(paste0("\n","><> Produce results output of ",settings$model.type," model for ",settings$assessment," ",settings$scenario," <><","\n"))
  
  
  
  #-----------------------------------------------------------
  # <><<><<><<><<><<><<>< Outputs ><>><>><>><>><>><>><>><>><>
  #-----------------------------------------------------------
  
  # run some mcmc convergence tests
  par.dat= data.frame(posteriors[params[c(1:7)]])
  
  geweke = coda::geweke.diag(data.frame(par.dat))
  pvalues <- 2*pnorm(-abs(geweke$z))
  pvalues
  
  heidle = coda::heidel.diag(data.frame(par.dat))
  
  # postrior means + 95% BCIs
  man.dat = data.frame(posteriors[params[8:10]],bmsyk=as.numeric(posteriors$SBmsy)/as.numeric(posteriors$K))
  
  # Depletion
  Depletion = posteriors$P[,c(1,n.years)]
  colnames(Depletion) = c(paste0("P",years[1]),paste0("P",years[n.years]))
  
  # Current stock status (Kobe posterior)
  H_Hmsy.cur = posteriors$HtoHmsy[,c(n.years)]
  B_Bmsy.cur = posteriors$BtoBmsy[,c(n.years)]
  
  
  # Prepare posterior quantaties
  man.dat = data.frame(man.dat,Depletion,B_Bmsy.cur,H_Hmsy.cur)
  
  results = t(cbind(apply(par.dat,2,quantile,c(0.025,0.5,0.975)))) 
  
  results = data.frame(Median = results[,2],LCI=results[,1],UCI=results[,3],Geweke.p=round(pvalues,3),Heidel.p = round(heidle[,3],3))
  
  ref.points = round(t(cbind(apply(man.dat,2,quantile,c(0.025,0.5,0.975)))),3)
  
  ref.points = data.frame(Median = ref.points[,2],LCI=ref.points[,1],UCI=ref.points[,3])
  # get number of parameters
  npar = length(par.dat)
  # number of years
  N=n.years
  
  #-------------------------------------------------------------------------
  # Save parameters, results table and current status posterior in csv files
  #-------------------------------------------------------------------------
  
  # Make standard results table with parameter estimates and reference points
  Table = rbind(data.frame(results)[c("K","r","psi","sigma2","m"),1:3],data.frame(ref.points))
  Table[4,] = round(sqrt((Table[4,])),3)
  rownames(Table)[4] = "sigma.proc"
  colnames(Table)  <- c("mu","lci","uci")
  
  #-----------------------------------------------
  # Stock trajectories
  #-----------------------------------------------
  #Bt = posteriors$SB
  #Ht = posteriors$H
  #Bt_Bmsy = posteriors$BtoBmsy
  #Ht_Hmsy = posteriors$HtoHmsy
  #Bt_K = posteriors$P
  catch.temp = matrix(rep(catch[,2],each=nrow(posteriors$SB)),ncol=nrow(jbinput$data$catch),nrow=nrow(posteriors$SB))
  
  yrdim = length(years)
  Stock_trj = array(data=NA,dim=c(yrdim,3,7),dimnames = list(years,c("mu","lci","uci"),c("B","F","BBmsy","FFmsy","BB0","procB","SPt")))
  for(i in 1:3){
    Stock_trj[,i,] =  cbind(t(apply(posteriors$SB[,1:yrdim],2,quantile,c(0.5,0.025,0.975)))[,i],
                            t(apply(posteriors$H[,1:yrdim],2,quantile,c(0.5,0.025,0.975)))[,i],
                            t(apply(posteriors$BtoBmsy[,1:yrdim],2,quantile,c(0.5,0.025,0.975)))[,i],
                            t(apply(posteriors$HtoHmsy[,1:yrdim],2,quantile,c(0.5,0.025,0.975)))[,i],
                            t(apply(posteriors$P[,1:yrdim],2,quantile,c(0.5,0.025,0.975)))[,i],
                            t(apply(posteriors$Proc.Dev[,1:yrdim],2,quantile,c(0.5,0.025,0.975)))[,i],
                            t(cbind(rep(0,3),apply(posteriors$SB[,-1]-posteriors$SB[,-ncol(posteriors$SB)]+catch.temp[,-ncol(posteriors$SB)],2,quantile,c(0.5,0.025,0.975)))[,1:yrdim])[,i]
    )
    
  }
  
  
  
  #------------------------
  # Production function
  #------------------------
  m.sp = median(posteriors$m)
  Bit = seq(1,median(posteriors$K),median(posteriors$K)/500)
  Cmsy = Bit*median(posteriors$Hmsy)
  B.sp = apply(posteriors$SB,2,mean)
  Hmsy.sp = median(posteriors$Hmsy)
  SB0.sp = median(posteriors$K)
  SP = Hmsy.sp/(1-1/m.sp)*Bit*(1-(Bit/SB0.sp)^(m.sp-1))
  Bmsy.sp = median(posteriors$SBmsy)
  MSY.sp = quantile(posteriors$SBmsy*posteriors$Hmsy,c(0.025,0.5,0.975))
  
  #------------------
  # Goodness-of-Fit
  #------------------
  DIC =round(mod$BUGSoutput$DIC,1)
  if(settings$Auxiliary) nA = ncol(AUXI)
  
  if(settings$CatchOnly==FALSE | settings$Auxiliary==TRUE ){
    
    # get residuals
    Resids = NULL
    StResid = NULL
    
    if(settings$CatchOnly==FALSE){
    for(i in 1:n.indices){
      Resids =rbind(Resids,log(CPUE[,i])-log(apply(posteriors$CPUE[,,i],2,quantile,c(0.5))))
    
      StResid =rbind(StResid,log(CPUE[,i]/apply(posteriors$CPUE[,,i],2,quantile,c(0.5)))/
                       apply(posteriors$TOE[,,i],2,quantile,c(0.5))+0.5*apply(posteriors$TOE[,,i],2,quantile,c(0.5)))
    }
    }
    if(settings$Auxiliary==TRUE){
      for(i in 1:nA){
        Resids =rbind(Resids,log(AUXI[,i])-log(apply(posteriors$Ahat[,,i],2,quantile,c(0.5))))
        StResid =rbind(StResid,log(AUXI[,i]/apply(posteriors$Ahat[,,i],2,quantile,c(0.5)))/
                         apply(posteriors$TAE[,,i],2,quantile,c(0.5))+0.5*apply(posteriors$TAE[,,i],2,quantile,c(0.5)))
        
       }
    }  
       
    
    Nobs =length(as.numeric(Resids)[is.na(as.numeric(Resids))==FALSE])
    DF = Nobs-npar
    RMSE = round(100*sqrt(sum(Resids^2,na.rm =TRUE)/DF),1)
    SDNR = round(sqrt(sum(StResid^2,na.rm =TRUE)/(Nobs-1)),2)
    Crit.value = (qchisq(.95, df=(Nobs-1))/(Nobs-1))^0.5
    # Produce statistice describing the Goodness of the Fit
  } else {
    
    Nobs =DF = RMSE = SDNR = Crit.value = NA
    
  }
  
  # Save Obs,Fit,Residuals
  jabba.res = PPC = NULL
  if(settings$CatchOnly==FALSE){
      
      Yr = years
      Yr = min(Yr):max(Yr)
      yr = Yr-min(years)+1
    
    for(i in 1:n.indices){
      exp.i = apply(posteriors$CPUE[,,i],2,quantile,c(0.5))[is.na(cpue[,i+1])==F]
      se.i = apply(posteriors$CPUE[,,i],2,sd)[is.na(cpue[,i+1])==F]
      expLCI.i = apply(posteriors$Ihat[,,i],2,quantile,c(0.025))[is.na(cpue[,i+1])==F]
      expUCI.i = apply(posteriors$Ihat[,,i],2,quantile,c(0.975))[is.na(cpue[,i+1])==F]
      obs.i = cpue[is.na(cpue[,i+1])==F,i+1]
      sigma.obs.i = (apply(posteriors$TOE[,,i],2,quantile,c(0.5)))[is.na(cpue[,i+1])==F]
      yr.i = Yr[is.na(cpue[,i+1])==F]
      jabba.res = rbind(jabba.res,data.frame(scenario=settings$scenario,name=names(cpue)[i+1],year=yr.i,obs=obs.i,obs.err=sigma.obs.i,hat=exp.i,hat.lci=expLCI.i,hat.uci=expUCI.i,residual=log(obs.i)-log(exp.i),retro.peels=peels))
    
      if(do.ppc){ # Posterior Predictive Check
        nmc = nrow(posteriors$K)
        mcmcfit = posteriors$Ihat[,is.na(cpue[,i+1])==F,i]
        mcmcppd = posteriors$CPUE[,is.na(cpue[,i+1])==F,i]

        observed = obs.i  
        observed = array(rep(observed,each=nmc),c(nmc,length(observed)))
        sigobs =  posteriors$TOE[,is.na(cpue[,i+1])==F,i]
        PPC = rbind(PPC,data.frame(name=names(cpue)[i+1],Dx=apply((log(observed)-log(mcmcfit))^2/sigobs,1,sum),   
                                   Dy = apply((log(mcmcppd)-log(mcmcfit))^2/sigobs,1,sum)))
      }
      
      
      
      }
  }
  
  if(settings$Auxiliary==TRUE){
      
      Yr = years
      Yr = min(Yr):max(Yr)
      yr = Yr-min(years)+1
      
      for(i in 1:nA){
      
      exp.i = apply(posteriors$Ahat[,,i],2,quantile,c(0.5))[is.na(auxiliary[,i+1])==F]
      se.i = apply(posteriors$Ahat[,,i],2,sd)[is.na(auxiliary[,i+1])==F]
      expLCI.i = apply(posteriors$Ahat[,,i],2,quantile,c(0.025))[is.na(auxiliary[,i+1])==F]
      expUCI.i = apply(posteriors$Ahat[,,i],2,quantile,c(0.975))[is.na(auxiliary[,i+1])==F]
      
      obs.i = auxiliary[is.na(auxiliary[,i+1])==F,i+1]
      sigma.obs.i = (apply(posteriors$TAE[,,i],2,quantile,c(0.5)))[is.na(auxiliary[,i+1])==F]
      
      yr.i = Yr[is.na(auxiliary[,i+1])==F]
      jabba.res = rbind(jabba.res,data.frame(scenario=settings$scenario,name=names(auxiliary)[i+1],year=yr.i,obs=obs.i,obs.err=sigma.obs.i,hat=exp.i,hat.lci=expLCI.i,hat.uci=expUCI.i,residual=log(obs.i)-log(exp.i),retro.peels=peels))

      if(do.ppc){ # Posterior Predictive Check
        nmc = nrow(posteriors$K)
        mcmcfit = posteriors$Ahat[,is.na(auxiliary[,i+1])==F,i]
        mcmcppd = posteriors$AUXI[,is.na(auxiliary[,i+1])==F,i]
        
        observed = obs.i  
        observed = array(rep(observed,each=nmc),c(nmc,length(observed)))
        sigobs =  posteriors$TAE[,is.na(auxiliary[,i+1])==F,i]
        PPC = rbind(PPC,data.frame(name=names(auxiliary)[i+1],Dx=apply((log(observed)-log(mcmcfit))^2/sigobs,1,sum),   
                                   Dy = apply((log(mcmcppd)-log(mcmcfit))^2/sigobs,1,sum)))
      }
      
  }
  }
  

  if(!jbinput$settings$Auxiliary){
  #----------------------------------
  # Predicted CPUE
  #----------------------------------
  cpue.hat = array(data=NA,dim=c(N,5,n.indices),list(years,c("mu","lci","uci","se","obserror"),names(cpue[,-1])))
  for(i in 1:n.indices){
    cpue.hat[,,i] = cbind(t(apply(posteriors$Ihat[,,i],2,quantile,c(0.5,0.025,0.975))),apply(log(posteriors$Ihat[,,i]),2,sd),(apply(posteriors$TOE[,,i],2,quantile,c(0.5))))
  }
  #------------------------------------
  # Posterior Predictive Distribution
  #------------------------------------
  
  cpue.ppd = array(data=NA,dim=c(N,5,n.indices),list(years,c("mu","lci","uci","se","obserror"),names(cpue[,-1])))
  for(i in 1:n.indices){
    cpue.ppd[,,i] = cbind(t(apply(posteriors$CPUE[,,i],2,quantile,c(0.5,0.025,0.975))),apply(log(posteriors$CPUE[,,i]),2,sd),(apply(posteriors$TOE[,,i],2,quantile,c(0.5))))
  }
}
  
  
  
    if(jbinput$settings$Auxiliary){
    #----------------------------------
    # Predicted CPUE
    #----------------------------------
    cpue.hat = array(data=NA,dim=c(N,5,n.indices+nA),list(years,c("mu","lci","uci","se","obserror"),names(cpue[,-1])))
    for(i in 1:n.indices){
      cpue.hat[,,i] = cbind(t(apply(posteriors$Ihat[,,i],2,quantile,c(0.5,0.025,0.975))),apply(log(posteriors$Ihat[,,i]),2,sd),(apply(posteriors$TOE[,,i],2,quantile,c(0.5))))
    }
    for(i in 1:nA){
    cpue.hat[,,n.indices+i] = cbind(t(apply(posteriors$Ahat[,,i],2,quantile,c(0.5,0.025,0.975))),apply(log(posteriors$Ahat[,,i]),2,sd),(apply(posteriors$TAE[,,i],2,quantile,c(0.5))))
    }
    #------------------------------------
    # Posterior Predictive Distribution
    #------------------------------------
    
    cpue.ppd = array(data=NA,dim=c(N,5,n.indices+nA),list(years,c("mu","lci","uci","se","obserror"),names(cpue[,-1])))
    for(i in 1:n.indices){
      cpue.ppd[,,i] = cbind(t(apply(posteriors$CPUE[,,i],2,quantile,c(0.5,0.025,0.975))),apply(log(posteriors$CPUE[,,i]),2,sd),(apply(posteriors$TOE[,,i],2,quantile,c(0.5))))
    }
    for(i in 1:nA){
      cpue.ppd[,,n.indices+i] = cbind(t(apply(posteriors$AUXI[,,i],2,quantile,c(0.5,0.025,0.975))),apply(log(posteriors$AUXI[,,i]),2,sd),(apply(posteriors$TAE[,,i],2,quantile,c(0.5))))
    }
  }
  
  
  
  #-----------------------------------
  # Note posteriors of key parameters
  #-----------------------------------
  sel.par = c(1,2,7,4,3,5)
  if(!settings$Auxiliary){
  out=data.frame(posteriors[params[sel.par]])
  } else {
    out=data.frame(posteriors[c(params[sel.par],"qA")])
  }
  outman = man.dat[,1:4]
  colnames(outman) = c("Fmsy","Bmsy","MSY","BmsyK")
  if(verbose)
    message(paste0("\n","\n",paste0("><> Scenario ", jbinput$settings$scenario,"_",jbinput$settings$model.type," completed in ",as.integer(save.time[3]/60)," min and ",round((save.time[3]/60-as.integer(save.time[3]/60))*100)," sec <><","\n")))
  
  
  
  #-------------------------------
  # summarize results in jabba list
  #-------------------------------
  jabba = list()
  jabba$assessment  = assessment
  jabba$scenario = scenario
  jabba$settings = c(jbinput$jagsdata,jbinput$settings)
  jabba$settings$mcmc = list(ni=ni,nt=nt,nb=nb,nc=nc,nsaved=nsaved)
  jabba$inputseries = list(cpue=cpue,se=se,catch=catch)
  jabba$pars=results
  jabba$estimates=Table
  jabba$yr = years
  jabba$catch = settings$TC
  
  if(settings$add.catch.CV ==TRUE){jabba$est.catch = data.frame(yr=years,mu=apply(posteriors$estC,2,quantile,c(0.5,0.025,0.975))[1,],
                                                                lci=apply(posteriors$estC,2,quantile,c(0.5,0.025,0.975))[2,],
                                                                uci=apply(posteriors$estC,2,quantile,c(0.5,0.025,0.975))[3,])} else {
                                                                  jabba$est.catch = "Require option: add.catch.cv = TRUE"
                                                                }
  jabba$cpue.hat=cpue.hat
  jabba$cpue.ppd=cpue.ppd
  jabba$PPC = NULL
  if(do.ppc) jabba$PPC = PPC 
  jabba$timeseries=Stock_trj
  jabba$refpts = data.frame(factor=assessment,level=settings$scenario,quant = c("hat","logse"), k=c(median(posteriors$K),sd(log(posteriors$K))),bmsy=c(median(posteriors$SBmsy),sd(log(posteriors$SBmsy))),
                            fmsy=c(median(posteriors$Hmsy),sd(log(posteriors$Hmsy))),msy=c(median(posteriors$SBmsy*posteriors$Hmsy),sd(log(posteriors$SBmsy*posteriors$Hmsy))))
  jabba$pfunc = data.frame(factor=assessment,level=scenario,SB_i=round(Bit,3),SP=round(SP,3),Hmsy=round(Hmsy.sp,4),r=round(Hmsy.sp*(m.sp-1)/(1-1/m.sp),4),m=round(m.sp,3),MSY=round(as.numeric(MSY.sp[2]),3),SB0=round(SB0.sp,3),Cmsy=round(Cmsy,3))
  
  if(settings$CatchOnly==FALSE | settings$Auxiliary==TRUE){
    jabba$diags = data.frame(factor=assessment,level=settings$scenario,name=jabba.res$name,year=jabba.res$year,season=1,obs=jabba.res$obs,hat=jabba.res$hat,hat.lci=jabba.res$hat.lci,hat.uci=jabba.res$hat.uci,residual=jabba.res$residual,retro.peels=jabba.res$retro.peels)
    # add hindcast qualifier
    if(jabba$diags$retro.peels[1]>0){
      jabba$diags$hindcast =  ifelse(jabba$diags$year%in%(max(jabba$yr)-(1:(jabba$diags$retro.peels[1]))+1),TRUE,FALSE)
    } else {
      jabba$diags$hindcast=FALSE
    }
  } else {
    jabba$diags = "Not Available for Catch-Only option"
    
  } 
    
  if(settings$CatchOnly==FALSE | settings$Auxiliary==TRUE){
    
    jabba$residuals = array(Resids,dim=c(nrow(Resids),ncol(Resids)),dimnames = list(unique(jabba.res$name),years))
    jabba$std.residuals = array(StResid,dim=c(nrow(StResid),ncol(StResid)),dimnames = list(unique(jabba.res$name),years))
    
 } else {
    jabba$residuals = "Not Available for Catch-Only option"
  }
  jabba$stats = data.frame(Stastistic = c("N","p","DF","SDNR","RMSE","DIC"),Value = c(Nobs,npar,DF,SDNR,RMSE,DIC))
  jabba$pars_posterior = out
  jabba$refpts_posterior = outman
  
  jabba$kobe = data.frame(factor=assessment,level=scenario,yr=years[N],stock=posteriors$BtoBmsy[,N],harvest=posteriors$HtoHmsy[,N],bk=posteriors$P[,N])
  
  # produce flqs with on step ahead bio
  if(jbinput$settings$add.catch.CV ==TRUE){
    flcatch = data.frame(age="all",year= c(years,max(years+1)),unit="unique",season="all",area="unique",iter=1,data= c(apply(posteriors$estC,2,quantile,0.5),NA),qname="catch")
  } else {
    flcatch = data.frame(age="all",year= c(years,max(years+1)),unit="unique",season="all",area="unique",iter=1,data= c(jbinput$jagsdata$TC,NA),qname="catch")
  }
  flqs = rbind( 
    data.frame(age="all",year= c(years,max(years+1)),unit="unique",season="all",area="unique",iter=1,data= apply(posteriors$SB,2,quantile,0.5),qname="biomass"),
    flcatch,
    data.frame(age="all",year= c(years,max(years+1)),unit="unique",season="all",area="unique",iter=1,data= apply(posteriors$BtoBmsy,2,quantile,0.5),qname="stock"),
    data.frame(age="all",year= c(years,max(years+1)),unit="unique",season="all",area="unique",iter=1,data= c(apply(posteriors$HtoHmsy,2,quantile,0.5),NA),qname="harvest")
  )  
  flqs$qname = factor(flqs$qname,levels=c("biomass","catch","stock","harvest"))
  
  jabba$flqs =  flqs
  
  # add b.ppdist
  if(jbinput$jagsdata$b.pr[3]==0){jabba$bppd = "No biomass prior used"} else {
    if(jbinput$jagsdata$b.pr[4]==0){jabba$bppd = posteriors$P[,which(years%in%jbinput$jagsdata$b.pr[3])]} 
    if(jbinput$jagsdata$b.pr[4]==1){jabba$bppd = posteriors$BtoBmsy[,which(years%in%jbinput$jagsdata$b.pr[3])]}
    if(jbinput$jagsdata$b.pr[4]==2){jabba$bppd = posteriors$HtoHmsy[,which(years%in%jbinput$jagsdata$b.pr[3])]}
  }  
  #jabba$prj.inits = prj.inits # New for external projections 
  
  
  
  jabba$kbtrj = NULL
  if(save.trj==TRUE){
    
    ## Create new full KOBE output
    
    kb = NULL
    for(i in 1:N){
      if(settings$add.catch.CV==FALSE){Ci = rep(jabba$catch[i],length(posteriors$K))} else {Ci = posteriors$estC[,i]}    
      kb = rbind(kb,data.frame(year=years[i],run=jabba$scenario,type="fit",iter=1:length(posteriors$K),
                               stock=posteriors$BtoBmsy[,i],harvest=posteriors$HtoHmsy[,i],B=posteriors$SB[,i],H=posteriors$H[,i],
                               Bdev=posteriors$Proc.Dev[,i],Catch=Ci,BB0=posteriors$P[,i]))
    }
    jabba$kbtrj = kb
  }
  
  #--------------------------------------------------------
  # Projections
  #-------------------------------------------------------
  
  if(jbinput$settings$projection==TRUE){
    if(verbose)
      message("\n","><> compiling Future Projections under fixed quota <><","\n")
    pyrs = jbinput$jagsdata$pyrs
    TACs = jbinput$jagsdata$TAC[pyrs,] 
    nTAC = length(TACs) 
    proj.yrs =  (years[n.years]+1):(years[n.years]+pyrs)
    posteriors$prBtoBmsy
    kbprj = NULL
    for(j in 1:nTAC){
      kb.temp = kb[kb$year==max(years),]
      kb.temp$run = paste0("C",TACs[j])
      for(i in 1:pyrs){
        kb.temp = rbind(kb.temp,data.frame(year=proj.yrs[i],run=paste0("C",TACs[j]),type="prj",iter=1:length(posteriors$K),
                                           stock=posteriors$prBtoBmsy[,i,j],harvest=posteriors$prHtoHmsy[,i,j],B=posteriors$prBtoBmsy[,i,j]*posteriors$SBmsy,H=posteriors$prHtoHmsy[,i,j]*posteriors$Hmsy,
                                           Bdev=0,Catch= TACs[j],BB0=posteriors$prP[,i,j]))
      }
      
      kbprj = rbind(kbprj,kb.temp) 
    }
    kbprj$run = factor(kbprj$run,levels=unique(kbprj$run))
    
    jabba$kbtrj = rbind( jabba$kbtrj,kbprj)
  } # End of projections compilation
  
  # Safe posteriors (Produces large object!)
  if(save.all==TRUE) jabba$posteriors = posteriors
  
  if(save.jabba==TRUE){
    save(jabba,file=paste0(output.dir,"/",settings$assessment,"_",settings$scenario,"_jabba.rdata"))
  }
  
  
  if(save.csvs==TRUE){
    # Save results
    write.csv(Stock_trj[,2,],paste0(output.dir,"/Stock_trj_",settings$assessment,"_",settings$scenario,".csv"))
    # Save model estimates and convergence p-values
    write.csv(data.frame(results),paste0(output.dir,"/Estimates_",settings$assessment,"_",settings$scenario,".csv"))
    write.csv(Table,paste0(output.dir,"/Results_",settings$assessment,"_",settings$scenario,".csv"))
    write.csv(jabba$stats,paste0(output.dir,"/GoodnessFit_",settings$assessment,"_",settings$scenario,".csv"))
  }
  
  
  
  return(jabba)
  
} # end of fit_jabba()
