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
#' @param save.jabba = FALSE
#' @param save.all = FALSE
#' @param output.dir path to save plot. default is getwd()
#' @return A result list containing estimates of model input, settings and results
#' @export
#' @examples
#' data(bet)
#' jbinput <- build_jabba()

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
                     save.all = FALSE,
                     save.jabba = FALSE,
                     save.csvs = FALSE,
                     save.prjkobe = FALSE,
                     output.dir = getwd()
){
  #write jabba model
  jabba2jags(jbinput)
  
  # mcmc saved
  nsaved = (ni-nb)/nt*nc
  # jabba model data
  jbd = jbinput$jagsdata
  # Initial starting values (new Eq)
  if(init.values==FALSE){
    inits = function(){list(K= rlnorm(1,log(jbd$K.pr[1])-0.5*0.3^2,0.3),r = rlnorm(1,log(jbd$r.pr[1]),jbd$r.pr[2]) ,q = runif(jbd$nq,min(jbd$I,na.rm=T)/max(jbd$TC,na.rm=T),mean(jbd$I,na.rm=T)/max(jbd$TC,na.rm=T)), isigma2.est=runif(1,20,100), itau2=runif(jbd$nvar,80,200))}
  }else {
    if(is.null(init.K))
      stop("\n","\n","><> Provide init.K guess for option init.values=TRUE  <><","\n","\n")
    if(is.null(init.r))
      stop("\n","\n","><> Provide init.r guess for option init.values=TRUE  <><","\n","\n")
    if(is.null(init.q))
      stop("\n","\n","><> Provide init.q vector guess for option init.values=TRUE  <><","\n","\n")
    if(length(init.q)!= jbinput$jagsdata$nq)
      stop("\n","\n","><> init.q vector must match length of estimable q's, length(unique(sets.q))   <><","\n","\n")
    inits = function(){ list(K= init.K,r=init.r,q = init.q, isigma2.est=runif(1,20,100), itau2=runif(jbd$nvar,80,200))}
  }

  out = output.dir
  if(file.exists(out)==FALSE) stop("\n","\n","><> output.dir does not exist <><","\n","\n")


  # retrospecitive peel
  years = jbinput$data$yr
  if(is.null(peels)) peels = 0
  if(peels > 0){
    jbd$I[(length(years)-peels+1) : length(years),]  = NA
  
  }
  jbinput$jagsdata$I = jbd$I # update
  # jabba model building conditions
  params = jbinput$settings$params


  ptm <- proc.time()

  mod <- R2jags::jags(jbd, inits,params,paste0(tempdir(),"/JABBA.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)  # adapt is burn-in

  proc.time() - ptm
  save.time = proc.time() - ptm

  # unpack
  settings= c(jbinput$data,jbinput$jagsdata,jbinput$settings)
  catch = settings$catch
  cpue= settings$cpue
  se = settings$se
  n.years = settings$N
  years = settings$yr
  assessment = settings$assessment
  scenario = settings$scenario
  CPUE = settings$I
  n.indices = settings$nI

  # if run with library(rjags)
  posteriors = mod$BUGSoutput$sims.list
  cat(paste0("\n","><> Produce results output of ",settings$model.type," model for ",settings$assessment," ",settings$scenario," <><","\n"))

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
  #Model  parameter
  apply(par.dat,2,quantile,c(0.025,0.5,0.975))

  man.dat = data.frame(posteriors[params[8:10]])
  #Management quantaties
  apply(man.dat,2,quantile,c(0.025,0.5,0.975))

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
  colnames(Table) <- c("mu","lci","uci")

  #-----------------------------------------------
  # Stock trajectories
  #-----------------------------------------------
  #Bt = posteriors$SB
  #Ht = posteriors$H
  #Bt_Bmsy = posteriors$BtoBmsy
  #Ht_Hmsy = posteriors$HtoHmsy
  #Bt_K = posteriors$P

  Stock_trj = array(data=NA,dim=c(ncol(posteriors$SB),3,6),dimnames = list(years,c("mu","lci","uci"),c("B","F","BBmsy","FFmsy","BB0","procB")))
  for(i in 1:3){
    Stock_trj[,i,] =  cbind(t(apply(posteriors$SB,2,quantile,c(0.5,0.025,0.975)))[,i],
                            t(apply(posteriors$H,2,quantile,c(0.5,0.025,0.975)))[,i],
                            t(apply(posteriors$BtoBmsy,2,quantile,c(0.5,0.025,0.975)))[,i],
                            t(apply(posteriors$HtoHmsy,2,quantile,c(0.5,0.025,0.975)))[,i],t(apply(posteriors$P,2,quantile,c(0.5,0.025,0.975)))[,i],t(apply(posteriors$Proc.Dev,2,quantile,c(0.5,0.025,0.975)))[,i])

  }

  #--------------------------------------------------------
  # Projections
  #-------------------------------------------------------
  
  if(jbinput$settings$projection==TRUE){
    cat("\n","><> compiling Future Projections under fixed quota <><","\n")
    pyrs = jbinput$jagsdata$pyrs
    TACs = jbinput$jagsdata$TAC[1,] 
    nTAC = length(TACs) 
    proj.yrs =  years[n.years]:(years[n.years]+pyrs)
    # Dims 1: saved MCMC,2: Years, 3:alternatic TACs, 4: P, H/Hmsy, B/Bmsy
    projections = array(NA,c(nsaved,length(proj.yrs),nTAC,3),dimnames = list(1:nsaved,proj.yrs,TACs,c("BB0","BBmsy","FFmsy")))
    
    for(i in 1:nTAC){
      projections[,,i,1] = cbind(posteriors$P[,(n.years):n.years],posteriors$prP[,,i])
    }
    
    for(i in 1:nTAC){
      projections[,,i,2] = cbind(posteriors$BtoBmsy[,(n.years):n.years],posteriors$prBtoBmsy[,,i])
    }
    for(i in 1:nTAC){
      projections[,,i,3] = cbind(posteriors$HtoHmsy[,(n.years):n.years],posteriors$prHtoHmsy[,,i])
    }
    
    
    Stock_prj = array(data=NA,dim=c(pyrs+1,3,length(TACs),3),dimnames = list(proj.yrs,c("mu","lci","uci"),TACs,c("BB0","BBmsy","FFmsy")))
    for(j in 1:length(TACs)){
      for(i in 1:3){
        Stock_prj[,i,j,] =  cbind(t(apply(projections[,,j,"BB0"],2,quantile,c(0.5,0.1,0.9)))[,i],
                                  t(apply(projections[,,j,"BBmsy"],2,quantile,c(0.5,0.1,0.95)))[,i],
                                  t(apply(projections[,,j,"FFmsy"],2,quantile,c(0.5,0.1,0.95)))[,i])
        
      }}
    
    
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

  if(settings$CatchOnly==FALSE){
    # get residuals
    Resids = NULL
    for(i in 1:n.indices){
      Resids =rbind(Resids,log(CPUE[,i])-log(apply(posteriors$CPUE[,,i],2,quantile,c(0.5))))
    }

    # Standardized Residuals
    StResid = NULL
    for(i in 1:n.indices){
      StResid =rbind(StResid,log(CPUE[,i]/apply(posteriors$CPUE[,,i],2,quantile,c(0.5)))/
                       apply(posteriors$TOE[,,i],2,quantile,c(0.5))+0.5*apply(posteriors$TOE[,,i],2,quantile,c(0.5)))
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
  jabba.res = NULL
  if(settings$CatchOnly==FALSE){
    for(i in 1:n.indices){

      Yr = years
      Yr = min(Yr):max(Yr)
      yr = Yr-min(years)+1

      exp.i = apply(posteriors$CPUE[,,i],2,quantile,c(0.5))[is.na(cpue[,i+1])==F]
      expLCI.i = apply(posteriors$CPUE[,,i],2,quantile,c(0.025))[is.na(cpue[,i+1])==F]
      expUCI.i = apply(posteriors$CPUE[,,i],2,quantile,c(0.975))[is.na(cpue[,i+1])==F]

      obs.i = cpue[is.na(cpue[,i+1])==F,i+1]
      sigma.obs.i = (apply(posteriors$TOE[,,i],2,quantile,c(0.5)))[is.na(cpue[,i+1])==F]

      yr.i = Yr[is.na(cpue[,i+1])==F]
      jabba.res = rbind(jabba.res,data.frame(scenario=settings$scenario,name=names(cpue)[i+1],year=yr.i,obs=obs.i,obs.err=sigma.obs.i,hat=exp.i,hat.lci=expLCI.i,hat.uci=expUCI.i,residual=log(obs.i)-log(exp.i),retro.peels=peels))
    }
  }

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



  #-----------------------------------
  # Note posteriors of key parameters
  #-----------------------------------
  sel.par = c(1,2,7,4,3,5)
  out=data.frame(posteriors[params[sel.par]])

  cat(paste0("\n","\n",paste0("><> Scenario ", jbinput$settings$scenario,"_",jbinput$settings$model.type," completed in ",as.integer(save.time[3]/60)," min and ",round((save.time[3]/60-as.integer(save.time[3]/60))*100)," sec <><","\n")))



  #-------------------------------
  # summarize results in jabba list
  #-------------------------------
  jabba = list()
  jabba$assessment  = assessment
  jabba$scenario = scenario
  jabba$settings = c(jbinput$jagsdata,jbinput$settings)
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
  jabba$timeseries=Stock_trj
  jabba$refpts = data.frame(factor=assessment,level=settings$scenario,quant = c("hat","logse"), k=c(median(posteriors$K),sd(log(posteriors$K))),bmsy=c(median(posteriors$SBmsy),sd(log(posteriors$SBmsy))),
                            fmsy=c(median(posteriors$Hmsy),sd(log(posteriors$Hmsy))),msy=c(median(posteriors$SBmsy*posteriors$Hmsy),sd(log(posteriors$SBmsy*posteriors$Hmsy))))
  jabba$pfunc = data.frame(factor=assessment,level=scenario,SB_i=round(Bit,3),SP=round(SP,3),Hmsy=round(Hmsy.sp,4),r=round(Hmsy.sp*(m.sp-1)/(1-1/m.sp),4),m=round(m.sp,3),MSY=round(as.numeric(MSY.sp[2]),3),SB0=round(SB0.sp,3),Cmsy=round(Cmsy,3))

  if(settings$CatchOnly==FALSE){
    jabba$diags = data.frame(factor=assessment,level=settings$scenario,name=jabba.res$name,year=jabba.res$year,season=1,obs=jabba.res$obs,hat=jabba.res$hat,hat.lci=jabba.res$hat.lci,hat.uci=jabba.res$hat.uci,residual=jabba.res$residual,retro.peels=jabba.res$retro.peels)
    jabba$residuals = array(Resids,dim=c(nrow(Resids),ncol(Resids)),dimnames = list(names(cpue)[-1],years))
    jabba$std.residuals = array(StResid,dim=c(nrow(StResid),ncol(StResid)),dimnames = list(names(cpue)[-1],years))


  } else {
    jabba$diags = "Not Available for Catch-Only option"
    jabba$residuals = "Not Available for Catch-Only option"
  }
  jabba$stats = data.frame(Stastistic = c("N","p","DF","SDNR","RMSE","DIC"),Value = c(Nobs,npar,DF,SDNR,RMSE,DIC))
  jabba$pars_posterior = out
  jabba$kobe = data.frame(factor=assessment,level=scenario,yr=years[N],stock=posteriors$BtoBmsy[,N],harvest=posteriors$HtoHmsy[,N])
  # add b.ppdist
  if(jbinput$jagsdata$b.pr[3]==0){jabba$bppd = "No biomass prior used"} else {
  if(jbinput$jagsdata$b.pr[4]==0){jabba$bppd = posteriors$P[,which(years%in%jbinput$jagsdata$b.pr[3])]} 
  if(jbinput$jagsdata$b.pr[4]==1){jabba$bppd = posteriors$BtoBmsy[,which(years%in%jbinput$jagsdata$b.pr[3])]}
  if(jbinput$jagsdata$b.pr[4]==2){jabba$bppd = posteriors$HtoHmsy[,which(years%in%jbinput$jagsdata$b.pr[3])]}
  }  
  if(jbinput$settings$projection==FALSE){ jabba$projections = "Required setting projection = TRUE in build_jabba()" } else{
    jabba$projections = Stock_prj
  }
  
  if(save.jabba==TRUE){
  save(jabba,file=paste0(output.dir,"/",settings$assessment,"_",settings$scenario,"_jabba.rdata"))
  }
  
  if(save.prjkobe==TRUE){
    prjkb = kobeJabbaProj(projections)
    save(prjkb,file=paste0(output.dir,"/",settings$assessment,"_",settings$scenario,"_prjkb.rdata"))
  }
  
  
  # Safe posteriors (Produces large object!)
  if(save.all==TRUE) save(posteriors,file=paste0(output.dir,"/",settings$assessment,"_",settings$scenario,"_posteriors.rdata"))
  
  
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
