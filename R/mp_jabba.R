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
#' @param quickmcmc option to run short mcmc
#' @param verbose option show cat comments and progress
#' @param par.quantile quantile(s) of parameter posterior, default median 0.5 
#' @param b.quantile quantile(s) of biomass posterior, default median 0.5
#' @return A result list containing estimates of model input, settings and results
#' @export
#' @examples
#' data(iccat)
#' jbinput <- build_jabba(catch=iccat$bet$catch,cpue=iccat$bet$cpue,se=iccat$bet$se,model.type="Fox",verbose=FALSE)
#' system.time(fit_jabba(jbinput,quickmcmc=TRUE,verbose=FALSE))
#' system.time(mp_jabba(jbinput))

mp_jabba = function(jbinput,
                     # MCMC settings
                     ni = 30000, # Number of iterations
                     nt = 5, # Steps saved
                     nb = 5000, # Burn-in
                     nc = 2, # number of chains
                     # init values
                     init.values = FALSE,
                     init.K = NULL,
                     init.r = NULL,
                     init.q = NULL,# vector,
                     par.quantile = 0.5,
                     b.quantile = 0.5,
                     quickmcmc = TRUE,
                     verbose=FALSE
){
  
  tmpath <- tempfile()
  dir.create(tmpath)
  if(!verbose) {
    progress.bar="none"
  } else {
    progress.bar="text"
  }
  #write jabba model
  jabba2jags(jbinput, tmpath)
  
  # mcmc saved
  nsaved = (ni-nb)/nt*nc
  # jabba model data
  jbd = jbinput$jagsdata
  
  if(quickmcmc==TRUE){
    ni = 6000
    nb = 1000
    nt = 2
    nc = 2
  }
  
  # a, b pars for beta prior
  ab = get_beta(max(min(jbinput$settings$psi.prior.raw[1],0.95),0.05),CV=0.05/max(min(jbinput$settings$psi.prior.raw[1],0.95),0.05))
  
  
  # Initial starting values (new Eq)
  if(init.values==FALSE){
    inits = function(){list(K= rlnorm(1,log(jbd$K.pr[1])-0.5*0.3^2,0.3),r = rlnorm(1,log(jbd$r.pr[1]),jbd$r.pr[2]) ,q = runif(jbd$nq,min(jbd$I,na.rm=T)/max(jbd$TC,na.rm=T),mean(jbd$I,na.rm=T)/max(jbd$TC,na.rm=T)),psi=rbeta(1,ab[1],ab[2]),isigma2.est=runif(1,20,100), itau2=runif(jbd$nvar,80,200))}
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
 
   years = jbinput$data$yr
   peels = 0

  jbinput$jagsdata$I = jbd$I # update
  # jabba model building conditions
  params = c("K","SBmsy","Hmsy","MSY","SB")
  
  
  ptm <- proc.time()
  
  mod <- R2jags::jags(jbd, inits,params,file.path(tmpath,"JABBA.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, quiet=!verbose, progress.bar = progress.bar)  # adapt is burn-in
  
  proc.time() - ptm
  save.time = proc.time() - ptm
  # if run with library(rjags)
  posteriors = mod$BUGSoutput$sims.list
  
  # run some mcmc convergence tests
  par.dat= data.frame(posteriors[params[c(1:4)]])
  names(par.dat) = c("K","Bmsy","Fmsy","MSY")
  geweke = coda::geweke.diag(data.frame(par.dat))
  pvalues <- 2*pnorm(-abs(geweke$z))
  heidle = coda::heidel.diag(data.frame(par.dat))
  
  # Refpoints
  jabba = list(
  B = data.frame(age="all",year= c(years,max(years+1)),unit="unique",season="all",area="unique",iter=1,data= apply(posteriors$SB,2,quantile,b.quantile)),
  refpts = apply(par.dat,2,quantile,par.quantile),
  convergence = data.frame(Geweke.p=round(pvalues,3),Heidel.p = round(heidle[,3],3))
  )
  
  return(jabba)
  
} # end of fit_jabba()
