#' Creates a data list used as input to JAGS
#'
#' Creates a data list with JABBA input and settings to be passed to fit_jabba()
#' @param catch  catch time series, requires data.frame(year, catch)
#' @param cpue cpue time series,  requires data.frame(year, cpue.1,cpue.2,...,cpue.N)
#' @param se optional log standard error (CV) time series,requires data.frame(year, se.1,se.2,...,se.N)
#' @param assessment = "example",
#' @param scenario = "s1",
#' @param model.type = c("Schaefer","Fox","Pella","Pella_m"),
#' @param add.catch.CV = c(TRUE,FALSE) option estimate catch with error
#' @param catch.cv  catch error on log-scale (default = 0.1)
#' @param catch.error can be random or directional under reporting "under"
#' @param auxiliary additional time series of either c(Z,Effort,B/B0,B/Bmsy,F/Fmsy)
#' @param auxiliary.se optional input standard errors on auxilary data 
#' @param auxiliary.type c("effort","z","bk","bbmsy","ffmsy") 
#' @param auxiliary.lag  lag option in years (default 1) for effort, z and ffmsy  
#' @param Plim = 0, # Set Plim = Blim/K where recruitment may become impaired (e.g. Plim = 0.25)
#' PRIORS
#' @param r.dist = c("lnorm","range"), # prior distribution for the intrinsic rate population increas 
#' @param r.prior = c(0.2,0.5), # prior(mu, lod.sd) for intrinsic rate of population increase   
#' @param K.dist = c("lnorm","range"), # prior distribution for unfished biomass  K = B0 
#' @param K.prior = NULL, # prior(mu,CV) for the unfished biomass K = B0
#' @param psi.dist = c("lnorm","beta"), # prior distribution for the initial biomass depletion B[1]/K
#' @param psi.prior = c(0.9,0.25), # prior(mu, CV) for the initial biomass depletion B[1]/K
#' @param b.prior = c(FALSE,0.3,NA,c("bk","bbmy","ffmsy")[1]), # depletion prior set as d.prior = c(mean,cv,yr,type=c("bk","bbmsy"))
#' @param BmsyK = 0.4, # Inflection point of the surplus production curve, requires Pella-Tomlinson (model = 3 | model 4)
#' @param shape.CV = 0.3, # CV of the shape m parameters, if estimated with Pella-Tomlinson (Model 4)
#' VARIANCE options
#' @param igamma = c(3,0.01), # prior for process error variance, default informative igamma ~ mean 0.07, CV 0.4
#' @param sets.q = 1:(ncol(cpue)-1), # assigns catchability q to different CPUE indices. Default is each index a seperate q
#' @param sigma.est = TRUE, # Estimate additional observation variance
#' @param sets.var = 1:(ncol(cpue)-1), # estimate individual additional variace
#' @param fixed.obsE  # Minimum fixed observation error
#' @param auxiliary.obsE # Fixed observation error for auxiliary data   
#' @param auxiliary.sigma # TRUE/FALSE 
#' @param qA.cv precision on lognormal prior for e.g. qA*Z (not applicable to effort)
#' @param sets.varA estimate individual additional variance 
#' @param sigma.proc =  TRUE, # TRUE: Estimate process error, else set to value
#' @param proc.dev.all = TRUE, # TRUE: All year, year = starting year
#' @param projection = FALSE, # Switch on by Projection = TRUE
#' @param TACs = NULL
#' @param TACint =  NULL, # default avg last 3 years
#' @param imp.yr = NULL, # default last year plus ONE
#' @param pyrs = NULL, # Set number of projections years
#' @param P_bound = c(0.02,1.3),  # Soft penalty bounds for b/k
#' @param sigmaobs_bound = 1, # Adds an upper bound to the observation variance
#' @param sigmaproc_bound = 0.2, # Adds an upper bound to the process variance
#' @param q_bounds c(10^-30,1000), # Defines lower and upper bounds for q
#' @param K_bounds c(0.01,10^10), # Defines lower and upper bounds for K
#' @param qA_bounds Defines lower and upper bounds for q of auxiliary data type effort 
#' @param harvest.label = c("Hmsy","Fmsy")[2], # choose label preference H/Hmsy versus Fmsy
#' @param catch.metric  "(t)" # Define catch input metric e.g. (tons) "000 t" 
#' @param verbose option show cat comments
#' @return List to be used as data input to JABBA JAGS model.
#' @importFrom grDevices rainbow
#' @importFrom stats dlnorm
#' @export

build_jabba <- function(
  catch = NULL,
  cpue = NULL,
  se =  NULL,
  auxiliary = NULL,
  auxiliary.se = NULL, 
  auxiliary.type = c("z","effort","bk","bbmsy","ffmsy")[1],
  assessment = "jabba",
  scenario = "jabba",
  model.type = c("Schaefer","Fox","Pella","Pella_m"),
  add.catch.CV = TRUE, # to match original assessment
  catch.cv = 0.1, # CV for catch error
  catch.error = c("random","under")[1], # 
  Plim = 0, # Set Plim = Blim/K where recruitment may become impaired (e.g. Plim = 0.25)
  r.dist = c("lnorm","range"), # prior distribution for the intrinsic rate population increas 
  r.prior = c(0.2,0.5), # prior(mu, lod.sd) for intrinsic rate of population increase   
  K.dist = c("lnorm","range"), # prior distribution for unfished biomass  K = B0 
  K.prior = NULL, # prior(mu,CV) for the unfished biomass K = B0
  psi.dist = c("lnorm","beta"), # prior distribution for the initial biomass depletion B[1]/K
  psi.prior = c(0.9,0.25), # depletionprior(mu, CV) for the initial biomass depletion B[1]/K
  b.prior = c(FALSE,0.3,NA,c("bk","bbmsy","ffmsy")[1]), # depletion prior set as b.prior = c(mean,cv,yr,type=c("bk","bbmsy","ffmsy))
  BmsyK = 0.4, # Inflection point of the surplus production curve, requires Pella-Tomlinson (model = 3 | model 4)
  shape.CV = 0.3, # CV of the shape m parameters, if estimated with Pella-Tomlinson (Model 4)
  sets.q = 1:(ncol(cpue)-1), # assigns catchability q to different CPUE indices. Default is each index a seperate q
  sigma.est = TRUE, # Estimate additional observation variance
  sets.var = 1:(ncol(cpue)-1), # estimate individual additional variance
  fixed.obsE = ifelse(is.null(se),0.1,0.001), # Minimum fixed observation error
  auxiliary.obsE  = ifelse(is.null(auxiliary.se),0.1,0.001),
  auxiliary.sigma = TRUE,
  auxiliary.lag = 1,  # lags between measure of effort, Z or ffmsy in years
  qA.cv = 0.1,  # adjusts precision qA for Z, F, ffmsy or bbmsy  
  sets.varA = 1:(ncol(auxiliary)-1), # estimate individual additional variance
  sigma.proc =  TRUE, # TRUE: Estimate observation error, else set to value
  proc.dev.all = TRUE, # TRUE: All year, year = starting year
  igamma = c(4,0.01), # informative mean 0.05, CV 0.4 from original paper
  projection = FALSE, # Switch on by Projection = TRUE
  TACs = NULL, # vector of fixed catches used for projections  
  TACint =  NULL, # default avg last 3 years
  imp.yr = NULL, # default last year plus ONE
  pyrs = NULL, # Set number of projections years
  # Penalties
  P_bound = c(0.02,1.3),  # Soft penalty bounds for b/k
  sigmaobs_bound = 1, # Adds an upper bound to the observation variance
  sigmaproc_bound = 0.2, # Adds an upper bound to the process variance
  q_bounds= c(10^-30,1000), # Defines lower and upper bounds for q
  K_bounds= c(0.01,10^10), # Defines lower and upper bounds for K
  qA_bounds = c(10^-30,1000), # Defines lower and upper bounds for q of auxiliary data type "effort"
  # Settings
  harvest.label = c("Hmsy","Fmsy")[2], # choose label preference H/Hmsy versus Fmsy
  catch.metric = "(t)", # Define catch input metric e.g. (tons) "000 t" etc
  verbose=TRUE
){
  
  
  # define model typue
  mod.names = model.type[1]
  model = which(c("Schaefer","Fox","Pella","Pella_m")%in%mod.names)
  # Set up output structure
  
  
  #-------------------------
  # Prepare input data
  #-------------------------
  if(verbose)
    message("\n","><> Prepare JABBA input data <><","\n","\n")
  
  
  if(is.null(catch)) stop("\n","\n","><> Catch Time series not provided <><","\n","\n")
  years = catch[,1]
  styr = min(years)
  endyr = max(years)
  n.years = length(years)
  if(length(which(styr:endyr%in%catch[,1]==FALSE))>0) stop("\n","\n","><> catch time series must have years in sequential order <><","\n","\n")
  
  if(is.null(se)) SE.I = FALSE
  if(is.null(cpue)){
    CatchOnly = TRUE
    if(verbose)
      message(paste0("\n","><> Running Catch-Only mode: CatchOnly = TRUE <><","\n","\n"))
    cpue= catch[,1:2]
    colnames(cpue) = c("Year","Dummy Index")
    cpue[,2] = NA
    cpue[1,2] = 1
    SE.I = FALSE
    #add.catch.CV =FALSE
    sigma.est = FALSE
    sets.q =1
    sets.var=1
    n.indices = 1
  } else {
    CatchOnly = FALSE
    
  }
  
 
  
  
  if(nrow(catch)!=nrow(cpue)) stop("\n","\n","><> cpue and catch differ in number of year <><","\n","\n")
  catches = names(catch)[2:ncol(catch)]
  n.catches = length(catches)
  styr.catch = min(catch[,1])
  styr.C = styr.catch-styr+1
  conv.catch = as.numeric(rbind(matrix(rep(NA,(styr.C-1)*n.catches),styr.C-1,n.catches),as.matrix(catch[,-1])))
  
  if(length(which(conv.catch<0.0001)>0)) {
    if(verbose)
      message(paste0("\n","><> Warnning: Replacing 0 Catch by small constant 0.0001 <><","\n","\n"))
    conv.catch[conv.catch< 0.001]=0
  }
  if(length(which(is.na(conv.catch)))>0) stop("\n","\n","><> Missing Catch values currently not permitted (NAs should be 0) <><","\n","\n")
  
  # Build Catch input
  Catch=matrix(conv.catch,nrow=n.years,ncol=n.catches)
  if(ncol(Catch)>1){
    if(verbose)
      message(paste0("\n","><> Aggrigating multiple catch colums (assumed by fleet) to a single total catch column <><","\n","\n"))
    Catch = apply(Catch,1,sum)
  }
  TC = as.numeric(Catch) # Total Catch
  
  # Catch CV option.
  if(add.catch.CV==TRUE){
    if(length(catch.cv)>1) {CV.C =catch.cv[,2]} else {CV.C = rep(catch.cv,length(TC))}
    if(verbose)
      message("\n","><> Assume Catch with error CV = ",mean(CV.C)," <><","\n","\n")
  } else {
    if(verbose)
      message("\n","><> Assume Catch to be known without error <><","\n","\n")
  }
  
  
  
  # Build CPUE input
  indices = names(cpue)[2:ncol(cpue)]
  n.indices = max(length(indices),1)
  styr.cpue = min(cpue[,1])
  styr.I = styr.cpue-styr+1
  # Convert input data to matrices for JAGS input
  conv.cpue = as.numeric(rbind(matrix(rep(NA,(styr.I-1)*n.indices),styr.I-1,n.indices),as.matrix(cpue[,-1])))
  CPUE=matrix(conv.cpue,nrow=n.years,ncol=n.indices)
  
  # Build Standard Error matrix
  if(is.null(se)){SE.I=FALSE} else {SE.I=TRUE}
  if(SE.I==FALSE){
    if(verbose)
      message("\n","><> SE.I=FALSE: Creating cpue.se dummy matrix <><","\n","\n")
    se = cpue
    conv.se = as.numeric(rbind(matrix(rep(NA,(styr.I-1)*n.indices),styr.I-1,n.indices),as.matrix(cpue[,-1])))
    se2 = matrix(ifelse(fixed.obsE>0,fixed.obsE^2,10^-10),n.years,n.indices)#/2
  } else{
    if(ncol(cpue)!=ncol(se) | nrow(cpue)!=nrow(se)) stop("\n","\n","><> SE and CPUE matrix must match <><","\n","\n")
    conv.se = as.numeric(rbind(matrix(rep(NA,(styr.I-1)*n.indices),styr.I-1,n.indices),as.matrix(se[,-1])))
    se2 = matrix(ifelse(is.na(conv.se),0.3^2,conv.se)^2,n.years,n.indices)+fixed.obsE^2#/2
  }
  
  # Adding auxiliary data
  
  if(is.null(auxiliary)){
    Auxiliary = FALSE
  } else {
    Auxiliary =TRUE
    if(nrow(catch)!=nrow(auxiliary)) stop("\n","\n","><> auxiliary data and catch differ in number of year <><","\n","\n")
    styr.aux = min(auxiliary[,1])
    styr.A = styr.aux-styr+1
    nA= ncol(auxiliary)-1 
    SE.A = TRUE
    if(is.null(auxiliary.se)) SE.A = FALSE
    if(SE.A==FALSE){
      if(verbose)
        message("\n","><> SE.A=FALSE: Creating auxiliary.se dummy matrix <><","\n","\n")
      A.se = auxiliary
      convA.se = as.numeric(rbind(matrix(rep(NA,(styr.A-1)*nA),styr.A-1,nA),as.matrix(auxiliary[,-1])))
      A.se2 = matrix(ifelse(auxiliary.obsE>0,auxiliary.obsE^2,10^-10),n.years,nA)#/2
    } else{
      if(ncol(auxiliary)!=ncol(auxiliary.se) | nrow(auxiliary)!=nrow(auxiliary.se)) stop("\n","\n","><> SE and CPUE matrix must match <><","\n","\n")
      convA.se = as.numeric(rbind(matrix(rep(NA,(styr.A-1)*nA),styr.A-1,nA),as.matrix(auxiliary.se[,-1])))
      A.se2 = matrix(ifelse(is.na(convA.se),0.3^2,convA.se)^2,n.years,nA)+auxiliary.obsE^2#/2
    }
    
    conv.aux = as.numeric(rbind(matrix(rep(NA,(styr.A-1)*nA),styr.A-1,nA),as.matrix(auxiliary[,-1])))
    AUXI=matrix(conv.aux,nrow=n.years,ncol=nA)
    if(verbose)
      message("\n","><> Fit Auxiliary data of type: ", auxiliary.type ,"<><","\n","\n")
  } 
  
  
  # Add depletion prior
  if((b.prior[1])==FALSE){
    b.yr = rep(0,n.years) # activate by setting one year to one
    b.pr = c(0.5,0.1,0,0) # not active
  } else {
    b.yr = rep(0,n.years)
    if(as.numeric(b.prior[3]) %in% years){
      b.pr = as.numeric(c(b.prior[1:2],years[which(years %in% b.prior[3])],which(c("bk","bbmsy","ffmsy")%in%b.prior[4])-1))
      b.yr[which(years %in% b.prior[3])] =1
      
    } else {
      b.yr[n.years] = 1
      b.pr = as.numeric(c(b.prior[1:2],years[n.years],which(c("bk","bbmsy","ffmsy")%in%b.prior[4])-1))
    }}
  
  
  #---------------------
  # Index color palette
  #---------------------
  jabba.colors = as.character(c('#e6194b', "#3cb44b", "#ffe119",
                                "#0082c8","#f58231", "#911eb4",
                                "#46f0f0", "#f032e6", "#d2f53c",
                                "#fabebe", "#008080","#e6beff", "#aa6e28",rainbow(10)))
  
  
  #----------------------------------------------------------
  # Get JABBA parameterization and suplus production function
  #----------------------------------------------------------
  
  # For Pella-Tomlinson
  getshape = function(model,BmsyK){
    if(model==3 | model==4){
      #-----------------------------------------------
      # find inflection point
      ishape = NULL
      # Find shape for  SBmsytoK
      ishape = seq(0.1,10,0.001)
      
      check.shape =((ishape)^(-1/(ishape-1))-BmsyK)^2
      #  Set shape (> 0, with 1.001 ~ Fox and 2 = Schaefer)
      shape =  ishape[check.shape==min(check.shape)]
    } else {shape=FALSE}
    return(shape)
  }
  shape = getshape(model=model,BmsyK=BmsyK)
  
  # Set shape m for Fox and Schaefer: Fox m ~1; Schaefer m =2
  if(shape==FALSE){
    if(model == 1){m=2} else {m = 1.001}}else{m=shape}
  
  #------------------------------------------------
  
  
  #----------------------------------------------------
  #----------------------------------------------------
  if(r.dist[1]=="range"){
    # initial range of r based on resilience (FishBase.org)
    if(length(r.prior)>1){ start.r = r.prior} else
      if(r.prior == "High") {
        start.r <- c(0.6,1.5)} else if(r.prior == "Medium") {
          start.r <- c(0.2,0.8)}    else if(r.prior == "Low") {
            start.r <- c(0.05,0.5)}  else { # i.e. res== "Very low"
              start.r <- c(0.015,0.1)}
    
    log.r = mean(log(start.r))
    sd.r = abs(log.r - log(start.r[1]))/2
    r.prior = c(exp(log.r),sd.r)
    CV.r = sqrt(exp(sd.r^2)-1)
  } else {
    log.r = log(r.prior[1])
    sd.r = r.prior[2]
    CV.r = sqrt(exp(sd.r^2)-1)
  }
  
  #----------------------------------------------------
  # Prepare K prior
  #----------------------------------------------------
  if(is.null(K.prior)) K.prior = c(8*max(Catch),1) # mean and CV
  if(length(K.prior)<2) K.prior[2] = 1
  
  
  if(K.dist[1]=="range"){
    log.K =  mean(log(K.prior))
    sd.K= abs(log.K - log(K.prior[1]))/2
    log.K = mean(log(K.prior))+0.5*sd.K^2 # Will be back-bias corrected in lognorm function
    CV.K = sqrt(exp(sd.K^2)-1)
  } else {
    
    CV.K = K.prior[2]
    sd.K=sqrt(log(CV.K^2+1))
    log.K = log(K.prior[1])#-0.5*sd.K^2
  }
  
  
  
  
  # Get input priors
  K.pr = plot_lnorm(exp(log.K),CV.K,Prior="K")
  r.pr = plot_lnorm(mu=exp(log.r),CV=CV.r,Prior="r")
  
  psi.dist = psi.dist[1]
  if(psi.dist=="beta"){
    psi.pr = get_beta(mu=psi.prior[1],CV=psi.prior[2],Min=0,Prior=paste0("Prior B(",years[1],")/K"))} else {
      psi.pr = plot_lnorm(mu=psi.prior[1],CV=psi.prior[2],Prior=paste0("Prior B(",years[1],")/K"))
    }
  if(model==4){
    dm = dlnorm((seq(0.001,5,0.1)),log(m),shape.CV)
    dm = dm/max(dm)
    bmsyk  = (seq(0.001,5,0.1))^(-1/(seq(0.001,5,0.1)-1))
  }
  
  
  # Note PRIORS and save input subfolder
  Priors =rbind(c(K.pr[1],CV.K),psi.prior,c(r.pr[1],CV.r))
  row.names(Priors) = c("K","Psi","r")
  colnames(Priors) = c("Mean","CV")
  Priors = data.frame(Priors)
  Priors$log.sd = sqrt(log(Priors[,2]^2+1))
  
  
  if(verbose) {
    message("\n","><> Model type:",model.type," <><","\n")
    if(model<4){message("\n","><> Shape m =",m ,"\n")} else {message("\n","><> Shape m is estmated with a mean",m,"and a CV",shape.CV,"\n")}
    message("\n","><> K prior mean =",K.pr[1],"and CV =", CV.K,"(log.sd = ",K.pr[2],")","\n")
    message("\n","><> r prior mean =",r.pr[1],"and CV =", CV.r,"(log.sd = ",r.pr[2],")","\n")
    message("\n","><> Psi (B1/K) prior mean =",psi.prior[1],"and CV =",psi.prior[2],"with",psi.dist,"destribution","\n")
  }
  
  
  
  #----------------------------------------------------------
  # Set up JABBA
  #----------------------------------------------------------
  if(verbose)
    message("\n","\n","\n","><> ALWAYS ENSURE to adjust default settings to your specific stock <><","\n","\n")
  # Plot MSY
  # remove scientific numbers
  options(scipen=999)
  
  #----------------------------------------------------------
  # Setup TAC projection
  #---------------------------------------------------------
  if(is.null(TACs)) TACs = round(seq(0.5,1.2,0.1)*TC[n.years],0)
  if(is.null(TACint)) TACint =  mean(c(TC[length(TC)-1],TC[length(TC)])) # avg last 3 years
  # default avg last 3 years
  if(is.null(imp.yr))imp.yr = years[n.years]+1 # default last year plus ONE
  if(is.null(pyrs)) pyrs = 10 # Set number of projections years
  
  if(projection==TRUE){
    nTAC = length(TACs)
    TAC = mat.or.vec(pyrs,nTAC)
    yr.last = max(years) # assessment year
    
    if(length(TACint)>1 & length(TACint) < imp.yr-yr.last-1){
      stop("too few intial catches for initial projections specified")}
    
    if(length(TACint)==1) TACint = rep(TACint,max(imp.yr-yr.last-1,1))
    
    for(i in 1:nTAC){
      if(imp.yr-yr.last-1>0){
        TAC[,i] = c(TACint,rep(TACs[i],pyrs-(imp.yr-yr.last-1)))} else{
          TAC[,i] = c(rep(TACs[i],pyrs-(imp.yr-yr.last-1)))
        }
    }
    
  }else{
    nTAC = 1
    TAC = TC[n.years] #
    pyrs = 1
  }
  
  
  
  
  #---------------------------------------------------------------
  # JABBA Schaefer/Fox Models 1-2, Pella 3, 4 estimate shape
  #---------------------------------------------------------------
  
  # Slope of hockey-stick
  slope.HS = ifelse(Plim==0,1/10^-10,1/Plim)
  nSel = 1 # setup for JABBA-SELECT version (in prep)
  nI = ncol(CPUE) # number of CPUE series
  stI = ifelse(proc.dev.all==TRUE,1, c(1:n.years)[is.na(apply(CPUE,1,mean,na.rm=TRUE))==FALSE][1]) #first year with CPUE
  
  # starting values
  nq = length(unique(sets.q))
  nvar = length(unique(sets.var))
  
  
  # JABBA input data
  surplus.dat = list(N=n.years, TC = TC,I=CPUE,SE2=se2,r.pr=r.pr,psi.pr=psi.pr,K.pr = K.pr,
                     nq=nq,nI = nI,nvar=nvar,sigma.fixed=ifelse(sigma.proc==TRUE,0,sigma.proc),
                     sets.var=sets.var, sets.q=sets.q,Plim=Plim,slope.HS=slope.HS,
                     nTAC=nTAC,pyrs=pyrs,TAC=TAC,igamma = igamma,stI=stI,pen.P = rep(0,n.years) ,pen.bk = rep(0,n.years),proc.pen=0,K.pen = 0,
                     obs.pen = rep(0,nvar),P_bound=P_bound,q_bounds=q_bounds,sigmaobs_bound=sigmaobs_bound,sigmaproc_bound=sigmaproc_bound,K_bounds=K_bounds,mu.m=m,b.yr=b.yr, b.pr = b.pr)
  
  
  
  # PARAMETERS TO MONITOR
  params <- c("K","r", "q", "psi","sigma2", "tau2","m","Hmsy","SBmsy", "MSY", "BtoBmsy","HtoHmsy","CPUE","Ihat","Proc.Dev","P","SB","H","prP","prBtoBmsy","prHtoHmsy","TOE")
  
  
  #-----------------------------------------------
  # If shape parameter is estimated (Model =4)
  if(model==4){
    surplus.dat$m.CV = shape.CV }
  #-----------------------------------------------
  # If Catch Estimation with CV is used
  if(add.catch.CV==TRUE){
    surplus.dat$CV.C = CV.C
    params = c(params,"estC")
    
  }
  
  if(Auxiliary){
    surplus.dat$A = AUXI
    surplus.dat$A.SE2 = A.se2
    surplus.dat$qA_bounds =qA_bounds 
    if(auxiliary.lag>0) surplus.dat$A.lag = auxiliary.lag
    surplus.dat$nA = nA
    params = c(params,c("Ahat","qA","AUXI","TAE"))
    if(auxiliary.sigma){
    surplus.dat$sets.varA = sets.varA 
    surplus.dat$nAvar = length(sets.varA) 
    params = c(params,c("eta2"))
    }
    if(auxiliary.type%in%c("z","f","ffmsy","bbmsy","bk")) surplus.dat$qA.cv = qA.cv 
  }
  
  
  #--------------------------
  # Capture Settings
  #--------------------------
  jbinput = list()
  jbinput$data = list()
  jbinput$jagsdata = list()
  jbinput$settings = list()
  jbinput$data$yr = years
  jbinput$data$catch = catch
  jbinput$data$cpue = cpue
  jbinput$data$se = se
  jbinput$data$auxiliary = auxiliary
  jbinput$data$auxiliary.se = auxiliary.se
  jbinput$jagsdata = surplus.dat
  jbinput$settings$params = params
  jbinput$settings$BmsyK = BmsyK 
  jbinput$settings$psi.dist = psi.dist
  jbinput$settings$psi.prior.raw = psi.prior
  jbinput$settings$SE.I = SE.I
  jbinput$settings$sigma.proc = sigma.proc
  jbinput$settings$sigma.est= sigma.est
  jbinput$settings$model.id = model
  jbinput$settings$model.type = mod.names
  jbinput$settings$add.catch.CV = add.catch.CV
  jbinput$settings$catch.cv = catch.cv
  jbinput$settings$catch.error = catch.error
  jbinput$settings$CatchOnly = CatchOnly
  jbinput$settings$proc.dev.all = proc.dev.all
  jbinput$settings$Auxiliary = Auxiliary
  jbinput$settings$auxiliary.type = auxiliary.type
  jbinput$settings$auxiliary.lag = auxiliary.lag
  jbinput$settings$auxiliary.sigma = auxiliary.sigma 
  jbinput$settings$projection = projection
  jbinput$settings$TAC.implementation = imp.yr
  jbinput$settings$catch.metric = catch.metric
  jbinput$settings$harvest.label = harvest.label
  jbinput$settings$assessment = assessment
  jbinput$settings$scenario = scenario
  jbinput$settings$cols = jabba.colors
  #capture.output(jbinput, file=paste0(output.dir,"/Settings.txt"))
  #-------------------------------------------------------------------
  # write JAGS MODEL
  
  #jabba2jags(jbinput)
  
  return(jbinput)
  
} # end of build_jabba()


