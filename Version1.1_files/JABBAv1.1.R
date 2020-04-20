##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><><
## Stock Assessment execution File for JABBA
## Developed by Henning Winker & Felipe Carvalho (Cape Town/Hawaii)  
##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><><

cat(paste0("\n","><>><>><>><>><>><>><>><>"))
cat(paste0("\n","><> Run Model ",Mod.names,"<><"))
cat(paste0("\n","><>><>><>><>><>><>><>><>","\n","\n"))
# setwd(paste(File))
dir.create(paste0(File,"/",assessment),showWarnings = FALSE)
dir.create(paste0(File,"/",assessment,"/",Scenario,"_",Mod.names),showWarnings = FALSE)
dir.create(paste0(File,"/",assessment,"/",Scenario,"_",Mod.names,"/Input"),showWarnings = FALSE)
input.dir = paste0(File,"/",assessment,"/",Scenario,"_",Mod.names,"/Input")

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Define objects to make sure they exist if not included in Prime file
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
if(exists("igamma")==FALSE) igamma = c(4,0.01)  # Generic process error prior
if(exists("BmsyK")==FALSE) BmsyK = 0.4  # JABBA default for Pella model
if(exists("Model")==FALSE){ model = 1; Mod.names = c("Schaefer")} # Run Schaefer if model is not specified
if(exists("proc.dev.all")==FALSE) proc.dev.all = FALSE # process error deviation if only catch is available  
if(exists("Plim")==FALSE) Plim = 0  # Standard non-compound model
if(exists("P_bound")==FALSE) P_bound = c(0.02,1)  # Soft penalty bounds for P 
if(exists("KOBE.plot")==FALSE) KOBE.plot = TRUE # Produces JABBA Kobe plot 
if(exists("KOBE.type")==FALSE) KOBE.type = c("ICCAT","IOTC")[2] # ICCAT uses 3 colors; IOTC 4 (incl. orange) 
if(exists("SP.plot")==FALSE) SP.plot = c("standard","phase")[2] # Produces standard or 'Kobe phase' SP plot  
if(exists("Biplot")==FALSE) Biplot= TRUE # Produces a "post-modern" biplot with buffer and target zones (Quinn & Collie 2005)
if(exists("save.trajectories")==FALSE) save.trajectories =FALSE # saves posteriors of P=B/K, B/Bmsy and H/Hmsy as .RData object 
if(exists("harvest.label")==FALSE) harvest.label = c("Hmsy","Fmsy")[2] # choose label preference H/Hmsy versus Fmsy
if(exists("CPUE.plot")==FALSE) CPUE.plot= TRUE # Runs state-tool to produce "alligned" multi-CPUE plot  
if(exists("catch.metric")==FALSE) catch.metric = "(t)" # Runs state-tool to produce "alligned" multi-CPUE plot  
if(exists("meanCPUE")==FALSE) meanCPUE = FALSE # Uses averaged CPUE from state-space tool instead of individual indices  
if(exists("Projection")==FALSE) Projection = FALSE # Use Projections: requires to define TACs vectors 
if(exists("save.projections")==FALSE) save.projections = FALSE# saves projection posteriors as .RData object 
if(exists("Reproduce.seed")==FALSE) Reproduce.seed = FALSE # If FALSE a random seed assigned to each run (default)
if(exists("TACint")==FALSE) TACint = mean(catch[nrow(catch)-3,2]:catch[nrow(catch),2]) # use mean catch from last years
if(exists("imp.yr")==FALSE) imp.yr = as.numeric(format(Sys.Date(), "%Y"))+1 # use next year from now
if(exists("st.value")==FALSE) st.value =FALSE

# Save entire posterior as .RData object
if(exists("save.all")==FALSE) save.all = FALSE #  
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>




#-------------------------
# Prepare input data
#-------------------------
cat(paste0("\n","><> Prepare input data <><","\n","\n"))
indices = names(cpue)[2:ncol(cpue)]
n.indices = max(length(indices),1)
catches = names(catch)[2:ncol(catch)]
n.catches = length(catches)

years=catch[,1]
styr = min(years)
endyr = max(years)
n.years = length(years)
styr.cpue = min(cpue[,1])
styr.I = styr.cpue-styr+1 


# Convert input data to matrices for JAGS input
conv.cpue = as.numeric(rbind(matrix(rep(NA,(styr.I-1)*n.indices),styr.I-1,n.indices),as.matrix(cpue[,-1])))
CPUE=matrix(conv.cpue,nrow=n.years,ncol=n.indices)


if(SE.I==FALSE){
  se = cpue  
  conv.se = as.numeric(rbind(matrix(rep(NA,(styr.I-1)*n.indices),styr.I-1,n.indices),as.matrix(cpue[,-1])))
  se2 = matrix(ifelse(fixed.obsE>0,fixed.obsE^2,10^-10),n.years,n.indices)#/2
} else{
  conv.se = as.numeric(rbind(matrix(rep(NA,(styr.I-1)*n.indices),styr.I-1,n.indices),as.matrix(se[,-1])))
  #conv.se = sqrt(conv.se^2+fixed.obsE^2) 
  se2 = matrix(ifelse(is.na(conv.se),0.3^2,conv.se)^2,n.years,n.indices)+fixed.obsE^2#/2
}

conv.catch = as.numeric(rbind(matrix(rep(NA,(styr.I-1)*n.catches),styr.I-1,n.catches),as.matrix(catch[,-1])))
Catch=matrix(conv.catch,nrow=n.years,ncol=n.catches)
Catch[is.na(Catch)] = 0 # Replace any NA by zero 

# Total Catch
TC = apply(Catch,1,sum)

#------------
# Plot Catch
#------------

cat(paste0("\n","><> Plot Catch in Input subfolder <><","\n","\n"))

Par = list(mfrow=c(1,1),mar = c(5, 5, 1, 1), mgp =c(3,1,0), tck = -0.02,cex=0.8)
png(file = paste0(input.dir,"/Catches_",assessment,".png"), width = 7, height = 5, 
    res = 200, units = "in")
par(Par)
plot(catch[,1],catch[,1],ylim=c(0,max(catch[,2:ncol(catch)],na.rm=TRUE)),ylab=paste0("Catch ",catch.metric),xlab="Year",type="n")
for(i in 2:ncol(catch)) lines(catch[,1],catch[,i],lty=(i-1),lwd=2)
legend("topright",paste(names(catch)[2:ncol(catch)]),lty=1:(ncol(catch)-1),bty="n")
dev.off()

#---------------------
# Index color palette
#---------------------
jabba.colors = as.character(rep(c('#e6194b', "#3cb44b", "#ffe119",
                                  "#0082c8","#f58231", "#911eb4",
                                  "#46f0f0", "#f032e6", "#d2f53c",
                                  "#fabebe", "#008080","#e6beff", "#aa6e28"),2))
#--------------------
# Set seed
#--------------------
if(Reproduce.seed==FALSE){ set.seed(ceiling(runif(1,min=0,max=1e6))) } else {set.seed(123)}

#---------------------------------------------------------------------------
# CPUE run State-Space model for averaging CPUE
#---------------------------------------------------------------------------
if(CPUE.plot==TRUE){ 
  cat(paste0("\n","><> Run State-Space CPUE averaging tool","\n"))
  #find first time-series with first CPUE
  q1.y = c(1:n.years)[is.na(apply(CPUE,1,mean,na.rm=TRUE))==FALSE][1] #first year with CPUE
  q1.I = which.max(CPUE[q1.y,])
  
  qs = c(q1.I,c(1:(ncol(cpue)-1))[-q1.I])
  
  
  sink("cpueAVG.jags")
  cat("
      model {
      
      # Prior specifications  
      eps <- 0.0000000000001 # small constant    
      
      iq[1] ~ dgamma(1000,1000)
      q[1] <-  pow(iq[1],-1)
      logq[1] <- log(1)
      for(i in 2:nI){
      iq[i] ~ dgamma(0.001,0.001)
      q[i] <- pow(iq[i],-1)
      logq[i] <-  log(q[i])
      }
      
      
      ")
  
  if(sigma.proc==TRUE){
    cat("
        # Process variance
        isigma2 <- isigma2.est 
        sigma2 <- pow(isigma2,-1)
        sigma <- sqrt(sigma2)
        fakesigma.fixed <- sigma.fixed # Prevent unused variable error msg    
        ",append=TRUE)  
  }else{ cat(" 
      isigma2 <- pow(sigma.fixed+eps,-2) 
             sigma2 <- pow(isigma2,-1)
             sigma <- sqrt(sigma2)
             
             ",append=TRUE)}
  
  if(sigma.est==TRUE){
    cat("
        # Obsevation variance
        # Observation error
        itau2~ dgamma(0.001,0.001)
        tau2 <- 1/itau2
        
        
        for(i in 1:nI)
        {
        for(t in 1:N)
        {
        var.obs[t,i] <- SE2[t,i]+tau2
        ivar.obs[t,i] <- 1/var.obs[t,i]
        # note total observation error (TOE)     
        TOE[t,i] <- sqrt(var.obs[t,i])
        
        }}
        ",append=TRUE)  
  }else{ cat(" 
      # Obsevation variance
             # Observation error
             itau2~ dgamma(2,2)
             tau2 <- 1/itau2
             
             
             for(i in 1:nI)
             {
             for(t in 1:N)
             {
             var.obs[t,i] <- SE2[t,i] # drop tau2
             fake.tau[t,i] <- tau2
             
             ivar.obs[t,i] <- 1/var.obs[t,i]
             # note total observation error (TOE)     
             TOE[t,i] <- sqrt(var.obs[t,i])
             
             }}
             
             ",append=TRUE)}
  
  # Run rest of code  
  cat("  
      # Process variance prior
      isigma2.est ~ dgamma(0.001,0.001)
      
      
      # Priors and constraints
      logY.est[1] ~ dnorm(logY1, 1)       # Prior for initial population size
      
      mean.r ~ dnorm(1, 0.001)             # Prior for mean growth rate
      
      # Likelihood
      # State process
      for (t in 1:(N-1)){
      r[t] ~ dnorm(mean.r, isigma2)
      logY.est[t+1] <- logY.est[t] + r[t] }
      
      # Observation process
      for (t in 1:N) {
      for(i in 1:nI){
      y[t,i] ~ dnorm(logY.est[t]+logq[i], ivar.obs[t,i])
      }}
      
      # Population sizes on real scale
      for (t in 1:N) {
      Y.est[t] <- exp(logY.est[t])
      }
      
  } 
      ",fill = TRUE)
  sink()
  
  
  
  q.init = 1
  mCPUE = as.matrix(CPUE[q1.y:n.years,qs])
  mSE2 = as.matrix(se2[q1.y:n.years,qs])
  if(n.indices>1) for(i in 2:n.indices){q.init[i] = mean(mCPUE[,i],na.rm=TRUE)/mean(mCPUE[,1],na.rm=TRUE)}
  # Bundle data
  jags.data <- list(y = log(mCPUE),SE2=mSE2, logY1 = log(mCPUE[1,1]), N = length(q1.y:n.years),nI=n.indices,sigma.fixed=ifelse(sigma.proc==TRUE,0,sigma.proc))
  
  # Initial values
  inits <- function(){list(isigma2.est=runif(1,20,100), itau2=runif(1,80,200), mean.r = rnorm(1),iq = 1/q.init)}
  
  # Parameters monitored
  parameters <- c("mean.r", "sigma","r", "Y.est","q")
  
  
  # Call JAGS from R (BRT 3 min)
  mod.cpue <- jags(jags.data, inits, parameters, "cpueAVG.jags", n.chains = nc, n.thin = max(nt,2), n.iter = max(ni/5,10000), n.burnin = nb/10)
  
  
  cat(paste0("\n","><> Plot State-Space CPUE fits  in Input subfolder <><","\n"))
  # get individual trends
  fitted <- lower <- upper <- NULL
  cpue.yrs = years[q1.y:n.years]
  
  for (t in 1:nrow(mCPUE)){
    fitted[t] <- median(mod.cpue$BUGSoutput$sims.list$Y.est[,t])
    lower[t] <- quantile(mod.cpue$BUGSoutput$sims.list$Y.est[,t], 0.025)
    upper[t] <- quantile(mod.cpue$BUGSoutput$sims.list$Y.est[,t], 0.975)}
  
  
  q.adj = apply(mod.cpue$BUGSoutput$sims.list$q,2,median)
  
  
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  png(file = paste0(input.dir,"/CPUE_",assessment,"_",Scenario,".png"), width = 5, height = 3.5, 
      res = 200, units = "in")
  par(Par)
  u.ylim = NULL
  for(i in 1:n.indices){ u.ylim = c(u.ylim,exp(log(mCPUE[,i]/q.adj[i])+1.96*sqrt(mSE2[,i])))}  
  ylim = c(0,max(u.ylim,na.rm=TRUE))
  plot(0, 0, ylim = ylim, xlim = range(cpue.yrs), ylab = "Expected CPUE", xlab = "Year", col = "black", type = "n")
  legend("topright",paste(indices),lwd=2,col=(jabba.colors)[1:n.indices],bty="n")
  polygon(x = c(cpue.yrs,rev(cpue.yrs)), y = c(lower,rev(upper)), col = "gray", border = "gray90")
  
  for(i in 1:n.indices)
  {
    shift = runif(1,-0.1,0.1)
    cols=jabba.colors[qs[i]]
    plotCI(cpue.yrs+shift,mCPUE[,i]/q.adj[i],ui=exp(log(mCPUE[,i]/q.adj[i])+1.96*sqrt(mSE2[,i])),li=exp(log(mCPUE[,i]/q.adj[i])-1.96*sqrt(mSE2[,i])),add=TRUE,col= cols,pt.bg = cols,pch=21,gap=0)
    lines(cpue.yrs+shift,mCPUE[,i]/q.adj[i], col = cols,lwd=2)
    points(cpue.yrs+shift,mCPUE[,i]/q.adj[i], bg = cols,pch=21)
  }
  lines(cpue.yrs,fitted,lwd=2)
  
  dev.off()
  
  logSE = apply(log(mod.cpue$BUGSoutput$sims.list$Y.est),2,sd)
  
  
  if(nrow(mCPUE)<n.years) {
    fitted = c(rep(NA,q1.y-1),fitted)
    logSE = c(rep(0.2,q1.y-1),logSE)
  }    
  avgCPUE = data.frame(Year=years,CPUE= fitted,logSE=logSE)
  
  write.csv(avgCPUE,paste0(input.dir,"/avgCPUE_",assessment,"_",Scenario,".csv"))
  
  if(meanCPUE==TRUE){
    cat(paste0("\n","><> Use average CPUE as input for JABBA <><","\n"))
    
    CPUE = as.matrix(avgCPUE[,2]) 
    cpue.check = cpue[,-1]
    cpue.check[is.na(cpue[,-1])]=0
    CPUE[,1] = ifelse(apply(cpue.check,1,sum)==0,rep(NA,length(CPUE[,1])),CPUE[,1])
    se2 =  as.matrix(avgCPUE[,3]^2)     
    n.indices=1
    indices = "All"
    sets.q =1
    sets.var =1
  }
  
  }


#--------------------------------------------------------------------------
# END of CPUE State-Space tool
#--------------------------------------------------------------------------


#-----------
# FUNCTIONS
#-----------
cat(paste0("\n","><> Prepare JABBA prior inputs <><","\n"))

#--------------------------------------------------
# Function to get beta prior parameters
#--------------------------------------------------
get_beta <- function(mu,CV,Min=0,Prior="x"){
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
  plot(x,pdf,type="l",xlim=range(x[pdf>0.01]),xlab=paste(Prior),ylab="",yaxt="n")
  polygon(c(x,rev(x)),c(rep(0,length(x)),rev(ifelse(pdf==Inf,100000,pdf))),col="grey")
  return(c(a,b))
}

#--------------------------------------------------
# Function to get gamma prior parameters
#--------------------------------------------------

get_gamma <- function(mu,CV,Prior="x"){
  a = seq(0.0001,1000,0.0001)
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
  plot(x,pdf,type="l",xlim=range(x[pdf>0.01]),xlab=paste(Prior),ylab="",yaxt="n")
  polygon(c(x,rev(x)),c(rep(0,length(x)),rev(ifelse(pdf==Inf,100000,pdf))),col="grey")
  return(c(a,b))
}


#--------------------------------------------------
# Function to get lognormal prior parameters
#--------------------------------------------------
plot_lnorm <- function(mu,CV,Prior="x"){
  sdev= sqrt(log(CV^2+1))
  rand.pr = rlnorm(1000,log(mu),sdev)
  x = seq(min(rand.pr),quantile(rand.pr,0.995),max(rand.pr/500))  
  pdf = dlnorm(x,log(mu),sdev)  
  plot(x,pdf,type="l",xlim=range(x),xlab=paste(Prior),ylab="",yaxt="n")
  polygon(c(x,rev(x)),c(rep(0,length(x)),rev(ifelse(pdf==Inf,100000,pdf))),col="grey")
  return(c(mu,sdev))
}

#------------------------------------
# Function kobeJabba for FLR
#------------------------------------
kobeJabba<-function(x,minyear=1){
  
  out=cbind(melt(x[,,2]),c(x[,,3]))
  names(out)=c("iter","year","stock","harvest")
  out$year=out$year+minyear-1
  out}

#-------------------------------------------------
# Function kobeJabbaProj for projections with FLR
#-------------------------------------------------
kobeJabbaProj<-function(x,minyear=1,tac=NULL){
  
  out=cbind(melt(x[,,,2]),c(x[,,,3]))
  names(out)=c("iter","year","tac","stock","harvest")
  out$year=out$year+minyear-1
  
  out}

#----------------------------------------------------
# Determine initial ranges for r
#----------------------------------------------------
if(r.dist=="range"){
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
if(K.dist=="range"){
  log.K = mean(log(K.prior))
  sd.K= abs(log.K - log(K.prior[1]))/2
  CV.K = sqrt(exp(sd.K^2)-1)
} else {
  
  log.K = log(K.prior[1])
  CV.K = K.prior[2]
  sd.K=sqrt(log(CV.K^2+1))
}


#----------------------------------------------------------
# Get JABBA parameterization and suplus production function
#----------------------------------------------------------

# For Pella-Tomlinson
if(Model==3){ 
  #-----------------------------------------------
  # find inflection point
  ishape = NULL
  # Find shape for  SBmsytoK 
  ishape = seq(0.1,10,0.001)
  
  check.shape =((ishape)^(-1/(ishape-1))-BmsyK)^2
  
  #  Set shape (> 0, with 1.001 ~ Fox and 2 = Schaefer)
  shape =  ishape[check.shape==min(check.shape)] 
} else {shape=FALSE}
#------------------------------------------------


# Set shape m for Fox and Schaefer: Fox m ~1; Schaefer m =2
if(shape==FALSE){
  if(Model == 1){m=2} else {m = 1.001}}else{m=shape}

cat(paste0("\n","><> Plot Prior distributions in Input subfolder  <><","\n"))

Par = list(mfrow=c(1,3),mai=c(0.5,0.1,0,.1),omi = c(0.1,0.2,0.1,0) + 0.1,mgp=c(2,1,0), tck = -0.02,cex=0.8)
png(file = paste0(input.dir,"/Priors_",assessment,"_",Scenario,".png"), width = 9, height = 3, 
    res = 200, units = "in")
par(Par)
K.pr = plot_lnorm(exp(log.K),CV.K,Prior="K")

if(psi.dist=="beta"){
  psi.pr = get_beta(mu=psi.prior[1],CV=psi.prior[2],Min=0,Prior=paste0("Prior B(",years[1],")/K"))} else {
    psi.pr = plot_lnorm(mu=psi.prior[1],CV=psi.prior[2],Prior=paste0("Prior B(",years[1],")/K"))  
  }


r.pr = plot_lnorm(mu=exp(log.r),CV=CV.r,Prior="r")
mtext(paste("Density"), side=2, outer=TRUE, at=0.5,line=1,cex=0.9)
dev.off() 


cat(paste0("\n","><> Plot assumed Surplus Production shape in Input subfolder  <><","\n"))


# Plot MSY
Par = list(mfrow=c(1,1),mai=c(0.6,0.3,0,.15),omi = c(0.1,0.2,0.2,0) + 0.1,mgp=c(2,1,0), tck = -0.02,cex=0.8)
png(file = paste0(input.dir,"/Production",assessment,"_",Scenario,".png"), width = 6, height = 5, 
    res = 200, units = "in")
par(Par)

# Get Bmsy/B0 as a fucntion of M 
Bmsy=(m)^(-1/(m-1))
P = seq(0.0001,1,0.001) 
SP = ifelse(P>Plim,r.pr[1]/(m-1)*P*(1-P^(m-1)),r.pr[1]/(m-1)*P*(1-P^(m-1))*4*P)
#if(is.null(refBmsy)==TRUE) refBmsy = Bmsy
plot(P,SP/max(SP),type="l",ylab="Relative Yield",xlab="B/B0",lwd=2)
mtext(paste("Relative Yield"), side=2, outer=TRUE, at=0.6,line=1,cex=0.9)
legend("topright",c("SPM"),col=c(1),lwd=2,bty="n")  
abline(v=Bmsy,lty=2)
dev.off()

# Note PRIORS and save input subfolder
Priors =rbind(K.pr,psi.prior,c(r.pr[1],CV.r))
row.names(Priors) = c("K","Psi","r")
colnames(Priors) = c("Mean","CV")                          
write.csv(Priors,paste0(input.dir,"/Priors",assessment,"_",Scenario,".csv"))


#----------------------------------------------------------
# Set up JABBA 
#----------------------------------------------------------
cat(paste0("\n","><> Set up JAGS input <><","\n"))
# Plot MSY
# remove scientific numbers
options(scipen=999)
#----------------------------------------------------------


# starting values
nq = length(unique(sets.q))
nvar = length(unique(sets.var))

#----------------------------------------------------------
# Setup TAC projection
#---------------------------------------------------------
if(Projection==TRUE){
  nTAC = length(TACs)
  TAC = mat.or.vec(pyrs,nTAC)
  yr.last = max(years) # assessment year  
  
  for(i in 1:nTAC){
    TAC[,i] = c(rep(TACint,imp.yr-yr.last-1),rep(TACs[i],pyrs-(imp.yr-yr.last-1)))  
  }
  
}else{
  nTAC = 1  
  TAC = TC[n.years] #  
  pyrs = 1
}

#---------------------------------------------------------------
# JABBA Schaefer/Fox Models 1-2, Pella 3
#---------------------------------------------------------------

# Slope of hockey-stick
slope.HS = ifelse(Plim==0,1/10^-10,1/Plim)


nSel = 1 # setup for JABBA-SELECT version (in prep)
nI = ncol(CPUE) # number of CPUE series
stI = ifelse(proc.dev.all==TRUE,1, c(1:n.years)[is.na(apply(CPUE,1,mean,na.rm=TRUE))==FALSE][1]) #first year with CPUE


# Initial starting values
inits <- function(){list(K= rlnorm(1,log.K,0.3),q = runif(nq,0.005,0.5), isigma2.est=runif(1,20,100), itau2=runif(nvar,80,200))}
# starting value option
if(st.value==TRUE){
inits <- function(){list(K= K.init,r=r.init,q = q.init, isigma2.est=runif(1,20,100), itau2=runif(nvar,80,200))}
}

# JABBA input data 
surplus.dat = list(N=n.years, TC = TC,I=CPUE,SE2=se2,m=m,r.pr=r.pr,psi.pr=psi.pr,K.pr = K.pr,
                   nq=nq,nI = nI,nvar=nvar,sigma.fixed=ifelse(sigma.proc==TRUE,0,sigma.proc),
                   sets.var=sets.var, sets.q=sets.q,pen.bk = rep(0,n.years),Plim=Plim,slope.HS=slope.HS,
                   nTAC=nTAC,pyrs=pyrs,TAC=TAC,igamma = igamma,stI=stI,TACint =TACint,P_bound=P_bound,proc.pen=0)


# JAGS model file
JABBA = "JABBA.jags"


# PARAMETERS TO MONITOR
params <- c("K","r", "q", "psi","sigma2", "tau2","m","Hmsy","SBmsy", "MSY", "BtoBmsy","HtoHmsy","CPUE","Proc.Dev","P","SB","prP","prBtoBmsy","prHtoHmsy","TOE")


cat(paste0("\n","><> RUN ",Mod.names," model for ",assessment," ",Scenario," in JAGS <><","\n","\n"))

# JAGS MODEL Standard
sink("JABBA.jags")
cat("
    
    model {
    
    # Prior specifications  
    eps <- 0.0000000000000000000000000000000001 # small constant    
    
    #Catchability coefficients
    for(i in 1:nq)
    {   
    q[i] ~ dunif(eps,100)    
    }  
    
    
    ")

if(psi.dist =="beta"){
  cat("
      # Beta Prior for Biomass depletion at the start (deteministic)
      psi ~ dbeta(psi.pr[1],psi.pr[2])
      ",append=TRUE)
} else {
  cat("
      # Lognormal for Biomass depletion at the start (deteministic)
      psi ~ dlnorm(log(psi.pr[1]),pow(psi.pr[2],-2)) #I(0.1,1.1)    
      ",append=TRUE)  
}

if(sigma.proc==TRUE){
  cat("
      # Process variance
      isigma2 <- isigma2.est 
      sigma2 <- pow(isigma2,-1)
      sigma <- sqrt(sigma2)
      fakesigma.fixed <- sigma.fixed # Prevent unused variable error msg    
      ",append=TRUE)  
}else{ cat(" 
      isigma2 <- pow(sigma.fixed+eps,-2) 
           sigma2 <- pow(isigma2,-1)
           sigma <- sqrt(sigma2)
           
           ",append=TRUE)}

if(sigma.est==TRUE){
  cat("
      # Obsevation variance
      for(i in 1:nvar)
      {
      # Observation error
      itau2[i]~ dgamma(0.001,0.001)
      tau2[i] <- 1/itau2[i]
      }
      
      for(i in 1:nI)
      {
      for(t in 1:N)
      {
      var.obs[t,i] <- SE2[t,i]+tau2[sets.var[i]] 
      ivar.obs[t,i] <- 1/var.obs[t,i]
      # note total observation error (TOE)     
      TOE[t,i] <- sqrt(var.obs[t,i]) # Total observation variance
      
      }}
      ",append=TRUE)  
}else{ cat(" 
      # Obsevation variance
           for(i in 1:nvar)
           {
           # Observation error
           itau2[i]~ dgamma(2,2)
           tau2[i] <- 1/itau2[i]
           }
           
           for(i in 1:nI)
           {
           for(t in 1:N)
           {
           var.obs[t,i] <- SE2[t,i] # drop tau2
           fake.tau[t,i] <- tau2[sets.var[i]]
           
           ivar.obs[t,i] <- 1/var.obs[t,i]
           # note total observation error (TOE)     
           TOE[t,i] <- sqrt(var.obs[t,i])
           
           }}
           
           ",append=TRUE)}

# Run rest of code  
cat("  
    # Process variance prior
    isigma2.est ~ dgamma(igamma[1],igamma[2])
    
    # Carrying Capacity SB0
    K ~ dlnorm(log(K.pr[1]),pow(K.pr[2], -2))
    
    # informative priors for Hmsy as a function of r
    r ~ dlnorm(log(r.pr[1]),pow(r.pr[2],-2))
    
    
    #Process equation
    Pmean[1] <- log(psi)
    iPV[1] <- ifelse(1<(stI),10000,isigma2) # inverse process variance
    P[1] ~ dlnorm(Pmean[1],iPV[1]) # set to small noise instead of isigma2
    penB[1]  <- ifelse(P[1]<P_bound[1],log(K*P[1])-log(K*P_bound[1]),ifelse(P[1]>P_bound[2],log(K*P[1])-log(K*P_bound[2]),0)) # penalty if Pmean is outside viable biomass
    
    # Process equation
    for (t in 2:N) 
    {
    Pmean[t] <- ifelse(P[t-1] > Plim,
    log(max(P[t-1] +  r/(m-1)*P[t-1]*(1-pow(P[t-1],m-1)) - TC[t-1]/K,0.005)),
    log(max(P[t-1] +  r/(m-1)*P[t-1]*(1-pow(P[t-1],m-1))*P[t-1]*slope.HS - TC[t-1]/K,0.005)))
    iPV[t] <- ifelse(t<(stI),10000,isigma2) # inverse process variance
    P[t] ~ dlnorm(Pmean[t],iPV[t])
    penB[t]  <- ifelse(P[t]<(P_bound[1]),log(K*P[t])-log(K*(P_bound[1])),ifelse(P[t]>P_bound[2],log(K*P[t])-log(K*(P_bound[2])),0)) # penalty if Pmean is outside viable biomass
    }
    
    # Process error deviation 
    for(t in 1:N){
    Proc.Dev[t] <- P[t]-exp(Pmean[t])} 
    
    # Enforce soft penalties on bounds for P
    for(t in 1:N){
    pen.bk[t] ~ dnorm(penB[t],1000) # enforce penalty with CV = 0.1
    }
    
    Hmsy <- r*pow(m-1,-1)*(1-1/m) 
    
    for (t in 1:N) 
    { 
    SB[t] <- K*P[t]    
    H[t] <- TC[t]/SB[t] 
    }
    
    # Observation equation in related to EB
    
    for(i in 1:nI)
    {
    for (t in 1:N) 
    { 
    Imean[t,i] <- log(q[sets.q[i]]*P[t]*K);
    I[t,i] ~ dlnorm(Imean[t,i],(ivar.obs[t,i]));
    CPUE[t,i] <- q[sets.q[i]]*P[t]*K
    }}
    
    
    #Management quantaties
    SBmsy_K <- (m)^(-1/(m-1))
    SBmsy <- SBmsy_K*K
    
    MSY <- SBmsy*Hmsy
    for (t in 1:N)
    {
    # use x y to put them towards the end of the alphabetically sorted  mcmc object
    #SP[t] <- pow(r.pella,-(m-1))*SB[t]*(1-pow(P[t],m-1))
    BtoBmsy[t] <- SB[t]/SBmsy
    HtoHmsy[t] <- H[t]/(Hmsy) 
    }
    
    
    # Enforce soft penalty on process deviance if sigma.proc > 0.2 
    proc.pen ~ dnorm(penProc,1000) # enforce penalty 
    penProc  <- ifelse(sigma>0.2,log(sigma)-log(0.2),0) 
    
    
    ", append=TRUE)

# PROJECTION
if(Projection==TRUE){
  cat("
      for(i in 1:nTAC){
      # Project first year into the future
      prPmean[1,i] <- ifelse(P[N] > Plim,
      log(max(P[N] +  Hmsy/(1-1/m)*P[N]*(1-pow(P[N],m-1)) - TACint/K,0.005)),
      log(max(P[N] +  Hmsy/(1-1/m)*P[N]*(1-pow(P[N],m-1))*4*P[N] - TACint/K,0.005)))
      prP[1,i] ~ dlnorm(prPmean[1,i],isigma2) 
      # Project all following years
      for(t in 2:pyrs){
      prPmean[t,i] <- ifelse(prP[t-1,i] > Plim,
      log(max(prP[t-1,i] +  Hmsy/(1-1/m)*prP[t-1,i]*(1-pow(prP[t-1,i],m-1)) - TAC[t-1,i]/K,0.001)),
      log(max(prP[t-1,i] +  Hmsy/(1-1/m)*prP[t-1,i]*(1-pow(prP[t-1,i],m-1))*slope.HS*prP[t-1,i] - TAC[t-1,i]/K,0.005)))
      # process error (as monte-carlo simular)
      prP[t,i] ~ dlnorm(prPmean[t,i],isigma2)}
      for(t in 1:pyrs){
      prB[t,i] <- prP[t,i]*K
      prH[t,i] <- TAC[t,i]/prB[t,i]
      prHtoHmsy[t,i] <- prH[t,i]/Hmsy
      prBtoBmsy[t,i] <- prB[t,i]/SBmsy
      }}  
      ",append=TRUE)} else {
        cat("
            #Prevent error for unused input    
            fakeTAC <-  TAC
            fakepyrs <- pyrs 
            fakenTAC <- nTAC    
            fakeTACint <- TACint
            prHtoHmsy <- 1
            prP <- 1 
            prBtoBmsy <- 1    
            ", append=TRUE)}  

cat("
    
} # END OF MODEL
    ",append=TRUE,fill = TRUE)
sink()


ptm <- proc.time()

mod <- jags(surplus.dat, inits,params,paste(JABBA), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)  # adapt is burn-in

proc.time() - ptm
save.time = proc.time() - ptm

cat(paste0("\n",paste0("><> Scenario ",Scenario,"_",Mod.names," completed in ",as.integer(save.time[3]/60)," min and ",round((save.time[3]/60-as.integer(save.time[3]/60))*100)," sec <><","\n")))

cat(paste0("\n","><> Produce results output of ",Mod.names," model for ",assessment," ",Scenario," <><","\n"))

# if run with library(rjags)
posteriors = mod$BUGSoutput$sims.list



#-----------------------------------------------------------
# <><<><<><<><<><<><<>< Outputs ><>><>><>><>><>><>><>><>><>
#-----------------------------------------------------------
output.dir = paste0(File,"/",assessment,"/",Scenario,"_",Mod.names,"/Output")
dir.create(output.dir, showWarnings = FALSE)

# run some mcmc convergence tests
par.dat= data.frame(posteriors[params[c(1:6)]])
geweke = geweke.diag(data.frame(par.dat))
pvalues <- 2*pnorm(-abs(geweke$z))
pvalues

heidle = heidel.diag(data.frame(par.dat))

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

results = round(t(cbind(apply(par.dat,2,quantile,c(0.025,0.5,0.975)))),3)

results = data.frame(Median = results[,2],LCI=results[,1],UCI=results[,3],Geweke.p=round(pvalues,3),Heidel.p = round(heidle[,3],3))

ref.points = round(t(cbind(apply(man.dat,2,quantile,c(0.025,0.5,0.975)))),3)

ref.points = data.frame(Median = ref.points[,2],LCI=ref.points[,1],UCI=ref.points[,3])

# get number of parameters
npar = length(par.dat)
# number of years
N=n.years



# Safe posteriors (Produces large object!)
if(save.all==TRUE) save(posteriors,file=paste0(output.dir,"/",assessment,"_",Mod.names,"_",Scenario,"_posteriors"))

#-------------------------------------------------------------------------
# Save parameters, results table and current status posterior in csv files
#-------------------------------------------------------------------------

# Save model estimates and convergence p-values
write.csv(data.frame(results),paste0(output.dir,"/Estimates_",assessment,"_",Scenario,".csv"))

# Make standard results table with parameter estimates and reference points
Table = rbind(data.frame(results)[c("K","r","psi","sigma2"),1:3],data.frame(ref.points))  
Table[4,] = round(sqrt((Table[4,])),3) 
rownames(Table)[4] = "sigma.proc"
write.csv(Table,paste0(output.dir,"/Results_",assessment,"_",Scenario,".csv"))
#Save posterior of recent assessment year (KOBE posterior)
write.csv(data.frame(BtoBmsy=B_Bmsy.cur,FtoFmsy=H_Hmsy.cur),paste0(output.dir,"/Status_posterior",assessment,".csv"))  


#----------------
# Total Landings
#----------------
Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)

png(file = paste0(output.dir,"/Landings_",assessment,"_",Scenario,".png"), width = 5, height = 3.5, 
    res = 200, units = "in")
par(Par)

cord.x <- c(years,rev(years))
y<-rep(0,length(years))
plot(years,(TC),type="l",ylim=c(0,max(TC)),lty=1,lwd=1.3,xlab="Year",ylab=paste0("Catch ('000 t",catch.metric,")"),main="")
polygon(cord.x,c(TC,rev(y)),col="gray",border=1,lty=1)
dev.off()

#------------------------------
# Plot Posteriors
#------------------------------
sel.par = c(1,2,7,4,3,5)

out=data.frame(posteriors[params[sel.par]])
if(nSel>1) out=out[,-c(3:(3+nSel-2))]

node_id = names(out)
#informative priors
Prs = as.matrix(cbind(K.pr,r.pr,c(0,0),psi.pr))

#Posteriors
Par = list(mfrow=c(round(length(node_id)/3+0.33,0),3),mai=c(0.4,0.1,0,.1),omi = c(0.3,0.5,0.1,0) + 0.1,mgp=c(1,0.1,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/Posteriors_",assessment,"_",Scenario,".png"),width  = 8, height = 2.5*round(length(node_id)/3,0), 
    res = 200, units = "in")
par(Par)


node_id = names(out)


#par(mfrow=c(4,2),oma=c(0,1,1,0), mar=c(4,4,1,1))

for(i in 1:length(node_id))
{
  
  post.par = as.numeric(unlist(out[paste(node_id[i])]))
  
  if(i==1){
    
    rpr = rlnorm(10000,log(K.pr[1]),K.pr[2]) 
    pdf = stats::density(post.par,adjust=2)  
    prior = dlnorm(sort(rpr),log(K.pr[1]),K.pr[2])   
    plot(pdf,type="l",ylim=range(prior,pdf$y),xlim=range(c(pdf$x,quantile(rpr,c(0.0001,0.95)))),yaxt="n",xlab="K",ylab="",xaxs="i",yaxs="i",main="")
    
    polygon(c(sort(rpr),rev(sort(rpr))),c(prior,rep(0,length(sort(rpr)))),col=gray(0.4,1))
    polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
    legend('topright',c("Prior","Posterior"),pch=22,pt.cex=1.5,pt.bg = c(grey(0.4,1),grey(0.8,0.6)),bty="n")
  }  
  
  
  if(i==2){
    
    rpr = rlnorm(10000,log(Prs[1,i]),Prs[2,i]) 
    pdf = stats::density(post.par,adjust=2) 
    prior = dlnorm(sort(rpr),log(Prs[1,i]),Prs[2,i])   
    plot(pdf$x,pdf$y,type="l",ylim=range(prior,pdf$y),xlim=range(c(post.par,quantile(rpr,c(0.0001,0.95)))),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i")
    
    polygon(c(sort(rpr),rev(sort(rpr))),c(prior,rep(0,length(sort(rpr)))),col=gray(0.4,1))
    polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
  }
  
  if(i==3){
    plot(1,1,type="n",xlim=range(0.5,2.5),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i")  
    abline(v=m,lwd=2)}
  
  
  
  if(i==4){
    if(psi.dist=="beta"){
      parm = fitdist(post.par[post.par<1 & post.par>0.01], "beta")$estimate
      rpr = rbeta(10000,(psi.pr[1]),psi.pr[2]) 
      pdf = stats::density(post.par,adjust=2)  
      prior = dbeta(sort(rpr),psi.pr[1],psi.pr[2])   
    } else {
      rpr = rlnorm(10000,log(psi.prior[1]),psi.prior[2]) 
      pdf = stats::density(post.par,adjust=2)  
      prior = dlnorm(sort(rpr),log(psi.prior[1]),psi.prior[2])}
    plot(pdf,type="l",ylim=range(quantile(c(prior,pdf$y,c(0,0.95)))),xlim=range(c(0.5,post.par,pdf$x,quantile(rpr,c(0.001,0.999)))),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i",main="")
    polygon(c(sort(rpr),rev(sort(rpr))),c(prior,rep(0,length(sort(rpr)))),col=gray(0.4,1))
    polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
    #legend('topright',c("Prior","Posterior"),pch=22,pt.cex=1.5,pt.bg = c(grey(0.4,1),grey(0.8,0.6)),bty="n")
  }        
  
  if(i>4){
    if(sigma.proc!=TRUE & i==length(node_id)) {
      plot(1,1,type="n",xlim=range(0,0.15^2),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i")  
      abline(v=sigma.proc^2,lwd=2)} else {
        
        pdf = stats::density(post.par,adjust=2)  
        plot(pdf,type="l",xlim=range(0,post.par),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i",main="")
        if(i==length(node_id)& igamma[1]>0.9){
          rpr = 1/rgamma(10000,igamma[1],igamma[2])
          prior = stats::density(rpr,adjust=2)
          polygon(c(prior$x,rev(prior$x)),c(prior$y,rep(0,length(prior$y))),col=gray(0.4,1))
        }
        
        polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
        #legend('topright',c("Posterior"),pch=22,pt.cex=1.5,pt.bg = c(grey(0.8,0.6)),bty="n")
      } }         
  
}
mtext(paste("Density"), side=2, outer=TRUE, at=0.5,line=1,cex=0.9)
dev.off()   



#-----------------------------
# MCMC chains of posteriors
#-----------------------------

Par = list(mfrow=c(round(length(node_id)/3+0.33,0),3),mai=c(0.4,0.1,0,.1),omi = c(0.3,0.5,0.1,0) + 0.1,mgp=c(1,0.1,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/MCMC_",assessment,"_",Scenario,".png"), width = 8, height = 2.5*round(length(node_id)/3,0), 
    res = 200, units = "in")
par(Par)
for(i in 1:length(node_id)){
  
  post.par = as.numeric(unlist(out[paste(node_id[i])]))
  plot(out[,i],xlab=paste(node_id[i]),ylab="",type="l",col=4)
  lines(rep(mean(out[,i]),length(out[,i])),col=2,lwd=2)   
}
dev.off()

#-----------------------------------------------------------
# <><<><<><<>< Produce JABBA Model Diagnostics ><>><>><>><>
#-----------------------------------------------------------
cat(paste0("\n","><> Producing JABBA Model Fit Diagnostics <><","\n"))


# extract predicted CPUE + CIs

N = n.years
series = 1:n.indices

check.yrs = apply(CPUE,1,sum,na.rm=TRUE)
cpue.yrs = years[check.yrs>0]

#CPUE FITS
Par = list(mfrow=c(round(n.indices/2+0.01,0),ifelse(n.indices==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.3,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/Fits_",assessment,"_",Scenario,".png"), width = 7, height = ifelse(n.indices==1,5,2.5)*round(n.indices/2+0.01,0), 
    res = 200, units = "in")
par(Par)


for(i in 1:n.indices){
  
  # set observed vs predicted CPUE
  #par(mfrow=c(1,1))
  Yr = years
  Yr = min(Yr):max(Yr)
  yr = Yr-min(years)+1
  
  fit = apply(posteriors$CPUE[,,i],2,quantile,c(0.025,0.5,0.975))
  mufit = mean(fit[2,])
  fit = fit/mufit
  cpue.i = CPUE[is.na(CPUE[,i])==F,i]
  yr.i = Yr[is.na(CPUE[,i])==F]
  se.i = sqrt(se2[is.na(CPUE[,i])==F,(i)])
  
  ylim = c(min(fit[,check.yrs>0]*0.8,exp(log(cpue.i)-1.96*se.i)/mufit), max(fit[,yr[check.yrs>0]]*1.4,exp(log(cpue.i)+1.96*se.i)/mufit))
  
  cord.x <- c(Yr,rev(Yr))
  cord.y <- c(fit[1,yr],rev(fit[3,yr]))
  
  # Plot Observed vs predicted CPUE
  plot(years,CPUE[,i],ylab="",xlab="",ylim=ylim,xlim=range(cpue.yrs),type='n',xaxt="n",yaxt="n")
  axis(1,labels=TRUE,cex=0.8)
  axis(2,labels=TRUE,cex=0.8)
  polygon(cord.x,cord.y,col=grey(0.5,0.5),border=0,lty=2)
  
  
  lines(Yr,fit[2,yr],lwd=2,col=1)
  if(SE.I ==TRUE | max(se2)>0.01){ plotCI(yr.i,cpue.i/mufit,ui=exp(log(cpue.i)+1.96*se.i)/mufit,li=exp(log(cpue.i)-1.96*se.i)/mufit,add=T,gap=0,pch=21,xaxt="n",yaxt="n")}else{
    points(yr.i,cpue.i/mufit,pch=21,xaxt="n",yaxt="n",bg="white")}
  
  legend('topright',paste(indices[i]),bty="n",y.intersp = -0.2,cex=0.9)
}

mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
mtext(paste("Normalized Index"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
dev.off()

#log CPUE FITS
Par = list(mfrow=c(round(n.indices/2+0.01,0),ifelse(n.indices==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.3,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/logFits_",assessment,"_",Scenario,".png"), width = 7, height = ifelse(n.indices==1,5,2.5)*round(n.indices/2+0.01,0), 
    res = 200, units = "in")
par(Par)


for(i in 1:n.indices){
  
  # set observed vs predicted CPUE
  #par(mfrow=c(1,1))
  Yr = years
  Yr = min(Yr):max(Yr)
  yr = Yr-min(years)+1
  
  fit = apply(posteriors$CPUE[,,i],2,quantile,c(0.025,0.5,0.975))
  mufit = mean(fit[2,])
  fit = fit/mufit
  cpue.i = CPUE[is.na(CPUE[,i])==F,i]
  yr.i = Yr[is.na(CPUE[,i])==F]
  se.i = sqrt(se2[is.na(CPUE[,i])==F,(i)])
  
  ylim = log(c(min(fit[,yr[is.na(CPUE[,i])==F]]*0.8,exp(log(cpue.i)-1.96*se.i)/mufit), max(fit[,yr[is.na(CPUE[,i])==F]]*1.3,exp(log(cpue.i)+1.96*se.i)/mufit)))
  
  cord.x <- c(Yr,rev(Yr))
  cord.y <- log(c(fit[1,yr],rev(fit[3,yr])))
  
  # Plot Observed vs predicted CPUE
  plot(years,CPUE[,i],ylab="",xlab="",ylim=ylim,xlim=range(yr.i),type='n',xaxt="n",yaxt="n")
  axis(1,labels=TRUE,cex=0.8)
  axis(2,labels=TRUE,cex=0.8)
  #polygon(cord.x,cord.y,col=grey(0.5,0.5),border=0,lty=2)
  
  
  lines(Yr,log(fit[2,yr]),lwd=2,col=4)
  if(SE.I ==TRUE | max(se2)>0.01){ plotCI(yr.i,log(cpue.i/mufit),ui=log(exp(log(cpue.i)+1.96*se.i)/mufit),li=log(exp(log(cpue.i)-1.96*se.i)/mufit),add=T,gap=0,pch=21,xaxt="n",yaxt="n")}else{
    points(yr.i,log(cpue.i/mufit),pch=21,xaxt="n",yaxt="n",bg="white")}
  legend('topright',paste(indices[i]),bty="n",y.intersp = -0.2,cex=0.8)
}

mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
mtext(paste("Log Index"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
dev.off()



# JABBA-residual plot
Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/Residuals_",assessment,"_",Scenario,".png"), width = 5, height = 3.5, 
    res = 200, units = "in")
par(Par)

Resids = NULL
for(i in 1:n.indices){
  Resids =rbind(Resids,log(CPUE[,i])-log(apply(posteriors$CPUE[,,i],2,quantile,c(0.5))))   
}

plot(Yr,Yr,type = "n",ylim=ifelse(rep(max(Resids,na.rm = T),2)>0.9,range(1.2*Resids,na.rm = T),range(c(-1.3,1.2))),xlim=range(cpue.yrs),ylab="log residuals",xlab="Year")
boxplot(Resids,add=TRUE,at=c(Yr),xaxt="n",col=grey(0.8,0.5),notch=FALSE,outline = FALSE)
abline(h=0,lty=2)
positions=runif(n.indices,-0.2,0.2)

for(i in 1:n.indices){
  for(t in 1:n.years){
    lines(rep((Yr+positions[i])[t],2),c(0,Resids[i,t]),col=jabba.colors[i])}
  points(Yr+positions[i],Resids[i,],col=1,pch=21,bg=jabba.colors[i])}
mean.res = apply(Resids,2,mean,na.rm =TRUE)
smooth.res = predict(loess(mean.res~Yr),data.frame(Yr=cpue.yrs))
lines(cpue.yrs,smooth.res,lwd=2)
DIC =round(mod$BUGSoutput$DIC,1)
# get degree of freedom
Nobs =length(as.numeric(Resids)[is.na(as.numeric(Resids))==FALSE])
DF = Nobs-npar

RMSE = round(100*sqrt(sum(Resids^2,na.rm =TRUE)/DF),1)

legend('topright',c(paste0("RMSE = ",RMSE,"%")),bty="n")
legend('bottomright',c(paste(indices),"Loess"),bty="n",col=1,pt.cex=1.1,cex=0.75,pch=c(rep(21,n.indices),-1),pt.bg=c(jabba.colors[series],1),lwd=c(rep(-1,n.indices),2))

dev.off()

#Save Residuals 
Res.CPUE = data.frame(Resids)
row.names(Res.CPUE) = indices   
colnames(Res.CPUE) = paste(Yr)
write.csv(Res.CPUE,paste0(output.dir,"/ResCPUE_",assessment,"_",Scenario,".csv"))

#---------------------------------------
# Stadardized Residuals
#--------------------------------------
Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/StandardizedResids_",assessment,"_",Scenario,".png"), width = 5, height = 3.5, 
    res = 200, units = "in")
par(Par)


# Standardized Residuals
StResid = NULL
for(i in 1:n.indices){
  StResid =rbind(StResid,log(CPUE[,i]/apply(posteriors$CPUE[,,i],2,quantile,c(0.5)))/
                   apply(posteriors$TOE[,,i],2,quantile,c(0.5))+0.5*apply(posteriors$TOE[,,i],2,quantile,c(0.5)))        
}

plot(Yr,Yr,type = "n",ylim=c(min(-1,-1.2*max(abs(StResid),na.rm = T)),max(1,1.2*max(abs(StResid),na.rm = T))),xlim=range(cpue.yrs),ylab="Standardized residuals",xlab="Year")
boxplot(StResid,add=TRUE,at=c(Yr),xaxt="n",col=grey(0.8,0.5),notch=FALSE,outline = FALSE)
abline(h=0,lty=2)
positions=runif(n.indices,-0.2,0.2)

for(i in 1:n.indices){
  for(t in 1:n.years){
    lines(rep((Yr+positions[i])[t],2),c(0,StResid[i,t]),col=jabba.colors[i])}
  points(Yr+positions[i],StResid[i,],col=1,pch=21,bg=jabba.colors[i])}
mean.res = apply(StResid,2,mean,na.rm =TRUE)
smooth.res = predict(loess(mean.res~Yr),data.frame(Yr=cpue.yrs))
lines(cpue.yrs,smooth.res,lwd=2)
DIC =round(mod$BUGSoutput$DIC,1)
SDNR = round(sqrt(sum(StResid^2,na.rm =TRUE)/(Nobs-1)),2)
Crit.value = (qchisq(.95, df=(Nobs-1))/(Nobs-1))^0.5
legend('topright',c(paste0("SDNR = ",SDNR,"(",round(Crit.value,2),")")),bty="n")
legend('bottomright',c(paste(indices),"Loess"),bty="n",col=1,cex=0.75,pt.cex=1.1,pch=c(rep(21,n.indices),-1),pt.bg=c(jabba.colors[series],1),lwd=c(rep(-1,n.indices),2))


dev.off()

#Save standardized Residuals 
StRes.CPUE = data.frame(StResid)
row.names(Res.CPUE) = indices   
colnames(Res.CPUE) = paste(Yr)
write.csv(Res.CPUE,paste0(output.dir,"/StResCPUE_",assessment,"_",Scenario,".csv"))

# Produce statistice describing the Goodness of the Fit
GOF = data.frame(Stastistic = c("N","p","DF","SDNR","RMSE","DIC"),Value = c(Nobs,npar,DF,SDNR,RMSE,DIC))
write.csv(GOF,paste0(output.dir,"/GOF_",assessment,"_",Scenario,".csv"))

#------------------------------
# Plot process error deviation
#------------------------------

proc.dev = apply(posteriors$Proc.Dev,2,quantile,c(0.025,0.5,0.975))

Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)

png(file = paste0(output.dir,"/ProcDev_",assessment,"_",Scenario,".png"), width = 5, height = 3.5, 
    res = 200, units = "in")
par(Par)

ylim = range(proc.dev)*1.1
cord.x <- c(years,rev(years))
cord.y <- c(proc.dev[1,],rev(proc.dev[3,]))
# Process Error
plot(years,proc.dev[2,],ylab="Process deviation P[t]",xlab="Year",ylim=ylim,type="n")
polygon(cord.x,cord.y,col='grey',border=0,lty=2)
lines(years,proc.dev[2,],lwd=2)
lines(years,rep(0,length(years)),lty=5)


dev.off()

#-----------------------------------------------------------
# <><<><<><<><<>< JABBA Management Plots ><>><>><>><>><>><>
#-----------------------------------------------------------
cat(paste0("\n","><> Producing Management Plots <><","\n"))

#-----------------------
# Plot Biomass B_t
#-----------------------

B_t = posteriors$SB
mu.B = apply(B_t,2,quantile,c(0.025,0.5,0.975))

Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/Biomass_",assessment,"_",Scenario,".png"), width = 5, height = 3.5, 
    res = 200, units = "in")
par(Par)

ylim = c(0, max(mu.B ))
cord.x <- c(years,rev(years))
cord.y <- c(mu.B [1,],rev(mu.B [3,]))

# B_t
plot(years,mu.B[2,],ylab=paste0("Biomass ",catch.metric),xlab="Year",ylim=ylim,type="n")
polygon(cord.x,cord.y,col='grey',border=0,lty=2)
lines(years,mu.B[2,],lwd=2,col=1)
lines(years,rep(mean(posteriors$SBmsy),length(years)),lty=5)
text((max(years)-min(years))/30+years[1],mean(posteriors$SBmsy)*1.11,expression(paste(B[MSY])))
dev.off()



#-----------------------
# Plot B/Bmsy and H/Hmsy
#-----------------------

HtoHmsy = posteriors$HtoHmsy
BtoBmsy = posteriors$BtoBmsy

mu.f = apply(HtoHmsy,2,quantile,c(0.025,0.5,0.975))
mu.b = apply(BtoBmsy,2,quantile,c(0.025,0.5,0.975))

f = HtoHmsy[,N]
b = BtoBmsy[,N]
Par = list(mfrow=c(1,2),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/TrendMSY_",assessment,"_",Scenario,".png"), width = 7, height = 3, 
    res = 200, units = "in")
par(Par)

ylim = c(0, max(mu.f))
cord.x <- c(years,rev(years))
cord.y <- c(mu.f[1,],rev(mu.f[3,]))

# H/Hmsy
plot(years,mu.f[2,],ylab=ifelse(harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))),xlab="Year",ylim=ylim,type="n")
polygon(cord.x,cord.y,col='grey',border=0,lty=2)
lines(years,mu.f[2,],lwd=2,col=1)
lines(years,rep(1,length(years)),lty=5)

ylim = c(0, max(mu.b,1.1))

cord.x <- c(years,rev(years))
cord.y <- c(mu.b[1,],rev(mu.b[3,]))
plot(years,mu.b[2,],ylab=expression(paste(B/B[MSY])),xlab="Year",ylim=ylim,type="n")
polygon(cord.x,cord.y,col='grey',border=0,lty=2)
lines(years,mu.b[2,],lwd=2,col=1)
lines(years,rep(1,length(years)),lty=5)
dev.off()

if(SP.plot!="phase"){
  #-----------------------------------------
  # Produce simple Production function plot
  #-----------------------------------------
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/SP_",assessment,"_",Scenario,".png"), width = 5, height = 4.5, 
      res = 200, units = "in")
  par(Par)
  Bit = seq(1,mean(posteriors$K),mean(posteriors$K)/500)
  SP = mean(posteriors$r)/(m-1)*Bit*(1-(Bit/median(posteriors$K))^(m-1))  
  B = apply(posteriors$SB,2,mean)
  MSY = quantile(posteriors$MSY,c(0.025,0.5,0.975)) 
  plot(Bit,SP,type = "n",ylim=c(0,max(c(max(Catch,na.rm=T),max(MSY*1.1)))),ylab=paste0("Surplus Production ",catch.metric),xlab="Biomass (t)")
  polygon(c(-10000,10^7,10^7,-10000),c(rep(MSY[1],2),rep(MSY[3],2)),border = 0,col=grey(0.5,0.4))
  lines(Bit,SP,col=2,lwd=2)
  lines(B,Catch,lty=1)
  points(B,Catch,cex=0.5,pch=4)
  sel.yr = c(1,round(quantile(1:N,0.7),0),N)
  points(B[sel.yr],Catch[sel.yr],col= 1,pch=c(22,21,24),bg="white",cex=1.5)
  abline(h=max(SP),col=4)
  sel.years = c(min(years),years[sel.yr[2]],max(years))
  lines(rep(median(posteriors$SBmsy),2),c(-1000,max(SP)),lty=2,col=2)
  legend('topright', 
         c(expression(B[MSY]),"MSY","Catch",paste(sel.years)), 
         lty=c(2,1,1,1,1,1),pch=c(-1,-1,4,22,21,24),pt.bg=c(0,0,0,rep("white",3)), 
         col=c(2,4,rep(1,4)),lwd=1,cex=0.9,pt.cex=c(-1,-1,0.5,rep(1.3,3)),bty="n")
  
  dev.off()
}

if(SP.plot=="phase"){
  #-----------------------------------------
  # Produce JABBA SP-phase plot
  #-----------------------------------------
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/SPphase_",assessment,"_",Scenario,".png"), width = 5, height = 4.5, 
      res = 200, units = "in")
  par(Par)
  Bit = seq(1,median(posteriors$K),median(posteriors$K)/500)
  Cmsy = Bit*median(posteriors$Hmsy)
  SP = median(posteriors$r)/(m-1)*Bit*(1-(Bit/median(posteriors$K))^(m-1))  
  B = apply(posteriors$SB,2,mean)
  MSY = quantile(posteriors$MSY,c(0.025,0.5,0.975)) 
  Bmsy.sp = median(posteriors$SBmsy)
  K.sp = median(posteriors$K)
  green.x = c(max(Bit,B),max(Bit,B),Bmsy.sp,Bmsy.sp,max(Bit))
  green.y = c(Bmsy.sp,0,0,max(SP),max(Cmsy))
  red.x = c(0,0,Bmsy.sp,Bmsy.sp,0)
  red.y = c(K.sp,0,max(SP),K.sp,K.sp)
  plot(Bit,SP,type = "n",ylim=c(0,max(c(max(Catch,na.rm=T)*1.05,max(MSY*1.1)))),xlim=c(0,max(Bit,B)),ylab=paste0("Surplus Production ",catch.metric),xlab="Biomass (t)",xaxs="i",yaxs="i")
  rect(0,0,K.sp*1.1,K.sp*1.1,col="green",border=0)
  rect(0,0,K.sp,K.sp,col="yellow",border=0)
  if(KOBE.type!="ICCAT") rect(0,max(SP),K.sp,K.sp,col="orange",border=0)
  polygon(green.x,green.y,border = 0,col="green")
  polygon(red.x,red.y,border = 0,col="red")
  
  ry.sp = Bit[Bit<=Bmsy.sp]
  for(i in 1:length(ry.sp)){
    
    lines(rep(Bit[i],2),c(Cmsy[i],SP[i]),col=ifelse(i %% 2== 0,"yellow","red"),lty=3)  
    #i = i+1
  }
  
  gy.sp = Bit[Bit>Bmsy.sp]
  for(i in (length(ry.sp)+1):length(Bit)){
    
    #lines(rep(Bit[i],2),c(max(SP),Cmsy[i]),col=ifelse(i %% 2== 0,ifelse(KOBE.type=="ICCAT","yellow","orange"),"green"),lty=3)  
    #i = i+1
  }
  
  
  polygon(c(-10000,10^7,10^7,-10000),c(rep(MSY[1],2),rep(MSY[3],2)),border = FALSE,col=rgb(0,0,1,0.4))
  lines(Bit,SP,col=4,lwd=2)
  lines(B,Catch,lty=1,lwd=1)
  points(B,Catch,cex=0.8,pch=16)
  lines(Bit,Cmsy,col=1,lwd=1,lty=2)
  N=n.years
  sel.yr = c(1,round(quantile(1:N,0.7),0),N)
  points(B[sel.yr],Catch[sel.yr],col= 1,pch=c(22,21,24),bg="white",cex=1.7)
  abline(h=max(SP),col=4,lty=5)
  sel.years =years[sel.yr]
  lines(rep(median(posteriors$SBmsy),2),c(-1000,max(SP)),lty=2,col=4)
  
  legend('topright', 
         c(expression(B[MSY]),"MSY","SP","Catch",paste(sel.years)), 
         lty=c(2,5,1,1,1,1,1),pch=c(-1,-1,-1,16,22,21,24),pt.bg=c(0,0,0,0,rep("white",3)), 
         col=c(4,4,4,rep(1,4)),lwd=c(1,1,2,1,1,1),cex=0.8,pt.cex=c(-1,-1,-1,0.5,rep(1.3,3)),bty="n")
  
  dev.off()
}


#----------------------
# Produce Kobe plot
#----------------------
# extract vectors BtoBmsy and FtoFmsy


if(KOBE.plot==TRUE){
  # prepare 
  
  # fit kernel function
  kernelF <- ci2d(b,f,nbins=151,factor=1.5,ci.levels=c(0.50,0.80,0.75,0.90,0.95),show="none",col=1,xlab= ifelse(harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))),ylab=expression(paste(B/B[MSY])))
  
  
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/Kobe_",assessment,"_",Scenario,".png"), width = 5, height = 4.5, 
      res = 200, units = "in")
  par(Par)
  
  #Create plot
  plot(1000,1000,type="b", xlim=c(0,2.5), ylim=c(0,max(apply(HtoHmsy,2,quantile,c(0.5)),quantile(f,0.85),2.)),lty=3,ylab=ifelse(harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))),xlab=expression(paste(B/B[MSY])),xaxs="i",yaxs="i")
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
  polygon(c(1,100,100,1),c(1,1,100,100),col=ifelse(KOBE.type=="ICCAT","yellow","orange"),border=0)
  polygon(c(0,1,1,0),c(1,1,100,100),col="red",border=0)
  
  polygon(kernelF$contours$"0.95",lty=2,border=NA,col="cornsilk4")
  polygon(kernelF$contours$"0.8",border=NA,lty=2,col="grey")
  polygon(kernelF$contours$"0.5",border=NA,lty=2,col="cornsilk2")
  points(mu.b[2,],mu.f[2,],pch=16,cex=1)
  
  
  lines(c1,c2,lty=3,lwd=0.7)
  lines(c2,c1,lty=3,lwd=0.7)
  lines(mu.b[2,],mu.f[2,], lty=1,lwd=1.)
  sel.yr = c(1,round(quantile(1:N,0.7),0),N)
  points(mu.b[2,sel.yr],mu.f[2,sel.yr],col=
           1,pch=c(22,21,24),bg="white",cex=1.9)
  
  # Get Propability
  Pr.green = sum(ifelse(b>1 & f<1,1,0))/length(b)*100
  Pr.red = sum(ifelse(b<1 & f>1,1,0))/length(b)*100
  
  if(KOBE.type=="ICCAT"){               
    Pr.yellow = (sum(ifelse(b<1 & f<1,1,0))+sum(ifelse(b>1 & f>1,1,0)))/length(b)*100} else {
      Pr.yellow = sum(ifelse(b<1 & f<1,1,0))/length(b)*100
      Pr.orange = sum(ifelse(b>1 & f>1,1,0))/length(b)*100
    }
  
  
  sel.years = c(years[sel.yr])
  ## Add legend
  if(KOBE.type=="ICCAT"){
    legend('topright', 
           c(paste(sel.years),"50% C.I.","80% C.I.","95% C.I.",paste0(round(c(Pr.red,Pr.yellow,Pr.green),1),"%")), 
           lty=c(1,1,1,rep(-1,7)),pch=c(22,21,24,rep(22,7)),pt.bg=c(rep("white",3),"cornsilk2","grey","cornsilk4","red","yellow","green"), 
           col=1,lwd=1.1,cex=0.9,pt.cex=c(rep(1.3,3),rep(1.7,3),rep(2.1,3)),bty="n")
  }else{
    legend('topright', 
           c(paste(sel.years),"50% C.I.","80% C.I.","95% C.I.",paste0(round(c(Pr.red,Pr.yellow,Pr.orange,Pr.green),1),"%")), 
           lty=c(1,1,1,rep(-1,8)),pch=c(22,21,24,rep(22,8)),pt.bg=c(rep("white",3),"cornsilk2","grey","cornsilk4","red","yellow","orange","green"), 
           col=1,lwd=1.1,cex=0.9,pt.cex=c(rep(1.3,3),rep(1.7,3),rep(2.2,4)),bty="n")  
    
  }
  dev.off()
}


if(Biplot==TRUE){
  
  #---------------------------------------------------------
  # Produce 'post-modern' biplot (see Quinn and Collie 2005)
  #---------------------------------------------------------
  
  # read ftarget,bthreshold
  ftarget<-0.8
  bthreshold<-0.2
  
  # fit kernel function
  kernelF <- ci2d(f,b,nbins=201,factor=2,ci.levels=c(0.50,0.80,0.75,0.90,0.95),show="none",col=1,ylab= ifelse(harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))),xlab=expression(paste(B/B[MSY])))
  
  
  Par = list(mfrow=c(1,1),mai=c(0.2,0.15,0,.15),omi = c(0.3,0.25,0.2,0) + 0.1, mgp =c(3,1,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/Biplot_",assessment,"_",Scenario,".png"), width = 5, height = 4.5, 
      res = 200, units = "in")
  par(Par)
  
  #Create plot
  plot(1000,1000,type="b", ylim=c(0,2.5), xlim=c(0,max(apply(HtoHmsy,2,quantile,c(0.5)),quantile(f,0.85),2.)),lty=3,xaxs="i",yaxs="i")
  
  # and fill areas using the polygon function
  fint = seq(0.001,100,0.01)
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
  lines(mu.f[2,],mu.b[2,], lty=1,lwd=1.)
  sel.yr = c(1,round(quantile(1:N,0.7),0),N)
  points(mu.f[2,sel.yr],mu.b[2,sel.yr],col=
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
      } else {
        Zone[i]<-Y
      }
    }}
  
  perGreen = round(length(Zone[Zone==0.15])/length(Zone)*100,1) 
  perYellow = round(length(Zone[Zone==0])/length(Zone)*100,1) 
  perRed = round(length(Zone[Zone==-0.15])/length(Zone)*100,1)
  
  mtext(expression(paste(B/B[MSY])), side=2, outer=TRUE, at=0.5,line=1,cex=0.9)
  mtext(ifelse(harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))), side=1, outer=TRUE, at=0.5,line=1,cex=0.9)
  
  text(0.65,2.4,paste0(perGreen,"%"))
  text(0.9,2.4,paste0(perYellow,"%"))
  text(1.2,2.4,paste0(perRed,"%"))
  
  dev.off()
  
}



#--------------------------------------
# Plot projections and safe posteriors
#--------------------------------------
if(Projection ==TRUE){
  cat(paste0("\n","><> Producing Future TAC Projections <><","\n"))
  
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/Projections_",assessment,"_",Scenario,".png"), width = 5, height = 3.5, 
      res = 200, units = "in")
  par(Par)
  
  proj.yrs =  years[n.years]:(years[n.years]+pyrs)
  # Dims 1: saved MCMC,2: Years, 3:alternatic TACs, 4: P, H/Hmsy, B/Bmsy
  projections = array(NA,c(nsaved,length(proj.yrs),nTAC,3))
  for(i in 1:nTAC){
    projections[,,i,1] = cbind(posteriors$P[,(n.years):n.years],posteriors$prP[,,i])
  }
  
  for(i in 1:nTAC){
    projections[,,i,2] = cbind(posteriors$BtoBmsy[,(n.years):n.years],posteriors$prBtoBmsy[,,i])
  }
  for(i in 1:nTAC){
    projections[,,i,3] = cbind(posteriors$HtoHmsy[,(n.years):n.years],posteriors$prHtoHmsy[,,i])
  }
  
  kjp = kobeJabbaProj(projections,proj.yrs[1])
  
  save(kjp,file=paste0(output.dir,"/",assessment,"_",Mod.names,"_",Scenario,"_projections"))
  
  # Change here for ICCAT Bmsy plot
  Traj = cbind(posteriors$P[,(n.years):n.years],posteriors$prP[,,nTAC])
  
  #plot(proj.yrs,apply(Traj,2,mean),ylim=c(0,1),xlim=c(min(proj.yrs),max(proj.yrs)+length(proj.yrs)*0.2),type="n",ylab="Biomass depletion (B/K)",xlab="Projection Years")
  plot(proj.yrs,apply(Traj,2,mean),ylim=c(0,1),xlim=c(min(proj.yrs),max(proj.yrs)),type="n",ylab="Biomass depletion (B/K)",xlab="Projection Years")
  
  cols = rev(seq(0.4,0.9,0.5/nTAC))
  plot.order = (1:nTAC)
  for(j in 1:(nTAC)){
    i =  plot.order[j]
    Traj = cbind(posteriors$P[,(n.years):n.years],posteriors$prP[,,i])
    polygon(c(proj.yrs,rev(proj.yrs)),c(apply(Traj,2,quantile,0.05),rev(apply(Traj,2,quantile,0.95))),col=grey(cols[i],1),border=NA)
  }
  for(j in 1:(nTAC)){
    i =  plot.order[j]
    Traj = cbind(posteriors$P[,(n.years):n.years],posteriors$prP[,,i])
    lines(proj.yrs,apply(Traj,2,median),col=rev(jabba.colors[j]),lwd=2)
  }
  lines(proj.yrs[1:(imp.yr-yr.last+1)],apply(Traj,2,median)[1:(imp.yr-yr.last+1)],col=1,lwd=2)
  BmsyK=(m)^(-1/(m-1))
  abline(h=BmsyK,lty=2,lwd=2)
  legend("topleft",paste(TACs,"(t)"),col=(jabba.colors[1:nTAC]),lwd=2,cex=0.8)   
  
  dev.off()    
  
}


if(save.trajectories==TRUE){
  cat(paste0("\n","><> Saving Posteriors of FRP trajectories <><","\n"))
  
  # FRP trajectories
  trajectories = array(NA,c(nsaved,n.years,3))
  trajectories[,,1] = posteriors$P 
  trajectories[,,2] = posteriors$BtoBmsy 
  trajectories[,,3] = posteriors$HtoHmsy
  
  kb=kobeJabba(trajectories,years[1])
  save(kb,file=paste0(output.dir,"/",assessment,"_",Mod.names,"_",Scenario,"_trajectories"))
  
}


cat(paste0("\n","\n","><> Scenario ",Mod.names,"_",Scenario," for ",assessment," - DONE! <><","\n"))





