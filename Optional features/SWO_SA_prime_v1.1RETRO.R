##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><><
## JABBA: Just Another Bayesian Biomass Assessment
## Input File for JABBA
## Developed by Henning Winker & Felipe Carvalho (Cape Town/Hawaii)  
##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  
# required packages
library(gplots)
library(coda)
library(rjags)
library(R2jags)
library("fitdistrplus")
library(reshape)

#----------------------------------------------------------------
# Setup working directories and output folder labels 
#-----------------------------------------------------------------
# Set Working directory file, where assessments are stored 
File = "C:/Work/Research/MS_JABBA/R1/final"
# Set working directory for JABBA R source code
JABBA.file = "C:/Work/Research/MS_JABBA/R1/final"
# JABBA version
version = "v1.1"
# Set Assessment file: assement folder within File that includes .csv input files
assessment = "SWO_SA" 
# add specifier for assessment (File names of outputs)

#----------------------------------------------------------------


#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Graphic, Output, Saving (.RData) settings 
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
KOBE.plot = TRUE # Produces JABBA Kobe plot 
KOBE.type = c("ICCAT","IOTC")[2] # ICCAT uses 3 colors; IOTC 4 (incl. orange) 
Biplot= TRUE # Produces a "post-modern" biplot with buffer and target zones (Quinn & Collie 2005)
SP.plot = c("standard","phase")[2] # Produces standard or 'Kobe phase' SP plot  
save.trajectories =TRUE # saves posteriors of P=B/K, B/Bmsy and H/Hmsy as .RData object 
harvest.label = c("Hmsy","Fmsy")[2] # choose label preference H/Hmsy versus Fmsy
CPUE.plot= FALSE # Runs state-tool to produce "alligned" multi-CPUE plot  
meanCPUE = FALSE # Uses averaged CPUE from state-space tool instead of individual indices  
Projection = FALSE # Use Projections: requires to define TACs vectors 
save.projections = FALSE # saves projection posteriors as .RData object 
catch.metric = "(t)" # Define catch input metric e.g. (tons) or "(000 t)" 
# Save entire posterior as .RData object
save.all = FALSE # (if TRUE, a very large R object of entire posterior is saved)  
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>


#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Optional: Note Scenarios
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# S3: Base-case Model with time blocks on ESP and JPN 
# Restrospective Analysis 6 years back
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# backward steps
nback = 6 
Scenarios = c("Ref",paste0("RS",1:nback)) 
B.rs = SP.rs = SP.max = Bx = MSY.rs = Bmsy.rs = Fmsy.rs = BtoBmsy.rs = FtoFmsy.rs =  NULL


# Execute multiple JABBA runs in loop 
for(s in 1:(nback+1)){
  Scenario = Scenarios[s] 
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  # Suplus Production model specifications
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  
  # Choose model type: 
  # 1: Schaefer
  # 2: Fox  
  # 3: Pella-Tomlinsson  
  
  Model = rep(3,nback+1)[s] 
  Mod.names = c("Schaefer","Fox","Pella")[Model]
  
  # Depensation opiton:
  # Set Plim = Blim/K where recruitment may become impaired (e.g. Plim = 0.25) 
  # Choose Plim = 0 to reduce to conventional Schaefer, Fox, Pella models 
  Plim = 0
  
  # Required specification for Pella-Tomlinson (Model = 3)
  BmsyK = 0.4 # Set Surplus Production curve inflection point
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>  
  #------------------------------------------------
  
  #--------------------------------------------------
  # Read csv files
  #--------------------------------------------------
  
  # Use SEs from csv file for abudance indices (TRUE/FALSE)
  SE.I = TRUE
  
  # Load assessment data
  catch = read.csv(paste0(File,"/",assessment,"/catch",assessment,".csv"))
  cpue = read.csv(paste0(File,"/",assessment,"/cpue",assessment,".csv"))#
  
  if(SE.I ==TRUE){
    se =  read.csv(paste0(File,"/",assessment,"/se",assessment,".csv"))
  }
  
  names(cpue)
  names(catch)
  
  
  #--------------------------------------------------
  # option to exclude CPUE time series or catch year
  #--------------------------------------------------
  
  # Set up Base-Case for SWO
  
  # Exclude Brazil I
  cpue = cpue[,-c(2)]
  se = se[,-c(2)]
  # Check
  names(cpue)
  ncol(catch)
  ncol(cpue)
  
  Ncol = ncol(cpue)
  Nrow = nrow(cpue)
  
  if(s>1){
  cpue[(Nrow+2-s):Nrow,2:Ncol] = matrix(NA,nrow=s-1,ncol=Ncol-1)  
  #se[(Nrow+2-s):Nrow,2:Ncol] = matrix(NA,nrow=s-1,ncol=Ncol-1)  
  }
   
  #------------------------------------------------------
  # Option use mean CPUE from state-space cpue averaging
  #-----------------------------------------------------
  meanCPUE = FALSE
  
  #------------------------------------------------
  # Prior for unfished biomass K
  #------------------------------------------------
  # The option are: 
  # a) Specify as a lognormal prior with mean and CV 
  # b) Specify as range to be converted into lognormal prior
  
  K.dist = c("lnorm","range")[1]
  
  # if lnorm use mean and CV; if range use lower,upper bound
  K.prior = c(200000,1) 
  
  #-----------------------------------------------------------
  # mean and CV and sd for Initial depletion level P1= SB/SB0
  #-----------------------------------------------------------
  # Set the initial depletion prior B1/K 
  # To be converted into a lognormal prior (with upper bound at 1.1)
  
  psi.dist= c("lnorm","beta")[1]
  # specify as mean and CV 
  psi.prior = c(1,0.25) 
  
  #--------------------------------------------------------------
  # Determine estimation for catchability q and observation error 
  #--------------------------------------------------------------
  # Assign q to CPUE
  sets.q = 1:(ncol(cpue)-1) 
  
  
  #----------------------------------------------------
  # Determine r prior
  #----------------------------------------------------
  # The option are: 
  # a) Specifying a lognormal prior 
  # b) Specifying a resiliance category after Froese et al. (2016; CMSY)
  # Resilience: "Very low", "Low", "Medium", High" (requires r.range = TRUE)
  
  # use range or mean/stdev
  r.dist = c("lnorm","range")[1] # Set to 1 for mean/stdev specifications
  
  r.prior = c(0.42,0.37) # as range with upper and lower bound of lnorm prior
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  # Observation Error
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  
  #To Estimate additional observation variance set sigma.add = TRUE
  sigma.est = TRUE
  
  # Series
  sets.var = 1:(ncol(cpue)-1) # estimate individual additional variace
  
  # As option for data-weighing
  # minimum fixed observation error for each variance set (optional choose 1 value for both)
  fixed.obsE = c(0.2) # Important if SE.I is not availble
  
  # Total observation error: TOE = sqrt(SE^2+sigma.est^2+fixed.obsE^2)
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  # Process Error
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  #Estimate set sigma.proc == True
  sigma.proc = TRUE
  # Determines if process error deviation are estimated for all years (TRUE)  
  # or only from the point the first abundance index becomes available (FALSE)
  proc.dev.all = FALSE 
  #------------------------------------------
  if(sigma.proc == TRUE){
    igamma = c(4,0.01) #specify inv-gamma parameters
    
    # Process error check
    gamma.check = 1/rgamma(1000,igamma[1],igamma[2])
    # check mean process error + CV
    mu.proc = sqrt(mean(gamma.check)); CV.proc = sd(sqrt(gamma.check))/mean(sqrt(gamma.check))
    
    # check CV
    round(c(mu.proc,CV.proc),3)
    quantile(sqrt(gamma.check),c(0.1,0.9))
  }else{
    sigma.proc = 0.07 #IF Fixed: typicallly 0.05-0.15 (see Ono et al. 2012)
  }
  #--------------------------------------------
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  # Optional: Do TAC Projections
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  Projection = FALSE # Switch on by Projection = TRUE 
  
  # Check final year catch 
  catch[nrow(catch),]
  
  # Set range for alternative TAC projections
  TACs = seq(10000,18000,1000) #example
  
  # Intermitted TAC to get to current year
  #TACint = mean(catch[nrow(catch)-3,2]:catch[nrow(catch),2]) # avg last 3 years
  TACint = 10058 # Catch for 2016
  # Set number of projections years
  pyrs = 10
  
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  # Execute model and produce output
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  
  # MCMC settings
  ni <- 30000 # Number of iterations
  nt <- 5 # Steps saved
  nb <- 5000 # Burn-in
  nc <- 2 # number of chains
  nsaved = (ni-nb)/nt*nc
  
   # Run model (BSPSPexe file must be in the same working directory)
  source(paste0(JABBA.file,"/JABBA",version,".R")) 
  
  # Note results
  Bx = cbind(Bx,Bit)
  SP.rs = cbind(SP.rs,SP)
  B.rs = cbind(B.rs,apply(posteriors$SB,2,median))
  MSY.rs = c(MSY.rs,median(posteriors$MSY))
  Bmsy.rs = c(Bmsy.rs,median(posteriors$SBmsy))
  Fmsy.rs = c(Fmsy.rs,median(posteriors$Hmsy))
  BtoBmsy.rs = cbind(BtoBmsy.rs,apply(posteriors$BtoBmsy,2,median))
  FtoFmsy.rs = cbind(FtoFmsy.rs,apply(posteriors$HtoHmsy,2,median))
  
  
  }# THE END
  
#><>><>><>><>><>><>><>><>
# Retrospect analysis
#><>><>><>><>><>><>><>><>

dir.create(paste0(File,"/",assessment,"/Retrospectives"),showWarnings = F)

cols = seq(0.1,0.7,0.7/nback)


#><> plot B
y = B.rs

# Retrospect analysis
Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(File,"/",assessment,"/Retrospectives/RetroB_",assessment,".png"), width = 5, height = 3.5, 
    res = 200, units = "in")
par(Par)

plot(years,years,type="n",ylim=c(0,max(y)),ylab="Biomass (t)",xlab="Year")
lines(years,y[,1],lwd=3,col=4)

for(i in 1:nback){
  lines(years[1:(n.years-i)],y[1:(n.years-i),i],col=grey(cols[i],1),lwd=1,lty=5)
  points(years[1:(n.years-i)],y[1:(n.years-i),i],col=grey(cols[i],1),pch=14+i,cex=0.5)
}

legend("topright",c("Reference",rev(paste(years[(n.years-nback):(n.years-1)]))),col=c(4,grey(cols,1)),bty="n",cex=0.8,pt.cex=0.7,pch=c(-1,seq(15,14+nback,1)),lwd=c(2,rep(1,nback)))
#abline(h=1,lty=2)
#mtext(paste("Year"), side=1, outer=T, at=0.5,line=1,cex=0.9)
dev.off()

#><> Plot BtoBmsy
y = BtoBmsy.rs


# Retrospect analysis
Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(File,"/",assessment,"/Retrospectives/RetroBBmsy_",assessment,".png"), width = 5, height = 3.5, 
    res = 200, units = "in")
par(Par)

plot(years,years,type="n",ylim=c(0,max(y)),ylab=expression(paste(B/B[MSY])),xlab="Year")
lines(years,y[,1],lwd=3,col=4)

for(i in 1:nback){
  lines(years[1:(n.years-i)],y[1:(n.years-i),i+1],col=grey(cols[i],1),lwd=1,lty=5)
  points(years[1:(n.years-i)],y[1:(n.years-i),i+1],col=grey(cols[i],1),pch=14+i,cex=0.5)
}

legend("topright",c("Reference",rev(paste(years[(n.years-nback):(n.years-1)]))),col=c(4,grey(cols,1)),bty="n",cex=0.8,pt.cex=0.7,pch=c(-1,seq(15,14+nback,1)),lwd=c(2,rep(1,nback)))
#abline(h=1,lty=2)
#mtext(paste("Year"), side=1, outer=T, at=0.5,line=1,cex=0.9)
dev.off()

#><> Plot FtoFmsy
y = FtoFmsy.rs


# Retrospect analysis
Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(File,"/",assessment,"/Retrospectives/RetroFFmsy_",assessment,".png"), width = 5, height = 3.5, 
    res = 200, units = "in")
par(Par)

plot(years,years,type="n",ylim=c(0,max(y)),ylab=expression(paste(F/F[MSY])),xlab="Year")
lines(years,y[,1],lwd=3,col=4)

for(i in 1:nback){
  lines(years[1:(n.years-i)],y[1:(n.years-i),i+1],col=grey(cols[i],1),lwd=1,lty=5)
  points(years[1:(n.years-i)],y[1:(n.years-i),i+1],col=grey(cols[i],1),pch=14+i,cex=0.5)
}

#legend("topright",c("Reference",rev(paste(years[(n.years-nback):(n.years-1)]))),col=c(4,grey(cols,1)),bty="n",cex=0.8,pt.cex=0.7,pch=c(-1,seq(15,14+nback,1)),lwd=c(2,rep(1,nback)))
#abline(h=1,lty=2)
#mtext(paste("Year"), side=1, outer=T, at=0.5,line=1,cex=0.9)
dev.off()





# ><> PLOT Surplusproduction
Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(File,"/",assessment,"/Retrospectives/RetroSP_",assessment,".png"), width = 5, height = 4, 
    res = 200, units = "in")
par(Par)

plot(years,years,type="n",ylim=c(0,max(SP.rs)),xlim=c(0,max(Bx)),ylab="Surplus Production (t)",xlab="Biomass (t)")

for(i in 1:nback){
  lines(Bx[,i+1],SP.rs[,i+1],col=grey(cols[i],1),lwd=1,lty=5)
  points(Bx[SP.rs[,i+1]==max(SP.rs[,i+1]),i+1],max(SP.rs[,i+1]),col=grey(cols[i],1),pch=14+i,cex=1.2)
  
}
lines(Bx[,1],SP.rs[,1],col=4,lwd=2)
points(Bx[SP.rs[,1]==max(SP.rs[,1]),1],max(SP.rs[,1]),col=4,cex=1.5,pch=16)
#lines(rep(Bx[SP.rs[,1]==max(SP.rs[,1]),1],2),c(0,max(SP.rs[,1])),lty=2)
#abline(h = max(SP.rs[,1]),lty=2)

legend("topright",c("Reference",rev(paste(years[(n.years-nback):(n.years-1)]))),col=c(4,grey(cols,1)),bty="n",cex=0.8,pt.cex=0.7,pch=c(-1,seq(15,14+nback,1)),lwd=c(2,rep(1,nback)))
#abline(h=1,lty=2)
#mtext(paste("Year"), side=1, outer=T, at=0.5,line=1,cex=0.9)
dev.off()

#Create matrix
# col: 
FRP.rho = c("B", "MSY", "Bmsy", "Fmsy", "BtoBmsy", "FtoFmsy")  
# row: 
yr.rho = c(2014:2009)
nyrs = length(years)
# First get numerator
rho.mat = mat.or.vec(nback,length(FRP.rho))
for(i in 1:nback) rho.mat[i,1] = (B.rs[nyrs-i,i+1]-B.rs[nyrs-i,1])/B.rs[nyrs-i,1] 
rho.mat[,2] = (MSY.rs[-1]-MSY.rs[1])/MSY.rs[1]
rho.mat[,3] = (Bmsy.rs[-1]-Bmsy.rs[1])/Bmsy.rs[1]
rho.mat[,4] = (Fmsy.rs[-1]-Fmsy.rs[1])/Fmsy.rs[1]
for(i in 1:nback) rho.mat[i,5] = (BtoBmsy.rs[nyrs-i,i+1]-BtoBmsy.rs[nyrs-i,1])/BtoBmsy.rs[nyrs-i,1] 
for(i in 1:nback) rho.mat[i,6] = (FtoFmsy.rs[nyrs-i,i+1]-FtoFmsy.rs[nyrs-i,1])/FtoFmsy.rs[nyrs-i,1] 

rho.mat= rbind(rho.mat,apply(rho.mat,2,mean))
rho.mat = round(rho.mat,3)

colnames(rho.mat) = paste(FRP.rho)
rownames(rho.mat) = paste(c(yr.rho,"Mean"))
write.csv(rho.mat,paste0(File,"/",assessment,"/Retrospectives/RHOmat_",assessment,".csv"))
