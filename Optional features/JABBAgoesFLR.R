
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# JABBA kobe projections
# Developed by Henning Winker & Laurie Kell (Madrid, 2017)
# Builds on FLR Kobe R package 
# https://github.com/flr/kobe/blob/master/R/kobe-jabba.R
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

#***********************************************
# This R code produces standardized output for *
# ICCAT scientific management advice from      *
# saved .RData Projections and Trajecties      *
#***********************************************

# Requires FLR Kobe package by Laurie Kell 
library(FLCore)
library(ggplotFL)
library(kobe)
require(plyr)
require(dplyr)
require(reshape2)
require(scales)

# Run for all Scenarios
for(s in 3:length(Scenarios)){
Scenario = Scenarios[s]

#-----------------------------------------------------------------------------
# repeat basic settings as in Prime
Model = c(3,3,3)[s] 
Mod.names = c("Schaefer","Fox","Pella")[Model]
#------------------------------------------------------------------------------

# sets up directories automatically
output.dir = paste0(File,"/",assessment,"/",Scenario,"_",Mod.names,"/Output")
dir.create(paste0(output.dir,"/FLRout"),showWarnings = FALSE)


# Load RData from runs
load(paste0(output.dir,"/",assessment,"_",Mod.names,"_",Scenario,"_projections"))
load(paste0(output.dir,"/",assessment,"_",Mod.names,"_",Scenario,"_trajectories"))

#><>><>><>><>><>><>><>><>><>><>><>
# FRL KOB package: Make Kobe plot
kbPrj=transform(kjp,tac=TACs[tac]) # Prepare projections object
#><>><>><>><>><>><>><>><>><>><>><>
# Make plot
Par = list(mfrow=c(1,1),mar = c(5, 5, 1, 1), mgp =c(3,1,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/FLRout/KobeFLR_",assessment,"_",Scenario,".png"), width = 6.5, height = 5.5, 
    res = 200, units = "in")
par(Par)
kobe:::kobePhaseMar3(subset(kb,year==2015),col = colorRampPalette(c("darkgrey"))
                     ,xlab = expression(B/B[MSY]),ylab = expression(F/F[MSY]))

dev.off()

#><>><>><>><>><>><>
# Projection Plot
#><>><>><>><>><>><>
stk=FLQuants(dlply(kbPrj,.(tac), with, as.FLQuant(data.frame(year=year,iter=iter,data=stock))))

# Plot
Par = list(mfrow=c(1,1),mar = c(5, 5, 1, 1), mgp =c(3,1,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/FLRout/ProjFLR_",assessment,"_",Scenario,".png"), width = 6, height = 4, 
    res = 200, units = "in")
par(Par)
p=plot(stk)$data
print(ggplot(p)+
         geom_hline(aes(yintercept=1))+
         geom_line(aes(year,`50%`,col=qname))+
         xlab("Year")+ylab(expression(B/B[MSY]))+ 
         labs(colour = "TAC")+coord_cartesian(ylim = c(0, 2)) +theme_classic())
dev.off()

  
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Produce Projection matrices for ICCAT advice
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

# Prepare matrices using library(kob)
   t.=ddply(kbPrj,.(year,tac), with, kobe:::smry(stock,harvest))
  
  # Compile projection matrices in tables 
  k2smTab=list()
  k2smTab[[1]]=cast(t.,tac~year,value="underFishing")
  k2smTab[[2]]=cast(t.,tac~year,value="underFished")
  k2smTab[[3]]=cast(t.,tac~year,value="green")
  
  
  nTAC = length(TACs)
  for(j in 1:3){
  # Define Projection matrix
  ypr = (as.numeric(format(Sys.Date(), "%Y"))+1):max(proj.yrs)
  npr = length(ypr)
  dy = min(ypr)-proj.yrs[1] 
  pjm = round(k2smTab[[j]][,-c(1:(dy+1))]*100,0)
  
  
  mat.names = c("/pjmF_","/pjmB_","/pjmKobe_")
  # Write table
  pjm.save = k2smTab[[j]][,-c(2:(dy+1))]
  pjm.save[,2:ncol(pjm.save)] = round(pjm.save[,2:ncol(pjm.save)]*100,1)
  write.csv(pjm.save,paste0(output.dir,paste(mat.names[j]),assessment,"_",Scenario,".csv"))
  
  
  op <- par(mar = rep(0, 4),mai = rep(0, 4),omi = rep(0, 4))
  png(file = paste0(output.dir,"/FLRout/",paste(mat.names[j]),assessment,"_",Scenario,".png"), width = 6.5, height = 6.5*nTAC/npr*0.4 , 
      res = 200, units = "in")
  par(op)
  plot(1,1,axes=FALSE,frame.plot=FALSE,xlab="",ylab="",type="n",ylim=c(0,nTAC+1),xlim=c(-1,npr+1))
  # first line
  rect(-1,1:nTAC+1,1,0:nTAC);  
  text(rep(0,nTAC+1),seq(0.5,nTAC+0.5,1),c(rev(paste(TACs)),"TAC | Year"),cex=0.8)
  
  # Set grey shading
  mat= as.matrix(k2smTab[[j]])[,-c(1:(dy))]
  mat[mat<0.5]=-1
  mat=(1-mat)*2
  mat[mat>1]=1
  
  rect(1:npr,rep(nTAC+1,npr),1:(npr+1),rep(nTAC,npr))
  text(c(seq(1.5,npr+0.5,1)),rep(nTAC+0.5,npr),c(paste(ypr)),cex=0.8)
  for(i in 1:nTAC) rect(1:npr,rep(nTAC+1-i,npr),2:(npr+1),rep(nTAC-i,npr),col=grey(mat[i,],0.5))
  for(i in 1:nTAC) text(c(seq(1.5,npr+0.5,1)),rep(nTAC+0.5-i,npr),paste(pjm[i,]),cex=0.8)
  dev.off()
  }
  
}# End of Scenario loop
 
# FLR done
 