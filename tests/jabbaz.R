
library(JABBA)
data("pol.sim")

# compute zage from catch-at-age
set.seed(1234)
z = zage(catch.n,ages=4:12)
z$qname = "Z"
# add some noise
z$data = z$data*rlnorm(nrow(z),0,0.1)

# Prepare df
df = rbind(catch,index,z)
  
inp = reshape2::dcast(df,year~qname,value.var="data",fun.aggregate=sum)
# JABBA data.frames
Index = inp[,c("year","Index")]
Index[Index==0] = NA
Index$Index[1:15] = NA  # exclude observations
Index$Index = Index$Index/mean(Index$Index,na.rm=T)
Catch = inp[,c("year","Catch")]
Z= inp[,c("year","Z")]
Z[Z==0] = NA
Z[1:15,-1] = NA # exclude observations

# Fit with Catch + Index: Simple Fox with r = Fmsy
jbI= build_jabba(catch=Catch,cpue=Index,
                  model.type = "Fox",scenario="Index",
                  r.prior = c(om$fmsy,0.3),
                  verbose=F,
                  psi.prior=c(0.7,0.3))

fI= fit_jabba(jbI,quickmcmc = T,verbose=F)


# Fit with Catch + Index + Effort = qZ 
jbIZ= build_jabba(catch=Catch,cpue=Index,auxiliary = Z,
                  model.type = "Fox",scenario="Index+Z",
                  r.prior = c(om$fmsy,0.3),
                  verbose=F,
                  auxiliary.lag = 4, # lag effect between impact and Z pop structure
                  auxiliary.type = "z", 
                  psi.prior=c(0.7,0.3))


fIZ= fit_jabba(jbIZ,quickmcmc = T,verbose=F)
jbplot_cpuefits(fIZ)
jbplot_residuals(fIZ)

jbpar(mfrow=c(2,2))
jbplot_ensemble(list(fI,fIZ),subplots = 1,add=T)
lines(data~year,om$stock,lwd=2,lty=2)
jbplot_ensemble(list(fI,fIZ),subplots = 2,add=T)
lines(data~year,om$harvest,lwd=2,lty=2)
jbplot_ensemble(list(fI,fIZ),subplots = 5,add=T)
jbplot_ensemble(list(fI,fIZ),subplots = 6,add=T)

## Project 15 years ahead
prj1 = fw_jabba(fIZ,quant="F",type="ratio",imp.yr =2,stochastic = T,imp.values=1,nyears = 15,AR1=T )

jbpar(mfrow=c(2,2))
jbplot_ensemble(prj1,subplots = 1,add=T)
lines(data~year,om$stock,lwd=2,lty=2)
jbplot_ensemble(prj1,subplots = 2,add=T)
lines(data~year,om$harvest,lwd=2,lty=2)
jbplot_ensemble(prj1,subplots = 5,add=T)
jbplot_ensemble(prj1,subplots = 6,add=T)
lines(data~year,om$catch,lwd=2,lty=2)


