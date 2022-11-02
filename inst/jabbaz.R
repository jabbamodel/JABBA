# 
# 
# #setwd("C:/Work/Research/Laurie/JabbaZ")
# #here()
# #load("C:/Work/Research/Laurie/JabbaZ/om.75.RData",verbose=T)
# library(FLCore)
# library(FLBRP)
# 
# library(JABBA)
# plot(om)
# 
# stk = iter(om,1)
# plot(stk)
# stk = window(stk,start=70,end=115)
# # priors
# fmsy = an(refpts(eq)["msy","harvest"])
# bmsy = an(refpts(eq)["msy","ssb"])
# b0 = an(refpts(eq)["virgin","ssb"])
# msy = an(refpts(eq)["msy","yield"])
# shape = bmsy/b0
# 
# 
# # Always convert Z to rate z =  1-exp(-Z)
# zobs =  function(flq,frange = "missing"){
#   out=NULL
#   CN = flq
#   if(missing(frange)){
#     Amin = which(apply(CN,1,median)==max(apply(CN,1,median))) 
#     Amax = nrow(CN)-1 
#    frange= Amin:Amax
#   }
#    
#   for(t in 1:ncol(CN)){
#     Ct= CN[,t]
#     df =   as.data.frame(log(Ct[frange,])) 
#     Z = -coef(lm(data~age,df))[[2]]
#     df = df[1,]
#     df$data = 1-exp(-Z)
#     df$age = "all"
#     out = rbind(out,df)
#   }
#   return(out)
# }
# 
# years = an(dimnames(stk)$year)
# 
# # Prepare df
# # Important don't use too many age-class because of lag
# # If I choose ages 5:15 the default lag 1 works - for 5:20 I need lag2+
# set.seed(1234)
# df = as.data.frame(FLQuants(Catch=catch(stk),
#     Index= vb(stk)*ar1rlnorm(0,years,sdlog=0.25),
#     Z1=zobs(catch.n(stk),frange=5:15)*ar1rlnorm(0,years,sdlog=0.1),
#     Z2=zobs(catch.n(stk),frange=5:15)*ar1rlnorm(0,years,sdlog=0.2)
#     ))
# 
# om = NULL
# om$stock = as.data.frame(ssb(stkfw)/bmsy)
# om$harvest = as.data.frame(fbar(stkfw)/fmsy)
# om$catch = as.data.frame(catch(stkfw))
# 
# om$fmsy = fmsy
# om$bmsy = bmsy
# om$msy = c(refpts(eq)["msy","yield"])
# 
# catch = Catch
# index = Index
# catch.n = as.data.frame(catch.n(stk))
# catch.n = catch.n[catch.n$year>=85,]
# 
# save(catch,index,catch.n,om,file="data/pol.sim.rdata")
# 
# 
# ggplot(df[df$qname%in%c("Z1","Z2"),],aes(year,data,color=qname))+geom_line()
# 
# 
# inp = reshape2::dcast(df,year~qname,value.var="data",fun.aggregate=sum)
# Index = inp[,c("year","Index")]
# Index[Index==0] = NA
# Index$Index[1:15] = NA  # exclude observations
# Index$Index = Index$Index/mean(Index$Index,na.rm=T)
# Catch = inp[,c("year","Catch")]
# Z= inp[,c("year","Z1","Z2")]
# Z.se = Z
# Z.se[,2:3] = 0.1
# 
# Z[Z==0] = NA
# Z[1:15,-1] = NA # exclude observations
# 
# # Fit with Catch + Index: Simple Fox with r = Fmsy
# jbinput= build_jabba(catch=Catch,cpue=Index,auxiliary=Z,
#                  auxiliary.se = Z.se,     
#                   model.type = "Fox",scenario="Index",
#                   r.prior = c(fmsy,0.2),
#                   auxiliary.type = "z",
#                   auxiliary.lag = 4,
#                   shape.CV = 0.2,
#                   auxiliary.obsE = 0.01,auxiliary.sigma = T,verbose=F,
#                   psi.prior=c(0.7,0.3))
# 
# fI= fit_jabba(jbinput,quickmcmc = T,verbose=F)
# jabba=fI
# jbplot_cpuefits(fI,index=1:2)
# jbplot_runstest(fI,index=1:2)
# jbplot_residuals(fI)
# 
# prj1 = fw_jabba(fI,quant="F",type="ratio",imp.yr =2,stochastic = T,imp.values=1,nyears = 15,AR1=T )
# 
# jbplot_ensemble(prj1)
# 
# # Compare: Black dashed line is TRUE pop dyn from OM
# stkfw = iter(om,1)
# jbpar(mfrow=c(2,2))
# jbplot_ensemble(prj1,subplots = 1,add=T)
# lines(data~year,data=as.data.frame(ssb(stkfw)/bmsy),lwd=2,lty=2)
# jbplot_ensemble(prj1,subplots = 2,add=T)
# lines(data~year,data=as.data.frame(fbar(stkfw)/fmsy),lwd=2,lty=2)
# jbplot_ensemble(prj1,subplots = 5,add=T)
# jbplot_ensemble(prj1,subplots = 6,add=T)
# lines(data~year,data=as.data.frame(catch(stkfw)),lwd=2,lty=2)
# 
# 
# 
# # Fit Catch + Index + Z, testing with lag = 4  
# jbIZ.lag4 = build_jabba(catch=Catch,cpue=Index,auxiliary = Z,
#                         model.type = "Fox",scenario="Ind+Z.lag4",
#                         r.prior = c(fmsy,0.3),
#                         auxiliary.type = "z",
#                         psi.prior=c(0.7,0.3),
#                         auxiliary.obsE = 0.1,
#                         auxiliary.sigma=TRUE,
#                         auxiliary.lag = 4, # Lag4
#                         qA.cv=0.1,
#                         verbose=F)
# 
# 
# # Laurie - how to objectively check lag of F to Z response in stock?
# error = ar1rlnorm(0,years,sdlog=0.)
# Z35 =(zobs(stk,frange=5:35)*error)
# Z20 =(zobs(stk,frange=5:20)*error)
# Z15 =(zobs(stk,frange=5:15)*error)
# Z10 =(zobs(stk,frange=5:10)*error)
# Zom = apply(1-exp(-z(stk)[2:29]),2,median) 
# invI= 1/(vb(stk)*error)
# invI = invI/mean(invI)*mean(Z35)
# plot(invI,Z35,Z10,Zom)
# 
# 
# # Fit Catch + Index + Z, testing with lag = 4  
# jbIZ.lag4 = build_jabba(catch=Catch,cpue=Index,auxiliary = Z,
#                    model.type = "Fox",scenario="Ind+Z.lag4",
#                    r.prior = c(fmsy,0.3),
#                    auxiliary.type = "z",
#                    psi.prior=c(0.7,0.3),
#                    auxiliary.obsE = 0.2,
#                    auxiliary.lag = 4, # Lag4
#                    verbose=F)
# 
# fIZ=  fit_jabba(jbIZ.lag4,quickmcmc = T,verbose=F)
# 
# # Fit only relative H ~ qZ, H is effort proxy with lag 4
# jbIE.lag4 = build_jabba(catch=Catch,cpue=Index,auxiliary = Z,
#                    model.type = "Fox",scenario="Ind+Eff.lag4",
#                    r.prior = c(fmsy,0.3),
#                    auxiliary.type = "effort",
#                    auxiliary.obsE = 0.2,
#                    auxiliary.lag = 4,
#                    verbose=F,
#                    psi.prior=c(0.7,0.3))
# 
# fIE=  fit_jabba(jbIE.lag4,quickmcmc = T,verbose=F)
# 
# 
# # Compare: Black dashed line is TRUE pop dyn from OM
# jbpar(mfrow=c(2,2))
# jbplot_ensemble(list(fI,fIZ,fIE),subplots = 1,add=T)
# lines(data~year,data=as.data.frame(ssb(stk)/bmsy),lwd=2,lty=2)
# jbplot_ensemble(list(fI,fIZ,fIE),subplots = 2,add=T)
# lines(data~year,data=as.data.frame(fbar(stk)/fmsy),lwd=2,lty=2)
# jbplot_ensemble(list(fI,fIZ,fIE),subplots = 5,add=T)
# jbplot_ensemble(list(fI,fIZ,fIE),subplots = 6,add=T)
# 
# # Check fits
# jbpar(mfrow=c(2,2))
# jbplot_cpuefits(fIZ,add=T)
# jbplot_cpuefits(fIE,add=T)
# 
# 
# 
# 
# # Compare Lags
# lags = 1:4
# jbs = lapply(lags,function(x){
#   jb = build_jabba(catch=Catch,cpue=Index,auxiliary = Z,
#                    model.type = "Fox",
#                    scenario=paste0("Ind+Z.lag",x),
#                    r.prior = c(fmsy,0.3),
#                    auxiliary.type = "z",
#                    psi.prior=c(0.7,0.3),
#                    auxiliary.obsE = 0.2,
#                    auxiliary.lag = x, # Lag2
#                    verbose=F)
#   fit_jabba(jb,quickmcmc = T,verbose=F)
# })
# 
# jbpar(mfrow=c(2,2))
# jbplot_ensemble(jbs,subplots = 1,add=T)
# lines(data~year,data=as.data.frame(ssb(stk)/bmsy),lwd=2,lty=2)
# jbplot_ensemble(jbs,subplots = 2,add=T)
# lines(data~year,data=as.data.frame(fbar(stk)/fmsy),lwd=2,lty=2)
# jbplot_ensemble(jbs,subplots = 5,add=T)
# jbplot_ensemble(jbs,subplots = 6,add=T)
# 
# jbpar(mfrow=c(2,2))
# for(i in 1:4){
#   jbplot_runstest(jbs[[i]],index=2,add=T)
#   legend("topleft",paste("Lag",i),bty="n") 
# }
# 
# jbpar(mfrow=c(2,2))
# for(i in 1:4){
#   jbplot_cpuefits(jbs[[i]],index=2,add=T)
#   legend("topleft",paste("Lag",i),bty="n") 
# }
