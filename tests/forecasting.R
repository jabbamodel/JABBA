

library(JABBA)
data(iccat)
jbinput0 <- build_jabba(catch=iccat$bet$catch,cpue=iccat$bet$cpue,se=iccat$bet$se,scenario = "Fox.prj",model.type="Fox",verbose=FALSE,
                        projection=TRUE,
                        TACs = seq(60,80,5)*1000,
                        TACint = 80000,
                        imp.yr = 2019
)
jbinput1 <- build_jabba(catch=iccat$bet$catch,cpue=iccat$bet$cpue,se=iccat$bet$se,scenario = "Fox",model.type="Fox",verbose=FALSE)
jbinput2 <- build_jabba(catch=iccat$bet$catch,cpue=iccat$bet$cpue,se=iccat$bet$se,scenario = "Schaefer",model.type="Schaefer",verbose=FALSE)
jbinput3 <- build_jabba(catch=iccat$bet$catch[10:68,],cpue=iccat$bet$cpue[10:68,],se=iccat$bet$se[10:68,],scenario = "Short",model.type="Fox",verbose=FALSE)

jabba0 = fit_jabba(jbinput0,quickmcmc=TRUE,verbose=FALSE)
jabba1 = fit_jabba(jbinput1,quickmcmc=TRUE,verbose=FALSE)
jabba2 = fit_jabba(jbinput2,quickmcmc=TRUE,verbose=FALSE)
jabba3 = fit_jabba(jbinput3,quickmcmc=TRUE,verbose=FALSE)


prjtac = jabba_fw(jabba1,values=seq(60,80,5)*1000,type="abs",status.quo = 80000,imp.yr = 2,stochastic = T)

prj1 = jabba_fw(jabba1,quant="F",type="ratio",imp.yr = 2,stochastic = T)
prj2 = jabba_fw(jabba2,quant="F",type="ratio",imp.yr = 2,stochastic = T)
joint = jabba_fw(list(jabba1,jabba2),quant="F",type="ratio",imp.yr = 2,stochastic = T)

jbpar(mfrow=c(3,2))
jbplot_ensemble(prj1,add=)

jbpar(mfrow=c(3,2))
jbplot_ensemble(list(jabba1,jabba3))

jbplot_summary(list(jabba1,jabba3))

