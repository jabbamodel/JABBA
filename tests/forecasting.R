

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

jbpar(mfrow=c(3,2))
jbplot_ensemble(list(jabba1,jabba2))




prjtac = jabba_fw(jabba1,values=seq(60,80,5)*1000,type="abs",status.quo = 80000,imp.yr = 2,stochastic = T)

prj1 = jabba_fw(jabba1,quant="F",type="ratio",imp.yr = 2,stochastic = T)


prj2 = jabba_fw(jabba2,quant="F",type="ratio",imp.yr = 2,stochastic = T)
joint = jabba_fw(list(jabba1,jabba2),quant="F",type="ratio",imp.yr = 2,stochastic = T)

jbpar(mfrow=c(3,2))
jbplot_ensemble(prj1,add=)

jbpar(mfrow=c(3,2))
jbplot_ensemble(list(jabba1,jabba2,jabba3))

jbplot_summary(list(jabba1,jabba3),cols =ss3col(3))

jbpar(mfrow=c(3,2))
jbplot_ensemble(prjtac)

jabba = list(jabba1,jabba2)

#-----------------------------------------
# AR1 test for ICCAT BUM 2018 (Blue Marlin)
#-----------------------------------------

bum = iccat$bum
jbinput = build_jabba(catch=bum$catch,cpue=bum$cpue,se=bum$se,
                     assessment="BUM",scenario = "basecase",
                     model.type = "Pella",
                     BmsyK = 0.36,
                     K.prior = c(50000,2), 
                     r.prior = c(0.098,0.18),
                     psi.prior = c(1,0.25),
                     fixed.obsE = c(0.25),
                     proc.dev.all = FALSE,
                     sigma.proc = 0.07,
                     verbose = F,rt=T,rho.prior = c(0.3,0.5))  

fit.bum = fit_jabba(jb.bum,quickmcmc = T)

jbplot_summary(fit.bum)
jbplot_procdev(fit.bum)


jbinput <- build_jabba(catch=iccat$bet$catch,cpue=iccat$bet$cpue,se=iccat$bet$se
                       ,scenario = "Fox",model.type="Fox",verbose=FALSE,
                       rt=TRUE)


