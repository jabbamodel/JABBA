library(JABBA)
whm = iccat$whm
# Simple default (Schaefer)
jb1 = build_jabba(catch=whm$catch,cpue=whm$cpue,se=NULL,assessment="WHM",
                       scenario = "Simple")
fit1 = fit_jabba(jb1,quickmcmc = T)
# Base case (not 100%)
jb2 = build_jabba(catch=whm$catch,cpue=whm$cpue,se=whm$se,assessment="WHM",
                  scenario = "BaseCase",model.type = "Pella",
                  r.prior = c(0.181,0.18),
                  BmsyK = 0.39,
                  igamma = c(0.001,0.001))

fit2 = fit_jabba(jb2,quickmcmc = T)

# Compare
fits = list(fit1,fit2)
jbplot_summary(fits)
jbplot_ensemble(fit1,joint=T,kbout=T)

# IOTC type
fw.io = fw_jabba(fit2,nyears=10,imp.yr=2,imp.values = seq(0.6,1.2,0.1),quant="Catch",type="ratio",nsq=3,stochastic = T)

jbplot_ensemble(fw.io)
jbpar(mfrow=c(2,2))
jbplot_ensemble(fw.io,add=T,subplots = 1)
jbplot_ensemble(fw.io,add=T,subplots = 2)
jbplot_ensemble(fw.io,add=T,subplots = 5)
jbplot_ensemble(fw.io,add=T,subplots = 6)

fw.ioar1 = fw_jabba(fit2,nyears=10,imp.yr=2,imp.values = seq(0.6,1.2,0.1),quant="Catch",type="ratio",nsq=3,stochastic = T,AR1=T)
jbpar(mfrow=c(2,2))
jbplot_ensemble(fw.ioar1,add=T,subplots = 1)
jbplot_ensemble(fw.ioar1,add=T,subplots = 2)
jbplot_ensemble(fw.ioar1,add=T,subplots = 5)
jbplot_ensemble(fw.ioar1,add=T,subplots = 6)

fw.ioar1.det = fw_jabba(fit2,nyears=10,imp.yr=2,imp.values = seq(0.6,1.2,0.1),quant="Catch",type="ratio",nsq=3,stochastic = F,AR1=T)
jbpar(mfrow=c(2,2))
jbplot_ensemble(fw.ioar1.det ,add=T,subplots = 1)
jbplot_ensemble(fw.ioar1.det ,add=T,subplots = 2)
jbplot_ensemble(fw.ioar1.det ,add=T,subplots = 5)
jbplot_ensemble(fw.ioar1.det ,add=T,subplots = 6)


# Ensemble 
fw.ioar1.det = fw_jabba(fits,nyears=10,imp.yr=2,imp.values = seq(0.6,1.2,0.1),quant="Catch",type="ratio",nsq=3,stochastic = F,AR1=T)
jbpar(mfrow=c(2,2))
jbplot_ensemble(fw.ioar1.det ,add=T,subplots = 1)
jbplot_ensemble(fw.ioar1.det ,add=T,subplots = 2)
jbplot_ensemble(fw.ioar1.det ,add=T,subplots = 5)
jbplot_ensemble(fw.ioar1.det ,add=T,subplots = 6)

# ICCAT Style
TACs = seq(300,800,50)

fw.iccat= fw_jabba(fits,nyears=10,imp.yr=2,initial = 475,imp.values = TACs,quant="Catch",type="abs",nsq=3,stochastic = F,AR1=T)

jbplot_ensemble(fw.iccat,legendcex = 0.4,xlim=c(2010,2027))





# Single Forecast for Base-Case model
fw1 = fw_jabba(fit2,nyears=10,imp.yr=2,imp.values = seq(0.7,1.2,0.1),quant="F",type="msy",stochastic = T)
jbpar(mfrow=c(3,2))
jbplot_ensemble(fw1)
# Zoom-in
jbplot_ensemble(fw1,xlim=c(2015,2027))

# Forecast with AR1 process error
fw1.ar1 = fw_jabba(fits[[2]],nyears=10,imp.yr=2,quant="F",type="msy",AR1=TRUE,stochastic = T)
# now compare
jbpar(mfrow=c(3,2))
for(i in 1:3){
jbplot_ensemble(fw1,subplots = i,add=T,xlim=c(2010,2028))
jbplot_ensemble(fw1.ar1,subplots = i,add=T,xlim=c(2010,2028))
}

# Do hindcast cross-validation
hc1 = hindcast_jabba(jb2,fit2,peels=1:5) 
hc2 = jbhcxval(hc1,AR1=T) # make forecasts with AR1
hc3 = jbhcxval(hc1,AR1=T,ndevs=3) # mean devs based on last 3 years
hc0 = jbhcxval(hc1,AR1=F) # check hc1 by redoing forecasts externally 




jbpar(mfrow=c(2,3))
for(i in 1:2){
jbplot_hcxval(hc1,index=c(8,11)[i],add=T,minyr = 2000,legend.add = F)
jbplot_hcxval(hc2,index=c(8,11)[i],add=T,minyr = 2000,legend.add = F)
jbplot_hcxval(hc3,index=c(8,11)[i],add=T,minyr = 2000,legend.add = F)
}

jbpar(mfrow=c(3,2))
jbplot_ensemble(list(random=hc1[[6]],ar1=hc3[[6]]))

jbplot_hcxval(hc1,index=11,add=T,minyr = 2000,legend.add = F)
hc=hc1


# 2-model ensemble forecast
fw2 = fw_jabba(fits,nyears=10,imp.yr=3,quant="F",type="msy")
jbplot_ensemble(fw2,xlim=c(1980,2028))
fw2.ar1 = fw_jabba(fits,nyears=10,imp.yr=3,quant="F",type="msy",AR1=TRUE)
jbplot_ensemble(fw2.ar1,xlim=c(1980,2028))
fw2.ar1.n5 = fw_jabba(fits[[1]],nyears=10,imp.yr=3,quant="F",type="msy",AR1=TRUE,ndevs=5)
jbplot_ensemble(fw2.ar1.n5,xlim=c(1980,2028))


