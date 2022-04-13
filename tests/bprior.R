

library(JABBA)

#><>><>><>><>><>><>><>><>><>><>><>
# Bigeye Tuna ICCAT
#><>><>><>><>><>><>><>><>><>><>><>

assessment = "BETiccat"
#------------------------------------------------------
# COMs with priors for b/k, b/bmsy and f/fmsy
#-------------------------------------------------------
data(iccat)
# get BET data
bet = iccat$bet
# Compile JABBA JAGS model and input object
jbinput = build_jabba(catch=bet$catch,cpue=bet$cpue,se=bet$se,assessment="bet",scenario = "cpue",model.type = "Fox")
# Fit JABBA (here mostly default value - careful)
fit = fit_jabba(jbinput,quickmcmc = TRUE) # quick run

# COM b/k prior in 2015
combk = build_jabba(catch=bet$catch,assessment="bet",scenario = "bk",model.type = "Fox",sigma.est = FALSE,fixed.obsE = 0.01,
                    b.prior = c(0.3,0.2,2015,"bk"))
fitbk = fit_jabba(combk,quickmcmc = TRUE)
# COM b/bmsy prior in 2015
combbmsy = build_jabba(catch=bet$catch,assessment="bet",scenario = "bbmsy",model.type = "Fox",sigma.est = FALSE,fixed.obsE = 0.01,
                    b.prior = c(0.7,0.2,2015,"bbmsy"))
fitbbmsy = fit_jabba(combbmsy,quickmcmc = TRUE)

# COM f/fmsy prior in 2015
comffmsy = build_jabba(catch=bet$catch,assessment="bet",scenario = "ffmsy",model.type = "Fox",sigma.est = FALSE,fixed.obsE = 0.01,
                       b.prior = c(1.5,0.2,2015,"ffmsy"))
fitffmsy = fit_jabba(comffmsy,quickmcmc = TRUE)

# Compare
jbplot_summary(list(fit,fitbk,fitbbmsy,fitffmsy))

# plot priors
jbpar(mfrow=c(1,3),plot.cex=0.8,labs=T)
jbplot_bprior(fitbk,add=T)
jbplot_bprior(fitbbmsy,add=T)
jbplot_bprior(fitffmsy,add=T)


