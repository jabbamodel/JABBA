
# Installation
# install.packages(devtools)
# devtools:: install_github("jabbamodel/JABBA")

library(JABBA)
File = "C:/Work/Research/GitHub/JABBApkg_testing/example" # LINK to your folder of choice here

#><>><>><>><>><>><>><>><>><>><>><>
# Bigeye Tuna ICCAT
#><>><>><>><>><>><>><>><>><>><>><>

assessment = "BETiccat"
output.dir = file.path(File,assessment)
dir.create(output.dir,showWarnings = F)
setwd(output.dir)
#------------------------------------------------------
# Simple example fit JABBA to Catch and CPUE with SEs
#-------------------------------------------------------
data(iccat)
# get BET data
bet = iccat$bet
# Compile JABBA JAGS model and input object
jbinput = build_jabba(catch=bet$catch,cpue=bet$cpue,se=bet$se,assessment=assessment,scenario = "Test",model.type = "Fox",sigma.est = FALSE,fixed.obsE = 0.01,igamma = c(0.001,0.001))
# Fit JABBA (here mostly default value - careful)
bet1 = fit_jabba(jbinput,save.jabba=TRUE,output.dir=output.dir,quickmcmc = TRUE) # quick run

# Make individual plots
jbplot_catch(bet1)
jbplot_catcherror(bet1)
jbplot_ppdist(bet1)
jbplot_mcmc(bet1)
jbplot_residuals(bet1)
jbplot_cpuefits(bet1)
jbplot_runstest(bet1)
jbplot_logfits(bet1)
jbplot_procdev(bet1)

# Try to improve runs test diagnostics by changing the variance settings
jbinput = build_jabba(catch=bet$catch,cpue=bet$cpue,se=NULL,assessment=assessment,scenario = "Ref",model.type = "Fox",sigma.est = TRUE,fixed.obsE = 0.1,igamma = c(0.001,0.001),psi.prior = c(1,0.1))
bet2 = fit_jabba(jbinput,save.jabba=TRUE,output.dir=output.dir)
# Check residual diags
jbplot_cpuefits(bet2)
jbplot_runstest(bet2)
jbplot_logfits(bet2)
# Improved
refinput = jbinput # Note as reference input 

# status summary
par(mfrow=c(3,2),mar = c(3.5, 3.5, 0.5, 0.1))
jbplot_trj(bet2,type="B",add=T)
jbplot_trj(bet2,type="F",add=T)
jbplot_trj(bet2,type="BBmsy",add=T)
jbplot_trj(bet2,type="FFmsy",add=T)
jbplot_spphase(bet2,add=T)
jbplot_kobe(bet2,add=T)

# Write all as png
jabba_plots(jabba=bet2,output.dir = output.dir)


#------------------------------------------------------
# Estimate shape m as function of Bmsy/K
#-------------------------------------------------------

# Compile JABBA JAGS model and input object
jbinput = build_jabba(catch=bet$catch,cpue=bet$cpue,se=NULL,assessment=assessment,scenario = "Est.Shape",model.type = "Pella_m",BmsyK=0.37,sigma.est = TRUE,fixed.obsE = 0.1,igamma = c(0.001,0.001),psi.prior = c(1,0.1))
bet3 = fit_jabba(jbinput,save.jabba=TRUE,output.dir=output.dir)

jbplot_ppdist(bet3) # check shape prior & posterior dist
# Compare
par(mfrow=c(2,2))
jbplot_trj(bet2,type="BBmsy",add=T)
jbplot_trj(bet3,type="BBmsy",add=T)
jbplot_kobe(bet2,add=T)
jbplot_kobe(bet3,add=T)



#------------------------------------------------------
# Catch-Only with biomass prior in 2010 as type B/Bmsy
#------------------------------------------------------
# Compile JABBA JAGS model and input object for Catch Only
# Add biomass prior based on B/Bmsy guestimate
jbinput = build_jabba(catch=bet$catch,model.type = "Fox",assessment=assessment,scenario =  "CatchOnly" ,b.prior=c(0.7,0.2,2010,"bbmsy"),psi.prior = c(1,0.1))
# Fit JABBA
bet4 = fit_jabba(jbinput,save.jabba=TRUE,output.dir=output.dir)

# Check depletion prior vs posterior
jbplot_bprior(bet4)
# Compare
par(mfrow=c(3,2))
jbplot_trj(bet2,type="BBmsy",add=T)
jbplot_trj(bet2,type="FFmsy",add=T)
jbplot_trj(bet3,type="BBmsy",add=T)
jbplot_trj(bet3,type="FFmsy",add=T)
jbplot_trj(bet4,type="BBmsy",add=T)
jbplot_trj(bet4,type="FFmsy",add=T)


#-------------------------------------------
# Make summary plot comparing the three scenarios
#-------------------------------------------
Scenarios = (c("Ref","Est_shape","CatchOnly")) # Scenarios to be loaded as Rdata objects
#  Check plot with CIs
jbplot_summary(assessment=assessment,scenarios = Scenarios,mod.path = output.dir)
# and without CIs
jbplot_summary(assessment=assessment,scenarios = Scenarios,plotCIs=FALSE)
# Check Base only
jbplot_summary(assessment=assessment,scenarios = Scenarios[1],prefix="SmryBase",as.png = F)
# Save comparison 
jbplot_summary(assessment=assessment,scenarios = Scenarios,prefix="Comp3runs",save.summary = T,as.png = T,output.dir = output.dir)


#----------------------------------------------------------------
# Conduct Retrospective Analysis and Hind-Cast Cross-Validation
#----------------------------------------------------------------
# Organize folders by creating a "retro" subfolder
retro.dir = file.path(output.dir,"retro")
dir.create(retro.dir,showWarnings = F)

# Run hindcasts for reference model (set plotall = TRUE if you want to save plots from all runs)
hc = jabba_hindcast(refinput,save.hc=T,plotall=T,output.dir = retro.dir,peels = 0:7)

# Retro Analysis Summary plot
jbplot_retro(hc,as.png = F,single.plots = F,output.dir = retro.dir)
# Save plot and note Mohn's rho statistic
mohnsrho = jbplot_retro(hc,as.png = T,single.plots = F,output.dir = retro.dir)
# Zoom-in
mohnsrho = jbplot_retro(hc,as.png = F,single.plots = F,output.dir = retro.dir,Xlim=c(2000,2014))
# eval mohnsrho
mohnsrho

# Do Hindcast Cross-Validation (hcxval) 
# show multiplot
jbplot_hcxval(hc,single.plots = F,as.png = F,output.dir=retro.dir)
# Zoom-in
jbplot_hcxval(hc,single.plots = F,as.png = F,output.dir=retro.dir,minyr=2000)
# save as png and note summary stats 
mase = jbplot_hcxval(hc,single.plots = F,as.png = TRUE,output.dir=retro.dir)
#check stats
mase

#---------------------------
# Run REF with projections
#---------------------------
jbinput = build_jabba(catch=bet$catch,cpue=bet$cpue,se=NULL,assessment=assessment,scenario = "Ref",
                      model.type = "Fox",sigma.est = TRUE,fixed.obsE = 0.1,igamma = c(0.001,0.001),psi.prior = c(1,0.1),
                      projection=TRUE,TACs=seq(45000,90000,5000))


betprj = fit_jabba(jbinput,output.dir=output.dir,save.csvs = T)
# plot with CIs (80% for projections)
jbplot_prj(betprj,type="BBmsy")
jbplot_prj(betprj,type="BB0")

# or without CIs (80% for projections) 
jbplot_prj(betprj,type="FFmsy", CIs=FALSE)


#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Original South Atlantic Swordfish example here
# Winker et al. (2018). JABBA: Just Another Biomass Assessment
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
swos = iccat$swos

assessment = "SWOSiccat"
scenario = "Base"
output.dir = file.path(File,assessment)
dir.create(output.dir,showWarnings = F)
setwd(output.dir)

# Compile JABBA JAGS model and input object
jbinput = build_jabba(catch=swos$catch,cpue=swos$cpue,se=swos$se,assessment=assessment,scenario = scenario,
                      model.type = "Pella",
                      BmsyK = 0.4,
                      r.prior=c(0.42,0.37),
                      K.prior = c(250000,1),
                      psi.prior = c(1,0.25),
                      fixed.obsE = 0.2,
                      add.catch.CV = FALSE,
                      proc.dev.all = FALSE, 
                      igamma=c(4,0.01),
                      P_bound = c(0.02,1.1))
                      
# fit JABBA
swos = fit_jabba(jbinput,save.jabba=TRUE,output.dir=output.dir)

# Plot all
jabba_plots(jabba=swos,output.dir = output.dir)

# Organize folders by creating a "retro" subfolder
retro.dir = file.path(output.dir,"retro")
dir.create(retro.dir,showWarnings = F)

# Run hindcasts
hc = jabba_hindcast(jbinput,save.hc=T,plotall=F,output.dir = retro.dir,peels = 0:7)

# Retro Analysis Summary plot
jbplot_retro(hc,as.png = F,single.plots = F,output.dir = retro.dir)
# Zoom-in
mohnsrho = jbplot_retro(hc,as.png = F,single.plots = F,output.dir = retro.dir,Xlim=c(1990,2015))

# Save plot and note Mohn's rho statistic
mohnsrho = jbplot_retro(hc,as.png = T,single.plots = F,output.dir = retro.dir)
# eval mohnsrho
mohnsrho

mohnsrho["rho.mu",]
