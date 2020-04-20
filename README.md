# JABBA: Just Another Bayesian Biomass Assessment
The materials in this repository present the stock assessment tool ‘Just Another Bayesian Biomass Assessment’ JABBA. The motivation for developing JABBA was to provide a user-friendly R to JAGS (Plummer) interface for fitting generalized Bayesian State-Space SPMs with the aim to generate reproducible stock status estimates and diagnostics. Building on recent advances in optimizing the fitting procedures through the development of Bayesian state-space modelling approaches, JABBA originates from a continuous development process of a Bayesian State-Space SPM tool that has been applied and tested in many assessments across oceans. JABBA was conceived in the Hawaiian Summer of 2015 as a collaboration between young researchers from South Africa and the Pacific Islands Fisheries Science Center (NOAA) in Honolulu, HI USA. The goal was to provide a bridge between age-structured and biomass dynamic models, which are still widely used. JABBA runs quickly and by default generates many useful plots and diagnosic tools for stock assessments.

Inbuilt JABBA features include:

+ Integrated state-space tool for averaging multiple CPUE series (+SE) for optional use in assessments
+ Automatic fitting of multiple CPUE time series and associated standard errors
+ Fox, Schaefer or Pella Tomlinson production function (optional as input Bmsy/K)
+ Kobe-type biplot plotting functions 
+ Forecasting for alternative TACs 
+ Residual and MCMC diagnostics 
+ Estimating or fixing the process variance
+ Optional estimation additional observation variance for individual or grouped CPUE time series
+ Easy implementation of time-block changes in selectivity
+ New: Estimating Catch with Error 
+ New: Estimating shape with prior
+ New: Catch-Only option with additional relative biomass priors
+ New: Inbuilt retrospective and hindcasting run and plotting options 

## Installing JABBA as R package

`library(devtools)` <br>
`install_github("jabbamodel/JABBA")`

## Test-drive JABBA

`library(JABBA)`
`data(iccat)`<br>

`# Compile JABBA JAGS model and input object for bigeye tuna (bet)`<br>
`jbinput = build_jabba(catch=iccat$bet$catch,cpue=iccat$bet$cpue,se=iccat$bet$se,assessment="BET",scenario = "TestRun",model.type = "Fox",sigma.est = FALSE,fixed.obsE = 0.01)`

`# Fit JABBA (here mostly default value - careful)`<br>
`bet1 = fit_jabba(jbinput)`

`# Make individual plots` <br>
`jbplot_catcherror(bet1)` <br>
`jbplot_ppdist(bet1)`     <br>
`jbplot_cpuefits(bet1)`<br>
`jbplot_logfits(bet1)`<br>

`# Plot Status Summary` <br>
`par(mfrow=c(3,2),mar = c(3.5, 3.5, 0.5, 0.1))` <br>
`jbplot_trj(bet1,type="B",add=T)` <br>
`jbplot_trj(bet1,type="F",add=T)` <br>
`jbplot_trj(bet1,type="BBmsy",add=T)`<br>
`jbplot_trj(bet1,type="FFmsy",add=T)`<br>
`jbplot_spphase(bet1,add=T)`<br>
`jbplot_kobe(bet1,add=T)`<br>
`# Test run end`

More detailed examples of JABBA applications to ICCAT stock can be found in [Examples_iccat.R](https://github.com/jabbamodel/JABBA/blob/master/Example/Examples_iccat.R).

An updated vignette is currently under development...

**Reference**  
[Winker, H., Carvalho, F., Kapur, M. (2018) <U>JABBA: Just Another Bayesian Biomass Assessment.</U> *Fisheries 
Research* **204**: 275-288.](https://www.sciencedirect.com/science/article/pii/S0165783618300845)   


--------------------------------------------------------------------------------

#### [JABBAbeta](https://github.com/Henning-Winker/JABBAbeta) GitHub repository
JABBA development version for testing new JABBA features and stock assessment examples 
