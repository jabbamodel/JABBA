# JABBA
JABBA (Just Another Bayesian Biomass Assessment) is a generalized Bayesian State-Space Surplus Production Model and represents a further development of the modeling framework applied in the 2016 ICCAT South Atlantic blue shark assessment (Carvalho and Winker 2016) and 2017 North Pacific blue shark assessment, the 2017 ICCAT Mediterranean Albacore assessment and the 2017 ICCAT North and South Atlantic shortfin mako assessments. JABBA was conceived in the Summer of 2015 as a collaboration between the South Africa Department of Agriculture, Forestry and Fisheries and the Pacific Islands Fisheries Science Center (NOAA) in Honolulu, HI USA. The goal was to provide a bridge between age-structured and biomass dynamic models, which are still widely used in data-limited fisheries and in developing countries on a routine basis. JABBA runs quickly and by default generates many useful plots and diagnosic tools for stock assessments.

Inbuilt JABBA features include:

+ Integrated state-space tool for averaging multiple CPUE series (+SE) for optional use in assessments (Fig. 1)
+ Automatic fitting of multiple CPUE time series and associated standard errors (Fig.2)
+ Fox, Schaefer or Pella Tomlinson production function (optional as input Bmsy/K)
+ Kobe-type biplot plotting functions 
+ Forecasting for alternative TACs 
+ Residual (Fig.5) and MCMC diagnostics 
+ Estimating or fixing the process variance
+ Optional estimation additional observation variance for individual or grouped CPUE time series
+ Easy implementation of time-block changes in selectivity

The materials in this repository present the stock assessment tool ‘Just Another Bayesian Biomass Assessment’ JABBA. The motivation for developing JABBA was to provide a user-friendly R to JAGS (Plummer) interface for fitting generalized Bayesian State-Space SPMs with the aim to generate reproducible stock status estimates and diagnostics. Building on recent advances in optimizing the fitting procedures through the development of Bayesian state-space modelling approaches (Meyer and Millar 1999, Millar and Meyer 2000, Thorson et al. 2014), JABBA originates from a continuous development process of a Bayesian State-Space SPM tool that has been applied and tested in assessments of South Atlantic blue shark (Ref 2015), North Pacific blue shark (Ref 2016), Mediterranean albacore tuna (Ref 2017) and North and South Atlantic shortfin mako shark (Ref 2017).

Reference
Winker, H., Carvalho, F., Kapur, M. JABBA: Just Another Bayesian Biomass Assessment. (In review at Fisheries 
Research). 

<B>A self-contained R package of JABBA is forthcoming.</b>
