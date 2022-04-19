#' hindcast_jabba()
#'
#' Wrapper to fit retrospectives 
#' @param jbinput object from build_jabba()
#' @param jbfit fitted model from fit_jabba
#' MCMC settings
#' @param ni number of iterations
#' @param nt thinning interval of saved iterations
#' @param nb burn-in
#' @param nc number of mcmc chains
#' @param quickmcmc Reduces MCMC iters for hindcasting reference run
#' Initial values
#' @param init.values if TRUE init values from fit are used 
#' @param peels sequence of retrospective peels default 1:5
#' @param verbose if FALSE run silent
#' @return hc containing estimates of key joint results from all hindcast run 
#' @export
#' @examples
#' data(iccat)
#' whm = iccat$whm
#' # ICCAT white marlin setup
#' jb = build_jabba(catch=whm$catch,cpue=whm$cpue,se=whm$se,assessment="WHM",scenario = "BaseCase",model.type = "Pella",r.prior = c(0.181,0.18),BmsyK = 0.39,igamma = c(0.001,0.001))
#' fit = fit_jabba(jb,quickmcmc=TRUE,verbose=TRUE)
#' hc = hindcast_jabba(jbinput=jb,fit=fit,peels=1:5)
#' jbplot_retro(hc)
#' jbplot_hcxval(hc,index=c(8,11))
#' hc.ar1 = jbhcxval(hc,AR1=TRUE) # do hindcasting with AR1
#' jbplot_hcxval(hc.ar1,index=c(8,11))

hindcast_jabba = function(jbinput,fit,
                          # MCMC settings
                          ni = NULL, # Number of iterations
                          nt = NULL, # Steps saved
                          nb = NULL, # Burn-in
                          nc = NULL, # number of chains
                          quickmcmc = TRUE,
                          # init values
                          init.values = TRUE,
                          peels = 1:5, # retro peel option
                          verbose=FALSE){
 
  
  runs = as.list(peels)
  if(is.null(ni) & !quickmcmc) ni = fit$settings$mcmc$ni
  if(is.null(nt) & !quickmcmc) ni = fit$settings$mcmc$nt
  if(is.null(nb) & !quickmcmc) ni = fit$settings$mcmc$nb
  if(is.null(nc) & !quickmcmc) ni = fit$settings$mcmc$nc
  
    if(init.values){
      Kin = fit$pars[1,1]
      rin = fit$pars[2,1]
      qin = fit$pars[(3:(2+fit$settings$nq)),1]
    } else {
      Kin = NULL
      rin = NULL
      qin = NULL
    }
    
    fithc = lapply(runs,function(x){
                       fit_jabba(jbinput,save.trj = TRUE,
                      ni = ni, # Number of iterations
                      nt = nt, # Steps saved
                      nb = nb, # Burn-in
                      nc = nc, # number of chains
                      init.values = init.values,
                      init.K = Kin,
                      init.r = rin,
                      init.q = qin,
                      quickmcmc = T,# vector
                      peels = x,verbose=verbose) # retro peel option
    })
    
    
    retro = c(list(fit),fithc)
    
    names(retro) = c(fit$scenario,paste0("-",max(fit$yr)-peels+1))
    # assign scenario names 
    retro = Map(function(x,y){
      x$settings$scenario = y
      x
    },x=retro,y=as.list(names(retro)))
    
    
    return(retro)
} #end of retro function
