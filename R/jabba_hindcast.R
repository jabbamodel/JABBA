#' jabba hindcasting function
#'
#' Wrapper to coduct histcasts for retrospective analysis and cross-validation
#' @param jbinput List of input variables as output by build_jabba()
#' MCMC settings
#' @param ni number of iterations
#' @param nt thinning interval of saved iterations
#' @param nb burn-in
#' @param nc number of mcmc chains
#' Initial values
#' @param init.values = FALSE,
#' @param init.K = NULL,
#' @param init.r = NULL,
#' @param init.q = NULL,# vector
#' @param peels sequence of retro spective peels default 0:5
#' @param save.jabba Save individual JABBA run outputs
#' @param output.dir path to save plot. default is getwd()
#' @param save.hc Save hindcast list output as .rdata
#' @param plotall if TRUE makes jabba_plots() for each run    
#' @param speedup Reduces MCMC after setting runs 2+ inits to first "full" reference run
#' @return hc containing estimates of key joint results from all hindcast run 
#' @export
jabba_hindcast = function(jbinput,
                          # MCMC settings
                          ni = 30000, # Number of iterations
                          nt = 5, # Steps saved
                          nb = 5000, # Burn-in
                          nc = 2, # number of chains
                          # init values
                          init.values = FALSE,
                          init.K = NULL,
                          init.r = NULL,
                          init.q = NULL,# vector
                          peels = 0:5, # retro peel option
                          save.jabba = FALSE,
                          output.dir = getwd(),
                          save.hc = FALSE,
                          plotall = FALSE,
                          speedup = TRUE){
  # hindcast object define object  
  hc = list(scenario = jbinput$settings$scenario, yr=jbinput$data$yr,catch=jbinput$jagsdata$TC,peels=peels,timeseries = NULL,refpts=NULL,pfunc=NULL,diags=NULL,settings=NULL)
  hc$settings$cols = jbinput$settings$cols
  hc$settings$harvest = jbinput$settings$harvest.label
  hc$settings$catch.metric = jbinput$settings$catch.metric
  
  Scenario = jbinput$settings$scenario
  for(i in 1:length(peels)){
    jbinput$settings$scenario = peels[i]
    if(i == 1 | speedup==FALSE){
      mci = ni
      mct = nt
      mcb = nb
      mcc = nc
      Kin = init.K
      rin = init.r
      qin = init.q
    } else {
      mci = 11000
      mct = 2
      mcb = 1000
      mcc = 2
      init.values=TRUE
    } 
    if(peels[i]%in%peels[2]){
      
      Kin = fithc$pars[1,1]
      rin = fithc$pars[2,1]
      qin = fithc$pars[(3:(2+fithc$settings$nq)),1]
    }
    
    fithc = fit_jabba(jbinput,save.jabba=save.jabba,output.dir=output.dir,
                      ni = mci, # Number of iterations
                      nt = mct, # Steps saved
                      nb = mcb, # Burn-in
                      nc = mcc, # number of chains
                      init.values = init.values,
                      init.K = Kin,
                      init.r = rin,
                      init.q = qin,# vector
                      peels = peels[i]) # retro peel option
    
    hc$timeseries$mu = rbind(hc$timeseries$mu,data.frame(factor=fithc$diags[1,1],level=fithc$diags[1,2],fithc$timeseries[,"mu",])) 
    hc$timeseries$lci = rbind(hc$timeseries$lci,data.frame(factor=fithc$diags[1,1],level=fithc$diags[1,2],fithc$timeseries[,"lci",])) 
    hc$timeseries$uci = rbind(hc$timeseries$uci,data.frame(factor=fithc$diags[1,1],level=fithc$diags[1,2],fithc$timeseries[,"uci",])) 
    hc$diags = rbind(hc$diags,fithc$diags)
    hc$refpts= rbind(hc$refpts,fithc$refpts[1,])
    hc$pfunc= rbind(hc$pfunc ,fithc$pfunc)
    
    if(plotall==TRUE){
      jabba_plots(fithc,output.dir = output.dir)
    }
    
  } # end of loop
  if(save.hc==TRUE){
    save(hc,file=paste0(output.dir,"/hc_",Scenario,".rdata"))  
  }
  return(hc)
} #end of hindcast function 
