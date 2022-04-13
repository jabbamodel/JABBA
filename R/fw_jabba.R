#{{{
#' jabba_fw()
#
#' External forward projections in JABBA 
#'
#' @param jabba objects from fit_jabba() or list of fit_jabba() objects  
#' @param quant quantity to forecast  c("Catch","F")
#' \itemize{
#'   \item Catch  
#'   \item F 
#' }    
#' @param values vector Catch or F values provide as absolute or ratios   
#' @param type choose the type of values provided  c("ratio","msy","abs")[1]
#' \itemize{
#'   \item ratio (relative to status quo Catch or F)   
#'   \item msy (relative to Fmsy or MSY)
#'   \item abs (input values are taken as absolute values)
#' }    
#' @param stochastic if FALSE, process error sigma.proc is set to zero 
#' @param AR1 if TRUE, projection account auto correlation in the process devs 
#' @param rho if AR1 = TRUE, the autocorrelation coefficient is estimated from the proc devs
#' @param sigma.proc option to specify the process error other than the posterior estimate
#' @param run option to assign a scenario name other than specified in build_jabba()
#' @param thin option to thin the posterior at rates > 1
#' @return data.frame of kobe posterior model + forecast scenarios
#' @export

#{{{
jabba_fw <- function(jabba,nyears = 10, imp.yr = 2,
                     status.quo = NULL,
                     values = seq(0.6,1.2,0.1),
                     quant = c("Catch","F")[1],
                     type = c("ratio","msy","abs")[1],
                     stochastic = c(TRUE, FALSE)[1],
                     AR1 = c(TRUE, FALSE)[2],
                     sigma.proc = NULL,
                     rho = NULL,
                     run = NULL,
                     thin=1){
  
  #TODO add hcr options (e.g. ICES style)
  if(!is.null(jabba$settings)){
    jabba=list(jabba)
    if(is.null(run)) run = jabba[[1]]$scenario
  }
  
  
  if(is.null(jabba[[1]]$kbtrj)) stop("Use option fit_jabba(...,save.trj=TRUE) to enable forecasting")
  if(is.null(run)) run = "Joint"
  thinning = seq(1,length(jabba[[1]]$pars_posterior$K),thin)
  year= jabba[[1]]$yr
  yrend = max( year)
  n = length(year)
  pars = do.call(rbind,lapply(jabba,function(x){
    y = x$pars_posterior[thinning,]
    y
  }))
  
  rho.est = mean(do.call(cbind,lapply(jabba,function(x){
    y =   as.numeric(x$timeseries[,"mu","procB"])[which(apply(x$settings$I,1,sum,na.rm=T)>0)[1]:length(year)]
    stats::cor(y [-length(y )],y[-1])
  })))
  
  trj = do.call(rbind,lapply(jabba,function(x){
    y = x$kbtrj
    y = y[y$iter%in%thinning,]
    y$run = run
    y
  }))
  
  if(quant=="Catch"){
    if(is.null(status.quo)){
      status.quo= mean(do.call(c,lapply(jabba,function(x){
        mean(x$catch[(n-2):n])})))
    }}
  if(quant=="F"){
    if(is.null(status.quo)){
      status.quo= mean(do.call(c,lapply(jabba,function(x){
        mean(x$timeseries[(n-2):n,1,2])})))
    }}
  
  
  
  inits = data.frame(trj[trj$year==yrend,])
  
  
  # Get quants
  k = pars[["K"]]
  r = pars[["r"]]
  m = pars[["m"]]
  if(is.null(sigma.proc)){
    sigma = sqrt(pars$sigma2)} else {
      sigma = sigma.proc
    }
  iters = length(k) 
  
  if(AR1){
    if(is.null(rho)) rho = rho.est
  }
  if(!AR1){
    rho=0  
  }
  if(stochastic){
    if(is.null(sigma.proc)) sigma=sqrt(pars$sigma2)
  }
  if(!stochastic){
    sigma = 0 
  }
  
  # get management quants
  fmsy = r/(m-1)*(1-1/m)
  shape = (m)^(-1/(m-1))
  bmsy= shape*k
  msy= bmsy*fmsy 
  pyears = c(year[n],year[n]+(1:nyears))  
  iters = length(k) 
  P = devs = B = H =  matrix(NA,ncol = nyears+1, nrow = iters)
  C = matrix(NA,ncol = nyears+1, nrow = iters)
  C[,1] =inits$Catch
  devs[,1]= inits$Bdev
  for(i in 2:ncol(P)){
    if(rho>0) devs[,i] = rho * devs[,i - 1] + sqrt(1 - rho^2) * rnorm(iters,0,sigma)
    if(rho==0) devs[,i] = rnorm(iters,0,sigma)
    
  }
  
  P[,1] = inits$BB0
  B[,1] = inits$B
  H[,1] = inits$H
  
  
  if(quant=="Catch"){
    C[,2:(imp.yr+1)][] = status.quo
  }
  if(quant=="F"){
    H[,2:(imp.yr+1)][] = status.quo
  }
  
  
  # Create list of forecast values
  if(type=="ratio"){
    vals = values*status.quo
    runs = paste0(round(100*values,0),"%")
  }
  if(type=="msy"){
    if(quant=="F") vals = values*median(fmsy)
    if(quant=="Catch") vals = values*median(msy)
    runs = paste0(round(100*values,0),"%")
  }
  if(type=="abs"){
    vals =  values 
    if(quant=="F") runs = paste0("F",vals)
    if(quant=="Catch"){
      runs = paste0("C",vals)
      if(max(vals)>=100000) runs = paste0("C",round(vals/1000,1))
    }
  } 
  kb = NULL
  
  fw.ls = split(data.frame(runs,vals),seq(length(runs)))
  names(fw.ls) = runs
  
  
  
  kobe = do.call(rbind,lapply(fw.ls,function(x){ 
    
    if(quant=="Catch")C[,(imp.yr+1):ncol(C)][] = x[[2]]
    if(quant=="F") H[,(imp.yr+1):ncol(C)][] = x[[2]]
    
    kb = data.frame(year=pyears[1],run=x[[1]],type="fit",iter=1:iters,
                    stock=B[,1]/bmsy,harvest=H[,1]/fmsy,B=B[,1],H=H[,1],
                    Bdev=devs[,1],Catch=C[,1],BB0=P[,1])
    
    for(i in 2:ncol(P)){
      P[,i] <- pmax((P[,i-1]+  r/(m-1)*P[,i-1]*(1-P[,i-1]^(m-1)) - C[,i-1]/k),0.005)*exp(devs[,i])
      B[,i] = P[,i]*k
      if(quant=="Catch"){
        H[,i] = C[,i]/B[,i]#,fmax*median(fmsy)) # fmax constraint
      } 
      if(quant=="F") C[,i] = H[,i]*B[,i]
      kb = rbind(kb,data.frame(year=pyears[i],run=x[[1]],type="prj",iter=1:iters,
                               stock=B[,i]/bmsy,harvest=H[,i]/fmsy,B=B[,i],H=H[,i],
                               Bdev=devs[,i],Catch=C[,i],BB0=P[,i]))
    }
    kb
  }))
  rownames(kobe) = 1:nrow(kobe)
  kb = rbind(trj,kobe)
  
  kb$run = factor(kb$run,levels=unique(kb$run)) 
  
  
  
  return(kb)
  
} 
# }}} End of function