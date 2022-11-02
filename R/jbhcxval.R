#' jbhcxval()
#' 
#' additional hindcast options with external foreward projections 
#'
#' @param hindcasts object (list of models) from hindcast_jabba() 
#' @param stochastic if FALSE, process error sigma.proc is set to zero 
#' @param AR1 if TRUE, projection account auto correlation in the process devs 
#' @param rho if AR1 = TRUE, the autocorrelation coefficient is estimated from the proc devs
#' @param sigma.proc option to specify the process error other than the posterior estimate
#' @param ndevs number years on the tail to set initial proc.error for forecasting  
#' @param run option to assign a scenario name other than specified in build_jabba()
#' @return data.frame of kobe posterior model + forecast scenarios
#' @export
#' @importFrom utils tail
#' @examples
#' data(iccat)
#' whm <- iccat$whm
#' 
#' # ICCAT white marlin setup
#' jb <- build_jabba(catch=whm$catch, 
#'                   cpue=whm$cpue,
#'                   se=whm$se,
#'                   assessment="WHM",
#'                   scenario = "BaseCase",
#'                   model.type = "Pella",
#'                   r.prior = c(0.181,0.18),
#'                   BmsyK = 0.39,
#'                   igamma = c(0.001,0.001))
#' fit <- fit_jabba(jb,quickmcmc=TRUE,verbose=TRUE)
#' 
#' # Hindcast
#' hc <- hindcast_jabba(jbinput=jb,fit=fit,peels=1:5)
#' jbplot_retro(hc)
#' jbplot_hcxval(hc,index=c(8,11))
#' hc.ar1 <- jbhcxval(hc,AR1=TRUE) # do hindcasting with AR1
#' jbplot_hcxval(hc.ar1,index=c(8,11))

jbhcxval <- function(hindcasts,
                     stochastic = c(TRUE, FALSE)[1],
                     AR1 = c(TRUE, FALSE)[1],
                     sigma.proc = NULL,
                     rho = NULL,
                     ndevs=1,
                     run = NULL){

peels = do.call(c,lapply(hindcasts,function(x){
  x$diags$retro.peels[1]
}))
peels = as.numeric(peels[peels>0])

# Cut internal forecasts
hc = lapply(hindcasts[-1],function(x){
og = x
nyears=length(tail(x$yr,x$diags$retro.peels[1]))
x$yr =  x$yr[x$yr%in%tail(x$yr,x$diags$retro.peels[1])==FALSE]
x$kbtrj = x$kbtrj[x$kbtrj$year%in%x$yr,]
fwtrj = fw_jabba(x,nyears=x$diags$retro.peels[1],
                  imp.yr = NULL,
                  quant = "Catch",
                  initial = x$catch[-c(1:length(x$yr))],
                  imp.values = 1,
                  type="abs",
                  stochastic = stochastic ,
                  AR1 = AR1,
                  sigma.proc = sigma.proc,
                  rho = rho,ndevs=ndevs)

fc =fwtrj[fwtrj$run ==unique(fwtrj$run)[2],]
fc = fc[fc$year!=min(fc$year),]
og$kbtrj =rbind(fwtrj[fwtrj$run ==unique(fwtrj$run)[1],],fc)
og
})

# Update forecast time-series
x= hc[[1]]
y = as.list(peels)[[1]]

hcts = Map(function(x,y){
ts = x$timeseries
ny = tail(1:length(x$yr),y)

for(i in 1:length(ny)){
posteriors = x$kbtrj[x$kbtrj$year==x$yr[ny[i]],]  
x$timeseries[ny[i],,1:6] =  cbind(quantile(posteriors$B,c(0.5,0.025,0.975)),
                       quantile(posteriors$H,c(0.5,0.025,0.975)),
                       quantile(posteriors$stock,c(0.5,0.025,0.975)),
                       quantile(posteriors$harvest,c(0.5,0.025,0.975)),
                       quantile(posteriors$BB0,c(0.5,0.025,0.975)),
                       quantile(posteriors$Bdev,c(0.5,0.025,0.975)))
}                       
x
},x=hc,y=as.list(peels))                       

# Update forecast CPUE
hcI = lapply(hcts, function(x){
qs = as.matrix(x$pars_posterior[,grep("q",names(x$pars_posterior))])
sets.q = x$settings$sets.q
nq = length(sets.q)
diags=x$diags
idxs = unique(diags$name)
for(j in 1:nq){
  if(tail(diags[diags$name==idxs[j],]$hindcast,1)==TRUE){
    sub = diags[diags$name==idxs[j]&diags$hindcast,]
    nhc = nrow(sub) 
    for(i in 1:nhc){
    hat = c(quantile(x$kbtrj[x$kbtrj$year==sub$year[i],]$B*qs[,sets.q[j]],
          c(0.5,0.025,0.975)))         
    diags[diags$name==idxs[j]&diags$hindcast&diags$year==sub$year[i],
          c("hat","hat.lci","hat.uci")] = hat
    } # end i
    }} # end j
   x$diags = diags
   x
})
out = c(hindcasts[1],hcI)
return(out)
}
#}}}
