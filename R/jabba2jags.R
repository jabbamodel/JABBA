#' Function to write JABBA model in JAGS
#'
#' Writes JAGS code for JABBA into a temporary directory
#' @param jbinput JABBA input object from build_jabba()
#' @export

jabba2jags = function(jbinput, dir){
  
  # JAGS MODEL Standard
  sink(file.path(dir,"JABBA.jags"))
  cat("

    model {

    # Prior specifications
    eps <- 0.0000000000000000000000000000000001 # small constant


    # Process variance prior
    isigma2.est ~ dgamma(igamma[1],igamma[2])

    # Carrying Capacity SB0
    K ~ dlnorm(log(K.pr[1]),pow(K.pr[2], -2))

    # informative priors for Hmsy as a function of r
    r ~ dlnorm(log(r.pr[1]),pow(r.pr[2],-2))

    ")
  
 
  if(jbinput$settings$q1.dist=="range"){
    cat("
    
    q[1] ~ dunif(q1.pr[1],q1.pr[2])
    
    ",append=TRUE)} else {
      cat("
    
    q[1] ~ dlnorm(log(q1.pr[1]),pow(q1.pr[2],-2))
    
    ",append=TRUE)
    }

   
  
  if(jbinput$jagsdata$nq>1){
    cat("
    
  #Catchability coefficients
  for(i in 2:nq)
  {
    q[i] ~ dunif(q_bounds[1],q_bounds[2])
  }
  
  ",append=TRUE)
  } else {
    cat("
    
    nqq <- nq
    bq1 <- q_bounds[1]
    bq2 <- q_bounds[2]
    
    ",append=TRUE)
  }
    
 
  
 
   if(jbinput$settings$model.id==4){
    cat("
      # Shape m prior
      m ~ dlnorm(log(mu.m),pow(m.CV,-2))
      ",append=TRUE)
  }else{ cat("
      m <- mu.m
    ",append=TRUE)}
  
  if(jbinput$settings$psi.dist =="beta"){
    cat("
      # Beta Prior for Biomass depletion at the start (deteministic)
      psi ~ dbeta(psi.pr[1],psi.pr[2])
      ",append=TRUE)
  } else {
    cat("
      # Lognormal for Biomass depletion at the start (deteministic)
      psi ~ dlnorm(log(psi.pr[1]),pow(psi.pr[2],-2)) #I(0.1,1.1)
      ",append=TRUE)
  }
  
  if(jbinput$settings$sigma.proc==TRUE){
    cat("
      # Process variance
      isigma2 <- isigma2.est
      sigma2 <- pow(isigma2,-1)
      sigma <- sqrt(sigma2)
      fakesigma.fixed <- sigma.fixed # Prevent unused variable error msg
      ",append=TRUE)
  }else{ cat("
      isigma2 <- pow(sigma.fixed+eps,-2)
           sigma2 <- pow(isigma2,-1)
           sigma <- sqrt(sigma2)

           ",append=TRUE)}
  
  if(jbinput$settings$sigma.est==TRUE){
    cat("
      # Observation variance
      for(i in 1:nvar)
      {
      # Observation error
      itau2[i]~ dgamma(0.001,0.001)
      tau2[i] <- 1/itau2[i]
      }

      for(i in 1:nI)
      {
      for(t in 1:N)
      {
      var.obs[t,i] <- SE2[t,i]+tau2[sets.var[i]]
      ivar.obs[t,i] <- 1/var.obs[t,i]
      # note total observation error (TOE)
      TOE[t,i] <- sqrt(var.obs[t,i]) # Total observation variance

      }}
      ",append=TRUE)
  }else{ cat("
      # Observation variance
           for(i in 1:nvar)
           {
           # Observation error
           itau2[i]~ dgamma(4,0.01)
           tau2[i] <- 1/itau2[i]
           }

           for(i in 1:nI)
           {
           for(t in 1:N)
           {
           var.obs[t,i] <- SE2[t,i] # drop tau2
           fake.tau[t,i] <- tau2[sets.var[i]]

           ivar.obs[t,i] <- 1/var.obs[t,i]
           # note total observation error (TOE)
           TOE[t,i] <- sqrt(var.obs[t,i])

           }}

           ",append=TRUE)}
  
  # Run standard JABBA
  if(jbinput$settings$add.catch.CV==FALSE){
    cat("
  for(t in 1:N){
  estC[t] <- TC[t]
  }
",append=TRUE)} else if(jbinput$settings$catch.error=="random"){
  cat("
  for(t in 1:N){
      estC[t] ~ dlnorm(log(TC[t]),pow(CV.C[t],-2))
  }
  

",append=TRUE)} else if (jbinput$settings$catch.error=="under"){
  cat("
  for(t in 1:N){
  cdev[t] ~ dnorm(0,pow(CV.C[t],-2))
  cpos[t] <- sqrt(pow(cdev[t],2))
  estC[t] <- exp(log(TC[t])+cpos[t]) 
  }
",append=TRUE)}  
  
  
  
  
  proc.mode = jbinput$settings$proc.mode
  if(is.null(proc.mode)) proc.mode = ifelse(isTRUE(jbinput$settings$proc.dyn),"dyn","const")

  ebsb.f = function(p, c) paste0("ebsbm[",c,",1] + (ebsbm[",c,",2]-ebsbm[",c,",1])*(1-exp(-ebsbm[",c,",3]*(min(max(",p,",ebsbm[",c,",4]),ebsbm[",c,",5])-ebsbm[",c,",4])))/(1-exp(-ebsbm[",c,",3]*(ebsbm[",c,",5]-ebsbm[",c,",4])))")
  ebsb.sub = function(x){
    if(isTRUE(jbinput$settings$index.ebsb.on)){
      x = gsub("EBSBIDX", "ebsbi[t,i]", x, fixed=TRUE)
      # rho_c(P): EB/SB of curve c; ebsbi: curve of index i; ebsbr: effective EB/SB of the removals (catch-share weighted harmonic mean)
      gsub("EBSBCODE", paste0("for(t in 1:N){ for(c in 1:(ncur+1)){ rhoc[t,c] <- ", ebsb.f("P[t]","c"), "; irhoc[t,c] <- 1/rhoc[t,c] }\n",
        "    for(i in 1:nI){ ebsbi[t,i] <- rhoc[t,ebsbci[i]] }\n",
        "    ebsbr[t] <- 1/inprod(ebsbcw[t,1:ncur], irhoc[t,1:ncur]) }"), x, fixed=TRUE)
    } else gsub("EBSBCODE", "", gsub("*EBSBIDX", "", x, fixed=TRUE), fixed=TRUE)
  }
  pex.txt = if(isTRUE(jbinput$settings$proc.xfix.on)) "pe.xfix[t]" else "BtoBmsy[t]"
  if(proc.mode=="floor"){
    cat("

    # Recruitment-driven process error floor + estimated additive process error
    #   sd_rec(x) = sig.msy * x^-k,  x = B/Bmsy truncated to pe.xlim
    #   total process variance = sd_rec(x)^2 + sigma2 (sigma2: estimated additive part)
    # iprocdyn rescales isigma2 so that isigma2*iprocdyn = 1/(sigma2 + sd_rec^2)
    ",append=TRUE)
    if(isTRUE(jbinput$settings$proc.k.est)){
      cat("
    # k (power of B/Bmsy) estimated with truncated normal prior; sig.msy fixed
    kpe ~ dnorm(pe.kpr[1], pow(pe.kpr[2],-2)) T(pe.kpr[3], pe.kpr[4])
    ",append=TRUE)
    } else {
      cat("
    kpe <- pe.floor[2]
    ",append=TRUE)
    }
    cat(gsub("PEX1", if(isTRUE(jbinput$settings$proc.xfix.on)) "pe.xfix[1]" else "psi/SBmsy_K", "
    # first-year state gets the same total process variance (depletion from psi, avoids a cycle through P[1])
    iproc1 <- sigma2/(sigma2 + pow(pe.floor[1]*pow(min(max(PEX1, pe.xlim[1]), pe.xlim[2]), -kpe), 2))
    "),append=TRUE)
    cat(gsub("PEX", pex.txt, "
    for(t in 1:N){
      xpe[t]      <- min(max(PEX, pe.xlim[1]), pe.xlim[2])
      sig2rec[t]  <- pow(pe.floor[1]*pow(xpe[t], -kpe), 2)
      iprocdyn[t] <- sigma2/(sigma2 + sig2rec[t])
      procsd[t]   <- sqrt(sigma2 + sig2rec[t])
    }

    "),append=TRUE)
  } else if(proc.mode=="dyn"){
    cat("
    iproc1 <- 1

    # Dynamic process-error scaling from pesim workflow
    # relsig = sigma.proc(B/Bmsy) / sigma.proc(Bmsy)
    # iprocdyn is an inverse-variance multiplier:
    #   iprocdyn = 1 / relsig^2
    #
    # Function:
    #   raw(x) = a_pe + (1 - a_pe) / (1 + (x / c_pe)^d_pe)
    #   relsig(x) = raw(x) / raw(1)
    #
    # Therefore relsig = 1 at B/Bmsy = 1 and iprocdyn = 1 at B/Bmsy = 1.
   
   # Emperical ices data
   #    a     c     d 
   # 0.128 0.223 1.326 
   a_pe <- 0.128
   c_pe <- 0.223
   d_pe <- 1.326 
   

    relsig.ref <- a_pe + (1 - a_pe) / (1 + pow(1 / c_pe, d_pe))

    for(t in 1:(N)){
      relsig.raw[t] <- a_pe + (1 - a_pe) / (1 + pow(BtoBmsy[t] / c_pe, d_pe))
      relsig[t]     <- relsig.raw[t] / relsig.ref
      iprocdyn[t]   <- pow(relsig[t], -2)
      procsd[t]     <- sigma*relsig[t]
    }

    ",append=TRUE)} else {
      cat("

    for(t in 1:(N+1)){
      iprocdyn[t] <- 1
    }
    iproc1 <- 1

    ",append=TRUE)
    }
  
  
  # Low-biomass production correction: expected log deviation lowb[t] = cpe * min(log(B/Bmsy), 0)
  lowb.mode = jbinput$settings$proc.lowb.mode
  if(is.null(lowb.mode)) lowb.mode = "none"
  if(lowb.mode=="none"){
    cat("
    for(t in 1:(N+1)){ lowb[t] <- 0 }
    ",append=TRUE)
  } else {
    if(lowb.mode=="fixed") cat("
    cpe <- pe.lowb[1]
    ",append=TRUE)
    if(lowb.mode=="est") cat("
    cpe ~ dnorm(pe.lowb[1], pow(pe.lowb[2],-2)) T(0, pe.lowb[3])
    ",append=TRUE)
    cat("
    # Realised production at low B falls below the deterministic curve (age structure + recruitment variability)
    for(t in 1:(N+1)){ lowb[t] <- cpe*min(log(max(BtoBmsy[t], pe.lowbx)), 0) }
    ",append=TRUE)
  }
  
if(jbinput$settings$catch.error!="fdev"){
  cat("

    #Process equation
    Pmean[1] <- log(psi)
    iPV[1] <- ifelse(1<(stI),10000,isigma2) # inverse process variance
    P[1] ~ dlnorm(Pmean[1],iPV[1]*ifelse(1<(stI),1,iproc1)) # set to small noise instead of isigma2
    penB[1]  <- ifelse(P[1]<P_bound[1],log(K*P[1])-log(K*P_bound[1]),ifelse(P[1]>P_bound[2],log(K*P[1])-log(K*P_bound[2]),0)) # penalty if Pmean is outside viable biomass
    penBK[1] <- 0

    # Process equation
    # Production is driven by lagged depletion Pl[t] = P[lagidx[t]] (lagidx[t] = max(t-proc.lag,1));
    # proc.lag = 0 gives the standard model (Pl = P[t-1] at t-1)
    for (t in 2:(N+1))
    {
    Pl[t-1] <- P[lagidx[t-1]]
    Pmean[t] <- ifelse(Pl[t-1] > Plim,
    log(max(P[t-1] +  r/(m-1)*Pl[t-1]*(1-pow(Pl[t-1],m-1)) - estC[t-1]/K,0.001)),
    log(max(P[t-1] +  r/(m-1)*Pl[t-1]*(1-pow(Pl[t-1],m-1))*Pl[t-1]*slope.HS - estC[t-1]/K,0.001))) + lowb[t-1]
    iPV[t] <- ifelse(t<(stI),10000,isigma2) # inverse process variance
    P[t] ~ dlnorm(Pmean[t],iPV[t]*ifelse(t<(stI),1,iprocdyn[t-1])) # dynamic process error only from the first index year
    }
    for (t in 2:N)
    {
    penB[t]  <- ifelse(P[t]<(P_bound[1]),log(K*P[t])-log(K*(P_bound[1])),ifelse(P[t]>P_bound[2],log(K*P[t])-log(K*(P_bound[2])),0)) # penalty if Pmean is outside viable biomass
    # Depletion prior
    penBK[t] <- ifelse(b.yr[t] < 1,0,log(ifelse(b.pr[4]<1,P[t],ifelse(b.pr[4]>1,HtoHmsy[t],BtoBmsy[t])))-log(b.pr[1]))
    }


    # Process error deviation
    for(t in 1:N){
    Proc.Dev[t] <- log(P[t]*K)-log(exp(Pmean[t])*K)}

    # Enforce soft penalties on bounds for P
    for(t in 1:N){
    pen.P[t] ~ dnorm(penB[t],1000) # enforce penalty with CV = 0.1
    pen.bk[t] ~ dnorm(penBK[t],pow(b.pr[2],-2))
    }

    Hmsy <- r*pow(m-1,-1)*(1-1/m)
    SBmsy_K <- (m)^(-1/(m-1))
    SBmsy <- SBmsy_K*K
    MSY <- SBmsy*Hmsy
    
    for (t in 1:N)
    {
    
    H[t] <- estC[t]/SB[t]
    HtoHmsy[t] <- H[t]/(Hmsy)
    }


    
    
    for (t in 1:(N+1)) # One step ahead biomass
    {
    SB[t] <- K*P[t]
    BtoBmsy[t] <- SB[t]/SBmsy
    }

    ",append=TRUE)}
 
  if(jbinput$settings$catch.error=="fdev"){ ### NOT ACTIVE
    cat("
    
      frho ~ dnorm(0.5,pow(0.5,-2))T(0.01,0.9)
      Hmean ~ dlnorm(log(0.05),0.0001) 
      #df ~  dnorm(0,0.001)
      #rho ~ dbeta(rho.pr[1],rho.pr[2])
      #rho <- rho1/1.1
      isigf2 ~ dgamma(0.001,0.001)
      sigf2 <- pow(isigf2,-1)
      sigf <- sqrt(sigf2)
      fdev[1] ~ dnorm(0,isigf2)
      for(t in 2:(N)){
      fdev[t] ~ dnorm(frho*fdev[t-1],isigf2) # Add AR1 on r
      }
      for(t in 1:(N)){
      H[t] <- Hmean*exp(fdev[t]-0.5*sigf2)
      }
      
    #Process equation
    Pmean[1] <- log(psi)
    iPV[1] <- ifelse(1<(stI),10000,isigma2) # inverse process variance
    P[1] ~ dlnorm(Pmean[1],iPV[1]*ifelse(1<(stI),1,iproc1)) # set to small noise instead of isigma2
    penB[1]  <- ifelse(P[1]<P_bound[1],log(K*P[1])-log(K*P_bound[1]),ifelse(P[1]>P_bound[2],log(K*P[1])-log(K*P_bound[2]),0)) # penalty if Pmean is outside viable biomass
    penBK[1] <- 0
    
    # Process equation
    for (t in 2:(N+1))
    {
    Pmean[t] <- ifelse(P[t-1] > Plim,
    log(max(P[t-1] +  r/(m-1)*P[t-1]*(1-pow(P[t-1],m-1)) - H[t-1]*P[t-1],0.001)),
    log(max(P[t-1] +  r/(m-1)*P[t-1]*(1-pow(P[t-1],m-1))*P[t-1]*slope.HS - H[t-1]*P[t-1],0.001)))
    iPV[t] <- ifelse(t<(stI),10000,isigma2*iprocdyn[t-1]) # inverse process variance
    P[t] ~ dlnorm(Pmean[t],iPV[t])
    }
    for (t in 2:N)
    {
    penB[t]  <- ifelse(P[t]<(P_bound[1]),log(K*P[t])-log(K*(P_bound[1])),ifelse(P[t]>P_bound[2],log(K*P[t])-log(K*(P_bound[2])),0)) # penalty if Pmean is outside viable biomass
    # Depletion prior
    penBK[t] <- ifelse(b.yr[t] < 1,0,log(ifelse(b.pr[4]<1,P[t],ifelse(b.pr[4]>1,HtoHmsy[t],BtoBmsy[t])))-log(b.pr[1]))
    }

    # Process error deviation
    for(t in 1:N){
    Proc.Dev[t] <- log(P[t]*K)-log(exp(Pmean[t])*K)}

    # Enforce soft penalties on bounds for P
    #for(t in 1:N){
    #pen.P[t] ~ dnorm(penB[t],1000) # enforce penalty with CV = 0.1
    #pen.bk[t] ~ dnorm(penBK[t],pow(b.pr[2],-2))
    #}

    Hmsy <- r*pow(m-1,-1)*(1-1/m)
    SBmsy_K <- (m)^(-1/(m-1))
    SBmsy <- SBmsy_K*K
    MSY <- SBmsy*Hmsy
    
    
    
    for (t in 1:N)
    {
    meanC[t] <- log(H[t]*P[t]*K)
    TC[t] ~ dlnorm(meanC[t],pow(CV.C[t],-2))
    estC[t] <- exp(meanC[t])
    
    #H[t] <- estC[t]/SB[t]
    HtoHmsy[t] <- H[t]/(Hmsy)
    }

    for (t in 1:(N+1)) # One step ahead biomass
    {
    SB[t] <- K*P[t]
    BtoBmsy[t] <- SB[t]/SBmsy
    }

    ",append=TRUE)}
  
  if(jbinput$settings$CatchOnly==FALSE){
    cat(ebsb.sub("
    # Observation equation in related to EB

    for(i in 1:nI)
    {
    for (t in 1:N)
    {
    Imean[t,i] <- log(q[sets.q[i]]*P[t]*K*EBSBIDX);
    I[t,i] ~ dlnorm(Imean[t,i],(ivar.obs[t,i]));
    CPUE[t,i] ~ dlnorm(Imean[t,i],(ivar.obs[t,i]))
    Ihat[t,i]  <- exp(Imean[t,i])

    }}
    EBSBCODE
  "),append=TRUE)}
  
    if(jbinput$settings$CatchOnly==TRUE){
      cat("
    # Observation equation in related to EB
    one <- sets.q
    for(i in 1:nI)
    {
    for (t in 1:N)
    {
    Imean[t,i] <- log(pow(P[1],-1)*P[t]);
    I[t,i] ~ dlnorm(Imean[t,i],(ivar.obs[t,i]));
    CPUE[t,i] ~ dlnorm(Imean[t,i],(ivar.obs[t,i]))   ####q[[i]]*P[t]*SB0*EBtoSB[t,i]
    Ihat[t,i]  <- exp(Imean[t,i])

    }}

  ",append=TRUE)} 
    
  cat("
      
    
   
    # Enforce soft penalty on K if < K_bounds >
    K.pen ~ dnorm(penK,1000) # enforce penalty
    penK  <- ifelse(K<(K_bounds[1]),log(K)-log(K_bounds[1]),ifelse(K>K_bounds[2],log(K)-log(K_bounds[2]),0)) # penalty if Pmean is outside viable biomass


    # Enforce soft penalty on process deviance if sigma.proc > 0.2
    proc.pen ~ dnorm(penProc,1000) # enforce penalty
    penProc  <- ifelse(sigma>sigmaproc_bound,log(sigma)-log(sigmaproc_bound),0)


    # Enforce soft penalty on observation error if sigma.obs > sigma_bound
    for(i in 1:nvar){
    obs.pen[i] ~ dnorm(penObs[i],1000) # enforce penalty
    penObs[i]  <- ifelse(pow(tau2[i],0.5)>sigmaobs_bound,log(pow(tau2[i],0.5))-log(sigmaobs_bound),0)
    }

    ", append=TRUE)
  # Add Auxiliary observation likelihood
  if(jbinput$settings$Auxiliary==TRUE){

    if(jbinput$settings$auxiliary.sigma==TRUE){
      cat("
      # Observation variance
      for(i in 1:nAvar)
      {
      # Observation error
      ieta2[i]~ dgamma(0.001,0.001)
      eta2[i] <- 1/ieta2[i]
      }

      for(i in 1:nA)
      {
      for(t in 1:N)
      {
      varA.obs[t,i] <- A.SE2[t,i]+eta2[sets.varA[i]]
      ivarA.obs[t,i] <- 1/varA.obs[t,i]
      # note total auxiliary observation error (TAE)
      TAE[t,i] <- sqrt(varA.obs[t,i]) # Total observation variance
      }}
      ",append=TRUE) }
      else{ cat("
      # Observation variance
          
         
          
      for(i in 1:nA)
      {
      for(t in 1:N)
      {
      varA.obs[t,i] <- A.SE2[t,i] # Fixed only
      ivarA.obs[t,i] <- 1/varA.obs[t,i]
      # note total auxiliary observation error (TAE)
      TAE[t,i] <- sqrt(varA.obs[t,i]) # Total observation variance
      }}

           ",append=TRUE)}
    
      if(jbinput$settings$auxiliary.type=="effort"){
      cat("
  
      for(i in 1:nA){ 
      qA[i] ~ dunif(qA_bounds[1],qA_bounds[2])} 
      for(t in 1:N){
      Ax[t] <- H[t]}
      ",append=TRUE)}
      
    if(jbinput$settings$auxiliary.type=="z"){
      cat("
      for(i in 1:nA){
      qA[i] ~ dlnorm(0,pow(qA.cv,-2))T(qA_bounds[1],qA_bounds[2])
      }  
      for(t in 1:N){
      #Ax[t] <- -log(max(1-H[t],0.01))-log(1-max(1-Hmsy,0.01))
      Ax[t] <- max(1-(1-Hmsy)*(1-H[t]),0.001)
      }
      ",append=TRUE)}
    
      if(jbinput$settings$auxiliary.type=="f"){
      cat("
      for(i in 1:nA){
      qA[i] ~ dlnorm(0,pow(qA.cv,-2))T(qA_bounds[1],qA_bounds[2])
      }  
      for(t in 1:N){
      Ax[t] <- -log(1-min(H[t],0.99))
      }
      ",append=TRUE)}
     
    if(jbinput$settings$auxiliary.type=="ffmsy"){
      cat("
      for(i in 1:nA){
      qA[i] ~ dlnorm(0,pow(qA.cv,-2))T(qA_bounds[1],qA_bounds[2])
      }  
      for(t in 1:N){
      Ax[t] <- HtoHmsy[t]
      }
      ",append=TRUE)}
    
    if(jbinput$settings$auxiliary.type=="bbmsy"){
      cat("
      for(i in 1:nA){
      qA[i] ~ dlnorm(0,pow(qA.cv,-2))T(qA_bounds[1],qA_bounds[2])
      }  
      for(t in 1:N){
      Ax[t] <- BtoBmsy[t]
      }
      ",append=TRUE)}
       
    if(jbinput$settings$auxiliary.type=="bk"){
      cat("
      for(i in 1:nA){
      qA[i] ~ dlnorm(0,pow(qA.cv,-2))T(qA_bounds[1],qA_bounds[2])
      }  
      for(t in 1:N){
      Ax[t] <- P[t]
      }
      ",append=TRUE)}
    
    if(jbinput$settings$auxiliary.lag>0){
      cat("
      for(i in 1:nA){ 
      for(t in 1:A.lag){
      
      #Amean[t,i] <- log(qA[i])+log(mean(Ax[1:t]))
      Amean[t,i] <- log(qA[i])+log(mean(Ax[1]))
      
      A[t,i] ~ dlnorm(Amean[t,i],ivarA.obs[t,i])
      AUXI[t,i] ~ dlnorm(Amean[t,i],ivarA.obs[t,i]) 
      Ahat[t,i]  <- exp(Amean[t,i])}}
      
      for(i in 1:nA){
      for(t in (A.lag+1):N){
      Amean[t,i] <- log(qA[i])+log(mean(Ax[(t-A.lag):(t-1)]))
      #Amean[t,i] <- log(qA[i])+log(mean(Ax[(t-A.lag)]))
      
      A[t,i] ~ dlnorm(Amean[t,i],ivarA.obs[t,i])
      AUXI[t,i] ~ dlnorm(Amean[t,i],ivarA.obs[t,i]) 
      Ahat[t,i]  <- exp(Amean[t,i])
      }}
  
    ",append=TRUE)}
    # No Lag
    if(jbinput$settings$auxiliary.lag==0){
     cat("
      for(i in 1:nA){
      for(t in 1:N){
      Amean[t,i] <- log(qA[i])+log(mean(Ax[t]))
      A[t,i] ~ dlnorm(Amean[t,i],ivarA.obs[t,i])
      AUXI[t,i] ~ dlnorm(Amean[t,i],ivarA.obs[t,i]) 
      Ahat[t,i]  <- exp(Amean[t,i])
      }}
      
      ",append=TRUE)
    }
  } # End of auxiliary
  
  # PROJECTION
  if(jbinput$settings$projection==TRUE){
    # Low-B correction for projected states
    if(lowb.mode=="none"){
      cat("
      for(i in 1:nTAC){ for(t in 1:pyrs){ prlowb[t,i] <- 0 }}
      ",append=TRUE)
    } else {
      cat("
      for(i in 1:nTAC){
      prlowb[1,i] <- lowb[N]
      for(t in 2:pyrs){ prlowb[t,i] <- cpe*min(log(max(prBtoBmsy[t-1,i], pe.lowbx)), 0) }}
      ",append=TRUE)
    }
    # Process error multiplier for projected states (prdyn), consistent with the historical period
    if(proc.mode=="floor"){
      cat("
      for(i in 1:nTAC){
      prx[1,i] <- min(max(BtoBmsy[N], pe.xlim[1]), pe.xlim[2])
      for(t in 2:pyrs){ prx[t,i] <- min(max(prBtoBmsy[t-1,i], pe.xlim[1]), pe.xlim[2]) }
      for(t in 1:pyrs){ prdyn[t,i] <- sigma2/(sigma2 + pow(pe.floor[1]*pow(prx[t,i], -kpe), 2)) }}
      ",append=TRUE)
    } else if(proc.mode=="dyn"){
      cat("
      for(i in 1:nTAC){
      prdyn[1,i] <- iprocdyn[N]
      for(t in 2:pyrs){
      prdyn[t,i] <- pow((a_pe + (1 - a_pe) / (1 + pow(prBtoBmsy[t-1,i] / c_pe, d_pe)))/relsig.ref, -2) }}
      ",append=TRUE)
    } else {
      cat("
      for(i in 1:nTAC){ for(t in 1:pyrs){ prdyn[t,i] <- 1 }}
      ",append=TRUE)
    }
    cat("
      for(i in 1:nTAC){
      # Project first year into the future
      prPmean[1,i] <- ifelse(P[N] > Plim,
      log(max(P[N] +  Hmsy/(1-1/m)*P[N]*(1-pow(P[N],m-1)) - TC[N]/K,0.005)),
      log(max(P[N] +  Hmsy/(1-1/m)*P[N]*(1-pow(P[N],m-1))*4*P[N] - TC[N]/K,0.005))) + prlowb[1,i]
      prP[1,i] ~ dlnorm(prPmean[1,i],isigma2*prdyn[1,i])
      # Project all following years
      for(t in 2:pyrs){
      prPmean[t,i] <- ifelse(prP[t-1,i] > Plim,
      log(max(prP[t-1,i] +  Hmsy/(1-1/m)*prP[t-1,i]*(1-pow(prP[t-1,i],m-1)) - TAC[t-1,i]/K,0.005)),
      log(max(prP[t-1,i] +  Hmsy/(1-1/m)*prP[t-1,i]*(1-pow(prP[t-1,i],m-1))*slope.HS*prP[t-1,i] - TAC[t-1,i]/K,0.005))) + prlowb[t,i]
      # process error (as monte-carlo simular)
      prP[t,i] ~ dlnorm(prPmean[t,i],isigma2*prdyn[t,i])}
      for(t in 1:pyrs){
      prB[t,i] <- prP[t,i]*K
      prH[t,i] <- TAC[t,i]/prB[t,i]
      prHtoHmsy[t,i] <- prH[t,i]/Hmsy
      prBtoBmsy[t,i] <- prB[t,i]/SBmsy
      }}
      ",append=TRUE)} else {
        cat("
            #Prevent error for unused input
            fakeTAC <-  TAC
            fakepyrs <- pyrs
            fakenTAC <- nTAC
            prHtoHmsy <- 1
            prP <- 1
            prBtoBmsy <- 1
            ", append=TRUE)}
  
  cat("
} # END OF MODEL
    ",append=TRUE,fill = TRUE)
  sink()
  if(isTRUE(jbinput$settings$proc.biascor)){
    jf = file.path(dir,"JABBA.jags"); txt = readLines(jf)
    a1 = "P[t] ~ dlnorm(Pmean[t],iPV[t]*ifelse(t<(stI),1,iprocdyn[t-1]))"
    txt = gsub(a1, "P[t] ~ dlnorm(Pmean[t]-0.5/(iPV[t]*ifelse(t<(stI),1,iprocdyn[t-1])),iPV[t]*ifelse(t<(stI),1,iprocdyn[t-1]))", txt, fixed=TRUE)
    txt = gsub("P[t] ~ dlnorm(Pmean[t],iPV[t])", "P[t] ~ dlnorm(Pmean[t]-0.5/iPV[t],iPV[t])", txt, fixed=TRUE)
    writeLines(txt, jf)
  }
  if(jbinput$settings$projection==TRUE && isTRUE(jbinput$settings$proc.lag>0)){
    # lagged production in projections: first year uses P[lagidx[N]] as the historical step; later years use the
    # combined historical and projected series allP
    jf = file.path(dir,"JABBA.jags"); txt = readLines(jf)
    txt = gsub("Hmsy/(1-1/m)*P[N]*(1-pow(P[N],m-1))", "Hmsy/(1-1/m)*P[lagidx[N]]*(1-pow(P[lagidx[N]],m-1))", txt, fixed=TRUE)
    txt = gsub("*4*P[N] - TC[N]", "*4*P[lagidx[N]] - TC[N]", txt, fixed=TRUE)
    txt = gsub("prPmean[1,i] <- ifelse(P[N] > Plim", "prPmean[1,i] <- ifelse(P[lagidx[N]] > Plim", txt, fixed=TRUE)
    txt = gsub("Hmsy/(1-1/m)*prP[t-1,i]*(1-pow(prP[t-1,i],m-1))", "Hmsy/(1-1/m)*prPl[t,i]*(1-pow(prPl[t,i],m-1))", txt, fixed=TRUE)
    txt = gsub("*slope.HS*prP[t-1,i]", "*slope.HS*prPl[t,i]", txt, fixed=TRUE)
    txt = gsub("prPmean[t,i] <- ifelse(prP[t-1,i] > Plim", "prPmean[t,i] <- ifelse(prPl[t,i] > Plim", txt, fixed=TRUE)
    txt = gsub("# process error (as monte-carlo simular)", "prPl[t,i] <- allP[prlagidx[t],i]\n      # process error (as monte-carlo simular)", txt, fixed=TRUE)
    txt = gsub("# Project first year into the future", "for(j in 1:N){ allP[j,i] <- P[j] }\n      for(t in 1:pyrs){ allP[N+t,i] <- prP[t,i] }\n      # Project first year into the future", txt, fixed=TRUE)
    writeLines(txt, jf)
  }
  if(isTRUE(jbinput$settings$ebsb.removal)){
    jf = file.path(dir,"JABBA.jags"); txt = readLines(jf)
    txt = gsub("estC[t-1]/K", "estC[t-1]/(K*ebsbr[t-1])", txt, fixed=TRUE)
    # Selex: harvest rate on exploitable biomass, h = C/EB = C/(SB*EB/SB); it removes the same fraction of SB.
    # hmsy = MSY/EBmsy equals Hmsy of the SB production function, and MSY in catch units is SBmsy*Hmsy*(EB/SB)(Pmsy)
    # reference allocation (ebsbref) for MSY, Hmsy and projections
    txt = gsub("H[t] <- estC[t]/SB[t]", "H[t] <- estC[t]/(SB[t]*ebsbr[t])", txt, fixed=TRUE)
    txt = gsub("MSY <- SBmsy*Hmsy", paste0("for(c in 1:ncur){ irhom[c] <- 1/(", ebsb.f("SBmsy_K","c"), ") }\n    ebsbr.msy <- 1/inprod(ebsbref[1:ncur], irhom[1:ncur])\n    MSY <- SBmsy*Hmsy*ebsbr.msy"), txt, fixed=TRUE)
    # projections: removals and harvest rate on EB
    txt = gsub("TC[N]/K", "TC[N]/(K*ebsbr[N])", txt, fixed=TRUE)
    txt = gsub("TAC[t-1,i]/K", "TAC[t-1,i]/(K*prebsb[t-1,i])", txt, fixed=TRUE)
    txt = gsub("prH[t,i] <- TAC[t,i]/prB[t,i]", paste0("for(c in 1:ncur){ iprho[t,i,c] <- 1/(", ebsb.f("prP[t,i]","c"), ") }\n      prebsb[t,i] <- 1/inprod(ebsbref[1:ncur], iprho[t,i,1:ncur])\n      prH[t,i] <- TAC[t,i]/(prB[t,i]*prebsb[t,i])"), txt, fixed=TRUE)
    writeLines(txt, jf)
  }
} # #END jabba JAGS