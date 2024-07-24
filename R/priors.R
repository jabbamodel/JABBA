
#' jbio()
#'
#' Creates a list for JABBA-Select stock object as input for the ASEM
#' @param amax maximum age
#' @param amin minimum age
#' @param nsexes if nsexes = 2, values are in order female, male
#' @param Loo asymptotic length
#' @param k brody growth coefficient
#' @param t0 theoretical age at zero length (vb growth)
#' @param aW parameter of length-weigth Wa = aW*La^bW
#' @param bW parameter of length-weigth Wa = aW*La^bW
#' @param mat c(lm,dm,0=len or 1 = age) logistic maturity ogive parameters
#' \itemize{
#'   \item m50 - age or length
#'   \item m95 - age or length
#'   \item P(mat) = 1/(1+exp(-log(19)*(x-m50)/(m95-x)))
#'   \item x - age or length
#' }
#' @param M Natural mortality (vectors of M at age can be adjusted after stk is created)
#' @param h steepness of BevHolt SSR
#' @param fec number of offspring per female - for sharks and rays
#' @param empty produces an empty list structure based on amin:amax if TRUE
#' @return jbio List of life history vectors and matrices
#' @export

jbio <- function(amin = 0,amax = 10,nsexes=2,nseas=1,Loo = c(80,65),k=c(0.2,0.25),t0=c(-0.5,-0.7),
                 aW = 0.01,bW=3.04,mat=c(40,45,0),M=c(0.3,0.4),h=0.7,fec=NA,empty=FALSE){
  stk = NULL
  # object
  age = seq(amin,amax,1/nseas)
  if(nseas>1) age = age[age>0]
  nage = length(age)
  stk$nsexes = nsexes
  stk$nseas = nseas
  stk$age = age
  stk$La = matrix(NA,nrow=nage,ncol=2)
  stk$Wa = matrix(NA,nrow=nage,ncol=2)
  stk$Mat = matrix(NA,nrow=nage,ncol=2)
  stk$M = matrix(NA,nrow=nage,ncol=2)
  stk$h = NA
  if(empty==F){
    if(NA%in%c(Loo[1],k[1],t0[1])==FALSE){
      stk$La[,1] = Loo[1]*(1-exp(-k[1]*(age-t0[1])))
      if(nsexes>1)
        stk$La[,2] = Loo[length(Loo)]*(1-exp(-k[length(k)]*(age-t0[length(t0)])))
      if(nsexes==1)
        stk$La[,2] = stk$La[,1]
      
      stk$La[stk$La<0] = 0.01
      
    }
    if(NA%in%c(aW[1],bW[1])==FALSE){
      stk$Wa[,1] = aW[1]*stk$La[,1]^bW[1]
      if(nsexes>1)
        stk$Wa[,2] = aW[length(aW)]*stk$La[,2]^bW[length(bW)]
      
      if(nsexes==1)
        stk$Wa[,2] = stk$Wa[,1]
    }
    
    if(NA%in%c(mat)==FALSE){
      if(mat[3]==0){x = stk$La[,1]} else {x= stk$age}
      stk$Mat = 1/(1+exp(-log(19)*(x-mat[1])/(mat[2]-mat[1])))}
    if(NA%in%c(M[1])==FALSE){
      stk$M[,1] = rep(M[1]/nseas,nage)
      if(nsexes>1)
        stk$M[,2] = rep(M[1]/nseas,nage)
      if(nsexes==1)
        stk$M[,2] = stk$M[,1]
    }
    if(is.na(h)==FALSE) stk$h = h
    if(is.na(fec)==FALSE){stk$fec =  stk$Mat*fec} else {stk$fec=NULL}
  }
  return(stk)
}

#' spr0() 
#'
#' Function to compute spr0
#' @param jbio object list from jbio()
#' @return spr0
#' @param byage if TRUE it return spr0_at_age
#' @return spr0 value(s)
#' @export
jbspr0 <- function(bio=jbio(),byage=FALSE){
survivors = exp(-cumsum(bio$M[,1]))
survivors[-1] = survivors [-length(survivors)]
survivors[1]=1
expZ=exp(-bio$M[length(bio$M[,1]),1])
survivors[length(bio$M[,1])]=survivors[length(bio$M[,1])]*(-1.0/(expZ-1.0))
fec=bio$Mat*bio$Wa[,1]
out = fec * survivors
if(!byage) out = sum(out) 
return(out)
}

#' jbleslie
#'
#' Function to compute r, generation time (gt) based on female life history 
#' @param bio list of history parameters 
#' \itemize{
#'   \item if fecundity is provided (e.g. pups per females) it will be used directly  
#'   \item if fecundity is NULL, a BevHolt as a function is used to compute the annual reproductive rate
#' }
#' @return list r, gt
#' @export
#' @examples
#' #Blond ray
#' # Use fecundity: 30 per female
#' bio.fec = jbio(amin = 0,amax = 20,Loo = c(154.7),k=c(0.129),t0=c(-1),aW = 0.01,bW=3.04,mat=c(5.5,6,1),M=c(0.26),fec=32,empty=FALSE)
#' jbleslie(bio.fec)
#' # Use steepness: h = 0.5
#' bio.h = jbio(amin = 0,amax = 20,Loo = c(154.7),k=c(0.129),t0=c(-1),aW = 0.01,bW=3.04,mat=c(5.5,6,1),M=c(0.26),h=0.5,empty=FALSE)
#' jbleslie(bio.h)
#' #
#' # Thornback ray
#' # Use fecundity: 100 per female
#' bio.fec = jbio(amin = 0,amax = 18,Loo = c(117.6),k=c(0.16),t0=c(-1),aW = 0.01,bW=3.04,mat=c(5.3,5.8,1),M=c(0.32),fec=100,empty=FALSE)
#' jbleslie(bio.fec)
#' # Use steepness: h = 0.5
#' bio.h = jbio(amin = 0,amax = 18,Loo = c(117.6),k=c(0.16),t0=c(-1),aW = 0.01,bW=3.04,mat=c(5.3,5.8,1),M=c(0.32),h=0.5,empty=FALSE)
#' jbleslie(bio.h)
#' 
jbleslie <- function(bio){ 
  nage = length(bio$age)
  age = bio$age
  if(is.null(bio$fec)){
    spr0 = jbspr0(bio,byage=FALSE)
    spr0_a = jbspr0(bio,byage=TRUE)
  # Reproductive output Rs for bonyfish
  rs = 4*bio$h/(spr0*(1-bio$h))
  wm = bio$Wa[,1]*bio$Mat
  mx = rs*wm
  GT = sum(age*spr0_a)/spr0
  } else {
    lx = exp(-cumsum(bio$M[,1]))
    lx[-1] = lx [-length(lx)]
    lx[1]=1
    expZ=exp(-bio$M[length(bio$M[,1]),1])
    lx[length(bio$M[,1])]=lx[length(bio$M[,1])]*(-1.0/(expZ-1.0)) 
    mx = lx*bio$fec/2
    R0 = sum(mx) # Net reproductive rate
    GT = sum(age*mx)/R0 #generation time
  }
  

  # Make Leslie matrix 
    L=mat.or.vec(nage,nage)
    L[1,] =   mx
    #fill rest of Matrix with Survival
    for(i  in 2:nage) L[i,(i-1)] = exp(-bio$M[i,1]) 
    # Plus group
    L[nage,nage] = exp(-bio$M[nage,1]) 
    # Net reproductive rate
    r=log(as.numeric(eigen(L)$values[1]))

  
  # return intrinsic rate of population increase r and generation GT
  return(list(r=r,gt=GT))
}



