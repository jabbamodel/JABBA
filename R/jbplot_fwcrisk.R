#{{{
#' jbplot_Crisk()
#'
#' Plots plots JABBA ensemble models + projections - joint or by run  
#' 
#' @param prjc objects from fit_jabba(),jabba_fw(), list of fit_jabba() or fit_jabba()$kbtrj    
#' @param subplots option choose from subplots 1:7 
#' \itemize{
#'   \item 1: Risk B < Bmsy  
#'   \item 2: Risk F > Fmsy 
#'   \item 3: Risk B < Bfrac (e.g. Blim, MSST)
#' }    
#' @param riskout if TRUE, produces the kb data.frame as output 
#' @param bfrac biomass fraction of Bmsy or B0 (subplot 8), default 0.5Bmsy
#' @param bref biomass fraction reference options c("bmsy","b0")
#' @param years TODO DOCUMENTATION
#' @param ylabs yaxis labels for quants
#' @param ylab.bref option to only specify BBfrac plot ylab 
#' final year of values to show for each model. By default it is set to the
#' @param xlab xaxis label
#' @param plot TODO DOCUMENTATION
#' @param as.png TODO DOCUMENTATION
#' @param col Optional vector of colors to be used for lines. Input NULL
#' @param pch Optional vector of plot character values
#' @param lty Optional vector of line types
#' @param lwd Optional vector of line widths
#' @param tickEndYr TRUE/FALSE switch to turn on/off extra axis mark at final
#' year in timeseries plots.
#' @param ylimAdj Multiplier for ylim parameter. Allows additional white space
#' @param xlim = NULL range of years
#' @param xaxs Choice of xaxs parameter (see ?par for more info)
#' @param yaxs Choice of yaxs parameter (see ?par for more info)
#' @param xylabs TODO DOCUMENTATION
#' @param type Type parameter passed to points (default 'o' overplots points on
#' top of lines)
#' @param legend Add a legend?
#' @param legendlabels Optional vector of labels to include in legend.
#' @param legend.loc Location of legend. Either a string like "topleft" or a vector
#' of two numeric values representing the fraction of the maximum in the x and y
#' dimensions, respectively. See ?legend for more info on the string options.
#' @param legendorder Optional vector of model numbers that can be used to have
#' the legend display the model names in an order that is different than that
#' which is represented in the summary input object.
#' @param legendncol Number of columns for the legend.
#' @param legendcex Defeult=1 Allows to adjust legend cex
#' @param legendsp Space between legend labels
#' @param pwidth Width of plot
#' @param pheight Height of plot
#' @param punits Units for PNG file
#' @param res Resolution for PNG file
#' @param ptsize Point size for PNG file
#' @param cex.main Character expansion for plot titles
#' @param plotdir Directory where PNG or PDF files will be written. By default
#' it will be the directory where the model was run.
#' @param filenameprefix Additional text to append to PNG or PDF file names.
#' It will be separated from default name by an underscore.
#' @param par list of graphics parameter values passed to par() function
#' @param verbose Report progress to R GUI?
#' @param shadecol uncertainty shading of hcxval horizon
#' @param shadealpha Transparency adjustment used to make default shadecol
#' @param new Create new empty plot window
#' @param add surpresses par() to create multiplot figs
#' @param run name for single models or joint ensembles
#' @param single.plots TODO DOCUMENTATION
#' @param fmax TODO DOCUMENTATION
#' @author Mostly adopted from ss3diags::SSplotEnsemble
#' @importFrom grDevices graphics.off adjustcolor png dev.off
#' @importFrom stats aggregate 
#' @importFrom graphics lines abline axis box
#' @export
#' @examples
#' data(iccat)
#' bet = iccat$bet 
#' # Fit Fox and Schaefer
#' jb1 <- build_jabba(catch=bet$catch,cpue=bet$cpue,se=bet$se,scenario = "Fox",model.type="Fox")
#' fit1 = fit_jabba(jb1,quickmcmc=TRUE,verbose=TRUE)
#' # Compare
#' jbplot_ensemble(fit1)
#' # Do forecast
#' prjc = fw_jabba(list(fit1,fit2),quant="Catch",type="abs",imp.values = seq(60,100,1)*1000)
#' jbplot_Crisk(prjc,subplots=1)
#' jbplot_Crisk(prjc,subplots=2)
#' jbplot_Crisk(prjc,subplots=3) # MSST as (1-0.5)*Bmsy with bfrac =0.5
#' jbplot_Crisk(prjc,subplots=3, bfrac=0.2,bref="b0") # Blim = 0.2 B0
#' jbplot_Crisk(prjc,subplots=3, bfrac=0.2,bref="b0",ylabs=expression(B/B[lim])) # Blim = 0.2 B0

jbplot_Crisk<- function(prjc,
                        subplots=2,
                        bfrac = 0.5,
                        bref = c("bmsy","b0")[1],
                        years = NULL,
                        riskout = FALSE,
                        ylabs = NULL,
                        ylab.bref = NULL,
                        xlab="Catch",
                        plot=TRUE,as.png=FALSE,
                        col=NULL, 
                        pch=NULL, lty=1, lwd=1.5,
                        tickEndYr=FALSE,
                        xlim=NULL, ylimAdj=1.05,
                        xaxs="i", yaxs="i",
                        xylabs=TRUE,
                        type="l", 
                        legend=TRUE, legendlabels="default", legend.loc="topleft",
                        legendorder="default",legendncol=1,legendcex=0.7,legendsp=0.8,
                        pwidth=6.5,pheight=5.0,punits="in",res=300,ptsize=10,cex.main=1,
                        plotdir=NULL,
                        filenameprefix="",
                        par=list(mar=c(5,4,1,1)+.1),
                        verbose=FALSE,
                        shadecol = NULL, shadealpha=0.3,new=TRUE,
                        add=FALSE,
                        single.plots = add,
                        fmax=5.0
                        ){ # plot different fits to a single index of abundance
     
      
   kb = prjc[prjc$type=="prj",]
   endyr = max(kb$year)
   if(is.null(years)) years = unique(kb$year)[-ifelse(length(unique(kb$year))>4,2,1)]
   kb = kb[kb$year%in%years,]
   
   kb = addBfrac(kb)$kb
   kb$riskBmsy = ifelse(kb$stock<1,1,0)
   kb$riskBfrac = ifelse(kb$BBfrac<1,1,0)
   kb$riskFmsy = ifelse(kb$harvest>1,1,0)
   
   
  if(as.png==TRUE){
    add=FALSE
  }
   if(!add) graphics.off()
   
  if(is.null(ylabs)){
    if(!is.null(bfrac) & bref[1]=="bmsy") bflab = bquote(Risk ~ B < ~(.(bfrac)~B[MSY])) 
    if(!is.null(bfrac) & bref[1]=="b0") bflab = bquote(Risk ~ B < ~(.(bfrac)~B[0]))
    if(is.null(bfrac)) bflab = "Not defined"
    if(!is.null(ylab.bref)) bflab = ylab.bref
    ylab.default = TRUE
    ylabs =  c(expression("Risk B <" ~ B[MSY]),expression("Risk F > " ~ F[MSY]),bflab)
  } else {
    ylab.default = FALSE
  }
  
  refquants=c("riskBmsy","riskFmsy","riskBfrac")
  
  kb = addBfrac(kb,bfrac=bfrac,bref=bref)$kb
  
  if(subplots==1) kbs= aggregate(cbind(riskBmsy)~Catch+year,kb,mean)
  if(subplots==2) kbs= aggregate(cbind(riskFmsy)~Catch+year,kb,mean)
  if(subplots==3) kbs= aggregate(cbind(riskBfrac)~Catch+year,kb,mean)
  
  n             <- length(years)
  if(is.null(col)){
  col = ss3col(n,1)
  shadecol <- ss3col(n,shadealpha)
  } else {
    shadecol= adjustcolor(col, alpha.f=shadealpha)
  }
  
  quants =  refquants[subplots]
 
  pngfun <- function(file){
    
    
    # if extra text requested, add it before extention in file name
    file <- paste0(filenameprefix, file)
    # open png file
    png(filename=file.path(plotdir,file),
        width=pwidth,height=pheight,units=punits,res=res,pointsize=ptsize)
    # change graphics parameters to input value
    par(par)
  }
  
  print=FALSE
  if(as.png) print <- TRUE
  if(as.png & is.null(plotdir)) plotdir = getwd()
  
  if(!single.plots){
  Par = list(mfrow=c(1,1),mai=c(0.45,0.49,0.1,.15),omi = c(0.15,0.15,0.1,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
  #if(as.png==TRUE){png(filename = paste0(plotdir,"/Risk_",jbs$assessment,".png"), width = 7, height = 6,
  if(as.png==TRUE){png(filename = paste0(plotdir,"/Risk_","jbs$assessment",".png"), width = 7, height = 6,
                       res = 200, units = "in")}
  par(Par)
  }
  
  plot_quants <- function(quant="riskFmsy"){  
   
    if(single.plots){ 
    if(as.png) print <- TRUE
    if(as.png & is.null(plotdir))
      stop("to print PNG files, you must supply a directory as 'plotdir'")
    }
    #-------------------------------------------------------------
    # plot function
    #-------------------------------------------------------------
    # get stuff from summary output (minimized)
    y = kbs[,quant]   
    C = kbs$Catch
    exp      <- y
    models <- 1:n    
    nlines <- length(models) 
    runs = years
    
   
    # if line stuff is shorter than number of lines, recycle as needed
    if(length(col) < nlines) col <- rep(col,nlines)[1:nlines]
    if(length(pch) < nlines) pch <- rep(pch,nlines)[1:nlines]
    if(length(lty) < nlines) lty <- rep(lty,nlines)[1:nlines]
    if(length(lwd) < nlines) lwd <- rep(lwd,nlines)[1:nlines]
    
    
    # subfunction to add legend
    legendfun <- function(legendlabels,cumulative=FALSE) {
      if(cumulative){
        legend.loc="topleft"
      }
      if(is.numeric(legend.loc)) {
        Usr <- par()$usr
        legend.loc <- list(x = Usr[1] + legend.loc[1] * (Usr[2] - Usr[1]),
                          y = Usr[3] + legend.loc[2] * (Usr[4] - Usr[3]))
      }
      
      # if type input is "l" then turn off points on top of lines in legend
      legend.pch <- pch
      if(type=="l"){
        legend.pch <- rep(NA,length(pch))
      }
      legend(legend.loc, legend=legendlabels[legendorder],
             col=col[legendorder], lty=lty[legendorder],seg.len = 2,
             lwd=lwd[legendorder], pch=legend.pch[legendorder], bty="n", ncol=legendncol,pt.cex=0.7,cex=legendcex,y.intersp = legendsp)
    }
    
    
    if(!is.expression(legendlabels[1]) &&
       legendlabels[1]=="default") legendlabels <- runs
    if(legendorder[1]=="default") legendorder <- 1:(nlines)
    
    
    
    # open new window if requested
    if(single.plots){
    if(plot & as.png==FALSE){
      if(!add) par(par)
      
    } else {
      
      if(!add) par(par)
    }
    }
    
   if(is.null(xlim)) xlim = c(max(min(C)),max(C)) 
    xmin = min(xlim)
    ylim <- c(0,1*ylimAdj)
    if(ylab.default){
    ylab = ylabs[which(refquants%in%quant)]} else {
    ylab = ylabs[which(quants%in%quant)]  
    }
    
    
    plot(0, type = "n", xlim = xlim, yaxs = yaxs, 
         ylim = ylim, xlab = xlab, ylab = ifelse(xylabs,ylab,""), axes = FALSE,cex.lab=0.9)
    
    
    
    for(iline in 1:nlines){
      catch <- kbs[kbs$year==runs[iline],]$Catch  
      lines(catch,exp[kbs$year == runs[iline]],col=col[iline],pch=pch[iline],lty=lty[iline],lwd=lwd[iline],type="l")
      }  
    abline(h=c(0.2,0.5,0.8),lty=c(3,2,3))
    
    
    if(legend){
      # add legend if requested
      
      legendfun(paste0(legendlabels))
    }
    #axis(1, at=c(min(xmin,min(yr)):max(endyrvec)))
    axis(1,cex.axis=0.8)
    
    if(tickEndYr) axis(1, at=max(xlim[2]),cex.axis=0.8)
    
    axis(2,at=seq(0,1,0.1),cex.axis=0.8)
    box()
  } # End of plot_quant function  
  legend.temp = legend  
  # Do plotting
  if(plot){ 
    # subplots
    for(s in 1:length(subplots)){
    if(print & single.plots){
    quant=quants[s]
    par(par)
    pngfun(paste0("ModelComp_",quant,".png",sep=""))
    plot_quants(quant)
    dev.off()
    }
    }
    # subplots
    for(s in 1:length(subplots)){
      if(verbose) cat(paste0("\n","Plot Comparison of ",quants[s],"\n"))
      
    if(!add & single.plots)par(par)
    quant=quants[s]
    plot_quants(quant)   
    }
    if(riskout) return(kbs)
  } # endplot

} 
#}}} end of jbplot_ensemble()
#-----------------------------------------------------------------------------------------

