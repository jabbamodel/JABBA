#{{{
#' jbplot_ensemble()
#'
#' Plots plots JABBA ensemble models + projections - joint or by run  
#' 
#' @param kb objects from fit_jabba(),jabba_fw(), list of fit_jabba() or fit_jabba()$kbtrj    
#' @param subplots option to c("stock","harvest","B","H","Bdev","Catch","BB0") 
#' \itemize{
#'   \item stock (B/Bmsy)  
#'   \item harvest (F/Fmsy) 
#'   \item B (Biomass) 
#'   \item H (Harvest rate)
#'   \item Bdev (Process Deviations)
#'   \item Catch
#'   \item BB0 (B/K) 
#' }    
#' @param joint if true it creates a joint ensemble from list of fit_jabba()
#' @param kbout if TRUE, produces the kb data.frame as output 
#' @param ylabs yaxis labels for quants
#' final year of values to show for each model. By default it is set to the
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
#' @param type Type parameter passed to points (default 'o' overplots points on
#' top of lines)
#' @param legend Add a legend?
#' @param legendlabels Optional vector of labels to include in legend.
#' @param legendloc Location of legend. Either a string like "topleft" or a vector
#' of two numeric values representing the fraction of the maximum in the x and y
#' dimensions, respectively. See ?legend for more info on the string options.
#' @param legendorder Optional vector of model numbers that can be used to have
#' the legend display the model names in an order that is different than that
#' which is represented in the summary input object.
#' @param legendncol Number of columns for the legend.
#' @param legendcex=1 Allows to adjust legend cex
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
#' @author Mostly adopted from ss3diags::SSplotEnsemble
#' @export
jbplot_ensemble<- function(kb,
                        subplots=c("B","H","stock","harvest","Bdev","Catch","BB0")[1:6],
                        joint=FALSE,
                        quantiles = c(0.025,0.975),
                        kbout = FALSE,
                        ylabs = NULL,
                        plot=TRUE,print=FALSE,png=print,
                        col=NULL, 
                        pch=NULL, lty=1, lwd=2,
                        tickEndYr=FALSE,
                        xlim=NULL, ylimAdj=1.05,
                        xaxs="i", yaxs="i",
                        xylabs=TRUE,
                        type="l", uncertainty=TRUE, 
                        legend=TRUE, legendlabels="default", legendloc="topleft",
                        legendorder="default",legendncol=1,legendcex=0.7,legendsp=0.8,
                        pwidth=6.5,pheight=5.0,punits="in",res=300,ptsize=10,cex.main=1,
                        plotdir=NULL,
                        filenameprefix="",
                        par=list(mar=c(5,4,1,1)+.1),
                        verbose=FALSE,
                        shadecol = NULL, shadealpha=0.3,new=TRUE,
                        add=TRUE,
                        run=NULL,
                        fmax=5.0
                        ){ # plot different fits to a single index of abundance
     
    if(!is.null(kb$settings) & length(unique(kb$kbtrj$run)>1)){ 
      kb = kb$kbtrj
     }      
     
      # if a list of fit_jabba() is provided
      if(class(kb)=="list"){   
      
      if(!is.null(kb$settings)){
        kb = list(kb)
        if(!is.null(run)) names(kb) = run
      }
      
        if(is.null(names(kb)[1])){
          run.ls = do.call(c,lapply(kb,function(x){
          x$scenario
          })) } else {
          run.ls = names(kb)
          }
          run.ls = as.list(run.ls)
    
          
          kb = do.call(rbind,Map(function(x,y){
          z = x$kbtrj
          z$run = y
          z
          },x=kb,y=run.ls))
          
          
          if(joint & !is.null(run)) kb$run = run
          if(joint & is.null(run)) kb$run = "Joint"
      }
  
   if(!is.null(xlim)) kb =kb[kb$year<=xlim[2] & kb$year>=xlim[1],]
  
   # Contraint on F/Fmsy
   kb$harvest[kb$type=="prj"] = pmin(kb[kb$type=="prj",]$harvest,fmax)
   kb$H[kb$type=="prj"]= pmin(fmax*median(kb[kb$type=="prj",]$H/kb[kb$type=="prj",]$harvest),kb[kb$type=="prj",]$H)
   
  if(print==TRUE | png==TRUE){
    add=FALSE
  }
   if(!add) graphics.off()
   
  if(is.null(ylabs)){
    ylab.default = TRUE
    ylabs =  c(expression(B/B[MSY]),expression(F/F[MSY]),"Biomass (t)","Fishing Mortality","Process Deviation","Catch (t)",expression(B/B[0]))
  } else {
    ylab.default = FALSE
  }
  
  C
  
  refquants=c("stock","harvest","B","H","Bdev","Catch","BB0")
  
  kbs = aggregate(cbind(stock,harvest,B,H,Bdev,Catch,BB0)~year+run,kb,
                   quantile, c(0.5,quantiles))
  
  n             <- length(unique(kbs$run))
  startyrs      <- min(kbs$year)
  endyrs        <- max(kbs$year)
  years         <- unique(kbs$year)
  
  col = ss3col(n,1)
  shadecol <- ss3col(n,shadealpha)
  quants = subplots
 
  pngfun <- function(file){
    
    
    # if extra text requested, add it before extention in file name
    file <- paste0(filenameprefix, file)
    # open png file
    png(filename=file.path(plotdir,file),
        width=pwidth,height=pheight,units=punits,res=res,pointsize=ptsize)
    # change graphics parameters to input value
    par(par)
  }
  
  
  if(png) print <- TRUE
  if(png & is.null(plotdir)) plotdir = getwd()
  
  
  
  
  plot_quants <- function(quant="Bdev"){  
    
    if(png) print <- TRUE
    if(png & is.null(plotdir))
      stop("to print PNG files, you must supply a directory as 'plotdir'")
    #-------------------------------------------------------------
    # plot function
    #-------------------------------------------------------------
    # get stuff from summary output (minimized)
    y = kbs[,quant]   
    Yr = kbs$year
    exp      <- y[,1]
    lower    <- y[,2]
    upper <- y[,3]
    models <- 1:n    
    nlines <- length(models) 
    runs = unique(kb$run)[models]
    
   
    # if line stuff is shorter than number of lines, recycle as needed
    if(length(col) < nlines) col <- rep(col,nlines)[1:nlines]
    if(length(pch) < nlines) pch <- rep(pch,nlines)[1:nlines]
    if(length(lty) < nlines) lty <- rep(lty,nlines)[1:nlines]
    if(length(lwd) < nlines) lwd <- rep(lwd,nlines)[1:nlines]
    
    
    # subfunction to add legend
    legendfun <- function(legendlabels,cumulative=FALSE) {
      if(cumulative){
        legendloc="topleft"
      }
      if(is.numeric(legendloc)) {
        Usr <- par()$usr
        legendloc <- list(x = Usr[1] + legendloc[1] * (Usr[2] - Usr[1]),
                          y = Usr[3] + legendloc[2] * (Usr[4] - Usr[3]))
      }
      
      # if type input is "l" then turn off points on top of lines in legend
      legend.pch <- pch
      if(type=="l"){
        legend.pch <- rep(NA,length(pch))
      }
      legend(legendloc, legend=legendlabels[legendorder],
             col=col[legendorder], lty=lty[legendorder],seg.len = 2,
             lwd=lwd[legendorder], pch=legend.pch[legendorder], bty="n", ncol=legendncol,pt.cex=0.7,cex=legendcex,y.intersp = legendsp)
    }
    
    
    if(!is.expression(legendlabels[1]) &&
       legendlabels[1]=="default") legendlabels <- runs
    if(legendorder[1]=="default") legendorder <- 1:(nlines)
    
    
    
    # open new window if requested
    if(plot & png==FALSE){
      if(!add) par(par)
      
    } else {
      
      if(!add) par(par)
    }
    
   if(is.null(xlim)) xlim = c(max(min(years)),max(years)) 
    xmin = min(xlim)
    ylim <- c(0,max(ifelse(uncertainty,max(upper[Yr>=xmin])*ylimAdj, ylimAdj*max(exp[Yr>=xmin])*1.05)))
    if(quant=="Bdev") ylim <- c(-max(ifelse(uncertainty,max(c(0.2,upper[Yr>=xmin],abs(lower[Yr>=xmin])))*ylimAdj, ylimAdj*max(abs(exp[Yr>=xmin]))*1.05)),max(0.2,ifelse(uncertainty,max(c(upper[Yr>=xmin],abs(lower[Yr>=xmin])))*ylimAdj, ylimAdj*max(abs(exp[Yr>=xmin]))*1.05)))
    
      
    if(ylab.default){
    ylab = ylabs[which(refquants%in%quant)]} else {
    ylab = ylabs[which(subplots%in%quant)]  
    }
    
    
    plot(0, type = "n", xlim = xlim, yaxs = yaxs, 
         ylim = ylim, xlab = ifelse(xylabs,"Year",""), ylab = ifelse(xylabs,ylab,""), axes = FALSE,cex.lab=0.9)
    
    if(uncertainty){
    for(iline in nlines:1){
    yr <- kbs[kbs$run==runs[iline],]$year  
    if(quant%in%c("B","stock","harvest","H","Bdev","BB0","Catch")){  
       polygon(c(yr,rev(yr)),c(lower[kbs$run == runs[iline]],rev(upper[kbs$run == runs[iline]])),col=shadecol[iline],border=shadecol)
    } else {
      adj <- 0.2*iline/nlines - 0.1
      arrows(x0=yr+adj, y0=lower[kbs$run == runs[iline]],
      x1=yr+adj, y1=upper[kbs$run == runs[iline]],
      length=0.02, angle=90, code=3, col=col[iline])
    }}
    }
    
    for(iline in 1:nlines){
      yr <- kbs[kbs$run==runs[iline],]$year  
      if(quant%in%c("B","stock","harvest","H","Bdev","BB0","Catch")){
        lines(yr,exp[kbs$run == runs[iline]],col=col[iline],pch=pch[iline],lty=lty[iline],lwd=lwd[iline],type="l")
      } else {
        points(yr,exp[kbs$run == runs[iline]],col=col[iline],pch=16,cex=0.8)
      } 
   
    }  
    if(quant == "stock") abline(h=1,lty=2)
    if(quant == "harvest") abline(h=1,lty=2)
    if(quant == "Bdev") abline(h=0,lty=2)
    
    if(legend){
      # add legend if requested
      
      legendfun(legendlabels)
    }
    
    #axis(1, at=c(min(xmin,min(yr)):max(endyrvec)))
    axis(1,cex.axis=0.8)
    
    if(tickEndYr) axis(1, at=max(xlim[2]),cex.axis=0.8)
    
    axis(2,cex.axis=0.8)
    box()
  } # End of plot_quant function  
  legend.temp = legend  
  # Do plotting
  if(plot){ 
    # subplots
    for(s in 1:length(subplots)){
    if(print){
    quant=subplots[s]
    par(par)
    pngfun(paste0("ModelComp_",quant,".png",sep=""))
    plot_quants(quant)
    dev.off()
    }
    }
    # subplots
    for(s in 1:length(subplots)){
      if(verbose) cat(paste0("\n","Plot Comparison of ",subplots[s],"\n"))
    if(subplots[s]!="Index"){  
    if(!add)par(par)
    quant=subplots[s]
    plot_quants(quant)   
    }else{  
       nfleets=length(unique(summaryoutput$indices$Fleet))
      
        for(fi in 1:nfleets){
        legend=F
        if(fi%in%legendindex) legend=TRUE
        indexfleets = unique(summaryoutput$indices$Fleet)[fi] 
        if(!add)par(par)
        plot_index(indexfleets)   
        legend = legend.temp 
      } # End of Fleet Loop
    }
    }  
    if(kbout) return(kb)
  } # endplot

} 
#}}} end of jbplot_ensemble()
#-----------------------------------------------------------------------------------------

