% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jabba_hindcast.R
\name{jabba_hindcast}
\alias{jabba_hindcast}
\title{DEPRECIATED}
\usage{
jabba_hindcast(
  jbinput,
  ni = 30000,
  nt = 5,
  nb = 5000,
  nc = 2,
  init.values = FALSE,
  init.K = NULL,
  init.r = NULL,
  init.q = NULL,
  peels = 0:5,
  save.jabba = FALSE,
  output.dir = getwd(),
  save.hc = FALSE,
  plotall = FALSE,
  speedup = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{jbinput}{List of input variables as output by build_jabba()
MCMC settings}

\item{ni}{number of iterations}

\item{nt}{thinning interval of saved iterations}

\item{nb}{burn-in}

\item{nc}{number of mcmc chains
Initial values}

\item{init.values}{= FALSE,}

\item{init.K}{= NULL,}

\item{init.r}{= NULL,}

\item{init.q}{= NULL,# vector}

\item{peels}{sequence of retro spective peels default 0:5}

\item{save.jabba}{Save individual JABBA run outputs}

\item{output.dir}{path to save plot. default is getwd()}

\item{save.hc}{Save hindcast list output as .rdata}

\item{plotall}{if TRUE makes jabba_plots() for each run}

\item{speedup}{Reduces MCMC after setting runs 2+ inits to first "full" reference run}

\item{verbose}{if FALSE run silent}
}
\value{
hc containing estimates of key joint results from all hindcast run
}
\description{
PLEASE USE hindcast_jabba() as the new wrapper
}
