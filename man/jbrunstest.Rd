% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jabba_utils.R
\name{jbrunstest}
\alias{jbrunstest}
\title{JABBA runs test}
\usage{
jbrunstest(jabba, index = NULL, mixing = "less")
}
\arguments{
\item{jabba}{output list from fit_jabba}

\item{index}{option to plot specific indices (numeric & in order)}

\item{mixing}{c("less","greater","two.sided"). Default "less" is checking for positive autocorrelation only}
}
\description{
Residual diagnostics with runs test p-value
}
\examples{
data(iccat)
bet= iccat$bet
jb = build_jabba(catch=bet$catch,cpue=bet$cpue,se=bet$se,assessment="BET",scenario = "Ref",model.type = "Pella",igamma = c(0.001,0.001),verbose=FALSE)
fit = fit_jabba(jb,quickmcmc=TRUE,verbose=FALSE)
jbrunstest(fit)
jbrunstest(fit,index=2)
jbplot_runstest(fit,verbose=FALSE)
}
