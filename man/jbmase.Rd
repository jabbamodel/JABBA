% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jabba_utils.R
\name{jbmase}
\alias{jbmase}
\title{jbmase()}
\usage{
jbmase(hc, naive.min = 0.1, index = NULL, residuals = FALSE, verbose = TRUE)
}
\arguments{
\item{hc}{object list of hindcasts from hindcast_jabba() or jbhcxval()}

\item{naive.min}{minimum MASE denominator (naive predictions) for MASE.adj (default = 0.1)}

\item{index}{option to compute for specific indices (numeric & in order)}

\item{residuals}{if TRUE, outputs individual prediction and naive residuals}

\item{verbose}{if FALSE run silent}
}
\value{
hc containing estimates of key joint results from all hindcast run
}
\description{
Computes Mean Absolute Scaled Errors as a measure of prediction skill
}
