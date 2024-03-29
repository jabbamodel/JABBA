% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jabba_utils.R
\name{zage}
\alias{zage}
\title{zage()}
\usage{
zage(ca, ages = "missing")
}
\arguments{
\item{ca}{data.frame input with column names year, age, data}

\item{ages}{for which the z slope is taken}
}
\value{
data.frame with data = z
}
\description{
Computes the Z[t] from catch-at-age data
}
