% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dependentLCM_fit.r
\name{getStart_domains}
\alias{getStart_domains}
\title{Choose starting domain values.
Note that the initial choice of theta (although set) is unimportant since we recalculate theta before applying it in MCMC/CPP.}
\usage{
getStart_domains(mat, classes, hparams)
}
\arguments{
\item{mat}{matrix. Raw data.}

\item{classes}{integer vector. The class of each observation.}

\item{hparams}{list. List of hyperparameters}
}
\description{
Choose starting domain values.
Note that the initial choice of theta (although set) is unimportant since we recalculate theta before applying it in MCMC/CPP.
}
\keyword{internal}
