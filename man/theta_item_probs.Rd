% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/describeDLCM.r
\name{theta_item_probs}
\alias{theta_item_probs}
\title{Calculate the probabilities of different patterns}
\usage{
theta_item_probs(
  items,
  this_sim,
  itrs = NULL,
  merge_itrs = TRUE,
  classes = NULL
)
}
\arguments{
\item{items}{integer vector. Which items should we include? (indexed starting at 1)}

\item{this_sim}{a mcmc simulation. Output from dependentLCM_fit()}

\item{itrs}{integer vector. Which iterations from this_sim should we include?}

\item{merge_itrs}{TRUE if you want average probabilities across simulations. FALSE if you want to examine each iteration individually}

\item{classes}{Which classes are we evaluating? By default all classes}
}
\value{
Probabilities for each pattern of items
}
\description{
Useful when getting probilities which cross domains.
}
