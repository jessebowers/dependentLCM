% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/describeDLCM.r
\name{dependence_intensity}
\alias{dependence_intensity}
\title{For each domain domain, compare the probabilities under local dependence (DLCM) versus under local independence (traditional LCM). Describes the amount of dependence captured by our DLCM model.}
\usage{
dependence_intensity(thetas_avg, dlcm, items_ids = NULL, class_pi = NULL)
}
\arguments{
\item{thetas_avg}{dataframe. Describes average response probabilities for different domains. As returned from dlcm.summary()}

\item{dlcm}{list. Fitted domain LCM}

\item{items_ids}{integerVector. Optional. Which domains do we want to process?}

\item{class_pi}{numericVector. Optional. Prior probability of being in each class. Used to aggregate certain terms across classes.}
}
\value{
Returns two objects.
"ratio_dfs" has one element per domain. For each domain we examine each response pattern and compare the response probability under conditional dependence versus conditional independence.
\itemize{
\item{"class"}{= Which class does this apply to? NA indicates this term is aggregated across classes.}
\item{"items_id"}{= Unique identifier for the items in the domain. Matches dlcm$mcmc$domains$items_id.}
\item{"prob_dependent"}{= Probability under local dependence. Same as input 'probs'.}
\item{"prob_marginal"}{= Probability of this response pattern under local independence.}
\item{"odds_ratio"}{= Compares prob_dependent and prob_marginal. Values far from '1' show strong dependence.}
\item{"odds_ratio_str"}{= odds_ratio as a fraction}
\item{"diff_prob"}{= gives the risk difference: prob_dependent minus prob_marginal.}
}
"kl_df" measures the KL divergence for each domain and class. The higher the KL divergence the greater the dependence.
\itemize{
\item{"class"}{= Which class does this apply to? NA indicates this term is aggregated across classes.}
\item{"items_id"}{= Unique identifier for the items in the domain. Matches dlcm$mcmc$domains$items_id.}
\item{"kl_divergence"}{= Calculates the Kullback-Leibler divergence. Assumes the dependent model is true. On log scale, gets the expected likelihood ratio under dependence versus independence. Zero indicates independence. The higher the number the greate the dependence.}
\item{"kl_maximum"}{= The KL divergence value under perfect dependence.}
\item{"kl_ratio"}{= The KL divergence scaled to 0 (independence) to 1 (perfect dependence).}
}
}
\description{
For each domain domain, compare the probabilities under local dependence (DLCM) versus under local independence (traditional LCM). Describes the amount of dependence captured by our DLCM model.
}
