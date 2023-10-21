#' @include utilities.r
NULL

##############
############## SUMMARY
##############


#' Get some summary statistics from fitted dependent LCM
#' @param dlcm list. Fitted dependent LCM
#' @param nwarmup integer. The first n iterations are considered warmup iterations
#' @return List with:
#' \itemize{
#' \item{"thetas_avg"}{= average response probabilities for different patterns}
#' \item{"thetas_avg_mode"}{= same information as thetas_avg, but just for the most common domain structure}
#' \item{"domain_items"}{= Which items are commonly grouped together in the same domain}
#' \item{"domain_items_all"}{= What are the most common domain structures? The column 'domains_merged' describes the domain structure as a string listing each item. Bars "|" separate domains, and commas "," separate items within a domain. For (partially) heterogeneous DLCMs, pluses "+" separate groups of classes (class2domain's) which may have different domain structures. Column 'n' gives the number of post-warmup iterations this domain structure appeared in.}
#' \item{"mode_domains"}{= Describes the single most common domain structure. Each row indicates one domain.}
#' \item{"domain_accept"}{= How often our metropolis step accepts its proposal. See dependentLCM_fit > domains_accept for details.}
#' \item{"classes"}{= For each subject this vector gives the most common class that observation belongs to.}
#' \item{"classes_cnts"}{= For each subject, in what number of iterations was this subject in each class?}
#' \item{"classes_pi"}{=The average prior probability of each class.}
#' \item{"dependence_intensity_dfs"}{=For each domain, how strong is the dependence? Compares response probabilities under local dependence vs local independence. See dependence_intensity().}
#' \item{"waic_summary"}{=See function 'dlcm.get_waic()'.}
#' \item{"waic_df"}{=See function 'dlcm.get_waic()'.}
#' }
#' @export
dlcm.summary <- function(dlcm, nwarmup=NULL, waic_method="agg") {
  
  if (is.null(nwarmup)) {
    nwarmup = 0
  }
  
  { # Summarize most common domain structures
    
    # Per domain
    domain_items <- (
      dlcm$mcmc$domains 
      %>% dplyr::filter(pattern_id==0, itr > nwarmup) 
      %>% dplyr::group_by(class2domain, items) 
      %>% dplyr::filter(class == min(class)) 
      %>% dplyr::summarize(nitems=max(nitems), n=dplyr::n(), .groups="keep") 
      %>% dplyr::arrange(-n)
    )
    
    # Per domains structure
    domain_items_all_raw <- (
      dlcm$mcmc$domains_merged
      %>% dplyr::group_by(itr) 
      # %>% dplyr::arrange(class2domain) # should already be ordered 
      %>% dplyr::summarize(domains_merged=paste0(domains_merged, collapse="+"))
    )
    domain_items_all <- (
      domain_items_all_raw
      %>% dplyr::filter(itr > nwarmup)
      %>% dplyr::group_by(domains_merged) 
      %>% dplyr::summarize(n=dplyr::n())
      %>% dplyr::mutate(perc=n/sum(n))
      %>% dplyr::arrange(-n)
    )
    
    # single most common domain structure
    first_mode_domain_itr <- unlist(domain_items_all_raw[which(domain_items_all_raw$domains_merged == unclass(domain_items_all[1,"domains_merged"]))[1],"itr"])
    mode_domains <- unique(dlcm$mcmc$domains %>% dplyr::filter(itr==first_mode_domain_itr, pattern_id==0) %>% .[,c("class", "class2domain", "items_id", "items")]) %>% dplyr::mutate(is_one=1)
    
    rm(domain_items_all_raw)
  }
  
  { # Response probabilities (thetas)
    
    # identify each domain and response pattern
    response_patterns <- dlcm$mcmc$response_patterns
    response_patterns$pattern_str <- sapply(response_patterns$pattern, paste0, collapse=", ")
    
    
    thetas_avg <- (
      dlcm$mcmc$domains
      %>% dplyr::filter(itr > nwarmup) 
      %>% dplyr::group_by(class, items_id, pattern_id) 
      %>% dplyr::summarize(
        items = dplyr::first(items)
        , n=dplyr::n()
        , prob=mean(prob)
        , .groups="drop")
      %>% dplyr::left_join(
        y=response_patterns[,c("items_id", "pattern_id", "pattern_str")]
        , by=c("items_id", "pattern_id")
      ) # lookup item_value
      %>% dplyr::left_join(
        y=mode_domains[,c("class", "items_id", "is_one")]
        , by=c("class", "items_id")
      ) # identify which domains are from the most common domain structure
      %>% dplyr::rename(is_mode=is_one, item_value=pattern_str)
      %>% .[, c("class", "items_id", "pattern_id", "items", "item_value", "is_mode", "n", "prob")]
    )
    
    dependence_intensity_dfs <- dependence_intensity(
      thetas_avg=thetas_avg
      , dlcm=dlcm
      , items_ids = unique(mode_domains$items_id)
    )


    thetas_avg_mode <- reshape2::dcast(
      data=thetas_avg %>% dplyr::filter(is_mode==TRUE)
      , formula = items_id + pattern_id + items + item_value ~ class
      , value.var = "prob"
    )
  }
  
  
  # Summarize Classes
  if (length(dlcm$mcmc$class_counts)>0) {
    classes_cnts <- dlcm$mcmc$class_counts # use existing summary
  } else if (dlcm$hparams$nclass>1) {
    # build summary
    itrs <- get_itrs_helper(nsaved=dlcm$hparams$save_itrs["classes"], nwarmup=nwarmup, nitr=dlcm$hparams$nitr)
    
    classes_cnts <- apply(
      dlcm$mcmc$classes[,itrs,drop=FALSE], 1
      , function(iclasses, class_levels) {
        return(table(factor(iclasses, levels=class_levels)))
      }
      , class_levels = seq(from=0, to=dlcm$hparams$nclass-1)
    )
    rownames(classes_cnts) <- paste0("class", rownames(classes_cnts))
    colnames(classes_cnts) <- paste0("obs", seq_len(ncol(classes_cnts)))
    
  } else {
    classes_cnts <- NA
  }
  
  if (!identical(classes_cnts, NA)) {
    classes <- apply(classes_cnts, 2, which.max) - 1
  } else {
    classes <- NA
  }
  
  
  waic <- dlcm.get_waic(
    dlcm
    , itrs=get_itrs_helper(nsaved=dlcm$mcmc$nitrLik, nwarmup=nwarmup, nitr=dlcm$hparams$nitr)
    , method=waic_method)
  names(waic) <- paste0("waic_", names(waic))
  
  itrs <- get_itrs_helper(nsaved=ncol(dlcm$mcmc$class_pi), nwarmup=nwarmup, nitr=dlcm$hparams$nitr)
  class_pi = rowMeans(dlcm$mcmc$class_pi[,itrs, drop=FALSE])
  
  return(c(
    list("thetas_avg"=thetas_avg, "domain_items"=domain_items, "domain_items_all"=domain_items_all, "classes"=classes, "class_pi"=class_pi, "thetas_avg_mode"=thetas_avg_mode, "mode_domains"=mode_domains, "first_mode_domain_itr"=first_mode_domain_itr, "classes_cnts"=classes_cnts, "dependence_intensity_dfs"=dependence_intensity_dfs)
    , waic
  ))
}

#' Generate a range of iterations which are both saved and occur after warmup.
#' @keywords internal
get_itrs_helper <- function(nsaved, nwarmup, nitr) {
  itrs <- seq(
    from=max(
      1 # take all saved itrs
      , nwarmup-(nitr-nsaved) # warmup itrs only
      )
    , to=nsaved # indexed based on saved location, not based on raw iteration
    )
  
  # itrs <- itrs + (nitr-nsaved) # convert to raw iteration, if desired
  
  return(itrs)
}


##############
############## WAIC
##############


#' Calculate likelihood and WAIC
#' @param dlcm Dependent Latent class model
#' @param method string from c("agg", "raw"). 
#' @param itrs integer vector. Which iterations should be include in calculation? Ignored if "agg" is chosen.
#' @return List with a 'summary' and 'df' value. The 'summary' is a vector with the following summary statitics:
#' \itemize{
#' \item{"nparams_avg"}{= The average number of model parameters per iteration.}
#' \item{"logLik_avg"}{= The average logLikelihood per iteration,}
#' \item{"lppd"}{= The LPPD (Log PointWise Predictive Density)}
#' \item{"waic_nparams1"}{= The effective number of parameters of the model according to WAIC method 1. Describes model complexity and tendency to overfit.}
#' \item{"waic_nparams2"}{= The effective number of parameters of the model according to WAIC method 2. Describes model complexity and tendency to overfit.}
#' \item{"waic1"}{= The WAIC (Watanabe-Akaike Information Criterion) using method 1. Approximates the overall (out of sample) goodness of fit of the model. Equals '-2 lppd + 2 waic_nparams.'}
#' \item{"waic2"}{= The WAIC (Watanabe-Akaike Information Criterion) using method 2. Approximates the overall (out of sample) goodness of fit of the model. Equals '-2 lppd + 2 waic_nparams.}
#' }
#' For information on wAIC methods 1/2 see Gelman 2014 (DOI https://doi.org/10.1007/s11222-013-9416-2). The 'df' value has a row for each iteration and reports the logLikelihood and raw number of parameters.
#' @export
dlcm.get_waic <- function(dlcm, method="agg", itrs=NULL) {
  if (method=="agg") {
    return(dlcm.get_waic_fromagg(dlcm=dlcm))
  } else if (method=="raw") {
    return(dlcm.get_waic_fromraw(dlcm=dlcm, itrs=itrs))
  }
}

#' Calculate likelihood and WAIC from aggregated information. See dlcm.get_waic() for details.
#' Warning: Only works if class_loglik_collapsed is saved, as in: dependentLCM_fit(..., save_itrs=c(class_loglik_collapsed=Inf, ...))
#' @keywords internal
dlcm.get_waic_fromagg <- function(dlcm) {
  
  itrs <- as.numeric(substring(names(dlcm$mcmc$itrLogLik), 4))
  likelihoods <- data.frame(
    itr = itrs
    , logLik = dlcm$mcmc$itrLogLik
    , nparameters = 0
  )
  
  theta_parameters <- (
    dlcm$mcmc$domains 
    %>% dplyr::filter(itr %in% itrs) 
    %>% dplyr::group_by(itr) 
    %>% dplyr::summarize(nparams = dplyr::n() - length(unique(domain)), .groups="keep"))
  likelihoods$nparameters <- (
    unlist(theta_parameters[match(likelihoods$itr, theta_parameters$itr), "nparams"])
    + dlcm$hparams$nclass-1 # class pis
    + length(unique(dlcm$hparams$class2domain)) * (dlcm$hparams$nitems-1) # domains
  )
  
  summary <- c(
    nparams_avg = mean(likelihoods$nparameters)
    , logLik_avg = mean(dlcm$mcmc$itrLogLik)
    , lppd = sum(log(dlcm$mcmc$obsLik / dlcm$mcmc$nitrLik))
    , waic_nparams1 = 2*sum(log(dlcm$mcmc$obsLik / dlcm$mcmc$nitrLik) - dlcm$mcmc$obsLogLik / dlcm$mcmc$nitrLik)
    , waic_nparams2 = sum(dlcm$mcmc$obsLogLik2 / (dlcm$mcmc$nitrLik-1) - dlcm$mcmc$nitrLik/(dlcm$mcmc$nitrLik-1)*(dlcm$mcmc$obsLogLik / dlcm$mcmc$nitrLik)^2)
  )
  summary["waic1"] <- -2 * (summary["lppd"] - summary["waic_nparams1"])
  summary["waic2"] <- -2 * (summary["lppd"] - summary["waic_nparams2"])
  summary["aic"] <- -2*summary["logLik_avg"] + 2*summary["nparams_avg"]
  summary["nitrs_used"] <- dlcm$mcmc$nitrLik
  
  return(list(
    summary=summary
    , df=likelihoods
  ))
}

#' Calculate likelihood and WAIC from raw class_loglik. See dlcm.get_waic() for details.
#' Warning: Only works if class_loglik is saved, as in: dependentLCM_fit(..., save_itrs=c(class_loglik=Inf, ...))
#' @keywords internal
dlcm.get_waic_fromraw <- function(dlcm, itrs=NULL) {
  
  if (is.null(itrs)) {
    # default in all iterations
    itrs_str <- seq_len(dlcm$hparams$nitr)
    itrs <- as.numeric(substring(names(dlcm$mcmc$itrLogLik), 4))
  } else {
    itrs_str <- paste0("itr", itrs)
  }
  
  
  #
  # Initialize
  #
  
  nitr <- length(itrs)
  likelihoods <- data.frame(
    matrix(0, nrow=nitr, ncol=3)
  )
  colnames(likelihoods) <- c("itr", "logLik", "nparameters")
  likelihoods$itr <- itrs
  
  summary <- c()
  
  #
  # Parameters, raw
  #
  
  theta_parameters <- (
    dlcm$mcmc$domains 
    %>% dplyr::filter(itr %in% itrs) 
    %>% dplyr::group_by(itr) 
    %>% dplyr::summarize(nparams = dplyr::n() - length(unique(domain)), .groups="keep"))
  likelihoods$nparameters <- (
    unlist(theta_parameters[match(likelihoods$itr, theta_parameters$itr), "nparams"])
    + dlcm$hparams$nclass-1 # class pis
    + length(unique(dlcm$hparams$class2domain)) * (dlcm$hparams$nitems-1) # domains
  )
  summary["nparams_avg"] <- mean(likelihoods$nparameters)
  
  #
  # Likelihood and WAIC
  #
  
  # Calculate likelihood for each observation in each iteration
  obsLogLiks <- sweep(dlcm$mcmc$class_loglik[,,itrs_str, drop=FALSE]
                      , c(1,3) # Include pi. Repeat for each observation (2)
                      , log(dlcm$mcmc$class_pi[,itrs_str, drop=FALSE]), "+")
  obsLogLiks <- apply(obsLogLiks, c(2,3), expSumLog)
  
  likelihoods$logLik <- colSums(obsLogLiks, 1) # sum across observations
  
  # Same as: LaplacesDemon::WAIC
  summary["logLik_avg"] <- mean(likelihoods$logLik)
  summary["lppd"] <- sum(apply(obsLogLiks, 1, expSumLog) - log(dim(obsLogLiks)[2])) 
  summary["waic_nparams1"] <- 2 * sum(apply(obsLogLiks, 1, expSumLog) - log(dim(obsLogLiks)[2])
                                      - apply(obsLogLiks, 1, mean))
  summary["waic_nparams2"] <- sum(apply(obsLogLiks, 1, var))
  summary["waic1"] <- -2 * (summary["lppd"] - summary["waic_nparams1"])
  summary["waic2"] <- -2 * (summary["lppd"] - summary["waic_nparams2"])
  summary["aic"] <- -2*summary["logLik_avg"] + 2*summary["nparams_avg"]
  summary["nitrs_used"] <- nrow(likelihoods)
  
  return(list(
    summary=summary
    , df=likelihoods
  ))
}


##############
############## Response Probabilities
##############


#' Calculate the probabilities of different patterns
#' @param items integer vector. Which items should we include? (indexed starting at 1)
#' @param this_sim a mcmc simulation. Output from dependentLCM_fit()
#' @param itrs integer vector. Which iterations from this_sim should we include?
#' @param merge_itrs TRUE if you want average probabilities across simulations. FALSE if you want to examine each iteration individually
#' @param classes Which classes are we evaluating? By default all classes
#' @description Useful when getting probilities which cross domains.
#' @returns Probabilities for each pattern of items
#' @export
theta_item_probs <- function(items, this_sim, itrs=NULL, merge_itrs=TRUE, classes=NULL) {
  
  if (is.null(itrs)) {
    itrs <- seq_len(this_sim$hparams$nitr)
  }
  
  class_filter <- TRUE
  if (!is.null(classes)
      & !identical(classes, seq_len(this_sim$hparams$nclass)-1)
  ) {
    this_sim$mcmc$class_pi[-(classes+1),] <- 0
    this_sim$mcmc$class_pi <- sweep(
      this_sim$mcmc$class_pi
      , 2
      , colSums(this_sim$mcmc$class_pi)
      , '/'
    )
    class_filter = (this_sim$mcmc$domains$class %in% classes)
  }
  
  
  afilter <- which(
    (this_sim$mcmc$domains$itr %in% itrs)
    & class_filter
  )
  nitr <- length(itrs)
  nitems <- length(items)
  items_colnames <- colnames(this_sim$mcmc$domains)[this_sim$hparams$domain_item_cols[items]]
  
  # Reduce thetas down to relevant patterns only (merging redundant patterns as necessary)
  thetas_agg <- aggregate(
    formula(paste0("prob ~ ", paste(c(items_colnames, "class", "itr"), collapse=" + ")))
    , this_sim$mcmc$domains[afilter, ]
    , sum
  )
  thetas_agg <- thetas_agg[!(rowSums(thetas_agg[, seq_len(nitems), drop=FALSE]) == -nitems), ] # Remove unrelated rows
  thetas_agg$prob_log <- log(thetas_agg$prob)
  thetas_agg$class_pi <- this_sim$mcmc$class_pi[cbind(thetas_agg$class+1, thetas_agg$itr)]
  
  all_patterns <- expand.grid(
    lapply(this_sim$hparams$item_nlevels[items], function(n) 0:(n-1))
  )
  colnames(all_patterns) <- items_colnames
  
  prob_helper <- function(ipattern) {
    # Purpose: Calculate the probability for this pattern
    # Implict: thetas_agg, this_sim, nitr, nitems
    imatch_pattern <- rowSums(
      sweep(thetas_agg[,seq_len(nitems), drop=FALSE], 2, ipattern, "==") # Matches pattern
      | thetas_agg[,seq_len(nitems), drop=FALSE] < 0 # Or is undefined
    ) == nitems
    iprobs <- (
      thetas_agg 
      %>% dplyr::filter(imatch_pattern) 
      %>% dplyr::group_by(itr, class) 
      %>% dplyr::summarize(prob=exp(sum(prob_log)), class_pi=dplyr::first(class_pi), .groups="keep"))
    iprobs$prob_pi <- iprobs$prob * iprobs$class_pi
    return(iprobs)
  }
  
  if (merge_itrs==TRUE) {
    # Get average across iterations
    out <- all_patterns
    out$prob <- apply(all_patterns, 1
                      , FUN = function(ipattern) {
                        iprobs <- prob_helper(ipattern)
                        return(sum(iprobs$prob_pi) / nitr)
                      })
  } else {
    # Return probabilities for each individual iteration
    out <- list(
      all_patterns=all_patterns
      , probs = t(apply(all_patterns, 1, FUN = function(ipattern) {
        iprobs <- prob_helper(ipattern)
        out <- (iprobs %>% dplyr::group_by(itr) %>% dplyr::summarize(probs=sum(prob_pi), .groups="keep")) %>% .$probs
        return(out)
      }))
    )
  }
  
  return(out)
}


#' For each domain domain, compare the probabilities under local dependence (DLCM) versus under local independence (traditional LCM). Describes the amount of dependence captured by our DLCM model.
#' @param thetas_avg dataframe. Describes average response probabilities for different domains. As returned from dlcm.summary()
#' @param dlcm list. Fitted dependent LCM
#' @param items_ids integerVector Which domains do we want to process?
#' @return DataFrames with:
#' \itemize{
#' \item{"prob_dependent"}{= Probability under local dependence. Same as input 'probs'.}
#' \item{"prob_marginal"}{= Probability of this response pattern under local independence.}
#' \item{"odds_ratio"}{= Compares prob_dependent and prob_marginal. Values far from '1' show strong dependence.}
#' \item{"odds_ratio_str"}{= odds_ratio as a fraction}
#' }
#' @export
dependence_intensity <- function(thetas_avg, dlcm, items_ids=NULL) {
  
  response_patterns <- dlcm$mcmc$response_patterns
  response_patterns$nitems <- sapply(response_patterns$items, length) 
  
  if (identical(items_ids, NULL)) {
    items_ids=unique(thetas_avg$items_id)
  }
  
  thetas_marginal_ratio <- (
    thetas_avg
    %>% dplyr::left_join(
      y=response_patterns
      , by=c("items_id", "pattern_id")
    )
    %>% dplyr::filter(nitems>1, items_id %in% items_ids)
    %>% dplyr::group_by(class, items_id)
    %>% dplyr::summarize(
      pattern_id=list(pattern_id)
      , patterns=list(do.call(rbind, pattern))
      , probs=list(prob), .groups="drop")
  )
  thetas_marginal_ratio$ratios <- mapply(
    function(values, probs, ...) {
      out <- dependence_intensity_one(values, probs)
      ids <- do.call(cbind.data.frame, list(...))
      out <- cbind.data.frame(ids, out)
      return(out)
    }
    , values=thetas_marginal_ratio$patterns
    , probs=thetas_marginal_ratio$probs
    # identifiers:
    , class=thetas_marginal_ratio$class
    , items_id=thetas_marginal_ratio$items_id
    , pattern_id=thetas_marginal_ratio$pattern_id
    , SIMPLIFY = FALSE
  )
  
  out <- split(
    thetas_marginal_ratio
    , f=thetas_marginal_ratio$items_id
  )
  out <- lapply(
    out
    , function(x) do.call(rbind.data.frame, x$ratios)
  )
  
  return(out)
}


#' Takes a single domain. Compare the probabilites under the given probabilities versus under local independence.
#' @param values characterMatrix One row for each response pattern to this domain. One column per item. Cell contains the value of a given item under the given response pattern.
#' @param probs numericVector. For a single class, probability of returning this response pattern.
#' @inherit dependence_intensity return
#' @keywords internal
dependence_intensity_one <- function(values, probs) {
  
  class(values) <- "character"
  
  item_prob_marginal <- lapply(
    seq_len(ncol(values))
    , function(iitem) {
      iprobs <- data.frame(value=values[,iitem], prob=probs)
      iprob_marginal <- iprobs %>% dplyr::group_by(value) %>% dplyr::summarize(prob=sum(prob))
      iprob_marginal <- setNames(iprob_marginal$prob, iprob_marginal$value)
      return(iprob_marginal)
    }
  )
  
  prob_marginal <- mapply(
    function(item_values, item_probs) {
      item_probs[item_values]
    }
    , item_values = unclass(as.data.frame(values))
    , item_probs = item_prob_marginal
  )
  prob_marginal <- apply(prob_marginal, 1, prod)
  
  odds_ratio <- (
    ( probs / (1-probs) )
    / (prob_marginal / (1-prob_marginal) )
  )
  
  odds_ratio_str <- sapply(
    odds_ratio
    , function(x) {
      if (x >= 1) {
        out_str <- format(x, digits=1, nsmall=1)
      } else {
        out_str <- paste0("1/", format(1/x, digits=1, nsmall=1), collapse="")
      }
    }
  )
  
  
  df_out <- data.frame(
    values
    , prob_dependent=probs
    , prob_marginal=prob_marginal
    , odds_ratio=odds_ratio
    , odds_ratio_str=odds_ratio_str
  )
  
  return(df_out)
}


##############
############## MISC
##############


#' What is the prior probability of this choice of domains?
#' Only correct up to a normalizing constant. Ignores identifiability restrictions.
#' @param x IntegerVector. The number of items in each domain. Ok to omit 0's
#' @param ndomains Integer. The total number of domains (including empty domains)
#' @param specific_items Boolean. 
#' If FALSE, we look for any domain which produces domains of this size regardless of what specific items they contain.
#' IF TRUE, we fix which items are in which domain, and calculate the probability of grouping these specific items together.
#' @param log Boolean. If TRUE give probability in log scale.
#' @param domain_theta_prior_type String. What domain prior are we using? One of c("permissive", "restrictive")
#' @param x_npatterns IntegerVector. For domain_theta_prior_type="restrictive". For each domain x, how many possible response patterns are there to that domain?
#' @export
ldomain_prior <- function(x, ndomains, specific_items=FALSE, log=TRUE, domain_theta_prior_type="permissive", x_npatterns=NULL) {
  
  lpattern_adjustment = 0 # for "restrictive"
  if (domain_theta_prior_type=="restrictive") {
    if (length(x) != length(x_npatterns)) {
      stop("ldomain_prior: For domain_theta_prior_type==restrictive we require x_npatterns to be set.")
    }
    lpattern_adjustment = sum(lfactorial(x_npatterns - 1))
  }
  
  if (specific_items==FALSE) {
    lprob <- (
      lfactorial(sum(x)) - sum(sapply(x, lfactorial)) # unique groups
      - sum(sapply(table(x[x>0]), lfactorial)) # non-unique groups
      + lfactorial(ndomains) - lfactorial(ndomains-sum(x > 0)) # permuting groups
      - sum(x) * log(ndomains) # denominator
      - lpattern_adjustment
    )
  } else { # specific_items==TRUE
    lprob <- (
      lfactorial(ndomains) - lfactorial(ndomains-sum(x > 0)) # permuting groups
      - sum(x) * log(ndomains) # denominator
      - lpattern_adjustment
    )
  }
  if (log==FALSE) {
    return(exp(lprob))
  }
  return(lprob)
}


#' For each iteration calculate the probability that a given observation is in each class.
#' Warning: Only works if class_loglik is saved, as in: dependentLCM_fit(..., save_itrs=c(class_loglik=Inf, ...))
#' @param dlcm Dependent latent class model outptut from dependentLCM_fit
#' @export
get_class_probs <- function(dlcm) {
  obsLogLiks <- sweep(dlcm$mcmc$class_loglik[,,, drop=FALSE]
                      , c(1,3) # Include pi. Repeat for each observation (2)
                      , log(dlcm$mcmc$class_pi[,, drop=FALSE]), "+")
  obsLogLiks_sum <- apply(obsLogLiks, c(2,3), expSumLog)
  obsLogLiks <- sweep(obsLogLiks
                      , c(2,3)
                      , obsLogLiks_sum
                      , "-"
  )
  return(exp(obsLogLiks))
}

#' For each iteration, calculate the prior probabilities of each parameter and conditional probability of our responses (marginalized over class)
#' Assumes domain prior of "permissive".
#' Might require save_itrs=c(all=0, classes=Inf) in DLCM model.
#' @inheritParams get_jointLikelihood_obs
#' @export
get_jointLikelihood <- function(dlcm, method="itrLogLik") {
  joint_loglikelihood_prior <- get_jointLikelihood_priors(dlcm)
  joint_loglikelihood_obs <- get_jointLikelihood_obs(dlcm, method)
  
  joint_loglikelihood <- dplyr::left_join(
    x=joint_loglikelihood_prior
    , y=joint_loglikelihood_obs
    , by="itr"
  )
  joint_loglikelihood$ltotal <- joint_loglikelihood$lpriors + joint_loglikelihood$obs_lprob
  return(joint_loglikelihood)
}


#' Calculates probability of each response collapsed on classes, then aggregated for each iteration
#' @param dlcm Dependent latent class model outptut from dependentLCM_fit
#' @param method Defines what variable from dlcm$mcmc is used to calculate resposne probabilities. One of c("itrLogLik", "class_loglik", "").
#' @keywords internal
get_jointLikelihood_obs <- function(dlcm, method) {
  
  if (method=="class_loglik") {
    
    # iterations with class_loglik set
    itrs <- seq(to=dlcm$hparams$nitr, length=dlcm$hparams$save_itrs["class_loglik"])
    
    obsLogLiks <- sweep(dlcm$mcmc$class_loglik[,,, drop=FALSE]
                        , c(1,3) # Include pi. Repeat for each observation (2)
                        , log(dlcm$mcmc$class_pi[,itrs , drop=FALSE]), "+")
    obsLogLiks <- apply(obsLogLiks, c(2,3), expSumLog)
    obsLogLiks_agg <- apply(obsLogLiks, 2, sum)
    
  } else if (method=="itrLogLik") {
    # probability of each response conditional on classes
    
    itrs <- seq(to=dlcm$hparams$nitr, length=dlcm$hparams$save_itrs["agg_loglik"])
    obsLogLiks_agg <- dlcm$mcmc$itrLogLik
  } else {
    # keep blank
    itrs = c(NA)
    obsLogLiks_agg = c(NA)
  }
  
  obsLogLiks_df <- data.frame(
    itr=itrs
    , obs_lprob=obsLogLiks_agg
  )
  return(obsLogLiks_df)
}

#' For each iteration, calculate the prior probabilities of each parameter
#' @inheritParams get_jointLikelihood_obs
#' @keywords internal
get_jointLikelihood_priors <- function(dlcm) {
  
  duplicated_domains <- duplicated(dlcm$mcmc$domains[, c("itr", "domain", "class2domain")]) # Want only one row per domain
  domain_prior <- (
    dlcm$mcmc$domains[!duplicated_domains,!(colnames(dlcm$mcmc$domains)%in% c("items"))]
    %>% dplyr::group_by(itr, class2domain)
    %>% dplyr::summarize(
      nitems_list=list(nitems)
      , domain_lprior=ldomain_prior(x=nitems, ndomains=dlcm$hparams$ndomains, specific_items=TRUE, log=TRUE)
      , .groups="keep")
    %>% dplyr::group_by(itr)
    %>% dplyr::summarize(domain_lprior=sum(domain_lprior), .groups="drop")
  )
  
  pi_prior <- log(apply(dlcm$mcmc$class_pi, 2, gtools::ddirichlet, alpha=dlcm$hparams$classPi_alpha*rep(1, dlcm$hparams$nclass)))
  
  thetas_prior <- (
    dlcm$mcmc$domains
    %>% dplyr::group_by(itr, domain)
    %>% dplyr::summarize(
      theta_lprior=LaplacesDemon::ddirichlet(x=prob, alpha=dlcm$hparams$theta_alpha*rep(1, dplyr::n()), log=TRUE)
      , .groups="drop"
    )
  )
  thetas_prior_agg <- (
    thetas_prior
    %>% dplyr::group_by(itr)
    %>% dplyr::summarize(theta_lprior=sum(theta_lprior))
  )
  
  joint_loglikelihood <- data.frame(
    itr = domain_prior$itr
    , domain_lprior=domain_prior$domain_lprior
    , pi_lprior=pi_prior
    , thetas_lprior = thetas_prior_agg$theta_lprior
  )
  joint_loglikelihood$lpriors_total <- rowSums(joint_loglikelihood[,c("domain_lprior", "pi_lprior", "thetas_lprior")])
  domain_strs <- (
    dlcm$mcmc$domains_merged
    %>% dplyr::group_by(itr) 
    # %>% dplyr::arrange(class2domain) # should already be ordered 
    %>% dplyr::summarize(domains_merged=paste0(domains_merged, collapse="+"))
    %>% .[,"domains_merged"]
  )
  joint_loglikelihood$domain <- domain_strs$domains_merged
  
  return(joint_loglikelihood)
}