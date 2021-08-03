# library(klaR) # for kmodes
# library(Rcpp)
# library(abind)
# library(dplyr)

#' dependentLCM: Dependent Latent Class Model
#'
#' @docType package
#' @name dependentLCM
#' @useDynLib dependentLCM, .registration = TRUE
NULL

NCLASS = 2
CLASSPI_ALPHA = 2
DOMAIN_ALPHA_RATIO = 2
THETA_ALPHA = 2
DOMAIN_PROPOSAL_EMPTY = 0.3
DOMAIN_PROPOSAL_SWAP = 0.3

#' Fits a bayesian dependent LCM model
#' @param nitr integer. Number of iterations to run the bayes MCMC
#' @inheritParams getStart_hparams
#' @inheritParams getStart_bayes_params
#' @export
dependentLCM_fit <- function(
  nitr
  # Data
  , df=NULL, mat=NULL
  # Hyperparameters
  ,nclass=NCLASS, ndomains=NULL, class2domain=NULL, classPi_alpha=CLASSPI_ALPHA, domain_alpha=NULL, domain_maxitems=NULL, theta_alpha=THETA_ALPHA, domain_proposal_empty=DOMAIN_PROPOSAL_EMPTY, domain_proposal_swap=DOMAIN_PROPOSAL_SWAP, domain_nproposals=NULL
  # Bayes parameters
  , class_pi = NULL, classes = NULL, thetas = NULL
  # Misc
  , class_init_method = "kmodes") {
  
  datetimes <- c(Sys.time())
  
  #
  # Set starting values
  #

  if (is.null(mat)) {
    mat <- getStart_matrix(df)
  }

  hparams <- getStart_hparams(
    df=mat
    # Hyperparameters
    ,nclass=nclass, ndomains=ndomains, class2domain=class2domain, classPi_alpha=classPi_alpha, domain_alpha=domain_alpha, domain_maxitems=domain_maxitems, theta_alpha=theta_alpha, domain_proposal_empty=domain_proposal_empty, domain_proposal_swap=domain_proposal_swap, domain_nproposals=domain_nproposals
  )

  bayesparams <- getStart_bayes_params(
    mat=mat, hparams=hparams
    , class_pi = class_pi, classes = classes, thetas = thetas
    # Misc
    , class_init_method = class_init_method
  )
  
  all_params <- list(
    mat = mat
    , hparams = hparams
    , bayesparams = bayesparams
  )

  if (check_params(all_params)) {
    stop()
  }

  #
  # Run MCMC
  #
  
  datetimes <- c(datetimes, Sys.time())
  dlcm = dependentLCM_fit_cpp(x=(all_params$mat),hparams_list=all_params$hparams, params_list=all_params$bayesparams, nitr=nitr)
  datetimes <- c(datetimes, Sys.time())
  
  #
  # Post Processing
  #

  # Clean up output
  dlcm$thetas_id <- do.call(cbind, dlcm$thetas_id)
  rownames(dlcm$thetas_id) <- c("itr", "class", "domain", "pattern_id")
  dlcm$thetas_patterns <- do.call(cbind, dlcm$thetas_patterns)
  rownames(dlcm$thetas_patterns) <- paste0("item_", 0:(nrow(dlcm$thetas_patterns)-1))
  dlcm$thetas_probs <- unlist(dlcm$thetas_probs)

  # Supplemental
  thetas_items <- apply(
    dlcm$thetas_patterns > -1
    , 2, function(x) paste0(c("", which(x)-1, ""), collapse=",")
  )
  thetas_nitems <- colSums(dlcm$thetas_patterns > -1)
  thetas_class2domain <- all_params$hparams$class2domain[dlcm$thetas_id["class",,drop=TRUE]+1]
  thetas_domain_id <- as.data.frame(t(rbind(dlcm$thetas_id, thetas_class2domain)))
  thetas_domain_id <- do.call(paste, thetas_domain_id) # paste rows together, faster than apply

  # Merge thetas attributes
  dlcm$thetas <- data.frame(
    t(dlcm$thetas_id)
    , class2domain = thetas_class2domain
    , prob=dlcm$thetas_probs
    , nitems = thetas_nitems
    , items = thetas_items
    , t(dlcm$thetas_patterns)
  )
  all_params$hparams$theta_item_cols <- grep("^item_[0-9]+", colnames(dlcm$thetas))
  dlcm$thetas_accept <- do.call(function(...) abind::abind(..., along=3), dlcm$thetas_accept)

  
  dlcm$domains_merged <- dplyr::summarize(dplyr::group_by(dplyr::filter(dlcm$thetas, pattern_id==0), itr, class2domain), domains_merged=paste(sort(unique(items)), collapse="|"), .groups="keep") # dlcm$domains_merged <- dlcm$thetas %>% dplyr::filter(pattern_id==0) %>% dplyr::group_by(itr, class2domain) %>% dplyr::summarize(domains_merged=paste(sort(unique(items)), collapse="|"), .groups="keep") 

  dlcm$thetas_id <- NULL
  dlcm$thetas_patterns <- NULL
  dlcm$thetas_probs <- NULL
  dlcm$nclass2domain <- NULL
  
  datetimes <- c(datetimes, Sys.time())
  dlcm$runtimes <- c(diff(datetimes), total=tail(datetimes,1)-datetimes[1])
  names(dlcm$runtimes)[1:3] <- c("pre", "mcmc", "post")
  
  return(list(hparams=all_params$hparams, mcmc=dlcm))
}


#' Generate hyperparameters
#' Generate a list of hyperparameters, adding default values when necessary
#' @param df dataframe. Raw data you are fitting with all assumed to be factors
#' @param nitems integer. Number of columns of df.
#' @param nclass integer. Number of subject latent classes
#' @param ndomains integer. Number of item domains
#' @param class2domain integer vector of length nclass. Classes with same value have same domains.
#' @param classPi_alpha numeric vector. Bayes hyperparameter giving prior for Pi.
#' @param domain_alpha numeric. Bayes hyperparameter giving prior for domains.
#' @param domain_maxitems iinteger. Maximum number of items which can be in a doamin.
#' @param theta_alpha numeric. Bayes hyperparemter giving the prior for theta/probabilties.
#' @param domain_proposal_empty numeric. Sets how often the domain metropolis proposal function pairs a nonempty domain with an empty domain.
#' @param domain_proposal_swap numeric. Sets how often the domain metropolis proposal function swaps items between domains.
#' @param domain_nproposals. Sets how many times the domain metropolis propsal function is called each iteration
#' @keywords internal
getStart_hparams <- function(
  df=NULL, nitems=NULL
  # Hyperparameters
  ,nclass=NCLASS, ndomains=NULL, class2domain=NULL, classPi_alpha=CLASSPI_ALPHA, domain_alpha=NULL, domain_maxitems=NULL, theta_alpha=THETA_ALPHA, domain_proposal_empty=DOMAIN_PROPOSAL_EMPTY, domain_proposal_swap=DOMAIN_PROPOSAL_SWAP, domain_nproposals=NULL
) {
  # Purpose: Add default hyperparameters

  if (is.null(nitems)) {
    nitems <- dim(df)[2]
  }
  item_nlevels <- apply(df, 2, max)+1

  # ndomains
  if (is.null(ndomains)){
    ndomains <- 100*nitems # Many more domains than items
  }

  # class2domain
  if (is.null(class2domain)) {
    class2domain <- rep(0, nclass) # all classes use same domains
  }
  class2domain <- as.integer(factor(class2domain))-1 # Make sequential starting at 0
  nclass2domain <- length(unique(class2domain))

  # classPi_alpha
  if (length(classPi_alpha)==1) {
    classPi_alpha <- rep(classPi_alpha, nclass)
  }

  # domain_alpha
  if (is.null(domain_alpha)) {
    domain_alpha <- nitems * DOMAIN_ALPHA_RATIO
  }

  # domain_maxitems
  if (is.null(domain_maxitems)) {
    domain_maxitems <- nitems # No restrictions
  }

  if (is.null(domain_nproposals)) {
    domain_nproposals <- nitems
  }

  # theta_alpha, no action taken

  return(list(
    nclass=nclass, ndomains=ndomains, class2domain=class2domain, classPi_alpha=classPi_alpha, domain_alpha=domain_alpha, domain_maxitems=domain_maxitems, theta_alpha=theta_alpha, nitems = nitems, item_nlevels = item_nlevels, nclass2domain = nclass2domain, domain_proposal_empty=domain_proposal_empty, domain_proposal_swap=domain_proposal_swap, domain_nproposals=domain_nproposals
  ))
}

#' Convert dataframe of factors to matrix of integers (indexed starting at 0)
#' @keywords internal
getStart_matrix <- function(df) {
  # Purpose: Convert to factor in integer form (0:k-1) and remove NA rows

  df <- as.data.frame(df)
  for (i in 1:ncol(df)) {
    df[, i] <- as.integer(factor(df[, i])) - 1
  }
  df <- as.matrix(df)
  df <- na.omit(df)
  mode(df) <- "integer"

  return(df)
}

#' Generates initial values for bayes parameters
#' @param mat matrix. The data you wish to analyze
#' @param hparams list. List of hyperparameters from getStart_hparams.
#' @param class_pi numeric vector size nclass. Initial condition for bayes parameter pi.
#' @param classes integer vector size nrow(df). Initial condition for subject classes.
#' @param thetas list. Initial values for thetas/probabilites.
#' @param class_init_method string. Decides how 'classes' is defaulted if NULL. One of "kmodes" or "random"
#' @keywords internal
getStart_bayes_params <- function(
  mat, hparams
  , class_pi = NULL, classes = NULL, thetas = NULL
  # Misc
  , class_init_method = "kmodes"
) {

  # classes
  if (is.null(classes)) {
    if (class_init_method=="kmodes") {
      classes <- getStart_class_kmodes(mat, hparams)
    }
    if (class_init_method=="random") {
      classes <- getStart_class_random(mat, hparams)
    }
  }

  # class_pi
  class_pi <- table(classes)+hparams$classPi_alpha
  class_pi <- class_pi / sum(class_pi)

  # thetas
  if (is.null(thetas)) {
    thetas <- getStart_thetas(mat, classes, hparams)
  } # Future maybe allow other types of inputs

  return(list(
    class_pi = class_pi, classes = classes, thetas = thetas
  ))
}


#' Choose class of each observation entirely at random
#' @param mat matrix. Raw data.
#' @param hparams list. List of hyperparameters
#' @keywords internal
getStart_class_random <- function(mat, hparams) {
  class_pi <- hparams$classPi_alpha / sum(hparams$classPi_alpha)
  nobs <- dim(mat)[1]
  classes <- sample(0:(hparams$nclass-1), size=nobs, replace=TRUE, prob=class_pi)

  # Ensure every class shows up atleast once:
  ids <- sample.int(nobs, hparams$nclass)
  classes[ids] <- 0:(hparams$nclass-1)

  return(classes)
}

#' Choose class of each observation using kmodes
#' @param mat matrix. Raw data.
#' @param hparams list. List of hyperparameters
#' @keywords internal
getStart_class_kmodes <- function(mat, hparams, iter.max=2, ...) {
  return(klaR::kmodes(mat, hparams$nclass, iter.max=iter.max, ...)$cluster-1)
}

#' Choose starting theta values.
#' @param mat matrix. Raw data.
#' @param classes integer vector. The class of each observation.
#' @param hparams list. List of hyperparameters
#' @keywords internal
getStart_thetas <- function(mat, classes, hparams) {
  # FUTURE: Maybe switch to internal C gibbs here? Keep this for now for validation
  thetas_list = list()
  for (iclass in 0:(hparams$nclass-1)) {
    thetas_list[[iclass+1]] <- lapply(
      1:hparams$nitems
      , function(icol) {
        ix = mat[classes==iclass, icol]
        ithetas = as.vector((table(ix) + hparams$theta_alpha))
        ithetas <- ithetas/sum(ithetas)
        iout = list(
          items=c(icol-1)
          # , item_nlevels = hparams$item_nlevels[icol] # now inferred
          , thetas = ithetas
        )
        return(iout)
      }
    )
  }
  return(thetas_list)
}

#' Check parameters for internal consistency
#' @param all_params list. List of parameters and hyperparameters
#' @keywords internal
check_params <- function(all_params) {
  check_out = FALSE

  if (all_params$hparams$nclass != length(all_params$hparams$class2domain)) {
    warning("nclass~class2domain mismatch")
    check_out = TRUE
  }

  names(all_params$bayesparams)

  return(check_out)
}


##############
############## UTILITIES
##############


#' Convert numeric vector x to percentages.
#' @keywords internal
get_perc <- function(x) {
  x / sum(x)
}

#' Calculate mode of x
#' @keywords internal
getmode <- function(x) {
  xcounts <- table(x)
  mode_name <- names(xcounts)[which.max(xcounts)]
  return(mode_name)
}

#' Reset all sinks. Sink is used to to write to file.
#' @keywords internal
sink.reset <- function(){
  for(i in seq_len(sink.number())){
    sink(NULL)
  }
}

#' Wrapper for dplyr::`%>%` 
#' @keywords internal
`%>%` <- function(...) {
  dplyr::`%>%`(...)
}


##############
############## SUMMARY
##############


#' Get some summary statistics from fitted dependent LCM
#' @param dlcm list. Fitted dependent LCM
#' @param nwarmup integer. The first n iterations are considered warmup iterations
#' @return List with:
#' thetas_avg = average response probabilities for different patterns
#' domain_items = Which items are commonly grouped together in the same domain
#' domain_nitems = How often the domains are of size k.
#' domain_accept = How often our metropolis step accepts its proposal
#' @export
dlcm.summary <- function(dlcm, nwarmup=1000) {

  dlcm$thetas$item_value <- apply(
    out$thetas[, out$theta_item_cols]
    , 1, function(x) paste0(x[x>-1], collapse=", ")
  )
  thetas_avg <- dlcm$thetas %>% dplyr::filter(itr > nwarmup) %>% dplyr::group_by(class, items, item_value) %>% dplyr::summarize(n=dplyr::n(), prob=mean(prob))

  domain_items <- rbind(
    dlcm$domains_merged %>% dplyr::group_by(class2domain, items=domains_merged) %>% dplyr::summarize(nitems=NA, n=dplyr::n(), all_items=TRUE) %>% dplyr::arrange(-n)
    , dlcm$thetas %>% dplyr::filter(pattern_id==0, itr > nwarmup) %>% dplyr::group_by(class2domain, items) %>% dplyr::summarize(nitems=max(nitems), n=dplyr::n(), all_items=FALSE) %>% dplyr::arrange(-n)
  )
  domain_nitems <- table((dlcm$thetas %>% dplyr::filter(pattern_id==0, itr > nwarmup))$nitems)

  # domain_accept <- apply(
  #   dlcm$thetas_accept[, , -(1:nwarmup)], 1
  #   , function(x) list(table(x))
  # )
  # domain_accept <- do.call(bind_rows, domain_accept)
  domain_accept <- table(dlcm$thetas_accept[, , -(1:nwarmup)])

  return(list(
    "thetas_avg"=thetas_avg, "domain_items"=domain_items, "domain_nitems"=domain_nitems, "domain_accept"=domain_accept
  ))
}


##############
############## SIMULATION
##############


#' Simulate random latent class model data
#' @param n integer. Number of observations
#' @param pis numeric vector. How often an observation is in class k.
#' @param thetas response probabilities. If a matrix of size Kxlength(pi) then the response probabilites of k bernouli items across length(pi) classes. Otherwise a list of length(pi) (per class) containing lists of length K (per item), containing vectors of response probabilities of each item.
#' @examples
#' \donotrun{
#' set.seed(4)
#' sim_list1 <- simulate_lcm(
#'   n = 1000
#' , pis = c(0.5, 0.5) # 2 classes equally likely
#'   , thetas = cbind(
#'     c(rep(0.2, 10), rep(0.8, 10)) # Response probabilities for 20 items in first class
#'     , c(rep(0.8, 10), rep(0.2, 10))
#' ))
#'
#' sims_list2 <- simulate_lcm(
#'  n=1000
#' , pis=c(0.5, 0.5) # 2 classes equally likely
#' , thetas_list = list(
#'   c( rep(list(c(0.8,0.2)), 20)
#'      , rep(list(c(0.6,0.3,0.05,0.05)), 2)
#'      , rep(list(c(rep(0.3,3), rep(0.1/5,5))), 1)
#'   ) # Response probabilities for 23 items in first class
#'   , c( rep(list(c(0.5,0.5)), 20)
#'        , rep(list(c(0.3,0.6,0.05,0.05)), 2)
#'        , rep(list(c(rep(0.4,2), rep(0.2/6,6))), 1)
#'   )
#' ))
#' }
#' @export
simulate_lcm <- function(n, pis, thetas) {

  nclasses <- length(pis)
  classes <- apply(rmultinom(n, 1, pis), 2, function(x) which(x==1))


  if (class(thetas) == "matrix") {
    responses <- sapply(
      classes
      , function(iclass) {runif(nrow(thetas)) > thetas[, iclass]}
    )+0
    responses <- t(responses)
  } else if (class(thetas) == "list") {
    nitems <- length(thetas[[1]])
    responses <- matrix(NA, nrow=n, ncol=nitems)
    for (iclass in 1:nclasses) {
      ifilter <- classes==iclass
      for (iitem in 1:nitems) {
        responses[ifilter, iitem] = apply(
          rmultinom(sum(ifilter), 1, thetas[[iclass]][[iitem]])
          , 2, function(x) which(x==1))
      }
    }
  }

  return(list(classes=classes, responses=responses))
}
