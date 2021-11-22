#' dependentLCM: Dependent Latent Class Model
#'
#' @docType package
#' @name dependentLCM
#' @useDynLib dependentLCM, .registration = TRUE
#' @importFrom dplyr "%>%"
NULL

NCLASS = 2
CLASSPI_ALPHA = 1
THETA_ALPHA = 1
CLASS_INIT_METHOD = "kmodes"
DOMAIN_PROPOSAL_EMPTY = 0.3
DOMAIN_PROPOSAL_SWAP = 0.2
STEPS_ACTIVE = c("thetas"=TRUE, "domains"=TRUE, "class_pi"=TRUE, "classes"=TRUE, "identifiable"=TRUE)
THETA_ALPHA_FUNNAMES = c("constant") # c("log", "average")
THETA_ALPHA_FUNNAME = THETA_ALPHA_FUNNAMES[1]
DOMAIN_MAXITEMS = 10
CLASS2DOMAIN_FUNS = list(
  "HOMO" = function(nclass) rep(0, nclass)
  , "HET" = function(nclass) seq(0, nclass-1)
)
CLASS2DOMAIN = "HOMO"

#' Fits a bayesian dependent LCM model
#' @param nitr integer. Number of iterations to run the bayes MCMC
#' @param cleanup boolean. TRUE to delete some redundant data from output
#' @inheritParams getStart_hparams
#' @inheritParams getStart_bayes_params
#' @export
dependentLCM_fit <- function(
  nitr
  # Data
  , df=NULL, mat=NULL
  # Hyperparameters
  ,nclass=NCLASS, ndomains=NULL, class2domain=CLASS2DOMAIN, classPi_alpha=CLASSPI_ALPHA, domain_maxitems=NULL, theta_alpha=THETA_ALPHA, domain_proposal_empty=DOMAIN_PROPOSAL_EMPTY, domain_proposal_swap=DOMAIN_PROPOSAL_SWAP, domain_nproposals=NULL, steps_active = STEPS_ACTIVE, theta_alpha_funname = THETA_ALPHA_FUNNAME
  # Bayes parameters
  , class_pi = NULL, classes = NULL, domains = NULL
  # Misc
  , class_init_method = CLASS_INIT_METHOD
  , warmup_settings = "default"
  ) {
  
  args <- as.list(environment()) # c(as.list(environment()), list(...))
  datetimes <- c("fun_start"=Sys.time())
  
  #
  # Set starting values
  #
  
  if (identical(warmup_settings, "default")) {
    if (identical(class2domain, "HET") | identical(class2domain, "het")) {
      # Default warmup for HET
      warmup_settings <- list()
      warmup_settings$nitr <- round(median(c(100, 1000, nitr/20)))
      warmup_settings$class2domain <- "HOMO"
    } else {
      # Otherwise default is no warmup
      warmup_settings = NULL
    }
  }
  
  if (is.null(mat)) {
    mat <- getStart_matrix(df)
  }
  
  hparams <- getStart_hparams(
    df=mat
    # Hyperparameters
    ,nclass=nclass, ndomains=ndomains, class2domain=class2domain, classPi_alpha=classPi_alpha, domain_maxitems=domain_maxitems, theta_alpha=theta_alpha, domain_proposal_empty=domain_proposal_empty, domain_proposal_swap=domain_proposal_swap, domain_nproposals=domain_nproposals, steps_active=steps_active, theta_alpha_funname=theta_alpha_funname
  )
  
  if (is.null(warmup_settings)) {
    # no warmup
    bayesparams <- getStart_bayes_params(
      mat=mat, hparams=hparams
      , class_pi = class_pi, classes = classes, domains = domains
      # Misc
      , class_init_method = class_init_method
    )
    warmup_sims <- NULL
  } else {
    # warmup
    warmup_args <- args # as.list(match.call())[-1]
    warmup_args[names(warmup_settings)] <- warmup_settings
    warmup_args$warmup_settings <- NULL
    warmup_sims <- do.call(dependentLCM_fit, warmup_args)
    bayesparams <- dlcm2paramargs(warmup_sims)$bayesparams
  }
  
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
  
  datetimes <- c(datetimes, "mcmc_start"=Sys.time())
  dlcm = dependentLCM_fit_cpp(x=(all_params$mat),hparams_list=all_params$hparams, params_list=all_params$bayesparams, nitr=nitr)
  datetimes <- c(datetimes, "mcmc_end"=Sys.time())
  
  #
  # Post Processing
  #
  
  # Clean up output
  dlcm$domains_id <- do.call(cbind, dlcm$domains_id)
  rownames(dlcm$domains_id) <- c("itr", "class", "domain", "pattern_id", "items_id")
  mode(dlcm$domains_id) <- "integer"
  dlcm$domains_patterns <- do.call(cbind, dlcm$domains_patterns)
  rownames(dlcm$domains_patterns) <- paste0("item_", 1:(nrow(dlcm$domains_patterns)))
  mode(dlcm$domains_patterns) <- "integer"
  dlcm$domains_lprobs <- unlist(dlcm$domains_lprobs)
  dlcm$domains_accept <- do.call(function(...) abind::abind(..., along=3), dlcm$domains_accept)
  dlcm$class_loglik <- do.call(function(...) abind::abind(..., along=3), dlcm$class_loglik)
  
  # name
  dlcm$class_pi <- set_dimnames(dlcm$class_pi, c("class", "itr"))
  dlcm$classes <- set_dimnames(dlcm$classes, c("obs", "itr"))
  dlcm$domains_accept <- set_dimnames(dlcm$domains_accept, c(NULL, NULL, "itr"))
  dlcm$class_loglik <- set_dimnames(dlcm$class_loglik, c("class", "obs", "itr"))
  
  # Supplemental
  domains_items <- apply(
    dlcm$domains_patterns > -1
    , 2, function(x) paste0(c("", which(x), ""), collapse=",")
  )
  domains_nitems <- as.integer(colSums(dlcm$domains_patterns > -1))
  domains_class2domain <- all_params$hparams$class2domain[dlcm$domains_id["class",,drop=TRUE]+1]
  domains_domain_id <- as.data.frame(t(rbind(dlcm$domains_id, domains_class2domain)))
  domains_domain_id <- do.call(paste, domains_domain_id) # paste rows together, faster than apply
  dlcm$domains_id["itr",] <- dlcm$domains_id["itr",] + 1L # start at 1
  
  # Merge domains attributes
  dlcm$domains <- data.frame(
    t(dlcm$domains_id)
    , class2domain = domains_class2domain
    , prob=exp(dlcm$domains_lprobs)
    , nitems = domains_nitems
    , items = domains_items
    , t(dlcm$domains_patterns)
    , stringsAsFactors = FALSE
  )
  all_params$hparams$domain_item_cols <- grep("^item_[0-9]+", colnames(dlcm$domains))
  
  dlcm$domains_merged <- as.data.frame(
    dlcm$domains 
    %>% dplyr::filter(pattern_id==0) 
    %>% dplyr::group_by(itr, class2domain)
    %>% filter(class == min(class))
    %>% dplyr::arrange(-nitems, items)
    %>% dplyr::summarize(domains_merged=paste(items, collapse="|"), .groups="keep")
  )
  
  # Delete redundant info
  dlcm$domains_id <- NULL
  dlcm$domains_patterns <- NULL
  dlcm$domains_lprobs <- NULL
  dlcm$nclass2domain <- NULL
  
  datetimes <- c(datetimes, "fun_end"=Sys.time())
  dlcm$runtimes <- as.numeric(c(diff(datetimes), total=tail(datetimes,1)-datetimes[1]))
  names(dlcm$runtimes) <- c("pre", "mcmc", "post", "total")
  
  return(list(hparams=all_params$hparams, mcmc=dlcm, warmup=warmup_sims))
}


#' Generate hyperparameters
#' Generate a list of hyperparameters, adding default values when necessary
#' @param df dataframe. Raw data you are fitting with all assumed to be factors
#' @param nitems integer. Number of columns of df.
#' @param nclass integer. Number of subject latent classes
#' @param ndomains integer. Number of item domains
#' @param class2domain integer vector of length nclass. Classes with same value have same domains.
#' If "HOMO" then defaults to c(0,0,...). If "HET" then defaults to c(0, 1, 2, ...)
#' @param classPi_alpha numeric vector. Bayes hyperparameter giving prior for Pi.
#' @param domain_maxitems iinteger. Maximum number of items which can be in a domain. Default is 10 (beyond 10 runs slowly).
#' @param theta_alpha numeric. Bayes hyperparemter giving the prior for theta/probabilties.
#' @param domain_proposal_empty numeric. Sets how often the domain metropolis proposal function pairs a nonempty domain with an empty domain.
#' @param domain_proposal_swap numeric. Sets how often the domain metropolis proposal function swaps items between domains.
#' @param domain_nproposals Sets how many times the domain metropolis propsal function is called each iteration
#' @param steps_active Named boolean vector of what actions to take during mcmc. If mcmc is skipped then initial values are kept as fixed.
#' thetas=TRUE to do gibbs on response probabilities, domains=TRUE to do metropolis on domains, class_pi=TRUE to do gibbs on class prior, classes=TRUE to do gibbs on class membership, identifiable=TRUE to check generic identifiability conditions of domains
#' @param theta_alpha_funname string Decides the prior for theta as domains get merged. 
#' "constant" for theta~Dirichlet(theta_alpha * rep(1, npatterns)). Produces domains of 1-4 items but typically not larger.
#' "log" for theta~Dirichlet(theta_alpha * ln(sum(npattens_i))/ln(npatterns) rep(1, npatterns))
#' "average" [depreciated] for theta~Dirichlet(theta_alpha * sum(npattens_i)/npatterns rep(1, npatterns)). Unstable for large domains
#' @keywords internal
getStart_hparams <- function(
  df=NULL, nitems=NULL
  # Hyperparameters
  ,nclass=NCLASS, ndomains=NULL, class2domain=NULL, classPi_alpha=CLASSPI_ALPHA, domain_maxitems=NULL, theta_alpha=THETA_ALPHA, domain_proposal_empty=DOMAIN_PROPOSAL_EMPTY, domain_proposal_swap=DOMAIN_PROPOSAL_SWAP, domain_nproposals=NULL, steps_active=STEPS_ACTIVE, theta_alpha_funname = THETA_ALPHA_FUNNAME
) {
  # Purpose: Add default hyperparameters
  
  if (is.null(nitems)) {
    nitems <- dim(df)[2]
  }
  item_nlevels <- apply(df, 2, max)+1
  
  # ndomains
  if (is.null(ndomains)){
    ndomains <- ndomains_singleton_mode(nitems, prop=2) # All singletons X2 more likely than next most likely choice of domains
  }
  
  # class2domain
  if ((typeof(class2domain) == "character") & length(class2domain)==1) {
    # convert string description into raw numbers we need
    class2domain = CLASS2DOMAIN_FUNS[[toupper(class2domain)]](nclass)
  }
  class2domain <- as.integer(factor(class2domain))-1 # Make sequential starting at 0
  nclass2domain <- length(unique(class2domain))

  # classPi_alpha
  if (length(classPi_alpha)==1) {
    classPi_alpha <- rep(classPi_alpha, nclass)
  }
  
  # domain_maxitems
  if (is.null(domain_maxitems)) {
    domain_maxitems <- min(nitems, DOMAIN_MAXITEMS)
  }
  
  # domain_nproposals
  if (is.null(domain_nproposals)) {
    domain_nproposals <- nitems
  }
  
  
  # steps_active (fill in missing values)
  steps_active_fn = STEPS_ACTIVE
  steps_active_fn[names(steps_active)] = steps_active
  
  return(list(
    nclass=nclass, ndomains=ndomains, class2domain=class2domain, classPi_alpha=classPi_alpha, domain_maxitems=domain_maxitems, theta_alpha=theta_alpha, nitems = nitems, item_nlevels = item_nlevels, nclass2domain = nclass2domain, domain_proposal_empty=domain_proposal_empty, domain_proposal_swap=domain_proposal_swap, domain_nproposals=domain_nproposals, steps_active = steps_active_fn, theta_alpha_funname = theta_alpha_funname
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
#' @param domains list. Initial values for domains/probabilites.
#' @param class_init_method string. Decides how 'classes' is defaulted if NULL. One of "kmodes" or "random" or "random_centers"
#' @keywords internal
getStart_bayes_params <- function(
  mat, hparams
  , class_pi = NULL, classes = NULL, domains = NULL
  # Misc
  , class_init_method = CLASS_INIT_METHOD
) {
  
  # classes
  if (is.null(classes)) {
    if (class_init_method=="kmodes") {
      classes <- getStart_class_kmodes(mat, hparams)
    }
    if (class_init_method=="random") {
      classes <- getStart_class_random(mat, hparams)
    }
    if (class_init_method=="random_centers") {
      classes <- getStart_class_random_centers(mat, hparams)
    }
  } else {
    classes = as.integer(factor(classes))-1
  }
  
  # class_pi
  if (is.null(class_pi)) {
    class_pi <- table(classes)+hparams$classPi_alpha
    class_pi <- class_pi / sum(class_pi)
  }
  
  # domains
  if (is.null(domains)) {
    domains <- getStart_domains(mat, classes, hparams)
  } # Future maybe allow other types of inputs
  
  return(list(
    class_pi = class_pi, classes = classes, domains = domains
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
getStart_class_kmodes <- function(mat, hparams, ...) {
  for (i in seq_len(3)) { 
    # Repeated attempts because klaR sometimes errors randomly
    tryCatch({
      clusters <- klaR::kmodes(mat, hparams$nclass, ...)$cluster-1
      break
    }
    , error=function(cond) {
      message("Warning: Issue with class_init_method = 'kmodes'. Retrying. If issue persists try class_init_method = 'random_centers'.\n Reference getStart_class_kmodes:  ", cond)
    }
    )
  }
  return(clusters)
}

#' Choose nclass random centers and put observations in their nearest center
#' Designed so that the chosen centers are as far from each other as possible
#' @param mat matrix. Raw data.
#' @param hparams list. List of hyperparameters
#' @param weight_prob bool. True if more common patterns should be more likely to be chosen as centers. False if all patterns should be equally likely.
#' @keywords internal
getStart_class_random_centers <- function(mat, hparams, weight_prob=TRUE) {
  nobs <- dim(mat)[1]
  nclass <- hparams$nclass
  attempts <- 3
  
  if (weight_prob==TRUE) {
    # More common patterns more likely to be chosen
    possible_centers <- mat # lazy copy (fast)
  } else {
    # All patterns equally likely to be chosen
    possible_centers <- unique(mat)
  }
  centers <- matrix(data=NA, nrow=nclass, ncol(mat))
  distance <- matrix(data=NA, nrow=nobs, ncol=nclass)
  
  # Get centers
  centers[1,] <-possible_centers[sample.int(n=nrow(possible_centers)
                                            , size=1),]
  distance[,1] = rowSums(sweep(mat, 2, centers[1,], FUN="=="))
  for (i in seq_len(nclass-1)+1) {
    # choose the center farthest from the current centers
    min_dist <- apply(distance, 1, min, na.rm=TRUE)
    center_inds <- which(min_dist==min(min_dist))
    centers[i,] <- possible_centers[sample(center_inds, size=1),]
    distance[,i] = rowSums(sweep(mat, 2, centers[i,], FUN="=="))
  }
  
  # Associate observation to nearest class
  classes <- apply(distance, 1, FUN=function(idist) sample(which(idist==min(idist)),1))
  
  return(classes-1)
}

#' Choose starting domain values.
#' Note that the initial choice of theta (although set) is unimportant since we recalculate theta before applying it in MCMC/CPP.
#' @param mat matrix. Raw data.
#' @param classes integer vector. The class of each observation.
#' @param hparams list. List of hyperparameters
#' @keywords internal
getStart_domains <- function(mat, classes, hparams) {
  # FUTURE: Maybe switch to internal C gibbs here? Keep this for now for validation
  domains_list = list()
  for (iclass in 0:(hparams$nclass-1)) {
    domains_list[[iclass+1]] <- lapply(
      1:hparams$nitems
      , function(icol) {
        ix = c(mat[classes==iclass, icol] # actual data
               , seq_len(hparams$item_nlevels[icol])-1 # +1 of each category (avoid zeros)
               )
        ithetas = as.vector(table(ix)) / length(ix)
        iout = list(
          items=c(icol-1)
          , thetas = ithetas
        )
        return(iout)
      }
    )
  }
  return(domains_list)
}

#' Check parameters for internal consistency
#' @param all_params list. List of parameters and hyperparameters
#' @keywords internal
check_params <- function(all_params) {
  is_problem = FALSE
  
  if (all_params$hparams$nclass != length(all_params$hparams$class2domain)) {
    warning("nclass~class2domain mismatch")
    is_problem = TRUE
  }
  
  if (all_params$hparams$ndomains < all_params$hparams$nitems) {
    warning("Must have more domains than items") # Fewer domains theoretically ok, but not implemented in code
    is_problem = TRUE
  }
  
  if (!setequal(names(all_params$hparams$steps_active), names(STEPS_ACTIVE))) {
    warning("steps_active invalid")
    is_problem = TRUE
  }
  
  if (max(all_params$bayesparams$classes) > all_params$hparams$nclass) {
    warning("classes > nclass")
    is_problem = TRUE
  }
  
  return(is_problem)
}

#' Take dependent latent class model (dlcm) output (i.e. from dependentLCM_fit)
#' and convert it into Hyper/Bayes parameter arguments to put into dependentLCM_fit.
#' Used namely for nesting dependentLCM_fit.
#' @param dlcm Dependent latent class model from dependentLCM_fit()
#' @param iter which iteration to pull values from (by default last itration)
#' @keywords internal
dlcm2paramargs <- function(dlcm, iter=NULL) {
  iter <- dlcm$mcmc$maxitr
  nitems <- dlcm$hparams$nitems
  
  domains_mat <- dlcm$mcmc$domains %>% filter(itr == iter)
  domains_mat$items <- (
    # apply(id2pattern(domains_mat$items_id, rep(2, nitems))==1, 2, which)
    unlist(apply(domains_mat[,dlcm$hparams$domain_item_cols] >= 0, 1, function(x) list(which(x)-1)), recursive = FALSE)
  )
  
  domains_list <- (domains_mat 
                   %>% group_by(class, domain)
                   %>% arrange(pattern_id) 
                   %>% dplyr::group_map(function(.x, .y) list(
                     items = unname(unlist(.x[1,"items"]))
                     , thetas=.x$prob
                     , name=.y # for filering in next step
                   )))
  domains_list <- lapply(
    seq_len(dlcm$hparams$nclass)-1
    , function(iclass) {
      domains_list[sapply(domains_list, function(x) x$name$class==iclass)]
    }
  )
  
  bayesparams <- list(
    class_pi = unname(dlcm$mcmc$class_pi[,iter])
    , classes = unname(dlcm$mcmc$classes[,iter])
    , domains =  domains_list
  )
  
  all_params <- list(
    hparams = dlcm$hparams
    , bayesparams = bayesparams
  )
  
  return(all_params)
}

#' What is the prior probability of this choice of domains?
#' Ignores identifiability restrictions
#' @param x IntegerVector. The number of items in each domain. Ok to omit 0's
#' @param D Integer. The total number of domains (including empty domains)
#' @param specific_items Boolean. 
#' If FALSE, we look for any domain which produces domains of this size regardless of what specific items they contain.
#' IF TRUE, we fix which items are in which domain, and calculate the probability of grouping these specific items together.
#' @param log Boolean. If TRUE give probability in log scale.
#' Assumes x are logged values. Calculates sum(e^x) and then converts back to log scale
#' @param x numeric vector in log scale
#' @export
ldomain_prior <- function(x, ndomains, specific_items=FALSE, log=TRUE) {
  
  if (specific_items==FALSE) {
    lprob <- (
      lfactorial(sum(x)) - sum(sapply(x, lfactorial)) # unique groups
      - sum(sapply(table(x[x>0]), lfactorial)) # non-unique groups
      + lfactorial(ndomains) - lfactorial(ndomains-sum(x > 0)) # permuting groups
      - sum(x) * log(ndomains) # denominator
    )
  } else { # specific_items==TRUE
    lprob <- (
      lfactorial(ndomains) - lfactorial(ndomains-sum(x > 0)) # permuting groups
      - sum(x) * log(ndomains) # denominator
    )
  }
  if (log==FALSE) {
    return(exp(lprob))
  }
  return(lprob)
}


#' How many domains (ndomains) we we need before 'all singleton domains' are prop-times more likely than 'one 2-item domain with rest singletons'? Using the prior only.
#' Solves: ldomain_prior(rep(1, nitems), ndomains, specific_items=FALSE, log=TRUE) - ldomain_prior(c(2,rep(1, nitems-2)), ndomains, specific_items=FALSE, log=TRUE) = log(prop)
#' @param nitems integer. Number of items in the data.
#' @param prop float. prop=1 forces singleton domains to be mode. prop>1 makes singleton domains increasingly frequent.
#' @export
ndomains_singleton_mode <- function(nitems, prop=1) {
  ceiling(nitems+prop*nitems*(nitems-1)/2-1)
}

##############
############## UTILITIES
##############

CAST_FUN <- list("integer"=as.integer, "double"=as.numeric)

#' Calculate mode (most common value)
#' @keywords internal
getMode <- function(x) {
  xtype <- typeof(x)
  xcounts <- table(x)
  xmode <- names(xcounts)[which.max(xcounts)]
  xmode <- CAST_FUN[[xtype]](xmode) 
  return(xmode)
}


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

#' Robust estimate for variance. Inspired by MCD but fast approximation instead
#' @param x numeric vector. Data we wish to calculate variance of
#' @param alpha percentage. Remove the alpha most extreme points before calculating variance.
#' @keywords internal
robust_var <- function(x, alpha) {
  if (alpha <= 0) {
    # Not robust
    return(var(x))
  }
  
  center <- median(x)
  center_diff <- abs(x-center)
  v <- var(x[center_diff < quantile(center_diff, 1-alpha)])
  return(v) 
}

#' Same as log(sum(exp(x))), but adjusted to improve precision.
#' Assumes x are logged values. Calculates sum(e^x) and then converts back to log scale
#' @param x numeric vector in log scale
#' @keywords internal
expSumLog <- function(x) {
  # Handle precision: log(e^a+e^b+e^c) = log(a * (1+e^(b-a)+e^(c-1))) = log(a) + log(1+e^(b-a)+e^(c-1)))
  xmax <- max(x)
  return(xmax + log(sum(exp(x-xmax))))
}

#' @name set_dimnames
#' @title set_dimnames
#' @description Names each axis of your array by naming each row/column axisName#
#' @param xarray The array you wish to name the axis of
#' @param axis_names The name of each axis of this array
#' @keywords internal
set_dimnames <- function(xarray, axis_names) {
  dimnames(xarray) <- lapply(
    seq_along(axis_names)
    , function(iaxis) {
      iaxis_name <- axis_names[iaxis]
      if (is.null(iaxis_name)) {
        return(null)
      }
      paste0(
        iaxis_name
        , seq_len(dim(xarray)[iaxis])
      )
    }
  )
  return(xarray)
}

`%>%` <- magrittr::`%>%`

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
dlcm.summary <- function(dlcm, nwarmup=NULL) {
  
  if (is.null(nwarmup)) {
    nwarmup = min(1000, floor(dlcm$mcmc$maxitr/2))
  }
  
  dlcm$mcmc$domains$item_value <- apply(
    dlcm$mcmc$domains[, dlcm$hparams$domain_item_cols]
    , 1, function(x) paste0(x[x>-1], collapse=", ")
  )
  thetas_avg <- dlcm$mcmc$domains %>% dplyr::filter(itr > nwarmup) %>% dplyr::group_by(class, items, item_value) %>% dplyr::summarize(n=dplyr::n(), prob=mean(prob), .groups="keep")
  
  domain_items_all <- dlcm$mcmc$domains_merged %>% dplyr::filter(itr > nwarmup) %>% dplyr::group_by(class2domain, items=domains_merged) %>% dplyr::summarize(n=dplyr::n(), .groups="keep")
  domain_items_all <- domain_items_all %>% group_by(class2domain) %>% mutate(perc = n / sum(n)) %>% dplyr::arrange(-perc)
  domain_items <- dlcm$mcmc$domains %>% dplyr::filter(pattern_id==0, itr > nwarmup) %>% dplyr::group_by(class2domain, items) %>% filter(class == min(class)) %>% dplyr::summarize(nitems=max(nitems), n=dplyr::n(), .groups="keep") %>% dplyr::arrange(-n)
  
  # domain_nitems <- table((dlcm$mcmc$domains %>% dplyr::filter(pattern_id==0, itr > nwarmup))$nitems)
  
  domain_accept <- table(dlcm$mcmc$domains_accept[, , -(1:nwarmup)])
  
  waic <- dlcm.get_waic(dlcm, itrs=nwarmup:dlcm$mcmc$maxitr)
  names(waic) <- paste0("waic_", names(waic))
  
  classes <- unname(apply(dlcm$mcmc$classes[,-(1:nwarmup)], 1, getMode))
  class_pi = rowMeans(dlcm$mcmc$class_pi[,-(1:nwarmup), drop=FALSE])
  
  return(c(
    list("thetas_avg"=thetas_avg, "domain_items"=domain_items, "domain_items_all"=domain_items_all, "domain_accept"=domain_accept, "classes"=classes, "class_pi"=class_pi)
    , waic
    ))
}

#' Calculate likelihood and WAIC
#' @param this_sim. Dependent Latent class model
#' @param itrs integer vector. Which iterations should be include in calculation?
#' @export
dlcm.get_waic <- function(this_sim, itrs=NULL) {
  
  alpha <- 0 # No robust MCD variance estimate in WAIC
  
  if (is.null(itrs)) {
    # default in all iterations
    itrs = 1:this_sim$mcmc$maxitr
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
    this_sim$mcmc$domains 
    %>% filter(itr %in% itrs) 
    %>% group_by(itr) 
    %>% summarize(nparams = n() - length(unique(domain)), .groups="keep"))
  likelihoods$nparameters <- (
    unlist(theta_parameters[match(likelihoods$itr, theta_parameters$itr), "nparams"])
    + this_sim$hparams$nclass-1 # class pis
    + length(unique(this_sim$hparams$class2domain)) * (this_sim$hparams$nitems-1) # domains
  )
  summary["nparams_avg"] <- mean(likelihoods$nparameters)
  
  #
  # Likelihood and WAIC
  #
  
  # Calculate likelihood for each observation in each iteration
  obsLogLiks <- sweep(this_sim$mcmc$class_loglik[,,itrs, drop=FALSE]
                      , c(1,3) # Include pi. Repeat for each observation (2)
                      , log(this_sim$mcmc$class_pi[,itrs, drop=FALSE]), "+")
  obsLogLiks <- apply(obsLogLiks, c(2,3), expSumLog)
  
  likelihoods$logLik <- colSums(obsLogLiks, 1) # sum across observations
  
  # Same as: LaplacesDemon::WAIC
  summary["logLik_avg"] <- mean(likelihoods$logLik)
  summary["lppd"] <- sum(apply(obsLogLiks, 1, expSumLog) - log(dim(obsLogLiks)[2])) 
  summary["waic_nparams1"] <- 2 * sum(apply(obsLogLiks, 1, expSumLog) - log(dim(obsLogLiks)[2])
                                      - apply(obsLogLiks, 1, mean))
  summary["waic_nparams2"] <- sum(apply(obsLogLiks, 1, robust_var, alpha=alpha))
  summary["waic1"] <- -2 * (summary["lppd"] - summary["waic_nparams1"])
  summary["waic2"] <- -2 * (summary["lppd"] - summary["waic_nparams2"])
  summary["aic"] <- -2*summary["logLik_avg"] + 2*summary["nparams_avg"]
  
  return(list(
    summary=summary
    , df=likelihoods
  ))
}

#' Calculate the probabilities of different patterns
#' @param items integer vector. Which items should we include? (indexed starting at 1)
#' @param this_sim a mcmc simulation. Output from dependentLCM_fit()
#' @param itrs integer vector. Which iterations from this_sim should we include?
#' @param merge_itrs TRUE if you want average probabilities across simulations. FALSE if you want to examine each iteration individually
#' @description Useful when getting probilities which cross domains.
#' @returns Probabilities for each pattern of items
#' @export
theta_item_probs <- function(items, this_sim, itrs, merge_itrs=TRUE) {
  
  itr_filter <- which(this_sim$mcmc$domains$itr %in% itrs)
  nitr <- length(itrs)
  nitems <- length(items)
  items_colnames <- colnames(this_sim$mcmc$domains)[this_sim$hparams$domain_item_cols[items]]
  
  # Reduce thetas down to relevant patterns only (merging redundant patterns as necessary)
  thetas_agg <- aggregate(
    formula(paste0("prob ~ ", paste(c(items_colnames, "class", "itr"), collapse=" + ")))
    , this_sim$mcmc$domains[itr_filter, ]
    , sum
  )
  thetas_agg <- thetas_agg[!(rowSums(thetas_agg[, 1:nitems, drop=FALSE]) == -nitems), ] # Remove unrelated rows
  thetas_agg$prob_log <- log(thetas_agg$prob)
  
  all_patterns <- expand.grid(
    lapply(this_sim$hparams$item_nlevels[items], function(n) 0:(n-1))
  )
  colnames(all_patterns) <- items_colnames
  
  prob_helper <- function(ipattern) {
    # Purpose: Calculate the probability for this pattern
    # Implict: thetas_agg, this_sim, nitr, nitems
    imatch_pattern <- rowSums(
      sweep(thetas_agg[,1:nitems, drop=FALSE], 2, ipattern, "==") # Matches pattern
      | thetas_agg[,1:nitems, drop=FALSE] < 0 # Or is undefined
    ) == nitems
    iprobs <- thetas_agg %>% filter(imatch_pattern) %>% group_by(itr, class) %>% summarize(prob=exp(sum(prob_log)), .groups="keep")
    iprobs$pi <- this_sim$mcm$class_pi[cbind(iprobs$class+1, iprobs$itr)]
    iprobs$prob_pi <- iprobs$prob * iprobs$pi
    return(iprobs)
  }
  
  if (merge_itrs==TRUE) {
    # Get average across iterations
    all_patterns$prob <- apply(all_patterns, 1
                               , FUN = function(ipattern) {
                                 iprobs <- prob_helper(ipattern)
                                 return(sum(iprobs$prob * iprobs$pi) / nitr)
                               })
    out <- all_patterns
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

#' For each iteration calculate the probability that a given observation is in each class.
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
#' , thetas = list(
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
    responses <- responses - 1
  }
  
  return(list(classes=classes, responses=responses))
}

#' Take the parameters of a simulated lcm (simulate_lcm)
#' and produce the probability of seeing this observation
#' @inheritParams simulate_lcm
#' @param xobs numeric vector. One vector observation
#' @note Helper function for sim_logprobs. Currently only supports thetas in list form.
#' @keywords internal
sim_logprob_one <- function(xobs, pis, thetas) {
  nclass <-length(pis)
  class_logprobs <- sapply(
    1:nclass
    , function(iclass){
      sum(log(mapply(
        function(theta, xobs) theta[xobs+1]
        , theta=thetas[[iclass]]
        , xobs=xobs
      )))}
  )
  class_logprobs <- class_logprobs + log(pis)
  return(expSumLog(class_logprobs))
}

#' Take the parameters of a simulated lcm (simulate_lcm)
#' and produce the probability of seeing each observation given
#' @inheritParams simulate_lcm
#' @param x numeric matrix. Data with one row per observation.
#' @note Currently only supports thetas in list form. FUTURE extend to include matrix thetas as well.
#' @export
sim_logprobs <- function(x, pis, thetas) {
  return(apply(x, 1, sim_logprob_one, pis=pis, thetas=thetas))
}
