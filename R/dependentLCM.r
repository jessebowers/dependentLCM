#' dependentLCM: Dependent Latent Class Model
#' @description 
#' Latent Class Models (LCMs) are used to cluster multivariate categorical data (e.g. group participants based on survey responses). Traditional LCMs assume a property called conditional independence. This assumption can be restrictive, leading to model misspecification and overparameterization. To combat this problem, we developed a novel Bayesian model called a Dependent Latent Class Model (DLCM), which permits conditional dependence. Compared to traditional LCMs, DLCMs are effective in applications with time series, overlapping items, and structural zeroes.
#' 
#' The primary function is dependentLCM_fit. 
#' 
#' Bowers, J., & Culpepper, S. (2022). Dependent Latent Class Models (Version 1). arXiv. \url{https://doi.org/10.48550/ARXIV.2205.08677}
#' 
#' \url{https://arxiv.org/abs/2205.08677}
#' 
#' \url{https://github.com/jessebowers/dependentLCM}
#' 
#' @docType package
#' @name dependentLCM
#' @useDynLib dependentLCM, .registration = TRUE
#' @importFrom dplyr "%>%"
#' @importFrom Rcpp evalCpp
NULL

# Default values
NCLASS = 2
CLASSPI_ALPHA = 1
THETA_ALPHA = 1
CLASS_INIT_METHODS = c("random_centers_polar", "random_centers", "random", "kmodes")
DOMAIN_PROPOSAL_EMPTY = 0.3
STEPS_ACTIVE = c("thetas"=TRUE, "domains"=TRUE, "class_pi"=TRUE, "classes"=TRUE, "identifiable"=TRUE, "classLikelihood"=TRUE)
DOMAIN_MAXITEMS = 10
CLASS2DOMAIN_FUNS = list(
  "HOMO" = function(nclass) rep(0, nclass)
  , "HET" = function(nclass) seq(0, nclass-1)
)
CLASS2DOMAINS = names(CLASS2DOMAIN_FUNS)

#' Fits a bayesian dependent LCM model
#' @inheritParams getStart_hparams
#' @inheritParams getStart_bayes_params
#' @inheritParams doWarmup
#' @return Returns a list with three items. hparams describes the hyperparmeters chosen by the model. These match the arguments input into the function plus any default values. mcmc describes the simulated (fitted) bayes parameter values. warmup describes extra/special warmup iterations using the 'warmup_settings' argument.
#' 
#' The mcmc item includes the following values. mcmc does NOT discard any warmup iterations.
#' \itemize{
#' \item{"class_pi"}{=A matrix describing the prior class probabilities of each class in each MCMC iteration. It contains one row per class one column per iteration, and each cell contains the prior probability of a subject being in that class.}
#' \item{"classes"}{=A matrix describing what class each subject belongs to. It contains one row per subject, one column per MCMC iteration, and each cell identifies the class of that subject in that iteration.}
#' \item{"nitr"}{=Gives the total number of iterations processed including warmup iterations}
#' \item{"domains_accept"}{=Each MCMC iteration, we attempt to change the domain structure #domain_nproposals times with a metropolis step. This 3-dimensional array describes whether the proposed change was accepted or rejected. A value of -2 indicates the proposal was rejected due to identifiability or max item constraints. A value of -1 or 2 indicates the proposed change is (equivalent to) no change. A value of 0 indicates the proposal was rejected. A value of 1 indicates the proposal was accepted.. The dimensions of the array are as follows. The third dimension has one slice per iteration. The first dimension gives one slice per 'domain_nproposals'. The second dimension has between 1 and nclass slices. A homogeneous DLCM has 1 slice, a heterogeneous DLCM has nclass slices, and a partially heterogeneous DLCM may be in-between.}
#' \item{"class_loglik"}{=A 3-dimensional array. The first dimension has one slice per class. The second dimension has one slice per subject/observation. The third dimension has one slice per MCMC iteration. Each cell contains the log likelihood that we would observe the response pattern given by this subject, if the subject was in this class, based on the parameters in this iteration.}
#' \item{"troubleshooting"}{=Used for investigationg bugs. Should be empty.}
#' \item{".Random.seed"}{=The state of the random seed before this function was executed. If you set the seed to this state and run again you should get the same result.}
#' \item{"domains"}{=Dataframe with one row for iteration X domain X pattern. For an MCMC iteration, it identifies what domains there are, and for each class what's the probability of getting a given response pattern to a given domain. It contains the following columns:
#' \itemize{
#' \item{"itr"}{=What MCMC iteration is this?}
#' \item{"class"}{=What class are we calculating the response probabilities for?}
#' \item{"domain"}{=ID. What domain is this?}
#' \item{"pattern_id"}{=Integer uniquely identifying a response pattern to thie items in this domain. We will calculate the probability of this response pattern. pattern_id needs to be paired with a items_id to make sense of it.}
#' \item{"items_id"}{=Integer uniquely identifying what items (what set of items) are in this domain.}
#' \item{"class2domain"}{=For heterogeneous and partially heterogeneous DLCMs, this identifies which group of latent classes this domain belongs to.}
#' \item{"prob"}{=What's the probability of getting this response pattern?}
#' \item{"nitems"}{=How many items are in this domain?}
#' \item{items"}{=String listing the items in this domain. Function of items_id.}
#' \item{item_#"}{=For each item #, gives the specific value of that item in this response pattern. A value of -1 indicates this item is not in this domain. item_# is a function of (items_id, pattern_id).}
#' }}
#' \item{"domains_merged"}{=Dataframe describing what domains were chosen for each MCMC iteration. There is one row for MCMC iteratation X class2domain. For (partially) heterogeneous DLCMs, class2domain allows different classes to have different domains. The string column domains_merged describes the domains with vertical bars "|" separating domains, and commas "," separating items wtihin a given domain.}
#' \item{"runtimes"}{=Describes how long the function ran in seconds. pre records the seconds of preprocessing before starting MCMC. mcmc describes how long it took to run the mcmc iterations. post describes how long it took to aggregate/transform the data after the MCMC is completed. Total gives the total time. There are 'secret' troubleshooting steps to get the runtime of specific C++ functions executed by this algorithm (these are stored under 'troubleshooting').}
#' }
#' @examples
#' \dontrun{
#' library(dependentLCM)
#' 
#' # Get Data
#' library(pks)
#' data(probability, package="pks")
#' xdf <- na.omit(probability[,c("b101", "b102", "b103", "b104", "b105", "b106", "b107", "b108", "b109", "b110", "b111", "b112", "b201", "b202", "b203", "b204", "b205", "b206", "b207", "b208", "b209", "b210", "b211", "b212")])
#' 
#' # Run Model
#' set.seed(4)
#' dlcm <- dependentLCM_fit(
#'   nitr = 6000
#'   , df=xdf
#'   , nclass=3
#' )
#' dlcm$summary <- dlcm.summary(dlcm, nwarmup=1000)
#' 
#' 
#' # Class of each observation
#' dlcm$summary$classes
#' 
#' # Which items are grouped together because they show local dependence? See ?dlcm.summary. We call groups of locally dependent items 'domains'.
#' dlcm$summary$domain_items_all
#' 
#' # Average response probabilities
#' dlcm$summary$thetas_avg_mode
#' }
#' @export
dependentLCM_fit <- function(
    nitr
    # Data
    , df=NULL, mat=NULL
    # Hyperparameters
    ,nclass=NCLASS, ndomains=NULL, class2domain=CLASS2DOMAINS[1], classPi_alpha=CLASSPI_ALPHA, domain_maxitems=NULL, theta_alpha=THETA_ALPHA, domain_proposal_empty=DOMAIN_PROPOSAL_EMPTY, domain_nproposals=NULL, steps_active = STEPS_ACTIVE
    # Bayes parameters
    , class_pi = NULL, classes = NULL, domains = NULL
    # Misc
    , class_init_method = CLASS_INIT_METHODS[1]
    , warmup_settings = "default", warmup_dlcm=NULL
    , collapse_classes_itrs = 0
) {
  
  #
  # Save starting info
  #
  
  args <- as.list(environment()) # c(as.list(environment()), list(...))
  datetimes <- c("fun_start"=Sys.time())
  
  # Save seed for reproducibility
  if (!exists(".Random.seed")) {
    set.seed(NULL)
  }
  .Random.seed_start <- .Random.seed
  
  #
  # Set initial parameter values
  #
  
  if (is.null(mat)) {
    mat <- getStart_matrix(df)
  }
  
  hparams <- getStart_hparams(
    nitr=nitr, df=mat
    # Hyperparameters
    ,nclass=nclass, ndomains=ndomains, class2domain=class2domain, classPi_alpha=classPi_alpha, domain_maxitems=domain_maxitems, theta_alpha=theta_alpha, domain_proposal_empty=domain_proposal_empty, domain_nproposals=domain_nproposals, steps_active=steps_active, collapse_classes_itrs=collapse_classes_itrs
  )
  
  warmup_dlcm <- doWarmup(
    args=args, warmup_settings=warmup_settings, warmup_dlcm=warmup_dlcm
  )
  
  bayesparams <- getStart_bayes_params(
    mat=mat, hparams=hparams
    , class_pi = class_pi, classes = classes, domains = domains
    # Misc
    , class_init_method = class_init_method
    , warmup_dlcm = warmup_dlcm
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
  
  datetimes <- c(datetimes, "mcmc_start"=Sys.time())
  dlcm = dependentLCM_fit_cpp(x=(all_params$mat),hparams_list=all_params$hparams, params_list=all_params$bayesparams)
  datetimes <- c(datetimes, "mcmc_end"=Sys.time())
  
  #
  # Post Processing
  #
  
  # Clean up output
  dlcm$domains_id <- do.call(cbind, dlcm$domains_id)
  rownames(dlcm$domains_id) <- c("itr", "class", "domain", "pattern_id", "items_id")
  mode(dlcm$domains_id) <- "integer"
  dlcm$domains_lprobs <- unlist(dlcm$domains_lprobs)
  dlcm$domains_accept <- do.call(function(...) abind::abind(..., along=3), dlcm$domains_accept)
  dlcm$class_loglik <- do.call(function(...) abind::abind(..., along=3), dlcm$class_loglik)
  dlcm$class_loglik_collapsed <- do.call(function(...) abind::abind(..., along=3), dlcm$class_loglik_collapsed)
  
  # name
  dlcm$class_pi <- set_dimnames(dlcm$class_pi, c("class", "itr"))
  dlcm$classes <- set_dimnames(dlcm$classes, c("obs", "itr"))
  dlcm$domains_accept <- set_dimnames(dlcm$domains_accept, c(NULL, NULL, "itr"))
  dlcm$class_loglik <- set_dimnames(dlcm$class_loglik, c("class", "obs", "itr"))
  # dlcm$class_loglik_collapsed <- set_dimnames(dlcm$class_loglik_collapsed, c("class", "obs", "itr"))
  
  # domains_patterns 
  dlcm$domains_patterns <- itemid2patterns(dlcm$domains_id["pattern_id",], dlcm$domains_id["items_id",], hparams[["item_nlevels"]])
  rownames(dlcm$domains_patterns) <- paste0("item_", seq_len(nrow(dlcm$domains_patterns)))
  mode(dlcm$domains_patterns) <- "integer"
  
  # Supplemental
  domains_items <- get_which_strs(dlcm$domains_patterns > -1)
  domains_nitems <- as.integer(colSums(dlcm$domains_patterns > -1))
  domains_class2domain <- all_params$hparams$class2domain[dlcm$domains_id["class",,drop=TRUE]+1]
  dlcm$domains_id["itr",] <- dlcm$domains_id["itr",] + 1L # start at 1
  dlcm$.Random.seed <- .Random.seed_start
  
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
    %>% dplyr::filter(class == min(class))
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
  
  return(list(hparams=all_params$hparams, mcmc=dlcm, warmup=warmup_dlcm))
}


#' Generate hyperparameters
#' Generate a list of hyperparameters, adding default values when necessary
#' @param nitr integer. Number of iterations to run the bayes MCMC
#' @param df dataframe. Raw data you are fitting with all assumed to be factors. Assumes no missing data (no NAs).
#' @param nitems integer. Number of columns of df.
#' @param nclass integer. Number of subject latent classes
#' @param ndomains integer. Number of item domains
#' @param class2domain integer vector of length nclass. Classes with same value have same domains.
#' If "HOMO" then defaults to c(0,0,...). If "HET" then defaults to c(0, 1, 2, ...)
#' @param classPi_alpha numeric vector. Bayes hyperparameter giving prior for Pi.
#' @param domain_maxitems iinteger. Maximum number of items which can be in a domain. Default is 10 (beyond 10 runs slowly).
#' @param theta_alpha numeric. Bayes hyperparemter giving the prior for theta/probabilties.
#' @param domain_proposal_empty numeric. Sets how often the domain metropolis proposal function pairs a nonempty domain with an empty domain.
#' @param domain_nproposals Sets how many times the domain metropolis propsal function is called each iteration
#' @param steps_active Named boolean vector of what actions to take during mcmc. If mcmc is skipped then initial values are kept as fixed.
#' thetas=TRUE to do gibbs on response probabilities, domains=TRUE to do metropolis on domains, class_pi=TRUE to do gibbs on class prior, classes=TRUE to do gibbs on class membership, identifiable=TRUE to check generic identifiability conditions of domains, classLikelihood=TRUE to get the likelihood that each observation is in each class.
#' @param collapse_classes_itrs integer. Which iterations should used collasped Gibbs when sampling class membership? This allows collapsing on class prior and response probability when evaluating class. The first collapse_classes_itrs (integer) iterations are collapsed.
#' @keywords internal
getStart_hparams <- function(
    nitr, df=NULL, nitems=NULL
    # Hyperparameters
    ,nclass=NCLASS, ndomains=NULL, class2domain=NULL, classPi_alpha=CLASSPI_ALPHA, domain_maxitems=NULL, theta_alpha=THETA_ALPHA, domain_proposal_empty=DOMAIN_PROPOSAL_EMPTY, domain_nproposals=NULL, steps_active=STEPS_ACTIVE, collapse_classes_itrs=0
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
    nitr=nitr, nclass=nclass, ndomains=ndomains, class2domain=class2domain, classPi_alpha=classPi_alpha, domain_maxitems=domain_maxitems, theta_alpha=theta_alpha, nitems = nitems, item_nlevels = item_nlevels, nclass2domain = nclass2domain, domain_proposal_empty=domain_proposal_empty, domain_nproposals=domain_nproposals, steps_active = steps_active_fn
    , theta_alpha_funname = "constant", collapse_classes_itrs = collapse_classes_itrs
  ))
}

#' Convert dataframe of factors to matrix of integers (indexed starting at 0)
#' @keywords internal
getStart_matrix <- function(df) {
  # Purpose: Convert to factor in integer form (0:k-1) and remove NA rows
  
  df <- as.data.frame(df)
  for (i in seq_len(ncol(df))) {
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
#' @param domains list. Initial values for domains/probabilites. Should be of the form:
#' domains = list(
#' class1 = list(
#'   domain1 = list(
#'     items = c(...) # items in this domain
#'     , thetas = c(...) # probabilities of each response to these items ordered per dependentLCM::id2pattern()
#'   )
#'   , domain2 = list(items=c(...), thetas=c(...))
#'   , ...
#' )
#' , class2 = list(
#'   domain1 = list(items=c(...), thetas=c(...))
#'   , domain2 = list(items=c(...), thetas=c(...))
#'   , ...
#' )
#' , ...
#' )
#' @param class_init_method string. Decides how 'classes' is defaulted if NULL. One of "kmodes" or "random" or "random_centers", "random_centers_polar"
#' @param warmup_dlcm list. A past DLCM fit (from dependentLCM_fit). The last iteration of the DLCM is used as the Bayes parameters.
#' @keywords internal
getStart_bayes_params <- function(
    mat, hparams
    , class_pi = NULL, classes = NULL, domains = NULL
    # Misc
    , class_init_method = CLASS_INIT_METHODS[1]
    , warmup_dlcm = NULL
) {
  
  # warmup_dlcm takes precidence
  # maybe in the future allow other parameters to overwrite??
  if (!is.null(warmup_dlcm)) {
    bayesparams <- dlcm2paramargs(warmup_dlcm)$bayesparams
    return(bayesparams)
  }
  
  # classes
  if (is.null(classes)) {
    if (hparams$nclass == 1) {
      classes <- rep(0, dim(mat)[1])
    } else if (identical(class_init_method,"kmodes")) {
      classes <- getStart_class_kmodes(mat, hparams)$classes
    } else if (identical(class_init_method,"random")) {
      classes <- getStart_class_random(mat, hparams)$classes
    } else if (identical(class_init_method,"random_centers")) {
      classes <- getStart_class_random_centers(mat, hparams, isWeighted=FALSE, isPolar=FALSE)$classes
    }  else if (identical(class_init_method,"random_centers_polar")) {
      classes <- getStart_class_random_centers(mat, hparams, isWeighted=FALSE, isPolar=TRUE)$classes
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

#' Does *EXTRA* MCMC warmup step if applicable
#' @param args list All arguments given to dependentLCM_fit
#' @param warmup_settings list. Optionally allows for an extra warmup cycle using different parameters vs the arguments given for the main process. Any parameters in the given list will overwrite the arguments of the same name in the warmup cycle.
#' @param warmup_dlcm list. A past DLCM fit (from dependentLCM_fit). The last iteration of the DLCM is used as the Bayes parameters. Takes precidence over warmup_settings.
#' @keywords internal
doWarmup <- function(args, warmup_settings=NULL, warmup_dlcm=NULL) {
  
  if (!is.null(warmup_dlcm)) {
    return(warmup_dlcm) # takes precedence
  }
  
  # get default warmup_settings
  class2domain <- args[["class2domain"]]
  nitr <- args[["nitr"]]
  if (identical(warmup_settings, "default")) {
    if (identical(class2domain, "HET") | identical(class2domain, "het") | (length(unique(class2domain))>1)) {
      # HET and partially HET default
      warmup_settings <- list()
      warmup_settings$nitr <- round(median(c(100, 1000, nitr/20)))
      warmup_settings$class2domain <- "HOMO"
    } else {
      # Other situations have no default
      warmup_settings = NULL
    }
  }
  
  if (is.null(warmup_settings)) {
    return(NULL) # No warmup to do
  }
  
  
  # use args anywhere warmup_settings is not specified
  warmup_settings <- c(
    warmup_settings
    , args[-which(names(args) %in% c(names(warmup_settings), "warmup_settings", "warmup_dlcm"))]
  )
  
  warmup_dlcm <- do.call(dependentLCM_fit, warmup_settings)
  
  return(warmup_dlcm)
}


#' Choose class of each observation entirely at random
#' @param mat matrix. Raw data.
#' @param hparams list. List of hyperparameters
#' @keywords internal
getStart_class_random <- function(mat, hparams) {
  class_pi <- hparams$classPi_alpha / sum(hparams$classPi_alpha)
  nobs <- dim(mat)[1]
  classes <- sample.int(hparams$nclass, size=nobs, replace=TRUE, prob=class_pi) - 1
  
  # Ensure every class shows up atleast once:
  ids <- sample.int(n=nobs, size=hparams$nclass, replace=FALSE)
  classes[ids] <- 0:(hparams$nclass-1)
  
  return(list(classes=classes))
}

#' Choose class of each observation using kmodes
#' @param mat matrix. Raw data.
#' @param hparams list. List of hyperparameters
#' @keywords internal
getStart_class_kmodes <- function(mat, hparams, ...) {
  success <- FALSE
  for (i in seq_len(3)) { 
    # Repeated attempts because klaR sometimes errors randomly
    tryCatch({
      kModesFit <- klaR::kmodes(mat, hparams$nclass, ...)
      out <- list(
        classes = kModesFit$cluster - 1
        , centers = kModesFit$modes
      )
      success <- TRUE
      break
    }
    , error=function(cond) {
      message("Warning: Issue with class_init_method = 'kmodes'. Retrying. If issue persists try class_init_method = 'random_centers'.\n Reference getStart_class_kmodes:  ", cond)
    }
    )
  }
  if (success == FALSE) {
    message("Warning: Issue with class_init_method = 'kmodes'. Using class_init_method = 'random_centers' instead.")
    out <- getStart_class_random_centers(mat, hparams, isWeighted=FALSE, isPolar=FALSE)
  }
  
  return(out)
}

#' Choose nclass random centers and put observations in their nearest center
#' Designed so that the chosen centers are as far from each other as possible
#' @param mat matrix. Raw data.
#' @param hparams list. List of hyperparameters
#' @param isWeighted bool. True if more common patterns should be more likely to be chosen as centers. False if all patterns should be equally likely.
#' @param isPolar bool. True to force centers to be as far from each other as possible
#' @param attempts int. Number of tries to get unique centers
#' @keywords internal
getStart_class_random_centers <- function(mat, hparams, isWeighted, isPolar) {
  nobs <- dim(mat)[1]
  nclass <- hparams$nclass
  
  if (isWeighted==FALSE) {
    # All patterns equally likely to be chosen, each pattern appears once
    centers_filter <- which(!duplicated(mat))
  } else {
    # More common patterns more likely to be chosen, patterns appear repeatedly
    centers_filter <- seq_len(nrow(mat))
  }
  
  distance <- matrix(data=NA, nrow=nobs, ncol=nclass)
  if (isPolar==TRUE) {
    # Have centers be as far as possible from each other
    
    # Get centers
    centers <- matrix(data=NA, nrow=nclass, ncol(mat))
    centers[1,] <-mat[sample.integers(x=centers_filter, size=1),]
    distance[,1] <- rowSums(sweep(mat, 2, centers[1,], FUN="!="))
    for (i in seq_len(nclass-1)+1) {
      # choose the center farthest from the current centers
      min_dist <- apply(distance[centers_filter,], 1, min, na.rm=TRUE)
      center_inds <- centers_filter[which(min_dist==max(min_dist))]
      centers[i,] <- mat[sample.integers(x=center_inds, size=1),]
      distance[,i] = rowSums(sweep(mat, 2, centers[i,], FUN="!="))
    }
  } else {
    # Centers entirely at random
    centers <-unique(mat[sample.integers(x=centers_filter, size=nclass),])
    last_processed = 0
    for (j in seq_len(nclass)) {
      
      for (i in last_processed+seq_len(nrow(centers)-last_processed)) {
        distance[,i] = rowSums(sweep(mat, 2, centers[i,], FUN="!="))
      }
      
      last_processed <- nrow(centers)
      if (last_processed==nclass) {
        break # Done. We have nclass unique centers
      }
      
      # Add more centers since we do not have enough
      center_inds <- centers_filter[apply(distance[centers_filter,], 1, min, na.rm=TRUE) > 0]
      
      centers <- unique(rbind(
        centers
        , mat[sample.integers(x=center_inds, size=nclass-nrow(centers)),]
      ))
    }
  }
  
  # Associate observation to nearest class
  classes <- apply(distance, 1, FUN=function(idist) sample.integers(x=which(idist==min(idist)), size=1))
  classes <- classes - 1 # start at 0
  
  return(list(classes=classes, centers=centers))
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
      seq_len(hparams$nitems)
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
  iter <- dlcm$hparams$nitr
  nitems <- dlcm$hparams$nitems
  
  domains_mat <- dlcm$mcmc$domains %>% dplyr::filter(itr == iter)
  domains_mat$items <- (
    # apply(id2pattern(domains_mat$items_id, rep(2, nitems))==1, 2, which)
    unlist(apply(domains_mat[,dlcm$hparams$domain_item_cols] >= 0, 1, function(x) list(which(x)-1)), recursive = FALSE)
  )
  
  domains_list <- (domains_mat 
                   %>% dplyr::group_by(class, domain)
                   %>% dplyr::arrange(pattern_id) 
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
#' @param ndomains Integer. The total number of domains (including empty domains)
#' @param specific_items Boolean. 
#' If FALSE, we look for any domain which produces domains of this size regardless of what specific items they contain.
#' IF TRUE, we fix which items are in which domain, and calculate the probability of grouping these specific items together.
#' @param log Boolean. If TRUE give probability in log scale.
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
#' @param prop float. prop=1 forces singleton domains to be mode. prop>1 makes singleton domains increasingly frequent. Default is prop=2.
#' @export
ndomains_singleton_mode <- function(nitems, prop=2) {
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


#' @name sample.df
#' @title sample.df
#' @description Samples rows from a dataframe or matrix
#' @param x The dataframe or matrix we are sampling from
#' @param size Number of rows to sample
#' @param ... Other options per sample.int()
#' @export
sample.df <- function(x, size, ...) {
  return(x[sample.int(n=nrow(x), size=size, ...),])
}

#' @name sample.integers
#' @title sample.integers
#' @description Samples values from an provided integer vector.
#' Safer version of sample() which has bad behavior if length(x)==1.
#' @param x integer vector we are sampling from
#' @param size Number of rows to sample
#' @param ... Other options per sample.int()
#' @keywords internal
sample.integers <- function(x, size, ...) {
  return(x[sample.int(n=length(x), size=size, ...)])
}

`%>%` <- magrittr::`%>%`

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
#' }
#' @export
dlcm.summary <- function(dlcm, nwarmup=NULL) {
  
  if (is.null(nwarmup)) {
    nwarmup = floor(min(c(1000, dlcm$hparams$nitr/2)))
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
    mode_domains <- unique(dlcm$mcmc$domains %>% dplyr::filter(itr==first_mode_domain_itr, pattern_id==0) %>% .[,c("class2domain", "items_id", "items")]) %>% dplyr::mutate(is_one=1)
    
    rm(domain_items_all_raw)
  }
  
  { # Response probabilities (thetas)
    dlcm$mcmc$domains$item_value <- apply(
      dlcm$mcmc$domains[, dlcm$hparams$domain_item_cols]
      , 1, function(x) paste0(x[x>-1], collapse=", ")
    )
    dlcm$mcmc$domains <- (
      dlcm$mcmc$domains
      %>% dplyr::left_join(mode_domains[,c("class2domain", "items_id", "is_one")], by=c("class2domain", "items_id"))
      %>% dplyr::mutate(is_mode= !is.na(is_one), is_one=NULL) # convert is_one to is_mode
    )
    
    thetas_avg <- (
      dlcm$mcmc$domains 
      %>% dplyr::filter(itr > nwarmup) 
      %>% dplyr::group_by(class, items_id, pattern_id) 
      %>% dplyr::summarize(
        items = dplyr::first(items)
        , item_value = dplyr::first(item_value)
        , is_mode = dplyr::first(is_mode)
        , n=dplyr::n()
        , prob=mean(prob)
        , .groups="keep")
    )
    
    thetas_avg_mode <- reshape2::dcast(
      data=thetas_avg %>% dplyr::filter(is_mode==TRUE)
      , formula = items_id + pattern_id + items + item_value ~ class
      , value.var = "prob"
    )
  }
  
  { # Summarize Classes
    
    classes_cnts <- apply(
      dlcm$mcmc$classes[,-seq_len(nwarmup)], 1
      , function(iclasses, class_vec_default) {
        table_cnts <- table(iclasses)
        out <- class_vec_default
        out[names(table_cnts)] <- table_cnts
        return(out)
      }
      , class_vec_default = setNames(rep(0, dlcm$hparams$nclass), paste0(seq_len(dlcm$hparams$nclass)-1))
    )
    rownames(classes_cnts) <- paste0("class", rownames(classes_cnts))
    
    classes <- apply(classes_cnts, 2, which.max) - 1 # unname(apply(dlcm$mcmc$classes[,-seq_len(nwarmup)], 1, getMode))
  }
  
  # domain_nitems <- table((dlcm$mcmc$domains %>% dplyr::filter(pattern_id==0, itr > nwarmup))$nitems)
  
  domain_accept <- table(dlcm$mcmc$domains_accept[, , -seq_len(nwarmup)])
  
  waic <- dlcm.get_waic(dlcm, itrs=nwarmup:dlcm$hparams$nitr)
  names(waic) <- paste0("waic_", names(waic))
  
  class_pi = rowMeans(dlcm$mcmc$class_pi[,-seq_len(nwarmup), drop=FALSE])
  
  return(c(
    list("thetas_avg"=thetas_avg, "domain_items"=domain_items, "domain_items_all"=domain_items_all, "domain_accept"=domain_accept, "classes"=classes, "class_pi"=class_pi, "thetas_avg_mode"=thetas_avg_mode, "mode_domains"=mode_domains, "first_mode_domain_itr"=first_mode_domain_itr, "classes_cnts"=classes_cnts)
    , waic
  ))
}

#' Calculate likelihood and WAIC
#' @param this_sim Dependent Latent class model
#' @param itrs integer vector. Which iterations should be include in calculation?
#' @export
dlcm.get_waic <- function(this_sim, itrs=NULL) {
  
  alpha <- 0 # No robust MCD variance estimate in WAIC
  
  if (is.null(itrs)) {
    # default in all iterations
    itrs = seq_len(this_sim$hparams$nitr)
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
    %>% dplyr::filter(itr %in% itrs) 
    %>% dplyr::group_by(itr) 
    %>% dplyr::summarize(nparams = dplyr::n() - length(unique(domain)), .groups="keep"))
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
#' \dontrun{
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
    for (iclass in seq_len(nclasses)) {
      ifilter <- classes==iclass
      for (iitem in seq_len(nitems)) {
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
    seq_len(nclass)
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


##############
############## LABEL SWAPPING
##############

#' For each iteration, identify which class labels are swapped and should be relabeled.
#' See identify_swaps() for procedure.
#' @param sim_classes Integer matrix describing what class each observation belongs to in each iteration. There should be one row per observation and one column per iteration as in dlcm$mcmc$classes.
#' @param classes_mode Integer vector describing the most likely (mode) class of each observation.
#' @param maxitr integer. How many attempts to make.
#' @param nclass integer. Number of classes supported by the DLCM
#' @keywords internal
identify_swaps_helper <- function(sim_classes, classes_mode, maxitr, nclass) {
  
  class_levels <- paste0(seq_len(nclass)-1)
  classes_mode <- factor(classes_mode, levels=class_levels)
  counts <- apply(
    sim_classes
    , 2, function(x) table(
      classes_mode
      , factor(x, levels=class_levels)
    )
  )
  
  e1071.matchClasses <- apply(
    counts, 2
    , function(icounts) {
      icounts <- matrix(icounts, nrow=nclass, ncol=nclass)
      as.integer(e1071::matchClasses(icounts, nclass, method="exact", verbose=FALSE))
    }
  )
  rownames(e1071.matchClasses) <- paste0("class_matches", seq_len(nrow(e1071.matchClasses)))
  e1071.matchClasses_id <- factor(apply(e1071.matchClasses-1, 2, paste0, collapse=","))
  identify_out <- data.frame(
    itr = seq_along(e1071.matchClasses_id)
    , valueIsOldClassName_positionIsNewClassIndex_string = e1071.matchClasses_id
    , t(e1071.matchClasses-1)
  )
  
  return(identify_out)
}

#' Applies label swapping based on an output from identify_swaps_helper.
#' @param dlcm An dependent latent class model from dependentLCM_fit().
#' @param identify_out Output from identify_swaps_helper.
#' @param nwarmup integer. Which iterations are warmup iterations?
#' @keywords internal
apply_swaps_helper <- function(
    dlcm # Single chain simulation
    , identify_out # Information on swaps. One row per iteration giving 1) ichain, 2) itr (assumes ordered), 3) valueIsOldClassName_positionIsNewClassIndex_string 4) class_matches
    , nwarmup
) {
  
  if (is.null(dlcm$label_swapping)) {
    dlcm$label_swapping <- list() # initialize
  }
  
  idcols <- c(1,2)
  unique_valueIsOldClassName_positionIsNewClassIndex_strings <- unique(identify_out$valueIsOldClassName_positionIsNewClassIndex_string)
  any_change <- FALSE
  
  for (istr in unique_valueIsOldClassName_positionIsNewClassIndex_strings) {
    
    itrs_swap <- ((identify_out %>% .[,"valueIsOldClassName_positionIsNewClassIndex_string"]) == istr)
    valueIsOldClassName_positionIsNewClassIndex <- unlist(identify_out %>% dplyr::filter(valueIsOldClassName_positionIsNewClassIndex_string==istr) %>% .[1,-idcols]) # setNames(stats_ordered_by_oldclass, nm=originalClasses)[paste0(valueIsOldClassName_positionIsNewClassIndex)] -> stats_ordered_by_newclass
    valueIsOldClassIndex_positionIsNewClassIndex <- valueIsOldClassName_positionIsNewClassIndex+1 # stats_ordered_by_oldclass[valueIsOldClassIndex_positionIsNewClassIndex] -> stats_ordered_by_newclass
    valueIsNewClassName_positionIsOldClassIndex <- order(valueIsOldClassName_positionIsNewClassIndex)-1 # valueIsNewClassName_positionIsOldClassIndex[old_class_id+1] -> new_class_id
    
    # # Examples to understand valueIsOldClassName_positionIsNewClassIndex vs valueIsOldClassIndex_positionIsNewClassIndex vs valueIsNewClassName_positionIsOldClassIndex above
    # class_levels <- (seq_len(dlcm$hparams$nclass)-1)
    # iitr <- which(itrs_swap)[1] # or other iteration applicable to this loop...
    # classes_mode <- factor(dlcm$label_swapping$classes, levels=class_levels)
    # classes_iitr <- factor(dlcm$mcmc$classes[,iitr], levels=class_levels)
    # counts <- table(classes_mode, classes_iitr)
    # counts # original, bad
    # counts[,paste0(valueIsOldClassName_positionIsNewClassIndex)] # fixed
    # counts[,valueIsOldClassIndex_positionIsNewClassIndex] # fixed
    # old_classes <- c(0,0,1,1,2,2)
    # valueIsNewClassName_positionIsOldClassIndex[old_classes+1] # new class IDs for old_classes observations
    
    if (identical(valueIsNewClassName_positionIsOldClassIndex, seq_along(valueIsNewClassName_positionIsOldClassIndex)-1)) {
      next # Nothing to reorder
    } else {
      any_change <- TRUE
    }
    
    dlcm$mcmc$class_pi[,itrs_swap] <- dlcm$mcmc$class_pi[valueIsOldClassIndex_positionIsNewClassIndex,itrs_swap]
    
    dlcm$label_swapping$valueIsNewClassName_positionIsOldClassIndex[, itrs_swap] <- valueIsNewClassName_positionIsOldClassIndex[dlcm$label_swapping$valueIsNewClassName_positionIsOldClassIndex[,itrs_swap]+1]
    
    ifilter <- matrix(
      data = itrs_swap
      , nrow=nrow(dlcm$mcmc$classes)
      , ncol=ncol(dlcm$mcmc$classes)
      , byrow=TRUE
    )
    dlcm$mcmc$classes[ifilter] <- valueIsNewClassName_positionIsOldClassIndex[dlcm$mcmc$classes[ifilter]+1]
    
    ifilter <- dlcm$mcmc$domains$itr %in% which(itrs_swap)
    dlcm$mcmc$domains$class[ifilter] <- valueIsNewClassName_positionIsOldClassIndex[dlcm$mcmc$domains$class[ifilter]+1]
  }
  
  if (any_change==TRUE) {
    dlcm$label_swapping$classes <- unname(apply(dlcm$mcmc$classes[,-seq_len(nwarmup)], 1, getMode))
  }
  
  return(dlcm)
}

#' Identifies and fixes label swapping among classes. 
#' For each iteration, it relabels class to most align with the overall most common (mode) class of each observation.
#' Warning01: Best for homogeneous DLCMs. Potentially fine for heterogenous DLCMs (label swapping allowed, but not informed by domains here). Not appropriate for *partially* heterogeneous DLCMs.
#' Warning02: May not be appropriate if a single class is too small and therefore never shows up as the mode for any observation.
#' Warning03: Relabels warmup iterations which may be unnecessary.
#' @param dlcm An dependent latent class model from dependentLCM_fit().
#' @param nwarmup integer. Which iterations are warmup iterations?
#' @param initial_target_classes integer vector with one value per observation. What set of class labels should we try to mirror, initially? If blank defaults to the most common class for each observation.
#' @param maxitr integer. How many attempts should we make to correct label swapping?
#' @return Returns an updated DLCM with corrected classes. Summary information (e.g. dlcm.summary()) is not modified. 
#' Information on what labels were swapped is given in output$label_swapping.
#' \itemize{
#' \item{"classes"}{= Integer vector. The most common (mode) class of each observation after relabeling. This serves as the target for relabeling (done iteratively).}
#' \item{"valueIsNewClassName_positionIsOldClassIndex"}{= Dataframe with one row per MCMC iteration describing what label swapping was done in that iteration (if any). Each iteration has a vector indexed by the original classes. For each old class index, we provide the new class name.}
#' \item{"valueIsOldClassName_positionIsNewClassIndex"}{= Dataframe with one row per MCMC iteration describing what label swapping was done in that iteration (if any). Each iteration has a vector indexed by the new classes. For each new class index, we provide the old class name.}
#' \item{"any_swaps"}{= Boolean. Did we do any label swapping? Same as nitrs_with_changes>0}
#' \item{"nitrs_attempted"}{= Integer. How many iterations did we run?}
#' \item{"nitrs_with_changes"}{= Integer. In how many iterations did we actually change class labels?}
#' }
#' @export
fix_class_label_switching <- function(dlcm, nwarmup, initial_target_classes=NULL, maxitr=5) {
  
  if (is.null(initial_target_classes)) {
    initial_target_classes <- unname(apply(dlcm$mcmc$classes[,-seq_len(nwarmup)], 1, getMode))
  }
  nclass <- dlcm$hparams$nclass
  dlcm$label_swapping <- list(
    classes = initial_target_classes
    , valueIsNewClassName_positionIsOldClassIndex = matrix(data=seq_len(nclass)-1, nrow=nclass, ncol=ncol(dlcm$mcmc$classes)) # i.e. valueIsNewClassName_positionIsOldClassIndex
  )
  dlcm_copy <- dlcm
  
  # Find and correct for label switching
  nitrs_with_changes <- 0
  for (jitr in seq_len(maxitr)) {
    
    identify_out <- identify_swaps_helper(sim_classes=dlcm_copy$mcmc$classes, classes_mode=dlcm_copy$label_swapping$classes, maxitr=maxitr, nclass=dlcm_copy$hparams$nclass)
    dlcm_copy <- apply_swaps_helper(dlcm=dlcm_copy, identify_out=identify_out, nwarmup = nwarmup)
    
    if (length(unique(identify_out[-seq_len(nwarmup),"valueIsOldClassName_positionIsNewClassIndex_string"])) <= 1) {
      break # Done
    }
    
    nitrs_with_changes <- nitrs_with_changes + 1
  }
  
  identify_out_fn <- identify_out
  if (jitr>1) {
    # Put swaps in terms of the original dlcm
    id_cols <- c(1,2)
    identify_out_fn[,-id_cols] <- t(apply(dlcm_copy$label_swapping$valueIsNewClassName_positionIsOldClassIndex, 2, order))
    identify_out_fn$valueIsOldClassName_positionIsNewClassIndex_string <- apply(identify_out_fn[,-id_cols], 1, paste0, collapse=",")
  }
  
  valueIsNewClassName_positionIsOldClassIndex_df <- data.frame(
    itr=seq_len(ncol(dlcm_copy$label_swapping$valueIsNewClassName_positionIsOldClassIndex))
    , valueIsNewClassName_positionIsOldClassIndex_string = apply(dlcm_copy$label_swapping$valueIsNewClassName_positionIsOldClassIndex, 2, paste0, collapse=",")
    , t(dlcm_copy$label_swapping$valueIsNewClassName_positionIsOldClassIndex)
    , row.names=paste0("itr", seq_len(ncol(dlcm_copy$label_swapping$valueIsNewClassName_positionIsOldClassIndex)))
  )
  dlcm_copy$label_swapping$valueIsOldClassName_positionIsNewClassIndex <- identify_out_fn
  dlcm_copy$label_swapping$any_swaps <- (nitrs_with_changes>0)
  dlcm_copy$label_swapping$nitrs_attempted <- jitr
  dlcm_copy$label_swapping$nitrs_with_changes <- nitrs_with_changes
  dlcm_copy$label_swapping$valueIsNewClassName_positionIsOldClassIndex <- valueIsNewClassName_positionIsOldClassIndex_df
  dlcm_copy$mcmc$domains$class_original <- dlcm$mcmc$domains$class
  dlcm_copy$mcmc$classes_original <- dlcm$mcmc$classes
  return(dlcm_copy)
}

