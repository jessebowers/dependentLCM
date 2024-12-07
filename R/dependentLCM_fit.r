#' @include utilities.r
NULL

#' dependentLCM: Domain Latent Class Model
#' @description 
#' Latent Class Models (LCMs) are used to cluster multivariate categorical data (e.g. group participants based on survey responses). Traditional LCMs assume a property called conditional independence. This assumption can be patternadjusted, leading to model misspecification and overparameterization. To combat this problem, we developed a novel Bayesian model called a Domain Latent Class Model (DLCM), which permits conditional dependence. Compared to traditional LCMs, DLCMs are effective in applications with time series, overlapping items, and structural zeroes.
#' 
#' The primary function is dependentLCM_fit. 
#' 
#' Jesse Bowers. Steve Culpepper. "Domain Latent Class Models." Bayesian Anal. Advance Publication 1 - 28, 2024. \url{https://doi.org/10.1214/24-BA1433}
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
DOMAIN_THETA_PRIOR_TYPE = c("bucket", "patternadjusted", "niave")
CLASS_INIT_METHODS = c("random_centers_polar", "random_centers", "random", "kmodes")
DOMAIN_PROPOSAL_EMPTY = 0.3
STEPS_ACTIVE = c("thetas"=TRUE, "domains"=TRUE, "class_pi"=TRUE, "classes"=TRUE, "identifiable"=TRUE, "likelihood"=TRUE, "class_collapse"=FALSE)
SAVE_ITRS <- c(domains_accept=0, class_loglik=0, class_loglik_collapsed=0, agg_loglik=Inf, classes=1, class_counts=Inf, all=Inf)
DOMAIN_MAXITEMS = 10
CLASS2DOMAIN_FUNS = list(
  "HOMO" = function(nclass) rep(0, nclass)
  , "HET" = function(nclass) seq(0, nclass-1)
  , "TRAD" = function(nclass) rep(0, nclass) # this term essentially ignored since domains are off
)
CLASS2DOMAINS = names(CLASS2DOMAIN_FUNS)

#' Fits a bayesian Domain LCM model
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
#' \item{"domains"}{=Dataframe with one row for iteration X domain X pattern. Each domain (items_id) indicates a group of items, and each pattern (pattern_id) indicates the response pattern of those items. To look up information on these join to dataframe response_patterns on columns items_id and pattern_id. For each MCMC iteration, the 'domains' dataframe gives the response probability for each class and response pattern. It contains the following columns:
#' \itemize{
#' \item{"itr"}{=What MCMC iteration is this?}
#' \item{"class"}{=What class are we calculating the response probabilities for?}
#' \item{"domain"}{=ID. What domain is this? For most purposes you should use items_id instead.}
#' \item{"pattern_id"}{=Integer uniquely identifying a response pattern to thie items in this domain. We will calculate the probability of this response pattern. pattern_id needs to be paired with a items_id to make sense of it.}
#' \item{"items_id"}{=Integer uniquely identifying what items (what set of items) are in this domain.}
#' \item{"class2domain"}{=For heterogeneous and partially heterogeneous DLCMs, this identifies which group of latent classes this domain belongs to.}
#' \item{"prob"}{=What's the probability of getting this response pattern?}
#' }}
#' \item{"response_patterns"}{=Dataframe with one row for each domain X pattern. This identifies the items in the domain, and the vaues of the respone pattern. It has the following columns:
#' \itemize{
#' \item{"items_id"}{=Integer uniquely identifying what items (what set of items) are in this domain.}
#' \item{"pattern_id"}{=Integer uniquely identifying a response pattern to thie items in this domain. We will calculate the probability of this response pattern. pattern_id needs to be paired with a items_id to make sense of it.}
#' \item{"items"}{=Vector containing the items in this domain. Function of items_id.}
#' \item{"pattern"}{=Vector containing the values in the response pattern. Function of with (items_id,pattern_id).}
#' \item{"items_id_first"}{=Boolean. Filtering on TRUE gives unique domains/items_id.}
#' \item{"nitems"}{=How many items are in this domain?}
#' \item{"item_#"}{=For each item #, gives the specific value of that item in this response pattern. A value of -1 indicates this item is not in this domain. item_# is a function of (items_id, pattern_id).}
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
#'   , save_itrs=c(all=6000-1000) # warmup of 1000 used here
#'   , df=xdf
#'   , nclass=3
#'   , class2domain = "HOMO"
#' )
#' dlcm$summary <- dlcm.summary(dlcm)
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
    ,nclass=NCLASS, ndomains=NULL, class2domain=CLASS2DOMAINS[1], classPi_alpha=CLASSPI_ALPHA, domain_maxitems=NULL, theta_alpha=THETA_ALPHA, domain_proposal_empty=DOMAIN_PROPOSAL_EMPTY, domain_nproposals=NULL, domain_theta_prior_type=DOMAIN_THETA_PRIOR_TYPE[1]
    # Bayes parameters
    , class_pi = NULL, classes = NULL, domains = NULL
    # Misc
    , steps_active = STEPS_ACTIVE, save_itrs=SAVE_ITRS
    , class_init_method = CLASS_INIT_METHODS[1]
    , warmup_settings = "default", warmup_dlcm=NULL
    , domainPriorKnowledgeLog=NULL, print_itrs=Inf
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
    ,nclass=nclass, ndomains=ndomains, class2domain=class2domain, classPi_alpha=classPi_alpha, domain_maxitems=domain_maxitems, theta_alpha=theta_alpha, domain_proposal_empty=domain_proposal_empty, domain_nproposals=domain_nproposals, steps_active=steps_active, save_itrs=save_itrs, domain_theta_prior_type=domain_theta_prior_type, domainPriorKnowledgeLog=domainPriorKnowledgeLog, print_itrs=print_itrs
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
  mcmc = dependentLCM_fit_cpp(x=(all_params$mat),hparams_list=all_params$hparams, params_list=all_params$bayesparams)
  datetimes <- c(datetimes, "mcmc_end"=Sys.time())
  
  #
  # Post Processing
  #
  
  #
  # Assemble mcmc$domains
  #
  
  # Clean up output
  mcmc$domains_id <- do.call(cbind, mcmc$domains_id)
  rownames(mcmc$domains_id) <- c("itr", "class", "domain", "pattern_id", "items_id")
  mode(mcmc$domains_id) <- "integer"
  mcmc$domains_lprobs <- unlist(mcmc$domains_lprobs)
  
  # response_patterns
  domains_id_unique_ids <- which(!duplicated(mcmc$domains_id[c("items_id", "pattern_id"), ], MARGIN=2, fromLast=TRUE))
  response_patterns <- as.data.frame(t(mcmc$domains_id[c("items_id", "pattern_id"), domains_id_unique_ids]))
  {
    idomains_patterns <- itemid2patterns(response_patterns$pattern_id, response_patterns$items_id, hparams[["item_nlevels"]])
    rownames(idomains_patterns) <- paste0("item_", seq_len(nrow(idomains_patterns)))
    mode(idomains_patterns) <- "integer"
    all_params$hparams$domain_item_cols <- rownames(idomains_patterns)
    response_patterns <- cbind.data.frame(
      response_patterns
      , t(idomains_patterns)
    )
  }
  response_patterns$items <- apply(
    response_patterns[,all_params$hparams$domain_item_cols]
    , 2
    , function(x) which(x>-1)
    , simplify = FALSE
  )
  response_patterns$pattern <- apply(
    response_patterns[,all_params$hparams$domain_item_cols]
    , 2, function(x) x[x>-1]
    , simplify = FALSE
  )
  response_patterns$items_id_first <- (response_patterns$pattern_id==0)
  response_patterns$nitems <- sapply(response_patterns$items, length)
  mcmc$response_patterns <- response_patterns
  rm(response_patterns)
  
  # Supplemental
  domains_class2domain <- all_params$hparams$class2domain[mcmc$domains_id["class",,drop=TRUE]+1]
  mcmc$domains_id["itr",] <- mcmc$domains_id["itr",] + 1L # start at 1
  mcmc$.Random.seed <- .Random.seed_start
  
  # Merge domains attributes
  mcmc$domains <- data.frame(
    t(mcmc$domains_id)
    , class2domain = domains_class2domain
    , prob=exp(mcmc$domains_lprobs)
    , stringsAsFactors = FALSE
  )
  
  # Delete redundant info
  mcmc$domains_id <- NULL
  mcmc$domains_lprobs <- NULL
  mcmc$nclass2domain <- NULL
  
  #
  # Assemble other
  #
  
  mcmc$domains_merged <- as.data.frame(
    mcmc$domains 
    %>% dplyr::filter(pattern_id==0) 
    %>% dplyr::group_by(itr, class2domain)
    %>% dplyr::filter(class == min(class))
    %>% dplyr::arrange(-nitems, items)
    %>% dplyr::summarize(domains_merged=paste(items, collapse="|"), .groups="keep")
  )
  
  if (mcmc$nitrLik > 0) {
    names(mcmc$itrLogLik) <- paste0("itr", seq(to=hparams$nitr, length.out=mcmc$nitrLik))
    names(mcmc$obsLik) <- paste0("obs", seq_along(mcmc$obsLik))
    names(mcmc$obsLogLik) <- paste0("obs", seq_along(mcmc$obsLogLik))
    names(mcmc$obsLogLik2) <- paste0("obs", seq_along(mcmc$obsLogLik2))
  }
  if (length(mcmc$domains_accept)>0) {
    mcmc$domains_accept <- do.call(function(...) abind::abind(..., along=3), mcmc$domains_accept)
    mcmc$domains_accept <- set_dimnames(mcmc$domains_accept, c("proposal", "class2domain", NA))
    dimnames(mcmc$domains_accept)[[3]] <- paste0("itr", seq(to=hparams$nitr, length.out=dim(mcmc$domains_accept)[3]))
  }
  if (length(mcmc$class_loglik)>0) {
    mcmc$class_loglik <- do.call(function(...) abind::abind(..., along=3), mcmc$class_loglik)
    mcmc$class_loglik <- set_dimnames(mcmc$class_loglik, c("class", "obs", NA))
    dimnames(mcmc$class_loglik)[[3]] <- paste0("itr", seq(to=hparams$nitr, length.out=dim(mcmc$class_loglik)[3]))
  }
  if (length(mcmc$class_loglik_collapsed)>0) {
    mcmc$class_loglik_collapsed <- do.call(function(...) abind::abind(..., along=3), mcmc$class_loglik_collapsed)
    mcmc$class_loglik_collapsed <- set_dimnames(mcmc$class_loglik_collapsed, c("class", "obs", NA))
    dimnames(mcmc$class_loglik_collapsed)[[3]] <- paste0("itr", seq(to=hparams$nitr, length.out=dim(mcmc$class_loglik_collapsed)[3]))
  }
  
  # name
  mcmc$classes <- set_dimnames(mcmc$classes, c("obs", "itr"))
  mcmc$class_pi <- set_dimnames(mcmc$class_pi, c("class", "itr"))
  mcmc$class_counts <- set_dimnames(mcmc$class_counts, c("class", "obs"))
  
  
  datetimes <- c(datetimes, "fun_end"=Sys.time())
  mcmc$runtimes <- as.numeric(c(diff(datetimes), total=tail(datetimes,1)-datetimes[1]))
  names(mcmc$runtimes) <- c("pre", "mcmc", "post", "total")
  
  return(list(hparams=all_params$hparams, mcmc=mcmc, warmup=warmup_dlcm))
}


#' Generate hyperparameters
#' Generate a list of hyperparameters, adding default values when necessary
#' @param nitr integer. Number of iterations to run the bayes MCMC
#' @param df dataframe. The data you wish to analyze. Should describe nominal data. Factor columns work best, but we make an honest effort to interpret any values given (see getStart_matrix()). This will be converted into an integer matrix (see 'matrix' argument). This 'matrix' argument takes precedence over the 'df' argument; only one of these two arguments is needed. Assumes no missing data (NA rows will be dropped).
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
#' \itemize{
#' \item{"thetas=TRUE"}{ to do gibbs on response probabilities}
#' \item{"domains=TRUE"}{ to do metropolis on domains}
#' \item{"class_pi=TRUE"}{ to do gibbs on class prior}
#' \item{"classes=TRUE"}{ to do gibbs on class membership}
#' \item{"identifiable=TRUE"}{ to check generic identifiability conditions of domains}
#' \item{"likelihood=TRUE"}{ to get the likelihood that each observation is in each class}
#' \item{"class_collapse=TRUE"}{ to collapse on theta when sampling classes}
#' }
#' @param save_itrs Named numeric vector. Gives the maximum number of MCMC iterations which should be saved. Can set components to infinity to get data from every iteration.
#' \itemize{
#' \item{"domains_accept=#}{ for saving domain metropolis accept/reject choices.}
#' \item{"class_loglik=#}{ saves the probability that a response is observed conditional on each class}
#' \item{"class_loglik_collapsed=#}{ as class_loglik but after collapsing on response probabilities}
#' \item{"agg_loglik=#}{ saves aggregate likelihood of each iteration}
#' \item{"classes=#}{ saves the class vector of each iteration.}
#' \item{"class_counts=#}{ saves the number of times an observation is in each class.}
#' \item{"all=#}{ saves everything. If lower, takes precidence over other terms.}
#' }
#' @param domain_theta_prior_type string. Defines what sort of prior is used for domains and theta. One of the following.
#' \itemize{
#' \item{"bucket=}{ recommended/default. has moderate regularization on domains. Domain uses "balls in buckets" prior, and theta uses dirichlet prior.}
#' \item{"patternadjusted=}{ has strong regularization on domains. As bucket, but domain prior is adjusted further to cancel out any theta prior.}
#' \item{"niave=}{ no regularization on domains. Bad outcomes result. For demonstration purposes only. Assumes all domains are equally likely.}
#' }
#' @param domainPriorKnowledgeLog Numeric. An nitem*nitem upper triangular matrix. Values of zero indicate no prior knowledge of the domain structure. Values greater than zero indicate that this pair of items should be more likely to be in the same domain. Values less than zero indicate that this pair of items should be less likely to be in the same domain.
#' @param print_itrs Integer. Optional. Prints progress by printing every print_itrs iterations.
#' @keywords internal
getStart_hparams <- function(
    nitr, df=NULL, nitems=NULL
    # Hyperparameters
    ,nclass=NCLASS, ndomains=NULL, class2domain=NULL, classPi_alpha=CLASSPI_ALPHA, domain_maxitems=NULL, theta_alpha=THETA_ALPHA, domain_proposal_empty=DOMAIN_PROPOSAL_EMPTY, domain_nproposals=NULL, steps_active=STEPS_ACTIVE, save_itrs=SAVE_ITRS, domain_theta_prior_type=DOMAIN_THETA_PRIOR_TYPE[1], domainPriorKnowledgeLog=NULL, print_itrs=Inf
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
  domains_off_override = FALSE
  if ((typeof(class2domain) == "character") & length(class2domain)==1) {
    class2domain <- toupper(class2domain)
    if (class2domain=="TRAD") {
      # Do traditional LCM without domains
      domains_off_override <- TRUE
    }
    
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
  if (steps_active_fn["classes"]==FALSE) {
    # if not doing classes, cannot collapse on classes
    steps_active_fn["class_collapse"] <- FALSE
  }
  if (steps_active_fn["classes"] & !steps_active_fn["class_collapse"]) {
    # if doing un-collapsed classes, classLikelihood will be calculated
    steps_active_fn["likelihood"] <- TRUE
  }
  if (domains_off_override) {
    steps_active_fn["domains"] <- FALSE
  }
  
  # steps_active (fill in missing values)
  save_itrs_fn = SAVE_ITRS
  save_itrs_fn[names(save_itrs)] = save_itrs
  save_itrs_fn[save_itrs_fn>nitr] = nitr
  save_itrs_fn[save_itrs_fn>save_itrs_fn["all"]] <- save_itrs_fn["all"]
  if (steps_active_fn["class_collapse"]==FALSE) {
    save_itrs_fn["class_loglik_collapsed"] <- 0
  }
  if (steps_active_fn["likelihood"]==FALSE) {
    save_itrs_fn["class_loglik"] <- 0
    save_itrs_fn["agg_loglik"] <- 0
  }
  
  if (identical(domainPriorKnowledgeLog, NULL)) {
    domainPriorKnowledgeLog = matrix(data=0, nrow=nitems, ncol=nitems)
  }
  
  return(list(
    nitr=nitr, nclass=nclass, ndomains=ndomains, class2domain=class2domain, classPi_alpha=classPi_alpha, domain_maxitems=domain_maxitems, theta_alpha=theta_alpha, nitems = nitems, item_nlevels = item_nlevels, nclass2domain = nclass2domain, domain_proposal_empty=domain_proposal_empty, domain_nproposals=domain_nproposals, steps_active = steps_active_fn, save_itrs=save_itrs_fn
    , domain_theta_prior_type = domain_theta_prior_type, domainPriorKnowledgeLog=domainPriorKnowledgeLog, print_itrs=print_itrs
  ))
}

#' Converts a dataframe of factors to matrix of integers (indexed starting at 0).
#' If dependentLCM_fit() is executed with the df argument, the mat argument is generated using this function.
#' @param df dataframe. The data you wish to analyze. Should describe nominal data. Factor columns work best, but we make an honest effort to interpret any values given.
#' @export
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
#' @param mat matrix. The data you wish to analyze. Should be an integer matrix describing nominal data. Each column should take values 0, 1, 2, ..., max-for-this-column. Assumptions: i) In each column the maximum value is achieved for atleast one observation. ii) There are no NAs.
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
      
      warmup_settings$save_itrs <- SAVE_ITRS
      warmup_settings$save_itrs[SAVE_ITRS>0] <- 1 # only archive last iteration warmup by default
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
  
  if (!setequal(names(all_params$hparams$save_itrs), names(SAVE_ITRS))) {
    warning("save_itrs invalid")
    is_problem = TRUE
  }
  
  if (max(all_params$bayesparams$classes) > all_params$hparams$nclass) {
    warning("classes > nclass")
    is_problem = TRUE
  }
  
  if (length(all_params$bayesparams$domains) != all_params$hparams$nclass) {
    warning("domains should have one element per class")
    is_problem = TRUE
  }
  
  domain_items <- lapply(
    all_params$bayesparams$domains
    , function(iclass_domains) sort(unname(unlist(lapply(iclass_domains, function(jdomain) jdomain$items))))
  )
  if (!all(sapply(domain_items, function(items) identical(items, seq_len(all_params$hparams$nitems)-1)))) {
    warning("items in domains do not match items in data")
    is_problem = TRUE
  }
  
  if (!(all_params$hparams$domain_theta_prior_type %in% DOMAIN_THETA_PRIOR_TYPE)) {
    warning("domain_theta_prior_type not valid")
    is_problem = TRUE
  }
  
  if ((all_params$hparams$domain_theta_prior_type == "patternadjusted") & (all_params$hparams$theta_alpha!=1)) {
    warning("domain_theta_prior_type=patternadjusted requires theta_alpha=1")
    is_problem = TRUE
  }
  
  domainPriorKnowledgeLog = all_params$hparams$domainPriorKnowledgeLog
  if (nrow(domainPriorKnowledgeLog) != all_params$hparams$nitems) {
    warning("domainPriorKnowledgeLog should be an nitems*nitems matrix")
    is_problem = TRUE
  } else if (ncol(domainPriorKnowledgeLog) != all_params$hparams$nitems) {
    warning("domainPriorKnowledgeLog should be an nitems*nitems matrix")
    is_problem = TRUE
  } else if (sum(domainPriorKnowledgeLog[lower.tri(domainPriorKnowledgeLog, diag=TRUE)] != 0, na.rm = TRUE)>0) {
    warning("domainPriorKnowledgeLog should be zero except in upper triangle")
    is_problem = TRUE
  } else if (max(domainPriorKnowledgeLog, na.rm = TRUE) > log(all_params$hparams$ndomains - all_params$hparams$nitems + 1)) {
    warning("domainPriorKnowledgeLog larger than recommended")
    # do not set is_problem
  }
  
  return(is_problem)
}

#' Take domain latent class model (dlcm) output (i.e. from dependentLCM_fit)
#' and convert it into Hyper/Bayes parameter arguments to put into dependentLCM_fit.
#' Used namely for nesting dependentLCM_fit.
#' @param dlcm Domain latent class model from dependentLCM_fit()
#' @param iter_diff Zero indicates last iteration. One indicates the previous iteration. Etc.
#' @keywords internal
dlcm2paramargs <- function(dlcm, iter_diff=0) {
  nitems <- dlcm$hparams$nitems
  
  domains_mat <- dlcm$mcmc$domains %>% dplyr::filter(itr == dlcm$hparams$nitr-iter_diff)
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
    class_pi = unname(dlcm$mcmc$class_pi[,dim(dlcm$mcmc$class_pi)[2]]-iter_diff)
    , classes = unname(dlcm$mcmc$classes[,dim(dlcm$mcmc$classes)[2]]-iter_diff)
    , domains =  domains_list
  )
  
  all_params <- list(
    hparams = dlcm$hparams
    , bayesparams = bayesparams
  )
  
  return(all_params)
}

#' How many domains (ndomains) we we need before 'all singleton domains' are prop-times more likely than 'one 2-item domain with rest singletons'? Using the prior only.
#' Solves: ldomain_prior(rep(1, nitems), ndomains, specific_items=FALSE, log=TRUE) - ldomain_prior(c(2,rep(1, nitems-2)), ndomains, specific_items=FALSE, log=TRUE) = log(prop)
#' @param nitems integer. Number of items in the data.
#' @param prop float. prop=1 forces singleton domains to be mode. prop>1 makes singleton domains increasingly frequent. Default is prop=2.
#' @keywords internal
ndomains_singleton_mode <- function(nitems, prop=2) {
  ceiling(nitems+prop*nitems*(nitems-1)/2-1)
}

