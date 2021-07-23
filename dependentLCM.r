library(klaR)
library(Rcpp)
library(abind)
library(dplyr)
library(tidyverse)

CLASSPI_ALPHA_DEFAULT = 2
DOMAIN_ALPHA_RATIO = 2
THETA_ALPHA_DEFAULT = 2
DOMAIN_PROPOSAL_RATIO_DEFAULT = 0.5

dependentLCM_fit <- function(nitr, ...) {
  
  all_params <- get_start(...)
  
  if (check_params(all_params)) {
    stop()
  }
  
  out = dependentLCM_fit_cpp(x=(all_params$mat),hparams_list=all_params$hparams, params_list=all_params$bayesparams, nitr=nitr)
  
  # Clean up output
  out$thetas_id <- do.call(cbind, out$thetas_id)
  rownames(out$thetas_id) <- c("itr", "class", "domain")
  out$thetas_patterns <- do.call(cbind, out$thetas_patterns)
  rownames(out$thetas_patterns) <- paste0("item_", 0:(nrow(out$thetas_patterns)-1))
  out$thetas_probs <- unlist(out$thetas_probs)
  
  # Supplemental
  thetas_nitems <- colSums(out$thetas_patterns > -1)
  thetas_class2domain <- all_params$hparams$class2domain[out$thetas_id["class",,drop=TRUE]+1]
  thetas_domain_id <- as.data.frame(t(rbind(out$thetas_id, thetas_class2domain)))
  thetas_domain_id <- do.call(paste, thetas_domain_id) # paste rows together, faster than apply
  thetas_domain_filter <- ave(seq_len(ncol(out$thetas_id)), thetas_domain_id, FUN = rank)
  thetas_domain_filter <- (thetas_domain_filter == 1)
  
  # Merge thetas attributes
  out$thetas <- data.frame(
    t(out$thetas_id)
    , class2domain = thetas_class2domain
    , prob=out$thetas_probs
    , nitems = thetas_nitems
    , domain_filter = thetas_domain_filter
    , t(out$thetas_patterns) 
  )
  theta_item_cols <- grep("^item_[0-9]+", colnames(out$thetas))
  out$thetas_accept <- do.call(function(...) abind(..., along=3), out$thetas_accept)
  
  out$thetas_id <- NULL
  out$thetas_patterns <- NULL
  out$thetas_probs <- NULL
  out$nclass2domain <- NULL
  
  return(c(out, all_params$hparams, list(theta_item_cols=theta_item_cols)))
}

get_start <- function(
  df=NULL, mat=NULL
  # Hyperparameters
  ,nclass=2, ndomains=NULL, class2domain=NULL, classPi_alpha=CLASSPI_ALPHA_DEFAULT, domain_alpha=NULL, domain_nitems=NULL, theta_alpha=THETA_ALPHA_DEFAULT, domain_proposal_ratio=DOMAIN_PROPOSAL_RATIO_DEFAULT
  # Bayes parameters
  , class_pi = NULL, classes = NULL,  domains_pi = NULL, thetas = NULL
  # Misc
  , class_init_method = "kmodes") {
  
  if (is.null(mat)) {
    mat <- get_start.matrix(df)
  }
  
  hparams <- get_start.hparams(
    df=mat
    # Hyperparameters
    ,nclass=nclass, ndomains=ndomains, class2domain=class2domain, classPi_alpha=classPi_alpha, domain_alpha=domain_alpha, domain_nitems=domain_nitems, theta_alpha=theta_alpha, domain_proposal_ratio=domain_proposal_ratio
  )
  
  bayesparams <- get_start.bayes_params(
    mat=mat, hparams=hparams
    , class_pi = class_pi, classes = classes,  domains_pi = domains_pi, thetas = thetas
    # Misc
    , class_init_method = class_init_method
  )
  
  return(list(
    mat = mat
    , hparams = hparams
    , bayesparams = bayesparams
  ))
}

get_start.hparams <- function(
  df=NULL, nitems=NULL
  # Hyperparameters
  ,nclass=2, ndomains=NULL, class2domain=NULL, classPi_alpha=CLASSPI_ALPHA_DEFAULT, domain_alpha=NULL, domain_nitems=NULL, theta_alpha=THETA_ALPHA_DEFAULT, domain_proposal_ratio=DOMAIN_PROPOSAL_RATIO_DEFAULT
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
    domain_alpha <- nitems/ndomains * DOMAIN_ALPHA_RATIO
  }
  if (length(domain_alpha)==1) {
    domain_alpha <- rep(domain_alpha, ndomains)
  }
  if (is.null(dim(domain_alpha))) {
    # Convert vector to matrix of appropriate size
    domain_alpha <- do.call(rbind, rep(list(domain_alpha), nclass2domain))
  }
  
  # domain_nitems
  if (is.null(domain_nitems)) {
    domain_nitems <- nitems # No restrictions
  }
  
  # theta_alpha, no action taken
  
  return(list(
    nclass=nclass, ndomains=ndomains, class2domain=class2domain, classPi_alpha=classPi_alpha, domain_alpha=domain_alpha, domain_nitems=domain_nitems, theta_alpha=theta_alpha, nitems = nitems, item_nlevels = item_nlevels, nclass2domain = nclass2domain, domain_proposal_ratio=domain_proposal_ratio
  ))
}

get_start.matrix <- function(df) {
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

get_start.bayes_params <- function(
  mat, hparams
  , class_pi = NULL, classes = NULL,  domains_pi = NULL, thetas = NULL
  # Misc
  , class_init_method = "kmodes"
) {
  
  # domains_pi
  if (is.null(domains_pi)) {
    domains_pi = hparams$domain_alpha / sum(hparams$domain_alpha)
  }
  
  # classes
  if (is.null(classes)) {
    if (class_init_method=="kmodes") {
      classes <- get_start.class_kmodes(mat, hparams)
    }
    if (class_init_method=="random") {
      classes <- get_start.class_random(mat, hparams)
    }
  }
  
  # class_pi
  class_pi <- table(classes)+hparams$classPi_alpha
  class_pi <- class_pi / sum(class_pi)
  
  # thetas
  if (is.null(thetas)) {
    thetas <- get_start.thetas(mat, classes, hparams)
  } # Future maybe allow other types of inputs
  
  return(list(
    class_pi = class_pi, classes = classes,  domains_pi = domains_pi, thetas = thetas
  ))
}



get_start.class_random <- function(mat, hparams) {
  class_pi <- hparams$classPi_alpha / sum(hparams$classPi_alpha)
  nobs <- dim(mat)[1]
  classes <- sample(0:(hparams$nclass-1), size=nobs, replace=TRUE, prob=class_pi)
  
  # Ensure every class shows up atleast once:
  ids <- sample.int(nobs, hparams$nclass)
  classes[ids] <- 0:(hparams$nclass-1)
  
  return(classes)
}

get_start.class_kmodes <- function(mat, hparams, iter.max=2, ...) {
  return(kmodes(mat, hparams$nclass, iter.max=iter.max, ...)$cluster-1)
}

get_start.thetas <- function(mat, classes, hparams) {
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
############## SIMULATION
##############

generate <- function(n, pis, thetas) {
  
  classes <- apply(rmultinom(n, 1, pis), 2, function(x) which(x==1))
  
  responses <- sapply(
    classes
    , function(iclass) {runif(nrow(thetas)) > thetas[, iclass]}
  )+0
  
  posterior <- lapply(
    1:ncol(thetas)
    , function(iclass){
      apply(responses, 2
            , function(iresp) {
              prod(c(thetas[which(iresp==1), iclass], (1-thetas)[which(iresp==0), iclass]))
            }
      )
    }
  )
  posterior <- do.call(cbind, posterior)
  posterior <- sweep(posterior, 2, pis, "*")
  posterior <- posterior / rowSums(posterior)
  
  return(list(
    classes = classes
    , responses = t(responses)
  ))
}

dlcm.thetas_items <- function(dlcm) {
  return(data.frame(
    items = apply(
      dlcm$thetas[, dlcm$theta_item_cols] > -1
      , 1, function(x) paste0(which(x)-1, collapse=", ")
    )
    , item_value = apply(
      dlcm$thetas[, dlcm$theta_item_cols]
      , 1, function(x) paste0(x[x>-1], collapse=", ")
    )
    , stringsAsFactors = FALSE
  ))
}

get_perc <- function(x) {
  x / sum(x)
}
getmode <- function(x) {
  xcounts <- table(x)
  mode_name <- names(xcounts)[which.max(xcounts)]
  return(mode_name)
}

dlcm.summary <- function(out, nwarmup=1000) {
  out$thetas[, c("items", "item_value")] <- dlcm.thetas_items(out)
  
  thetas_avg <- out$thetas %>% filter(itr > nwarmup) %>% group_by(class, items, item_value) %>% summarize(n=n(), prob=mean(prob))
  domain_items <- out$thetas %>% filter(domain_filter==TRUE, itr > nwarmup) %>% group_by(items) %>% summarize(nitems=max(nitems), n=n())
  domain_nitems <- table((out$thetas %>% filter(itr > nwarmup, domain_filter==TRUE))$nitems)
  domain_accept <- apply(
    out$thetas_accept[, , -(1:nwarmup)], 1
    , function(x) list(table(x))
  )
  domain_accept <- do.call(bind_rows, domain_accept)
  
  return(list(
    "thetas_avg"=thetas_avg, "domain_items"=domain_items, "domain_nitems"=domain_nitems, "domain_accept"=domain_accept
  ))
}