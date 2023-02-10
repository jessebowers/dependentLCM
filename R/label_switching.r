#' @include utilities.r
NULL

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
#' @return Returns an updated DLCM with corrected classes. Only modifies dlcm$mcmc[c('domains', 'classes', 'class_pi')]. Other information is not modified, including summary information (e.g. dlcm.summary()).
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
  if (dlcm_copy$label_swapping$any_swaps==TRUE) {
    dlcm_copy$mcmc$domains$class_original <- dlcm$mcmc$domains$class
    dlcm_copy$mcmc$classes_original <- dlcm$mcmc$classes
  }
  return(dlcm_copy)
}

