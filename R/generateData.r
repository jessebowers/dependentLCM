#' @include utilities.r
NULL

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
  
  
  if ("matrix" %in% class(thetas)) {
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
