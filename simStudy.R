sink.reset <- function(){
  for(i in seq_len(sink.number())){
    sink(NULL)
  }
}
sink.reset()
outfile_path <- "C:/Users/drakethrice/Desktop/analysis-output.txt" # ; close(file(outfile_path, open="w")); sink(outfile_path)

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

generate2 <- function(n, pis, thetas_list) {
  nclasses <- length(pis)
  nitems <- length(thetas_list[[1]])
  
  classes <- apply(rmultinom(n, 1, pis), 2, function(x) which(x==1))
  
  responses <- matrix(NA, nrow=n, ncol=nitems)
  for (iclass in 1:nclasses) {
    ifilter <- classes==iclass
    for (iitem in 1:nitems) {
      responses[ifilter, iitem] = apply(
        rmultinom(sum(ifilter), 1, thetas_list[[iclass]][[iitem]])
        , 2, function(x) which(x==1))
    }
  }
  
  return(list(
    classes=classes
    , responses=responses
  ))
}

sim_list <- generate2(n=1000
          , pis=c(0.5, 0.5)
          , thetas_list = list(
            c( rep(list(c(0.8,0.2)), 20)
               , rep(list(c(0.6,0.3,0.05,0.05)), 2)
               , rep(list(c(rep(0.3,3), rep(0.1/5,5))), 1)
            )
            , c( rep(list(c(0.5,0.5)), 20)
                 , rep(list(c(0.3,0.6,0.05,0.05)), 2)
                 , rep(list(c(rep(0.4,2), rep(0.2/6,6))), 1)
            )
          ))
pattern_lookups <- list(
  items2 = expand.grid(c(0,1), c(0,1))[c(4,1,2,3), ]
  , items3 = expand.grid(c(0,1), c(0,1), c(0,1))[c(8,1,2:7), ]
)
sim_list$repsonses2 <- cbind(
  sim_list$responses[, 1:20]
  , pattern_lookups[["items2"]][sim_list$responses[, 21], ]
  , pattern_lookups[["items2"]][sim_list$responses[, 22], ]
  , pattern_lookups[["items3"]][sim_list$responses[, 23], ]
)



# system.time(
#   {
#     nitr <- 6000
#     close(file(outfile_path, open="w")); sink(outfile_path)
#     rm(out)
#     set.seed(4)
#     sim_list <- generate(
#       n = 1000
#       , pis = c(0.5, 0.5)
#       , thetas = cbind(
#         c(rep(0.2, 10), rep(0.8, 10))
#         , c(rep(0.8, 10), rep(0.2, 10))
#       ))
#     set.seed(5)
#     out <- dependentLCM_fit(df = sim_list$responses, nclass=2, class2domain = c(0,0), nitr= nitr, class_init_method = "kmodes", domain_proposal_ratio=0.5)
#     out_summary <- dlcm.summary(out)
#     sink.reset()
#   }
# )

# system.time(
#   {
#     nitr <- 6000
#     close(file(outfile_path, open="w")); sink(outfile_path)
#     rm(out)
#     set.seed(4)
#     sim_list <- generate(
#       n = 1000
#       , pis = c(1)
#       , thetas = cbind(
#         c(rep(0.2, 10), rep(0.8, 10))
#       ))
#     set.seed(5)
#     out <- dependentLCM_fit(df = sim_list$responses, nclass=1, class2domain = c(0), nitr= nitr, class_init_method = "kmodes", domain_proposal_ratio=0.5)
#     out_summary <- dlcm.summary(out)
#     sink.reset()
#   }
# )

# system.time(
#   {
#     nitr <- 100
#     close(file(outfile_path, open="w")); sink(outfile_path)
#     rm(out)
#     set.seed(4)
#     sim_list <- generate(
#       n = 200 # 1000
#       , pis = c(1)
#       , thetas = cbind(
#         c(rep(0.2, 1), rep(0.8, 1), rep(0.9,20))
#       ))
#     sim_list$responses2 <- cbind(
#       do.call(cbind, rep(list(sim_list$responses[, 1]), 10))
#       , do.call(cbind, rep(list(sim_list$responses[, 2]), 10))
#     )
#     sim_list$responses3 <- sim_list$responses2 - sim_list$responses[, -c(1,2)] * (2 * sim_list$responses2 - 1)
# 
#     set.seed(4)
#     out <- dependentLCM_fit(df = sim_list$responses3, nclass=1, class2domain = c(0), nitr= nitr, class_init_method = "kmodes"
#                             , domain_maxitems=2, domain_proposal_swap=0.3
#     )
#     out_summary <- dlcm.summary(out)
#     sink.reset()
#     }
# )


system.time(
  {
    nitr <- 6000
    close(file(outfile_path, open="w")); sink(outfile_path)
    rm(out, out_summary)
    set.seed(4)
    sim_list <- generate2(n=1000
                          , pis=c(0.5, 0.5)
                          , thetas_list = list(
                            c( rep(list(c(0.8,0.2)), 20)
                               , rep(list(c(0.6,0.3,0.05,0.05)), 2)
                               , rep(list(c(rep(0.3,3), rep(0.1/5,5))), 1)
                            )
                            , c( rep(list(c(0.5,0.5)), 20)
                                 , rep(list(c(0.3,0.6,0.05,0.05)), 2)
                                 , rep(list(c(rep(0.4,2), rep(0.2/6,6))), 1)
                            )
                          ))
    pattern_lookups <- list(
      items2 = expand.grid(c(0,1), c(0,1))[c(4,1,2,3), ]
      , items3 = expand.grid(c(0,1), c(0,1), c(0,1))[c(8,1,2:7), ]
    )
    sim_list$responses2 <- cbind(
      sim_list$responses[, 1:20]
      , pattern_lookups[["items2"]][sim_list$responses[, 21], ]
      , pattern_lookups[["items2"]][sim_list$responses[, 22], ]
      , pattern_lookups[["items3"]][sim_list$responses[, 23], ]
    )
    
    set.seed(4)
    out <- dependentLCM_fit(df = sim_list$responses2, nclass=2, class2domain = c(0, 0), nitr= nitr
                            , class_init_method = "kmodes"
                            , domain_proposal_swap=0.3
    )
    out_summary <- dlcm.summary(out)
    sink.reset()
  }
)

