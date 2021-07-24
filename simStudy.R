sink.reset <- function(){
  for(i in seq_len(sink.number())){
    sink(NULL)
  }
}
sink.reset()
outfile_path <- "C:/Users/drakethrice/Desktop/analysis-output.txt" # ; close(file(outfile_path, open="w")); sink(outfile_path)


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

system.time(
  {
    nitr <- 6000
    close(file(outfile_path, open="w")); sink(outfile_path)
    rm(out)
    set.seed(4)
    sim_list <- generate(
      n = 1000
      , pis = c(1)
      , thetas = cbind(
        c(rep(0.2, 1), rep(0.8, 1), rep(0.9,20))
      ))
    sim_list$responses2 <- cbind(
      do.call(cbind, rep(list(sim_list$responses[, 1]), 10))
      , do.call(cbind, rep(list(sim_list$responses[, 2]), 10))
    )
    sim_list$responses3 <- sim_list$responses2 - sim_list$responses[, -c(1,2)] * (2 * sim_list$responses2 - 1)
    
    set.seed(4)
    out <- dependentLCM_fit(df = sim_list$responses3, nclass=1, class2domain = c(0), nitr= nitr, class_init_method = "kmodes", domain_proposal_ratio=0.5
                            , domain_maxitems=2
    )
    out_summary <- dlcm.summary(out)
    sink.reset()
  }
)

