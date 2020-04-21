#' Initialize simulation log
#'
#' @importFrom digest digest
#'
#' @export

initialize_sim <- function(sim_sample_size = 30,
                           sim_n_t = c(400, 800, 1600),
                           sim_m = 3,
                           sim_n_dep = c(1, 2, 4),
                           sim_q_emiss = 5,
                           sim_gamma_var = 1,
                           sim_emiss_var = 1,
                           sim_noisiness = c(0.03, 0.09, 0.15),
                           sim_overlapping = c("low", "moderate", "high"),

                           # Simulation parameters:
                           sim_iter = 5000,
                           sim_burnin = 1750,
                           sim_repetitions = 250) {

    # Generate log with scenarios:
    scenarios_log <- expand.grid(sim_sample_size,
                                 sim_n_t,
                                 sim_m,
                                 sim_n_dep,
                                 sim_q_emiss,
                                 sim_gamma_var,
                                 sim_emiss_var,
                                 sim_noisiness,
                                 sim_overlapping,
                                 sim_iter,
                                 sim_burnin,
                                 sim_repetitions)

    colnames(scenarios_log) <- c("sample_size",
                             "n_t",
                             "m",
                             "n_dep",
                             "q_emiss",
                             "gamma_var",
                             "emiss_var",
                             "noisiness",
                             "overlapping",
                             "iter",
                             "burnin",
                             "repetitions")

    scenarios_log <- cbind(scenarios_log,
                           scenario_uid = apply(scenarios_log, 1, digest))
    scenarios_log <- cbind(scenarios_log,
                           uid = apply(scenarios, 1, digest))

    # Save rds file
    saveRDS(scenarios_log, "scenarios_log.rds")

}
