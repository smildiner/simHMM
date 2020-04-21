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
                           sim_repetitions = 1:250) {

    # Generate log with scenarios:
    scenarios_log <- expand.grid("sample_size" = sim_sample_size,
                                 "n_t" = sim_n_t,
                                 "m" = sim_m,
                                 "n_dep" = sim_n_dep,
                                 "q_emiss" = sim_q_emiss,
                                 "gamma_var" = sim_gamma_var,
                                 "emiss_var" = sim_emiss_var,
                                 "noisiness" = sim_noisiness,
                                 "overlapping" = sim_overlapping,
                                 "iter" = sim_iter,
                                 "burnin" = sim_burnin,
                                 "repetitions" = sim_repetitions)

    # Add scenario id
    scenarios_log <- cbind(scenarios_log,
                           scenario_uid = apply(scenarios_log, 1, digest))

    # Add simulation id
    scenarios_log <- cbind(scenarios_log,
                           uid = apply(scenarios_log, 1, digest))

    # Save rds file
    if ("inputs" %in% dir()) {
        saveRDS(scenarios_log, "inputs/scenarios_log.rds")
    } else {
        dir.create("inputs")
        saveRDS(scenarios_log, "inputs/scenarios_log.rds")
    }

}
