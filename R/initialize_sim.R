#' Initialize simulation log
#'
#' @importFrom digest digest
#' @importFrom dplyr group_by sample_frac
#' @importFrom magrittr %>%
#'
#' @export

init_log <- function(sim_sample_size = 30,
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
                     sim_burnin = 1500,
                     sim_repetitions = 1:250,

                     # Save complete results:
                     save_frac_res = 0.02,
                     seed = 123) {

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
                                 "repetitions" = sim_repetitions,
                                 stringsAsFactors = FALSE)

    # Add scenario id
    scenarios_log <- cbind.data.frame(scenarios_log,
                           scenario_uid = apply(scenarios_log[,-which(names(scenarios_log) == "repetitions")], 1, digest),
                           stringsAsFactors = FALSE)

    # Add simulation id
    scenarios_log <- cbind.data.frame(scenarios_log,
                           uid = apply(scenarios_log, 1, digest),
                           stringsAsFactors = FALSE)

    # Specify which simulation will be saved
    set.seed(seed)
    save_res_idx <- scenarios_log %>% group_by(scenario_uid) %>% sample_frac(size = save_frac_res)
    scenarios_log[["save_all"]] <- scenarios_log[["uid"]] %in% save_res_idx[["uid"]]

    # Save rds file
    if ("inputs" %in% dir()) {
        saveRDS(scenarios_log, "inputs/scenarios_log.rds")
    } else {
        dir.create("inputs")
        dir.create("outputs")
        dir.create("outputs/complete_results")
        dir.create("outputs/results")
        saveRDS(scenarios_log, "inputs/scenarios_log.rds")
    }

}
