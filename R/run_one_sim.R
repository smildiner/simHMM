#' Run one simulation
#'
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @importFrom purrr map
#'
#' @export

run_one_sim <- function(uid, seed){

    # Get simulation parameters
    sim_pars <- get_pars(uid = uid)

    # # Return model parameters
    # return(list(sample_size = scenario[["sample_size"]],
    #             n_t = scenario[["n_t"]],
    #             m = scenario[["m"]],
    #             n_dep = scenario[["n_dep"]],
    #             q_emiss = rep(scenario[["q_emiss"]], scenario[["n_dep"]]),
    #             gamma_var = scenario[["gamma_var"]],
    #             emiss_var = rep(scenario[["emiss_var"]], scenario[["n_dep"]]),
    #             noisiness = scenario[["noisiness"]],
    #             overlapping = as.character(scenario[["overlapping"]]),
    #             iter = scenario[["iter"]],
    #             burnin = scenario[["burnin"]],
    #             repetitions = scenario[["repetitions"]],
    #             scenario_uid = as.character(scenario[["scenario_uid"]]),
    #             uid = as.character(scenario[["uid"]]),
    #             save_all = scenario[["save_all"]],
    #             gamma_sim = gamma_sim,
    #             emiss_sim = emiss_sim))

    # Simulate data
    sim_data <- sim_mHMM(
        # Number of observations for each person
        n_t = sim_pars[["n_t"]],
        # Number of persons
        n = sim_pars[["sample_size"]],
        # Number of categories per dependent variable
        q_emiss = sim_pars[["q_emiss"]],
        # Type of emission distributions
        data_distr = "categorical",
        # Number of states
        m = sim_pars[["m"]],
        # Number of emission distributions
        n_dep = sim_pars[["n_dep"]],
        # Transition probabilities
        gamma = sim_pars[["gamma_sim"]],
        # Emission distribution means + var
        emiss_distr = emiss[1:n_dep],
        # Between-subject variance for TPM
        var_gamma = gammaVar,
        # Between-subject variance for emission distributions
        var_emiss = rep(emissVar, n_dep),
        # Additional arguments
        return_ind_par = TRUE
    )

    # Fit mHMMbayes model

    # Get MAP estimates

    # Get evaluation metrics

    # Save results

}
