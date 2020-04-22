#' Run one simulation
#'
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @importFrom purrr map
#'
#' @export

run_one_sim <- function(uid, seed){

    # Get simulation parameters
    model_pars <- get_pars(uid = uid)

    # # Return model parameters
    #             iter = scenario[["iter"]],
    #             burnin = scenario[["burnin"]],
    #             repetitions = scenario[["repetitions"]],
    #             scenario_uid = as.character(scenario[["scenario_uid"]]),
    #             uid = as.character(scenario[["uid"]]),
    #             save_all = scenario[["save_all"]],

    # Simulate data
    sim_data <- sim_mHMM(
        # Number of observations for each person
        n_t = model_pars[["n_t"]],
        # Number of persons
        n = model_pars[["sample_size"]],
        # Type of emission distributions
        data_distr = "categorical",
        # Number of states
        m = model_pars[["m"]],
        # Number of emission distributions
        n_dep = model_pars[["n_dep"]],
        # Number of categories per dependent variable
        q_emiss = model_pars[["q_emiss"]],
        # Transition probabilities
        gamma = model_pars[["gamma_sim"]],
        # Emission distribution means + var
        emiss_distr = model_pars[["emiss_sim"]],
        # Between-subject variance for TPM
        var_gamma = model_pars[["gamma_var"]],
        # Between-subject variance for emission distributions
        var_emiss = model_pars[["emiss_var"]],
        # Additional arguments
        return_ind_par = TRUE
    )

    # Fit mHMMbayes model
    model_output <- fit_mHMM(
        # Number of states
        m = model_pars[["m"]],
        # Number of emission distributions
        n_dep = model_pars[["n_dep"]],
        # Number of categories per dependent variable
        q_emiss = model_pars[["q_emiss"]],
        # Transition probabilities
        gamma = model_pars[["gamma_sim"]],
        # Emission distribution means + var
        emiss = model_pars[["emiss_sim"]],
        # Number of iterations
        iter = model_pars[["iter"]],
        # Burn-in iterations
        burnin = model_pars[["burnin"]],
        # Starting values for the transition probabilities
        start_gamma = NULL,
        # Starting values for the emission distributions
        start_emiss = NULL,
        # Simulated data
        data_sim = sim_data
    )

    # Get MAP estimates

    # Get evaluation metrics

    # Save results: add the actual outcomes
    if (model_pars[["save_all"]]) {
        complete_data <- list("sim_data" = sim_data,
             "output" = model_output)
        saveRDS(object = complete_data, file = paste0("outputs/complete_results/",model_pars[["uid"]],".rds"))
    } else {
        saveRDS(object = out, file = paste0("outputs/",as.character(uidTarget),".rds"))
    }

}
