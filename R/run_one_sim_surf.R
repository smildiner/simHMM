#' Run one simulation in SurfSara Lisa
#'
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @importFrom purrr map
#'
#' @export

run_one_sim_surf <- function(pars, light = FALSE, save_subj_data = TRUE){

    # Legend
    #   pars[1] = sample_size
    #   pars[2] = n_t
    #   pars[3] = n_dep
    #   pars[4] = noisiness
    #   pars[5] = overlapping
    #   pars[6] = iter
    #   pars[7] = burnin
    #   pars[8] = repetitions
    #   pars[9] = scenario_uid
    #   pars[10] = uid
    #   pars[11:17] = .Random.seed

    exe_time <- system.time({

        # Put in the right format
        pars[9] <- deparse(substitute(pars[9]))
        pars[10] <- deparse(substitute(pars[10]))

        # Set L'Ecuyer random seed
        RNGkind("L'Ecuyer-CMRG")
        set.seed(42)
        .Random.seed <- matrix(pars[11:17], nrow = 1)

        # Store the current state of the stream of RNG
        seed_state <- list(state = .Random.seed,
                           kind = RNGkind())

        # Get simulation parameters
        model_pars <- get_pars_surf(pars)

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
            data_sim = sim_data,
            # Fit mHMM with lower memory use
            light = light,
            # Save subject level results
            save_subj_data = save_subj_data
        )

        # Add between subject variace to the output
        if(light == FALSE) {
            model_output <- c(model_output, get_var_bar(model_output))
        }

        # Get MAP estimates
        map_out <- MAP(model_output)

        # Get credibility intervals
        cci_out <- get_cci(model_output)
    })

    # Add scenario and iteration info
    out <- list(seed = seed_state,
                scenario_uid = model_pars[["scenario_uid"]],
                uid = model_pars[["uid"]],
                time = exe_time[[3]],
                truth = sim_data,
                map = map_out,
                cci = cci_out)

    # Save results: add the actual outcomes
    saveRDS(object = out, file = paste0(model_pars[["uid"]],".rds"))

}