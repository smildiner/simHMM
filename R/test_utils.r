#' Run one simulation in SurfSara Lisa
#'
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @importFrom purrr map
#'
#' @export

test_run_one_sim_surf <- function(pars, light = FALSE, save_subj_data = TRUE, baseline = FALSE, convergence = FALSE, save_path = FALSE, progress_line = FALSE){

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
    #   pars[11:17] = .Random.seed for simulation
    #   pars[18] = new/old
    #   pars[19] = good/bad
    ##   pars[18:24] = .Random.seed for convergence run

    exe_time <- system.time({

        # Put in the right format
        pars[9] <- pars[9]
        pars[10] <- pars[10]

        # Set L'Ecuyer random seed
        RNGkind("L'Ecuyer-CMRG")
        set.seed(42)
        .Random.seed <<- as.integer(matrix(as.numeric(pars[11:17]), nrow = 1))

        # Store the current state of the stream of RNG
        seed_state <- list(state = .Random.seed,
                           kind = RNGkind())

        # Get simulation parameters
        model_pars <- test_get_pars_surf(pars, baseline, as.character(pars[18]))

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

        # if(convergence == TRUE){
        #     # Set L'Ecuyer random seed
        #     RNGkind("L'Ecuyer-CMRG")
        #     set.seed(42)
        #     .Random.seed <<- as.integer(matrix(as.numeric(pars[18:24]), nrow = 1))
        #
        #     # Store the current state of the stream of RNG
        #     seed_convergence <- list(state = .Random.seed,
        #                              kind = RNGkind())
        # }

        # Fit mHMMbayes model
        model_output <- test_fit_mHMM(
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
            # Save local decoding
            save_path = save_path,
            # Save subject level results
            save_subj_data = save_subj_data,
            progress_line = progress_line,
            good_bad = as.character(pars[19])
        )

        # Add empirical between subject variance to the output
        if(light == FALSE) {
            model_output <- c(model_output, get_var_bar(model_output))
        }

        # Save local decoding?
        if(save_path == TRUE){
            local_decode <- lapply(model_output[["sample_path"]], function(e) t(apply(e[,(model_pars[["burnin"]]+1):model_pars[["iter"]]],1,function(r) {
                r <- factor(r, levels = 1:model_pars[["m"]], labels = paste0("S",1:model_pars[["m"]]))
                table(r)/length(r)
            })))
        }

        # Get MAP estimates
        map_out <- MAP(model_output)

        # Get credibility intervals
        cci_out <- get_cci(model_output)
    })

    # Add scenario and iteration info

    out <- list(seed = seed_state,
                # seed = seed_state,
                scenario_uid = model_pars[["scenario_uid"]],
                uid = model_pars[["uid"]],
                time = exe_time[[3]],
                truth = sim_data,
                map = map_out,
                cci = cci_out,
                output = model_output)

    if(save_path == TRUE){
        out <- c(out, ldecoding = list(local_decode))
    }

    return(out)

}


#' Fit mHMM using mHMMfast
#'
#' @export

test_fit_mHMM <- function(m,
                          n_dep,
                          q_emiss,
                          gamma,
                          emiss,
                          iter,
                          burnin,
                          start_gamma = NULL,
                          start_emiss = NULL,
                          data_sim,
                          light = FALSE,
                          save_path = FALSE,
                          save_subj_data = TRUE,
                          progress_line = FALSE,
                          good_bad = "good") {

    # Set starting values
    # gamma_start
    if (is.null(start_gamma)) {

        if(good_bad == "good"){
            start_gamma <- diag(c(runif(1, 0.85, 0.95),runif(1, 0.65, 0.75),runif(1, 0.45, 0.55)))
            # start_gamma <- diag(c(1,0,0))
            start_gamma[lower.tri(start_gamma) | upper.tri(start_gamma)] <- (1 - diag(start_gamma)) / (m - 1)
        } else {
            start_gamma <- diag(runif(1, 0.7, 0.9), m)
            start_gamma[lower.tri(start_gamma) | upper.tri(start_gamma)] <- (1 - diag(start_gamma)) / (m - 1)
        }

    }

    # emiss_start
    if (is.null(start_emiss)) {

        start_emiss <- emiss
        # Set seed to simulate datasets
        for (q in 1:n_dep) {

            start_emiss[[q]] <- int_to_prob(prob_to_int(emiss[[q]]) + matrix(runif(m*(q_emiss[q]-1), -1, 1), byrow = T, nrow = m, ncol = (q_emiss[q] - 1)))

        }

    }

    # Run the model
    ti <- Sys.time()
    if(light == FALSE){
        out <- simHMM::mHMMfast(s_data = data_sim$obs,
                                gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                                # xx = xx_vec,
                                start_val = c(list(start_gamma), start_emiss),
                                mcmc = list(J = iter, burn_in = burnin),
                                return_path = save_path,
                                show_progress = progress_line)
    } else {
        out <- simHMM::mHMMlight(s_data = data_sim$obs,
                                 gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                                 # xx = xx_vec,
                                 start_val = c(list(start_gamma), start_emiss),
                                 mcmc = list(J = iter, burn_in = burnin),
                                 return_path = save_path,
                                 show_progress = progress_line,
                                 save_subj_data = save_subj_data)
    }
    out[["time"]] <- Sys.time() - ti

    # Return model output
    return(out)

}

#' Load simulation parameters
#'
#' @export

test_get_pars_surf <- function(pars, baseline = FALSE, new_old = "new") {

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
    #   pars[18] = new or old
    #   pars[19] = good or bad

    # Set marker for noisiness
    if(as.numeric(pars[4]) == 0.03){
        noise <- "low"
    } else if(as.numeric(pars[4]) == 0.09){
        noise <- "moderate"
    } else if(as.numeric(pars[4]) == 0.15){
        noise <- "high"
    }


    # Specify the correct gamma, emiss and eps_str
    if(baseline == TRUE) {
        gamma_sim <- matrix(c(0.9, 0.05, 0.05,
                              0.1, 0.7, 0.2,
                              0.2, 0.3, 0.5), ncol = 3, byrow = TRUE)

        gamma_var <- matrix(c(1.12, 1.12,
                              0.16, 0.16,
                              0.13, 0.13), nrow = 3, byrow = TRUE)

    } else {
        gamma_sim <- matrix(c(0.9, 0.05, 0.05,
                              0.1, 0.7, 0.2,
                              0.2, 0.3, 0.5), ncol = 3, byrow = TRUE)

        gamma_var <- matrix(c(1.12, 1.12,
                              0.16, 0.16,
                              0.13, 0.13), nrow = 3, byrow = TRUE)

        # gamma_sim <- matrix(c(0.96, 0.02, 0.02,
        #                       0.03, 0.94, 0.03,
        #                       0.04, 0.04, 0.92), ncol = 3, byrow = TRUE)
    }

    if(new_old == "new") {
        gamma_sim <- matrix(c(0.9, 0.05, 0.05,
                              0.1, 0.7, 0.2,
                              0.2, 0.3, 0.5), ncol = 3, byrow = TRUE)

        gamma_var <- matrix(c(1.12, 1.12,
                              0.16, 0.16,
                              0.13, 0.13), nrow = 3, byrow = TRUE)

        emiss_sim <- list("low" = list(matrix(c(0.95, 0.02, 0.01, 0.01, 0.01,
                                                0.02, 0.48, 0.48, 0.01, 0.01,
                                                0.02, 0.01, 0.01, 0.48, 0.48), nrow = 3, ncol = 5, byrow = T),

                                       matrix(c(0.02, 0.01, 0.01, 0.48, 0.48,
                                                0.95, 0.02, 0.01, 0.01, 0.01,
                                                0.02, 0.48, 0.48, 0.01, 0.01), nrow = 3, ncol = 5, byrow = T),

                                       matrix(c(0.02, 0.48, 0.48, 0.01, 0.01,
                                                0.02, 0.01, 0.01, 0.48, 0.48,
                                                0.95, 0.02, 0.01, 0.01, 0.01), nrow = 3, ncol = 5, byrow = T),

                                       matrix(c(0.02, 0.01, 0.48, 0.01, 0.48,
                                                0.95, 0.02, 0.01, 0.01, 0.01,
                                                0.02, 0.48, 0.01, 0.48, 0.01), nrow = 3, ncol = 5, byrow = T)),

                          "moderate" = list(matrix(c(0.95, 0.02, 0.01, 0.01, 0.01,
                                                     0.02, 0.48, 0.48, 0.01, 0.01,
                                                     0.02, 0.01, 0.48, 0.48, 0.01), nrow = 3, ncol = 5, byrow = T),

                                            matrix(c(0.02, 0.01, 0.01, 0.48, 0.48,
                                                     0.95, 0.02, 0.01, 0.01, 0.01,
                                                     0.02, 0.01, 0.48, 0.48, 0.01), nrow = 3, ncol = 5, byrow = T),

                                            matrix(c(0.02, 0.48, 0.01, 0.01, 0.48,
                                                     0.02, 0.01, 0.01, 0.48, 0.48,
                                                     0.95, 0.02, 0.01, 0.01, 0.01), nrow = 3, ncol = 5, byrow = T),

                                            matrix(c(0.02, 0.01, 0.48, 0.01, 0.48,
                                                     0.95, 0.02, 0.01, 0.01, 0.01,
                                                     0.02, 0.48, 0.01, 0.01, 0.48), nrow = 3, ncol = 5, byrow = T)),

                          "high" = list(matrix(c(0.02, 0.95, 0.01, 0.01, 0.01,
                                                 0.02, 0.48, 0.48, 0.01, 0.01,
                                                 0.02, 0.01, 0.48, 0.48, 0.01), nrow = 3, ncol = 5, byrow = T),

                                        matrix(c(0.02, 0.01, 0.01, 0.48, 0.48,
                                                 0.02, 0.01, 0.01, 0.01, 0.95,
                                                 0.02, 0.01, 0.48, 0.48, 0.01), nrow = 3, ncol = 5, byrow = T),

                                        matrix(c(0.02, 0.48, 0.01, 0.01, 0.48,
                                                 0.02, 0.01, 0.01, 0.48, 0.48,
                                                 0.02, 0.01, 0.01, 0.95, 0.01), nrow = 3, ncol = 5, byrow = T),

                                        matrix(c(0.02, 0.01, 0.48, 0.01, 0.48,
                                                 0.02, 0.01, 0.95, 0.01, 0.01,
                                                 0.02, 0.48, 0.01, 0.01, 0.48), nrow = 3, ncol = 5, byrow = T)))

        eps_str <- list("low" = list(matrix(c(-4*1, 1, 1, 1, 1,
                                              1/1.5, -3/2*1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5,
                                              1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5, -3/2*1/1.5), nrow = 3, ncol = 5, byrow = T),

                                     matrix(c(1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5, -3/2*1/1.5,
                                              -4*1, 1, 1, 1, 1,
                                              1/1.5, -3/2*1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5), nrow = 3, ncol = 5, byrow = T),

                                     matrix(c(1/1.5, -3/2*1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5,
                                              1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5, -3/2*1/1.5,
                                              -4*1, 1, 1, 1, 1), nrow = 3, ncol = 5, byrow = T),

                                     matrix(c(1/1.5, 1/1.5, -3/2*1/1.5, 1/1.5, -3/2*1/1.5,
                                              -4*1, 1, 1, 1, 1,
                                              1/1.5, -3/2*1/1.5, 1/1.5, -3/2*1/1.5, 1/1.5), nrow = 3, ncol = 5, byrow = T)),

                        "moderate" = list(matrix(c(-4*1, 1, 1, 1, 1,
                                                   1/1.5, -3/2*1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5,
                                                   1/1.5, 1/1.5, -3/2*1/1.5, -3/2*1/1.5, 1/1.5), nrow = 3, ncol = 5, byrow = T),

                                          matrix(c(1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5, -3/2*1/1.5,
                                                   -4*1, 1, 1, 1, 1,
                                                   1/1.5, 1/1.5, -3/2*1/1.5, -3/2*1/1.5, 1/1.5), nrow = 3, ncol = 5, byrow = T),

                                          matrix(c(1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5,
                                                   1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5, -3/2*1/1.5,
                                                   -4*1, 1, 1, 1, 1), nrow = 3, ncol = 5, byrow = T),

                                          matrix(c(1/1.5, 1/1.5, -3/2*1/1.5, 1/1.5, -3/2*1/1.5,
                                                   -4*1, 1, 1, 1, 1,
                                                   1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5), nrow = 3, ncol = 5, byrow = T)),

                        "high" = list(matrix(c(1, -4*1, 1, 1, 1,
                                               1/1.5, -3/2*1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5,
                                               1/1.5, 1/1.5, -3/2*1/1.5, -3/2*1/1.5, 1/1.5), nrow = 3, ncol = 5, byrow = T),

                                      matrix(c(1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5, -3/2*1/1.5,
                                               1, 1, 1, 1, -4*1,
                                               1/1.5, 1/1.5, -3/2*1/1.5, -3/2*1/1.5, 1/1.5), nrow = 3, ncol = 5, byrow = T),

                                      matrix(c(1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5,
                                               1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5, -3/2*1/1.5,
                                               1, 1, 1, -4*1, 1), nrow = 3, ncol = 5, byrow = T),

                                      matrix(c(1/1.5, 1/1.5, -3/2*1/1.5, 1/1.5, -3/2*1/1.5,
                                               1, 1, -4*1, 1, 1,
                                               1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5), nrow = 3, ncol = 5, byrow = T)))

        emiss_var <- list("low" = list("low" = list(matrix(c(1.01, 1.01, 1.01, 1.01,
                                                             0.11, 0.11, 0.11, 0.11,
                                                             0.11, 0.11, 0.11, 0.11), nrow = 3, ncol = 4, byrow = T),

                                                    matrix(c(0.11, 0.11, 0.11, 0.11,
                                                             1.01, 1.01, 1.01, 1.01,
                                                             0.11, 0.11, 0.11, 0.11), nrow = 3, ncol = 4, byrow = T),

                                                    matrix(c(0.11, 0.11, 0.11, 0.11,
                                                             0.11, 0.11, 0.11, 0.11,
                                                             1.01, 1.01, 1.01, 1.01), nrow = 3, ncol = 4, byrow = T),

                                                    matrix(c(0.11, 0.11, 0.11, 0.11,
                                                             1.01, 1.01, 1.01, 1.01,
                                                             0.11, 0.11, 0.11, 0.11), nrow = 3, ncol = 4, byrow = T)),

                                       "moderate" = list(matrix(c(0.60, 0.60, 0.60, 0.60,
                                                                  0.14, 0.14, 0.14, 0.14,
                                                                  0.14, 0.14, 0.14, 0.14), nrow = 3, ncol = 4, byrow = T),

                                                         matrix(c(0.14, 0.14, 0.14, 0.14,
                                                                  0.60, 0.60, 0.60, 0.60,
                                                                  0.14, 0.14, 0.14, 0.14), nrow = 3, ncol = 4, byrow = T),

                                                         matrix(c(0.14, 0.14, 0.14, 0.14,
                                                                  0.14, 0.14, 0.14, 0.14,
                                                                  0.60, 0.60, 0.60, 0.60), nrow = 3, ncol = 4, byrow = T),

                                                         matrix(c(0.14, 0.14, 0.14, 0.14,
                                                                  0.60, 0.60, 0.60, 0.60,
                                                                  0.14, 0.14, 0.14, 0.14), nrow = 3, ncol = 4, byrow = T)),

                                       "high" = list(matrix(c(0.82, 0.82, 0.82, 0.82,
                                                              0.17, 0.17, 0.17, 0.17,
                                                              0.17, 0.17, 0.17, 0.17), nrow = 3, ncol = 4, byrow = T),

                                                     matrix(c(0.17, 0.17, 0.17, 0.17,
                                                              0.82, 0.82, 0.82, 0.82,
                                                              0.17, 0.17, 0.17, 0.17), nrow = 3, ncol = 4, byrow = T),

                                                     matrix(c(0.17, 0.17, 0.17, 0.17,
                                                              0.17, 0.17, 0.17, 0.17,
                                                              0.82, 0.82, 0.82, 0.82), nrow = 3, ncol = 4, byrow = T),

                                                     matrix(c(0.17, 0.17, 0.17, 0.17,
                                                              0.82, 0.82, 0.82, 0.82,
                                                              0.17, 0.17, 0.17, 0.17), nrow = 3, ncol = 4, byrow = T))),


                          "moderate" = list("low" = list(matrix(c(1.01, 1.01, 1.01, 1.01,
                                                                  0.11, 0.11, 0.11, 0.11,
                                                                  0.11, 0.11, 0.11, 0.11), nrow = 3, ncol = 4, byrow = T),

                                                         matrix(c(0.11, 0.11, 0.11, 0.11,
                                                                  1.01, 1.01, 1.01, 1.01,
                                                                  0.11, 0.11, 0.11, 0.11), nrow = 3, ncol = 4, byrow = T),

                                                         matrix(c(0.11, 0.11, 0.11, 0.11,
                                                                  0.11, 0.11, 0.11, 0.11,
                                                                  1.01, 1.01, 1.01, 1.01), nrow = 3, ncol = 4, byrow = T),

                                                         matrix(c(0.11, 0.11, 0.11, 0.11,
                                                                  1.01, 1.01, 1.01, 1.01,
                                                                  0.11, 0.11, 0.11, 0.11), nrow = 3, ncol = 4, byrow = T)),

                                            "moderate" = list(matrix(c(0.60, 0.60, 0.60, 0.60,
                                                                       0.14, 0.14, 0.14, 0.14,
                                                                       0.14, 0.14, 0.14, 0.14), nrow = 3, ncol = 4, byrow = T),

                                                              matrix(c(0.14, 0.14, 0.14, 0.14,
                                                                       0.60, 0.60, 0.60, 0.60,
                                                                       0.14, 0.14, 0.14, 0.14), nrow = 3, ncol = 4, byrow = T),

                                                              matrix(c(0.14, 0.14, 0.14, 0.14,
                                                                       0.14, 0.14, 0.14, 0.14,
                                                                       0.60, 0.60, 0.60, 0.60), nrow = 3, ncol = 4, byrow = T),

                                                              matrix(c(0.14, 0.14, 0.14, 0.14,
                                                                       0.60, 0.60, 0.60, 0.60,
                                                                       0.14, 0.14, 0.14, 0.14), nrow = 3, ncol = 4, byrow = T)),

                                            "high" = list(matrix(c(0.82, 0.82, 0.82, 0.82,
                                                                   0.17, 0.17, 0.17, 0.17,
                                                                   0.17, 0.17, 0.17, 0.17), nrow = 3, ncol = 4, byrow = T),

                                                          matrix(c(0.17, 0.17, 0.17, 0.17,
                                                                   0.82, 0.82, 0.82, 0.82,
                                                                   0.17, 0.17, 0.17, 0.17), nrow = 3, ncol = 4, byrow = T),

                                                          matrix(c(0.17, 0.17, 0.17, 0.17,
                                                                   0.17, 0.17, 0.17, 0.17,
                                                                   0.82, 0.82, 0.82, 0.82), nrow = 3, ncol = 4, byrow = T),

                                                          matrix(c(0.17, 0.17, 0.17, 0.17,
                                                                   0.82, 0.82, 0.82, 0.82,
                                                                   0.17, 0.17, 0.17, 0.17), nrow = 3, ncol = 4, byrow = T))),


                          "high" = list("low" = list(matrix(c(0.35, 0.35, 0.35, 0.35,
                                                              0.11, 0.11, 0.11, 0.11,
                                                              0.11, 0.11, 0.11, 0.11), nrow = 3, ncol = 4, byrow = T),

                                                     matrix(c(0.11, 0.11, 0.11, 0.11,
                                                              0.35, 0.35, 0.35, 0.35,
                                                              0.11, 0.11, 0.11, 0.11), nrow = 3, ncol = 4, byrow = T),

                                                     matrix(c(0.11, 0.11, 0.11, 0.11,
                                                              0.11, 0.11, 0.11, 0.11,
                                                              0.35, 0.35, 0.35, 0.35), nrow = 3, ncol = 4, byrow = T),

                                                     matrix(c(0.11, 0.11, 0.11, 0.11,
                                                              0.35, 0.35, 0.35, 0.35,
                                                              0.11, 0.11, 0.11, 0.11), nrow = 3, ncol = 4, byrow = T)),

                                        "moderate" = list(matrix(c(0.15, 0.15, 0.15, 0.15,
                                                                   0.14, 0.14, 0.14, 0.14,
                                                                   0.14, 0.14, 0.14, 0.14), nrow = 3, ncol = 4, byrow = T),

                                                          matrix(c(0.14, 0.14, 0.14, 0.14,
                                                                   0.15, 0.15, 0.15, 0.15,
                                                                   0.14, 0.14, 0.14, 0.14), nrow = 3, ncol = 4, byrow = T),

                                                          matrix(c(0.14, 0.14, 0.14, 0.14,
                                                                   0.14, 0.14, 0.14, 0.14,
                                                                   0.15, 0.15, 0.15, 0.15), nrow = 3, ncol = 4, byrow = T),

                                                          matrix(c(0.14, 0.14, 0.14, 0.14,
                                                                   0.15, 0.15, 0.15, 0.15,
                                                                   0.14, 0.14, 0.14, 0.14), nrow = 3, ncol = 4, byrow = T)),

                                        "high" = list(matrix(c(0.18, 0.18, 0.18, 0.18,
                                                               0.17, 0.17, 0.17, 0.17,
                                                               0.17, 0.17, 0.17, 0.17), nrow = 3, ncol = 4, byrow = T),

                                                      matrix(c(0.17, 0.17, 0.17, 0.17,
                                                               0.18, 0.18, 0.18, 0.18,
                                                               0.17, 0.17, 0.17, 0.17), nrow = 3, ncol = 4, byrow = T),

                                                      matrix(c(0.17, 0.17, 0.17, 0.17,
                                                               0.17, 0.17, 0.17, 0.17,
                                                               0.18, 0.18, 0.18, 0.18), nrow = 3, ncol = 4, byrow = T),

                                                      matrix(c(0.17, 0.17, 0.17, 0.17,
                                                               0.18, 0.18, 0.18, 0.18,
                                                               0.17, 0.17, 0.17, 0.17), nrow = 3, ncol = 4, byrow = T))))

        emiss_var <- emiss_var[[as.character(pars[5])]][[noise]]

    } else {
        gamma_sim <- matrix(c(0.96, 0.02, 0.02,
                              0.03, 0.94, 0.03,
                              0.04, 0.04, 0.92), ncol = 3, byrow = TRUE)
        gamma_var <- 1

        emiss_sim <- list("low" = list(matrix(c(0.96, 0.01, 0.01, 0.01, 0.01,
                                                0.02, 0.47, 0.47, 0.02, 0.02,
                                                0.02, 0.02, 0.02, 0.47, 0.47), nrow = 3, ncol = 5, byrow = T),

                                       matrix(c(0.02, 0.47, 0.47, 0.02, 0.02,
                                                0.96, 0.01, 0.01, 0.01, 0.01,
                                                0.02, 0.02, 0.02, 0.47, 0.47), nrow = 3, ncol = 5, byrow = T),

                                       matrix(c(0.01, 0.96, 0.01, 0.01, 0.01,
                                                0.47, 0.02, 0.47, 0.02, 0.02,
                                                0.02, 0.02, 0.02, 0.47, 0.47), nrow = 3, ncol = 5, byrow = T),

                                       matrix(c(0.47, 0.02, 0.47, 0.02, 0.02,
                                                0.02, 0.02, 0.02, 0.47, 0.47,
                                                0.01, 0.96, 0.01, 0.01, 0.01), nrow = 3, ncol = 5, byrow = T),

                                       matrix(c(0.01, 0.01, 0.96, 0.01, 0.01,
                                                0.47, 0.47, 0.02, 0.02, 0.02,
                                                0.02, 0.02, 0.02, 0.47, 0.47), nrow = 3, ncol = 5, byrow = T),

                                       matrix(c(0.01, 0.01, 0.01, 0.96, 0.01,
                                                0.02, 0.47, 0.47, 0.02, 0.02,
                                                0.47, 0.02, 0.02, 0.02, 0.47), nrow = 3, ncol = 5, byrow = T),

                                       matrix(c(0.02, 0.47, 0.47, 0.02, 0.02,
                                                0.47, 0.02, 0.02, 0.02, 0.47,
                                                0.01, 0.01, 0.01, 0.96, 0.01), nrow = 3, ncol = 5, byrow = T),

                                       matrix(c(0.01, 0.01, 0.01, 0.01, 0.96,
                                                0.02, 0.47, 0.47, 0.02, 0.02,
                                                0.47, 0.02, 0.02, 0.47, 0.02), nrow = 3, ncol = 5, byrow = T)),

                          "moderate" = list(matrix(c(0.96, 0.01, 0.01, 0.01, 0.01,
                                                     0.02, 0.47, 0.47, 0.02, 0.02,
                                                     0.02, 0.02, 0.47, 0.47, 0.02), nrow = 3, ncol = 5, byrow = T),

                                            matrix(c(0.02, 0.47, 0.47, 0.02, 0.02,
                                                     0.96, 0.01, 0.01, 0.01, 0.01,
                                                     0.02, 0.47, 0.02, 0.02, 0.47), nrow = 3, ncol = 5, byrow = T),

                                            matrix(c(0.01, 0.96, 0.01, 0.01, 0.01,
                                                     0.47, 0.02, 0.47, 0.02, 0.02,
                                                     0.02, 0.47, 0.02, 0.02, 0.47), nrow = 3, ncol = 5, byrow = T),

                                            matrix(c(0.47, 0.02, 0.02, 0.02, 0.47,
                                                     0.02, 0.02, 0.02, 0.47, 0.47,
                                                     0.01, 0.96, 0.01, 0.01, 0.01), nrow = 3, ncol = 5, byrow = T),

                                            matrix(c(0.01, 0.01, 0.96, 0.01, 0.01,
                                                     0.47, 0.02, 0.47, 0.02, 0.02,
                                                     0.02, 0.47, 0.02, 0.02, 0.47), nrow = 3, ncol = 5, byrow = T),

                                            matrix(c(0.01, 0.01, 0.01, 0.96, 0.01,
                                                     0.02, 0.47, 0.02, 0.02, 0.47,
                                                     0.47, 0.02, 0.02, 0.02, 0.47), nrow = 3, ncol = 5, byrow = T),

                                            matrix(c(0.47, 0.47, 0.02, 0.02, 0.02,
                                                     0.47, 0.02, 0.02, 0.47, 0.02,
                                                     0.01, 0.01, 0.01, 0.01, 0.96), nrow = 3, ncol = 5, byrow = T),

                                            matrix(c(0.02, 0.47, 0.47, 0.02, 0.02,
                                                     0.01, 0.01, 0.01, 0.01, 0.96,
                                                     0.47, 0.02, 0.02, 0.02, 0.47), nrow = 3, ncol = 5, byrow = T)),

                          "high" = list(matrix(c(0.01, 0.96, 0.01, 0.01, 0.01,
                                                 0.02, 0.47, 0.47, 0.02, 0.02,
                                                 0.02, 0.47, 0.47, 0.02, 0.02), nrow = 3, ncol = 5, byrow = T),

                                        matrix(c(0.02, 0.47, 0.02, 0.02, 0.47,
                                                 0.01, 0.01, 0.01, 0.01, 0.96,
                                                 0.02, 0.47, 0.02, 0.02, 0.47), nrow = 3, ncol = 5, byrow = T),

                                        matrix(c(0.01, 0.01, 0.96, 0.01, 0.01,
                                                 0.47, 0.02, 0.47, 0.02, 0.02,
                                                 0.47, 0.02, 0.47, 0.02, 0.02), nrow = 3, ncol = 5, byrow = T),

                                        matrix(c(0.02, 0.02, 0.02, 0.47, 0.47,
                                                 0.02, 0.02, 0.02, 0.47, 0.47,
                                                 0.01, 0.01, 0.01, 0.96, 0.01), nrow = 3, ncol = 5, byrow = T),

                                        matrix(c(0.02, 0.47, 0.02, 0.47, 0.02,
                                                 0.01, 0.01, 0.01, 0.96, 0.01,
                                                 0.02, 0.47, 0.02, 0.47, 0.02), nrow = 3, ncol = 5, byrow = T),

                                        matrix(c(0.47, 0.02, 0.47, 0.02, 0.02,
                                                 0.47, 0.02, 0.47, 0.02, 0.02,
                                                 0.96, 0.01, 0.01, 0.01, 0.01), nrow = 3, ncol = 5, byrow = T),

                                        matrix(c(0.47, 0.47, 0.02, 0.02, 0.02,
                                                 0.47, 0.47, 0.02, 0.02, 0.02,
                                                 0.01, 0.96, 0.01, 0.01, 0.01), nrow = 3, ncol = 5, byrow = T),

                                        matrix(c(0.47, 0.02, 0.47, 0.02, 0.02,
                                                 0.01, 0.01, 0.96, 0.01, 0.01,
                                                 0.47, 0.02, 0.47, 0.02, 0.02), nrow = 3, ncol = 5, byrow = T)))

        eps_str <- list("low" = list(matrix(c(-4*1, 1, 1, 1, 1,
                                              1/1.5, -3/2*1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5,
                                              1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5, -3/2*1/1.5), nrow = 3, ncol = 5, byrow = T),

                                     matrix(c(1/1.5, -3/2*1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5,
                                              -4*1, 1, 1, 1, 1,
                                              1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5, -3/2*1/1.5), nrow = 3, ncol = 5, byrow = T),

                                     matrix(c(1, -4*1, 1, 1, 1,
                                              -3/2*1/1.5, 1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5,
                                              1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5, -3/2*1/1.5), nrow = 3, ncol = 5, byrow = T),

                                     matrix(c(-3/2*1/1.5, 1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5,
                                              1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5, -3/2*1/1.5,
                                              1, -4*1, 1, 1, 1), nrow = 3, ncol = 5, byrow = T),

                                     matrix(c(1, 1, -4*1, 1, 1,
                                              -3/2*1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5, 1/1.5,
                                              1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5, -3/2*1/1.5), nrow = 3, ncol = 5, byrow = T),

                                     matrix(c(1, 1, 1, -4*1, 1,
                                              1/1.5, -3/2*1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5,
                                              -3/2*1/1.5, 1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5), nrow = 3, ncol = 5, byrow = T),

                                     matrix(c(1/1.5, -3/2*1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5,
                                              -3/2*1/1.5, 1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5,
                                              1, 1, 1, -4*1, 1), nrow = 3, ncol = 5, byrow = T),

                                     matrix(c(1, 1, 1, 1, -4*1,
                                              1/1.5, -3/2*1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5,
                                              -3/2*1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5, 1/1.5), nrow = 3, ncol = 5, byrow = T)),

                        "moderate" = list(matrix(c(-4*1, 1, 1, 1, 1,
                                                   1/1.5, -3/2*1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5,
                                                   1/1.5, 1/1.5, -3/2*1/1.5, -3/2*1/1.5, 1/1.5), nrow = 3, ncol = 5, byrow = T),

                                          matrix(c(1/1.5, -3/2*1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5,
                                                   -4*1, 1, 1, 1, 1,
                                                   1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5), nrow = 3, ncol = 5, byrow = T),

                                          matrix(c(1, -4*1, 1, 1, 1,
                                                   -3/2*1/1.5, 1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5,
                                                   1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5), nrow = 3, ncol = 5, byrow = T),

                                          matrix(c(-3/2*1/1.5, 1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5,
                                                   1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5, -3/2*1/1.5,
                                                   1, -4*1, 1, 1, 1), nrow = 3, ncol = 5, byrow = T),

                                          matrix(c(1, 1, -4*1, 1, 1,
                                                   -3/2*1/1.5, 1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5,
                                                   1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5), nrow = 3, ncol = 5, byrow = T),

                                          matrix(c(1, 1, 1, -4*1, 1,
                                                   1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5,
                                                   -3/2*1/1.5, 1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5), nrow = 3, ncol = 5, byrow = T),

                                          matrix(c(-3/2*1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5, 1/1.5,
                                                   -3/2*1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5, 1/1.5,
                                                   1, 1, 1, 1, -4*1), nrow = 3, ncol = 5, byrow = T),

                                          matrix(c(1/1.5, -3/2*1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5,
                                                   1, 1, 1, 1, -4*1,
                                                   -3/2*1/1.5, 1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5), nrow = 3, ncol = 5, byrow = T)),

                        "high" = list(matrix(c(1, -4*1, 1, 1, 1,
                                               1/1.5, -3/2*1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5,
                                               1/1.5, -3/2*1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5), nrow = 3, ncol = 5, byrow = T),

                                      matrix(c(1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5,
                                               1, 1, 1, 1, -4*1,
                                               1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5), nrow = 3, ncol = 5, byrow = T),

                                      matrix(c(1, 1, -4*1, 1, 1,
                                               -3/2*1/1.5, 1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5,
                                               -3/2*1/1.5, 1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5), nrow = 3, ncol = 5, byrow = T),

                                      matrix(c(1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5, -3/2*1/1.5,
                                               1/1.5, 1/1.5, 1/1.5, -3/2*1/1.5, -3/2*1/1.5,
                                               1, 1, 1, -4*1, 1), nrow = 3, ncol = 5, byrow = T),

                                      matrix(c(1/1.5, -3/2*1/1.5, 1/1.5, -3/2*1/1.5, 1/1.5,
                                               1, 1, 1, -4*1, 1,
                                               1/1.5, -3/2*1/1.5, 1/1.5, -3/2*1/1.5, 1/1.5), nrow = 3, ncol = 5, byrow = T),

                                      matrix(c(-3/2*1/1.5, 1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5,
                                               -3/2*1/1.5, 1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5,
                                               -4*1, 1, 1, 1, 1), nrow = 3, ncol = 5, byrow = T),

                                      matrix(c(-3/2*1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5, 1/1.5,
                                               -3/2*1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5, 1/1.5,
                                               1, -4*1, 1, 1, 1), nrow = 3, ncol = 5, byrow = T),

                                      matrix(c(-3/2*1/1.5, 1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5,
                                               1, 1, -4*1, 1, 1,
                                               -3/2*1/1.5, 1/1.5, -3/2*1/1.5, 1/1.5, 1/1.5), nrow = 3, ncol = 5, byrow = T)))

        emiss_var <- rep(1, as.numeric(pars[3]))
    }

    # Set the correct ones
    emiss_sim <- emiss_sim[[as.character(pars[5])]]
    eps_str <- eps_str[[as.character(pars[5])]]

    # Set emission distribution with noisiness
    eps <- as.numeric(pars[4])
    emiss_sim <- lapply(seq_along(emiss_sim),function(e, emiss_sim, eps_str, eps) {emiss_sim[[e]] + eps_str[[e]]*eps}, emiss_sim, eps_str, eps)

    # Return model parameters
    return(list(sample_size  = as.numeric(pars[1]),
                n_t          = as.numeric(pars[2]),
                m            = 3,
                n_dep        = as.numeric(pars[3]),
                q_emiss      = rep(5, as.numeric(pars[3])),
                # gamma_var    = 1,
                # emiss_var    = rep(1, as.numeric(pars[3])),
                gamma_var    = gamma_var,
                emiss_var    = emiss_var,
                noisiness    = as.numeric(pars[4]),
                overlapping  = as.character(pars[5]),
                iter         = as.numeric(pars[6]),
                burnin       = as.numeric(pars[7]),
                repetitions  = as.numeric(pars[8]),
                scenario_uid = pars[9],
                uid          = pars[10],
                save_all     = FALSE,
                gamma_sim    = gamma_sim,
                emiss_sim    = emiss_sim[1:as.numeric(pars[3])]))

}
