#' Fit mHMM using mHMMfast
#'
#' @export

fit_mHMM <- function(m,
                     n_dep,
                     q_emiss,
                     gamma,
                     emiss,
                     iter,
                     burnin,
                     start_gamma = NULL,
                     start_emiss = NULL,
                     data_sim) {

    # Set starting values
    # gamma_start
    if (is.null(start_gamma)) {

        start_gamma <- diag(runif(1, 0.7, 0.9), m)
        start_gamma[lower.tri(start_gamma) | upper.tri(start_gamma)] <- (1 - diag(start_gamma)) / (m - 1)

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
    out <- simHMM::mHMMfast(s_data = data_sim$obs,
                    gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                    # xx = xx_vec,
                    start_val = c(list(start_gamma), start_emiss),
                    mcmc = list(J = iter, burn_in = burnin),
                    return_path = FALSE,
                    show_progress = FALSE)
    out[["time"]] <- Sys.time() - ti

    # Return model output
    return(out)

}
