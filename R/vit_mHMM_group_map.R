#' Obtain hidden state sequence for each subject using the Viterbi
#' algorithm
#' @export

vit_mHMM_group_map <- function(object, s_data, data_distr){
    if(!(data_distr %in% c("categorical","continuous","poisson"))){
        stop("The data distribution is not categorical, continuous or poisson.")
    }
    id         <- unique(s_data[,1])
    n_subj     <- length(id)
    if(length(object$PD_subj) != n_subj){
        stop("s_data used should be identical to the data used for creating the object in mHMM.
         The number of subjects in the datasets are not the same.")
    }
    m           <- length(object[["gamma_prob_bar"]][["median"]])
    n_vary      <- table(s_data[,1])
    max_n       <- max(n_vary)
    state_seq   <- matrix(NA_integer_,ncol = n_subj, nrow = max_n)
    state_probs <- vector("list", length = n_subj)
    for(s in 1:n_subj){
        state_probs[[s]] <- matrix(NA_real_, ncol = m, nrow = n_vary[s])
    }

    # input      <- object$input
    n_dep      <- ncol(s_data) - 1
    m          <- sqrt(length(object$gamma_prob_bar[["median"]]))
    if(data_distr == "categorical"){
        q_emiss <- sapply(object$emiss_prob_bar, function(q) (length(q[["median"]])/m) )
    }

    # if(is.null(burn_in)){
    #     burn_in  <- input$burn_in
    # }
    # J          <- input$J

    # if (burn_in >= (J-1)){
    #     stop(paste("The specified burn in period should be at least 2 points smaller
    #            compared to the number of iterations J, J =", J))
    # }

    if(data_distr == "categorical"){
        est_emiss  <- rep(list(lapply(q_emiss, dif_matrix, rows = m)), n_subj)
        start <- c(0, q_emiss * m)
        for(i in 1:n_subj){
            for(j in 1:n_dep){
                est_emiss[[i]][[j]][] <- matrix(object$PD_subj[[i]][["median"]][(sum(start[1:j]) + 1) : sum(start[1:(j+1)])],
                                                byrow = TRUE, ncol = q_emiss[j], nrow = m)
            }
        }
        est_gamma <- lapply(object$gamma_int_subj, function(s) int_to_prob(matrix(s[["median"]], nrow = m, byrow = TRUE)) )
        for(s in 1:n_subj){
            emiss <- est_emiss[[s]]
            gamma    <- est_gamma[[s]]
            probs    <- cat_mult_fw_r_to_cpp(x = as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep),
                                             m = m, emiss = emiss, n_dep = n_dep, gamma = gamma, delta = NULL)[[1]]
            # $forward_p
            state_seq[1:n_vary[s], s] <- apply(probs, 2, which.max)
            state_probs[[s]] <- t(probs)
        }
    } else if(data_distr == "continuous"){
        est_emiss  <- rep(list(rep(list(matrix(NA_real_,nrow = m, ncol = 2)),n_dep)), n_subj)
        for(s in 1:n_subj){
            for(q in 1:n_dep){
                # est_emiss[[s]][[q]][] <- matrix(object$PD_subj[[s]][["median"]][(1+(q-1)*m) : (m+(q-1)*m)],
                #                                 byrow = TRUE, nrow = m, ncol = 1)
                est_emiss[[s]][[q]][] <- matrix(c(object$PD_subj[[s]][["median"]][((q-1) * m + 1) : ((q-1) * m + m)],
                                                  object$PD_subj[[s]][["median"]][(n_dep * m + (q-1) * m + 1) : (n_dep * m + (q-1) * m + m)]),
                                                ncol = 2, nrow = m, byrow = FALSE)

            }
        }
        est_gamma <- lapply(object$gamma_int_subj, function(s) int_to_prob(matrix(s[["median"]], nrow = m, byrow = TRUE)) )
        for(s in 1:n_subj){
            emiss   <- est_emiss[[s]]
            gamma   <- est_gamma[[s]]
            probs   <- cont_mult_fw_r_to_cpp(x = as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep),
                                             m = m, emiss = emiss, n_dep = n_dep, gamma = gamma)[[1]]
            state_seq[1:n_vary[s], s] <- apply(probs, 2, which.max)
            state_probs[[s]] <- t(probs)
        }
    } else if(data_distr == "poisson"){
        est_emiss  <- rep(list(rep(list(matrix(NA_real_,nrow = m, ncol = 1)),n_dep)), n_subj)
        for(s in 1:n_subj){
            for(q in 1:n_dep){
                est_emiss[[s]][[q]][] <- matrix(exp(object$emiss_mu_bar[[q]][["median"]]),
                                                byrow = TRUE, nrow = m, ncol = 1)
            }
        }
        est_gamma <- int_to_prob(matrix(object$gamma_int_bar[["median"]], nrow = m, byrow = TRUE) )
        for(s in 1:n_subj){
            emiss   <- est_emiss[[s]]
            gamma   <- est_gamma
            probs   <- pois_mult_fw_r_to_cpp(x = as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep),
                                             m = m, emiss = emiss, n_dep = n_dep, gamma = gamma)[[1]]
            state_seq[1:n_vary[s], s] <- apply(probs, 2, which.max)
            state_probs[[s]] <- t(probs)
        }
    }
    colnames(state_seq) <- paste0("subj_", id)
    names(state_probs) <- paste0("subj_",id)
    return(list("state_seq" = state_seq, "state_probs" = state_probs))
}

