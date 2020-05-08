#' Get between subject variance
#'
#' @export
#'

get_var_bar <- function(x){

    # Save number of hidden states
    m <- x[["input"]][["m"]]
    J <- x[["input"]][["J"]]
    n_dep <- x[["input"]][["n_dep"]]
    q_emiss <- x[["input"]][["q_emiss"]]

    # Get between subject variance
    gamma_prob_var <- do.call(rbind,lapply(seq(J), function(j) {
        apply(do.call(rbind,lapply(x[["PD_subj"]], "[",j,paste0("S", rep(1:m, each = m), "toS", rep(1:m, m)))),2,var)
    }))
    gamma_int_var <- do.call(rbind,lapply(seq(J), function(j) {
        apply(do.call(rbind,lapply(x[["gamma_int_subj"]], "[",j,)),2,var)
    }))

    # Emiss
    emiss_prob_var <- lapply(seq(n_dep), function(q) {
        do.call(rbind, lapply(seq(J), function(j) {
            apply(do.call(rbind,lapply(x[["PD_subj"]], "[",j,paste("q", q, "_emiss", rep(1:q_emiss[q], m), "_S", rep(1:m, each = q_emiss[q]), sep = ""))),2,var)
        }))
    })
    emiss_int_var <- lapply(seq(n_dep), function(q) {
        do.call(rbind, lapply(seq(J), function(j) {
            apply(do.call(rbind,lapply(lapply(x[["emiss_int_subj"]], "[[",q), "[",j,)),2,var)
        }))
    })

    return(list(gamma_prob_var = gamma_prob_var, emiss_prob_var = emiss_prob_var,
                gamma_int_var = gamma_int_var, emiss_int_var = emiss_int_var))

}
