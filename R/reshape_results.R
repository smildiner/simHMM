# Classes and method for sleepsimR results

# From Jasper Ginn; see: <https://github.com/JasperHG90/sleepsimReval/blob/master/R/results.R>


#' Postprocessing utility function for parameter estimates
#'
#' @param z list. contains the MAP values of the parameter estimates.
#' @param m integer. Number of hidden states
#'
#' @return data frame in which each parameter estimate has been placed in its own column.
#'
#' @export
postprocess_param_est <- function(z, m) {
    # Create names
    nams <- paste0("state", 1:m)
    # n dep
    n_dep <- length(z)
    # Add to each
    for(idx in seq_along(z)) {
        names(z[[idx]]$mean) <- nams
        names(z[[idx]]$median) <- nams
        names(z[[idx]]$SE) <- nams
    }
    # To data frame
    df <- as.data.frame(z)
    # Replace periods by underscores
    colnames(df) <- gsub("\\.", "_", colnames(df))
    # Return
    return(df)
}

#' Postprocess utility function for gamma_prob_bar
#'
#' @param z list. contains the MAP values of the transition probabilities.
#' @param m integer. Number of hidden states
#'
#' @return data frame in which each parameter estimate has been placed in its own column.
#'
#' @export
postprocess_gamma_int <- function(z, m) {
    # Number of values is equal to m x (m-1)
    # Make names
    nams <- c()
    for(idx_col in 1:m) {
        for(idx_row in 1:m) {
            nams <- c(nams, paste0("S",idx_col, "toS", idx_row))
        }
    }
    # Subset mean
    out <- c(z$median, z$SE)
    # Ignore SE --> names
    names(out) <- c(paste0(nams, "_median"), paste0(nams, "_SE"))
    # data frame and return
    return(data.frame(out))
}

#' Postprocess credible intervals
#'
#' @param z list. contains the credible intervals of the transition probabilities.
#' @param m integer. Number of hidden states
#'
#' @return data frame in which each CCI has been placed in its own column. (with _lower and _upper appended).
#'
#' @export
postprocess_ci <- function(z, m) {
    # Create names
    out_mp <- vector("list", length(z))
    for(lst_idx in seq_along(z)) {
        tmp <- z[[lst_idx]]
        nm <- names(z)[lst_idx]
        if(nm == "gamma_prob_bar") {
            # Make names
            nams <- c()
            for(idx_col in 1:m) {
                for(idx_row in 1:m) {
                    for(rngnm in c("lower", "upper")) {
                        nams <- c(nams, paste0("gamma_prob_bar_S",idx_col, "toS", idx_row, "_", rngnm))
                    }
                }
            }
            # Add names
            names(tmp) <- nams
            # To data frame
            out_mp[[lst_idx]] <- as.data.frame(tmp)
        } else {
            for(var_idx in seq_along(tmp)) {
                nms <- c()
                for(ele_idx in 1:(length(tmp[[var_idx]])/2)) {
                    nms <- c(nms, c(paste0("state", ele_idx,"_lower"), paste0("state", ele_idx,"_upper")))
                }
                names(tmp[[var_idx]]) <- nms
            }
            tmpdf <- as.data.frame(tmp)
            colnames(tmpdf) <- gsub("\\.", "_", colnames(tmpdf))
            out_mp[[lst_idx]] <-tmpdf
        }
    }
    # Cbind
    return(
        do.call(cbind.data.frame, out_mp)
    )
}

#' Postprocess subject-specific metrics
#'
#' @param z list. contains the MAP values of the parameter estimates for each subject.
#' @param m integer. Number of hidden states
#'
#' @return data frame with rows equal to the number of subjects and columns equal to the number of states x dependent variables x
#'
#' @export
postprocess_subject_specific <- function(z, m) {
    # Number of subjects
    subjs <- paste0("subject_", 1:length(z))
    # Depvar/state names
    ## MODIFY TO ACCEPT VARIABLE NUMBER OF DEPENDENT VARIABLES
    depvars <- c("dep_1", "EOG_median_theta", "EOG_min_beta")
    # States
    states <- c("state_1", "state_2", "state_3")
    # Expand grid
    grid <- expand.grid(states, depvars)
    # Paste
    statenames <- paste0(grid$Var2, "_", grid$Var1)
    # For each subject, make data frame
    subj_out <- vector("list", length(z))
    names(subj_out) <- subjs
    for(idx in seq_along(subjs)) {
        subj_est <- z[[idx]]
        for(est_idx in seq_along(subj_est)) {
            names(subj_est[[est_idx]]) <- statenames
        }
        # Bind
        subj_est <- do.call(cbind.data.frame, subj_est)
        # Replace periods
        colnames(subj_est) <- gsub("\\.", "_", colnames(subj_est))
        # Add to results
        subj_out[[idx]] <- subj_est
    }
    # Bind
    b <- do.call(rbind.data.frame, subj_out)
    # Unname rows
    row.names(b) <- 1:nrow(b)
    # Return
    return(b)
}
