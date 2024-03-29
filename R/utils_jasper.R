
## Utility functions

#' Convert a list of probabilities into an m x m transition probability matrix
#'
#' @param x an mHMM object
#'
#' @return m x m transition probability matrix, where m is the number of hidden states
#' @export

# get_subject_tpm <- function(x, ...) {
#     UseMethod("get_subject_tpm", x)
# }

get_subject_tpm <- function(x) {
    # Select probs
    p <- apply(x$gamma_prob_bar,2, mean)
    # Select states
    m <- x$input$m
    # Create m x m matrix and return
    matrix(
        p,
        ncol = m,
        nrow = m,
        byrow = TRUE
    )
}

#' Burn function for model output
#'
#' @param x an mHMM object
#'
#' @return mHMM_cont object for which the burn-in samples have been removed
#'         for each parameter.
#' @export

# burn <- function(x, ...) {
#     UseMethod("burn", x)
# }

burn <-  function(x) {
    # Number of burn_in samples
    burn_in <- x$input$burn_in
    J <- x$input$J
    # For each element (and nested elements), remove the burn-in samples
    inp <- x$input
    x$input <- NULL
    # Get data types for each
    dtypes <- vapply(x, function(x) mode(x), "string")
    for(idx in seq_along(x)) {
        # If character, pass
        if(mode(x[[idx]]) == "character") {
            next
            # If label switching or state orders, pass
        } else if(names(x)[idx] == "label_switch") {
            next
        } else if(names(x)[idx] == "state_orders") {
            next
        } else if(names(x)[idx] == "gamma_naccept") {
            next
        } else if(names(x)[idx] == "gamma_paccept") {
            next
        } else if(names(x)[idx] == "emiss_naccept") {
            next
        } else if(names(x)[idx] == "emiss_paccept") {
            next
        } else if(names(x)[idx] == "time") {
            next
        } else if (mode(x[[idx]]) == "list") {
            for(subj_idx in seq_along(x[[idx]])) {
                if(mode(x[[idx]][[subj_idx]]) == "list") {
                    for(n_dep_idx in seq_along(mode(x[[idx]][[subj_idx]]))) {
                        x[[idx]][[subj_idx]][[n_dep_idx]] <- x[[idx]][[subj_idx]][[n_dep_idx]][(burn_in+1):J,]
                    }
                } else {
                    x[[idx]][[subj_idx]] <- x[[idx]][[subj_idx]][(burn_in+1):J,]
                }
            }
        } else if (mode(x[[idx]]) == "numeric") {
            x[[idx]] <- x[[idx]][(burn_in+1):J,]
        }
    }
    # Create new object and return
    x$input <- inp
    class(x) <- "mHMM"
    # Return
    return(x)
}

#' Retrieve MAP estimates for parameters
#'
#' @param x an mHMM object
#'
#' @return List. Maximum a Posteriori (MAP) estimates for each parameter.
#'         names of the elements are identical of the names of the input
#'         parameters
#' @export

# MAP <- function(x, ...) {
#     UseMethod("MAP", x)
# }

MAP <- function(x) {
    # Remove burn-in samples
    feelthebern <- burn(x)
    # Remove input
    feelthebern$input <- NULL
    # Get data types for each
    dtypes <- vapply(feelthebern, function(x) mode(x), "string")
    # Remove character types
    feelthebern <- feelthebern[!dtypes == "character"]
    # For each, collect MAP
    map_out <- vector("list", length(feelthebern))
    # Names
    names(map_out) <- names(feelthebern)
    for(param_idx in seq_along(feelthebern)) {
        if(names(feelthebern)[param_idx] %in% c("label_switch", "input", "state_orders","sample_path","time")) {
            next
        }
        # if numeric, compute MAP
        if (class(feelthebern[[param_idx]]) == "numeric") {
            map_out[[param_idx]][["mean"]] <- unname(mean(feelthebern[[param_idx]]))
            map_out[[param_idx]][["median"]] <- unname(median(feelthebern[[param_idx]]))
            map_out[[param_idx]][["SE"]] <- unname(sd(feelthebern[[param_idx]]))
        } else if(mode(feelthebern[[param_idx]]) == "numeric") {
            map_out[[param_idx]][["mean"]] <- unname(apply(feelthebern[[param_idx]], 2, mean))
            map_out[[param_idx]][["median"]] <- unname(apply(feelthebern[[param_idx]], 2, median))
            map_out[[param_idx]][["SE"]] <- unname(apply(feelthebern[[param_idx]], 2, sd))
        } else {
            if(mode(feelthebern[[param_idx]]) == "list") {
                for(n_subj in seq_along(feelthebern[[param_idx]])) {
                    if(mode(feelthebern[[param_idx]][[n_subj]]) == "list") {
                        map_out[[param_idx]][[n_subj]] <- lapply(feelthebern[[param_idx]][[n_subj]], function(x) {
                            list(
                                "mean" = unname(apply(x, 2, mean)),
                                "median" = unname(apply(x, 2, median)),
                                "SE" = unname(apply(x, 2, sd))
                            )
                        })
                    } else {
                        map_out[[param_idx]][[n_subj]] <- list(
                            "mean" = unname(apply(feelthebern[[param_idx]][[n_subj]], 2, mean)),
                            "median" = unname(apply(feelthebern[[param_idx]][[n_subj]], 2, median)),
                            "SE" = unname(apply(feelthebern[[param_idx]][[n_subj]], 2, sd))
                        )
                    }
                }
            }
        }
    }
    # Add time
    map_out[[which(names(feelthebern) == "time")]] <- feelthebern[[which(names(feelthebern) == "time")]]
    # # Remove label switch
    # map_out$label_switch <- NULL
    # Return
    return(map_out)
}


#' Compute upper and lower values of the 95\% credible interval
#'
#' @param x posterior values for a parameter
#' @param type either one of '0.95' for 95\% CI or '0.99' for 99\% CI
#'
#' @return numeric vector. Element 1 is the lower 95\% CI, element 2 is the upper 95\% CI.
#'
#' @export
credible_interval <- function(x, type=c("0.95", "0.99")) {
    type <- match.arg(type)
    if(type == "0.95") {
        return(unname(apply(x, 2, function(y) quantile(y, c(0.025, 0.975)))))
    } else {
        return(unname(apply(x, 2, function(y) quantile(y, c(0.005, 0.995)))))
    }
}

#' Retrieve CCI estimates for parameters
#'
#' @param x an mHMM object
#'
#' @return List. Confidence credibility intervals (CCIs) estimates for each parameter.
#'         names of the elements are identical of the names of the input
#'         parameters
#' @export

# get_cci <- function(x, ...) {
#     UseMethod("get_cci", x)
# }

get_cci <- function(x) {
    # Remove burn-in samples
    feelthebern <- burn(x)
    # Remove input
    feelthebern$input <- NULL
    # Get data types for each
    dtypes <- vapply(feelthebern, function(x) mode(x), "string")
    # Remove character types
    feelthebern <- feelthebern[!dtypes == "character"]
    # For each, collect MAP
    cci_out <- vector("list", length(feelthebern))
    # Names
    names(cci_out) <- names(feelthebern)
    for(param_idx in seq_along(feelthebern)) {
        if(names(feelthebern)[param_idx] %in% c("label_switch", "input", "state_orders","sample_path","time")) {
            next
        }
        # if numeric, compute MAP
        if (class(feelthebern[[param_idx]]) == "numeric") {
            cci_out[[param_idx]][["CCI_95"]] <- unname(as.vector(quantile(feelthebern[[param_idx]], c(0.0275,0.975))))
            cci_out[[param_idx]][["CCI_99"]] <- unname(as.vector(quantile(feelthebern[[param_idx]], c(0.005, 0.995))))
        } else if(mode(feelthebern[[param_idx]]) == "numeric") {
            cci_out[[param_idx]][["CCI_95"]] <- unname(as.vector(credible_interval(feelthebern[[param_idx]], "0.95")))
            cci_out[[param_idx]][["CCI_99"]] <- unname(as.vector(credible_interval(feelthebern[[param_idx]], "0.99")))
        } else {
            if(mode(feelthebern[[param_idx]]) == "list") {
                for(n_subj in seq_along(feelthebern[[param_idx]])) {
                    if(mode(feelthebern[[param_idx]][[n_subj]]) == "list") {
                        cci_out[[param_idx]][[n_subj]] <- lapply(feelthebern[[param_idx]][[n_subj]], function(x) {
                            list(
                                "CCI_95" = unname(as.vector(credible_interval(x, "0.95"))),
                                "CCI_99" = unname(as.vector(credible_interval(x, "0.99")))
                            )
                        })
                    } else {
                        cci_out[[param_idx]][[n_subj]] <- list(
                            "CCI_95" = unname(as.vector(credible_interval(feelthebern[[param_idx]][[n_subj]], "0.95"))),
                            "CCI_99" = unname(as.vector(credible_interval(feelthebern[[param_idx]][[n_subj]], "0.99")))
                        )
                    }
                }
            }
        }
    }
    # Add time
    cci_out[[which(names(feelthebern) == "time")]] <- feelthebern[[which(names(feelthebern) == "time")]]
    # # Remove label switch
    # cci_out$label_switch <- NULL
    # Return
    return(cci_out)
}
