## Evaluation metrics (% bias etc) go here
# See https://cran.r-project.org/web/packages/rsimsum/vignettes/A-introduction.html
# and https://github.com/JasperHG90/sleepsimReval/blob/master/R/sim_metrics.R

#' Compute the parameter bias of a vector estimated parameters versus the ground-truth value of that parameter (taken from Jasper)
#'
#' @param true_param_value numeric. Value of the ground-truth parameter.
#' @param simulated_param_values k-length numeric vector. k >= 1 and holds the parameter values of the estimated parameters.
#'
#' @return numeric vector. Contains two elements: (1) average bias of simulated values versus the ground-truth value and (2) MCMC SE of the bias value
#'
#' @details This function computes the percentage bias by using the signed mean difference.
#'
#' @seealso Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
#'
#' @export
bias <- function(true_param_value, simulated_param_values) {
    bias <- mean(simulated_param_values - true_param_value)
    bmcse <- bias_MCMC_SE(simulated_param_values)
    ret <- c(bias, bmcse)
    names(ret) <- c("bias", "MCMC_SE")
    return(ret)
}

#' MCMC Standard Error of the bias value
#'
#' @param x k-length numeric vector. k >= 1 and holds the parameter values of the estimated parameters.
#'
#' @seealso Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
#'
#' @return Bias MCMC SE
bias_MCMC_SE <- function(x) {
    nsim <- length(x)
    a <- (1 / (nsim - 1))
    b <- sum((x - mean(x))^2)
    return(
        sqrt((a * b)/(nsim))
    )
}


#' Compute the parameter bias of a vector estimated parameters versus the ground-truth value of that parameter (modified from Jasper)
#'
#' @param true_param_value numeric. Value of the ground-truth parameter.
#' @param simulated_param_values k-length numeric vector. k >= 1 and holds the parameter values of the estimated parameters.
#'
#' @return numeric vector. Contains two elements: (1) average bias of simulated values versus the ground-truth value and (2) MCMC SE of the bias value
#'
#' @details This function computes the percentage bias by using the signed mean difference.
#'
#' @seealso Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
#'
#' @export
rel_bias <- function(true_param_value, simulated_param_values) {
    bias <- mean((simulated_param_values - true_param_value)/abs(true_param_value))
    bmcse <- rel_bias_MCMC_SE(simulated_param_values, true_param_value)
    ret <- c(bias, bmcse)
    names(ret) <- c("rel_bias", "rel_MCMC_SE")
    return(ret)
}

#' MCMC Standard Error of the bias value
#'
#' @param x k-length numeric vector. k >= 1 and holds the parameter values of the estimated parameters.
#'
#' @seealso Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
#'
#' @return Bias MCMC SE
rel_bias_MCMC_SE <- function(x, y) {
    nsim <- length(x)
    a <- (1 / (nsim - 1))
    b <- sum(((x - mean(x))/abs(y))^2)
    return(
        sqrt((a * b)/(nsim))
    )
}


#' Compute the empirical SE
#'
#' @param x numeric vector. Simulated parameter estimates
#'
#' @return numeric vector with two values: (1) empirical standard error and (2) MCMC SE of this value
#'
#' @seealso Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
#'
#' @export
empirical_SE <- function(x) {
    ESE <- sqrt((1/(length(x) - 1)) * sum((x - mean(x))^2))
    MCMCSE <- empirical_MCMC_SE(ESE, length(x))
    ret <- c(ESE, MCMCSE)
    names(ret) <- c("empirical_se", "MCMC_SE")
    return(ret)
}

#' Compute the MCMC SE
#'
#' @param x numeric vector. Simulated parameter estimates.
#'
#' @seealso Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
#'
#' @return numeric scalar. MCMC standard error.
empirical_MCMC_SE <- function(emp_se, nsim) {
    emp_se / sqrt(2*nsim - 1)
}

#' Compute coverage
#'
#' @param CI list. Each element contains the upper and lower values of the 95\% CI of an interation.
#' @param true_param_value numeric scalar. True value of the parameter.
#'
#' @return numeric vector with two values: (1) Probability that the 95\% CI contains the true value and (2) MCMC Standard Error
#'
#' @seealso Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
#'
#' @export
coverage <- function(CI, true_param_value) {
    cvr <- mean(vapply(CI, function(x) (x[1] <= true_param_value) & (true_param_value <= x[2]), 0))
    cvrSE <- coverage_MCMC_SE(cvr, length(CI))
    ret <- c(cvr, cvrSE)
    names(ret) <- c("coverage", "MCMC_SE")
    return(ret)
}

#' Coverage MCMC Standard Error
#'
#' @param coverage scalar. Coverage of the scenario under investigation.
#' @param nsim scalar. Number of iterations used in the scenario.
#'
#' @seealso Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
#'
#' @return MCMC standard error of coverage.
coverage_MCMC_SE <- function(coverage, nsim) {
    sqrt((coverage * (1-coverage))/(nsim))
}

#' Mean Squared Error (MSE)
#'
#' @param simulated_param_values k-length numeric vector. k >= 1 and holds the parameter values of the estimated parameters.
#' @param true_param_value numeric. Value of the ground-truth parameter.
#'
#' @return Mean-Squared Error
#'
#' @seealso Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
#'
#' @export
MSE <- function(simulated_param_values, true_param_value) {
    MSE_est <- mean((simulated_param_values - true_param_value)^2)
    MSE_est_mcmc_se <- MSE_MCMC_SE(MSE_est, simulated_param_values,
                                   true_param_value, length(simulated_param_values))
    ret <- c(MSE_est, MSE_est_mcmc_se)
    names(ret) <- c("MSE", "MCMC_SE")
    return(ret)
}

#' MSE MCMC Standard Error
#'
#' @param MSE Mean-Squared Error. \link[sleepsimReval]{MSE}.
#' @param estimates k-length numeric vector. k >= 1 and holds the parameter values of the estimated parameters.
#' @param true_value numeric. Value of the ground-truth parameter.
#' @param nsim scalar. Number of iterations used in the scenario.
#'
#' @seealso Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
#'
#' @return MCMC Standard Error of MSE.
MSE_MCMC_SE <- function(MSE, estimates, true_value, nsim) {
    a <- sum(((estimates - true_value)^2 - MSE)^2)
    b <- nsim * (nsim - 1)
    return(sqrt(a/b))
}

#' Compute the average model SE
#'
#' @param SE the SD of the posterior distribution for each parameter in the simulation scenario.
#'
#' @return average model SE
#'
#' @seealso Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
#'
#' @details Note that the SE here is not computed around the point estimate. Rather, it is taken from the posterior distribution of the parameter. We input SE (standard error), the formula needs the estimated variance. Hence, the SE is squared in the formula below.
#'
#' @export
average_model_SE <- function(SE) {
    n_sim <- length(SE)
    modSE <- sqrt(1 / n_sim * sum(SE^2))
    modSE_mcmc_se <- average_model_SE_mcmc_se(SE, length(SE), modSE)
    ret <- c(modSE, modSE_mcmc_se)
    names(ret) <- c("modSE", "MCMC_SE")
    return(ret)
}

#' Compute the MCMC SE for the average model SE
#'
#' @param SE the SD of the posterior distribution for each parameter in the simulation scenario.
#' @param nsim scalar. Number of iterations used in the scenario.
#'
#' @seealso Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
#'
#' @return average model SE MCMC SE
average_model_SE_mcmc_se <- function(SE, n_sim, modSE) {
    avg_varest <- mean(SE^2)
    num_a <- (1/(n_sim-1)) * sum((SE^2 - avg_varest) ^ 2)
    den_a <- 4 * n_sim * modSE^2
    return(
        sqrt(num_a / den_a)
    )
}

#' Wrapper function that computes all simulation metrics
#'
#' This wrapper function is compatible with dplyr and can be called in a dplyr chain
#'  of arguments.
#'
#' @param true_values vector of length n. Contains true population values.
#' @param simulated_values vector of length n. Contains simulated values.
#' @param lower_cci vector of length n. Lower 95\% CCI.
#' @param lower_cci vector of length n. Upper 95\% CCI.
#' @param SE vector of length n. Standard Error of the posterior distribution.
#' @param compute_multimodal boolean. If TRUE, this function will compute a test on the simulated parameter values to check if the distribution of parameter estimates is multimodal. See \link[multimode]{modetest}.
#'
#' @importFrom multimode modetest
#'
#' @return data frame containing:
#' \describe{
#'   \item{bias}{percent bias of the simulation scenario, computed as a percentage relative to the true value.}
#'   \item{bias_mcmc_se}{MCMC standard error of bias estimate, computed as a percentage relative to the bias estimate.}
#'   \item{empirical_se}{empirical standard error computed from the simulated values.}
#'   \item{empirical_se_mcmc_se}{MCMC standard error of the empirical SE.}
#'   \item{modSE}{model standard error computed from the simulated values.}
#'   \item{modSE_mcmc_se}{MCMC standard error of the model SE.}
#'   \item{MSE}{Mean-Squared Error of the simulated values.}
#'   \item{MSE_mcmc_se}{MCMC standard error of the simulated values.}
#'   \item{coverage}{Coverage given the 95\% CCI.}
#'   \item{coverage_mcmc_se}{MCMC standard error of coverage.}
#'   \item{bias_corr_coverage}{Bias-adjusted coverage. Instead of using the true population value, use the mean of the simulated values. Useful to check whether poor coverage is the result of bias. See 'See also' for reference.}
#'   \item{bias_corr_coverage_mcmc_se}{MCMC standard error of bias-adjusted coverage.}
#'   \item{multimodal}{p-value of the \link[multimode]{modetest} used to check for multimodal distributions.}
#' }
#'
#' @seealso Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
#'
#' @examples
#' \dontrun{
#' tst <- data %>%
#'    group_by(scenario_id) %>%
#'    # This is how you use the function
#'    do(summarize_simulation_metrics(.$true_values, .$simulated_values
#'                                    .$lower_cci, .$upper_cci, FALSE))
#' }
#'
#' @export
summarize_simulation_scenario <- function(true_values, simulated_values,
                                          lower_cci, upper_cci, SE = NULL,
                                          compute_multimodal = FALSE) {
    out_bias <- unname(bias(true_values[1], simulated_values))
    out_MSE <- unname(MSE(simulated_values, true_values[1]))
    out_ESE <- unname(empirical_SE(simulated_values))
    # Compute coverage
    cci <- map2(lower_cci, upper_cci, function(x,y) c(x, y))
    out_coverage <- unname(coverage(cci, true_values[1]))
    # Compute bias-corrected coverage
    out_coverage_bc <- unname(coverage(cci, mean(simulated_values)))
    df <- data.frame(
        "bias" = out_bias[1],
        "bias_mcmc_se" = out_bias[2],
        "empirical_se" = out_ESE[1],
        "empirical_se_mcmc_se" = out_ESE[2],
        "MSE" = out_MSE[1],
        "MSE_mcmc_se" = out_MSE[2],
        "coverage" = out_coverage[1],
        "coverage_mcmc_se" = out_coverage[2],
        "bias_corr_coverage" = out_coverage_bc[1],
        "bias_corr_coverage_mcmc_se" = out_coverage_bc[2]
    )
    if(compute_multimodal) {
        df$multimodal <- modetest(simulated_values, mod0 = 1)$p.value
    }
    if(!is.null(SE)) {
        out_modSE <- unname(average_model_SE(SE))
        df$modSE = out_modSE[1]
        df$modSE_mcmc_se = out_modSE[2]
    }
    return(df)
}
