#----------------------#
#  Auxiliary functions :
#----------------------#

#' @keywords internal
# Calculates the probabilities of observing each state at each point in time given
# the observations of all dependent variables, used for the forward probabilities
# Based on Zuchini 2016.
all1 <- function(x, emiss, n_dep, data_distr){
    inp <- rep(list(NULL), n_dep)
    if(data_distr == "categorical"){
        for(q in 1:n_dep){
            inp[[q]] <- t(emiss[[q]][,x[,q]])
        }
    } else if (data_distr == "continuous"){
        for(q in 1:n_dep){
            inp[[q]] <- outer(x[,q], Y = emiss[[q]][,1], FUN = dnorm, sd = rep(sqrt(emiss[[q]][,2]), each = dim(x)[1]))
        }
    } else if (data_distr == "poisson") {
        for(q in 1:n_dep){
            inp[[q]] <- outer(x[,q], emiss[[q]][,1], FUN = dpois)
        }
    }
    allprobs <- Reduce("*", inp)
    return(allprobs)
}

#' @keywords internal
# Could maybe made external
# Calculates the forward probabilities, used for sampling the state sequence
# Based on Zuchini 2016.
cat_mult_fw_r_to_cpp <- function(x, m, emiss, n_dep, gamma, delta = NULL){
    if(is.null(delta)) {
        delta <- solve(t(diag(m) - gamma + 1), rep(1, m))
    }
    n        <- dim(x)[1]
    allprobs <- all1(x = x, emiss = emiss, n_dep = n_dep, data_distr = "categorical")
    out <- cat_mult_fw_cpp(allprobs = allprobs, gamma = gamma, m = m, n = n, delta = delta)
    return(out)
}

#' @keywords internal
# Could maybe made external
# Calculates the forward probabilities, used for sampling the state sequence
# Based on Zuchini 2016.
cont_mult_fw_r_to_cpp <- function(x, m, emiss, n_dep, gamma, delta = NULL){
    if(is.null(delta)) {
        delta <- solve(t(diag(m) - gamma + 1), rep(1, m))
    }
    n        <- dim(x)[1]
    allprobs <- all1(x = x, emiss = emiss, n_dep = n_dep, data_distr = "continuous")
    out <- cat_mult_fw_cpp(allprobs = allprobs, gamma = gamma, m = m, n = n, delta = delta)
    return(out)
}

#' @keywords internal
# Could maybe made external
# Calculates the forward probabilities, used for sampling the state sequence
# Based on Zuchini 2016.
pois_mult_fw_r_to_cpp <- function(x, m, emiss, n_dep, gamma, delta = NULL){
    if(is.null(delta)) {
        delta <- solve(t(diag(m) - gamma + 1), rep(1, m))
    }
    n        <- dim(x)[1]
    allprobs <- all1(x = x, emiss = emiss, n_dep = n_dep, data_distr = "poisson")
    out <- cat_mult_fw_cpp(allprobs = allprobs, gamma = gamma, m = m, n = n, delta = delta)
    return(out)
}

#' @keywords internal
# Obtain fractional log likelihood for multinomial intercept only model, bassed on P. Rossi 2004
llmnl_int_frac <- function(beta, Obs, n_cat, pooled_likel, w, wgt){
    return((1 - w) * llmnl_int(beta = beta, Obs = Obs, n_cat = n_cat) + w * wgt * pooled_likel)
}

#' @keywords internal
# Obtain mnl -Expected[Hessian]  for intercept only model, bassed on P.Rossi 2004
mnlHess_int <- function(int, Obs, n_cat){
    n_Obs 	<- length(Obs)
    betas   <- matrix(c(0, int), byrow = T, ncol = n_cat)
    prob    <- exp(betas) / sum(exp(betas))
    Hess    <- (diag(x = prob[-1], nrow = n_cat-1) - outer(prob[-1],prob[-1])) * n_Obs
    return(Hess)
}

#' @keywords internal
# one run of the random walk metropolis sampler for an intercept only multinomial distribution
# this means no covariates at the lower/time level
mnl_RW_once <- function(int1, Obs, n_cat, mu_int_bar1, V_int1, scalar, candcov1) {
    # obtain likelihood and transition prob with the parameters sampled in the previous iteration and current sampled state sequence
    oldloglike	 		<- llmnl_int(beta = int1, Obs = Obs, n_cat = n_cat)
    oldpostlike	 		<- oldloglike + dmvnorm(int1, mu_int_bar1, V_int1, log = TRUE)
    probold				  <- int_to_prob(int1)

    # obtain new parameters for gamma from proposal distribution plus new likelihood
    int_new		 		  <- int1 + rmvnorm(1, rep(0, (n_cat - 1)), scalar^2 * candcov1, method = "svd")
    newloglike	 		<- llmnl_int(beta = int_new, Obs = Obs, n_cat = n_cat)
    newpostlike	 		<- newloglike + dmvnorm(int_new, mu_int_bar1, V_int1, log = TRUE)
    probnew				  <- int_to_prob(int_new)

    # determine to use the updated or current (previous iteration) gamma values of the parameters
    acc 				   <- min(log(1), (newpostlike - oldpostlike))
    if(acc < log(1)) {
        unif         <- log(runif(1))
    } else {
        unif         <- log(1)
    }
    if (unif <= acc) {
        draw_int		<- int_new
        accept			<- 1
        prob			  <- probnew
    } else {
        draw_int		<- int1
        accept			<- 0
        prob			  <- probold
    }
    return(list(draw_int = draw_int, accept = accept, prob = prob))
}

#------------------------------#
#  Changed auxiliary functions :
#------------------------------#

#' @keywords internal
# simple functions used in mHMM
dif_matrix <- function(rows, cols){
    return(matrix(, ncol = cols, nrow = rows))
}

#' @keywords internal
nested_list <- function(n_dep, m){
    return(rep(list(vector("list", n_dep)),m))
}

#' @keywords internal
dif_vector <- function(x){
    return(numeric(x))
}

#' @keywords internal
is.whole <- function(x) {
    return(is.numeric(x) && floor(x) == x)
}

#' @keywords internal
is.mHMM <- function(x) {
    inherits(x, "mHMM")
}

#' @keywords internal
is.mHMM_gamma <- function(x) {
    inherits(x, "mHMM_gamma")
}

#' @keywords internal
hms <- function(t){
    paste(formatC(t %/% (60*60) %% 24, width = 2, format = "d", flag = "0"),
          formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0"),
          formatC(t %% 60, width = 2, format = "d", flag = "0"),
          sep = ":")
}

#' @keywords internal
# computes probabilities from intercepts
#' @export
int_to_prob <- function(int1) {
    if(is.matrix(int1)){
        prob1 <- matrix(nrow = nrow(int1), ncol = ncol(int1) + 1)
        for(r in 1:nrow(int1)){
            exp_int1 	<- matrix(exp(c(0, int1[r,])), nrow  = 1)
            prob1[r,] <- exp_int1 / as.vector(exp_int1 %*% c(rep(1, (dim(exp_int1)[2]))))
        }
    } else {
        exp_int1 	<- matrix(exp(c(0, int1)), nrow  = 1)
        prob1 		<- exp_int1 / as.vector(exp_int1 %*% c(rep(1, (dim(exp_int1)[2]))))
    }
    return(round(prob1,4))
}

#' @keywords internal
# computes intercepts from probabilities, per row of input matrix
# first catagory is reference catagory
#' @export
prob_to_int <- function(prob1){
    prob1 <- prob1 + 0.00001
    b0 <- matrix(NA, nrow(prob1), ncol(prob1)-1)
    sum_exp <- numeric(nrow(prob1))
    for(r in 1:nrow(prob1)){
        sum_exp[r] <- (1/prob1[r,1]) - 1
        for(cr in 2:ncol(prob1)){
            #for every b0 except the first collumn (e.g. b012 <- log(y12/y11-y12))
            b0[r,(cr-1)] <- log(prob1[r,cr]*(1+sum_exp[r]))
        }
    }
    return(round(b0,4))
}
