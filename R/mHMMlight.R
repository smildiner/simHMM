#' Get estimations with mHMMlight
#'
# @param parameter
#'
#' @return #' @return \code{mHMM} returns an object of class \code{mHMM}, which has
#'   \code{print} and \code{summary} methods to see the results.
#'
# Package imports
# not sure if all functions given below for packages are actually still used, check!
#' @importFrom mvtnorm dmvnorm rmvnorm dmvt rmvt
#' @importFrom MCMCpack rdirichlet rwish
#' @importFrom stats optim rnorm runif median
# for RCpp
#' @useDynLib simHMM, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @export

mHMMlight <- function(s_data, gen, xx = NULL, start_val, mcmc, return_path = FALSE, print_iter, show_progress = TRUE,
                     gamma_hyp_prior = NULL, emiss_hyp_prior = NULL, gamma_sampler = NULL, emiss_sampler = NULL,
                     subj_data = FALSE){

    if(!missing(print_iter)){
        warning("The argument print_iter is deprecated; please use show_progress instead to show the progress of the algorithm.")
    }
    # Initialize data -----------------------------------
    # dependent variable(s), sample size, dimensions gamma and conditional distribuiton
    n_dep			 <- gen$n_dep
    dep_labels <- colnames(s_data[,2:(n_dep+1)])
    id         <- unique(s_data[,1])
    n_subj     <- length(id)
    subj_data  <- rep(list(NULL), n_subj)
    if(sum(sapply(s_data, is.factor)) > 0 ){
        stop("Your data contains factorial variables, which cannot be used as input in the function mHMM. All variables have to be numerical.")
    }
    for(s in 1:n_subj){
        subj_data[[s]]$y <- as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep)
    }
    ypooled    <- n <- NULL
    n_vary     <- numeric(n_subj)
    m          <- gen$m
    q_emiss 		 <- gen$q_emiss
    emiss_int_mle <- rep(list(NULL), n_dep)
    emiss_mhess   <- rep(list(NULL), n_dep)
    for(q in 1:n_dep){
        emiss_int_mle[[q]] <- matrix(, m, (q_emiss[q] - 1))
        emiss_mhess[[q]] <- matrix(, (q_emiss[q] - 1) * m, (q_emiss[q] - 1))
    }
    for(s in 1:n_subj){
        ypooled   <- rbind(ypooled, subj_data[[s]]$y)
        n         <- dim(subj_data[[s]]$y)[1]
        n_vary[s] <- n
        subj_data[[s]]	<- c(subj_data[[s]], n = n, list(gamma_converge = numeric(m), gamma_int_mle = matrix(, m, (m - 1)),
                                                        gamma_mhess = matrix(, (m - 1) * m, (m - 1)), emiss_converge =
                                                            rep(list(numeric(m)), n_dep), emiss_int_mle = emiss_int_mle, emiss_mhess = emiss_mhess))
    }
    n_total 		<- dim(ypooled)[1]

    # covariates
    n_dep1 <- 1 + n_dep
    nx <- numeric(n_dep1)
    if (is.null(xx)){
        xx <- rep(list(matrix(1, ncol = 1, nrow = n_subj)), n_dep1)
        nx[] <- 1
    } else {
        if(!is.list(xx) | length(xx) != n_dep1){
            stop("If xx is specified, xx should be a list, with the number of elements equal to the number of dependent variables + 1")
        }
        for(i in 1:n_dep1){
            if (is.null(xx[[i]])){
                xx[[i]] <- matrix(1, ncol = 1, nrow = n_subj)
                nx[i] <- 1
            } else {
                nx[i] <- ncol(xx[[i]])
                if (sum(xx[[i]][,1] != 1)){
                    stop("If xx is specified, the first column in each element of xx has to represent the intercept. That is, a column that only consists of the value 1")
                }
                if(nx[i] > 1){
                    for(j in 2:nx[i]){
                        if(is.factor(xx[[i]][,j])){
                            stop("Factors currently cannot be used as covariates, see help file for alternatives")
                        }
                        if((length(unique(xx[[i]][,j])) == 2) & (sum(xx[[i]][,j] != 0 & xx[[i]][,j] !=1) > 0)){
                            stop("Dichotomous covariates in xx need to be coded as 0 / 1 variables. That is, only conisting of the values 0 and 1")
                        }
                        if(length(unique(xx[[i]][,j])) > 2){
                            xx[[i]][,j] <- xx[[i]][,j] - mean(xx[[i]][,j])
                        }
                    }
                }
            }
        }
    }

    # Initialize mcmc argumetns
    J 				<- mcmc$J
    burn_in			<- mcmc$burn_in


    # Initalize priors and hyper priors --------------------------------
    # Initialize gamma sampler
    if(is.null(gamma_sampler)) {
        gamma_int_mle0  <- rep(0, m - 1)
        gamma_scalar    <- 2.93 / sqrt(m - 1)
        gamma_w         <- .1
    } else {
        gamma_int_mle0  <- gamma_sampler$gamma_int_mle0
        gamma_scalar    <- gamma_sampler$gamma_scalar
        gamma_w         <- gamma_sampler$gamma_w
    }

    # Initialize emiss sampler
    if(is.null(emiss_sampler)){
        emiss_int_mle0 <- rep(list(NULL), n_dep)
        emiss_scalar	<- rep(list(NULL), n_dep)
        for(q in 1:n_dep){
            emiss_int_mle0[[q]] <- rep(0, q_emiss[q] - 1)
            emiss_scalar[[q]] 	<- 2.93 / sqrt(q_emiss[q] - 1)
        }
        emiss_w		<- .1
    } else {
        emiss_int_mle0	<- emiss_sampler$emiss_int_mle0
        emiss_scalar 	<- emiss_sampler$emiss_scalar
        emiss_w    		<- emiss_sampler$emiss_w
    }

    # Initialize Gamma hyper prior
    if(is.null(gamma_hyp_prior)){
        gamma_mu0	  <- rep(list(matrix(0,nrow = nx[1], ncol = m - 1)), m)
        gamma_K0			<- diag(1, nx[1])
        gamma_nu			<- 3 + m - 1
        gamma_V			  <- gamma_nu * diag(m - 1)
    } else {
        ###### BUILD in a warning / check if gamma_mu0 is a matrix when given, with  nrows equal to the number of covariates
        gamma_mu0			<- gamma_hyp_prior$gamma_mu0
        gamma_K0			<- gamma_hyp_prior$gamma_K0
        gamma_nu			<- gamma_hyp_prior$gamma_nu
        gamma_V			  <- gamma_hyp_prior$gamma_V
    }

    # Initialize Pr hyper prior
    # emiss_mu0: for each dependent variable, emiss_mu0 is a list, with one element for each state.
    # Each element is a matrix, with number of rows equal to the number of covariates (with the intercept being one cov),
    # and the number of columns equal to q_emiss[q] - 1.
    emiss_mu0	  <- rep(list(vector("list", m)), n_dep)
    emiss_nu	    <- rep(list(NULL), n_dep)
    emiss_V	    <- rep(list(NULL), n_dep)
    emiss_K0     <- rep(list(NULL), n_dep)
    if(is.null(emiss_hyp_prior)){
        for(q in 1:n_dep){
            for(i in 1:m){
                emiss_mu0[[q]][[i]]		<- matrix(0, ncol = q_emiss[q] - 1, nrow = nx[1 + q])
            }
            emiss_nu[[q]]		<- 3 + q_emiss[q] - 1
            emiss_V[[q]]		  <- emiss_nu[[q]] * diag(q_emiss[q] - 1)
            emiss_K0[[q]]		<- diag(1, nx[1 + q])
        }
    } else {
        for(q in 1:n_dep){
            # emiss_hyp_prior[[q]]$emiss_mu0 has to contain a list with lenght equal to m, and each list contains matrix with number of rows equal to number of covariates for that dep. var.
            # stil build in a CHECK, with warning / stop / switch to default prior
            emiss_mu0[[q]]	 <- emiss_hyp_prior$emiss_mu0[[q]]
            emiss_nu[[q]]	 <- emiss_hyp_prior$emiss_nu[[q]]
            emiss_V[[q]]		 <- emiss_hyp_prior$emiss_V[[q]]
            emiss_K0[[q]]	 <- diag(emiss_hyp_prior$emiss_K0, nx[1 + q])
        }
    }


    # Define objects used to store data in mcmc algorithm, not returned ----------------------------
    # overall
    c <- llk <- numeric(1)
    sample_path <- lapply(n_vary, dif_matrix, cols = J)
    trans <- rep(list(vector("list", m)), n_subj)

    # gamma
    gamma_int_mle_pooled <- gamma_pooled_ll <- vector("list", m)
    gamma_c_int <- rep(list(matrix(, n_subj, (m-1))), m)
    gamma_mu_int_bar <- gamma_V_int <- vector("list", m)
    gamma_mu_prob_bar <- rep(list(numeric(m)), m)
    gamma_naccept <- matrix(0, n_subj, m)

    # emiss
    cond_y <- lapply(rep(n_dep, n_subj), nested_list, m = m)
    emiss_int_mle_pooled <- emiss_pooled_ll <- rep(list(vector("list", n_dep)), m)
    emiss_c_int <- rep(list(lapply(q_emiss - 1, dif_matrix, rows = n_subj)), m)
    emiss_mu_int_bar <- emiss_V_int <- rep(list(vector("list", n_dep)), m)
    emiss_mu_prob_bar <- rep(list(lapply(q_emiss, dif_vector)), m)
    emiss_naccept <- rep(list(matrix(0, n_subj, m)), n_dep)


    # Define objects that are returned from mcmc algorithm ----------------------------
    # Define object for subject specific posterior density, put start values on first row
    if(length(start_val) != n_dep + 1){
        stop("The number of elements in the list start_val should be equal to 1 + the number of dependent variables,
             and should not contain nested lists (i.e., lists within lists)")
    }
    PD 					  <- matrix(, nrow = J, ncol = sum(m * q_emiss) + m * m + 1)
    PD_emiss_names   <- paste("q", 1, "_emiss", rep(1:q_emiss[1], m), "_S", rep(1:m, each = q_emiss[1]), sep = "")
    if(n_dep > 1){
        for(q in 2:n_dep){
            PD_emiss_names <- c(PD_emiss_names, paste("q", q, "_emiss", rep(1:q_emiss[q], m), "_S", rep(1:m, each = q_emiss[q]), sep = ""))
        }
    }
    colnames(PD) 	<- c(PD_emiss_names, paste("S", rep(1:m, each = m), "toS", rep(1:m, m), sep = ""), "LL")

    PD[1, ((sum(m * q_emiss) + 1)) :((sum(m * q_emiss) + m * m))] <- unlist(sapply(start_val, t))[1:(m*m)]
    PD[1, 1:((sum(m * q_emiss)))] <- unlist(sapply(start_val, t))[(m*m + 1): (m*m + sum(m * q_emiss))]

    if(subj_data){
        PD_subj				<- rep(list(PD), n_subj)
    }

    # Define object for population posterior density (probabilities and regression coefficients parameterization )
    gamma_prob_bar		<- matrix(, nrow = J, ncol = (m * m))
    colnames(gamma_prob_bar) <- paste("S", rep(1:m, each = m), "toS", rep(1:m, m), sep = "")
    gamma_prob_bar[1,] <- PD[1,(sum(m*q_emiss) + 1):(sum(m * q_emiss) + m * m)]
    emiss_prob_bar			<- lapply(q_emiss * m, dif_matrix, rows = J)
    names(emiss_prob_bar) <- dep_labels
    for(q in 1:n_dep){
        colnames(emiss_prob_bar[[q]]) <- paste("Emiss", rep(1:q_emiss[q], m), "_S", rep(1:m, each = q_emiss[q]), sep = "")
        start <- c(0, q_emiss * m)
        emiss_prob_bar[[q]][1,] <- PD[1,(sum(start[1:q]) + 1):(sum(start[1:q]) + (m * q_emiss[q]))]
    }
    gamma_int_bar				<- matrix(, nrow = J, ncol = ((m-1) * m))
    colnames(gamma_int_bar) <- paste("int_S", rep(1:m, each = m-1), "toS", rep(2:m, m), sep = "")
    gamma_int_bar[1,] <- as.vector(prob_to_int(matrix(gamma_prob_bar[1,], byrow = TRUE, ncol = m, nrow = m)))

    gamma_V_int_bar <- matrix(, nrow = J, ncol = ((m-1) * (m-1) * m))
    colnames(gamma_V_int_bar) <- paste("var_int_S", rep(1:m, each = (m-1)*(m-1)), "toS", rep(2:m, each=m-1), "_with_", "int_S", rep(1:m, each = (m-1)*(m-1)), "toS", rep(2:m, m), sep = "")
    gamma_V_int_bar[1,] <- unlist(lapply(gamma_V, function(e) as.vector(t(e))))

    if(nx[1] > 1){
        gamma_cov_bar				<- matrix(, nrow = J, ncol = ((m-1) * m) * (nx[1] - 1))
        colnames(gamma_cov_bar) <- paste( paste("cov", rep(1 : (nx[1] - 1),each = nx[1]-1), "_", sep = ""), "S", rep(1:m, each = (m-1) * (nx[1] - 1)), "toS", rep(2:m, m * (nx[1] - 1)), sep = "")
        gamma_cov_bar[1,] <- 0
    } else{
        gamma_cov_bar <- "No covariates where used to predict the transition probability matrix"
    }
    emiss_int_bar			<- lapply((q_emiss-1) * m, dif_matrix, rows = J)
    names(emiss_int_bar) <- dep_labels

    emiss_V_int_bar <- lapply((q_emiss-1) * (q_emiss-1) * m, dif_matrix, rows = J)
    names(emiss_V_int_bar) <- dep_labels

    for(q in 1:n_dep){
        colnames(emiss_int_bar[[q]]) <-  paste("int_Emiss", rep(2:q_emiss[q], m), "_S", rep(1:m, each = q_emiss[q] - 1), sep = "")
        emiss_int_bar[[q]][1,] <- as.vector(prob_to_int(matrix(emiss_prob_bar[[q]][1,], byrow = TRUE, ncol = q_emiss[q], nrow = m)))

        emiss_V_int_bar[[q]] <- matrix(, nrow = J, ncol = ((q_emiss[q]-1) * (q_emiss[q]-1) * m))
        colnames(emiss_V_int_bar[[q]]) <- paste("var_int_Emiss", rep(2:q_emiss[q], each = (q_emiss[q]-1)),"_with_Emiss",rep(2:q_emiss[q], (q_emiss[q]-1)), "_S", rep(1:m, each = (q_emiss[q] - 1)*(q_emiss[q] - 1)), sep = "")
    }

    if(sum(nx[-1]) > n_dep){
        emiss_cov_bar			<- lapply((q_emiss-1) * m * (nx[-1] - 1 ), dif_matrix, rows = J)
        names(emiss_cov_bar) <- dep_labels
        for(q in 1:n_dep){
            if(nx[1 + q] > 1){
                colnames(emiss_cov_bar[[q]]) <-  paste( paste("cov", rep(1 : (nx[1+q] - 1),each = nx[1+q]-1), "_", sep = ""), "emiss", rep(2:q_emiss[q], m * (nx[1 + q] - 1)), "_S", rep(1:m, each = (q_emiss[q] - 1) * (nx[1 + q] - 1)), sep = "")
                emiss_cov_bar[[q]][1,] <- 0
            } else {
                emiss_cov_bar[[q]] <- "No covariates where used to predict the emission probabilities for this outcome"
            }
        }
    } else{
        emiss_cov_bar <- "No covariates where used to predict the emission probabilities"
    }

    # Define object for subject specific posterior density (regression coefficients parameterization )
    if(subj_data){
        gamma_int_subj			<- rep(list(gamma_int_bar), n_subj)
        emiss_int_subj			<- rep(list(emiss_int_bar), n_subj)
    }

    # Put starting values in place for fist run forward algorithm
    emiss_sep 			<- vector("list", n_dep)
    for(q in 1:n_dep){
        start <- c(0, q_emiss * m)
        emiss_sep[[q]] <- matrix(PD[1,(sum(start[1:q]) + 1):(sum(start[1:q]) + (m * q_emiss[q]))], byrow = TRUE, ncol = q_emiss[q], nrow = m)
    }
    emiss				  <- rep(list(emiss_sep), n_subj)
    gamma 			<- rep(list(matrix(PD[1,(sum(m*q_emiss) + 1):(sum(m * q_emiss) + m * m)], byrow = TRUE, ncol = m)), n_subj)
    delta 			<- rep(list(solve(t(diag(m) - gamma[[1]] + 1), rep(1, m))), n_subj)


    # Start analysis --------------------------------------------
    # Run the MCMC algorithm
    itime <- proc.time()[3]
    if(show_progress == TRUE){
        cat("Progress of the Bayesian mHMM algorithm:", "\n")
        pb <- utils::txtProgressBar(min = 2, max = J, style = 3)
    }
    for (iter in 2 : J){

        # For each subject, obtain sampled state sequence with subject individual parameters ----------
        for(s in 1:n_subj){
            # Run forward algorithm, obtain subject specific forward proababilities and log likelihood
            forward				<- cat_mult_fw_r_to_cpp(x = subj_data[[s]]$y, m = m, emiss = emiss[[s]], gamma = gamma[[s]], n_dep = n_dep, delta=NULL)
            alpha         <- forward[[1]]
            c             <- max(forward[[2]][, subj_data[[s]]$n])
            llk           <- c + log(sum(exp(forward[[2]][, subj_data[[s]]$n] - c)))

            if(subj_data){
                PD_subj[[s]][iter, sum(m * q_emiss) + m * m + 1] <- llk
            }

            # Using the forward probabilites, sample the state sequence in a backward manner.
            # In addition, saves state transitions in trans, and conditional observations within states in cond_y
            trans[[s]]					                  <- vector("list", m)
            sample_path[[s]][n_vary[[s]], iter] 	<- sample(1:m, 1, prob = c(alpha[, n_vary[[s]]]))
            for(t in (subj_data[[s]]$n - 1):1){
                sample_path[[s]][t,iter] 	              <- sample(1:m, 1, prob = (alpha[, t] * gamma[[s]][,sample_path[[s]][t + 1, iter]]))
                trans[[s]][[sample_path[[s]][t,iter]]]	<- c(trans[[s]][[sample_path[[s]][t, iter]]], sample_path[[s]][t + 1, iter])
            }
            for (i in 1:m){
                trans[[s]][[i]] <- c(trans[[s]][[i]], 1:m)
                trans[[s]][[i]] <- rev(trans[[s]][[i]])
                for(q in 1:n_dep){
                    cond_y[[s]][[i]][[q]] <- c(subj_data[[s]]$y[sample_path[[s]][, iter] == i, q], 1:q_emiss[q])
                }
            }
        }

        # The remainder of the mcmc algorithm is state specific
        for(i in 1:m){

            # Obtain MLE of the covariance matrices and log likelihood of gamma and emiss at subject and population level -----------------
            # used to scale the propasal distribution of the RW Metropolis sampler

            # population level, transition matrix
            trans_pooled			  <- factor(c(unlist(sapply(trans, "[[", i)), c(1:m)))
            gamma_mle_pooled		<- optim(gamma_int_mle0, llmnl_int, Obs = trans_pooled,
                                       n_cat = m, method = "BFGS", hessian = T,
                                       control = list(fnscale = -1))
            gamma_int_mle_pooled[[i]]  <- gamma_mle_pooled$par
            gamma_pooled_ll[[i]]			<- gamma_mle_pooled$value

            # population level, conditional probabilities, seperate for each dependent variable
            for(q in 1:n_dep){
                cond_y_pooled					      <- numeric()
                ### MOET OOK ECHT BETER KUNNEN, eerst # cond_y_pooled				<- unlist(sapply(cond_y, "[[", m))
                for(s in 1:n_subj){
                    cond_y_pooled             <- c(cond_y_pooled, cond_y[[s]][[i]][[q]])
                }
                emiss_mle_pooled		<- optim(emiss_int_mle0[[q]], llmnl_int, Obs = c(cond_y_pooled, c(1:q_emiss[q])),
                                           n_cat = q_emiss[q], method = "BFGS", hessian = T,
                                           control = list(fnscale = -1))
                emiss_int_mle_pooled[[i]][[q]]  <- emiss_mle_pooled$par
                emiss_pooled_ll[[i]][[q]]				<- emiss_mle_pooled$value
            }

            # subject level
            for (s in 1:n_subj){
                wgt 				<- subj_data[[s]]$n / n_total

                # subject level, transition matrix
                gamma_out					<- optim(gamma_int_mle_pooled[[i]], llmnl_int_frac, Obs = c(trans[[s]][[i]], c(1:m)),
                                       n_cat = m, pooled_likel = gamma_pooled_ll[[i]], w = gamma_w, wgt = wgt,
                                       method="BFGS", hessian = TRUE, control = list(fnscale = -1))
                if(gamma_out$convergence == 0){
                    subj_data[[s]]$gamma_converge[i] <- 1
                    subj_data[[s]]$gamma_int_mle[i,] <- gamma_out$par
                    subj_data[[s]]$gamma_mhess[(1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1)), ]	<-
                        mnlHess_int(int = gamma_out$par, Obs = c(trans[[s]][[i]], c(1:m)), n_cat =  m)
                } else {
                    subj_data[[s]]$gamma_converge[i] <- 0
                    subj_data[[s]]$gamma_int_mle[i,] <- rep(0, m - 1)
                    subj_data[[s]]$gamma_mhess[(1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1)), ]	<- diag(m-1)
                }
                # if this is first iteration, use MLE for current values RW metropolis sampler
                if (iter == 2){
                    gamma_c_int[[i]][s,]		<- gamma_out$par
                }

                # subject level, conditional probabilities, seperate for each dependent variable
                for(q in 1:n_dep){
                    emiss_out				<- optim(emiss_int_mle_pooled[[i]][[q]], llmnl_int_frac, Obs = c(cond_y[[s]][[i]][[q]], c(1:q_emiss[q])),
                                          n_cat = q_emiss[q], pooled_likel = emiss_pooled_ll[[i]][[q]],
                                          w = emiss_w, wgt = wgt, method = "BFGS", hessian = TRUE, control = list(fnscale = -1))
                    if(emiss_out$convergence == 0){
                        subj_data[[s]]$emiss_converge[[q]][i]	 <- 1
                        subj_data[[s]]$emiss_int_mle[[q]][i,] <- emiss_out$par
                        subj_data[[s]]$emiss_mhess[[q]][(1 + (i - 1) * (q_emiss[q] - 1)):((q_emiss[q] - 1) + (i - 1) * (q_emiss[q] - 1)), ]		<-
                            mnlHess_int(int = subj_data[[s]]$emiss_int_mle[[q]][i,], Obs = c(cond_y[[s]][[i]][[q]], c(1:q_emiss[q])), n_cat =  q_emiss[q])
                    } else {
                        subj_data[[s]]$emiss_converge[[q]][i]	 <- 0
                        subj_data[[s]]$emiss_int_mle[[q]][i,]  <- rep(0, q_emiss[q] - 1)
                        subj_data[[s]]$emiss_mhess[[q]][(1 + (i - 1) * (q_emiss[q] - 1)):((q_emiss[q] - 1) + (i - 1) * (q_emiss[q] - 1)), ]	<- diag(q_emiss[q] - 1)
                    }
                    # if this is first iteration, use MLE for current values RW metropolis sampler
                    if (iter == 2){
                        emiss_c_int[[i]][[q]][s,]	<- emiss_out$par
                    }
                }
            }


            # Sample pouplaton values for gamma and conditional probabilities using Gibbs sampler -----------
            # gamma_mu0_n and gamma_mu_int_bar are matrices, with the number of rows equal to the number of covariates, and ncol equal to number of intercepts estimated
            gamma_mu0_n           <- solve(t(xx[[1]]) %*% xx[[1]] + gamma_K0)  %*% (t(xx[[1]]) %*% gamma_c_int[[i]] + gamma_K0 %*% gamma_mu0[[i]])
            gamma_V_n             <- gamma_V + t(gamma_c_int[[i]] - xx[[1]] %*% gamma_mu0_n) %*% (gamma_c_int[[i]] - xx[[1]] %*% gamma_mu0_n) + t(gamma_mu0_n - gamma_mu0[[i]]) %*% gamma_K0 %*% (gamma_mu0_n - gamma_mu0[[i]])
            gamma_V_int[[i]]      <- solve(rwish(S = solve(gamma_V_n), v = gamma_nu + n_subj))
            gamma_mu_int_bar[[i]] <- gamma_mu0_n + solve(chol(t(xx[[1]]) %*% xx[[1]] + gamma_K0)) %*% matrix(rnorm((m - 1) * nx[1]), nrow = nx[1]) %*% t(solve(chol(solve(gamma_V_int[[i]]))))
            gamma_exp_int				  <- matrix(exp(c(0, gamma_mu_int_bar[[i]][1,] )), nrow  = 1)
            gamma_mu_prob_bar[[i]] 	<- gamma_exp_int / as.vector(gamma_exp_int %*% c(rep(1,(m))))

            for(q in 1:n_dep){
                emiss_mu0_n                 <- solve(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]) %*% (t(xx[[1 + q]]) %*% emiss_c_int[[i]][[q]] + emiss_K0[[q]] %*% emiss_mu0[[q]][[i]])
                emiss_V_n                   <- emiss_V[[q]] + t(emiss_c_int[[i]][[q]] - xx[[1 + q]] %*% emiss_mu0_n) %*% (emiss_c_int[[i]][[q]] - xx[[1 + q]] %*% emiss_mu0_n) + t(emiss_mu0_n - emiss_mu0[[q]][[i]]) %*% emiss_K0[[q]] %*% (emiss_mu0_n - emiss_mu0[[q]][[i]])
                emiss_V_int[[i]][[q]]       <- solve(rwish(S = solve(emiss_V_n), v = emiss_nu[[q]] + n_subj))
                emiss_mu_int_bar[[i]][[q]]	 <- emiss_mu0_n + solve(chol(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]])) %*% matrix(rnorm((q_emiss[q] - 1) * nx[1 + q]), nrow = nx[1 + q]) %*% t(solve(chol(solve(emiss_V_int[[i]][[q]]))))
                emiss_exp_int				       <- matrix(exp(c(0, emiss_mu_int_bar[[i]][[q]][1, ])), nrow  = 1)
                emiss_mu_prob_bar[[i]][[q]] <- emiss_exp_int / as.vector(emiss_exp_int %*% c(rep(1, (q_emiss[q]))))
            }

            # Sample subject values for gamma and conditional probabilities using RW Metropolis sampler -----------
            for (s in 1:n_subj){
                gamma_candcov_comb 			<- chol2inv(chol(subj_data[[s]]$gamma_mhess[(1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1)), ] + chol2inv(chol(gamma_V_int[[i]]))))
                gamma_RWout					    <- mnl_RW_once(int1 = gamma_c_int[[i]][s,], Obs = trans[[s]][[i]], n_cat = m, mu_int_bar1 = c(t(gamma_mu_int_bar[[i]]) %*% xx[[1]][s,]), V_int1 = gamma_V_int[[i]], scalar = gamma_scalar, candcov1 = gamma_candcov_comb)
                gamma[[s]][i,]  	<- gamma_RWout$prob
                if(subj_data){
                    PD_subj[[s]][iter, c((sum(m * q_emiss) + 1 + (i - 1) * m):(sum(m * q_emiss) + (i - 1) * m + m))] <- gamma_RWout$prob
                }
                gamma_naccept[s, i]			<- gamma_naccept[s, i] + gamma_RWout$accept
                gamma_c_int[[i]][s,]		<- gamma_RWout$draw_int

                if(subj_data){
                    gamma_int_subj[[s]][iter, (1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1))] <- gamma_c_int[[i]][s,]
                }

                start <- c(0, q_emiss * m)
                for(q in 1:n_dep){
                    emiss_candcov_comb		     <- chol2inv(chol(subj_data[[s]]$emiss_mhess[[q]][(1 + (i - 1) * (q_emiss[q] - 1)):((q_emiss[q] - 1) + (i - 1) * (q_emiss[q] - 1)), ] + chol2inv(chol(emiss_V_int[[i]][[q]]))))
                    emiss_RWout				       <- mnl_RW_once(int1 = emiss_c_int[[i]][[q]][s,], Obs = cond_y[[s]][[i]][[q]], n_cat = q_emiss[q], mu_int_bar1 = c(t(emiss_mu_int_bar[[i]][[q]]) %*% xx[[1 + q]][s,]), V_int1 = emiss_V_int[[i]][[q]], scalar = emiss_scalar[[q]], candcov1 = emiss_candcov_comb)
                    emiss[[s]][[q]][i,]		   <- emiss_RWout$prob
                    if(subj_data){
                        PD_subj[[s]][iter, (sum(start[1:q]) + 1 + (i - 1) * q_emiss[q]):(sum(start[1:q]) + (i - 1) * q_emiss[q] + q_emiss[q])] <- emiss_RWout$prob
                    }
                    emiss_naccept[[q]][s, i]	 <- emiss_naccept[[q]][s, i] + emiss_RWout$accept
                    emiss_c_int[[i]][[q]][s,] <- emiss_RWout$draw_int

                    if(subj_data){
                        emiss_int_subj[[s]][[q]][iter, (1 + (i - 1) * (q_emiss[q] - 1)) : ((q_emiss[q] - 1) + (i - 1) * (q_emiss[q] - 1))]	<- emiss_c_int[[i]][[q]][s,]
                    }
                }

                if(i == m){
                    delta[[s]] 		<- solve(t(diag(m) - gamma[[s]] + 1), rep(1, m))
                }
            }
        }


        # End of 1 MCMC iteration, save output values --------
        gamma_int_bar[iter, ]				   	<- unlist(lapply(gamma_mu_int_bar, "[",1,))
        gamma_V_int_bar[iter, ] <- unlist(lapply(gamma_V_int, function(e) as.vector(t(e))))
        if(nx[1] > 1){
            gamma_cov_bar[iter, ]      	<- unlist(lapply(gamma_mu_int_bar, "[",-1,))
        }
        gamma_prob_bar[iter,]			<- unlist(gamma_mu_prob_bar)
        for(q in 1:n_dep){
            emiss_int_bar[[q]][iter, ]	<- as.vector(unlist(lapply(
                lapply(emiss_mu_int_bar, "[[", q), "[",1,)
            ))
            if(nx[1+q] > 1){
                emiss_cov_bar[[q]][iter, ]  <- as.vector(unlist(lapply(
                    lapply(emiss_mu_int_bar, "[[", q), "[",-1,)
                ))
            }
            emiss_prob_bar[[q]][iter,]	<- as.vector(unlist(sapply(emiss_mu_prob_bar, "[[", q)))
            emiss_V_int_bar[[q]][iter,] <- unlist(lapply(emiss_V_int, function(e) as.vector(t(e[[q]]))))
        }
        if(show_progress == TRUE){
            utils::setTxtProgressBar(pb, iter)
        }
    }
    if(show_progress == TRUE){
        close(pb)
    }

    # End of function, return output values --------
    ctime = proc.time()[3]
    # message(paste("Total time elapsed (hh:mm:ss):", hms(ctime-itime)))

    out <- list(input = list(m = m, n_dep = n_dep, q_emiss = q_emiss, J = J,
                             burn_in = burn_in, n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                gamma_int_bar = gamma_int_bar,
                gamma_V_int_bar = gamma_V_int_bar,
                gamma_cov_bar = gamma_cov_bar,
                emiss_int_bar = emiss_int_bar,
                emiss_V_int_bar = emiss_V_int_bar,
                emiss_cov_bar = emiss_cov_bar,
                gamma_prob_bar = gamma_prob_bar, emiss_prob_bar = emiss_prob_bar,
                gamma_naccept = gamma_naccept, emiss_naccept = emiss_naccept)

    if(return_path == TRUE){
        out <- c(out,list(sample_path = sample_path))
    }

    if(subj_data){
        out <- c(out, list(PD_subj = PD_subj, gamma_int_subj = gamma_int_subj, emiss_int_subj = emiss_int_subj))
    }

    class(out) <- append(class(out), "mHMM")
    return(out)
    }
