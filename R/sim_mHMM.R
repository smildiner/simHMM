#' @export

# Function created by Emmeke Aarts: https://github.com/emmekeaarts/mHMMbayes/tree/develop

# write test for multivariate categorical data
# write tests for continuous data
# write stops for correctly specifying things like n_dep, data_distr, beta for multivariate data
# write stop to check that number of beta's specified corresponds to number of covariates
#   specified in xx (besides checking that both of them are present)

sim_mHMM <- function(n_t, n, data_distr = 'categorical', m, n_dep = 1,
                     start_state = NULL, q_emiss = NULL, gamma, emiss_distr, xx_vec = NULL, beta = NULL,
                     var_gamma = 0.1, var_emiss = NULL, return_ind_par = FALSE){

    #############
    # Inbuild checks for correct specification of parameters ---------------------
    #############

    if (dim(gamma)[1] != m){
        stop(paste("The transiton probability matrix gamma should be a", m, "by", m, "matrix."))
    }
    if (dim(gamma)[2] != m){
        stop(paste("The transiton probability matrix gamma should be a", m, "by", m, "matrix."))
    }
    if(!isTRUE(all.equal(apply(gamma,1,sum), rep(1,m)))){
        stop("The elements in each row of the transition probability matrix gamma should sum up to 1")
    }
    if(!is.list(emiss_distr)){
        stop("The format of emiss_distr should be a list with", n_dep, "elements.")
    }
    if(length(emiss_distr) != n_dep){
        stop("The number of dependent variables specified in n_dep and the number of elements specified in the list emiss_distr should be equal")
    }
    if(data_distr == "categorical" & length(q_emiss) != n_dep){
        stop("The lenght of q_emiss specifying the number of output categories for each of the number of dependent variables should equal the number of dependent variables specified in n_dep")
    }
    for(i in 1:n_dep){
        if (dim(emiss_distr[[i]])[1] != m){
            stop(paste("The number of rows of emission distribution matrix in element", i, "should be
                       equal to the number of states, which is", m, "."))
        }
        if(data_distr == 'categorical'){
            if (dim(emiss_distr[[i]])[2] != q_emiss[i]){
                stop(paste("The number of columns of the emission distribution matrix should be
                           equal to the number of observable categories, which is", q_emiss[i], ". See emission distribution in element", i, "."))
            }
            if(!isTRUE(all.equal(apply(emiss_distr[[i]], 1, sum), rep(1, m)))){
                stop("The elements in each row of the emission distribution matrix should sum up to 1, see emission distribution in element", i, ".")
            }
        }
        # if(data_distr == 'continuous'){
        #
        # }
    }

    if((is.null(xx_vec) & !is.null(beta)) | (!is.null(xx_vec) & is.null(beta))){
        stop("Either only xx_vec or only beta is specified. Please specify both 1) the values for the covariate
             in xx_vec and 2) the values of the regression parameters in beta, to allow correct simulation of the
             data.")
    }
    if(!is.null(xx_vec)){
        if((!is.null(xx_vec[[1]]) & is.null(beta[[1]])) |
           (!is.null(xx_vec[[2]]) & is.null(beta[[2]]))){
            stop("Either only xx_vec or only beta is specified in one of the elements.
                 Please specify both 1) the values for the covariate in xx_vec and 2)
                 the values of the regression parameters in beta if either one is not
                 empty, to allow correct simulation of the data.")
        }
        }
    if(!is.null(beta)){
        # extend to all 1 + n_dep
        if((!is.null(beta[[1]]) & is.null(xx_vec[[1]])) |
           (!is.null(beta[[2]]) & is.null(xx_vec[[2]]))){
            stop("Either only xx_vec or only beta is specified in one of the elements.
                 Please specify both 1) the values for the covariate in xx_vec and 2)
                 the values of the regression parameters in beta if either one is not
                 empty, to allow correct simulation of the data.")
        }
        }
    if(!is.null(xx_vec)){
        # extend to all 1 + n_dep
        if((!is.null(xx_vec[[1]]) & length(xx_vec[[1]]) != n) |
           (!is.null(xx_vec[[2]]) & length(xx_vec[[2]]) != n)){
            stop("The length of the vectors in xx_vec should be equal to the number of subjects to be simulated,
                 set in n, if (the element in) xx_vec is not set to NULL.")
        }
        }
    if (!is.null(beta)){
        if (!is.null(beta[[1]])){
            if ((dim(beta[[1]])[1] != (m)) | (dim(beta[[1]])[2] != (m-1))){
                stop(paste("The first element of beta to predict the transiton probability matrix gamma should be a m (", m, " ) by m - 1 (", m - 1, ") matrix."))
            }
        }
        if (!is.null(beta[[2]]) & data_distr == 'categorical'){
            # extend to all 1 + n_dep and continuous
            if((dim(beta[[2]])[1] != (m)) | (dim(beta[[2]])[2] != (q_emiss[1]-1))){
                stop(paste("The second element of beta to predict the emission distribution should be a m (", m, ") by q_emiss - 1 (", q_emiss[1] - 1, ") matrix."))
            }
        }
    }
    if(is.null(xx_vec)){
        xx_vec <- rep(list(NULL), n_dep + 1)
        for(i in 1:(n_dep + 1)){
            xx_vec[[i]] <- rep(1,n)
        }
    } else {
        for(i in 1:(n_dep + 1)){
            if(is.null(xx_vec[[i]])) {
                xx_vec[[i]] <- rep(1,n)
            }
        }
    }
    if(is.null(beta)){
        beta <- rep(list(NULL), n_dep + 1)
        beta[[1]] <- matrix(0, ncol = m - 1, nrow = m)
        for(i in 2:(n_dep + 1)){
            if(data_distr == 'categorical'){
                beta[[i]] <- matrix(0, ncol = q_emiss[i-1] - 1, nrow = m)
            } else if (data_distr == 'continuous'){
                beta[[i]] <- matrix(0, ncol = 1, nrow = m)
            }
        }
    } else {
        if(is.null(beta[[1]])) {
            beta[[1]] <- matrix(0, ncol = m - 1, nrow = m)
        }
        for (i in 2:(n_dep + 1)){
            if (is.null(beta[[i]])) {
                if(data_distr == 'categorical'){
                    beta[[i]] <- matrix(0, ncol = q_emiss[i-1] - 1, nrow = m)
                } else if (data_distr == 'continuous'){
                    beta[[i]] <- matrix(0, ncol = 1, nrow = m)
                }
            }
        }
    }

    if(data_distr == 'continuous'){
        if(n == 1){
            var_gamma <- 0
            var_emiss <- rep(0, n_dep)
        }

        if(is.null(var_emiss)){
            var_emiss <- rep(0.1, n_dep)
        } else if(length(var_emiss) != n_dep){
            stop("The lenght of var_emiss specifying variance between subjects in each of the number of dependent variables should equal the number of dependent variables specified in n_dep")
        }
    }

    # UPDATED CHECKS: else if(length(var_emiss)==){old check}, if(is.list(var_emiss)){new check}
    if(data_distr == 'categorical'){

        # If only 1 subject
        if(n == 1){
            var_gamma <- matrix(rep(0, m*(m-1)),nrow = m, byrow = TRUE)
            var_emiss <- rep(list(NULL), n_dep)
            for(i in 1:n_dep){
                var_emiss[[i]] <- matrix(rep(0, m*(q_emiss[i]-1)),nrow = m, byrow = TRUE)
            }
        }

        # If a single value of var_gamma specified, use for all categories
        if(length(var_gamma) == 1){
            var_gamma <- matrix(rep(var_gamma, m*(m-1)),nrow = m, byrow = TRUE)
        } else if(is.matrix(var_gamma)){
            if (dim(var_gamma)[1] != m){
                stop(paste("The betweem-subject variance matrix for the transition distribution should be a", m, "by", m-1, "matrix."))
            }
            if (dim(var_gamma)[2] != m-1){
                stop(paste("The betweem-subject variance matrix for the transition distribution should be a", m, "by", m-1, "matrix."))
            }
        }

        # If a single value of var_emiss specified, use for all categories and n_dep
        if(is.null(var_emiss)){
            var_emiss <- rep(list(NULL), n_dep)
            for(i in 1:n_dep){
                var_emiss[[i]] <- matrix(rep(0.1, m*(q_emiss[i]-1)),nrow = m, byrow = TRUE)
            }
        } else if(is.numeric(var_emiss) & length(var_emiss) == n_dep){
            arg_var_emiss <- var_emiss
            var_emiss <- rep(list(NULL), n_dep)
            for(i in 1:n_dep){
                var_emiss[[i]] <- matrix(rep(arg_var_emiss[i], m*(q_emiss[i]-1)),nrow = m, byrow = TRUE)
            }
        } else if(is.list(var_emiss)){
            for(i in 1:n_dep){
                if(dim(var_emiss[[i]])[2] != q_emiss[i]-1){
                    stop(paste("The number of columns of the between-subject variance for the emission distribution should be
                           equal to the number of observable categories minus one, which is", q_emiss[i], ". See emission distribution in element", i, "."))
                }
            }
        } else if(length(var_emiss) != n_dep){
            stop("The lenght of var_emiss specifying variance between subjects in each of the number of dependent variables should equal the number of dependent variables specified in n_dep. Note that var_emiss can either by a list of matrices, or a numeric vector.")
        }
    }

    #############
    # Simulating the data ---------------------
    #############

    states <- matrix(ncol = 2, nrow = n_t*n)
    states[,1] <- rep(1:n, each = n_t)
    obs <- matrix(ncol = 1 + n_dep, nrow = n_t*n)
    obs[,1] <- rep(1:n, each = n_t)
    sub_gamma <- rep(list(NULL), n)
    sub_emiss <- rep(list(vector("list", n_dep)), n)
    mnl_gamma <- prob_to_int(gamma)
    if(data_distr == "categorical"){
        mnl_emiss <- rep(list(NULL), n_dep)
        for(i in 1:n_dep){
            mnl_emiss[[i]] <- prob_to_int(emiss_distr[[i]])
        }
    }
    for(j in 1:n){
        # sub_gamma[[j]] <- int_to_prob(mnl_gamma + xx_vec[[1]][j] * beta[[1]] +
        #                                   rnorm(n = m * (m-1), mean = 0, sd = sqrt(var_gamma)))
        sub_gamma[[j]] <- int_to_prob(mnl_gamma + xx_vec[[1]][j] * beta[[1]] +
                                          rnorm(n = m * (m-1), mean = 0, sd = sqrt(as.numeric(var_gamma))))
        for(i in 1:n_dep){
            if(data_distr == "categorical"){
                # sub_emiss[[j]][[i]] <- int_to_prob(mnl_emiss[[i]] + xx_vec[[1+i]][j] * beta[[1+i]] +
                #                                        rnorm(n = m * (q_emiss[i]-1), mean = 0, sd = sqrt(var_emiss[i])))
                sub_emiss[[j]][[i]] <- int_to_prob(mnl_emiss[[i]] + xx_vec[[1+i]][j] * beta[[1+i]] +
                                                       rnorm(n = m * (q_emiss[i]-1), mean = 0, sd = sqrt(as.numeric(var_emiss[[i]]))))
            } else if(data_distr == "continuous"){
                sub_emiss[[j]][[i]] <- emiss_distr[[i]]
                sub_emiss[[j]][[i]][,1] <- emiss_distr[[i]][,1] +  xx_vec[[1+i]][j] * beta[[1+i]] +
                    rnorm(n = m, mean = 0, sd = sqrt(var_emiss[i]))
            }
        }

        if(n_t != 0){
            init <- solve(t(diag(m) - sub_gamma[[j]] + 1), rep(1, m))
            if (is.null(start_state)){
                states[((j-1) * n_t + 1), 2] <- sample(x = 1:m, size = 1, prob = init)
            } else {
                states[((j-1) * n_t + 1), 2] <- start_state
            }
            if(data_distr == "categorical"){
                for(i in 1:n_dep){
                    obs[((j-1) * n_t + 1), (1+i)] <- sample(x = 1:q_emiss[i], size = 1, prob = sub_emiss[[j]][[i]][states[((j-1) * n_t + 1), 2],])
                }
            } else if (data_distr == "continuous"){
                for(i in 1:n_dep){
                    obs[((j-1) * n_t + 1), (1+i)] <- rnorm(1, mean = sub_emiss[[j]][[i]][states[((j-1) * n_t + 1), 2],1], sd = sqrt(sub_emiss[[j]][[i]][states[((j-1) * n_t + 1), 2],2]))
                }
            }
            for(t in 2:n_t){
                states[((j-1) * n_t + t), 2] <- sample(x = 1:m, size = 1, prob = sub_gamma[[j]][states[((j-1) * n_t + t - 1), 2],])
                if(data_distr == "categorical"){
                    for(i in 1:n_dep){
                        obs[((j-1) * n_t + t), (1+i)] <- sample(x = 1:q_emiss[i], size = 1, prob = sub_emiss[[j]][[i]][states[((j-1) * n_t + t), 2],])
                    }
                } else if (data_distr == "continuous"){
                    for(i in 1:n_dep){
                        obs[((j-1) * n_t + t), (1+i)] <- rnorm(1, mean = sub_emiss[[j]][[i]][states[((j-1) * n_t + t), 2],1], sd = sqrt(sub_emiss[[j]][[i]][states[((j-1) * n_t + t), 2],2]))
                    }
                }
            }
        }
    }

    #############
    # Returning output  ---------------------
    #############
    colnames(states) <- c("subj", "state")
    colnames(obs)    <- c("subj", paste("observation", 1:n_dep))
    if (return_ind_par == FALSE & n_t != 0){
        return(list(states = states, obs = obs))
    } else if (return_ind_par == TRUE & n_t != 0){
        return(list(states = states, obs = obs, subject_gamma = sub_gamma, subject_emiss = sub_emiss))
    } else if (n_t == 0){
        return(list(subject_gamma = sub_gamma, subject_emiss = sub_emiss))
    }
    }
