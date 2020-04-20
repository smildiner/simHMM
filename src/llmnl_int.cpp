#include <Rcpp.h>
using namespace Rcpp;

// Evaluate loglikelihood for intercept only MNL, bassed on P. Rossi 2004
//' @keywords internal
// [[Rcpp::export(rng = false)]]
double llmnl_int(NumericVector beta, IntegerVector Obs, int n_cat) {

    int n_Obs = Obs.size();
    //2: Calculate log sum only once:
    // double expBetas_log_sum = log(sum(exp(betas)));
    double expBetas_log_sum = 1.0; // std::exp(0)
    for (int i = 1; i < n_cat; i++) {
        expBetas_log_sum += std::exp(beta[i-1]);
    };
    expBetas_log_sum = std::log(expBetas_log_sum);

    double ll_sum = 0;
    //3: Use n_Obs, to avoid calling Xby.size() every time
    for (int i = 0; i < n_Obs; i++) {
        if(Obs[i] == 1L) continue;
        ll_sum += beta[Obs[i]-2L];
    };
    //4: Use that we know denom is the same for all I:
    ll_sum = ll_sum - expBetas_log_sum * n_Obs;
    return ll_sum;
}
