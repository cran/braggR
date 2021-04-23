#include <Rcpp.h>
#include <random>
#include <math.h>
using namespace std;
using namespace Rcpp;
static const double pi = 3.141592653589793238462643383280;

bool valid_parameters(double delta, double gamma, double rho, double N) {
  bool check_bool = ( ((N*gamma*gamma/(N-1) - gamma/(N-1)) < rho) && (rho < gamma) && (gamma <= delta));
  return check_bool;
}

double logposterior(NumericVector P,
                    double rho, double gamma, double delta, double mu0,
                    double sum_P, double sum_P2, double sum_P_outer, double sum_P_outer_plus,
                    double N, double a, double b){
  double ret_log = -std::numeric_limits<double>::infinity();

  if(valid_parameters(delta, gamma, rho, N)){
    double mu = mu0/sqrt(1-gamma);
    double S3 = sum_P2 - 2*mu*sum_P + N*mu*mu;
    double S4 = S3 + 2*sum_P_outer - 2*mu*sum_P_outer_plus + N*(N-1)*mu*mu;
    double eta1 = (delta + (N-1)*rho)/(1-gamma);
    double eta2 = (delta-rho)/(1-gamma);
    ret_log = -N/2*log(2*pi) - 0.5*log(eta1) - (N-1)/2*log(eta2) - (N*S3-S4)/(2*N*eta2) - S4/(2*N*eta1); // Likelihood
    ret_log = ret_log - (0.5+a)*log(eta1) - (0.5+b)*log(eta2) - (2.5+a)*log(1-gamma); // Prior
  }
  return ret_log;
}

// [[Rcpp::export(.sample_parameters_cpp)]]
NumericMatrix sample_parameters_cpp(NumericVector p,
                                    double p0 = 0.5,
                                    double alpha = -1, double beta = -1,
                                    double a = 0.5, double b = 0.5,
                                    int num_sample = 100000,
                                    int seed = 1){
  int N = p.size();
  bool beta_distr = ((alpha>0) && (beta>0));

  NumericVector P = qnorm(p, 0.0, 1.0);
  double mu0;
  double sum_P = sum(P);
  double sum_P2 = sum(P*P);
  double sum_P_outer = 0;
  double sum_P_outer_plus = 0;
  for(int iter1 = 0; iter1 < N; iter1++){
    for(int iter2 = 0; iter2 < iter1; iter2++){
      sum_P_outer += P[iter1] * P[iter2];
      sum_P_outer_plus += P[iter1] + P[iter2];
    }
  }

  default_random_engine generator;
  generator.seed(seed);
  exponential_distribution<double> exp_RV(1.0);
  uniform_real_distribution<double> uni_RV(0.0,1.0);
  normal_distribution<double> normal_RV(0.0,1.0);
  gamma_distribution<double> gamma_RV_alpha(alpha,1.0);
  gamma_distribution<double> gamma_RV_beta(beta,1.0);

  NumericMatrix samples(num_sample, 5);
  NumericMatrix::Column p_post = samples(_,0);
  NumericMatrix::Column rhos = samples(_,1);
  NumericMatrix::Column gammas = samples(_,2);
  NumericMatrix::Column deltas = samples(_,3);
  NumericMatrix::Column p0s = samples(_,4);

  // Starting points:
  if(beta_distr){
    p0s[0] = 0.5;
  } else {
    p0s[0] = p0;
    mu0 = R::qnorm5(p0, 0.0, 1.0, 1, 0);
  }
  rhos[0] = 0;
  gammas[0] = 1/((double)N+1);
  deltas[0] = 2;

  double width = 10; // Should not be too small
  double L_delta, H_delta, L_gamma, H_gamma, L_rho, H_rho, u_log_lik, param_log_lik, y;
  double rho_cand, delta_cand, gamma_cand;
  double mu_cond, delta_cond;
  double X_cond, Y_cond;

  for(int iter = 1;  iter < num_sample; iter++){

    if(beta_distr){
      // Sample p0:
      X_cond = gamma_RV_alpha(generator);
      Y_cond = gamma_RV_beta(generator);
      p0s[iter] = X_cond / (X_cond + Y_cond);
      mu0 = R::qnorm5(p0s[iter], 0.0, 1.0, 1, 0);
    } else {
      p0s[iter] = p0;
    }

    // Sample rho, gamma, delta
    u_log_lik = logposterior(P, rhos[iter-1], gammas[iter-1], deltas[iter-1],
                             mu0, sum_P, sum_P2, sum_P_outer, sum_P_outer_plus, N,
                             a, b);
    y = u_log_lik - exp_RV(generator);

    L_delta = deltas[iter-1] - uni_RV(generator)*width;
    H_delta = L_delta + width;
    L_gamma = 0.0;
    H_gamma = 1.0;
    L_rho = 0.0;
    H_rho = 1.0;

    param_log_lik = -std::numeric_limits<double>::infinity();
    while(param_log_lik < y){
      rho_cand = L_rho + uni_RV(generator)*(H_rho-L_rho);           //runif(1, L_rho, H_rho)
      gamma_cand = L_gamma + uni_RV(generator)*(H_gamma-L_gamma);   //runif(1, L_gamma, H_gamma)
      delta_cand =  L_delta + uni_RV(generator)*(H_delta-L_delta);  //runif(1, L_delta, H_delta)

      param_log_lik = logposterior(P, rho_cand, gamma_cand, delta_cand,
                                   mu0, sum_P, sum_P2, sum_P_outer, sum_P_outer_plus, N,
                                   a, b);

      // Fix dimension if needed:
      if(rho_cand < rhos[iter-1]){
	        L_rho = rho_cand;
      } else {
	        H_rho = rho_cand;
      }

      if(delta_cand < deltas[iter-1]){
	        L_delta = delta_cand;
      } else {
	        H_delta = delta_cand;
      }

      if(gamma_cand < gammas[iter-1]){
	        L_gamma = gamma_cand;
      } else {
	        H_gamma = gamma_cand;
      }
    }

    rhos[iter] = rho_cand;
    deltas[iter] = delta_cand;
    gammas[iter] = gamma_cand;

    mu_cond = mu0 +  sqrt(1-gamma_cand)*gamma_cand / (delta_cand + (N-1)*rho_cand) * (sum_P - N*mu0 / sqrt(1-gamma_cand));
    delta_cond = 1 - gamma_cand*gamma_cand*N / (delta_cand + (N-1)*rho_cand);
    p_post[iter] = 1-R::pnorm(0, mu_cond, sqrt(delta_cond), 1, 0);
  }

  return samples;
}
