// Estimates more similar parameters to LB-SPR
// with more similar inputs
//
// Requires:
//   1) estimates of VB growth parameters, including CV on Linf
//   2) maturity at size
//   3) length composition of catch
// 
//  Estimates:
//   1) selectivity at length
//   2) annual F as fixed effect
//   3) annual R as random effect 

// #define TMB_LIB_INIT R_init_LIME
#include <TMB.hpp>
#include <iostream>

// dlnorm
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

// transformation to ensure correlation parameters are between -1 and 1
template <class Type>
Type rho_trans(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // data
  DATA_MATRIX(dat_growth); // 2 columns: age, length
  DATA_INTEGER(n_islands); // number of islands = 3
  DATA_IVECTOR(island); // 1, 2, or 3

  // parameters
  PARAMETER_VECTOR(mu); // mean vonB pars: Linf, logK, t0
  PARAMETER_VECTOR(log_sigma_mu); // variance of mu
  PARAMETER_VECTOR(rho_untrans); // untransformed correlation of mu (21, 31, 32)
  PARAMETER(log_CV_L); // coefficient of variation in the growth curve (shared all islands)
  PARAMETER_MATRIX(emat); // island-specific growth par deviations from mean, random effects

  using namespace density;
  Type jnll = 0;

  // transformations
  vector<Type> sigma_mu = exp(log_sigma_mu.array()).matrix(); // ensures sds are positive
  vector<Type> rho(3);
  for(int i=0; i<3; i++) rho(i) = rho_trans(rho_untrans(i)); // ensures correlations are between -1 and 1
  Type CV_L = exp(log_CV_L); // ensures CV is positive

  vector<Type> theta_global(3);
  theta_global(0) = mu(0);
  theta_global(1) = exp(mu(1));
  theta_global(2) = mu(2);

  // construct Sigma
  matrix<Type> Sigma(3,3);
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      Sigma(i,j) = sigma_mu(i) * sigma_mu(j);
    }
  }
  Sigma(1,0) *= rho(0);
  Sigma(2,0) *= rho(1);
  Sigma(2,1) *= rho(2);
  Sigma(0,1) *= rho(0);
  Sigma(0,2) *= rho(1);
  Sigma(1,2) *= rho(2);
  MVNORM_t<Type> your_dmnorm(Sigma);

  // likelihood of island-specific growth par deviations from mean
  matrix<Type> theta_untrans(3,n_islands); // island-specific growth pars
  for(int j=0; j<n_islands; j++){
    jnll += your_dmnorm(emat.col(j));
    for(int i=0; i<3; i++) theta_untrans(i,j) = mu(i) + emat(i,j);
  }

  // island-specific growth pars, transformed to Linf, vbk, t0
  matrix<Type> theta(3,n_islands);
  for(int j=0; j<n_islands; j++){
    theta(0,j) = theta_untrans(0,j);
    theta(1,j) = exp(theta_untrans(1,j));
    theta(2,j) = theta_untrans(2,j);
  }

  // likelihood of growth data, using island-specific growth pars
  int n_g = dat_growth.rows();
  vector<Type> pred_len(n_g);
  for(int g=0; g<n_g; g++){
    pred_len(g) = theta(0,island(g)-1)*(1-exp(-theta(1,island(g)-1)*(dat_growth(g,0)-theta(2,island(g)-1))));
    jnll -= dnorm(dat_growth(g,1), pred_len(g), pred_len(g)*CV_L, true);
  }

  ADREPORT(theta);
  ADREPORT(theta_global);
  ADREPORT(CV_L);

  REPORT(Sigma);
  REPORT(mu);
  REPORT(emat);
  REPORT(log_sigma_mu);
  REPORT(pred_len);
  REPORT(jnll); 

  return(jnll);
}
