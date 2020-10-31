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
  // ========== Inputs ============================

    // Indices
    DATA_INTEGER(n_t); // total number of years
    DATA_INTEGER(n_lb); // number of length bins
    DATA_INTEGER(n_fl); // number of fleets
    DATA_INTEGER(n_a); // number of ages

    // Data in likelihood
    DATA_ARRAY(LF_tlf); // length composition
    DATA_MATRIX(n_LF_ft); // number of independent observations by fleet annually, between 1 and C_t
    DATA_MATRIX(I_ft); // CPUE for each year
    // DATA_INTEGER(I_type); // if I_type=0, no index data, if I_type=1, numbers, if I_type=2, biomass    
    DATA_MATRIX(C_ft); // catch each year
    DATA_INTEGER(C_type); // if C_type=0, no catch data, if C_type=1, numbers, if C_type=2, biomass
    DATA_VECTOR(ML_t); // mean length each year
    DATA_IVECTOR(ML_years); // for each year, 1 = have length data, 0 = no data

    // Known values
    DATA_VECTOR(ages); // ages
    DATA_VECTOR(match_ages); // annual ages
    // DATA_VECTOR(L_a); // length-at-age
    // DATA_VECTOR(W_a); // weight-at-age
    DATA_VECTOR(lw_pars); // length-weight params (a, b)
    DATA_VECTOR(Mat_l); // maturity at length
    DATA_SCALAR(M); // natural mortality
    DATA_SCALAR(h); // steepness
    // DATA_VECTOR(Mat_a); // maturity-at-age
    DATA_VECTOR(lbhighs); // upper length bins
    DATA_VECTOR(lbmids); // mid length bins

    // penalties
    DATA_INTEGER(Fpen); // penalty on annual fishing mortality
    DATA_INTEGER(SigRpen); // penalty on sigmaR
    DATA_VECTOR(SigRprior); //mean and standard deviation for sigmaR prior

    // option for fixed time series for selectivity
    DATA_IVECTOR(selex_type_f); //0 =fixed, 1=estimated
    DATA_MATRIX(vals_selex_ft);

    // option for likelihood distribution for length comps
    DATA_INTEGER(LFdist); // 0=multinomial, 1=dirichlet-multinomial

    // option for shorter time-step than years
    DATA_IVECTOR(S_yrs); // matching each time step with a year
    DATA_INTEGER(n_s); // number of time steps within a year
    DATA_INTEGER(n_y); // number of years

    //mirror options
    DATA_INTEGER(mirror_theta);
    DATA_INTEGER(mirror_q);

    // estimate total F instead of by fleet
    DATA_INTEGER(est_totalF);
    DATA_VECTOR(prop_f);

    // new growth data
    DATA_MATRIX(dat_growth); // 2 columns: age, length
    DATA_INTEGER(n_islands); // number of islands = 3
    DATA_IVECTOR(island); // 1, 2, or 3

  // ======== Parameters =================================
    // Fixed, estimated parameters
    PARAMETER_MATRIX(log_F_ft);  // fishing mortality by fleet
    PARAMETER_VECTOR(log_q_f); // catachability by fleet associated with index
    PARAMETER(beta); // equilibrium recruitment
    PARAMETER(log_sigma_R); // recruitment variation
    PARAMETER_VECTOR(log_S50_f); // Length at 50% selectivity
    PARAMETER_VECTOR(log_Sdelta_f); // log(L95 - L50)
    PARAMETER(log_sigma_F) // fishing mortality standard deviation
    PARAMETER(log_sigma_C); // observation error - catch
    PARAMETER(log_sigma_I); // observation error - index
    PARAMETER(log_CV_L); // coefficient of variation in the age-length curve
    PARAMETER_VECTOR(log_theta); // dirichlet-multinomial parameter

    // estimate growth inside model
    PARAMETER_VECTOR(growth_mu); // mean vonB pars: Linf, logK, t0
    PARAMETER_VECTOR(log_sigma_mu); // variance of mu
    PARAMETER_VECTOR(rho_untrans); // untransformed correlation of mu (21, 31, 32)
    PARAMETER_MATRIX(emat); // island-specific growth par deviations from mean, random effects

    // Random effects
    PARAMETER_VECTOR(Nu_input); // temporal variation in recruitment


  // ============ Global values ============================

  using namespace density;
  Type jnll=0;
  vector<Type> jnll_comp(8);
  jnll_comp.setZero();

  // ======= Transform parameters =========================

  // variation terms
  Type sigma_F = exp(log_sigma_F);
  Type sigma_R = exp(log_sigma_R);
  Type sigma_C = exp(log_sigma_C);
  Type sigma_I = exp(log_sigma_I);
  Type CV_L = exp(log_CV_L);

  // transformations
  vector<Type> sigma_mu = exp(log_sigma_mu.array()).matrix(); // ensures sds are positive
  vector<Type> rho(3);
  for(int i=0; i<3; i++) rho(i) = rho_trans(rho_untrans(i)); // ensures correlations are between -1 and 1
  Type linf = growth_mu(0);
  Type vbk = exp(growth_mu(1));
  Type t0 = growth_mu(2);

  // selectivity
  vector<Type> S50_f(n_fl);
  //vector<Type> Sdelta_f(n_fl);
  vector<Type> S95_f(n_fl);
  for(int f=0;f<n_fl;f++){
    S50_f(f) = exp(log_S50_f(f));
    S95_f(f) = S50_f(f) + exp(log_Sdelta_f(f));
  }

  // dirichlet-multinomial parameter
  vector<Type> theta(n_fl);
  for(int f=0;f<n_fl;f++){
    if(mirror_theta!=1) theta(f) = exp(log_theta(f));
    if(mirror_theta==1) theta(f) = exp(log_theta(0));
  }

  // catchability coefficient
  vector<Type> q_f(n_fl);
  for(int f=0;f<n_fl;f++){
    if(mirror_q!=1) q_f(f) = exp(log_q_f(f));
    if(mirror_q==1) q_f(f) = exp(log_q_f(0));
  } 

  // ========= Probability of being length at age  =============================
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
    jnll_comp(7) += your_dmnorm(emat.col(j));
    for(int i=0; i<3; i++) theta_untrans(i,j) = growth_mu(i) + emat(i,j);
  }

  // island-specific growth pars, transformed to Linf, vbk, t0
  matrix<Type> growthpars(3,n_islands);
  for(int j=0; j<n_islands; j++){
    growthpars(0,j) = theta_untrans(0,j);
    growthpars(1,j) = exp(theta_untrans(1,j));
    growthpars(2,j) = theta_untrans(2,j);
  }

  // likelihood of growth data, using island-specific growth pars
  int n_g = dat_growth.rows();
  vector<Type> pred_len(n_g);
  for(int g=0; g<n_g; g++){
    pred_len(g) = growthpars(0,island(g)-1)*(1-exp(-growthpars(1,island(g)-1)*(dat_growth(g,0)-growthpars(2,island(g)-1))));
    jnll_comp(7) -= dnorm(dat_growth(g,1), pred_len(g), pred_len(g)*CV_L, true);
  }

  // Calculate expected length and weight at age
  // 3rd column of growthpars has Little Cayman-specific growth pars
  vector<Type> L_a(n_a);
  vector<Type> W_a(n_a);
  for(int a=0; a<n_a; a++){
    L_a(a) = growthpars(0,2)*(1-exp(-growthpars(1,2)*(ages(a) - growthpars(2,2))));
    W_a(a) = lw_pars(0)*pow(L_a(a),lw_pars(1));
  }

  /////probability of being in a length bin given INPUT age
  matrix<Type> plba(n_a,n_lb);
  Type sum_sublast = 0;
  for(int a=0;a<n_a;a++){
    for(int l=0;l<n_lb;l++){
      if(l==0){
        plba(a,l) = pnorm(lbhighs(l), L_a(a), L_a(a)*CV_L);
        sum_sublast += plba(a,l);
      }
      if(l>=1){
        if(l<(n_lb-1)){
          plba(a,l) = pnorm(lbhighs(l), L_a(a), L_a(a)*CV_L) - pnorm(lbhighs(l-1), L_a(a), L_a(a)*CV_L);
          sum_sublast += plba(a,l);
        }
        if(l==(n_lb-1)) plba(a,l) = Type(1.0) - sum_sublast;
      }
    }
    sum_sublast = 0;
  }

  // ========= By fleet  =============================

  // Transform vectors
  matrix<Type> F_fy(n_fl,n_y); // number of total time steps
  // int index;
  for(int f=0;f<n_fl;f++){
    for(int y=0;y<n_y;y++){
        if(est_totalF==0) F_fy(f,y) = exp(log_F_ft(f,y));
        if(est_totalF==1) F_fy(f,y) = exp(log_F_ft(0,y)) * prop_f(f);
    }
  }

  // redo seasonal to include option to estimate annual or seasonal F (mirror_F)
  //with high enough sample size, could estimate seasonal F
  matrix<Type> F_ft(n_fl,n_t); // annually
  int tmp;
  for(int f=0;f<n_fl;f++){
    for(int t=0;t<n_t;t++){
      tmp = S_yrs(t) - 1;
      // match fishing mortality in the total number of time steps to the annual parameter F
      F_ft(f,t) = F_fy(f,tmp)/n_s;
    }
  }

  // total F
  vector<Type> F_t(n_t);
  F_t.setZero();
  for(int f=0;f<n_fl;f++){
    for(int t=0;t<n_t;t++){
      F_t(t) += F_ft(f,t);
    }
  }
  vector<Type> F_y(n_y);
  F_y.setZero();
  for(int f=0;f<n_fl;f++){
    for(int y=0;y<n_y;y++){
      F_y(y) += F_fy(f,y);
    }
  }

  // selectivity at length
  matrix<Type> S_fl(n_fl,n_lb);
  S_fl.setZero();
  for(int f=0;f<n_fl;f++){
    for(int l=0;l<n_lb;l++){
      if(selex_type_f(f)==0) S_fl(f,l) = vals_selex_ft(f,l);
      if(selex_type_f(f)==1) S_fl(f,l) = 1 / (1 + exp(-log(Type(19))*(lbmids(l) - S50_f(f))/(S95_f(f) - S50_f(f))));
    }   
  }

  // selectivity at age
  matrix<Type> S_fa(n_fl,n_a);
  S_fa.setZero();
  vector<Type> sub_plba(n_lb);
  sub_plba.setZero();
  for(int f=0;f<n_fl;f++){
    for(int a=0;a<n_a;a++){
      sub_plba = plba.row(a);
      for(int l=0;l<n_lb;l++){
        S_fa(f,a) += sub_plba(l)*S_fl(f,l);
      }
      sub_plba.setZero();
    }     
  }

  // maturity at age
  vector<Type> Mat_a(n_a);
  Mat_a.setZero();
  for(int a=0; a<n_a; a++){
    sub_plba = plba.row(a);
    for(int l=0; l<n_lb; l++){
      Mat_a(a) += sub_plba(l) * Mat_l(l);
    }
    sub_plba.setZero();
  }         

  // ============ equilibrium spawning biomass ===============
  // sum up relative spawning biomass at the beginning of each year
  // exp(beta)*Type(n_s) - estimates of beta were biased low due to spreading annual recruitment throughout the year
  // -M * Type(n_s) use annual estimate of M, not seasonal
  int n_a2;
  n_a2 = match_ages.size();
  Type SB0 = 0;
  for(int a=0;a<n_a;a++){
    // for(int a2=0;a2<n_a2;a2++){
      SB0 += (exp(beta) * exp(-M * Type(a))) * W_a(a) * Mat_a(a);
    // }
    
  }
  
  // // ============ joint F rate including selectivity ===========================

  // calculate fishing mortality at age over time by fleet
  // includes fishing pressure by fleet, selectivity, and catchability by fleet
  //catchability could be sinusoidal by quarter 
  array<Type> F_atf(n_a,n_t,n_fl);
  for(int f=0;f<n_fl;f++){
    for(int t=0;t<n_t;t++){
      for(int a=0;a<n_a;a++){
        F_atf(a,t,f) = F_ft(f,t) * S_fa(f,a);
      }
    }
  }

  // combine total fishing mortality at age over time across fleets
  matrix<Type> F_ta(n_t,n_a);
  F_ta.setZero();
  for(int f=0;f<n_fl;f++){
    for(int t=0;t<n_t;t++){
      for(int a=0;a<n_a;a++){
        F_ta(t,a) += F_atf(a,t,f);     
      }      
    }
  }

  // =========================== Population dynamics ===============================  

  // ============ initialize =============  
  vector<Type> R_t(n_t); //recruitment
  R_t.setZero();
  Type R0 = exp(beta);
  R_t(0) = exp(beta) * exp(Nu_input(0) - pow(sigma_R,2)/Type(2));
  // R_t(0) = ((4 * h * exp(beta) * SB_t(0)) / (SB0 * (1-h) + SB_t(0) * (5*h-1))) * exp(Nu_input(S_yrs(0)-1) - pow(sigma_R,2)/Type(2));
  // R_t(0) = exp(beta) * exp(-pow(sigma_R,2)/Type(2));

  //over time by age
  matrix<Type> N_ta(n_t,n_a); // abundance
  matrix<Type> SB_ta(n_t,n_a); // spawning biomass
  matrix<Type> TB_ta(n_t,n_a); // total biomass
  matrix<Type> SN_ta(n_t,n_a); // number of spawners at age a, time t
  N_ta.setZero();
  SB_ta.setZero();
  TB_ta.setZero();
  SN_ta.setZero(); 

  //over time
  vector<Type> N_t(n_t); // abundance
  vector<Type> SB_t(n_t); //spawning biomass
  vector<Type> TB_t(n_t); //total biomass
  vector<Type> SN_t(n_t); //number of spawners
  N_t.setZero();
  SB_t.setZero();
  TB_t.setZero();
  SN_t.setZero();  

  N_ta(0,0) = R_t(0);
  // N_ta(0,1) = R0 * exp(-0.5*pow(sigma_R,2)) * exp(-M - F_ta(0,0));
  N_ta(0,1) = R0 * exp(-M - F_ta(0,0));
  for(int a=0;a<n_a;a++){
    // Population abundance
    // if(a==0){
    //   N_ta(0,a) = R_t(0);
    // }
    // if((a>=1) & (a<(n_a-1))){
    //   N_ta(0,a) = N_ta(0,a-1) * exp(-M - F_ta(0,a-1));
    // }
    if((a>=2) & (a<(n_a-1))){
      N_ta(0,a) = N_ta(0,a-1) * exp(-M - F_ta(0,a-1));
    }    
    if(a==(n_a-1)){
      N_ta(0,a) = (N_ta(0,a-1) * exp(-M - F_ta(0,a-1))) / (1 - exp(-M - F_ta(0,a-1)));
    }

    // Spawning biomass
    SB_ta(0,a) = N_ta(0,a) * Mat_a(a) * W_a(a);

    // Spawning numbers
    SN_ta(0,a) = N_ta(0,a) * Mat_a(a);    

    // Total biomass
    TB_ta(0,a) = N_ta(0,a) * W_a(a);

    //Annual values
    if(a>0) N_t(0) += N_ta(0,a);
    SB_t(0) += SB_ta(0,a);
    SN_t(0) += SN_ta(0,a);
    TB_t(0) += TB_ta(0,a);
  }

  //catch by time, age, fleet
  array<Type> Cn_taf(n_t,n_a,n_fl);
  array<Type> Cw_taf(n_t,n_a,n_fl);
  Cn_taf.setZero();
  Cw_taf.setZero();

  //catch by time, age
  matrix<Type> Cn_ta(n_t,n_a); //catch in numbers
  matrix<Type> Cw_ta(n_t,n_a); //catch in biomass
  Cn_ta.setZero();
  Cw_ta.setZero();

  //catch by fleet, time
  matrix<Type> Cn_ft(n_fl,n_t);
  matrix<Type> Cw_ft(n_fl,n_t);
  Cn_ft.setZero();
  Cw_ft.setZero();

  //catch by time
  vector<Type> Cn_t_hat(n_t); //catch in numbers
  vector<Type> Cw_t_hat(n_t); //catch in biomass
  Cn_t_hat.setZero();
  Cw_t_hat.setZero();

  //year 1
  for(int f=0;f<n_fl;f++){
    for(int a=0;a<n_a;a++){
      Cn_taf(0,a,f) = N_ta(0,a) * (Type(1.0) - exp(-M - F_atf(a,0,f))) * (F_atf(a,0,f))/(M + F_atf(a,0,f));
      Cw_taf(0,a,f) = W_a(a) * Cn_taf(0,a,f);

      Cn_ta(0,a) += Cn_taf(0,a,f);
      Cw_ta(0,a) += Cw_taf(0,a,f);

      Cn_ft(f,0) += Cn_taf(0,a,f);
      Cw_ft(f,0) += Cw_taf(0,a,f);

      Cn_t_hat(0) += Cn_taf(0,a,f);
      Cw_t_hat(0) += Cw_taf(0,a,f);
    }
  }

  // ============ project forward in time =============  
  for(int t=1;t<n_t;t++){
    
    // Age-structured dynamics
    for(int a=1; a<n_a-1; a++) N_ta(t,a) = N_ta(t-1,a-1) * exp(-M - F_ta(t,a-1));
    N_ta(t,n_a-1) = (N_ta(t-1,n_a-2) * exp(-M - F_ta(t,n_a-2))) + (N_ta(t-1,n_a-1) * exp(-M - F_ta(t,n_a-1)));

    for(int a=1; a<n_a; a++){
      SB_ta(t,a) = N_ta(t,a)*Mat_a(a)*W_a(a); // Spawning biomass
      SN_ta(t,a) = N_ta(t,a)*Mat_a(a); // Spawning numbers 
      TB_ta(t,a) = N_ta(t,a)*W_a(a); // Total biomass

      // Annual totals across ages - recruits (age-0) contribute later
      N_t(t) += N_ta(t,a);
      SB_t(t) += SB_ta(t,a);
      TB_t(t) += TB_ta(t,a);
      SN_t(t) += SN_ta(t,a);
    }

    // Recruitment in year t is a function of SSB in year t, not t-1 (measure adults in year t, they produce age-0 right then)
    // but this shouldn't matter bc h = 1 (BH function isn't used)
    R_t(t) = ((4 * h * exp(beta) * SB_t(t)) / (SB0 * (1-h) + SB_t(t) * (5*h-1))) * exp(Nu_input(S_yrs(t)-1) - pow(sigma_R,2)/Type(2));
    if(S_yrs(t) == S_yrs(t-1)) R_t(t) = 0;
    N_ta(t,0) = R_t(t);
    N_t(t) += N_ta(t,0);
    TB_t(t) += TB_ta(t,0);
  }

  //catch by fleet
  for(int f=0;f<n_fl;f++){
    for(int t=1;t<n_t;t++){
      for(int a=0;a<n_a;a++){
        Cn_taf(t,a,f) = N_ta(t-1,a) * (Type(1.0) - exp(-M - F_atf(a,t,f))) * (F_atf(a,t,f))/(M + F_atf(a,t,f));
        Cw_taf(t,a,f) = Cn_taf(t,a,f) * W_a(a);  

        Cn_ta(t,a) += Cn_taf(t,a,f);
        Cw_ta(t,a) += Cw_taf(t,a,f);

        Cn_ft(f,t) += Cn_taf(t,a,f);
        Cw_ft(f,t) += Cw_taf(t,a,f);

        Cn_t_hat(t) += Cn_taf(t,a,f);
        Cw_t_hat(t) += Cw_taf(t,a,f);
      }      
    }
  }

  // // ========== Length composition ================================

  // probability of being harvested at an age
  array<Type> page(n_t,n_a,n_fl);
  for(int f=0;f<n_fl;f++){
    for(int t=0;t<n_t;t++){
      for(int a=0;a<n_a;a++){
        page(t,a,f) = (N_ta(t,a)*S_fa(f,a))/(N_t(t));
      }
    }    
  }

  // probability of sampling given length bin
  matrix<Type> plb_init(n_t,n_lb);
  matrix<Type> page_temp(n_t,n_a);
  matrix<Type> plb_sums(n_fl,n_t);
  array<Type> plb(n_t,n_lb,n_fl);

  for(int f=0;f<n_fl;f++){
    //find page for fleet f saved as page_temp
    for(int t=0;t<n_t;t++){
      for(int a=0;a<n_a;a++){
          page_temp(t,a) = page(t,a,f);
      }
    }
    //for each fleet f:
    plb_init = page_temp*plba;
    for(int t=0;t<n_t;t++){
      plb_sums(f,t) = 0;
      for(int l=0;l<n_lb;l++){
        // if(plb_init(t,l)==0) plb_init(t,l) = 1e-20;
        plb_sums(f,t) += plb_init(t,l) + Type(1e-20);
      }
    }
    for(int t=0;t<n_t;t++){
      for(int l=0;l<n_lb;l++){
        plb(t,l,f) = plb_init(t,l)/plb_sums(f,t);
      }
    }    
  }

  // expected mean length each year
  matrix<Type> ML_ft_hat(n_fl,n_t);
  matrix<Type> Vul_pop(n_fl,n_t);
  Type temp = 0;
  Type temp2 = 0;
  for(int f=0;f<n_fl;f++){
    for(int t=0;t<n_t;t++){
      for(int a=0;a<n_a;a++){
        temp2 += N_ta(t,a)*S_fa(f,a);
      }
      Vul_pop(f,t) = temp2;
      temp2 = 0;
      for(int l=0;l<n_lb;l++){
         temp += Vul_pop(f,t)*plb(t,l,f)*lbmids(l);
      }
      ML_ft_hat(f,t) = temp/Vul_pop(f,t);
      temp = 0;
    }    
  }


  // // ========= spawning potential ratio ==============================
  vector<Type> SPR_t(n_t);
  SPR_t.setZero();
  vector<Type> SPRequil_t(n_t);
  SPRequil_t.setZero();  

  //no relation to last year
  matrix<Type> NFequil_ta(n_t,n_a);
  matrix<Type> NUequil_ta(n_t,n_a);
  vector<Type> SBF_t(n_t);
  vector<Type> SBU_t(n_t);
  SBF_t.setZero();
  SBU_t.setZero();
  for(int t=0;t<n_t;t++){
    for(int a=0;a<n_a;a++){
      if(a==0){
        // NFequil_ta(t,a) = 1;
        // NUequil_ta(t,a) = 1;
        NFequil_ta(t,a) = R0;
        NUequil_ta(t,a) = R0;        
      }
      if((a>=1) & (a<(n_a-1))){
        NFequil_ta(t,a) = NFequil_ta(t,a-1) * exp(-M - F_ta(t,a-1));
        NUequil_ta(t,a) = NUequil_ta(t,a-1) * exp(-M);
      }
      if(a==(n_a-1)){
        NFequil_ta(t,a) = (NFequil_ta(t,a-1) * exp(-M - F_ta(t,a-1))) / (1 - exp(-M - F_ta(t,a-1)));
        NUequil_ta(t,a) = (NUequil_ta(t,a-1) * exp(-M)) / (1-exp(-M));
      }
    }

    for(int a=0;a<n_a;a++){
      SBF_t(t) += NFequil_ta(t,a) * Mat_a(a) * W_a(a);
      SBU_t(t) += NUequil_ta(t,a) * Mat_a(a) * W_a(a);
    }

    SPR_t(t) = SB_t(t)/SBU_t(t);
    SPRequil_t(t) = SBF_t(t)/SBU_t(t);
  }


  // Depletion
  vector<Type> D_t(n_t);
  D_t.setZero();
  for(int t=0;t<n_t;t++){
     // D_t(t) = SB_t(t)/SB0; // SB0 calculated incorrectly on line 
    D_t(t) = SB_t(t)/SB_t(0);
  }

  // ========= Build likelihood ==============================

  // Likelihood contribution from observations
  matrix<Type> log_pL_ft(n_fl,n_t);
  log_pL_ft.setZero();
  matrix<Type> neff(n_fl,n_t);
  neff.setZero();
  Type checklc;
    vector<Type> LFprob(n_lb);
    vector<Type> LFraw(n_lb);
    vector<Type> prob(n_lb);
    vector<Type> sum1(n_t);
    vector<Type> sum2(n_t);
    sum1.setZero();
    sum2.setZero();

  for(int f=0;f<n_fl;f++){
    for(int t=0;t<n_t;t++){

      // for each new time step, set check if LC data exists to zero
      checklc = 0;

      // for each new time step, set observed and predicted data in each length bin to zero
      prob.setZero();
      LFraw.setZero();

      // read data for each length bin
      for(int l=0;l<n_lb;l++){
        prob(l) = plb(t,l,f);
        LFraw(l) = LF_tlf(t,l,f);
        checklc += LFraw(l);
      }

      // if lengths were observed this year, calculate likelihood
      if(checklc > 0){
        // multinomial
        if(LFdist==0){
          LFprob = n_LF_ft(f,t)*(LFraw/LFraw.sum());
          log_pL_ft(f,t) += dmultinom(LFprob, prob, true);
        }
        // dirichlet-multinomial
        if(LFdist==1){
          LFprob = LFraw/(LFraw.sum());
          for(int l=0;l<n_lb;l++){
            sum1(t) += lgamma(LFraw.sum()*LFprob(l)+1);
            sum2(t) += lgamma(LFraw.sum()*LFprob(l) + theta(f)*LFraw.sum()*prob(l)) - lgamma(theta(f)*LFraw.sum()*prob(l));
          }
          log_pL_ft(f,t) += lgamma(LFraw.sum()+1) - sum1(t) + lgamma(theta(f)*LFraw.sum()) - lgamma(LFraw.sum() + theta(f)*LFraw.sum()) + sum2(t);
        }
      }

    }
  }

  if(LFdist==1){
    for(int f=0;f<n_fl;f++){
      for(int t=0;t<n_t;t++){

        LFraw.setZero();
        checklc = 0;
        for(int l=0;l<n_lb;l++){
          LFraw(l) = LF_tlf(t,l,f);
          checklc += LFraw(l);
        }
        if(checklc > 0) neff(f,t) = (1 + theta(f)*LFraw.sum())/(1+theta(f));      
      }
    }
  }

  matrix<Type> I_ft_hat(n_fl,n_t);
  for(int f=0;f<n_fl;f++){
    for(int t=0;t<n_t;t++){
        // I_ft_hat(f,t) = q_f(f)*TB_t(t);
        I_ft_hat(f,t) = q_f(f)*SN_t(t); // index is numbers of spawners, not biomass and no immature fish
    }
  }

  matrix<Type> log_pI_ft(n_fl,n_t);
  log_pI_ft.setZero();
  for(int f=0;f<n_fl;f++){
    for(int t=0;t<n_t;t++){
        // probability of index at that sample
        if(I_ft(f,t)>0) log_pI_ft(t) += dlognorm( I_ft(f,t), log(I_ft_hat(f,t)), sigma_I, true);
    }
  }
  matrix<Type> log_index_resid = log(I_ft.array()) - log(I_ft_hat.array());
  ADREPORT(log_index_resid);

  matrix<Type> log_pC_ft(n_fl,n_t);
  log_pC_ft.setZero();
  for(int f=0;f<n_fl;f++){
    for(int t=0;t<n_t;t++){  
          // probability of index at that sample
        if(C_ft(f,t)>=0){
          if(C_type==1) log_pC_ft(f,t) += dlognorm( C_ft(f,t), log(Cn_ft(f,t)), sigma_C, true);
          if(C_type==2) log_pC_ft(f,t) += dlognorm( C_ft(f,t), log(Cw_ft(f,t)), sigma_C, true);
        }
    }
  }
  
  matrix<Type> log_pML_ft(n_fl,n_t);
  log_pML_ft.setZero();
  // for(int f=0;f<n_fl;f++){
  //   for(int t=0;t<n_t;t++){
  //         if(ML_ft(f,t)>0){
  //           log_pML_ft(f,t) += dnorm( ML_ft(f,t), ML_ft_hat(f,t), ML_ft_hat(f,t)*CV_L, true);
  //         }
  //   }
  // }
  matrix<Type> ML_resid(n_fl,n_t);
  for(int f=0;f<n_fl;f++){
    for(int t=0;t<n_t;t++){
      if(ML_years(t) == 1) ML_resid(f,t) = ML_t(t) - ML_ft_hat(f,t);
    }
  }
  ADREPORT(ML_resid);

  // ------ likelihood components ----------//
  // ============ Probability of random effects =============
  jnll_comp(0) = 0;
    for(int y=0;y<n_y;y++){
  // for(int y=0;y<n_y-1;y++){
      jnll_comp(0) -= dnorm(Nu_input(y), Type(0.0), sigma_R, true);
    }    


  // ============ Probability of data =============
  jnll_comp(1) = 0;
  jnll_comp(1) = Type(-1)*sum( log_pL_ft );
  jnll_comp(2) = 0;
  jnll_comp(2) = Type(-1)*sum( log_pML_ft ); 

  jnll_comp(3) = 0;
  jnll_comp(3) = Type(-1)*sum( log_pI_ft );
  jnll_comp(4) = 0;
  jnll_comp(4) = Type(-1)*sum( log_pC_ft );

    // F
    jnll_comp(5) = 0;
    if(Fpen==1){
        // for(int y=4;y<n_y;y++) jnll_comp(5) -= dnorm(F_y(y), F_y(y-1), sigma_F, true);
        for(int y=5;y<n_y;y++) jnll_comp(5) -= dnorm(F_y(y), F_y(y-1), sigma_F, true);
        // for(int t=1;t<n_t;t++) jnll_comp(5) -= dnorm(F_t(t), F_t(t-1), sigma_F, true);
    }

    // SigmaR
    Type sigrp;
    sigrp = 0;
    jnll_comp(6) = 0;
    sigrp = dlognorm(sigma_R, log(SigRprior(0)), SigRprior(1), true);
    if(SigRpen==1) jnll_comp(6) = Type(-1)*sigrp;
    // // penalize R(2012) and R(2010) relative to R(2011). Rec(2011) = R(12) = 12th Nu_input = Nu_input(11)
    // vector<Type> RecPen(2);
    // RecPen(0) = Type(0);
    // RecPen(1) = Nu_input(12)-Nu_input(11); // Rec(2012) - Rec(2011) should be negative. Add penalty if not.
    // jnll_comp(6) += Type(500)*pow(RecPen.maxCoeff(),2);
    // jnll_comp(6) += Type(500)*(Nu_input(12)-Nu_input(11));
    // jnll_comp(6) -= pow(max(Type(0), Nu_input(10)-Nu_input(11)),2);

    jnll = sum(jnll_comp);

    vector<Type> lN_t(n_t);
    vector<Type> lSB_t(n_t);
    vector<Type> lSN_t(n_t);    
    vector<Type> lTB_t(n_t);
    vector<Type> lR_t(n_t);
    vector<Type> lF_t(n_t);
    matrix<Type> lF_ft(n_fl,n_t);
    matrix<Type> lF_fy(n_fl,n_y);
    matrix<Type> lC_ft(n_fl,n_t);
    matrix<Type> lI_ft(n_fl,n_t);
    vector<Type> lD_t(n_t);
    for(int t=0;t<n_t;t++){
      lN_t(t) = log(N_t(t));
      lSB_t(t) = log(SB_t(t));
      lSN_t(t) = log(SN_t(t));      
      lTB_t(t) = log(TB_t(t));
      lR_t(t) = log(R_t(t));
      lF_t(t) = log(F_t(t));
      lD_t(t) = log(D_t(t));
    }
    for(int f=0;f<n_fl;f++){
      for(int t=0;t<n_t;t++){
        if(C_type==0) lC_ft(f,t) = log(Cw_ft(f,t));
        if(C_type==1) lC_ft(f,t) = log(Cn_ft(f,t));
        if(C_type==2) lC_ft(f,t) = log(Cw_ft(f,t));
        lI_ft(f,t) = log(I_ft_hat(f,t));
        lF_ft(f,t) = log(F_ft(f,t));
      }
      for(int y=0;y<n_y;y++){
        lF_fy(f,y) = log(F_fy(f,y));
      }
    }

    vector<Type> lF_y(n_y);
    for(int y=0;y<n_y;y++){
      lF_y(y) = log(F_y(y));
    }



  // // ============ Reporting section ======================================
  // Type R2011 = R_t(12)/R0; // 2011 recruitment / mean recruitment (get SE this way)
  Type muR = exp(beta - pow(sigma_R,2)/Type(2));
  Type R2011 = R_t(12) / muR;
  ADREPORT(R2011);
  ADREPORT(R0);
  ADREPORT(muR);

  ADREPORT( lN_t );
  ADREPORT( lR_t );
  ADREPORT( lI_ft );
  ADREPORT( lC_ft );
  ADREPORT( ML_ft_hat );
  ADREPORT( I_ft_hat );
  ADREPORT( SN_t );  
  ADREPORT( lSB_t );
  ADREPORT( lF_t );
  ADREPORT( lF_ft );
  ADREPORT( lF_y );
  ADREPORT( lF_fy);
  ADREPORT( lD_t );
  ADREPORT( SPR_t );
  ADREPORT( SPRequil_t );
  ADREPORT( S50_f );
  ADREPORT( S95_f );
  ADREPORT(vbk); 
  ADREPORT(linf);
  ADREPORT(t0);
  ADREPORT(growthpars);
  // ADREPORT( S_fl );

  // Parameters
  REPORT( q_f );
  REPORT( sigma_F );
  REPORT( beta );
  REPORT( sigma_R );
  REPORT( log_sigma_R );
  REPORT( S50_f );
  REPORT( S95_f );
  REPORT( S_fa );
  REPORT( S_fl );
  REPORT( sigma_C );
  REPORT( sigma_I );
  REPORT( CV_L );
  REPORT( SPR_t );
  REPORT( SPRequil_t );

  // Random effects
  REPORT( Nu_input );

   // State variables
  REPORT( R_t );
  REPORT( F_t );
  REPORT( F_y );
  REPORT( N_t );
  REPORT( F_ft );
  REPORT( prop_f );
  REPORT( F_fy );
  REPORT( F_ta );

  // Predicted quantities
  REPORT( ML_ft_hat );
  REPORT( I_ft_hat );
  REPORT(pred_len);
  // Derived quantities
  REPORT( Cn_t_hat );
  REPORT( Cw_t_hat );
  REPORT( Cn_ft );
  REPORT( Cw_ft );
  REPORT( SB_t );
  REPORT( SN_t );  
  REPORT( TB_t );
  REPORT( SB0 );
  REPORT(D_t);
  REPORT(N_ta);
  REPORT(Cn_ta);
  REPORT(plba);
  REPORT(page);
  REPORT(plb);
  ADREPORT(W_a);
  ADREPORT(L_a);
  ADREPORT(Mat_a);
  REPORT(M);
  REPORT(h);
  // growth pars
  REPORT(Sigma);
  REPORT(growth_mu);
  REPORT(emat);
  REPORT(log_sigma_mu);
    // Likelihoods
  REPORT(log_pC_ft);
  REPORT(log_pI_ft);
  REPORT(log_pL_ft);
  REPORT(neff);
  REPORT(theta);
  REPORT(log_pML_ft);
  REPORT(sigrp);
  REPORT(jnll_comp);
  REPORT(jnll); 
  return(jnll);
}
