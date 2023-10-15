// inspired from https://github.com/lionel68/spde_stan/blob/main/script/stan_spde.stan
data {
  int n;    // n obs
  int p;    // n par
  int n_knots;    // rows sparse matrix
  matrix[n, 2] lat_lon;
  vector[n] y;         //The response
  matrix[n, p] X;         //Design matrix for fixed effects
  int n_non_zero_M;
  int n_non_zero_A;
  matrix[3,2]L;
  matrix[n_knots, n_knots ] M0;     // SPDE matrices from INLA
  matrix[n_knots, n_knots ] M1;
  matrix[n_knots, n_knots ] M2;
  matrix[ n,n_knots] A;     //Matrix for interpolating points witin triangles
}
transformed data{
  // vector[n_knots] zeroes = rep_vector(0, n_knots);
  vector[n_non_zero_M] M0_v = csr_extract_w(M0);
  vector[n_non_zero_M] M1_v = csr_extract_w(M1);
  vector[n_non_zero_M] M2_v = csr_extract_w(M2);
  array[(n_knots+1)] int row_id_M = csr_extract_u(M0);
  array[n_non_zero_M] int col_id_M = csr_extract_v(M0);
  tuple(vector[n_non_zero_A],array[n_non_zero_A] int,array[(n+1)] int) A_sparse;
  A_sparse = csr_extract(A);
  
}

//========================
  // PARAMETERS
//========================
  parameters {
    vector[p] beta; 
    vector[2] tau_kappa_log;
    vector[n_knots] u;   // spatial random effect  
    real<lower=0> sigma;
  }


//========================
  // T.PARAMETERS
//========================
  transformed parameters {
    vector[3] llmda = L * tau_kappa_log;
    vector[3] lmda = exp(llmda);
    //------------------------------------------
    vector[n] eta;
    vector[n_knots] beta_s;
    {
    vector[n_knots] beta_s2;
    vector[n] B_delta;
    vector[n_non_zero_M] Q;
    vector[n] delta;
    Q = lmda[1] * M0_v + lmda[2] * M1_v + lmda[3] * M2_v;

    // beta_s = Q*u;
    beta_s = csr_matrix_times_vector(n_knots, n_knots, Q,col_id_M,row_id_M,u );
    
    beta_s2 =  u - beta_s;
    
    delta = csr_matrix_times_vector(n,n_knots,A_sparse.1,
                      A_sparse.2,A_sparse.3,beta_s2); // spatial

    eta = X * beta  +   delta;
    }
  }

//========================
  // MODEL
//========================
  model {
    //---------------------------------------------
    u  ~ normal( 0 , 2 );
    sigma ~ exponential(1);
    beta ~ normal(0,1);
    tau_kappa_log[1] ~ normal(-2.58, 0.5); // tau
    tau_kappa_log[2] ~ normal(1, 0.5); // kappa
    y ~  normal(eta, sigma);
  }

//========================
  // G.QUANTITIES
//========================
  generated quantities{
    real range = sqrt(8)/exp(tau_kappa_log[2]);
    real tau = exp(tau_kappa_log[1]);
    real kappa = exp(tau_kappa_log[2]);
    real denom = (pow(tau,2)) * 4*pi() * (pow(kappa,2));
        // see Lindgren et al. (2011) for this formula
    real sigma_spde = 1 / (sqrt(denom));
  }
