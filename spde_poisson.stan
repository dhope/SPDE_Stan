// inspired from https://github.com/lionel68/spde_stan/blob/main/script/stan_spde.stan
data {
  int n;    // n obs
  int p;    // n par
  int n_knots;    // rows sparse matrix
  // matrix[n, 2] lat_lon;
  array[n] int y;         //The response
  matrix[n, p] X;         //Design matrix for fixed effects
  int n_non_zero_M;
  int n_non_zero_A;
  matrix[3,2]L;
  matrix[n_knots, n_knots ] M0;     // SPDE matrices from INLA
  matrix[n_knots, n_knots ] M1;
  matrix[n_knots, n_knots ] M2;
  matrix[ n,n_knots] A;     //Matrix for interpolating points witin triangles
  vector[2] prior_mean_tau_kappa_log;
  vector[2] prior_sd_tau_kappa_log;
  matrix[557, p] pred_x;
  matrix[557, n_knots] pred_mat;
}
transformed data{
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
    // vector[ n_knots] mus;
  }


//========================
  // T.PARAMETERS
//========================
  transformed parameters {
    vector[3] llmda = L * tau_kappa_log;
    vector[3] lmda = exp(llmda);
    //------------------------------------------
    vector[n] eta;
    real log_lik_u;
    
    // matrix[n_knots, n_knots] Q;
    vector[n] delta;
     vector[n_non_zero_M] Qa;
     vector[n_knots] pen_v;
    {
     
     
      // matrix[n_knots, n_knots] Q;
      Qa = lmda[1] * M0_v + lmda[2] * M1_v + lmda[3] * M2_v;
      pen_v = csr_matrix_times_vector(n_knots, n_knots,
                      Qa,col_id_M,row_id_M,
                       u);

      // log_lik_u = -1 * dot_product(pen_v, u);
      log_lik_u = -1 * (to_row_vector(pen_v) * u);
      
     
        // Equivalent to:
      // log_lik_u2 = u' * Q * u;
      // Q = csr_to_dense_matrix(n_knots, n_knots, Qa, col_id_M, row_id_M);
      
      delta = csr_matrix_times_vector(n,n_knots,A_sparse.1,
                      A_sparse.2,A_sparse.3,u); // spatial
                      
      eta = X * beta + delta ;
      
    }
  }

//========================
  // MODEL
//========================
  model {
    //---------------------------------------------
    u  ~ std_normal();//student_t(3, 0, 3.2);
    //
    beta[1] ~ normal(log(11.5),0.2);
    // beta[2] ~ normal(0,2);

    tau_kappa_log[1] ~ normal(prior_mean_tau_kappa_log[1],
                               prior_sd_tau_kappa_log[1]); // tau
    tau_kappa_log[2] ~ normal(prior_mean_tau_kappa_log[2],
                            prior_sd_tau_kappa_log[2]); // kappa
    
    target += log_lik_u;
    target += poisson_log_lpmf(y | eta);
    
    

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
    
    vector[557] lambda_p = pred_x * beta + pred_mat * u;
    
    
    
    
  }
