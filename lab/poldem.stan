// Bollen's Political Democracy SEM with a mean structure, written
// independently of inlavanilla's log-posterior. Priors are INLAvaan's (and
// blavaan's) defaults, on the natural scale; Stan adds the Jacobians.
//
// The MCMC reference in lab/DESIGN.md (October 2026): CmdStan 2.37, 4 chains
// of 2,000 warm-up and 10,000 draws, adapt_delta = 0.95, R-hat 1.00, no
// divergences. Start the chains at the posterior mode: from Stan's default
// random inits, two of four chains stuck in a spurious mode (lam_x2 about
// -20, var(ind60) about 0.006) with lp__ 73 lower.
data {
  int<lower=1> N;
  matrix[N, 11] Y; // x1 x2 x3 y1 ... y8
}
parameters {
  vector[11] nu;
  real lam_x2;
  real lam_x3;
  real a;
  real b;
  real c;
  real beta_dem60_ind60;
  real beta_dem65_ind60;
  real beta_dem65_dem60;
  vector<lower=0>[3] sd_lat; // ind60, dem60, dem65
  vector<lower=0>[11] sd_res;
  vector<lower=-1, upper=1>[6] rho; // y1~~y5 y2~~y4 y2~~y6 y3~~y7 y4~~y8 y6~~y8
}
transformed parameters {
  cov_matrix[11] Sigma;
  {
    matrix[11, 3] Lambda = rep_matrix(0, 11, 3);
    matrix[3, 3] B = rep_matrix(0, 3, 3);
    matrix[11, 11] Theta = diag_matrix(square(sd_res));
    array[6, 2] int pr = {{4, 8}, {5, 7}, {5, 9}, {6, 10}, {7, 11}, {9, 11}};
    matrix[3, 3] A;
    Lambda[1, 1] = 1; Lambda[2, 1] = lam_x2; Lambda[3, 1] = lam_x3;
    Lambda[4, 2] = 1; Lambda[5, 2] = a; Lambda[6, 2] = b; Lambda[7, 2] = c;
    Lambda[8, 3] = 1; Lambda[9, 3] = a; Lambda[10, 3] = b; Lambda[11, 3] = c;
    B[2, 1] = beta_dem60_ind60;
    B[3, 1] = beta_dem65_ind60;
    B[3, 2] = beta_dem65_dem60;
    for (k in 1:6) {
      real v = rho[k] * sd_res[pr[k, 1]] * sd_res[pr[k, 2]];
      Theta[pr[k, 1], pr[k, 2]] = v;
      Theta[pr[k, 2], pr[k, 1]] = v;
    }
    A = inverse(diag_matrix(rep_vector(1, 3)) - B);
    Sigma = quad_form_sym(diag_matrix(square(sd_lat)), (Lambda * A)') + Theta;
    Sigma = 0.5 * (Sigma + Sigma');
  }
}
model {
  nu ~ normal(0, 32);
  {lam_x2, lam_x3, a, b, c} ~ normal(0, 10);
  {beta_dem60_ind60, beta_dem65_ind60, beta_dem65_dem60} ~ normal(0, 10);
  sd_lat ~ gamma(1, 0.5);
  sd_res ~ gamma(1, 0.5);
  (rho + 1) / 2 ~ beta(1, 1);
  for (n in 1:N) {
    Y[n]' ~ multi_normal(nu, Sigma);
  }
}
generated quantities {
  real logv_x2 = 2 * log(sd_res[2]);
  real logv_dem65 = 2 * log(sd_lat[3]);
}
