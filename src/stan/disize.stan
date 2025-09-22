functions {
  real partial_posterior(
    // Thread-specific ----
    array[,] int counts_slice,
    int start,
    int end,
    // Shared ----
    int n_obs,
    int n_fe,
    matrix fe_design,
    int n_re,
    vector re_design_x,
    array[] int re_design_j,
    array[] int re_design_p,
    array[] int re_id,
    array[] int batch_id,
    vector intercept,
    matrix raw_fe_coefs,
    matrix raw_re_coefs,
    vector fe_tau,
    matrix fe_lambda,
    vector re_tau,
    matrix re_lambda,
    vector sf,
    real iodisp
  ) {
    real log_prob = 0;
    vector[n_obs] log_mu;
    vector[n_obs] fe_effect;
    vector[n_obs] batch_effect = sf[batch_id];

    
    int n_slice_feats = end - start + 1; // # of features in slice
    for (i in 1 : n_slice_feats) {
      int feat_i = start + i - 1;

      // Priors
      // Half-cauchy prior over random-effects variance
      log_prob += cauchy_lpdf(re_lambda[ : , feat_i] | 0, 1);

      // Horseshoe prior over fixed-effects
      log_prob += cauchy_lpdf(fe_lambda[ : , feat_i] | 0, 1);
      log_prob += std_normal_lpdf(raw_fe_coefs[ : , feat_i]);

      // Normal prior over (raw) random-effects
      log_prob += std_normal_lpdf(raw_re_coefs[ : , feat_i]);
      
      // Effect from the experimental design
      log_mu = rep_vector(intercept[feat_i], n_obs);
      
      // Fixed-effects
      if (n_fe != 0) {
        log_mu += fe_design * (
          raw_fe_coefs[ : , feat_i] .* (fe_lambda[ : , feat_i] .* fe_tau)
        );
      }

      // Random-effects
      if (n_re != 0) {
        log_mu += csr_matrix_times_vector(
          n_obs,
          n_re,
          re_design_x,
          re_design_j,
          re_design_p,
          raw_re_coefs[ : , feat_i] .* (re_lambda[re_id, feat_i] .* re_tau[re_id])
        );
      }

      // Adjust for batch-effect
      log_mu += batch_effect;

      // Likelihood
      log_prob += neg_binomial_2_log_lpmf(counts_slice[i] | log_mu, iodisp);
    }
    return log_prob;
  }
}
data {
  // Dimensions ----
  int<lower=1> n_obs; // # of observations
  int<lower=1> n_feats; // # of features
  int<lower=0> n_fe; // # of fixed-effects per-gene
  int<lower=0> n_re; // # of random-effects per-gene
  int<lower=0> n_nz_re; // # of nonzero elements in the random-effects design matrix
  int<lower=0> n_re_terms; // # of random-effects terms
  int<lower=1> n_batches; // # of batches
  
  // Design Matrices ----
  // fixed-effects
  matrix[n_obs, n_fe] fe_design;

  // random-effects
  vector[n_nz_re] re_design_x;
  array[n_nz_re] int re_design_j;
  array[n_obs + 1] int re_design_p;

  // Index Variables ----
  array[n_re] int<lower=1, upper=n_re_terms> re_id; // random-effects term membership
  array[n_obs] int<lower=1, upper=n_batches> batch_id; // batch membership

  // Response ----
  array[n_feats, n_obs] int<lower=0> counts; // counts matrix

  // Configuration for Multi-threading ----
  int<lower=1> grainsize;
}
parameters {
  // Feature Expression ----
  vector[n_feats] intercept;
  matrix[n_fe, n_feats] raw_fe_coefs;
  matrix[n_re, n_feats] raw_re_coefs;

  // Batch Effect ----
  simplex[n_batches] raw_sf;

  // Shrinkage ----
  vector<lower=0>[n_fe] fe_tau;
  matrix<lower=0>[n_fe, n_feats] fe_lambda;

  vector<lower=0>[n_re_terms] re_tau;
  matrix<lower=0>[n_re_terms, n_feats] re_lambda;

  // Inverse Overdispersion ----
  real<lower=0> iodisp;
}
transformed parameters {
  // Size Factors ----
  vector[n_batches] sf = log(raw_sf) + log(n_batches);
}
model {
  // Parallel posterior eval ----
  target += reduce_sum(
    partial_posterior,
    counts,
    grainsize,
    // Shared ----
    n_obs,
    n_fe,
    fe_design,
    n_re,
    re_design_x,
    re_design_j,
    re_design_p,
    re_id,
    batch_id,
    intercept,
    raw_fe_coefs,
    raw_re_coefs,
    fe_tau,
    fe_lambda,
    re_tau,
    re_lambda,
    sf,
    iodisp
  );
}
