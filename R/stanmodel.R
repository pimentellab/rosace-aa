#' @importFrom cmdstanr write_stan_file
#'
NULL

WriteStanModel <- function(type) {
  if (type == "growth_pos_blosum_act") {
    stan_program <- "
    data {
      int<lower=0> T; // # of time points
      int<lower=0> V; // # of variants
      int<lower=0> P; // # of positions
      int<lower=0> M; // # of mean group
      int<lower=0> B; // # of blosum group
      array[V] int vMAPp;
      array[V] int vMAPm;
      array[V] int vMAPb;
      vector[T] t; // time
      array[V] vector[T] m; // normalized count
      vector[B] blosum_count; // Count of each BLOSUM group (including SYN)
      int<lower=0> P_syn; // # of synonymous positions
    }
    transformed data {
      vector[B-1] diagon = rep_vector(0, B-1);
      diagon[1] = 1.0;
      vector[B-2] off_diagonal = rep_vector(0, B-2);
      real<lower=0> off_diag_ss = 0.;
      for (i in 1:B-2) {
        off_diagonal[i] = (-1.0 / (B - 2) - off_diag_ss) / diagon[i];
        off_diag_ss = off_diag_ss + off_diagonal[i] * off_diagonal[i];
        diagon[i+1] = sqrt(1.0 - off_diag_ss);
      }
      matrix[B-1, B-2] nu_multiplier = rep_matrix(0.0, B-1, B-2);
      for (i in 1:B-2) {
        nu_multiplier[i][i] = diagon[i];
        for (j in i:B-2) {
          nu_multiplier[j+1][i] = off_diagonal[i];
        }
      }
      vector[B-1] w = blosum_count[:B-1];
      real count_mean = mean(w);
      for (i in 1:B-1) {
        w[i] = count_mean / w[i];
      }
      matrix[B-1, B-1] weight_multiplier = diag_matrix(w);
      nu_multiplier = weight_multiplier * nu_multiplier;
    }
    parameters {
      vector[P] phi; // Slope per position
      vector<lower=0>[P] sigma2; // Noise per position
      vector<lower=0,upper=1>[P-P_syn] rho_wo_syn; // Scaling factor per position
      vector[B-2] nu_raw;
      vector<lower=0>[M] epsilon2;
      vector[V] eta2; // std_normal for beta
      vector[V] b; // intercept
    }
    transformed parameters {
      vector[P] rho = rep_vector(.5, P);
      for (i in 1:P-P_syn) {
        rho[i] = rho_wo_syn[i];
      }
      vector[V] beta; // total slope per variant (beta + xi)
      vector[B-1] nu_wo_syn = nu_multiplier * nu_raw; // BLOSUM group effect (excl. synonymous)
      vector[B] nu = append_row(nu_wo_syn, 0); // nu with a 0 appended to the end
      for (v in 1:V) {
        beta[v] = phi[vMAPp[v]] + rho[vMAPp[v]] * nu[vMAPb[v]] + eta2[v] * sqrt(sigma2[vMAPp[v]]);
      }
    }
    model {
      phi ~ normal(0, 1);
      nu_raw ~ normal(0, 0.5);
      rho_wo_syn ~ beta(1.5, 1.5);
      sigma2 ~ inv_gamma(1, 1); 
      eta2 ~ normal(0, 1);
      epsilon2 ~ inv_gamma(1, 1);
      b ~ normal(0, 0.25);
      for (v in 1:V) {
        m[v] ~ normal(b[v] + beta[v] * t, sqrt(epsilon2[vMAPm[v]]));
      }
    }
    "
  } else if (type == "growth_pos_blosum_act_nosyn") {
    stan_program <- "
    data {
      int<lower=0> T; // # of time points
      int<lower=0> V; // # of variants
      int<lower=0> P; // # of positions
      int<lower=0> M; // # of mean group
      int<lower=0> B; // # of blosum group
      array[V] int vMAPp;
      array[V] int vMAPm;
      array[V] int vMAPb;
      vector[T] t; // time
      array[V] vector[T] m; // normalized count
      vector[B] blosum_count; // Count of each BLOSUM group (including SYN, no SYN this version)
    }
    transformed data {
      vector[B] diagon = rep_vector(0, B);
      diagon[1] = 1.0;
      vector[B-1] off_diagonal = rep_vector(0, B-1);
      real<lower=0> off_diag_ss = 0.;
      for (i in 1:B-1) {
        off_diagonal[i] = (-1.0 / (B - 1) - off_diag_ss) / diagon[i];
        off_diag_ss = off_diag_ss + off_diagonal[i] * off_diagonal[i];
        diagon[i+1] = sqrt(1.0 - off_diag_ss);
      }
      matrix[B, B-1] nu_multiplier = rep_matrix(0.0, B, B-1);
      for (i in 1:B-1) {
        nu_multiplier[i][i] = diagon[i];
        for (j in i:B-1) {
          nu_multiplier[j+1][i] = off_diagonal[i];
        }
      }
      vector[B] w = blosum_count[:B];
      real count_mean = mean(w);
      for (i in 1:B) {
        w[i] = count_mean / w[i];
      }
      matrix[B, B] weight_multiplier = diag_matrix(w);
      nu_multiplier = weight_multiplier * nu_multiplier;
    }
    parameters {
      vector[P] phi; // Slope per position
      vector<lower=0>[P] sigma2; // Noise per position
      vector<lower=0,upper=1>[P] rho; // Scaling factor per position
      vector[B-1] nu_raw;
      vector<lower=0>[M] epsilon2;
      vector[V] eta2; // std_normal for beta
      vector[V] b; // intercept
    }
    transformed parameters {
      vector[V] beta; // total slope per variant (beta + xi)
      vector[B] nu = nu_multiplier * nu_raw; // BLOSUM group effect (no synonymous)
      for (v in 1:V) {
        beta[v] = phi[vMAPp[v]] + rho[vMAPp[v]] * nu[vMAPb[v]] + eta2[v] * sqrt(sigma2[vMAPp[v]]);
      }
    }
    model {
      phi ~ normal(0, 1);
      nu_raw ~ normal(0, 0.5);
      rho ~ beta(1.5, 1.5);
      sigma2 ~ inv_gamma(1, 1); 
      eta2 ~ normal(0, 1);
      epsilon2 ~ inv_gamma(1, 1);
      b ~ normal(0, 0.25);
      for (v in 1:V) {
        m[v] ~ normal(b[v] + beta[v] * t, sqrt(epsilon2[vMAPm[v]]));
      }
    }
    "
  } else if (type == "growth_pos_blosum") {
    stan_program <- "
    data {
      int<lower=0> T; // # of time points
      int<lower=0> V; // # of variants
      int<lower=0> P; // # of positions
      int<lower=0> M; // # of mean group
      int<lower=0> B; // # of blosum group
      array[V] int vMAPp;
      array[V] int vMAPm;
      array[V] int vMAPb;
      vector[T] t; // time
      array[V] vector[T] m; // normalized count
      vector[B] blosum_count; // Count of each BLOSUM group (including SYN)
    }
    transformed data {
      vector[B-1] diagon = rep_vector(0, B-1);
      diagon[1] = 1.0;
      vector[B-2] off_diagonal = rep_vector(0, B-2);
      real<lower=0> off_diag_ss = 0.;
      for (i in 1:B-2) {
        off_diagonal[i] = (-1.0 / (B - 2) - off_diag_ss) / diagon[i];
        off_diag_ss = off_diag_ss + off_diagonal[i] * off_diagonal[i];
        diagon[i+1] = sqrt(1.0 - off_diag_ss);
      }
      matrix[B-1, B-2] nu_multiplier = rep_matrix(0.0, B-1, B-2);
      for (i in 1:B-2) {
        nu_multiplier[i][i] = diagon[i];
        for (j in i:B-2) {
          nu_multiplier[j+1][i] = off_diagonal[i];
        }
      }
      vector[B-1] w = blosum_count[:B-1];
      real count_mean = mean(w);
      for (i in 1:B-1) {
        w[i] = count_mean / w[i];
      }
      matrix[B-1, B-1] weight_multiplier = diag_matrix(w);
      nu_multiplier = weight_multiplier * nu_multiplier;
    }
    parameters {
      vector[P] phi; // slope per position
      vector[B-2] nu_raw;
      vector<lower=0>[P] sigma2;
      vector<lower=0>[M] epsilon2;
      vector[V] eta2; // std_normal for beta
      vector[V] b; // intercept
    }
    transformed parameters {
      vector[B-1] nu_wo_syn = nu_multiplier * nu_raw;
      vector[B] nu = append_row(nu_wo_syn, 0); // nu with a 0 appended to the end
      vector[V] beta; // slope per variants
      for (v in 1:V) {
        beta[v] = phi[vMAPp[v]] + nu[vMAPb[v]] + eta2[v] * sqrt(sigma2[vMAPp[v]]);
      }
    }
    model {
      phi ~ normal(0, 1);
      nu_raw ~ normal(0, 0.5);
      sigma2 ~ inv_gamma(1, 1);
      eta2 ~ normal(0, 1);
      epsilon2 ~ inv_gamma(1, 1);
      b ~ normal(0, 0.25);
      for (v in 1:V) {
        m[v] ~ normal(b[v] + beta[v] * t, sqrt(epsilon2[vMAPm[v]]));
      }
    }
    "
  } else if (type == "growth_pos_blosum_nosyn") {
    stan_program <- "
    data {
      int<lower=0> T; // # of time points
      int<lower=0> V; // # of variants
      int<lower=0> P; // # of positions
      int<lower=0> M; // # of mean group
      int<lower=0> B; // # of blosum group
      array[V] int vMAPp;
      array[V] int vMAPm;
      array[V] int vMAPb;
      vector[T] t; // time
      array[V] vector[T] m; // normalized count
      vector[B] blosum_count; // Count of each BLOSUM group (no SYN in this model)
    }
    transformed data {
      vector[B] diagon = rep_vector(0, B);
      diagon[1] = 1.0;
      vector[B-1] off_diagonal = rep_vector(0, B-1);
      real<lower=0> off_diag_ss = 0.;
      for (i in 1:B-1) {
        off_diagonal[i] = (-1.0 / (B - 1) - off_diag_ss) / diagon[i];
        off_diag_ss = off_diag_ss + off_diagonal[i] * off_diagonal[i];
        diagon[i+1] = sqrt(1.0 - off_diag_ss);
      }
      matrix[B, B-1] nu_multiplier = rep_matrix(0.0, B, B-1);
      for (i in 1:B-1) {
        nu_multiplier[i][i] = diagon[i];
        for (j in i:B-1) {
          nu_multiplier[j+1][i] = off_diagonal[i];
        }
      }
      vector[B] w = blosum_count[:B];
      real count_mean = mean(w);
      for (i in 1:B) {
        w[i] = count_mean / w[i];
      }
      matrix[B, B] weight_multiplier = diag_matrix(w);
      nu_multiplier = weight_multiplier * nu_multiplier;
    }
    parameters {
      vector[P] phi; // slope per position
      vector[B-1] nu_raw;
      vector<lower=0>[P] sigma2;
      vector<lower=0>[M] epsilon2;
      vector[V] eta2; // std_normal for beta
      vector[V] b; // intercept
    }
    transformed parameters {
      vector[B] nu = nu_multiplier * nu_raw;
      vector[V] beta; // slope per variants
      for (v in 1:V) {
        beta[v] = phi[vMAPp[v]] + nu[vMAPb[v]] + eta2[v] * sqrt(sigma2[vMAPp[v]]);
      }
    }
    model {
      phi ~ normal(0, 1);
      nu_raw ~ normal(0, 0.5);
      sigma2 ~ inv_gamma(1, 1);
      eta2 ~ normal(0, 1);
      epsilon2 ~ inv_gamma(1, 1);
      b ~ normal(0, 0.25);
      for (v in 1:V) {
        m[v] ~ normal(b[v] + beta[v] * t, sqrt(epsilon2[vMAPm[v]]));
      }
    }
    "
  } else if (type == "growth_pos") {
    stan_program <- "
    data {
      int<lower=0> T; // # of time points
      int<lower=0> V; // # of variants
      int<lower=0> P; // # of positions
      int<lower=0> M; // # of mean group
      array[V] int vMAPp;
      array[V] int vMAPm;
      vector[T] t; // time
      array[V] vector[T] m; // normalized count
    }
    parameters {
      vector[P] phi; // slope per position
      vector<lower=0>[P] sigma2; 
      vector<lower=0>[M] epsilon2;
      vector[V] eta2; // std_normal for beta
      vector[V] b; // intercept
    }
    transformed parameters {
      vector[V] beta; // slope per variants
      for (v in 1:V) {
        beta[v] = phi[vMAPp[v]] + eta2[v] * sqrt(sigma2[vMAPp[v]]);
      }
    }
    model {
      phi ~ normal(0, 1);
      sigma2 ~ inv_gamma(1, 1);
      eta2 ~ normal(0, 1);
      epsilon2 ~ inv_gamma(1, 1);
      b ~ normal(0, 0.25);
      for (v in 1:V) {
        m[v] ~ normal(b[v] + beta[v] * t, sqrt(epsilon2[vMAPm[v]]));
      }
    }
    "
  } else if (type == "growth_nopos") {
    stan_program <- "
    data {
      int<lower=0> T; // # of time points
      int<lower=0> V; // # of variants
      int<lower=0> M; // # of mean group
      array[V] int vMAPm;
      vector[T] t; // time
      array[V] vector[T] m; // normalized count
    }
    parameters {
      vector[V] beta;
      vector<lower=0>[M] epsilon2;
      vector[V] b; // intercept
    }
    model {
      beta ~ normal(0, 1);
      epsilon2 ~ inv_gamma(1, 1);
      b ~ normal(0, 0.25);
      for (v in 1:V) {
        m[v] ~ normal(b[v] + beta[v] * t, sqrt(epsilon2[vMAPm[v]]));
      }
    }
    "
  } else {
    stop("Invalid experiment type. Currently support 'growth' and 'binding'.")
  }

  file_pedantic <- cmdstanr::write_stan_file(code = stan_program,
                       dir = getOption("cmdstanr_write_stan_file_dir", tempdir()))

  return(file_pedantic)
}


