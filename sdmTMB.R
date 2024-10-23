library(ggplot2)
library(dplyr)
library(sdmTMB)
library(mgcv)
source("https://raw.githubusercontent.com/dill/SPDE-smoothing/master/supplementary/mgcv_spde_smooth.R")
source("spde_smooth.R") 
library(fmesher)
library(rstan) # for plot() method
options(mc.cores = parallel::detectCores()) # use rstan parallel processing

library(cmdstanr)


set.seed(123)
predictor_dat <- data.frame(
  X = runif(500), Y = runif(500),
  a1 = rnorm(500)
)

mesh_INLA <- fm_mesh_2d_inla(predictor_dat[,c("X", "Y")], cutoff = 0.1)
mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1, mesh = mesh_INLA)
# plot(mesh)
# mesh$mesh$n
sim_dat <- sdmTMB_simulate(
  formula = ~a1,
  data = predictor_dat,
  mesh = mesh,
  family = gaussian(),
  range = 0.3,
  phi = 0.2,
  sigma_O = 0.2,
  seed = 123,
  B = c(0.8, -0.4) # B0 = intercept, B1 = a1 slope
)

ggplot(sim_dat, aes(X, Y, colour = observed)) +
  geom_point() +
  scale_color_viridis_c()

fit <- sdmTMB(
  observed ~ a1,
  data = sim_dat,
  mesh = mesh,
  family = gaussian(),
  spatial = "on"
)
fit

fit <- sdmTMB(
  observed ~ a1,
  data = sim_dat,
  mesh = mesh,
  family = gaussian(),
  spatial = "on",
  priors = sdmTMBpriors(
    # location = vector of means; scale = vector of standard deviations:
    b = normal(location = c(0, 0), scale = c(5, 2)),
  )
)
fit


# grab the internal parameter list at estimated values:
pars <- sdmTMB::get_pars(fit)
# create a 'map' vector for TMB
# factor NA values cause TMB to fix or map the parameter at the starting value:
kappa_map <- factor(rep(NA, length(pars$ln_kappa)))

# rebuild model updating some elements:
fit_mle <- update(
  fit,
  control = sdmTMBcontrol(
    start = list(
      ln_kappa = pars$ln_kappa #<
    ),
    map = list(
      ln_kappa = kappa_map #<
    )
  ),
  do_fit = FALSE #<
)



fit_stan <- tmbstan::tmbstan(
  fit_mle$tmb_obj,
  iter = 1000, chains = 2,
  seed = 8217 # ensures repeatability
)


plot(fit_stan)
pars_plot <- c("b_j[1]", "b_j[2]", "ln_tau_O", "omega_s[1]")

bayesplot::mcmc_trace(fit_stan, pars = pars_plot)
bayesplot::mcmc_pairs(fit_stan, pars = pars_plot)


## My model

## Setup data for stan -----------------------------
smooth_spde <- smooth.construct.spde.test(predictor_dat,
                                          coords = c("X", "Y"),
                                          mesh_in = mesh_INLA, 
                                          knots = mesh_INLA$n )


# Add small values to matrices so that can combine them as sparse matricies in stan.
S <- smooth_spde$S
S[[1]][which((S[[3]]!=0)&S[[1]]==0)] <- .Machine$double.eps # set to non zero to allow all to be sparse
S[[2]][which((S[[3]]!=0)&S[[2]]==0)] <- .Machine$double.eps


stan_data <- list(
  n=nrow(sim_dat),#   int n;    // n obs
  p =2, # int p;    // n par
  n_knots = nrow(smooth_spde$S[[1]]),# int row_spar;    // rows sparse matrix
  lat_lon = as.matrix(sim_dat[,c("X", "Y")] ),
  n_non_zero_M = sum(S[[1]]!=0), # int n_non_zero_M;
  n_non_zero_A = sum(smooth_spde$A != 0), # int n_non_zero_A;
  L = smooth_spde$L,
  y = sim_dat$observed,# vector[n] y;         //The response
  X = model.matrix(fit_mle),# matrix[n, p] X;         //Design matrix for fixed effects
  M0 = (S[[1]]),# matrix[row_spar, col_spar] M0;     // SPDE matrices from INLA
  M1 = (S[[2]]),# matrix[row_spar, col_spar] M1;
  M2 = (S[[3]]),# matrix[row_spar, col_spar] M2;
  A = (smooth_spde$A),# matrix[n, col_spar] A;     //Matrix for interpolating points witin triangles
  lambda = 1,
  prior_mean_tau_kappa_log = c(0,1.858606),#log(c(3.603, 0.429)),
  prior_sd_tau_kappa_log = c(3,0.001)
)



# Build and Run Stan Model ----------------------------
library(cmdstanr)
mod_stan <- cmdstan_model("spde_sparse.stan")

samples <- mod_stan$sample(data = stan_data,
                           chains = 2,seed = 8217,
                           parallel_chains = 2,
                           iter_warmup = 500, 
                           iter_sampling = 500)
pars_plot2 <- c("beta[1]", "beta[2]", "tau_kappa_log[1]", "sigma_spde")
dp <- samples$draws(variables = pars_plot2)
bayesplot::mcmc_trace(dp)
bayesplot::mcmc_pairs(dp)



mesh_INLA <- fm_mesh_2d(pcod[,c("X", "Y")], cutoff = 10)
mesh <- make_mesh(
  pcod,
  xy_cols = c("X", "Y"),
  mesh = mesh_INLA
  
)
plot(mesh)
plot(mesh_INLA)


m <- sdmTMB(
  data = pcod,
  formula = present ~ depth_scaled + depth_scaled2,
  mesh = mesh,
  family = binomial(link = "logit"),
  spatial = "on"
)
m
AIC(m)


## Setup data for stan -----------------------------
smooth_spde <- smooth.construct.spde.test(pcod,
                                          coords = c("X", "Y"),
                                          mesh_in = mesh_INLA, 
                                          knots = mesh_INLA$n )


# Add small values to matrices so that can combine them as sparse matricies in stan.
S <- smooth_spde$S
S[[1]][which((S[[3]]!=0)&S[[1]]==0)] <- .Machine$double.eps # set to non zero to allow all to be sparse
S[[2]][which((S[[3]]!=0)&S[[2]]==0)] <- .Machine$double.eps


stan_data <- list(
  n=nrow(pcod),#   int n;    // n obs
  p =3, # int p;    // n par
  n_knots = nrow(smooth_spde$S[[1]]),# int row_spar;    // rows sparse matrix
  lat_lon = as.matrix(pcod[,c("X", "Y")] ),
  n_non_zero_M = sum(S[[1]]!=0), # int n_non_zero_M;
  n_non_zero_A = sum(smooth_spde$A != 0), # int n_non_zero_A;
  L = smooth_spde$L,
  y = pcod$present,# vector[n] y;         //The response
  X = model.matrix(m),# matrix[n, p] X;         //Design matrix for fixed effects
  M0 = (S[[1]]),# matrix[row_spar, col_spar] M0;     // SPDE matrices from INLA
  M1 = (S[[2]]),# matrix[row_spar, col_spar] M1;
  M2 = (S[[3]]),# matrix[row_spar, col_spar] M2;
  A = (smooth_spde$A),# matrix[n, col_spar] A;     //Matrix for interpolating points witin triangles
  lambda = 1
)

brms <- 
brms::brm( data = pcod,
           formula = present ~ depth_scaled + depth_scaled2 + s(X, Y, bs = 'tp'), family = "bernoulli", fit = F)

# Build and Run Stan Model ----------------------------
library(cmdstanr)
mod_stan <- cmdstan_model("spde_binom.stan")

samples <- mod_stan$sample(data = stan_data,
                           chains = 2,
                           parallel_chains = 2,
                           iter_warmup = 1000, 
                           iter_sampling = 1000)

dd <- samples$draws(variables = glue::glue("eta[{1:nrow(pcod)}]"))
# dda <- samples$draws(variables = c("tau", "kappa", "beta[1]", 'range', 'sigma'))

pcod$pred <- apply(dd, 3, mean) 
pcod$pred_sdmTMB <- predict(m)$est
ss <- samples$summary()

ss |> filter(stringr::str_detect(variable, "beta"))


tidy(m, conf.int = TRUE)

ss |> filter(stringr::str_detect(variable, "range|sigma_spde"))

tidy(m, "ran_pars", conf.int = TRUE)


ggplot(pcod, aes(pred_sdmTMB, pred)) + geom_point()



bayesplot::mcmc_trace(samples$draws(variables = c("tau_kappa_log","tau", "kappa",
                                                  "beta")))


set.seed(123)
samps <- sdmTMBextra::predict_mle_mcmc(m, mcmc_warmup = 100, mcmc_iter = 101)


pcod$mcmc <- samps

ggplot(pcod, aes(mcmc, pred)) + geom_point()

