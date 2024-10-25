# Attempt to reproduce excellent work by Miller et al. 2020
# for use in Stan.

# Code from 
# https://github.com/dill/SPDE-smoothing

# aral sea example (2d smoothing)
# data from the gamair package
# original source http://seawifs.gsfc.nasa.gov/

## LIBRARIES 
library(mgcv)
source("https://raw.githubusercontent.com/dill/SPDE-smoothing/master/supplementary/mgcv_spde_smooth.R")
source("spde_smooth.R") 
library(fmesher)
library(INLA)
library(gamair) # package containing the data 
# for plotting:
library(ggplot2)
library(gridExtra)
# set seed
set.seed(25853)

## DATA
data(aral)
# boundary of observation window 
data(aral.bnd)

aral <- aral |> 
  dplyr::filter(!is.na(chl))
# split boundary into segments for INLA 
bnd <- fm_as_segm(sf::st_as_sf(tibble::tibble(X=aral.bnd$lon, Y=aral.bnd$lat),coords = c("X", "Y"),
                               crs = 4326))
loc <- cbind(aral$lon, aral$lat)

# Build a mesh within the boundary
mesh <- fm_mesh_2d_inla(boundary=bnd,
                        max.edge=c(0.384, 0.5),
                        min.angle=c(30, 21),
                        max.n=c(48000, 16000), ## Safeguard against large meshes.
                        max.n.strict=c(128000, 128000), ## Don't build a huge mesh!
                        cutoff=0.01, ## Filter away adjacent points.
                        offset=c(0.1, 0.3)) ## Offset for extra boundaries, if needed.
class(mesh) <- 'inla.mesh'
#### FIT SPDE MODEL WITH MGCV ################################################# 
mod <- gam(chl ~ s(lon, lat, bs = "spde", k = mesh$n, xt = list(mesh = mesh)),
           data = aral,
           control =  gam.control(scalePenalty = FALSE),
           method = "REML")
aral$pred_gam <- predict(mod)

g_test <- readr::read_rds("gam_testing_object.rds")
names(g_test)
(g_test$smooth[[1]]$bs.dim)
g_test$U1 |> heatm()


s_spde <- smooth.construct.spde.smooth.spec(s(lon, lat, bs = "spde", k = mesh$n, xt = list(mesh = mesh)),
    data = aral)

gratia::penalty(s_spde)


## Setup data for stan -----------------------------
smooth_spde <- smooth.construct.spde.test(aral,coords = c("lon", "lat"), mesh_in = mesh, knots = mesh$n )


# Add small values to matrices so that can combine them as sparse matricies in stan.
S <- smooth_spde$S
S[[1]][which((S[[3]]!=0)&S[[1]]==0)] <- .Machine$double.eps # set to non zero to allow all to be sparse
S[[2]][which((S[[3]]!=0)&S[[2]]==0)] <- .Machine$double.eps


stan_data <- list(
  n=nrow(aral),#   int n;    // n obs
  p =1, # int p;    // n par
  n_knots = nrow(smooth_spde$S[[1]]),# int row_spar;    // rows sparse matrix
  lat_lon = as.matrix(aral[,c("lon", "lat")] ),
  n_non_zero_M = sum(S[[1]]!=0), # int n_non_zero_M;
  n_non_zero_A = sum(smooth_spde$A != 0), # int n_non_zero_A;
  L = smooth_spde$L,
  y = aral$chl,# vector[n] y;         //The response
  X = matrix(rep(1, nrow(aral)), ncol =1),# matrix[n, p] X;         //Design matrix for fixed effects
  M0 = (S[[1]]),# matrix[row_spar, col_spar] M0;     // SPDE matrices from INLA
  M1 = (S[[2]]),# matrix[row_spar, col_spar] M1;
  M2 = (S[[3]]),# matrix[row_spar, col_spar] M2;
  A = (smooth_spde$A),# matrix[n, col_spar] A;     //Matrix for interpolating points witin triangles
  lambda = 1,
  prior_mean_tau_kappa_log = log(rev(c(3.603, 0.429))),
  prior_sd_tau_kappa_log = c(3,3)
)



# Build and Run Stan Model ----------------------------
library(cmdstanr)
mod_stan <- cmdstan_model("spde_sparse.stan")

samples <- mod_stan$sample(data = stan_data,
                           chains = 4,
                           parallel_chains = 4,
                           iter_warmup = 1000, 
                           iter_sampling = 1000,
                           output_dir = 
                             "F:/TMP_STAN")


sum_ <- samples$summary()

d <- samples$draws()

shinystan::launch_shinystan((samples))

stanfit <- rstan::read_stan_csv(list.files(
  "F:/TMP_STAN", full.names = T
))
shinystan::launch_shinystan(stanfit)

## Extract  values for comparison --------------

s_mle <- mod_stan$optimize(data = stan_data, seed = 123, jacobian = F,algorithm = 'newton',
                           output_dir = 
                             "F:/TMP_STAN")

ss <- s_mle$summary()
dplyr::filter(ss, stringr::str_detect(variable, "^u\\[")) |> 
  dplyr::pull(estimate) |> plot(y=mod$coefficients)

dd <- samples$draws(variables = glue::glue("eta[{1:485}]"))

 
aral$pred <- apply(dd, 3, mean) 
aral$pred_sd <- apply(dd, 3, sd) 
aral$ss_u <- dplyr::filter(ss, stringr::str_detect(variable, "^eta")) |> 
  dplyr::pull(estimate)

ggplot(aral, aes(chl, pred)) +
  geom_point() +
  geom_linerange(aes(ymin = pred-pred_sd, ymax = pred+pred_sd))+
  geom_point(aes(y = pred_gam),
             colour = 'red',
             alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0)# +
  geom_point(aes(y = ss_u), colour = 'grey', 
             alpha = 0.2)

ggplot(aral, aes(pred_gam, pred)) +
  geom_point(alpha = 0.3) +
  geom_linerange(aes(ymin = pred-pred_sd, ymax = pred+pred_sd),
                 alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0)




bayesplot::mcmc_trace(samples$draws(variables = c("tau_kappa_log","tau", "kappa",
                                                  "beta")))

bayesplot::mcmc_trace(samples$draws(variables = "log_lik_u") )
bayesplot::mcmc_trace(samples$draws(variables = "lambda") )
bayesplot::mcmc_pairs(samples$draws(
  variables = c("tau_kappa_log[1]", "tau_kappa_log[2]",
                "sigma", "beta")) )

bayesplot::mcmc_areas_ridges(samples$draws(variables = "delta"))


aral |> 
  # sf::st_as_sf(c('lon', 'lat'), crs = 4326) |> 
  ggplot(aes(lon, lat, fill = pred)) +
  geom_raster()

aral |> 
  ggplot(aes(lon, lat, fill = chl)) +
  geom_raster()





aral |> 
  ggplot(aes(lon, lat, fill = pred_gam)) +
  geom_raster()




dda <- samples$draws(variables = c("tau", "kappa", "beta[1]", 'range', 'sigma_spde', "sigma"))

## get estimates for mgcv
kappa <- mod$sp[2]
tau <- mod$sp[1]
# compute correlation range (rho) and marginal variance (sigma)
rho <- sqrt(8) / kappa
# see Lindgren et al. (2011) for this formula
sigma <- 1 / (sqrt(tau^2 * 4*pi * kappa^2))


cat("mgcv:\n")
cat("kappa=", kappa, "\n")
cat("tau=", tau, "\n\n")
cat("rho=", rho, "\n\n")
cat("sigma=", sigma, "\n\n")


cat("Stan averages")
apply(dda, 3, mean)
apply(dda, 3, sd)
