
library(simr)
library(lme4)
library(tidyverse)

#### Helper functions ####

# Test multiple values for the numbers of clusters
power_extend <- function(fit, along, within, new_size, time_points, ...) {
  if (!missing(time_points)) {
    power_lst <- vector("list", length(time_points))
    for (j in seq_along(time_points)) {
      sim_m <- extend(fit, along = along, within = within, 
                      values = (1:time_points[j]) - 1)
      power_lst[[j]] <- powerSim(sim_m, ...)
    }
    setNames(power_lst, time_points)
  } else {
    power_lst <- vector("list", length(new_size))
    for (j in seq_along(new_size)) {
      sim_m <- extend(fit, along = along, within = within, n = new_size[j])
      power_lst[[j]] <- powerSim(sim_m, ...)
    }
    setNames(power_lst, new_size)
  }
}
# Compute proportions and CI
powerinterval <- function(object, alpha = object$alpha, level = 0.95, 
                          method = getSimrOption("binom"), ...) {
  x <- sum(object$pval < alpha, na.rm = TRUE)
  n <- object$n
  if (is.na(x) || is.na(n) || (n == 0)) {
    return(rep(NA, 3))
  }
  binom::binom.confint(x, n, level, method, ...)[c("mean", "lower", "upper")]
}
# Add a plotting function to get power curve
plot_power_extend <- function(object, target = .80, ...) {
  nsim <- object[[1]]$n
  df <- purrr::map_dfr(object, powerinterval, .id = "n", ...)
  df$n <- as.integer(df$n)
  ggplot(df, aes(x = n, y = mean, group = 1)) + 
    geom_hline(yintercept = target) + 
    geom_smooth(method = "glm", 
                formula = cbind(y * nsim, (1 - y) * nsim) ~ x, 
                method.args = list(family = binomial("probit"))) + 
    geom_pointrange(aes(ymin = lower, ymax = upper)) + 
    labs(y = "power") + 
    ylim(0, 1)
}


#### Model 1 ####
# with random slope for the users' habit strength

# Cluster ID (start with 10 clusters)
num_clus <- 10
id <- 1:num_clus
# Continuous lv-2 predictor
habit <- rnorm(num_clus)
# Force SD = 1
habit <- (habit - mean(habit)) / sd(habit) * (num_clus - 1) / num_clus
# Level-2 data
lv2_dat <- tibble(id, habit = habit)
# Cluster means of X, with ICC(X) = 0.2
nlike_cm <- rnorm(num_clus, mean = 0, sd = sqrt(0.2))
# Force SD = sqrt(0.2)
nlike_cm <- (nlike_cm - mean(nlike_cm)) / sd(nlike_cm) *
  sqrt(0.2) * (num_clus - 1) / num_clus
# Level-2 data
lv2_dat <- tibble(id, nlike_cm = nlike_cm, habit = habit)
# Expand each cluster to include more rows
clus_size <- 15 # Cluaster size
lv2_dat <- lv2_dat %>%
  slice(rep(1:n(), each = clus_size))
# Within-cluster component of X, with sigma^2(X) = 0.8
num_obs <- num_clus * clus_size
nlike_cmc <- rnorm(num_obs, mean = 0, sd = 1)
nlike_cmc <- nlike_cmc - ave(nlike_cmc, lv2_dat$id, FUN = mean)
# Force SD = sqrt(0.8)
nlike_cmc <- nlike_cmc / sd(nlike_cmc) * sqrt(0.8) * (num_obs - 1) / num_obs
# Expand each cluster to include more rows
(sim_dat <- mutate(lv2_dat, nlike_cmc = nlike_cmc))
#dataset-centered likes as well...
N_x <- length(id) * nj
x <- rnorm(N_x)
##standardize x:
like <- (x - mean(x)) / sd(x) * (N_x - 1) / N_x
sim_dat$like <- like
#### SIMULATION

gams2 <- c(0, 0.3, 0.3, 0.2)  # gamma_00, gamma_01, gamma_10 #fixed 
taus2 <- matrix(c(0.25, 0,
                         0, 0.10), nrow = 2)

sim_m2 <- makeLmer(y ~ habit + like + like:habit +
                     (like | id),
                   fixef = gams2, VarCorr = taus2, sigma = sigma,
                   data = sim_dat)

powers_m2 <- power_extend(sim_m2, along = "id", 
                          new_size = c(10, 20, 50, 100), 
                          nsim = 2000, 
                          test = simr::fixed("habit:like", 
                                             method = "anova"))



plot_power_extend(powers_m2) #POWER PLOT FOR SUPPLEMENT

r.squaredGLMM(sim_m2) #look at marginal rsquared

# gamma_00, gamma_01, gamma_02, gamma_10, gamma_20
gams <- c(0, 0.3, 0.2)
taus <- matrix(c(0.25, 0,
                 0, 0.10), nrow = 2)  # tau^2_0 = 0.25, tau^2_1 = 0.1
sigma <- 1
sim_m0 <- makeLmer(y ~ habit + like  +
                     (like | id),
                   fixef = gams, VarCorr = taus, sigma = sigma,
                   data = sim_dat)
r.squaredGLMM(sim_m0)
(r.squaredGLMM(sim_m2) - r.squaredGLMM(sim_m0)) / (1 - (r.squaredGLMM(sim_m2))) ## THE F2 ESTIMATE


#################################################################
###################NOT USED IN THE PAPER#########################
#################################################################
### FOR PARTICIPANT_MEANCENTERED
# gamma_00, gamma_01, gamma_02, gamma_10, gamma_20
gams <- c(0, 0.3, 0.3, 0.3, 0.2)
taus <- matrix(c(0.25, 0,
                 0, 0.10), nrow = 2)  # tau^2_0 = 0.25, tau^2_1 = 0.1
sigma <- 1

sim_m1 <- makeLmer(y ~ habit + nlike_cm + nlike_cmc + nlike_cmc:habit +
                     (nlike_cmc | id),
                   fixef = gams, VarCorr = taus, sigma = sigma,
                   data = sim_dat)

powers_m1 <- power_extend(sim_m1, along = "id", 
                          new_size = c(10, 20, 50, 100), 
                          nsim = 2000, 
                          test = simr::fixed("habit:nlike_cmc", 
                                             method = "anova"))

