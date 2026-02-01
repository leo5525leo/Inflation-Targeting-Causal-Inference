# Load packages 
library(tidyverse)
library(zoo)
library(gsynth)
library(dplyr)
library(ggplot2)

# Load data
df <- read.csv("inflation_targeting_simulated.csv")

# # Prepare data. Standardize variable types (factor, integer) and sort panel by country-year
df <- df %>%
  mutate(
    country = as.factor(country),
    year    = as.integer(year),
    it.policy = as.integer(it.policy),
    exchange.rate = factor(exchange.rate)  # Fixed / Floating
  ) %>%
  arrange(country, year)


# Compute 5-year rolling standard deviation of annual inflation (volatility measure)
df <- df %>%
  group_by(country) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    infl_vol_5y = rollapply(
      inflation, width = 5,
      FUN = sd, align = "right",
      fill = NA, na.rm = TRUE
    )
  ) %>%
  ungroup()


# Clean panel 
panel_df <- df %>%
  filter(year >= 1984) %>%     # remove rows that cannot be used
  arrange(country, year)


# Step 6: Run the GSCM (NO COVARIATES)
res_gscm <- gsynth(
  Y = "infl_vol_5y",
  D = "it.policy",
  data = panel_df,
  index = c("country", "year"),
  force = "unit",                         # country FE + year FE
  CV = TRUE,                                 # cross-validation to pick # of factors
  r = c(0, 5),                               # search range for latent factors
  se = TRUE,
  inference = "parametric",
  nboots = 500,
  min.T0 = 7.                                # require ≥7 pre-treatment years
)

summary(res_gscm)

# Event-study style ATT plot (treatment vs synthetic)
plot(res_gscm, "gap")

# Plot actual treated average vs estimated counterfactual average
plot(res_gscm, "counterfactual")

res_gscm$att.avg
#-0.8330652

est_avg <- res_gscm$est.avg; est_avg
att_hat   <- est_avg["ATT.avg", "Estimate"];att_hat
ci_lower  <- est_avg["ATT.avg", "CI.lower"];ci_lower
ci_upper  <- est_avg["ATT.avg", "CI.upper"];ci_upper




# Robustness
# ============================
# PLACEBO TEST (UNIT TEST)
# ============================
# Purpose: Randomly assign treatment to units that were never treated
# This helps check for false positives
# First run placebo for one sample.
# Second, build a placebo distribution. Difference combination of treated countries

# Placebo one time
# 1. Identify NEVER-treated countries
never_treated <- panel_df %>%
  group_by(country) %>%
  summarize(ever_treated = max(it.policy)) %>%
  filter(ever_treated == 0) %>%
  pull(country)

length(never_treated)  # should be 30

# 2. Randomly pick 50% for placebo treatment
set.seed(123)   # reproducible
placebo_treated_units <- sample(never_treated, length(never_treated)/2)

# 3. Construct placebo dataset
placebo_df <- panel_df %>%
  filter(country %in% never_treated) %>%   # Only never-treated countries
  mutate(
    fake_treat = ifelse(country %in% placebo_treated_units & year >= 2005, 1, 0)
  )

# Check counts
table(placebo_df$fake_treat, placebo_df$year == 2005)[,1]
length(unique(placebo_df$country[placebo_df$fake_treat == 1]))
length(unique(placebo_df$country[placebo_df$fake_treat == 0]))

# 4. Run GSCM with placebo treatment
res_placebo <- gsynth(
  Y = "infl_vol_5y",
  D = "fake_treat",
  data = placebo_df,
  index = c("country", "year"),
  force = "unit",
  CV = TRUE,
  r = c(0, 5),
  se = TRUE,
  inference = "parametric",
  nboots = 500,
  min.T0 = 7
)

# 5. Plot placebo results
res_placebo$att.avg
plot(res_placebo, type = "gap")            # ATT graph (should be ~0)
plot(res_placebo, type = "counterfactual") # actual vs synthetic

# This placebo test increases confidence in your GSCM model’s credibility. 
# It shows that when treatment is assigned randomly, 
# the estimator does not find large or systematic treatment effects.



# Distribution of placebo
# Create an empty container to store placebo ATT estimates
n_sim <- 300
placebo_att <- numeric(n_sim)

for (i in 1:n_sim) {
  # Randomly assign placebo treatment to 50% of never-treated units
  placebo_units <- sample(never_treated, length(never_treated)/2)
  
  # Construct the placebo dataset
  placebo_df <- panel_df %>%
    filter(country %in% never_treated) %>%  # Keep only never-treated units
    mutate(fake_treat = ifelse(country %in% placebo_units & year >= 2005, 1, 0))  # Assign placebo treatment from 2005
  
  # Run GSCM with placebo treatment
  res <- gsynth(
    Y = "infl_vol_5y",
    D = "fake_treat",
    data = placebo_df,
    index = c("country", "year"),
    force = "unit",  # or "two-way"
    CV = TRUE,
    r = c(0, 5),
    se = FALSE,
    min.T0 = 7
  )
  
  # Store the average ATT for this simulation
  placebo_att[i] <- mean(res$att.avg, na.rm = TRUE)
}

# Plot the distribution of placebo ATTs
# Choose x-limits that include both placebo draws and the real ATT
real_att <- -0.83  # your estimated ATT
x_min <- min(placebo_att, real_att) - 0.05
x_max <- max(placebo_att, real_att) + 0.05
hist(placebo_att,
     main = "Placebo Distribution of ATT",
     xlab = "Estimated ATT",
     col  = "lightgray",
     border = "white",
     xlim = c(x_min, x_max))
# Null line at 0
abline(v = 0, col = "red", lwd = 2)
# Mean placebo line
abline(v = mean(placebo_att), col = "blue", lwd = 2, lty = 2)
# Real ATT line
abline(v = real_att, col = "darkgreen", lwd = 2)
legend("topright",
       legend = c("Null (0)", "Placebo mean", "Real ATT"),
       col    = c("red", "blue", "darkgreen"),
       lty    = c(1, 2, 1),
       lwd    = 2,
       bty    = "n")




# ------------------------------------------------------------------------
# 5) Bias and variance (20 points):
# ------------------------------------------------------------------------

# ---------------------------
# 1. REAL ATT for volatility 
# ---------------------------

# Check file
truth <- read.csv("inflation_targeting_instructor.csv")
## 0 -  General Check
treated <- truth[truth$it.policy == 1, ]
# Number of treated observation 
n_treated <- nrow(treated)
# Number of treated and control in the Instructor 
table_by_country <- aggregate(it.policy ~ country, data = truth, FUN = max)
table(table_by_country$it.policy)
# Check the variable te.mult (Treatment effect multiplier)
# Correlation check
cor(truth$ite_inflation, truth$te.mult, use = "complete.obs")
# Simple ratio
subset <- truth[truth$te.mult != 0, ]
head(subset$ite_inflation / subset$te.mult)


# Rolling inflation
# Group by country and calculate 5-year rolling SD for both potential outcomes
truth <- truth %>%group_by(country) %>%arrange(year) %>%
  mutate(
    roll_vol_y1 = rollapply(inflation_y1, width = 5, FUN = sd, fill = NA, align = "right"),
    roll_vol_y0 = rollapply(inflation_y0, width = 5, FUN = sd, fill = NA, align = "right"),
    ite_rolling_vol = roll_vol_y1 - roll_vol_y0)
# Keep only treated observations
treated_roll <- truth %>% filter(it.policy == 1)
# Step 3: Compute true ATT based on rolling volatility
true_ATT_rolling_vol <- mean(treated_roll$ite_rolling_vol, na.rm = TRUE)
true_ATT_rolling_vol
# -0.3249706
# Use rolling inflation to match with define outcome in the project


# Step 5B – Compare estimated ATT vs true ATT
# (1) Estimated ATT from gsynth result
estimated_ATT <- mean(res_gscm$att.avg, na.rm = TRUE)
estimated_ATT
# (2) Bias
bias <- estimated_ATT - true_ATT_rolling_vol
bias
#-0.5080946

# Our estimated ATT of –0.83 under the GSCM model falls short of the true ATT of –0.32. 
# This suggests that, in this draw, the estimator underestimates the full effect of inflation targeting on reducing inflation volatility. 
# The bias of -0.5 may reflect either random sampling error or systematic bias, which we further investigate in the next steps.


# ---------------------------
# 1. Monte Carlo Using true DGP 
# (Changing the seed) Simple model NO COVARIATE
# ---------------------------

# === Set number of simulations ===
B <- 500  # You can increase to 500 later
att_estimates <- numeric(B)
covered       <- logical(B)   # TRUE/FALSE if CI contains true ATT
true_ATT      <- -0.325       # whatever your DGP’s true ATT is

# === Loop starts here ===
for (b in 1:B) {
  cat("Simulation", b, "\n")
  
  # Step 1: Generate a random seed and inject into the simulation script
  original_lines <- readLines("simulate_inflation_targeting_panel_LeeOcampos_export.R")
  new_seed <- sample(1:99999999, 1)
  modified_lines <- gsub("SEED <- .*", paste0("SEED <- ", new_seed), original_lines)
  writeLines(modified_lines, "temp_simulation.R")
  
  # Step 2: Run the modified simulation file
  source("temp_simulation.R")
  
  # Step 3: Load the new simulated dataset
  df <- as.data.frame(student_out)  # already created in the simulation script
  
  # Step 4: Compute 5-year rolling volatility
  df <- df %>%
    group_by(country) %>%
    arrange(year) %>%
    mutate(infl_vol = rollapply(inflation, width = 5, FUN = sd, fill = NA, align = "right")) %>%
    ungroup()%>%
    filter(year >= 1984)
  
  # Step 5: Run GSCM
  tryCatch({
    res_gscm <- gsynth(
      Y = "infl_vol",
      D = "it.policy",
      data = df,
      index = c("country", "year"),
      force = "unit",  # or "unit", adjust as needed
      se = TRUE,
      inference = "parametric",  
      CV = TRUE,
      min.T0 = 7 
    )
    
    est_avg  <- res_gscm$est.avg
    att_hat  <- est_avg["ATT.avg", "Estimate"]
    ci_lower <- est_avg["ATT.avg", "CI.lower"]
    ci_upper <- est_avg["ATT.avg", "CI.upper"]
    
    
    # Step 6: Store estimated ATT
    att_estimates[b] <- att_hat
    covered[b]       <- (true_ATT >= ci_lower & true_ATT <= ci_upper)
    
  }, error = function(e) {
    cat("Error in iteration", b, ":", conditionMessage(e), "\n")
    att_estimates[b] <- NA
    covered[b]       <- NA
  })
}

# === Analyze results ===
att_estimates_clean <- att_estimates[!is.na(att_estimates)]
covered_clean       <- covered[!is.na(covered)]

# Summary stats
mean_ATT <- mean(att_estimates_clean)
bias <- mean_ATT - (-0.325)  # Adjust to true ATT you’re using
variance <- var(att_estimates_clean)
rmse <- sqrt(mean((att_estimates_clean - (-0.325))^2))
coverage  <- mean(covered_clean)  # THIS is your 95% CI coverage

cat("\nSummary of Monte Carlo ATT estimates:\n")
cat("Mean ATT: ", mean_ATT, "\n")
cat("Bias:     ", bias, "\n")
cat("Variance: ", variance, "\n")
cat("RMSE:     ", rmse, "\n")
cat("95% CI coverage: ", coverage, "\n")


# === Plot the sampling distribution ===
hist(att_estimates_clean, breaks = 30, col = "skyblue", main = "Sampling Distribution of GSCM ATT",
     xlab = "Estimated ATT")
abline(v = -0.833, col = "purple", lwd = 2, lty = 2)  # Your original estimate
abline(v = -0.325, col = "red", lwd = 2, lty = 2)     # True ATT from DGP
abline(v = mean_ATT, col = "darkgreen", lwd = 2)      # Monte Carlo mean estimate
legend("topright", legend = c("First Estimate", "True ATT", "Mean Estimate"),
       col = c("purple", "red", "darkgreen"), lty = c(2, 2, 1), lwd = 2)

lower <- quantile(att_estimates_clean, 0.025)
upper <- quantile(att_estimates_clean, 0.975)
att_trimmed <- att_estimates_clean[att_estimates_clean >= lower & att_estimates_clean <= upper]

hist(att_trimmed, breaks = 30, col = "skyblue", main = "Sampling Distribution of GSCM ATT",
     xlab = "Estimated ATT")
abline(v = -0.833, col = "purple", lwd = 2, lty = 2)  #  Original estimate
abline(v = -0.325, col = "red", lwd = 2, lty = 2)     # True ATT from DGP
abline(v = mean_ATT, col = "darkgreen", lwd = 2)      # Monte Carlo mean estimate
legend("topright", legend = c("First Estimate", "True ATT", "Mean Estimate"),
       col = c("purple", "red", "darkgreen"), lty = c(2, 2, 1), lwd = 2)

# Bias and Variance Analysis 

# Bias and Variance Evaluation # Another option for interpretation (number and be difference When running the Monte Carlo)
# Our initial GSCM estimate of the ATT on 5-year rolling inflation volatility was –0.833, while the true ATT in the data-generating process (DGP) is –0.325. To understand whether this large effect was due to sampling 
# variation or systematic issues, we conducted 500 Monte Carlo simulations using new data draws from the same DGP. The average ATT estimate across these draws was –0.757, very close to our original estimate, 
# indicating that our initial result was not an outlier. However, this also revealed that the GSCM estimator is substantially biased downward in this setting, with a bias of –0.432 and an RMSE of 1.125. 
# Therefore, the discrepancy between our estimate and the true ATT was not due to sampling variance or bad luck, but rather due to systematic bias in the estimator under this DGP.


# ------------------------------------------------------------------------
# Part 6- Improve the model o change DGP:
# ------------------------------------------------------------------------

# Now the know what our model is bias, how to reduce this bias? two options
# 1 - Change estimator specification: Try model with covariate 
# 2 - Change the DGP
# Use Monte Carlo to evaluate


# -----------------------------------------------
# Model WITH COVARIATES
# -----------------------------------------------

# Adjust exchange rate
panel_df <- panel_df %>%
  mutate(exchange.rate.num = ifelse(exchange.rate == "Floating", 1, 0))

res_gscm_cov <- gsynth(
  Y = "infl_vol_5y",
  D = "it.policy",
  X = c("cbi", "trade.open"),
  data = panel_df,
  index = c("country", "year"),
  force = "unit",
  CV = TRUE,
  r = c(0, 5),
  se = TRUE,
  inference = "parametric",
  nboots = 500,
  min.T0 = 7
)

# Governance AND prior inflation and exchange rate is already capture in the unit FE
summary(res_gscm_cov)

# Event-study style ATT plot (treatment vs synthetic)
plot(res_gscm_cov, "gap")

# Plot actual treated average vs estimated counterfactual average
plot(res_gscm_cov, "counterfactual")

res_gscm_cov$att.avg
# -0.5854747

# This number is closer to true ATT=-0.3
# Now also run placebo to check robustness

# -----------------------------------------------
# Placebo Test with Covariates
# -----------------------------------------------

set.seed(12345)  # for reproducibility

# Step 1: Keep only never-treated countries from the original dataset
never_treated_countries <- panel_df %>%
  group_by(country) %>%
  summarize(ever_treated = max(it.policy)) %>%
  filter(ever_treated == 0) %>%
  pull(country)

placebo_df <- panel_df %>%
  filter(country %in% never_treated_countries)

# Step 2: Randomly assign 50% of these to placebo treatment starting in 2005
placebo_treated <- sample(unique(placebo_df$country), size = floor(length(unique(placebo_df$country)) / 2))

placebo_df <- placebo_df %>%
  mutate(
    fake_treat = ifelse(country %in% placebo_treated & year >= 2005, 1, 0)
  )

# Step 3: Run GSCM with placebo treatment and covariates
res_placebo_cov <- gsynth(
  Y = "infl_vol_5y", 
  D = "fake_treat",
  X = c("cbi", "trade.open"),
  data = placebo_df,
  index = c("country", "year"),
  force = "unit",
  CV = TRUE,
  r = c(0, 5),
  se = TRUE,
  inference = "parametric",
  nboots = 500,
  min.T0 = 7
)

# Step 4: Plot results
plot(res_placebo_cov, type = "gap")  # Event-study style
plot(res_placebo_cov, type = "counterfactual")  # Actual vs counterfactual

# Step 5: Inspect average estimated placebo effect
res_placebo_cov$att.avg

# To assess the risk of spurious findings, we ran a placebo test using only never-treated countries and randomly assigned a false inflation targeting treatment to half of them starting in 2005. 
# Under the GSCM model with covariates (central bank independence and trade openness), the estimated average treatment effect was just 0.027, close to zero. This indicates that our model does not detect artificial effects when no true policy change occurs — confirming that the earlier negative ATT estimates are not driven by overfitting or model artifacts.





# now also rum Monte Carlo with covariate to check bias
# ---------------------------
# 1. Monte Carlo with COVARIATE
# ---------------------------

B <- 500 # Number of simulations
att_estimates <- numeric(B)
num_failures <- 0  # Track errors

# === Monte Carlo Loop ===
for (b in 1:B) {
  cat("Simulation", b, "\n")
  
  # Step 1: Generate new seed and overwrite SEED in sim script
  original_lines <- readLines("simulate_inflation_targeting_panel_LeeOcampos_export.R")
  new_seed <- sample(1:99999999, 1)
  modified_lines <- gsub("SEED <- .*", paste0("SEED <- ", new_seed), original_lines)
  writeLines(modified_lines, "temp_simulation.R")
  
  # Step 2: Run simulation
  source("temp_simulation.R")
  
  # Step 3: Load data
  df <- as.data.frame(student_out)
  
  # Step 4: Rolling 5-year inflation volatility
  df <- df %>%
    group_by(country) %>%
    arrange(year) %>%
    mutate(infl_vol = rollapply(inflation, width = 5, FUN = sd, fill = NA, align = "right")) %>%
    ungroup() %>%
    filter(year >= 1984)
  
  # Step 5: GSCM with covariates
  tryCatch({
    res_gscm <- gsynth(
      Y = "infl_vol",
      D = "it.policy",
      X = c("cbi", "trade.open"),
      data = df,
      index = c("country", "year"),
      force = "two-way",
      CV = TRUE,
      r = c(0, 5),
      se = FALSE,
      min.T0 = 7
    )
    
    att_estimates[b] <- mean(res_gscm$att.avg, na.rm = TRUE)
    
  }, error = function(e) {
    cat("Error in iteration", b, ":", conditionMessage(e), "\n")
    att_estimates[b] <- NA
    num_failures <<- num_failures + 1
  })
}

# === Analyze results ===
att_clean <- na.omit(att_estimates)

true_ATT <- -0.325  # Adjust if needed
mean_ATT <- mean(att_clean)
bias <- mean_ATT - true_ATT
variance <- var(att_clean)
rmse <- sqrt(mean((att_clean - true_ATT)^2))

cat("\nMonte Carlo Summary:\n")
cat("Valid estimates: ", length(att_clean), " / ", B, "\n")
cat("Mean ATT:        ", mean_ATT, "\n")
cat("Bias:            ", bias, "\n")
cat("Variance:        ", variance, "\n")
cat("RMSE:            ", rmse, "\n")

# === Plot full distribution ===
hist(att_clean, breaks = 30, col = "skyblue",
     main = "Sampling Distribution of GSCM ATT (Two-way FE-With Covariates)",
     xlab = "Estimated ATT (Volatility)")
abline(v = -0.833, col = "purple", lwd = 2, lty = 2)
abline(v = true_ATT, col = "red", lwd = 2, lty = 2)
abline(v = mean_ATT, col = "darkgreen", lwd = 2)
legend("topright", legend = c("Your Estimate", "True ATT", "Mean Estimate"),
       col = c("purple", "red", "darkgreen"), lty = c(2, 2, 1), lwd = 2)

# === Plot trimmed (central 95%) distribution ===
lower <- quantile(att_clean, 0.025)
upper <- quantile(att_clean, 0.975)
att_trimmed <- att_clean[att_clean >= lower & att_clean <= upper]

hist(att_trimmed, breaks = 30, col = "skyblue",
     main = "Trimmed GSCM ATT Distribution (Central 95%)",
     xlab = "Estimated ATT (Volatility)")
abline(v = -0.833, col = "purple", lwd = 2, lty = 2)
abline(v = true_ATT, col = "red", lwd = 2, lty = 2)
abline(v = mean_ATT, col = "darkgreen", lwd = 2)
legend("topright", legend = c("First Estimate", "True ATT", "Mean Estimate"),
       col = c("purple", "red", "darkgreen"), lty = c(2, 2, 1), lwd = 2)


# The model with covariate is less bias that the simple model. Mean ATT = -0.681. Compare wiht the simple model (0.72)
# Rember true ATT = -0.32





# ---------------------------
# Change the DGP
# ---------------------------

# Now change the DGP and calculate Monte Carlo with this new code 
# ================================
# Monte Carlo Simulation with EXOGENOUS Adoption
# ================================

# Parameters
B <- 3  # number of Monte Carlo draws
att_estimates <- numeric(B)

for (b in 1:B) {
  cat("Simulation", b, "\n")
  
  # Step 1: Load and modify DGP to make adoption exogenous
  original_lines <- readLines("simulate_inflation_targeting_panel_LeeOcampos_export.R")
  new_seed <- sample(1:99999999, 1)
  modified_lines <- gsub("SEED <- .*", paste0("SEED <- ", new_seed), original_lines)
  modified_lines <- gsub("use_lag_infl\\s*=\\s*TRUE", "use_lag_infl = FALSE", modified_lines)
  writeLines(modified_lines, "temp_sim_exog.R")
  source("temp_sim_exog.R")  # this creates `student_out`
  
  # Step 2: Prepare data and compute rolling volatility
  df <- as.data.frame(student_out) %>%
    group_by(country) %>%
    arrange(year) %>%
    mutate(infl_vol = rollapply(inflation, width = 5, FUN = sd, fill = NA, align = "right")) %>%
    ungroup() %>%
    filter(year >= 1984)
  
  # Step 3: Estimate GSCM with covariates
  tryCatch({
    res_gscm <- gsynth(
      Y = "infl_vol",
      D = "it.policy",
      X = c("cbi", "trade.open"),
      data = df,
      index = c("country", "year"),
      force = "unit",
      CV = TRUE,
      r = c(0, 5),
      se = FALSE,
      min.T0 = 7
    )
    
    att_estimates[b] <- mean(res_gscm$att.avg, na.rm = TRUE)
    
  }, error = function(e) {
    cat("⚠️ Error in iteration", b, ":", conditionMessage(e), "\n")
    att_estimates[b] <- NA
  })
}

# Analyze Results
att_clean <- na.omit(att_estimates)
true_ATT <- -0.325  # update if needed
mean_ATT <- mean(att_clean)
bias <- mean_ATT - true_ATT
variance <- var(att_clean)
rmse <- sqrt(mean((att_clean - true_ATT)^2))

cat("\nMonte Carlo Results — Exogenous Adoption + Covariates:\n")
cat("Valid estimates: ", length(att_clean), "/", B, "\n")
cat("Mean ATT:        ", mean_ATT, "\n")
cat("Bias:            ", bias, "\n")
cat("Variance:        ", variance, "\n")
cat("RMSE:            ", rmse, "\n")

# Optional: Plot
hist(att_clean, breaks = 30, col = "skyblue",
     main = "Sampling Distribution (Exogenous Adoption)",
     xlab = "Estimated ATT")
abline(v = -0.585, col = "purple", lwd = 2, lty = 2)
abline(v = true_ATT, col = "red", lwd = 2, lty = 2)
abline(v = mean_ATT, col = "darkgreen", lwd = 2)
legend("topright", legend = c("First Estimate", "True ATT", "Mean Estimate"),
       col = c("purple", "red", "darkgreen"), lty = c(2, 2, 1), lwd = 2)

# === Plot trimmed (central 95%) distribution ===
lower <- quantile(att_clean, 0.025)
upper <- quantile(att_clean, 0.975)
att_trimmed <- att_clean[att_clean >= lower & att_clean <= upper]

hist(att_trimmed, breaks = 30, col = "skyblue",
     main = "Trimmed GSCM ATT Distribution (Central 95%)",
     xlab = "Estimated ATT (Volatility)")
abline(v = -0.585, col = "purple", lwd = 2, lty = 2)
abline(v = true_ATT, col = "red", lwd = 2, lty = 2)
abline(v = mean_ATT, col = "darkgreen", lwd = 2)
legend("topright", legend = c("First Estimate", "True ATT", "Mean Estimate"),
       col = c("purple", "red", "darkgreen"), lty = c(2, 2, 1), lwd = 2)


# This change is the DGP process make reduce bias. Mean ATT= -0.52. Much closer to true ATT=-0.32


# another option 
# ================================
# Monte Carlo Simulation with EXOGENOUS Adoption and stronger effect
# ================================

B <- 3  # Number of Monte Carlo simulations
att_estimates <- numeric(B)

for (b in 1:B) {
  cat("Simulation", b, "\n")
  
  # --- Step 1: Load base simulation script
  original_lines <- readLines("simulate_inflation_targeting_panel_LeeOcampos_export.R")
  
  # --- Step 2: Force exogenous adoption
  modified_lines <- gsub("use_lag_infl\\s*=\\s*TRUE", "use_lag_infl = FALSE", original_lines)
  
  # --- Step 3: Strengthen treatment effect (tau_it from -0.7 to -1.5)
  modified_lines <- gsub("VOL\\$tau_it.*?-0.7", "VOL$tau_it <- -1.5", modified_lines)
  
  # --- Step 4: Update random seed
  new_seed <- sample(1:99999999, 1)
  modified_lines <- gsub("SEED <- .*", paste0("SEED <- ", new_seed), modified_lines)
  
  # --- Step 5: Save modified script and run it
  writeLines(modified_lines, "temp_sim_strong_exog.R")
  source("temp_sim_strong_exog.R")  # creates student_out
  
  # --- Step 6: Compute volatility
  df <- as.data.frame(student_out) %>%
    group_by(country) %>%
    arrange(year) %>%
    mutate(infl_vol = rollapply(inflation, width = 5, FUN = sd, fill = NA, align = "right")) %>%
    ungroup() %>%
    filter(year >= 1984)
  
  # --- Step 7: Run GSCM
  tryCatch({
    res_gscm <- gsynth(
      Y = "infl_vol",
      D = "it.policy",
      X = c("cbi", "trade.open"),
      data = df,
      index = c("country", "year"),
      force = "unit",
      CV = TRUE,
      r = c(0, 5),
      se = FALSE,
      min.T0 = 7
    )
    
    att_estimates[b] <- mean(res_gscm$att.avg, na.rm = TRUE)
    
  }, error = function(e) {
    cat("Error in simulation", b, ":", conditionMessage(e), "\n")
    att_estimates[b] <- NA
  })
}

# === ANALYSIS ===
att_clean <- na.omit(att_estimates)

true_ATT <- -0.325
mean_ATT <- mean(att_clean)
bias <- mean_ATT - true_ATT
variance <- var(att_clean)
rmse <- sqrt(mean((att_clean - true_ATT)^2))

cat("\nMonte Carlo Results — Exogenous Adoption + Stronger Effect:\n")
cat("Valid estimates: ", length(att_clean), "/", B, "\n")
cat("Mean ATT:        ", mean_ATT, "\n")
cat("Bias:            ", bias, "\n")
cat("Variance:        ", variance, "\n")
cat("RMSE:            ", rmse, "\n")

# PLOT ===
hist(att_clean, breaks = 30, col = "skyblue",
     main = "Sampling Distribution (Strong Effect, Exogenous)",
     xlab = "Estimated ATT")
abline(v = true_ATT, col = "red", lwd = 2, lty = 2)
abline(v = mean_ATT, col = "darkgreen", lwd = 2)
legend("topright", legend = c("True ATT", "Mean Estimate"),
       col = c("red", "darkgreen"), lty = c(2, 1), lwd = 2)

# === Plot trimmed (central 95%) distribution ===
lower <- quantile(att_clean, 0.025)
upper <- quantile(att_clean, 0.975)
att_trimmed <- att_clean[att_clean >= lower & att_clean <= upper]

hist(att_trimmed, breaks = 30, col = "skyblue",
     main = "Trimmed GSCM ATT Distribution (Central 95%)",
     xlab = "Estimated ATT (Volatility)")
abline(v = -0.585, col = "purple", lwd = 2, lty = 2)
abline(v = true_ATT, col = "red", lwd = 2, lty = 2)
abline(v = mean_ATT, col = "darkgreen", lwd = 2)
legend("topright", legend = c("First Estimate", "True ATT", "Mean Estimate"),
       col = c("purple", "red", "darkgreen"), lty = c(2, 2, 1), lwd = 2)


# another option 
# ================================
# Monte Carlo Simulation with EXOGENOUS (prof suggestion)
# ================================

B <- 500  # Number of Monte Carlo simulations
att_estimates <- numeric(B)

for (b in 1:B) {
  cat("Simulation", b, "\n")
  
  # --- Step 1: Load base simulation script
  original_lines <- readLines("simulate_inflation_targeting_panel_LeeOcampos_export.R")
  
  # --- Step 2: Turn off lagged inflation in the hazard (make adoption more exogenous)
  modified_lines <- gsub("use_lag_infl\\s*=\\s*TRUE", "use_lag_infl = FALSE", original_lines)
  
  # --- Step 3: Remove institutional effect from hazard: adopt_inst = 0
  modified_lines <- gsub("adopt_inst\\s*=\\s*0\\.40","adopt_inst = 0.0", modified_lines)
  
  # --- Step 4: Update random seed
  new_seed <- sample(1:99999999, 1)
  modified_lines <- gsub("SEED <- .*", paste0("SEED <- ", new_seed), modified_lines)
  
  # --- Step 5: Save modified script and run it
  writeLines(modified_lines, "temp_sim_strong_exog.R")
  source("temp_sim_strong_exog.R")  # creates student_out
  
  # --- Step 6: Compute volatility
  df <- as.data.frame(student_out) %>%
    group_by(country) %>%
    arrange(year) %>%
    mutate(infl_vol = rollapply(inflation, width = 5, FUN = sd, fill = NA, align = "right")) %>%
    ungroup() %>%
    filter(year >= 1984)
  
  # --- Step 7: Run GSCM
  tryCatch({
    res_gscm <- gsynth(
      Y = "infl_vol",
      D = "it.policy",
      X = c("cbi", "trade.open"),
      data = df,
      index = c("country", "year"),
      force = "unit",
      CV = TRUE,
      r = c(0, 5),
      se = FALSE,
      min.T0 = 7
    )
    
    att_estimates[b] <- mean(res_gscm$att.avg, na.rm = TRUE)
    
  }, error = function(e) {
    cat("Error in simulation", b, ":", conditionMessage(e), "\n")
    att_estimates[b] <- NA
  })
}

# === ANALYSIS ===
att_clean <- na.omit(att_estimates)

true_ATT <- -0.325
mean_ATT <- mean(att_clean)
bias <- mean_ATT - true_ATT
variance <- var(att_clean)
rmse <- sqrt(mean((att_clean - true_ATT)^2))

cat("\nMonte Carlo Results — Exogenous Adoption + Stronger Effect:\n")
cat("Valid estimates: ", length(att_clean), "/", B, "\n")
cat("Mean ATT:        ", mean_ATT, "\n")
cat("Bias:            ", bias, "\n")
cat("Variance:        ", variance, "\n")
cat("RMSE:            ", rmse, "\n")

# PLOT ===
hist(att_clean, breaks = 30, col = "skyblue",
     main = "Sampling Distribution (Strong Effect, Exogenous)",
     xlab = "Estimated ATT")
abline(v = true_ATT, col = "red", lwd = 2, lty = 2)
abline(v = mean_ATT, col = "darkgreen", lwd = 2)
legend("topright", legend = c("True ATT", "Mean Estimate"),
       col = c("red", "darkgreen"), lty = c(2, 1), lwd = 2)

# === Plot trimmed (central 95%) distribution ===
lower <- quantile(att_clean, 0.025)
upper <- quantile(att_clean, 0.975)
att_trimmed <- att_clean[att_clean >= lower & att_clean <= upper]

hist(att_trimmed, breaks = 30, col = "skyblue",
     main = "Trimmed GSCM ATT Distribution (Central 95%)",
     xlab = "Estimated ATT (Volatility)")
abline(v = -0.833, col = "purple", lwd = 2, lty = 2)
abline(v = true_ATT, col = "red", lwd = 2, lty = 2)
abline(v = mean_ATT, col = "darkgreen", lwd = 2)
legend("topright", legend = c("First Estimate", "True ATT", "Mean Estimate"),
       col = c("purple", "red", "darkgreen"), lty = c(2, 2, 1), lwd = 2)





