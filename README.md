# Inflation Targeting and Inflation Volatility

Causal inference using Generalized Synthetic Control and Monte Carlo simulation.

---

##  Project Overview

This project studies whether **adopting an inflation-targeting (IT) regime reduces
inflation volatility**, or whether observed declines in volatility merely reflect
global disinflation trends since the 1980s.

Using a panel of **45 emerging economies (1980–2023)**, the analysis applies
modern causal inference methods designed for settings with
**heterogeneous units, staggered treatment adoption, and time-varying unobservables**.

- **Final paper / presentation (PDF):**  
  [Inflation-Targeting Regimes and Inflation Volatility](Causal_Inference_Project.pdf)


---

##  Research Question

- Does inflation targeting causally reduce inflation volatility?
- Or are post-adoption declines driven by global trends unrelated to IT?
- How do estimator bias and variance behave under different specifications?

---

##  Empirical Strategy

### Main Method
- **Generalized Synthetic Control Method (GSCM)**  
  (implemented via latent factor models)

Why GSCM:
- Parallel trends assumptions required by DiD are unlikely to hold
- Standard Synthetic Control is limited to single treated units
- GSCM accommodates:
  - Multiple treated units
  - Staggered adoption
  - Unobserved common shocks via latent factors

### Outcome
- Inflation volatility measured as a **5-year rolling standard deviation**
  of annual inflation.

### Covariates
- Central bank independence
- Trade openness
- Institutional quality
- Lagged inflation volatility
- GDP per capita
- Exchange rate regime

---

##  Monte Carlo Simulation

To evaluate estimator performance, the project conducts Monte Carlo simulations to study:

- Bias of GSCM under endogenous treatment adoption
- Variance and RMSE across specifications
- Coverage of bootstrap confidence intervals

Simulation results show that:
- Simple unit-FE GSCM overstates treatment effects
- Adding covariates reduces bias but does not eliminate it
- **Two-way fixed effects + covariates** substantially reduce bias and RMSE
- Making adoption exogenous in the DGP further improves estimator performance

---

