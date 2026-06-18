# Cost-Effectiveness of Lenacapavir vs TDF/FTC for HIV PrEP in MSM in Spain

A deterministic, compartmental HIV transmission model used to evaluate the
cost-effectiveness of subcutaneous lenacapavir (LEN) versus oral
tenofovir/emtricitabine (TDF/FTC) for pre-exposure prophylaxis (PrEP) among
men who have sex with men (MSM) in Spain.

## Model overview

The model tracks a population of ~893,000 MSM stratified by:

- **Age** (1–85 years)
- **Sexual risk group** (very-high, high, low, very-low)
- **Disease state** (susceptible → 7 chronic HIV stages → late-stage → AIDS)
- **PrEP status** (off / on)
- **ART status** (untreated / on treatment)

The state space is a five-dimensional array `C[t, a, s, r, p]` with weekly
time steps over a 20-year horizon. PrEP is available to the very-high and
high risk groups only; lenacapavir and TDF/FTC differ in their residual
infection multiplier (`RiskRed = 1 − efficacy`) and weekly PrEP cost.

Proportionate mixing across risk groups drives the force of infection. The
model is calibrated to 2019–2023 Spanish epidemiological data using Latin
Hypercube Sampling and local optimisation.

## Repository structure

### Core simulation

| File | Function | Purpose |
|------|----------|---------|
| `HIVSIM02exp05.R` | `HIVSIM()` | Transmission model used during **calibration** (simplified parameter interface) |
| `HIVSIM05.R` | `HIVSIMPREP()` | Transmission model used for **CEA** (full parameter passing, PrEP compartments) |
| `HIVSIM03.R`–`HIVSIM06.R` | various | Earlier/alternative model versions (not used in the main analysis) |
| `SettingArray.R` | — | Initialises compartment arrays by risk group |

Both `HIVSIM()` and `HIVSIMPREP()` implement the same epidemiological model.
`HIVSIMPREP()` accepts all parameters as explicit arguments (safe for
parallel cluster execution), while `HIVSIM()` reads some from the global
environment.

### Parameter sampling

| File | Purpose |
|------|---------|
| `invbeta.R` | Draw from Beta(mean, sd) — used in PSA |
| `invbeta_u.R` | Beta via inverse CDF of a uniform — used in LHS calibration |
| `invgamma.R` / `invgamma_u.R` | Gamma distribution (random / inverse-CDF) |
| `invlognorm.R` / `invlognorm_u.R` | Log-normal distribution |
| `rdirichlet_u.R` | Dirichlet via inverse CDF |
| `draw_params_from_u.R` | Master function: maps a vector of 28 uniform(0,1) values to all model parameters (used by LHS calibration) |

### Economic evaluation

| File | Purpose |
|------|---------|
| `QalyCalc.R` | Discounted QALY accumulation (3% annual rate, weekly steps) |
| `CostCalc.R` | Discounted cost accumulation |
| `QALYCost.R` | Aggregates simulation output → total QALYs + costs (no-PrEP scenario) |
| `QALYCostP.R` | Same, with PrEP compartment costs included |
| `DALYCost.R` | DALY-based alternative (not used in main analysis) |

`QALYCostP.R` reads `PrEPcost`, `Treatcost`, `AIDScost`, and QALY weights
(`WAIDS`, `WHIV`, `WM350`, `WTREAT`, `WSusceptible`) from the global
environment.

### Data processing utilities

| File | Purpose |
|------|---------|
| `make_slices_HIVSIM.R` | Index slices for extracting compartments from the output array |
| `rs_at.R` | Row-sum utility for compartment extraction by age |
| `sim_weekly_to_yearly.R` | Aggregates weekly simulation output to yearly summaries |

### Input data

| File | Contents |
|------|----------|
| `MortPobGen.csv` | Age-specific mortality rates (general population, Spain) |
| `Costes PREP.xlsx` | PrEP drug cost data |
| `Best_theta_Param` | Calibrated parameter vector (28 parameters, R binary) |

### Analysis scripts (R Markdown)

Run these in the order listed below.

#### 1. Calibration

| File | Purpose |
|------|---------|
| `PrEP_CEA_Calibration_LHS.Rmd` | Latin Hypercube Sampling calibration (5,000 runs). Produces `ResultParam` and initial parameter rankings. |
| `PrEP_CEA_Calibration_LHS_LocOpt.Rmd` | Local optimisation starting from the best LHS candidates. Produces `Best_theta_Param`. |

**Depends on:** `HIVSIM02exp05.R`, `HIVSIMCAL_LCA_LHS.R`,
`draw_params_from_u.R`, all `inv*_u.R` functions, `MortPobGen.csv`.

#### 2. Base-case cost-effectiveness analysis

| File | Purpose |
|------|---------|
| `PrEP_CEA_LHS.Rmd` | Runs the 20-year simulation for TDF/FTC and LEN at the calibrated point estimates. Computes base-case ICERs, incremental QALYs, and incremental costs. Produces summary tables and the CE plane figure. |

**Depends on:** `HIVSIM05.R` (for `HIVSIMPREP()`), `QALYCostP.R`,
`QalyCalc.R`, `CostCalc.R`, `Best_theta_Param`, `MortPobGen.csv`.

#### 3. Sensitivity analyses

| File | Purpose |
|------|---------|
| `PrEP_CEA_OWSA_LCA_LHS.Rmd` | One-way sensitivity analysis (±20% on ~29 parameters). Produces tornado plot data. |
| `PrEP_CEA_PSA_paired.Rmd` | **Paired** probabilistic sensitivity analysis. Both arms share common parameter draws per iteration; only `RiskRed` and `PrEPcost` differ. Produces CE plane, CEAC, and saved result matrices (`MAT2_LHS_R`, `MAT2LEN_R`). |

> **Note:** `PrEP_CEA_PSA_LHS.Rmd` is an earlier, **unpaired** version of the
> PSA where TDF and LEN arms drew parameters independently. This produced
> spurious variance in the incremental results. Use `PrEP_CEA_PSA_paired.Rmd`
> instead.

**Depends on:** same as step 2, plus `invbeta.R`, `invgamma.R`,
`DirichletReg`, `parallel`.

#### 4. Figures

| File | Purpose |
|------|---------|
| `PSA figures.Rmd` | Loads saved PSA matrices and produces publication-quality CE plane, CEAC, and combined TIFF figures. |

**Depends on:** `MAT2_LHS_R`, `MAT2LEN_R` (produced by step 3).

#### 5. Reporting

| File | Purpose |
|------|---------|
| `CHEERS_Report.Rmd` | Full manuscript-style report following the CHEERS 2022 checklist. |

## How to run

### Prerequisites

```r
install.packages(c(
  "ggplot2", "dplyr", "tidyr", "parallel",
  "DirichletReg", "BCEA", "lhs", "cowplot",
  "scales", "flextable", "officer"
))
```

R ≥ 4.1 is required (the code uses the base pipe `|>`).

### Execution order

1. **Open `PrEP-LEN.Rproj`** in RStudio.

2. **Calibration** (only needed once, or when changing model structure):
   - Knit `PrEP_CEA_Calibration_LHS.Rmd`
   - Knit `PrEP_CEA_Calibration_LHS_LocOpt.Rmd`
   - This produces `Best_theta_Param` (~15 min on 20 cores).

3. **Base-case CEA**:
   - Knit `PrEP_CEA_LHS.Rmd`

4. **Sensitivity analyses**:
   - Knit `PrEP_CEA_OWSA_LCA_LHS.Rmd` (OWSA; slow — runs ~60 simulations)
   - Knit `PrEP_CEA_PSA_paired.Rmd` (PSA; adjust `repeats` for desired
     sample size; 5 repeats × n_cores by default)

5. **Figures**:
   - Knit `PSA figures.Rmd`

### Parallelism

The calibration and PSA scripts use `parallel::makeCluster()` with
`detectCores()`. On a machine with *n* cores and `repeats = k`, the PSA
produces *n × k* iterations. For a full analysis, set `repeats` to at least
50 (giving ≥1,000 iterations on a 20-core machine).

## Key parameters

| Parameter | TDF/FTC | Lenacapavir | Source |
|-----------|---------|-------------|--------|
| Efficacy (point estimate) | 86% | 96% | Grant et al. (iPrEx); Bekker et al. (PURPOSE 1) |
| Residual risk (`RiskRed`) | 0.14 | 0.04 | 1 − efficacy |
| PSA SE for efficacy | 0.0973 | 0.0434 | Derived from trial 90%/95% CIs |
| Weekly PrEP cost | €26.94 | €898.82 | — |
| Weekly PrEP uptake (`W`) | 0.00192 | 0.00192 | SiPrEP programme data |
| Weekly PrEP discontinuation (`Woff`) | 0.0038 | 0.0038 | SiPrEP (~15%/year) |
| Discount rate | 3%/year | 3%/year | Spanish HTA guidelines |
| Time horizon | 20 years | 20 years | — |

## Data flow diagram

```
MortPobGen.csv ──┐
                 ├──► Calibration (LHS + local opt.) ──► Best_theta_Param
invbeta_u.R etc.─┘                                              │
                                                                 ▼
                                                     ┌──────────────────────┐
                                                     │  HIVSIMPREP()        │
                                                     │  (HIVSIM05.R)        │
                                                     │                      │
                                                     │  Shared parameters   │
                                                     │  drawn once per PSA  │
                                                     │  iteration           │
                                                     └────┬────────────┬────┘
                                                          │            │
                                                   RiskRed_TDF   RiskRed_LEN
                                                   PrEPcost=27   PrEPcost=899
                                                          │            │
                                                          ▼            ▼
                                                     QALYCostP()  QALYCostP()
                                                          │            │
                                                          └─────┬──────┘
                                                                │
                                                         ΔQALYs, ΔCosts
                                                                │
                                                       ┌────────┴────────┐
                                                       │                 │
                                                    CE plane           CEAC
                                                    (scatter)      (NMB-based)
```


