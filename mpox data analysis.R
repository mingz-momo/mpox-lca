## Data analysis

# ========================
# 0. Load required packages
# ========================
library(tidyverse)
library(stringr)
library(poLCA)
library(vcd)
library(FactoMineR)
library(factoextra)
library(missMDA)
library(scales)

# =========================
# 1. Data processing
# =========================

# Read control and case group data
control<-read.csv("data/control.csv", header = TRUE, sep = ',',
                  stringsAsFactors = FALSE, check.names = FALSE,
                  encoding = "UTF-8")
case<-read.csv("data/case.csv",header = TRUE, sep = ',',
                    stringsAsFactors = FALSE, check.names = FALSE,
                    encoding = "UTF-8")

# =========================
# 1.1. Data sorting
# =========================

# Subset variables of interest for control and case groups
control_select <- control[,c(1,9,12,13,14,19,20,21,22,24,25,26,27,28,29,30,31,
                           42,43,44,45,46,47,48,49,71,72,73,74,75,76,77,78,79,
                           80,81,82,83,84,85,87,88,89,90,91,98,105,112,119,
                           120,131)]
case_select <- case[,c(1,10,14,15,16,17,18,20,21,22,23,40,41,42,43, 44,45,46,
                       47,94,95,96,97,98,99,113,114,115,116,117,118,119,120,
                       121,122,123,124,125,126,127,129,130,131,132,133,143,153,
                       160,167,168,180)]

# Standardize variable names across datasets
colnames(control_select)[1:51]<-c("id","gender","mpox","HIV.status",
                                  "travel.3months","birth","marital.status",
                                  "education","occupation","income",
                                  "smallpox.vaccine","hypertension","diabetes",
                                  "hyperlipidemia","heart.disease",
                                  "liver.disease","cancer","syphilis",
                                  "chlamydia","gonorrhea","genital.herpes",
                                  "genital.mycoplasma","trichomonas.vaginalis",
                                  "genital.warts","lymphogranuloma.venereum",
                                  "orientation","anal.role","msm.freq",
                                  "msm.partners","msm.condom","msm.group",
                                  "Rush","poppers","capsule.zero","heroin",
                                  "ice","ecstasy.pills","K.powder","magou",
                                  "cannabis","hetero","hetero.freq",
                                  "hetero.condom","HIV.freq",
                                  "overseas.travel","domestic.travel",
                                  "msm.contact","hetero.contact",
                                  "overseas.party","domestic.party",
                                  "interview.date")
colnames(case_select)[1:51] <- c("id","mpox","gender","birth","marital.status",
                                 "education","occupation","travel.3months",
                                 "income","smallpox.vaccine","HIV.status",
                                 "syphilis","chlamydia","gonorrhea",
                                 "genital.herpes","genital.mycoplasma",
                                 "trichomonas.vaginalis","genital.warts",
                                 "lymphogranuloma.venereum","hypertension",
                                 "diabetes","hyperlipidemia","heart.disease",
                                 "liver.disease","cancer","orientation",
                                 "anal.role","msm.freq","msm.group",
                                 "msm.partners","msm.condom","Rush","poppers",
                                 "capsule.zero","heroin","ice","ecstasy.pills",
                                 "K.powder","magou","cannabis","hetero",
                                 "hetero.freq","hetero.partners",
                                 "hetero.condom","overseas.travel",
                                 "domestic.travel", "msm.contact",
                                 "hetero.contact","overseas.party",
                                 "domestic.party" ,"interview.date"
)

# Unifying the HIV codes of case and control groups
control_select$HIV.status[control_select$HIV.status == 3] <- 2 #In the initial data, the HIV-positive status was coded as 3 in the control group and 2 in the case group.

# Merge datasets into one analytical sample
combined_data <- dplyr::bind_rows(control_select, case_select)

# Derived variable construction
combined_data <- combined_data %>%
  mutate(
    
    # Merge highly related drugs
    rush.poppers = if_else(Rush == 1 | poppers == 1, 1, 0),
    
    # Other drug use indicator
    other.drugs = if_else(
      capsule.zero == 1 | heroin == 1 | ice == 1 |
        ecstasy.pills == 1 | K.powder == 1 |
        magou == 1 | cannabis == 1,
      1, 0
    )
  )

# STI composite variable (any STI infection)
combined_data$STI <- ifelse(
  combined_data$syphilis == 1 |
    combined_data$chlamydia == 1 |
    combined_data$gonorrhea == 1 |
    combined_data$genital.herpes == 1 |
    combined_data$genital.mycoplasma == 1 |
    combined_data$trichomonas.vaginalis == 1 |
    combined_data$genital.warts == 1 |
    combined_data$lymphogranuloma.venereum == 1,
  1, 0
)

# Age calculation and grouping
combined_data <- combined_data %>%
  mutate(
    # Standardize date format
    birth = gsub("/", "-", birth),
    interview.date = gsub("/", "-", interview.date),
    
    # Convert to Date type
    birth = as.Date(birth, format = "%Y-%m-%d"),
    interview.date = as.Date(interview.date, format = "%Y-%m-%d"),
    
    # Compute age
    age = as.numeric(difftime(interview.date, birth, units = "days")) / 365.25,
    
    # Categorize age
    age.group = case_when(
      age < 18 ~ "0",
      age >= 18 & age < 28 ~ "1",
      age >= 28 & age < 38 ~ "2",
      age >= 38 & age < 48 ~ "3",
      age >= 48 ~ "4",
      TRUE ~ NA_character_
    )
  )


# Remove redundant original variables after recoding
combined_data <- combined_data[
  , !(names(combined_data) %in% c("syphilis", "chlamydia", "gonorrhea",
                                  "genital.herpes", "genital.mycoplasma",
                                  "trichomonas.vaginalis", "genital.warts", 
                                  "lymphogranuloma.venereum", "Rush", "poppers",
                                  "capsule.zero", "heroin", "ice","ecstasy.pills",
                                  "K.powder", "magou", "cannabis",
                                  "id", "birth", "interview.date"))
]

# =========================
# 1.2. Quality control
# =========================

# Calculate missing rate per variable (Missing-data assessment)
missing_rate <- sapply(combined_data, function(x) mean(is.na(x)))

# Calculate dominant category proportion (Low-variance filtering)
max_prop <- sapply(combined_data, function(x) {
  if (is.numeric(x) && length(unique(x)) > 10) {
    return(NA)
  } else {
    tab <- table(as.factor(x))
    if (length(tab) == 0) return(NA)
    return(max(tab) / sum(tab))
  }
})

# Summarize QC metrics
flagged_vars <- data.frame(
  variable = names(combined_data),
  missing_rate = missing_rate,
  max_prop = max_prop,
  flag_missing = missing_rate > 0.2, # >20% missing
  flag_low_variance = max_prop > 0.90 # >90% single-level dominance
)

# Identify problematic variables
high_missing_vars <- flagged_vars %>%
  filter(flag_missing) %>%
  pull(variable)
low_variance_vars <- flagged_vars %>%
  filter(flag_low_variance) %>%
  pull(variable)

# Organize removal lists
var_lists <- list(
  high_missing_only = setdiff(high_missing_vars),
  low_variance_only = setdiff(low_variance_vars)
)

# Remove low-quality variables
combined_data <- combined_data %>%
  dplyr::select(-all_of(c(var_lists$high_missing_only,
                          var_lists$low_variance_only)))

# Multicollinearity check (Cramér's V)
data_factor_clean <- combined_data %>%
  dplyr::select(-mpox) %>%                   
  mutate(across(everything(), as.factor))

# Function for Cramér's V
cramers_v <- function(x, y) {
  tbl <- table(x, y, useNA = "no")    
  if (min(dim(tbl)) < 2) return(NA_real_)
  vcd::assocstats(tbl)$cramer
}

# Pairwise correlation matrix
vars <- names(data_factor_clean)
cramers_v_pairs <- expand_grid(var1 = vars, var2 = vars) %>%
  filter(var1 < var2) %>%
  mutate(
    cramer_v = map2_dbl(
      var1, var2,
      ~ cramers_v(data_factor_clean[[.x]], data_factor_clean[[.y]])
    )
  )

# Pairwise correlation matrix
high_corr_pairs <- cramers_v_pairs %>%
  filter(!is.na(cramer_v) & cramer_v > 0.7) %>% # Cramer's V >0.70
  arrange(desc(cramer_v))
print(high_corr_pairs)

# =========================
# 2. Latent Class Analysis (LCA)
# =========================

LCA_data_final<-combined_data

# Selected LCA indicators
lca_vars <- c("HIV.status","education","marital.status",
              "STI","orientation","anal.role",
              "msm.freq","msm.partners","msm.condom","msm.group",
              "domestic.travel",
              "rush.poppers","age.group")

# Build formula for poLCA
f <- as.formula(paste("cbind(", paste(lca_vars, collapse = ", "), ") ~ 1"))

# Convert to factors (required by poLCA)
LCA_data_final[lca_vars] <- lapply(LCA_data_final[lca_vars], function(x) {
  factor(x)   
})

# =========================
# 2.1. Fit multiple LCA models
# =========================

# Function to compute normalized entropy for LCA models
normalized_entropy <- function(posterior, eps = 1e-12) {
  posterior <- as.matrix(posterior)
  N <- nrow(posterior)
  K <- ncol(posterior)
  
  # If only one class, entropy is undefined
  if (K <= 1) return(NA_real_)
  
  # Avoid log(0) by adding a small constant
  p <- pmax(posterior, eps)
  
  # Raw entropy based on posterior class probabilities
  E_raw <- -sum(p * log(p))
  
  # Normalized entropy (scaled to 0-1)
  # 1 = perfect classification, 0 = complete uncertainty
  E_norm <- 1 - (E_raw / (N * log(K)))
  
  return(E_norm)
}

# Initialize an empty data frame to store model fit statistics
# k: number of latent classes
# AIC, BIC: information criteria for model comparison
# aBIC: sample-size adjusted Bayesian Information Criterion
# Entropy: classification quality measure
fit_stats <- data.frame(
  k = integer(),
  AIC = numeric(),
  BIC = numeric(),
  aBIC = numeric(),
  Entropy = numeric()
)

set.seed(123)  # Ensure reproducibility of model estimation
fit_list <- list()

# Fit latent class models with different numbers of classes (k = 1 to 6)
for (k in 1:6) {
  
  # Fit LCA model using poLCA
  # maxiter: maximum number of EM iterations
  # nrep: number of random starts to avoid local maxima
  # na.rm = FALSE ensures explicit handling of missing data
  fit <- poLCA(f, data = LCA_data_final, nclass = k,
               maxiter = 500, nrep = 10, na.rm = FALSE)
  
  # Store fitted model object for later inspection
  fit_list[[k]] <- fit
  
  # Compute sample-size adjusted BIC (aBIC)
  # Penalizes model complexity with an adjustment for sample size
  aBIC <- -2 * fit$llik + log((nrow(LCA_data_final) + 2) / 24) * fit$npar
  
  # Compute normalized entropy as a measure of classification quality
  ent01 <- normalized_entropy(fit$posterior)
  
  # Store model fit statistics for comparison across different k
  fit_stats <- rbind(fit_stats,
                     data.frame(k = k,
                                AIC = fit$aic,
                                BIC = fit$bic,
                                aBIC = aBIC,
                                Entropy = ent01))
  
  # Clear temporary model object to save memory
  rm(fit)
  gc()
}

# =========================
# 2.2. Class-specific proportion of mpox cases estimation
# =========================

# Extract 3-class LCA solution
fit3 <- fit_list[[3]] 

# Assign most likely latent class to each individual
LCA_data_final$class3 <- fit3$predclass

# Calculate proportion of mpox cases by latent class
LCA_data_final %>%
  group_by(class3) %>%
  summarise(
    n = n(),  # sample size per class
    mpox_cases = sum(mpox == 1, na.rm = TRUE),   # number of mpox cases
    mpox_rate = mpox_cases / n   # mpox cases proportion
  )

# =========================
# 2.3. Latent class profile extraction
# =========================

# Function to compute class-specific response distributions
# This function calculates the distribution (counts and proportions) of categorical responses within each latent class.
get_class_profile <- function(data, class_var, vars) {
  
  # Subset data to include only class membership variable and selected indicators
  data_subset <- data[, c(class_var, vars), drop = FALSE]
  
  # Remove observations with missing class assignment
  data_subset <- data_subset %>%
    dplyr::filter(!is.na(.data[[class_var]]))
  
  # Reshape data from wide to long format
  # Each row corresponds to one variable-level observation
  data_subset %>%
     tidyr::pivot_longer(
      cols = -1,  
      names_to = "variable",
      values_to = "level"
    ) %>%
    
    # Count frequency of each response category within each latent class
    dplyr::group_by(!!sym(class_var), variable, level) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
    
    # Compute within-class proportions for each category
    dplyr::mutate(prop = n / sum(n)) %>%
    
    # Remove grouping structure for downstream analysis
    dplyr::ungroup()
}

# Extract class profiles for 3-class model
profile3 <- get_class_profile(LCA_data_final, "class3", lca_vars)

# Save profile table for reporting
write.csv(profile3,"results/LCA-new.csv", row.names = FALSE)

# =========================
# 3. Vaccination prioritization analysis
# =========================

# This section evaluates different vaccine allocation strategies by estimating the proportion of mpox cases that would be covered under constrained vaccine supply scenarios.

# =========================
# 3.1 Risk scoring
# =========================

# Behavioral risk scoring

# Helper function: safely convert categorical variables to numeric values
to_num <- function(x) as.numeric(as.character(x))
df_behavior <- LCA_data_final %>%
  mutate(
    
    # Sexual behavior frequency (ordinal variable converted to numeric)
    # Missing values are treated as 0 (assumed no reported behavior)
    msm_freq_n = ifelse(is.na(msm.freq), 0, to_num(msm.freq)),
    
    # Number of sexual partners (treated as continuous/ordinal proxy of exposure intensity)
    # Missing values are assumed to represent no reported partners
    msm_partners_n = ifelse(is.na(msm.partners), 0, to_num(msm.partners)),
    
    # Group sex participation frequency (proxy for high-risk exposure contexts)
    # Missing values are imputed as 0 (no reported participation)
    msm_group_n = ifelse(is.na(msm.group), 0, to_num(msm.group)),
    
    # Condom use risk score (recoded into increasing risk levels)
    # Higher values indicate lower protection / higher behavioral risk
    condom_risk_n = case_when(
      msm.condom == "1" ~ 0,  # always condom use (lowest risk)
      msm.condom == "3" ~ 1,  # often use
      msm.condom == "2" ~ 2,  # occasionally use
      msm.condom == "4" ~ 3,  # never use (highest risk)
      is.na(msm.condom) ~ 0,  
      TRUE ~ 0
    ),
    
    # Drug use indicator
    # Binary coding: 1 = use, 0 = non-use or missing
    drug_use = case_when(
      rush.poppers == "1" ~ 1,
      rush.poppers == "0" ~ 0,
      is.na(rush.poppers) ~ 0,  
      TRUE ~ 0
    ),
    
    # Composite behavioral risk score
    # Equal-weight additive index as baseline specification
    score_behavior =
      1.0 * msm_freq_n +
      1.0 * msm_partners_n +
      1.0 * msm_group_n +
      1.0 * condom_risk_n +
      1.0 * drug_use
  )

# Clinical risk scoring

df_clinical <- df_behavior %>%
  mutate(
    
    # HIV status indicator
    HIV_pos = case_when(
      HIV.status == "1" ~ 1,        
      HIV.status == "2" ~ 0,  
      is.na(HIV.status) ~ 0,
      TRUE ~ 0
    ),
    
    # STI history indicator
    any_STI = case_when(
      STI == "1"  ~ 1,
      STI == "0" ~ 0,
      TRUE ~ 0   
    ),
    
    # Composite clinical risk score
    # Simple additive index assuming equal contribution of clinical risk factors
    score_clinical =
      1 * HIV_pos +
      1 * any_STI 
  )

# =========================
# 3.2 Case coverage evaluation
# =========================

# Sequence of vaccine supply scenarios
N_seq <- seq(0, 80, by = 5)

# Compute case coverage given a prioritization score
get_case_coverage <- function(df, score_var, N){
  
  # Number of people to vaccinate
  n_take <- min(N, nrow(df))
  
  # Select top-N highest risk individuals
  topN <- df %>%
    mutate(tiebreak = runif(nrow(df))) %>%                 
    arrange(desc(.data[[score_var]]), desc(tiebreak)) %>% 
    slice_head(n = n_take)
  
  # Proportion of total cases covered
  sum(topN$mpox == 1, na.rm = TRUE) / sum(df$mpox == 1, na.rm = TRUE)
}

# Bootstrap confidence intervals for coverage
bootstrap_coverage <- function(df, score_var, N, B = 1000){
  cov_boot <- replicate(B, {
    
    # resample individuals with replacement
    idx <- sample(seq_len(nrow(df)), replace = TRUE)
    df_b <- df[idx, ]
    
    # compute coverage
    get_case_coverage(df_b, score_var, N)
  })
  tibble(
    mean = mean(cov_boot, na.rm = TRUE),
    lwr  = quantile(cov_boot, 0.025, na.rm = TRUE),
    upr  = quantile(cov_boot, 0.975, na.rm = TRUE)
  )
}

# Generate full coverage curve across vaccine supply levels
coverage_curve <- function(df, score_var, strategy_name){
  purrr::map_dfr(N_seq, function(N){
    res <- bootstrap_coverage(df, score_var, N)
    tibble(strategy = strategy_name, N = N,
           coverage = res$mean, lwr = res$lwr, upr = res$upr)
  })
}

# =========================
# 3.3 Scenario-based strategies
# =========================

# Behavioral risk-based strategies

# Define behavioral components used in scoring
behavior_vars <- c("msm_freq_n", "msm_partners_n", "msm_group_n", "condom_risk_n", "drug_use")

# Pre-specified weighting scenarios
weight_scenarios <- list(
  equal = c(1, 1, 1, 1, 1),
  freq_focused = c(2, 1, 1, 1, 1),
  partners_focused = c(1, 2, 1, 1, 1),
  group_focused = c(1, 1, 2, 1, 1),
  condom_focused = c(1, 1, 1, 2, 1),
  drug_focused = c(1, 1, 1, 1, 2)
)

# Construct weighted behavior score
make_behavior_fixed <- function(df, weights){
  df %>%
    mutate(score_behavior_fixed =
             weights[1]*msm_freq_n +
             weights[2]*msm_partners_n +
             weights[3]*msm_group_n +
             weights[4]*condom_risk_n +
             weights[5]*drug_use)
}

# Evaluate all behavioral scenarios
curve_behavior_fixed_all <- purrr::map_dfr(names(weight_scenarios), function(scn){
  df_tmp <- make_behavior_fixed(df_behavior, weight_scenarios[[scn]])
  coverage_curve(df_tmp, "score_behavior_fixed", paste0("Behavioral risk-based (", scn, ")"))
})

# Clinical risk-based strategies

# Same framework applied to clinical risk indicators
clinical_vars <- c("HIV_pos", "any_STI")

# Grid search over HIV vs STI weighting
hiv_weights <- seq(0.3, 0.7, by = 0.1)

# Function to construct a weighted clinical risk score under a given HIV weighting scenario
make_clinical_fixed <- function(df, w_hiv){
  w_sti <- 1 - w_hiv
  df %>%
    mutate(score_clinical_fixed = w_hiv * HIV_pos + w_sti * any_STI)
}

# Evaluate clinical weighting scenarios
curve_clinical_fixed_all <- purrr::map_dfr(hiv_weights, function(w_hiv){
  df_tmp <- make_clinical_fixed(df_clinical, w_hiv)
  coverage_curve(df_tmp, "score_clinical_fixed",
                 paste0("Clinical risk-based (HIV=", w_hiv, ",STI=", round(1-w_hiv,1), ")"))
})

# =========================
# 3.4 LCA-informed prioritization strategy
# =========================

# Use posterior probability of high-risk class
vaccine_data <- LCA_data_final %>%
  mutate(score_lca = fit3$posterior[, 2])  # Class 3 has the highest proportion of mpox (in the dataset, Class 3 is assigned a value of 2)
curve_lca <- coverage_curve(vaccine_data, "score_lca", "LCA-informed")

# Combine all main strategies
curve_main_all <- dplyr::bind_rows(
  curve_behavior_fixed_all,
  curve_clinical_fixed_all,
  curve_lca
)

# Save the main results data
write.csv(curve_main_all, "results/main_results_coverage_curves.csv", row.names = FALSE)

# ============================================================
# 3.5 Spearman-derived weights
# ============================================================

# This function computes variable weights based on the absolute Spearman correlation between each risk factor and mpox outcome. 
# Weights are normalized to sum to a predefined total weight.

get_spearman_weights <- function(df, outcome, vars, total_weight){
  
  # Compute absolute Spearman correlation for each predictor
  w <- sapply(vars, function(v){
    
    # If variable not in dataset, assign zero weight
    if (!v %in% names(df)) return(0)
    
    # Compute correlation (robust to missing values)
    val <- suppressWarnings(abs(cor(df[[outcome]], df[[v]],
                                    use = "complete.obs",
                                    method = "spearman")))
    
    # Replace NA correlations with zero
    if (is.na(val)) val <- 0
    val
  })
  
  # If all correlations are zero, assign equal weights
  if (sum(w) == 0) w <- rep(1, length(vars))
  
  # Normalize weights to sum to total_weight
  w_norm <- w / sum(w) * total_weight
  names(w_norm) <- vars
  return(w_norm)
}

# Bootstrap evaluation of behavioral risk-based strategy
bootstrap_coverage_spearman_behavior <- function(df, N, B = 1000){
  cov_boot <- replicate(B, {
    
    # Bootstrap sampling
    idx <- sample(seq_len(nrow(df)), replace = TRUE)
    df_b <- df[idx, ]
    
    # Estimate Spearman-derived weights within bootstrap sample
    w <- get_spearman_weights(
      df = df_b,
      outcome = "mpox",
      vars = behavior_vars,
      total_weight = 5
    )
    
    # Construct weighted behavioral risk score
    df_b <- df_b %>%
      mutate(
        score_behavior_spear =
          w["msm_freq_n"]     * msm_freq_n +
          w["msm_partners_n"] * msm_partners_n +
          w["msm_group_n"]    * msm_group_n +
          w["condom_risk_n"]  * condom_risk_n +
          w["drug_use"]       * drug_use
      )
    
    # Compute case coverage for top-N prioritization
    get_case_coverage(df_b, "score_behavior_spear", N)
  })
  
  # Summarize bootstrap distribution
  tibble(
    mean = mean(cov_boot, na.rm = TRUE),
    lwr  = quantile(cov_boot, 0.025, na.rm = TRUE),
    upr  = quantile(cov_boot, 0.975, na.rm = TRUE)
  )
}

# Generate full coverage curve across vaccine supply levels (N)
curve_behavior_spear <- purrr::map_dfr(N_seq, function(N){
  res <- bootstrap_coverage_spearman_behavior(df_behavior, N, B = 1000)
  tibble(strategy="Behavioral risk-based (Spearman-derived weights)", N=N, coverage=res$mean, lwr=res$lwr, upr=res$upr)
})

# Clinical risk-based strategy (Spearman-derived weighted)

# Bootstrap evaluation of clinical risk-based strategy
bootstrap_coverage_spearman_clinical <- function(df, N, B = 1000){
  cov_boot <- replicate(B, {
    idx <- sample(seq_len(nrow(df)), replace = TRUE)
    df_b <- df[idx, ]
    
    # Estimate Spearman-derived weights within bootstrap sample
    w <- get_spearman_weights(df_b, "mpox", clinical_vars, total_weight = 2)
    
    # Construct weighted clinical risk score
    df_b <- df_b %>%
      mutate(
        score_clinical_spear = w["HIV_pos"] * HIV_pos + 
          w["any_STI"] * any_STI
      )
    
    # Compute case coverage for top-N prioritization
    get_case_coverage(df_b, "score_clinical_spear", N)
  })
  
  # Summarize bootstrap distribution
  tibble(
    mean = mean(cov_boot, na.rm = TRUE),
    lwr  = quantile(cov_boot, 0.025, na.rm = TRUE),
    upr  = quantile(cov_boot, 0.975, na.rm = TRUE)
  )
}

# Generate full coverage curve across vaccine supply levels (N)
curve_clinical_spear <- purrr::map_dfr(N_seq, function(N){
  res <- bootstrap_coverage_spearman_clinical(df_clinical, N, B = 1000)
  tibble(strategy="Clinical risk-based (Spearman-derived weights)", N=N, coverage=res$mean, lwr=res$lwr, upr=res$upr)
})


# Combine all strategies
curve_sa_all <- bind_rows(curve_behavior_spear, curve_clinical_spear, curve_lca)

# Save results for reproducibility
write.csv(curve_sa_all, "results/sensitivity_spearman_curves.csv", row.names = FALSE)

# =========================
# 4. Sensitivity Analysis
# =========================

# =========================
# 4.1 Latent Class Analysis (LCA)
# =========================

# Select LCA indicators (HIV excluded)
lca_vars_sens <- c("education","marital.status",
                   "STI","orientation","anal.role",
                   "msm.freq","msm.partners","msm.condom","msm.group",
                   "domestic.travel",
                   "rush.poppers","age.group")

# Build formula for poLCA
f_sens <- as.formula(paste("cbind(", paste(lca_vars_sens, collapse = ", "), ") ~ 1"))

# Initialize an empty data frame to store model fit statistics
fit_stats_sens <- data.frame(
  k = integer(),
  AIC = numeric(),
  BIC = numeric(),
  aBIC = numeric(),
  Entropy = numeric()
)
set.seed(123)  # Ensure reproducibility of model estimation
fit_list_sens <- list()

# Fit latent class models with different numbers of classes (k = 1 to 6)
for (k in 1:6) {
  
  # Fit LCA model using poLCA
  fit <- poLCA(f_sens, data = LCA_data_final, nclass = k,
               maxiter = 500, nrep = 10, na.rm = FALSE)
  
  # Store fitted model object for later inspection
  fit_list_sens[[k]] <- fit
  
  # Compute sample-size adjusted BIC (aBIC)
  aBIC <- -2 * fit$llik + log((nrow(LCA_data_final) + 2) / 24) * fit$npar
  
  # Compute normalized entropy as a measure of classification quality
  ent01 <- normalized_entropy(fit$posterior)
  
  # Store model fit statistics for comparison across different k
  fit_stats_sens <- rbind(fit_stats_sens,
                          data.frame(k = k,
                                     AIC = fit$aic,
                                     BIC = fit$bic,
                                     aBIC = aBIC,
                                     Entropy = ent01))
  rm(fit)  
  gc()  
}

# Extract 3-class LCA solution
fit_sens <- fit_list_sens[[3]]
LCA_data_sens<-LCA_data_final

# Assign most likely latent class to each individual
LCA_data_sens$class <- fit_sens$predclass

# Calculate proportion of mpox cases by latent class
LCA_data_sens %>%
  group_by(class) %>%
  summarise(
    n = n(),   # sample size per class
    mpox_cases = sum(mpox == 1, na.rm = TRUE),   # number of mpox cases
    mpox_rate = mpox_cases / n   # mpox cases proportion
  )

# Extract class profiles for 3-class model 
profile_sens <- get_class_profile(LCA_data_sens, "class", lca_vars)   # the reason for using lca_vars is to obtain a table that includes the HIV test results

# Save profile table for reporting
write.csv(profile_sens,"results/LCA-sensitivity.csv", row.names = FALSE)

# =========================
# 4.2 Vaccination prioritization analysis
# =========================

# This section evaluates different vaccine allocation strategies by estimating how many mpox cases would be covered under limited vaccine supply scenarios.

# Risk scoring

# Behavioral risk scoring

# Helper function: safely convert categorical variables to numeric values
to_num <- function(x) as.numeric(as.character(x))
df_behavior_sens <- LCA_data_sens %>%
  mutate(
    
    # Sexual behavior frequency (ordinal variable converted to numeric)
    # Missing values are treated as 0 (assumed no reported behavior)
    msm_freq_n     = ifelse(is.na(msm.freq), 0, to_num(msm.freq)),
    
    # Number of sexual partners (treated as continuous/ordinal proxy of exposure intensity)
    # Missing values are assumed to represent no reported partners
    msm_partners_n = ifelse(is.na(msm.partners), 0, to_num(msm.partners)),
    
    # Group sex participation frequency (proxy for high-risk exposure contexts)
    # Missing values are imputed as 0 (no reported participation)
    msm_group_n    = ifelse(is.na(msm.group), 0, to_num(msm.group)),
    
    # Condom use risk score (recoded into increasing risk levels)
    # Higher values indicate lower protection / higher behavioral risk
    condom_risk_n = case_when(
      msm.condom == "1" ~ 0,
      msm.condom == "3" ~ 1,
      msm.condom == "2" ~ 2,
      msm.condom == "4" ~ 3,
      is.na(msm.condom) ~ 0,  
      TRUE ~ 0
    ),
    
    # Drug use indicator
    # Binary coding: 1 = use, 0 = non-use or missing
    drug_use = case_when(
      rush.poppers == "1" ~ 1,
      rush.poppers == "0" ~ 0,
      is.na(rush.poppers) ~ 0,  
      TRUE ~ 0
    ),
    
    # Composite behavioral risk score
    # Equal-weight additive index as baseline specification
    score_behavior =
      1.0 * msm_freq_n +
      1.0 * msm_partners_n +
      1.0 * msm_group_n +
      1.0 * condom_risk_n +
      1.0 * drug_use
  )

# Clinical risk scoring

df_clinical_sens <- df_behavior_sens %>%
  mutate(
    
    # HIV status indicator
    HIV_pos = case_when(
      HIV.status == "1" ~ 1,        
      HIV.status == "2" ~ 0,  
      is.na(HIV.status) ~ 0,
      TRUE ~ 0
    ),
    
    # STI history indicator
    any_STI = case_when(
      STI == "1"  ~ 1,
      STI == "0" ~ 0,
      TRUE ~ 0   
    ),
    
    # Composite clinical risk score 
    # Simple additive index assuming equal contribution of clinical risk factors
    score_clinical =
      1 * HIV_pos +
      1  * any_STI 
  )

# Scenario-based strategies

# Evaluate all behavior scenarios
curve_behavior_fixed_sens <- purrr::map_dfr(names(weight_scenarios), function(scn){
  df_tmp <- make_behavior_fixed(df_behavior_sens, weight_scenarios[[scn]])
  coverage_curve(df_tmp, "score_behavior_fixed", paste0("Behavioral risk-based (", scn, ")"))
})

# Evaluate clinical weighting scenarios
curve_clinical_fixed_sens <- purrr::map_dfr(hiv_weights, function(w_hiv){
  df_tmp <- make_clinical_fixed(df_clinical_sens, w_hiv)
  coverage_curve(df_tmp, "score_clinical_fixed",
                 paste0("Clinical risk-based (HIV=", w_hiv, ",STI=", round(1-w_hiv,1), ")"))
})


# LCA-informed prioritization strategies

# Use posterior probability of high-risk class
vaccine_data_sens <- LCA_data_sens %>%
  mutate(score_lca = fit_sens$posterior[, 3])  # When the HIV status was removed from the LCA results, the mpox cases proportion of class 3 was the highest.
curve_lca_sens <- coverage_curve(vaccine_data_sens, "score_lca", "LCA-informed")

# Combine all strategies
curve_main_sens <- dplyr::bind_rows(
  curve_behavior_fixed_sens,
  curve_clinical_fixed_sens,
  curve_lca_sens
)

# Save the main results data
write.csv(curve_main_sens, "results/sensitivity_results_coverage_curves.csv", row.names = FALSE)

# Spearman-derived weights

# Behavioral risk-based strategy (Spearman-derived weights)

# Generate full coverage curve across vaccine supply levels (N)
curve_behavior_spear_sens <- purrr::map_dfr(N_seq, function(N){
  res <- bootstrap_coverage_spearman_behavior(df_behavior_sens, N, B = 1000)
  tibble(strategy="Behavioral risk-based (Spearman-derived weights)", N=N, coverage=res$mean, lwr=res$lwr, upr=res$upr)
})

# Clinical risk-based strategy (Spearman-derived weights)

# Generate full coverage curve across vaccine supply levels (N)
curve_clinical_spear_sens <- purrr::map_dfr(N_seq, function(N){
  res <- bootstrap_coverage_spearman_clinical(df_clinical_sens, N, B = 1000)
  tibble(strategy="Clinical risk-based (Spearman-derived weights)", N=N, coverage=res$mean, lwr=res$lwr, upr=res$upr)
})

# Combine all strategies
curve_sa_sens <- bind_rows(curve_behavior_spear_sens, curve_clinical_spear_sens, curve_lca_sens)

# Save results for reproducibility
write.csv(curve_sa_sens, "results/sensitivity_spearman_curves_results.csv", row.names = FALSE)