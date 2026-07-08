## Visualization

source("mpox data analysis.R")

# =========================
# 1. Latent Class Analysis (LCA)
# =========================

# Model fit visualization (Figure S1)
fit_stats_unique <- fit_stats %>% distinct(k, .keep_all = TRUE)
pdf("results/LCA_fit_metrics.pdf", width = 10, height = 6)

# Colors for plotting
col_aic <- "#8D2F25"
col_bic <- "#3E608D"
col_abic <- "#C47F1D"   
col_entropy <- "#72b063"
par(mar = c(5, 5, 4, 5))

# Plot AIC, BIC, aBIC
plot(fit_stats_unique$k, fit_stats_unique$AIC, type = "b",
     pch = 16, col = col_aic, lty = 1, lwd = 2,
     xlab = "Number of Classes (k)",
     ylab = "AIC / BIC / aBIC",
     ylim = range(c(fit_stats_unique$AIC, fit_stats_unique$BIC, fit_stats_unique$aBIC), na.rm = TRUE),
     xaxt = "n", main = "LCA Model Fit Indices by Number of Classes")
axis(1, at = fit_stats_unique$k, labels = fit_stats_unique$k)
lines(fit_stats_unique$k, fit_stats_unique$BIC, type = "b",
      pch = 17, col = col_bic, lty = 1, lwd = 2)
lines(fit_stats_unique$k, fit_stats_unique$aBIC, type = "b",
      pch = 18, col = col_abic, lty = 1, lwd = 2)

# Grid lines
abline(h = seq(floor(min(fit_stats_unique$AIC, na.rm = TRUE)/100)*100,
               ceiling(max(fit_stats_unique$BIC, na.rm = TRUE)/100)*100, by = 50),
       col = "lightgray", lty = "dotted")
abline(v = fit_stats_unique$k, col = "lightgray", lty = "dotted")

# Entropy (secondary axis)
ent_ok <- !is.na(fit_stats_unique$Entropy)
par(new = TRUE)
plot(fit_stats_unique$k[ent_ok], fit_stats_unique$Entropy[ent_ok], type = "b",
     pch = 15, col = col_entropy, lty = 2, lwd = 2,
     axes = FALSE, xlab = "", ylab = "",
     xlim = range(fit_stats_unique$k),
     ylim = range(fit_stats_unique$Entropy, na.rm = TRUE))
axis(side = 4)
mtext("Entropy", side = 4, line = 3)
pch_aic <- 16
pch_bic <- 17
pch_abic <- 18
pch_entropy <- 15

# Legend
legend("topright",
       legend = c("AIC", "BIC", "aBIC", "Entropy"),
       col = c(col_aic, col_bic, col_abic, col_entropy),
       pch = c(pch_aic, pch_bic, pch_abic, pch_entropy),
       lty = c(1, 1, 1, 2),
       lwd = 2,
       bty = "o",
       bg = "white",
       cex = 0.9)
dev.off()

# =========================
# 2. Mpox radar plot
# =========================

# Define numeric coding scheme for radar plot normalization
# (maps categorical responses into standardized 0-1 scale)
radar_codes <- list(
  
  # clinical characteristics
  "HIV.status" = c("1" = 1, "2" = 0),                 
  "STI" = c("0" = 0, "1" = 1),    
  
  # demographic characteristic
  "age.group" = c("1" = 0, "2" = 1/3, "3" = 2/3, "4" = 1),
  "education" = c("1" = 0, "2" = 1/3, "3" = 2/3, "4" = 1),
  "marital.status" = c("1" = 0, "3" = 1/3, "2" = 2/3, "4" = 1),
  "orientation" = c("1" = 0, "2" = 1, "4" = 0.5),
  
  # behavioral characteristic
  "anal.role" = c("1" = 0, "2" = 0.5, "3" = 1),
  "msm.condom" = c("1" = 0, "3" = 1/3, "2" = 2/3, "4" = 1),
  "msm.freq" = c("1" = 0, "2" = 0.5, "3" = 1),
  "msm.partners" = c("1" = 0, "2" = 0.5, "3" = 1),
  "msm.group" = c("1" = 0, "2" = 1/3, "3" = 2/3, "4" = 1),
  "rush.poppers" = c("0" = 0, "1" = 1),
  "domestic.travel" = c("2" = 0, "1" = 1)         
)

# Define plotting order of variables in radar chart
var_radar_order <- c(
  "HIV.status", "STI",
  "age.group", "education", "marital.status", "orientation", 
  "anal.role", "msm.condom", "msm.freq", "msm.partners", "msm.group", "rush.poppers",
  "domestic.travel"
)

# Convert LCA output into long format (class × variable × category)
lca_radar_data <- map_dfr(lca_vars, function(v) {
  LCA_data_final %>%
    filter(!is.na(class3)) %>%
    mutate(
      class3 = as.factor(class3),
      level = as.character(.data[[v]])
    ) %>%
    count(class3, level, name = "n") %>%
    
    # Compute within-class proportions
    group_by(class3) %>%
    mutate(
      prop = n / sum(n),
      variable = v
    ) %>%
    ungroup() %>%
    dplyr::select(class3, variable, level, n, prop)
})

# Re-label latent classes for interpretability
lca_radar_data <- lca_radar_data %>%
  mutate(class3 = as.character(class3)) %>%
  mutate(class3 = case_when(
    
    # Reorder class labels for consistency in visualization
    class3 == "2" ~ "3",
    class3 == "3" ~ "2",
    TRUE ~ class3
  )) %>%
  mutate(class3 = factor(class3))

# Compute weighted mean scores for radar plot
lca_radar_main <- lca_radar_data %>%
  mutate(
    class3 = as.numeric(as.character(class3)),
    class3 = paste0("Class", class3),
    variable = as.character(variable),
    level = as.character(level)
  ) %>%
  
  # Map categorical responses to numeric radar scores
  group_by(class3, variable) %>%
  mutate(
    code = radar_codes[[variable[1]]][level]
  ) %>%
  
  # Compute weighted mean score per variable per class
  summarise(
    weighted_mean = sum(as.numeric(code) * prop, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  
  # Convert to wide format (class × variable matrix)
  pivot_wider(names_from = variable, values_from = weighted_mean) %>%
  as.data.frame()

# Set class labels as row names
rownames(lca_radar_main) <- lca_radar_main$class3
lca_radar_main$class3 <- NULL

# Reorder variables for radar display
lca_radar_main <- lca_radar_main[, var_radar_order, drop = FALSE]

# Add min/max reference for radar scaling
radar_df <- rbind(
  max = rep(1, ncol(lca_radar_main)),
  min = rep(0, ncol(lca_radar_main)),
  lca_radar_main
)

# Mpox radar visualization (Figure 2A)
pdf("results/Radar_main.pdf", width = 10, height = 10)

# Convert to long format for ggplot2
plot_radar_all <- lca_radar_main %>%
  rownames_to_column("class3") %>%
  pivot_longer(
    cols = -class3,
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    variable = factor(variable, levels = var_radar_order),
    id = as.numeric(variable)
  )

# Compute mean profile across classes
mean_points_df <- plot_radar_all %>%
  group_by(variable, id) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

# Prepare axis labels for polar coordinates
label_df <- plot_radar_all %>%
  distinct(variable, id) %>%
  arrange(id) %>%
  mutate(
    angle = 90 - 360 * (id - 0.5) / n(),
    hjust = ifelse(angle < -90, 1, 0),
    angle = ifelse(angle < -90, angle + 180, angle)
  )

# Reorder values within variable (for consistent stacking)
plot_radar_all <- plot_radar_all %>%
  group_by(variable) %>%
  arrange(desc(value), .by_group = TRUE) %>%
  ungroup()

# Background grid lines
grid_df <- data.frame(y = c(0.25, 0.5, 0.75, 1.0))

# Spokes (radial axes)
spoke_df <- data.frame(id = seq_along(var_radar_order))

# Generate radar plot
ggplot(plot_radar_all, aes(x = factor(id), y = value, fill = class3)) +
  
  # Background radial grid
  geom_hline(
    data = grid_df,
    aes(yintercept = y),
    color = "grey80",
    linewidth = 0.4
  ) +
  
  # Class-specific bars
  geom_col(
    position = "identity",
    alpha = 1,
    width = 0.9,
    color = "white"
  ) +
  
  # Radial spokes
  geom_segment(
    data = spoke_df,
    aes(x = factor(id), xend = factor(id), y = 0, yend = 1.02),
    color = "grey35",
    linetype = "dashed",
    linewidth = 0.35,
    inherit.aes = FALSE
  ) +
  
  # Mean profile points
  geom_point(
    data = mean_points_df,
    aes(x = factor(id), y = mean_value),
    inherit.aes = FALSE,
    shape = 21,
    size = 1.5,
    fill = "black",
    color = "black"
  )+
  
  # Convert to polar coordinate system
  coord_polar(start = 0) +
  ylim(-0.35, 1.15) +
  
  # Variable labels
  geom_text(
    data = label_df,
    aes(x = factor(id), y = 1.08, label = variable, angle = angle, hjust = hjust),
    inherit.aes = FALSE,
    size = 3.6
  ) +
  
  # Manual class colors
  scale_fill_manual(values = c(
    "Class3" = "#F2724D",
    "Class2" = "#4683B4",
    "Class1" = "#67B3AD"
  )) +
  
  # Clean theme for publication-quality figure
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  )
plot_radar_all
dev.off()

# =========================
# 3. MCA analysis (Multiple Correspondence Analysis)
# =========================

# Prepare dataset for MCA
mca_df <- LCA_data_final

# Select variables used in MCA projection
mca_df <- mca_df %>%
  dplyr::select(
    HIV.status, STI,
    age.group, education, marital.status, orientation,
    anal.role, msm.condom, msm.freq, msm.partners, msm.group, rush.poppers,domestic.travel,
    class3, mpox
  )

# Reorder class labels for consistency across analysis
mca_df <- mca_df %>%
  mutate(class3 = as.character(class3)) %>%
  mutate(class3 = case_when(
    class3 == "2" ~ "3",   # class 2 was the group with the highest mpox cases proportion in original data
    class3 == "3" ~ "2",
    TRUE ~ class3
  )) %>%
  mutate(class3 = factor(class3))

# Convert all variables to categorical factors (required for MCA)
mca_df <- mca_df %>%
  mutate(across(everything(), as.factor))

# Relabel categorical variables for interpretability
mca_df <- mca_df %>%
  mutate(
    HIV.status = factor(HIV.status,
                        levels = c("1", "2"),
                        labels = c("HIV+", "HIV-")),
    
    STI = factor(STI,
                 levels = c("0", "1"),
                 labels = c("STI-", "STI+")),
    
    age.group = factor(age.group,
                       levels = c("1", "2", "3", "4"),
                       labels = c("18-28", "28-38", "38-48", "48-58")),
    
    education = factor(education,
                       levels = c("1", "2", "3", "4"),
                       labels = c("Junior_or_below",
                                  "Highschool_or_vocational",
                                  "Bachelor_or_associate",
                                  "Graduate_plus")),
    
    marital.status = factor(marital.status,
                            levels = c("1", "2", "3", "4"),
                            labels = c("Unmarried",
                                       "Married",
                                       "Cohabiting",
                                       "Divorced_or_widowed")),
    
    orientation = factor(orientation,
                         levels = c("1", "2", "4"),
                         labels = c("Gay", "Bisexual", "Not_sure")),
    
    domestic.travel = factor(domestic.travel,
                                levels = c("1", "2"),
                                labels = c("Local_travel", "Non_local_travel")),
    
    anal.role = factor(anal.role,
                       levels = c("1", "2", "3"),
                       labels = c("Insertive", "Versatile", "Receptive")),
    
    msm.condom = factor(msm.condom,
                        levels = c("1", "2", "3", "4"),
                        labels = c("Always", "Occasionally", "Often", "Never")),
    
    msm.freq = factor(msm.freq,
                      levels = c("1", "2", "3"),
                      labels = c("<=1_per_month", "2_4_per_month", ">=5_per_month")),
    
    msm.partners = factor(msm.partners,
                          levels = c("1", "2", "3"),
                          labels = c("<=1_partner", "2_3_partners", ">=4_partners")),
    
    msm.group = factor(msm.group,
                       levels = c("1", "2", "3", "4"),
                       labels = c("0_groupsex", "1_2_groupsex", "3_4_groupsex", ">=5_groupsex")),
    
    rush.poppers = factor(rush.poppers,
                          levels = c("0", "1"),
                          labels = c("No_poppers", "Yes_poppers")),
    
    class3 = factor(class3,
                    levels = c("1", "2", "3"),
                    labels = c("Class 1", "Class 2", "Class 3")),
    
    mpox = factor(mpox,
                  levels = c("1", "2"),
                  labels = c("Case", "Control"))
  )

# Estimate number of dimensions for MCA imputation
nb <- estim_ncpMCA(mca_df, ncp.max = 5)

# Impute missing categorical values using MCA
mca_imp <- imputeMCA(mca_df, ncp = nb$ncp)

# Perform MCA on completed dataset
res_mca <- MCA(mca_imp$completeObs,
               quali.sup = c(ncol(mca_df)-1, ncol(mca_df)),
               graph = FALSE)

# Eigenvalues (explained variance per dimension)
eig <- get_eigenvalue(res_mca)
print(eig)
eig[1:5, ]

# Variable contributions to dimensions
var_contrib <- res_mca$var$contrib

# Extract individual coordinates from MCA space
ind_coord <- as.data.frame(res_mca$ind$coord) %>%
  rownames_to_column("id")

# Combine MCA coordinates with class labels
plot_mca_df <- bind_cols(ind_coord, mca_df %>% dplyr::select(class3, mpox))

# MCA visualization (Figure 2B)
pdf("results/LCA-MCA.pdf", width = 10, height = 6)
p_class <- ggplot(plot_mca_df, aes(x = `Dim 1`, y = `Dim 2`, color = class3)) +
  
  # Reference axes
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  
  # Individual points in MCA space
  geom_point(size = 2, alpha = 0.75) +
  
  # Confidence ellipse for each latent class
  stat_ellipse(aes(fill = class3),
               geom = "polygon",
               alpha = 0.15,
               color = NA) +
  
  # Axis labels with explained variance
  labs(
    title = "Projection of LCA classes onto MCA risk space",
    x = paste0("Dim 1 (", round(eig[1, 2], 1), "%)"),
    y = paste0("Dim 2 (", round(eig[2, 2], 1), "%)"),
    color = "LCA class",
    fill = "LCA class"
  ) +
  
  # Class color scheme
  scale_color_manual(values = c(
    "Class 3" = "#F2724D",
    "Class 2" = "#4683B4",
    "Class 1" = "#67B3AD"
  )) +
  scale_fill_manual(values = c(
    "Class 3" = "#F2724D",
    "Class 2" = "#4683B4",
    "Class 1" = "#67B3AD"
  )) +
  
  # # Clean theme for publication-quality figure
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )
p_class
dev.off()

# =========================
# 4. Mpox case proportion by latent class
# =========================

plot_proportion <- LCA_data_final

# Reorder class labels consistently
plot_proportion <- plot_proportion %>%
  mutate(class3 = as.character(class3)) %>%
  mutate(class3 = case_when(
    class3 == "2" ~ "3",   # class 2 was the group with the highest mpox cases proportion in original data
    class3 == "3" ~ "2",
    TRUE ~ class3
  )) %>%
  mutate(class3 = factor(class3))

# Compute class-specific proportion of mpox cases
plot_proportion <- plot_proportion %>%
  mutate(
    class3 = factor(class3,
                    levels = c(1, 2, 3),
                    labels = c("Class 1", "Class 2", "Class 3"))
  ) %>%
  group_by(class3) %>%
  summarise(
    n = n(),   # total per class
    cases = sum(mpox == 1, na.rm = TRUE),   #number of cases
    rate = cases / n   # proportion of cases
  ) %>%
  ungroup()

# Visualization of mpox cases proportion (Figure 2B)
pdf("results/LCA-proportion.pdf", width = 15, height = 3)

p_proportion_bar <- ggplot(plot_proportion,
                           aes(x = class3, y = rate, fill = class3)) +
  
  # Bar plot of proportion
  geom_col(width = 0.65, color = "black") +
  
  # Add percentage labels
  geom_text(aes(label = percent(rate, accuracy = 0.1)),
            hjust = -0.2,
            size = 4) +
  
  # Format y-axis as percentage
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),          
    breaks = seq(0, 1, by = 0.2),  
    expand = c(0, 0)
  ) +
  
  # Class colors
  scale_fill_manual(values = c(
    "Class 3" = "#F2724D",
    "Class 2" = "#4683B4",
    "Class 1" = "#67B3AD"
  )) +
  
  # Horizontal bar layout
  coord_flip() +
  
  # Axis labels
  labs(
    x = "Latent class",
    y = "Mpox prevalence"
  ) +
  
  # Clean theme for publication-quality figure
  theme_classic() +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 11),
    legend.position = "none"
  )
p_proportion_bar
dev.off()

# =========================
# 5. Vaccination prioritization analysis
# =========================

# Scenario-based strategies

# Recode strategies into broader methodological groups for comparison
curve2 <- curve_main_all %>%
  mutate(group = case_when(
    str_detect(strategy, "^Behavioral risk-based") ~ "Behavioral risk-based",
    str_detect(strategy, "^Clinical risk-based")  ~ "Clinical risk-based",
    str_detect(strategy, "^LCA")            ~ "LCA-informed",
    TRUE ~ "Other"
  ))

# Compute uncertainty bands per strategy group
band <- curve2 %>%
  filter(group %in% c("Behavioral risk-based", "Clinical risk-based", "LCA-informed")) %>%
  group_by(group, N) %>%
  summarise(
    lwr_band = min(lwr, na.rm = TRUE),
    upr_band = max(upr, na.rm = TRUE),
    .groups = "drop"
  )

# Define fill colors for uncertainty ribbons by strategy group
group_fill <- c(
  "Behavioral risk-based" = "#deebf7",  
  "Clinical risk-based"  = "#e5f5e0",  
  "LCA-informed"   = "#fee0d2"   
)

# Behavioral risk-based strategies (gradient blues)
behavior_cols <- c("#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#084594")

# Clinical risk-based strategies (gradient greens)
clinical_cols <- c("#c7e9c0", "#a1d99b", "#74c476", "#41ab5d", "#238b45")

# LCA-informed strategy (single representative red)
lca_cols <- c("#de2d26")

# Combine all strategy-level colors into a single named vector
strategy_colors <- c(behavior_cols[1:6], clinical_cols[1:5], lca_cols)

# Assign colors to each strategy based on unique strategy labels in the dataset
names(strategy_colors) <- unique(curve2$strategy)

# Visualization of case coverage (Figure 3A)
pdf("results/main_results_coverage_curves_95CI.pdf", width = 10, height = 6)
p_main <- ggplot() +
  
  # uncertainty ribbons
  geom_ribbon(
    data = band,
    aes(x = N, ymin = lwr_band, ymax = upr_band, fill = group),
    alpha = 0.25,
    color = NA
  ) +
  
  # coverage curves
  geom_line(
    data = curve2,
    aes(x = N, y = coverage, color = strategy),
    linewidth = 1.1
  ) +
  geom_point(
    data = curve2,
    aes(x = N, y = coverage, color = strategy),
    size = 2.8
  ) +
  
  # color mapping
  scale_color_manual(values = strategy_colors) +
  scale_fill_manual(values = group_fill, guide = "none") +
  
  # percentage formatting
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Number of vaccine doses",
    y = "Case coverage",
    color = "Strategy"
  ) +
  
  # Clean theme for publication-quality figure
  theme_classic() +
  theme(
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor.y = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.7, "lines"),
    legend.key.width  = unit(1.2, "lines")
  )

print(p_main)
dev.off()

# Spearman-derived adaptive weights

# Visualization of spearman-derived analysis results (Figure 3B)
pdf("results/spearman_CI.pdf", width = 10, height = 6)
p_sa <- ggplot(
  curve_sa_all,
  aes(x = N, y = coverage, color = strategy, fill = strategy)
) +
  
  # 95% bootstrap CI
  geom_ribbon(
    aes(ymin = lwr, ymax = upr),
    alpha = 0.18,
    color = NA
  ) +
  
  # Mean coverage curve
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.8) +
  labs(
    x = "Number of vaccine doses",
    y = "Case coverage",
    color = "Strategy"
  ) +
  
  # Manual color mapping for strategies
  scale_color_manual(values = c(
    "Behavioral risk-based (Spearman-derived weights)" = "#6baed6",
    "Clinical risk-based (Spearman-derived weights)"  = "#74c476",
    "LCA-informed"                          = "#de2d26"
  )) +
  scale_fill_manual(values = c(
    "Behavioral risk-based (Spearman-derived weights)" = "#6baed6",
    "Clinical risk-based (Spearman-derived weights)"  = "#74c476",
    "LCA-informed"                          = "#de2d26"
  ), guide = "none") +
  
  # Format y-axis as percentage
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  
  # Clean theme for publication-quality figure
  theme_classic() +
  theme(
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor.y = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.7, "lines"),
    legend.key.width  = unit(1.2, "lines")
  )

print(p_sa)
dev.off()

# =========================
# 6. Sensitivity Analysis
# =========================

# Recode strategies into broader methodological groups for comparison
curve_sens <- curve_main_sens %>%
  mutate(group = case_when(
    str_detect(strategy, "^Behavioral risk-based") ~ "Behavioral risk-based",
    str_detect(strategy, "^Clinical risk-based")  ~ "Clinical risk-based",
    str_detect(strategy, "^LCA")            ~ "LCA-informed",
    TRUE ~ "Other"
  ))

# Compute uncertainty bands per strategy group
band_sens <- curve_sens %>%
  filter(group %in% c("Behavioral risk-based", "Clinical risk-based", "LCA-informed")) %>%
  group_by(group, N) %>%
  summarise(
    lwr_band = min(lwr, na.rm = TRUE),
    upr_band = max(upr, na.rm = TRUE),
    .groups = "drop"
  )

# Define fill colors for uncertainty ribbons by strategy group
group_fill <- c(
  "Behavioral risk-based" = "#deebf7",  
  "Clinical risk-based"  = "#e5f5e0",  
  "LCA-informed"   = "#fee0d2"   
)

# Behavioral risk-based strategies (gradient blues)
behavior_cols <- c("#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#084594")

# Clinical risk-based strategies (gradient greens)
clinical_cols <- c("#c7e9c0", "#a1d99b", "#74c476", "#41ab5d", "#238b45")

# LCA-informed strategy (single representative red)
lca_cols <- c("#de2d26")

# Combine all strategy-level colors into a single named vector
strategy_colors <- c(behavior_cols[1:6], clinical_cols[1:5], lca_cols)

# Assign colors to each strategy based on unique strategy labels in the dataset
names(strategy_colors) <- unique(curve_sens$strategy)

# Visualization of case coverage (Figure S2A)
pdf("results/sensitivity_results_coverage_curves_CI.pdf", width = 10, height = 6)
p_sens <- ggplot() +
  
  # Uncertainty ribbons
  geom_ribbon(
    data = band_sens,
    aes(x = N, ymin = lwr_band, ymax = upr_band, fill = group),
    alpha = 0.25,
    color = NA
  ) +
  
  # Coverage curves
  geom_line(
    data = curve_sens,
    aes(x = N, y = coverage, color = strategy),
    linewidth = 1.1
  ) +
  geom_point(
    data = curve_sens,
    aes(x = N, y = coverage, color = strategy),
    size = 2.8
  ) +
  
  # Color mapping
  scale_color_manual(values = strategy_colors) +
  scale_fill_manual(values = group_fill, guide = "none") +
  
  # Format y-axis as percentage
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Number of vaccine doses",
    y = "Case coverage",
    color = "Strategy"
  ) +
  
  # Clean theme for publication-quality figure
  theme_classic() +
  theme(
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor.y = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.7, "lines"),
    legend.key.width  = unit(1.2, "lines")
  )
print(p_sens)
dev.off()

# Spearman-derived adaptive weights

# Visualization of sensitivity analysis results (Figure S2B)
pdf("results/sensitivity_spearman_plot_CI.pdf", width = 10, height = 6)
p_sa_sens <- ggplot(
  curve_sa_sens,
  aes(x = N, y = coverage, color = strategy, fill = strategy)
) +
  
  # 95% bootstrap CI
  geom_ribbon(
    aes(ymin = lwr, ymax = upr),
    alpha = 0.08,
    color = NA
  ) +
  
  # Mean coverage curve
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.8) +
  labs(
    x = "Number of vaccine doses",
    y = "Case coverage",
    color = "Strategy"
  ) +
  
  # Manual color mapping for strategies
  scale_color_manual(values = c(
    "Behavioral risk-based (Spearman-derived weights)" = "#6baed6",
    "Clinical risk-based (Spearman-derived weights)"  = "#74c476",
    "LCA-informed"                          = "#de2d26"
  )) +
  scale_fill_manual(values = c(
    "Behavioral risk-based (Spearman-derived weights)" = "#6baed6",
    "Clinical risk-based (Spearman-derived weights)"  = "#74c476",
    "LCA-informed"                          = "#de2d26"
  ), guide = "none") +
  
  # Format y-axis as percentage
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  
  # Clean theme for publication-quality figure
  theme_classic() +
  theme(
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor.y = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.7, "lines"),
    legend.key.width  = unit(1.2, "lines")
  )
print(p_sa_sens)
dev.off()