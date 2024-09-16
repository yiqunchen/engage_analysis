library(readxl)
library(tidyverse)
library(reshape)
library(lme4)
library(nlme)
library(broom)
library(broom.mixed)
library(corrplot)
library(table1)

# i can check that 
# HVLTImmediateRecall/HVLTDelayedRecall 
# 
# Load data
igtd <- read_excel("~/Downloads/Engage_IGT_NAN_2021.xlsx")
wcstd <- read_excel("~/Downloads/Engage_WCST_NAN_2021.xlsx")
data <- read_excel("~/Downloads/Engage_Raw.xlsx", sheet = "DATASET")

# Filter data for specific timepoints and create 'week' variable
data <- data %>% 
  filter(timepoint %in% c("baseline", "wk 2", "wk 4", "wk 6", "wk 8", "wk 9")) %>%
  mutate(week = case_when(
    timepoint == "baseline" ~ 0,
    timepoint == "wk 2" ~ 2,
    timepoint == "wk 4" ~ 4,
    timepoint == "wk 6" ~ 6,
    timepoint == "wk 8" ~ 8,
    timepoint == "wk 9" ~ 9
  ))

# Merge datasets
merged1 <- merge(data, igtd, by.x = c("patient_id", "timepoint"), 
                 by.y = c("Patient_ID", "Interview_ID"), all = TRUE)

merged.f <- merge(merged1, wcstd, by.x = c("patient_id", "timepoint"), 
                  by.y = c("patient_id", "interview_id"), all = TRUE)

# Replace 9999 with NA
merged.f[merged.f == 9999] <- NA

# Arrange data for analysis
df_analysis <- merged.f %>%
  arrange(patient_id, timepoint)

# Function to add baseline values for specific columns
add_baseline_values <- function(df, column_name) {
  baseline_temp <- df %>% 
    filter(week == 0) %>%
    select(patient_id, all_of(column_name)) %>%
    rename_with(~ paste0("baseline_", .), all_of(column_name))
  
  df %>% left_join(baseline_temp, by = "patient_id")
}

# Compute the ratio of HVLTDelayedRecall to HVLTImmediateRecall at each timepoint
df_analysis <- df_analysis %>%
  mutate(hvlt_ratio = HVLTDelayedRecall / HVLTImmediateRecall)
# Columns to compute baseline values for
baseline_columns <- c("HVLTDelayedRecall", "HVLTImmediateRecall", "DSB", "Total_Money", 
                      "TOTALErr", "Stroop3", "Ham24tot", "whodastotpro", "hvlt_ratio")

# Apply the function for each baseline column
for (col in baseline_columns) {
  df_analysis <- add_baseline_values(df_analysis, col)
}


analyze_inclusion_exclusion_diff <- function(df_input, 
                                             analysis_rows = c("baseline_HVLTDelayedRecall",
                                                               "baseline_HVLTImmediateRecall",
                                                               "baseline_hvlt_ratio",
                                                               "baseline_DSB", 
                                                               "baseline_Total_Money", 
                                                               "baseline_TOTALErr",
                                                               "baseline_Stroop3", 
                                                               "baseline_Ham24tot", 
                                                               "baseline_whodastotpro")) {
  # Create a function to check if a patient has complete baseline information
  has_complete_info <- function(id) {
    patient_data <- df_input %>% 
      filter(patient_id == id)
    
    baseline_row <- patient_data %>% 
      filter(timepoint == "baseline")
    
    # Ensure that baseline_row has the required columns
    baseline_complete <- all(!is.na(baseline_row[analysis_rows]))
    
    return(baseline_complete)
  }
  
  # Group patients by completeness of baseline data
  patient_groups <- df_input %>%
    distinct(patient_id) %>%
    mutate(group = if_else(map_lgl(patient_id, has_complete_info), 
                           "Complete", "Incomplete"))
  
  # Join the groups back to the original data
  data_with_groups <- df_input %>%
    left_join(patient_groups, by = "patient_id")
  
  # Calculate averages, standard deviations, and counts for both groups at week 0 and week 9
  results <- data_with_groups %>%
    filter(week %in% c(0, 9)) %>%
    group_by(group, week) %>%
    summarise(
      count = n(),
      across(c(Ham24tot, whodastotpro, HVLTDelayedRecall, Stroop3, 
               mmsetotal, HVLTImmediateRecall, hvlt_ratio),
             list(avg = ~mean(., na.rm = TRUE), 
                  sd = ~sd(., na.rm = TRUE)),
             .names = "{.col}_{.fn}"),
      .groups = "drop"
    )
  
  # Calculate overall averages and standard deviations
  results_overall <- data_with_groups %>%
    filter(week %in% c(0, 9)) %>%
    group_by(week) %>%
    summarise(
      count = n(),
      across(c(Ham24tot, whodastotpro, HVLTDelayedRecall, Stroop3, 
               mmsetotal, HVLTImmediateRecall, hvlt_ratio),
             list(avg = ~mean(., na.rm = TRUE), 
                  sd = ~sd(., na.rm = TRUE)),
             .names = "{.col}_{.fn}"),
      .groups = "drop"
    )
  
  return(list(results = results,
              results_overall = results_overall,
              data_with_groups = data_with_groups))
}

# Example usage
supp_tab_1 <- analyze_inclusion_exclusion_diff(df_analysis)
supp_tab_1$results
supp_tab_1$results_overall

data_with_groups <- supp_tab_1$data_with_groups
data_with_groups$Treatment <- toupper(data_with_groups$Treatment)
# List of metrics to test
metrics <- c("Ham24tot", "whodastotpro", "HVLTDelayedRecall", "Stroop3", "mmsetotal",
             "hvlt_ratio", "HVLTImmediateRecall")

# Initialize an empty list to store results
results_list <- list()

# Perform t-tests for each metric at week 0 and week 9
for (w in c(0, 9)) {
  for (m in metrics) {
    if (w == 9 & m == 'mmsetotal') {
      cat("Skipping mmsetotal at week 9\n")
      next
    }
    # Filter data for the current week
    week_data <- supp_tab_1$data_with_groups %>% filter(week == w)
    
    # Perform t-test
    t_result <- tryCatch({
      t.test(week_data[[m]] ~ week_data$group)
    }, error = function(e) {
      cat("Error in t-test for", m, "at week", w, ":", conditionMessage(e), "\n")
      return(NULL)
    })
    
    # If t-test was successful, store results
    if (!is.null(t_result)) {
      tidy_result <- tidy(t_result) %>%
        mutate(
          week = w,
          metric = m
        )
      results_list[[length(results_list) + 1]] <- tidy_result
    }
  }
}

# Combine all results into a single dataframe
all_results <- bind_rows(results_list)


# pick out the participants who could be able to run for all these regressions
df_analysis_complete <- df_analysis %>% 
  drop_na(any_of(c("baseline_TOTALErr", "baseline_Stroop3",
                   "baseline_Total_Money", "baseline_DSB", 
                   "baseline_HVLTDelayedRecall",
                   "baseline_HVLTImmediateRecall",
                   "baseline_whodastotpro",
                   "baseline_Ham24tot")))

df_analysis_complete$Treatment_upper <- factor(toupper(df_analysis_complete$Treatment))
df_analysis_complete_supp_table_1 <- df_analysis_complete %>% 
  filter(week==9|week==0)
df_analysis_complete_week9 <- df_analysis_complete %>% filter(week==9)

# Function to run LMM for change scores
run_change_score_lmm <- function(df, outcome_var) {
  formula <- as.formula(paste0("I(", outcome_var, " - baseline_", outcome_var, ") ~ ",
                               "relevel(Treatment_upper,ref = 'PST') +",
                               "I(Ham24tot - baseline_Ham24tot) + ",
                               "I(whodastotpro - baseline_whodastotpro)"))
  
  lmm <- lme(fixed = formula,
             random = reStruct(~1 | patient_id, REML = TRUE),
             na.action = na.omit,
             data = df)
  
  list(summary = summary(lmm),
       tidy = tidy(lmm, conf.int = TRUE, effects = "fixed"))
}

# Ensure baseline variables exist
create_baseline <- function(df, var_name) {
  baseline_var <- paste0("baseline_", var_name)
  
  if (baseline_var %in% names(df)) {
    return(df)
  }
  
  baseline <- df %>% 
    filter(week == 0) %>%
    select(patient_id, !!sym(var_name)) %>%
    rename_with(~paste0("baseline_", .), -patient_id)
  
  df %>% left_join(baseline, by = "patient_id")
}

# List of all measures to analyze
measures <- c("Stroop3", "DSB", "Total_Money", "TOTALErr", "HVLTDelayedRecall",
              "hvlt_ratio")

# Ensure all baseline variables exist
for (measure in c(measures, "Ham24tot", "whodastotpro")) {
  df_analysis_complete_week9 <- create_baseline(df_analysis_complete_week9, measure)
}

# Run analysis for each measure
results <- list()

for (measure in measures) {
  results[[measure]] <- run_change_score_lmm(df_analysis_complete_week9, measure)
  
  cat("\n\nResults for", measure, ":\n")
  #print(results[[measure]]$summary)
  #print(results[[measure]]$tidy)
}

library(dplyr)
library(tidyr)
library(Hmisc)
library(corrplot)
library(kableExtra)

# Function to create baseline variable if it doesn't exist
create_baseline <- function(df, var_name) {
  baseline_var <- paste0("baseline_", var_name)
  
  if (baseline_var %in% names(df)) {
    return(df)
  }
  
  baseline <- df %>% 
    filter(week == 0) %>%
    select(patient_id, !!sym(var_name)) %>%
    rename_with(~paste0("baseline_", .), -patient_id)
  
  df %>% left_join(baseline, by = "patient_id")
}

# Function to compute change scores
compute_change_scores <- function(df, measures) {
  for (measure in measures) {
    change_var <- paste0(measure, "_change")
    df <- df %>%
      mutate(!!change_var := !!sym(measure) - !!sym(paste0("baseline_", measure)))
  }
  df
}

# List of measures
cognitive_measures <- c("Stroop3", "DSB", "Total_Money", "TOTALErr", "HVLTDelayedRecall",
                        "hvlt_ratio")
clinical_measures <- c("Ham24tot", "whodastotpro")
all_measures <- c(cognitive_measures, clinical_measures)

# Ensure all baseline variables exist and compute change scores
for (measure in all_measures) {
  df_analysis_complete_week9 <- create_baseline(df_analysis_complete_week9, measure)
}

df_analysis_complete_week9 <- compute_change_scores(df_analysis_complete_week9, all_measures)

# Function to compute correlation matrix with p-values
compute_correlation_matrix <- function(df, measures) {
  change_vars <- paste0(measures, "_change")
  cor_matrix <- Hmisc::rcorr(as.matrix(df[, change_vars]), type = "pearson")
  return(cor_matrix)
}

# Function to create formatted correlation table
create_correlation_table <- function(cor_matrix, measures, title) {
  change_vars <- paste0(measures, "_change")
  cor_table <- data.frame(
    Measure1 = rep(change_vars, each = length(change_vars)),
    Measure2 = rep(change_vars, times = length(change_vars)),
    Correlation = as.vector(cor_matrix$r),
    P_value = as.vector(cor_matrix$P)
  ) %>%
    filter(Measure1 < Measure2) %>%
    mutate(
      Correlation = round(Correlation, 3),
      P_value = round(P_value, 3),
      Significance = case_when(
        P_value < 0.001 ~ "***",
        P_value < 0.01 ~ "**",
        P_value < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  kable(cor_table, caption = paste("Correlation Table -", title)) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
    save_kable(file = paste0(title, "_correlation_table.html"))
  
  return(cor_table)
}

# Compute cross-correlations between cognitive and clinical measures
cross_cor <- Hmisc::rcorr(as.matrix(df_analysis_complete_week9[, c(paste0(cognitive_measures, "_change"), 
                                                            paste0(clinical_measures, "_change"))]), 
                   type = "pearson")


# Create cross-correlation table
cross_table <- data.frame(
  Cognitive_Measure = rep(paste0(cognitive_measures, "_change"), each = length(clinical_measures)),
  Clinical_Measure = rep(paste0(clinical_measures, "_change"), times = length(cognitive_measures)),
  Correlation = as.vector(cross_cor$r[paste0(cognitive_measures, "_change"), paste0(clinical_measures, "_change")]),
  P_value = as.vector(cross_cor$P[paste0(cognitive_measures, "_change"), paste0(clinical_measures, "_change")])
) %>%
  mutate(
    Correlation = round(Correlation, 3),
    P_value = round(P_value, 3),
    Significance = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01 ~ "**",
      P_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

library(ppcor)

# List of cognitive measures and WHO-DAS
measures <- c("HVLTDelayedRecall", "DSB", "Total_Money", "TOTALErr", "Stroop3", "WHO_DAS",
              "hvlt_ratio")

# Calculate change scores for all measures and HAM-D
# Function to run partial correlation
run_partial_correlation <- function(x, y, z) {
  result <- pcor.test(x, y, z)
  data.frame(
    measure = deparse(substitute(x)),
    estimate = result$estimate,
    p_value = result$p.value,
    lower_ci = result$conf.int[1],
    upper_ci = result$conf.int[2]
  )
}

# Run partial correlations
partial_list <- list()
for (measure in cognitive_measures){
  df_temp <- df_analysis_complete_week9 %>%
    dplyr::select(
      paste0(measure, "_change"),
      whodastotpro_change,
      Ham24tot_change
    ) %>%
    na.omit()
  partial_list[[measure]] <- pcor.test(
    df_temp[[paste0(measure, "_change")]],
    df_temp$whodastotpro_change,
    df_temp$Ham24tot_change
  )
}
partial_correlations <- map_dfr(
  cognitive_measures,
  ~run_partial_correlation(
    df_analysis_complete_week9[[paste0(.x, "_change")]],
    df_analysis_complete_week9$whodastotpro_change,
    df_analysis_complete_week9$Ham24tot_change
  )
)

# Print results
print(partial_correlations)


# whether changes in stroop 3 can be predicted by the baseline score
# stroops
lm_stroop_3_ham24 <- lme(Stroop3 ~ baseline_Stroop3+
                     relevel(Treatment_upper,ref = 'PST')+I(Ham24tot-baseline_Ham24tot) , 
                     random=reStruct(~1 | patient_id, REML=TRUE),
                     na.action=na.omit,
                   data=df_analysis_complete_week9)

summary(lm_stroop_3_ham24)
tidy(lm_stroop_3_ham24, conf.int = T,effects = "fixed")


lm_stroop_3_whodas <- lme(Stroop3 ~ baseline_Stroop3+
                    relevel(Treatment_upper,ref = 'PST')+I(whodastotpro-baseline_whodastotpro) , 
                    random=reStruct(~1 | patient_id, REML=TRUE),
                    na.action=na.omit,
                  data=df_analysis_complete_week9)

summary(lm_stroop_3_whodas)
tidy(lm_stroop_3_whodas, conf.int = T,effects = "fixed")

# HVLTDelayedRecall
lm_hvlt_ham <- lme(fixed=HVLTDelayedRecall ~ baseline_HVLTDelayedRecall+
      relevel(Treatment_upper,ref = 'PST')+
      I(Ham24tot-baseline_Ham24tot) , 
    random=reStruct(~1 | patient_id, REML=TRUE),
    data=df_analysis_complete_week9,
    na.action=na.omit)

summary(lm_hvlt_ham)
tidy(lm_hvlt_ham, conf.int = T)

lm_HVLT_who <- lme(fixed=HVLTDelayedRecall ~ baseline_HVLTDelayedRecall+
                           relevel(Treatment_upper,ref = 'PST')+
                             I(whodastotpro-baseline_whodastotpro) , 
                           random=reStruct(~1 | patient_id, REML=TRUE),
                         data=df_analysis_complete_week9,
                         na.action=na.omit)

summary(lm_HVLT_who)
tidy(lm_HVLT_who, conf.int = T,effects = "fixed")

# baseline_DSB
lm_dsb_ham <- lme(fixed=DSB ~ baseline_DSB+
                     relevel(Treatment_upper,ref = 'PST')+
                     I(Ham24tot-baseline_Ham24tot) , 
                   random=reStruct(~1 | patient_id, REML=TRUE),
                   data=df_analysis_complete_week9,
                   na.action=na.omit)

tidy(lm_dsb_ham, conf.int = T,effects = "fixed")

summary(lm_dsb_ham)

lm_dsb_who <- lme(fixed=DSB ~ baseline_DSB+
                     relevel(Treatment_upper,ref = 'PST')+
                     I(whodastotpro-baseline_whodastotpro) , 
                   random=reStruct(~1 | patient_id, REML=TRUE),
                   data=df_analysis_complete_week9,
                   na.action=na.omit)

summary(lm_dsb_who)
tidy(lm_dsb_who, conf.int = T,effects = "fixed")


# baseline_wcst
lm_wcst_ham <- lme(fixed=TOTALErr ~ baseline_TOTALErr+
                    relevel(Treatment_upper,ref = 'PST')+
                    I(Ham24tot-baseline_Ham24tot) , 
                  random=reStruct(~1 | patient_id, REML=TRUE),
                  data=df_analysis_complete_week9,
                  na.action=na.omit)

summary(lm_wcst_ham)
tidy(lm_wcst_ham, conf.int = T,effects = "fixed")


lm_wcst_who <- lme(fixed=TOTALErr ~ baseline_TOTALErr+
                    relevel(Treatment_upper,ref = 'PST')+
                    I(whodastotpro-baseline_whodastotpro) , 
                  random=reStruct(~1 | patient_id, REML=TRUE),
                  data=df_analysis_complete_week9,
                  na.action=na.omit)

summary(lm_wcst_who)
tidy(lm_wcst_who, conf.int = T,effects = "fixed")

# baseline_igt
lm_igt_ham <- lme(fixed=Total_Money ~ baseline_Total_Money+
                    relevel(Treatment_upper,ref = 'PST')+
                    I(Ham24tot-baseline_Ham24tot) , 
                  random=reStruct(~1 | patient_id, REML=TRUE),
                  data=df_analysis_complete_week9,
                  na.action=na.omit)

summary(lm_igt_ham)
tidy(lm_igt_ham, conf.int = T,effects = "fixed")

lm_igt_who <- lme(fixed=Total_Money ~ baseline_Total_Money+
                    relevel(Treatment_upper,ref = 'PST')+
                    I(whodastotpro-baseline_whodastotpro) , 
                  random=reStruct(~1 | patient_id, REML=TRUE),
                  data=df_analysis_complete_week9,
                  na.action=na.omit)

summary(lm_igt_who)
tidy(lm_igt_who, conf.int = T,effects = "fixed")


# HVLTDelayedRecall
lm_hvlt_ratio_ham <- lme(fixed=hvlt_ratio ~ baseline_hvlt_ratio+
                     relevel(Treatment_upper,ref = 'PST')+
                     I(Ham24tot-baseline_Ham24tot) , 
                   random=reStruct(~1 | patient_id, REML=TRUE),
                   data=df_analysis_complete_week9,
                   na.action=na.omit)

summary(lm_hvlt_ratio_ham)
tidy(lm_hvlt_ratio_ham, conf.int = T)

lm_hvlt_ratio_who <- lme(fixed=hvlt_ratio ~ baseline_hvlt_ratio+
                           relevel(Treatment_upper,ref = 'PST')+
                           I(whodastotpro-baseline_whodastotpro) , 
                         random=reStruct(~1 | patient_id, REML=TRUE),
                         data=df_analysis_complete_week9,
                         na.action=na.omit)

summary(lm_hvlt_ratio_who)
tidy(lm_hvlt_ratio_who, conf.int = T)


lm_HVLT_who <- lme(fixed=HVLTDelayedRecall ~ baseline_HVLTDelayedRecall+
                     relevel(Treatment_upper,ref = 'PST')+
                     I(whodastotpro-baseline_whodastotpro) , 
                   random=reStruct(~1 | patient_id, REML=TRUE),
                   data=df_analysis_complete_week9,
                   na.action=na.omit)

summary(lm_HVLT_who)
tidy(lm_HVLT_who, conf.int = T,effects = "fixed")


df_analysis$Treatment <- toupper(df_analysis$Treatment)
# HVLTDelayedRecall
lmm_adjusted_hvlt_baseline <- lme(fixed = whodastotpro ~ baseline_hvlt_ratio + week + 
                           Age + as.factor(Race) + as.factor(Gender) + Treatment,
                         random = ~1 + week | patient_id,
                         data = df_analysis,
                         na.action = na.omit)

print(summary(lmm_adjusted_hvlt_baseline))
print(tidy(lmm_adjusted_hvlt_baseline, conf.int = TRUE))

lmm_adjusted_hvlt_baseline_ham <- lme(fixed = Ham24tot ~ baseline_hvlt_ratio + week + 
                                    Age + as.factor(Race) + as.factor(Gender) + Treatment,
                                  random = ~1 + week | patient_id,
                                  data = df_analysis,
                                  na.action = na.omit)

print(summary(lmm_adjusted_hvlt_baseline_ham))
print(tidy(lmm_adjusted_hvlt_baseline_ham, conf.int = TRUE))


# List of baseline variables
baseline_vars <- c("baseline_hvlt_ratio", "baseline_HVLTDelayedRecall", "baseline_HVLTImmediateRecall",
                   "baseline_DSB", "baseline_Total_Money", "baseline_TOTALErr", "baseline_Stroop3")

# List of outcomes
outcomes <- c("ham_response", "ham_retention")

# Initialize an empty dataframe to store results
results_table <- data.frame()

# Loop through each baseline variable and outcome to run the glm models
for (outcome in outcomes) {
  for (baseline_var in baseline_vars) {
    
    # Define the formula for the GLM
    formula <- as.formula(paste(outcome, "~", baseline_var, "+ Age + as.factor(Race) + as.factor(Gender) + Treatment"))
    
    # Fit the GLM model
    model <- glm(formula, data = df_analysis_hamt_cat, family = binomial(), na.action = na.omit)
    
    # Extract the odds ratio, p-value, and confidence intervals for the baseline variable
    model_tidy <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(term == baseline_var) %>%
      mutate(outcome = outcome, baseline_var = baseline_var)
    
    # Append the results to the results_table
    results_table <- bind_rows(results_table, model_tidy)
  }
}

# Pivot the results so that each baseline variable has two sets of results: one for each outcome
results_table_pivot <- results_table %>%
  mutate(OR = round(estimate, 4),
         p.value = round(p.value,2),
         CI = paste0("[", round(conf.low, 4), ", ", round(conf.high, 4), "]")) %>%
  dplyr::select(baseline_var, outcome, OR, p.value, CI) %>%
  pivot_wider(names_from = outcome, values_from = c(OR, p.value, CI),
              names_glue = "{outcome}_{.value}") %>%
  dplyr::rename(`Ham Response OR` = ham_response_OR,
         `Ham Response P-value` = ham_response_p.value,
         `Ham Response 95% CI` = ham_response_CI,
         `Ham Retention OR` = ham_retention_OR,
         `Ham Retention P-value` = ham_retention_p.value,
         `Ham Retention 95% CI` = ham_retention_CI)


