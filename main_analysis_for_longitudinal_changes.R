library('readxl')
library("tidyverse")
# new packages!
library(reshape)
library(lme4)
library(nlme)
library(broom.mixed)
library(corrplot)

# Load data
igtd <- read_excel("~/Downloads/Engage_IGT_NAN_2021.xlsx")
wcstd <- read_excel("~/Downloads/Engage_WCST_NAN_2021.xlsx")
data <- read_excel("~/Downloads/Engage_Raw.xlsx", sheet = "DATASET")
# we first need to clean the data a little bit 
table(data$timepoint)
# for this project we are gonna keep only up to wk8

data <- data %>% 
  filter(timepoint %in% c("baseline","wk 2","wk 4","wk 6","wk 8","wk 9")) %>%
  mutate(week = case_when(
    timepoint == "baseline" ~ 0,
    timepoint == "wk 2" ~ 2,
    timepoint == "wk 4" ~ 4,
    timepoint == "wk 6" ~ 6,
    timepoint == "wk 8" ~ 8,
    timepoint == "wk 9" ~ 9
  ))


# Clean data ----
# Change to data frames
# Yiqun: this is not needed since read excel already converted the data into data frames
#data.df <- as.data.frame(data)
#igtd.df <- as.data.frame(igtd)
#wcstd.df <- as.data.frame(wcstd)

# Merge datasets 
merged1 <- merge(data, igtd, by.x = c("patient_id", "timepoint"), 
                 by.y = c("Patient_ID", "Interview_ID"), all = TRUE)

merged.f <- merge(merged1, wcstd, by.x = c("patient_id", "timepoint"), 
                  by.y = c("patient_id", "interview_id"), all = TRUE)

# Change any 9999 values to NA
merged.f[merged.f==9999] <- NA

# Delete MMSE (lots of missing data)
# we don't have to delete all of them since we are utilizing all the data points now!
# df_analysis <- subset(merged.f, select = -c(mmsetotal))

# we first probably want to fill in the missing data???
# visualize analysis
df_analysis <- merged.f %>%
  arrange(patient_id, timepoint)

# let's say that we are focused on HVLT delayed recall at baseline 
df_analysis$HVLTDelayedRecall
df_analysis$Ham24tot # let's say it's ham 24 tot
summary(df_analysis$Ham24tot)
summary(df_analysis$HVLTDelayedRecall)


visual_temp <- df_analysis %>% 
  dplyr::select(patient_id, HVLTDelayedRecall)

baseline_temp <- df_analysis %>% 
  filter(week == 0) %>%
  select(patient_id, HVLTDelayedRecall) %>%
  mutate(baseline_HVLTDelayedRecall = HVLTDelayedRecall) %>%
  select(-HVLTDelayedRecall)

head(baseline_temp)

df_analysis <- df_analysis %>% 
  left_join(baseline_temp, by = "patient_id")

df_analysis %>% select(patient_id, timepoint, baseline_HVLTDelayedRecall)
df_analysis$baseline_HVLTDelayedRecall

df_analysis <- df_analysis %>%
  mutate(Treatment = case_when(
    Treatment %in% c("Engage", "ENGAGE") ~ "ENGAGE",
    TRUE ~ as.character(Treatment)
  ))
# overall spaghetti plot
ggplot(df_analysis, 
       aes(x = week, y=Ham24tot, group=patient_id)) +
  geom_line(alpha=0.5) +
  xlab("Weeks post-baseline") + 
  ylab("Ham 24 total score") +
  theme_classic()

summary(df_analysis$baseline_HVLTDelayedRecall)
# numbers are for demo only!!
df_analysis <- df_analysis %>%
  mutate(
    bl_hvltDR_level = case_when(
      is.na(baseline_HVLTDelayedRecall) ~ 'missing',
      baseline_HVLTDelayedRecall <= 7 ~ 'low',
      (baseline_HVLTDelayedRecall > 7) & (baseline_HVLTDelayedRecall <=10) ~ 'medium',
      baseline_HVLTDelayedRecall > 10 ~'high'
    )
  )
# overall spaghetti plot, grouped by baseline_HVLTDelayedRecall?
ggplot(df_analysis, aes(x = week, y=Ham24tot, group=patient_id, colour=bl_hvltDR_level)) +
  geom_line(alpha=0.5) +
  xlab("Weeks post-baseline") + 
  ylab("Ham 24 total score") +
  theme_classic()


ggplot(df_analysis, aes(x = week, y=Ham24tot, group=patient_id)) +
  geom_line(alpha=0.5) +
  xlab("Weeks post-baseline") + 
  ylab("Ham 24 total score") +
  theme_classic()+
  facet_grid(~bl_hvltDR_level)

## look at correlations too 


df_analysis_wide <- reshape(df_analysis %>% 
                              select(c(patient_id, week, Ham24tot)), 
                            timevar = c("week"), idvar = "patient_id",
                        v.names = c("Ham24tot"),
                        direction = "wide")

cor_mat <- cor(df_analysis_wide[,c("Ham24tot.0", "Ham24tot.2", "Ham24tot.4",
                                   "Ham24tot.6","Ham24tot.8","Ham24tot.9")],
               use="pairwise.complete.obs")

cor_mat

table(df_analysis[!is.na(df_analysis$baseline_HVLTDelayedRecall)&
                    !(duplicated(df_analysis$patient_id)),'Gender'])

summary(df_analysis[!is.na(df_analysis$baseline_HVLTDelayedRecall)&
                      !(duplicated(df_analysis$patient_id)),'Age'])

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA")) 
corrplot(cor_mat, method="color", col=col(200),type="upper", 
         addCoef.col = "black",
         tl.col="black", tl.srt=45)

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
cognitive_measures <- c("Stroop3", "DSB", "Total_Money", "TOTALErr", "HVLTDelayedRecall")
for (measure in cognitive_measures){
  df_analysis <- create_baseline(df_analysis, measure)
}

# HVLTDelayedRecall
lmm_adjusted_hvlt <- lme(fixed = whodastotpro ~ baseline_HVLTDelayedRecall + week + Age + as.factor(Race) + as.factor(Gender) + Treatment,
                         random = ~1 + week | patient_id,
                         data = df_analysis,
                         na.action = na.omit)

print(summary(lmm_adjusted_hvlt))
print(tidy(lmm_adjusted_hvlt, conf.int = TRUE))

# DSB
lmm_adjusted_dsb <- lme(fixed = whodastotpro ~ baseline_DSB + week + Age + as.factor(Race) + as.factor(Gender) + Treatment,
                        random = ~1 + week | patient_id,
                        data = df_analysis,
                        na.action = na.omit)

print(summary(lmm_adjusted_dsb))
print(tidy(lmm_adjusted_dsb, conf.int = TRUE))

# Total_Money
lmm_adjusted_money <- lme(fixed = whodastotpro ~ baseline_Total_Money + week + Age + as.factor(Race) + as.factor(Gender) + Treatment,
                          random = ~1 + week | patient_id,
                          data = df_analysis,
                          na.action = na.omit)

print(summary(lmm_adjusted_money))
print(tidy(lmm_adjusted_money, conf.int = TRUE))

# TOTALErr
lmm_adjusted_totalerr <- lme(fixed = whodastotpro ~ baseline_TOTALErr + week + Age + as.factor(Race) + as.factor(Gender) + Treatment,
                             random = ~1 + week | patient_id,
                             data = df_analysis,
                             na.action = na.omit)

print(summary(lmm_adjusted_totalerr))
print(tidy(lmm_adjusted_totalerr, conf.int = TRUE))

# Stroop3
lmm_adjusted_stroop3 <- lme(fixed = whodastotpro ~ baseline_Stroop3 + week + Age + as.factor(Race) + as.factor(Gender) + Treatment,
                            random = ~1 + week | patient_id,
                            data = df_analysis,
                            na.action = na.omit)

print(summary(lmm_adjusted_stroop3))
print(tidy(lmm_adjusted_stroop3, conf.int = TRUE))


# List of cognitive measures
cognitive_measures <- c("HVLTDelayedRecall", "DSB", "Total_Money", "TOTALErr", "Stroop3")

# Function to run LMM for a given cognitive measure
run_lmm <- function(measure, data) {
  formula <- as.formula(paste("Ham24tot ~ baseline_", measure, 
                              " + week + Age + as.factor(Race) + as.factor(Gender) + Treatment", 
                              sep = ""))
  
  model <- lme(fixed = formula,
               random = ~1 + week | patient_id,
               data = data,
               na.action = na.omit)
  
  print(paste("Results for", measure))
  print(summary(model))
  print(tidy(model, conf.int = TRUE))
  
  return(model)
}

# Run LMM for each cognitive measure
lmm_results <- lapply(cognitive_measures, run_lmm, data = df_analysis)
summary(lmm_results$Stroop3)
summary(lmm_results$DSB)
summary(lmm_results$Total_Money)
summary(lmm_results$TOTALErr)
summary(lmm_results$HVLTDelayedRecall)


# Name the results list
names(lmm_results) <- cognitive_measures

