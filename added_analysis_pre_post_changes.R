library(readxl)
library(tidyverse)
library(reshape)
library(lme4)
library(nlme)
library(broom.mixed)
library(corrplot)
library(table1)

# Load data
igtd <- read_excel("~/Downloads/Engage_IGT_NAN_2021.xlsx")
wcstd <- read_excel("~/Downloads/Engage_WCST_NAN_2021.xlsx")
data <- read_excel("~/Downloads/Engage_Raw.xlsx", sheet = "DATASET")


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
head(df_analysis, 10)


baseline_temp <- df_analysis %>% 
  filter(week == 0) %>%
  select(patient_id, HVLTDelayedRecall) %>%
  mutate(baseline_HVLTDelayedRecall = HVLTDelayedRecall) %>%
  select(-HVLTDelayedRecall)

df_analysis <- df_analysis %>% 
  left_join(baseline_temp, by = "patient_id")

baseline_dsb <- df_analysis %>% 
  filter(week == 0) %>%
  select(patient_id, DSB) %>%
  mutate(baseline_DSB = DSB) %>%
  select(-DSB)

df_analysis <- df_analysis %>% 
  left_join(baseline_dsb, by = "patient_id")

#### 
baseline_Total_Money <- df_analysis %>% 
  filter(week == 0) %>%
  select(patient_id, Total_Money) %>%
  mutate(baseline_Total_Money = Total_Money) %>%
  select(-Total_Money)

df_analysis <- df_analysis %>% 
  left_join(baseline_Total_Money, by = "patient_id")

#### 
baseline_TOTALErr <- df_analysis %>% 
  filter(week == 0) %>%
  select(patient_id, TOTALErr) %>%
  mutate(baseline_TOTALErr = TOTALErr) %>%
  select(-TOTALErr)

df_analysis <- df_analysis %>% 
  left_join(baseline_TOTALErr, by = "patient_id")

#### 
baseline_Stroop3 <- df_analysis %>% 
  filter(week == 0) %>%
  select(patient_id, Stroop3) %>%
  mutate(baseline_Stroop3 = Stroop3) %>%
  select(-Stroop3)

df_analysis <- df_analysis %>% 
  left_join(baseline_Stroop3, by = "patient_id")

#### 
baseline_Ham24tot <- df_analysis %>% 
  filter(week == 0) %>%
  select(patient_id, Ham24tot) %>%
  mutate(baseline_Ham24tot = Ham24tot) %>%
  select(-Ham24tot)

df_analysis <- df_analysis %>% 
  left_join(baseline_Ham24tot, by = "patient_id")

#### 
baseline_whodastotpro <- df_analysis %>% 
  filter(week == 0) %>%
  select(patient_id, whodastotpro) %>%
  mutate(baseline_whodastotpro = whodastotpro) %>%
  select(-whodastotpro)

df_analysis <- df_analysis %>% 
  left_join(baseline_whodastotpro, by = "patient_id")


# pick out the participants who could be able to run for all these regressions
df_analysis_complete <- df_analysis %>% 
  drop_na(any_of(c("baseline_TOTALErr", "baseline_Stroop3",
                   "baseline_Total_Money", "baseline_DSB", 
                   "baseline_HVLTDelayedRecall",
                   "baseline_whodastotpro",
                   "baseline_Ham24tot")))
length(unique(df_analysis_complete$patient_id))



df_analysis_complete$Treatment_upper <- factor(toupper(df_analysis_complete$Treatment))
df_analysis_complete_week9 <- df_analysis_complete %>% filter(week==9)


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


