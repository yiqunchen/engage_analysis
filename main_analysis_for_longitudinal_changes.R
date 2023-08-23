
library('readxl')
library("tidyverse")
# new packages!
library(reshape)
library(lme4)
library(nlme)
library(broom.mixed)
library(corrplot)

# Load data
#data <- read_excel("~/Downloads/Engage_Data_NAN_2021.xlsx")
igtd <- read_excel("~/Downloads/Engage_IGT_NAN_2021.xlsx")
wcstd <- read_excel("~/Downloads/Engage_WCST_NAN_2021.xlsx")
data <- read_excel("~/Downloads/Engage_Raw.xlsx", sheet = "DATASET")
# we first need to clean the data a little bit 
data
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
head(df_analysis, 10)

# let's say that we are focused on HVLT delayed recall at baseline 
df_analysis$HVLTDelayedRecall
df_analysis$Ham24tot # let's say it's ham 24 tot
summary(df_analysis$Ham24tot)
summary(df_analysis$HVLTDelayedRecall)


visual_temp <- df_analysis %>% 
  select(patient_id, HVLTDelayedRecall)

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


#### restrict to a population of choice



lmm_unadjust_rand_int <- lme(fixed=Ham24tot ~ baseline_HVLTDelayedRecall+week, 
                               random=reStruct(~1 | patient_id, REML=TRUE), 
                               data=df_analysis,
                               subset= !is.na(baseline_HVLTDelayedRecall),
                               na.action=na.omit)

summary(lmm_unadjust_rand_int)
tidy(lmm_unadjust_rand_int, conf.int = TRUE)




lmm_unadjust_rand_slope <- lme(fixed=Ham24tot ~ baseline_HVLTDelayedRecall+week, 
                    random=reStruct(~1+week | patient_id, REML=TRUE), 
                    data=df_analysis,
                    subset= !is.na(baseline_HVLTDelayedRecall),
                    na.action=na.omit)

summary(lmm_unadjust_rand_slope)
tidy(lmm_unadjust_rand_slope, conf.int = TRUE)



lmm_adjusted_rand_slope <- lme(fixed=Ham24tot ~ baseline_HVLTDelayedRecall+week+
                                 Age+
                                 as.factor(Race)+as.factor(Gender)+Treatment, 
                             random=reStruct(~1+week | patient_id, REML=TRUE), 
                             data=df_analysis,
                             subset= !is.na(baseline_HVLTDelayedRecall),
                             na.action=na.omit)

summary(lmm_adjusted_rand_slope)
tidy(lmm_adjusted_rand_slope, conf.int = TRUE)


baseline_dsb <- df_analysis %>% 
  filter(week == 0) %>%
  select(patient_id, DSB) %>%
  mutate(baseline_DSB = DSB) %>%
  select(-DSB)

df_analysis <- df_analysis %>% 
  left_join(baseline_dsb, by = "patient_id")


lmm_unadjust_rand_slope <- lme(fixed=Ham24tot ~ baseline_DSB+week, 
                               random=reStruct(~1+week | patient_id, REML=TRUE), 
                               data=df_analysis,
                               subset= !is.na(baseline_DSB),
                               na.action=na.omit)

summary(lmm_unadjust_rand_slope)
tidy(lmm_unadjust_rand_slope, conf.int = TRUE)

lmm_adjusted_rand_slope <- lme(fixed=Ham24tot ~ baseline_DSB+week+
                                 Age+
                                 as.factor(Race)+as.factor(Gender)+Treatment, 
                               random=reStruct(~1+week | patient_id, REML=TRUE), 
                               data=df_analysis,
                               subset= !is.na(baseline_DSB),
                               na.action=na.omit)

summary(lmm_adjusted_rand_slope)
tidy(lmm_adjusted_rand_slope, conf.int = TRUE)

#### 
baseline_dsb <- df_analysis %>% 
  filter(week == 0) %>%
  select(patient_id, DSB) %>%
  mutate(baseline_DSB = DSB) %>%
  select(-DSB)

df_analysis <- df_analysis %>% 
  left_join(baseline_dsb, by = "patient_id")

lmm_unadjust_rand_slope <- lme(fixed=Ham24tot ~ baseline_DSB+week, 
                               random=reStruct(~1+week | patient_id, REML=TRUE), 
                               data=df_analysis,
                               subset= !is.na(baseline_DSB),
                               na.action=na.omit)

summary(lmm_unadjust_rand_slope)
tidy(lmm_unadjust_rand_slope, conf.int = TRUE)

lmm_adjusted_rand_slope <- lme(fixed=Ham24tot ~ baseline_DSB+week+
                                 Age+
                                 as.factor(Race)+as.factor(Gender)+Treatment, 
                               random=reStruct(~1+week | patient_id, REML=TRUE), 
                               data=df_analysis,
                               subset= !is.na(baseline_DSB),
                               na.action=na.omit)

summary(lmm_adjusted_rand_slope)
tidy(lmm_adjusted_rand_slope, conf.int = TRUE)

#### 
baseline_Total_Money <- df_analysis %>% 
  filter(week == 0) %>%
  select(patient_id, Total_Money) %>%
  mutate(baseline_Total_Money = Total_Money) %>%
  select(-Total_Money)

df_analysis <- df_analysis %>% 
  left_join(baseline_Total_Money, by = "patient_id")

lmm_unadjust_rand_slope <- lme(fixed=Ham24tot ~ baseline_Total_Money+week, 
                               random=reStruct(~1+week | patient_id, REML=TRUE), 
                               data=df_analysis,
                               subset= !is.na(baseline_Total_Money),
                               na.action=na.omit)

summary(lmm_unadjust_rand_slope)
tidy(lmm_unadjust_rand_slope, conf.int = TRUE)

lmm_adjusted_rand_slope <- lme(fixed=Ham24tot ~ baseline_Total_Money+week+
                                 Age+
                                 as.factor(Race)+as.factor(Gender)+Treatment, 
                               random=reStruct(~1+week | patient_id, REML=TRUE), 
                               data=df_analysis,
                               subset= !is.na(baseline_Total_Money),
                               na.action=na.omit)

summary(lmm_adjusted_rand_slope)
tidy(lmm_adjusted_rand_slope, conf.int = TRUE)

#### 
baseline_TOTALErr <- df_analysis %>% 
  filter(week == 0) %>%
  select(patient_id, TOTALErr) %>%
  mutate(baseline_TOTALErr = TOTALErr) %>%
  select(-TOTALErr)

df_analysis <- df_analysis %>% 
  left_join(baseline_TOTALErr, by = "patient_id")

lmm_unadjust_rand_slope <- lme(fixed=Ham24tot ~ baseline_TOTALErr+week, 
                               random=reStruct(~1+week | patient_id, REML=TRUE), 
                               data=df_analysis,
                               subset= !is.na(baseline_TOTALErr),
                               na.action=na.omit)

summary(lmm_unadjust_rand_slope)
tidy(lmm_unadjust_rand_slope, conf.int = TRUE)

lmm_adjusted_rand_slope <- lme(fixed=Ham24tot ~ baseline_TOTALErr+week+
                                 Age+
                                 as.factor(Race)+as.factor(Gender)+Treatment, 
                               random=reStruct(~1+week | patient_id, REML=TRUE), 
                               data=df_analysis,
                               subset= !is.na(baseline_TOTALErr),
                               na.action=na.omit)

summary(lmm_adjusted_rand_slope)
tidy(lmm_adjusted_rand_slope, conf.int = TRUE)


#### 
baseline_Stroop3 <- df_analysis %>% 
  filter(week == 0) %>%
  select(patient_id, Stroop3) %>%
  mutate(baseline_Stroop3 = Stroop3) %>%
  select(-Stroop3)

df_analysis <- df_analysis %>% 
  left_join(baseline_Stroop3, by = "patient_id")

lmm_unadjust_rand_slope <- lme(fixed=Ham24tot ~ baseline_Stroop3+week, 
                               random=reStruct(~1+week | patient_id, REML=TRUE), 
                               data=df_analysis,
                               subset= !is.na(baseline_Stroop3),
                               na.action=na.omit)

summary(lmm_unadjust_rand_slope)
tidy(lmm_unadjust_rand_slope, conf.int = TRUE)

lmm_adjusted_rand_slope <- lme(fixed=Ham24tot ~ baseline_Stroop3+week+
                                 Age+
                                 as.factor(Race)+as.factor(Gender)+Treatment, 
                               random=reStruct(~1+week | patient_id, REML=TRUE), 
                               data=df_analysis,
                               subset= !is.na(baseline_Stroop3),
                               na.action=na.omit)

summary(lmm_adjusted_rand_slope)
tidy(lmm_adjusted_rand_slope, conf.int = TRUE)


 