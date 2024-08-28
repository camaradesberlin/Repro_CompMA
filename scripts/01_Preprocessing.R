################################################################################
## SR-SAVI Reproducibility challenge
## Computational reproducibility of Ripley et al 2021
## Authors: Maria Economou, Torsten Rackoll
## Date: 27.08.2024
################################################################################

# Libraries ---------------------------------------------------------------


library(tidyverse)
library(readxl)
library(metafor)


# Load data ---------------------------------------------------------------


raw_data <- read_excel(here::here("data", "Copy of RIC_infarct volume data_updated.xlsx"))


# Inspect structure of data -----------------------------------------------

str(raw_data)
colnames(raw_data)

df <- raw_data %>%
  select(-c('Method of infarct measurement (latest method used)', 'Was infarct volume provided or calculated?',
         'If calculated, what figure was volume obtained from?', '95% CI', 'IQR')) %>%
  rename(
    "Author" = "First Author (Year)",
    "Condition" = "Control vs RIC",
    "Int.Label" = "Intervention Group",
    "Time.Inf" = "Latest timepoint (days)",
    "Unit" = "Unit of measurement",
    "Ref.vol" = "Reference point for infarct volume (% of total brain/hemisphere)",
    "Type.RIC" = "Type of RIC",
    "Time.cond" = "Pre-conditioning or post-conditioning?",
    "Time.RIC" = "Time of RIC onset relative to onset of stroke (hours)",
    "No.Limbs" = "Number of limbs used for RIC?",
    "Limbs" = "Which limbs were used for RIC?",
    "No.Session" = "Numberof RIC sessions",
    "No.Cycles" = "Number of RIC cycles",
    "Occlusion" = "Duration of occlusion (minutes)",
    "Reperfusion" = "Duration of reperfusion (minutes)",
    "Stroke.Type" = "Stroke type (permanent or reperfusion model)",
    "Model" = "Stroke type",
    "Occlusion.time" = "Occlusion time (min)"
  ) %>% 
  mutate(
    across(c(
      Anesthesia, Species, Condition, Ref.vol, Type.RIC,
      Time.cond, Limbs, Sex, Stroke.Type, Model), 
      as_factor),
    across(c(
      Mean, Median, Time.RIC),
      as.numeric),
    across(c(
      n, No.Limbs, No.Session, No.Cycles, Occlusion, Reperfusion, Occlusion.time), 
      as.integer)) %>% # warning: NA's introduced by coercion. 
     fill(RefID, Author, Species, Stroke.Type, Model, Occlusion.time
    )

# why does column Median contain text (same information as Int.Label) for the Ren et al 2008 study?

var <- c("Anesthesia", "Int.Label", "n", "Mean", "Median", "SD", "SEM",
         "Type.RIC", "Time.cond", "Time.RIC", "No.Limbs", "Limbs", "No.Session",
         "No.Cycles", "Occlusion", "Reperfusion")

check_studies <- df %>% 
  group_by(RefID) %>% 
  tally() %>% 
  filter(n > 2) %>% 
  select(RefID) %>% 
  left_join(df, by = "RefID")

studies_multi_contr <- check_studies %>% 
  group_by(RefID, Condition) %>% 
  tally() %>% 
  pivot_wider(names_from = Condition, values_from = n) %>% 
  filter(Control == 1) %>% 
  select(RefID) %>% 
  left_join(df, by = "RefID") %>% 
  ungroup() %>% 
  group_by(RefID) %>% 
  pivot_wider(names_from = Condition,
              values_from = c(n, Mean, Median, SD, SEM)) %>% 
  fill(n_Control, Mean_Control, Median_Control, SD_Control, SEM_Control)

df.wide <- df %>%
  group_by(RefID) %>%
  pivot_wider(names_from = Condition,
              values_from = all_of(var))

df.wide_es <- escalc(measure = "SMD",   # Specifies the type of effect size (Standardized Mean Difference)
                     m1i = Mean_RIC,         # Mean of treatment group
                     sd1i = SD_RIC,       # SD of treatment group
                     n1i = n_RIC,         # Sample size of treatment group
                     m2i = Mean_Control,         # Mean of control group
                     sd2i = SD_Control,       # SD of control group
                     n2i = n_Control,         # Sample size of control group
                     data = df_wide)

save(df.wide_es, file = here::here("data", "RIC_Ripley_wide_2024-08-27.RData"))
