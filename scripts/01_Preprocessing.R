################################################################################
## SR-SAVI Reproducibility challenge
## Computational reproducibility of Ripley et al 2021
## Authors: Maria Economou, Torsten Rackoll
## Date: 27.08.2024
################################################################################
library(tidyverse)
library(readxl)

#### Load data
df <- read_excel(here::here("data", "Copy of RIC_infarct volume data_updated.xlsx"))

### inspect structure of data
str(df)
colnames(df)

df <- df %>%
  select(-c('Method of infarct measurement (latest method used)', 'Was infarct volume provided or calculated?',
         'If calculated, what figure was volume obtained from?')) %>%
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
  ) %>%  mutate(across(c(
    Anesthesia, Species, Condition, Ref.vol, Type.RIC,
    Time.cond, Limbs, Sex, Stroke.Type, Model)
    , as_factor),
    across(c(
      Mean, Median, Time.RIC
    ), as.numeric
           ),
    across(c(
      n, No.Limbs, No.Session, No.Cycles, Occlusion, Reperfusion, Occlusion.time
    ), as.integer))) %>%# warning: NA's introduced by coercion. 
     fill(RefID, Author, Species
    )


df.wide <- df %>%
  group_by(RefID)
