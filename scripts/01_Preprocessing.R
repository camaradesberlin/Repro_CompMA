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

# str(raw_data)
# colnames(raw_data)

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
  fill(RefID, Author, Species, Stroke.Type, Model, Occlusion.time) %>% 
  mutate(
    across(
      c(Anesthesia, Species, Condition, Ref.vol, Type.RIC,
        Time.cond, Limbs, Sex, Stroke.Type, Model), 
      as_factor),
    # still need to deal with sample size (it's a range sometimes and can't be converted to numeric = NA)
    # now we extract the lower boundary as sample size
    n = ifelse(grepl("-", n), stringr::str_extract(n, "[^-]+"), n), 
    across(
      c(Mean, Median, Time.RIC, n),
      as.numeric)
  ) %>% 
  # fill variable of sex
  group_by(RefID) %>% 
  fill(Sex, .direction = "downup") %>% 
  ungroup() %>% 
  # harmonize sex variable
  mutate(Sex = case_when(Sex %in% "male" ~ "Male",
                         Sex %in% "50% male" ~ "Both",
                         .default = Sex)) %>% 
  mutate(Condition2 = Condition)

# column Median contain text (same information as Int.Label) for the Ren et al 2008 study?
# this is converted into NA with the as.numeric() above

var <- c("Anesthesia", "Int.Label", "n", "Mean", "Median", "SD", "SEM",
         "Type.RIC", "Time.cond", "Time.RIC", "No.Limbs", "Limbs", "No.Session",
         "No.Cycles", "Occlusion", "Reperfusion")


# Identify studies with more than one comparison (one control group vs one intervention group)
# to enable correct pivoting of data into wide format for meta-analysis

check_studies <- df %>% 
  group_by(RefID) %>% 
  tally() %>% 
  filter(n > 2) %>% 
  select(RefID) %>% 
  left_join(df, by = "RefID")

# Pivot studies with several experimental conditions and one control group
# do we need adjustment of control group values?
# from protocol: For studies using a single control group and multiple experimental groups, 
# control group sample sizes were split to avoid double counting control animals.

studies_one_contr <- check_studies %>% 
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
  fill(n_Control, Mean_Control, Median_Control, SD_Control, SEM_Control, .direction = "downup") %>% 
  filter(!Condition2 %in% "Control") %>%
  ungroup()

# Pivot studies with several experiment-control comparisons

study_264 <- check_studies %>% 
  filter(RefID %in% "264") %>% 
  mutate(group = case_when(grepl(" A", Int.Label) ~ "A",
                           grepl("B", Int.Label) ~ "B")) %>%
  pivot_wider(names_from = Condition,
              values_from = c(n, Mean, Median, SD, SEM)) %>% 
  group_by(group) %>% 
  fill(n_Control, Mean_Control, Median_Control, SD_Control, SEM_Control, .direction = "downup") %>% 
  filter(!Condition2 %in% "Control") %>%
  ungroup() %>% 
  select(-group)

study_138 <- check_studies %>% 
  filter(RefID %in% "138") %>% 
  mutate(Int.Label = case_when(grepl("12h|12 h", Int.Label) ~ "12h",
                           grepl("24", Int.Label) ~ "24h",
                           grepl("120", Int.Label) ~ "120h")) %>%
  pivot_wider(names_from = Condition,
              values_from = c(n, Mean, Median, SD, SEM)) %>% 
  group_by(Int.Label) %>% 
  fill(n_Control, Mean_Control, Median_Control, SD_Control, SEM_Control, .direction = "downup") %>% 
  filter(!Condition2 %in% "Control") %>%
  ungroup()
  
# for studies 157, 317, 324: not clear which control is used for which experimental condition from the labels

# protocol: Where experimental arms outnumbered control group animals, the review team made decisions on the most relevant experimental arms, 
# in consultation with stroke experts. 

problem <- check_studies %>% 
  filter(RefID %in% c(157, 317, 324))

# bind all multiple comparison studies together
studies_multi_comp <- rbind(study_264, study_138)

df_wide <- df %>% 
  filter(!RefID %in% check_studies$RefID) %>% 
  group_by(RefID) %>%
  pivot_wider(names_from = Condition,
              values_from = all_of(var)) %>% 
  select(-c(Condition2)) %>% 
  fill(everything(), .direction = "downup") %>% 
  distinct() %>% 
  ungroup %>% 
  bind_rows(studies_one_contr) %>% 
  bind_rows(studies_multi_comp) %>% 
  # need this to only keep comparisons with quantitative data
  filter(!grepl("Infarct volume not provided|Infarct volume not reported|No infarct volume reported", Comments)) %>% 
  # harmonize error measure - will we convert SEM to SD?
  mutate(SD_RIC = case_when(is.na(SD_RIC) & !is.na(SEM_RIC) ~ SEM_RIC * sqrt(n_RIC),
                            .default = SD_RIC),
         SD_Control = case_when(is.na(SD_Control) & !is.na(SEM_Control) ~ SEM_Control * sqrt(n_Control),
                            .default = SD_Control))

df.wide_es <- escalc(measure = "SMD",   # Specifies the type of effect size (Standardized Mean Difference)
                     m1i = Mean_RIC,         # Mean of treatment group
                     sd1i = SD_RIC,       # SD of treatment group
                     n1i = n_RIC,         # Sample size of treatment group
                     m2i = Mean_Control,         # Mean of control group
                     sd2i = SD_Control,       # SD of control group
                     n2i = n_Control,         # Sample size of control group
                     data = df_wide)
