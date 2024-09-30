################################################################################
## SR-SAVI Reproducibility challenge
## Computational reproducibility of Ripley et al 2021
## Authors: Maria Economou, Torsten Rackoll
## Date: 30.09.2024
################################################################################

# Libraries ---------------------------------------------------------------

pacman::p_load(tidyverse, readxl, meta, metafor, knitr, here)

# remotes::install_github("MathiasHarrer/dmetar")
library(dmetar)

# generate bib file for R packages
knitr::write_bib(c("metafor","dmetar","meta","tidyverse"), file = here::here("docs/R_references.bib"))

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
  # code N/A as NA
  mutate(Sex = case_when(Sex %in% "male" ~ "Male",
                         Sex %in% "50% male" ~ "Both",
                         .default = Sex)) %>% 
  # harmonize sex variable
  mutate(Occlusion.time = na_if(Occlusion.time, "N/A")) %>% 
  # need this to only keep comparisons with quantitative data
  filter(!grepl("Infarct volume not provided|Infarct volume not reported|No infarct volume reported", Comments)) %>%
  mutate(Condition2 = Condition)

# column Median contain text (same information as Int.Label) for the Ren et al 2008 study?
# this is converted into NA with the as.numeric() above

var <- c("Anesthesia", "Int.Label", "n", "Mean", "Median", "SD", "SEM",
         "Type.RIC", "Time.cond", "Time.RIC", "No.Limbs", "Limbs", "No.Session",
         "No.Cycles", "Occlusion", "Reperfusion")


df <- df %>% 
  # manually exclude cohorts that the authors excluded (red highlighted in suppl material)
  mutate(keep = case_when(RefID %in% 288 & !Int.Label %in% c("10/10 3 cycles","15/15 3 cycles","MCAO") ~ FALSE,
                          RefID %in% 286 & Int.Label %in% c("RIC-3-15","RIC-6-15") ~ FALSE,
                          RefID %in% 327 & Int.Label %in% "RIPC 30 min" ~ FALSE,
                          RefID %in% 264 & Int.Label %in% c("RIC A1","RIC A4","RIC A5","RIC A6","RIC B4") ~ FALSE,
                          .default = TRUE)) %>% 
  filter(keep %in% TRUE) %>% 
  select(-keep)

# Identify studies with more than one comparison (one control group vs one intervention group)
# to enable correct pivoting of data into wide format for meta-analysis

check_studies <- df %>% 
  group_by(RefID) %>% 
  tally() %>% 
  filter(n > 2) %>% 
  select(RefID) %>% 
  left_join(df, by = "RefID")

# Identify studies where # experimental arms > N control
# studies_manually_adjusted <- check_studies %>%
#   group_by(RefID, Condition) %>%
#   tally() %>%
#   filter(!Condition %in% "Control") %>%
#   rename(n_exp = n) %>%
#   left_join((check_studies %>%
#                filter(Condition %in% "Control") %>%
#                select(RefID, n))) %>%
#   filter(n_exp > n) %>%
#   pull(RefID) %>%
#   unique()


# Identify studies with several experimental conditions and one control group
# and pivot into wide format

studies_one_contr <- check_studies %>% 
  # remove unused control group from RefID 317
  filter(!(RefID %in% 317 & Int.Label %in% "Control")) %>% 
  group_by(RefID, Condition) %>% 
  tally() %>% 
  pivot_wider(names_from = Condition, values_from = n) %>% 
  filter(Control == 1) %>% 
  select(RefID) %>% 
  left_join(df, by = "RefID") %>% 
  ungroup() %>% 
  group_by(RefID) %>%
  # pivot_wider(names_from = Condition,
  #             values_from = c(n, Mean, Median, SD, SEM)) %>% 
  pivot_wider(names_from = Condition,
              values_from = c(n, Anesthesia, Mean, SD, SEM, Median)) %>% 
  fill(n_Control, Mean_Control, Median_Control, SD_Control, SEM_Control, Anesthesia_Control, .direction = "downup") %>% 
  filter(!Condition2 %in% "Control") %>%
  ungroup()

studies_one_contr_sampleAdj <- studies_one_contr %>% 
  group_by(RefID) %>% 
  tally() %>% 
  left_join(studies_one_contr, by = "RefID") %>% 
  relocate(n_comp = n, .before = n_Control) %>% 
  mutate(n_Control = n_Control/n_comp) %>% 
  select(-n_comp)
  

# Identify studies with several experiment-control comparisons
# and pivot into wide format

to_pivot <- check_studies %>%
  filter(!RefID %in% studies_one_contr$RefID) %>%
  pull(RefID) %>%
  unique()

study_138 <- check_studies %>% 
  filter(RefID %in% "138") %>% 
  mutate(Int.Label = case_when(grepl("12h|12 h", Int.Label) ~ "12h",
                               grepl("24", Int.Label) ~ "24h",
                               grepl("120", Int.Label) ~ "120h")) %>%
  pivot_wider(names_from = Condition,
              values_from = c(n, Anesthesia, Mean, Median, SD, SEM)) %>% 
  group_by(Int.Label) %>% 
  fill(n_Control, Mean_Control, Median_Control, SD_Control, SEM_Control, Anesthesia_Control, .direction = "downup") %>% 
  filter(!Condition2 %in% "Control") %>%
  ungroup()

study_157 <- check_studies %>% 
  filter(RefID %in% "157") %>% 
  filter(!grepl("MCAO \\(", Int.Label)) %>% 
  mutate(group = Anesthesia) %>% 
  pivot_wider(names_from = Condition,
              values_from = c(n, Anesthesia, Mean, Median, SD, SEM)) %>% 
  group_by(group) %>% 
  fill(n_Control, n_RIC, Mean_Control, Mean_RIC, Median_Control, Median_RIC, 
       SD_Control, SD_RIC, SEM_Control, SEM_RIC, Anesthesia_Control, Anesthesia_RIC, .direction = "downup") %>% 
  filter(!Condition2 %in% "Control") %>%
  ungroup() %>% 
  select(-group)
  
study_324 <- check_studies %>% 
  filter(RefID %in% "324") %>% 
  mutate(group = case_when(grepl("Rapid", Int.Label) ~ "Rapid",
                           grepl("12h", Int.Label) ~ "12h",
                           grepl("2d", Int.Label) ~ "2d")) %>%
  filter(!is.na(group)) %>% 
  pivot_wider(names_from = Condition,
              values_from = c(n, Anesthesia, Mean, Median, SD, SEM)) %>% 
  group_by(group) %>% 
  fill(n_Control, n_RIC, Mean_Control, Mean_RIC, Median_Control, Median_RIC, 
       SD_Control, SD_RIC, SEM_Control, SEM_RIC, Anesthesia_Control, Anesthesia_RIC, .direction = "downup") %>% 
  filter(!Condition2 %in% "Control") %>%
  ungroup()

study_324_sampleAdj <- study_324 %>% 
  group_by(RefID, group) %>% 
  tally() %>% 
  left_join(study_324, by = c("RefID", "group")) %>% 
  relocate(n_comp = n, .before = n_Control) %>% 
  mutate(n_Control = case_when(n_comp > 1 ~ n_Control/n_comp,
                               .default = n_Control)) %>% 
  select(-c(n_comp, group))


study_264 <- df %>%
  filter(RefID %in% "264") %>%
  mutate(group = case_when(grepl(" A", Int.Label) ~ "A",
                           grepl("B", Int.Label) ~ "B")) %>%
  pivot_wider(names_from = Condition,
              values_from = c(n, Anesthesia, Mean, Median, SD, SEM)) %>%
  group_by(group) %>%
  fill(n_Control, Mean_Control, Median_Control, SD_Control, SEM_Control, Anesthesia_Control, .direction = "downup") %>%
  filter(!Condition2 %in% "Control") %>%
  ungroup()

study_264_sampleAdj <- study_264 %>% 
  group_by(RefID, group) %>% 
  tally() %>% 
  left_join(study_264, by = c("RefID", "group")) %>% 
  relocate(n_comp = n, .before = n_Control) %>% 
  mutate(n_Control = case_when(n_comp > 1 ~ n_Control/n_comp,
                               .default = n_Control)) %>% 
  select(-c(n_comp, group))

# Merge fixed studies
check_studies_fixed <- rbind(studies_one_contr_sampleAdj, 
                             study_138, study_157, study_264_sampleAdj, 
                             study_324_sampleAdj)


df_wide <- df %>% 
  filter(!RefID %in% check_studies$RefID) %>% 
  # add studies with sampled comparisons
  group_by(RefID) %>%
  # pivot_wider(names_from = Condition,
  #             values_from = all_of(var)) %>% 
  pivot_wider(names_from = Condition,
              values_from = c(n, Anesthesia, Mean, Median, SD, SEM)) %>% 
  select(-c(Condition2, Int.Label)) %>% 
  fill(everything(), .direction = "downup") %>% 
  distinct() %>% 
  ungroup() %>% 
  # rename(Type.RIC = Type.RIC_RIC)
  # add studies with different combinations of experimental arms + control groups
  bind_rows(check_studies_fixed) %>% 
  # harmonize error measure - will we convert SEM to SD?
  mutate(SD_RIC = case_when(is.na(SD_RIC) & !is.na(SEM_RIC) ~ SEM_RIC * sqrt(n_RIC),
                            .default = SD_RIC),
         SD_Control = case_when(is.na(SD_Control) & !is.na(SEM_Control) ~ SEM_Control * sqrt(n_Control),
                            .default = SD_Control)) %>% 
  arrange(RefID)

# vars to check if complete
# type RIC

data.table::fwrite(df_wide, here::here("data","tidy","df_wide.csv"), row.names = F, col.names = T)

