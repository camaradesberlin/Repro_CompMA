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
  select(-'Method of infarct measurement (latest method used)', 'Was infarct volume provided or calculated?') %>%
  rename(
    
  )