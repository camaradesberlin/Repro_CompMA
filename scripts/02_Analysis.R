# libraries

# source data

df_wide <- data.table::fread(here::here("data","tidy","df_wide.csv"), header = T)

# Main Analysis -----------------------------------------------------------

# random effects inverse variance meta-analyses and presented with 95% confidence intervals.
# # I2 for heterogeneity
df.wide_es <- escalc(measure = "SMD",   # Specifies the type of effect size (Standardized Mean Difference)
                     m1i = Mean_RIC,         # Mean of treatment group
                     sd1i = SD_RIC,       # SD of treatment group
                     n1i = n_RIC,         # Sample size of treatment group
                     m2i = Mean_Control,         # Mean of control group
                     sd2i = SD_Control,       # SD of control group
                     n2i = n_Control,         # Sample size of control group
                     data = df_wide)



# Sub-group analyses ------------------------------------------------------



