
# Define fill colors for contour
col.contour = c("gray75", "gray85", "gray95")

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

## calculate random effect inverse variance meta-analysis
ma_results_main <- metagen(
  TE = yi,
  seTE = vi,
  studlab = Author,
  data = df.wide_es,
  sm = "SMD",
  fixed = FALSE,
  random = TRUE,
  title = "Random effect inverse variance meta-analysis of the effects on infarct volume"
)

summary(ma_results_main)

## assessment of funnel plot asymmetry
# Generate funnel plot
meta::funnel(ma_results_main, 
             yaxis = "invvar",
             contour = c(0.9, 0.95, 0.99),
             col.contour = col.contour)

# Add a legend
legend(x = 1.6, y = 50, 
       legend = c("p < 0.1", "p < 0.05", "p < 0.01"),
       fill = col.contour)

# Add a title
title("Contour-Enhanced Funnel Plot (Remote Ischemic Conditioning)")

# Eggers regression
ma_result_egger <- eggers.test(ma_results_main)
## test suggest presence of funnel plot asymmetry

# Sub-group analyses ------------------------------------------------------
# Species
species <- update(ma_results_main,
                  subgroup = Species)

# Sex
sex <- update(ma_results_main,
              subgroup = Sex)

# Number of limbs
no.limbs <- update(ma_results_main,
                   subgroup = No.Limbs_RIC)

# which limb
limb <- update(ma_results_main,
               subgroup = Limbs_RIC)

# type of RIC
type.RIC <- update(ma_results_main,
                   subgroup = Type.RIC_RIC)

# stroke model
stroke.model <- update(ma_results_main,
                       subgroup = Stroke.Type) 

# anaesthesia
anaesthesia <- update(ma_results_main,
                      subgroup = Anesthesia_RIC)

# information on randomisation and blinding is not in the data

## Meta-regression analyses
# number of sessions
no.session <- metareg(ma_results_main, ~No.Session_RIC) 
#error message but don't understand why; potentially try metaregression

# number of cycles
no.cycles <- metareg(ma_results_main, ~No.Cycles_RIC)

# duration of RIC occlusion
time.occlusion <-  metareg(ma_results_main, ~Occlusion)
time.occlusion.cond <- df.wide_es %>%
  group_by(Time.cond) %>%
  group_map(~ rma(yi = yi,
                  sei = vi,
                  data = df.wide_es,
                  mods = ~Occlusion))

# duration of RIC reperfusion
time.reperfusion <- metareg(ma_results_main, ~Reperfusion)
time.reperfusion.cond <- df.wide_es %>%
  group_by(Time.cond) %>%
  group_map(~ rma(yi = yi,
                  sei = vi,
                  data = df.wide_es,
                  mods = ~Reperfusion))

# time of RIC onset
RIC.onset <- metareg(ma_results_main, ~Time.RIC)
RIC.onset.cond <- df.wide_es %>%
  group_by(Time.cond) %>%
  group_map(~ rma(yi = yi,
                  sei = vi,
                  data = df.wide_es,
                  mods = ~Time.RIC))

# latest infarct measurement timepoint
Time.infarct <- metareg(ma_results_main, ~Time.Inf)