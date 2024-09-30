
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
                     data = df_wide) %>% 
  # remove study with Median (unused in MA)
  filter(!is.na(yi)) %>% 
  arrange(Author)

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


## Forest plot ----------------

pdf(file = here::here("figs","forestplot.pdf"), width = 9, height = 26)

meta::forest(ma_results_main, 
             sortvar = ma_results_main$data$Author,
             prediction = FALSE,
             print.I2 = TRUE,
             overall = TRUE,
             xlim = c(-10, 10),
             leftlabs = c("Author"))

dev.off()

pdf(file = here::here("figs","forestplot_byES.pdf"), width = 9, height = 26)

meta::forest(ma_results_main, 
             sortvar = TE,
             prediction = FALSE,
             print.I2 = TRUE,
             overall = TRUE,
             xlim = c(-10, 10),
             leftlabs = c("Author"))

dev.off()


## Funnel plot asymmetry ----------------

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
                   subgroup = No.Limbs)

# which limb
limb <- update(ma_results_main,
               subgroup = Limbs)

# type of RIC
type.RIC <- update(ma_results_main,
                   subgroup = Type.RIC)

# stroke model
stroke.model <- update(ma_results_main,
                       subgroup = Stroke.Type) 

# anaesthesia
anesthesia <- update(ma_results_main,
                      subgroup = Anesthesia_RIC)

# time.cond
time.cond <- update(ma_results_main,
                     subgroup = Time.cond)

# information on randomisation and blinding is not in the data

## Meta-regression analyses
# number of sessions
no.session <- metareg(ma_results_main, ~No.Session) 
#error message but don't understand why; potentially try metaregression

# number of cycles
no.cycles <- metareg(ma_results_main, ~No.Cycles)

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



# Small data --------------------------------------------------------------

main_i2 <- round(ma_results_main$I2*100, digits = 0)
main_i2_lower <- round(ma_results_main$lower.I2*100, digits = 0)
main_i2_upper <- round(ma_results_main$upper.I2*100, digits = 0)

main_smd <- round(ma_results_main$TE.random, digits = 2)
main_lower <- round(ma_results_main$lower.random, digits = 2)
main_upper <- round(ma_results_main$upper.random, digits = 2)

ripley_main_smd <- -2.19

es_change <-  abs((main_smd - ripley_main_smd) / ((main_smd + ripley_main_smd) / 2) * 100)

es_change <- round(es_change, digits = 1)

# smd_male <- sex$TE.random.w["Male"]
# smd_female <- sex$TE.random.w["Female"]


Rversion <- getRversion()
tidyverse_version <- utils::packageVersion("tidyverse")
metafor_version <- utils::packageVersion("metafor")
meta_version <- utils::packageVersion("meta")



# Save small data ---------------------------------------------------------

# saveRDS(list(es_change = es_change, main_i2 = main_i2, 
#              main_i2_lower = main_i2_lower, 
#              main_i2_upper = main_i2_upper,
#              main_smd = main_smd, main_lower = main_lower, 
#              main_upper = main_upper), "data/small_data.rds")

# save(ma_results_main, sex, species, no.limbs, no.cycles, time.occlusion,
#      time.occlusion.cond, time.reperfusion, 
#      time.reperfusion.cond, RIC.onset.cond, RIC.onset,
#      no.cycles, no.session, anesthesia, stroke.model, limb,
#      file = here::here("data","subgroup_models.Rdata"))

save(Rversion, tidyverse_version, metafor_version,
     meta_version, file = here::here("data","Renv.Rdata"))

save(es_change, main_i2,
     main_i2_upper, main_i2_lower,
     main_smd, main_lower, main_upper,
     file = here::here("data","small_data.Rdata"))

