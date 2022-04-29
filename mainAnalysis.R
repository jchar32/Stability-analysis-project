################################################################################
#                             Main Analysis Script
################################################################################

# Overview:
  # - Script acts on a dataset comprised on 30 participants over 7 gait variables
  #  and several demographic variables
  # - Data are cleaned and evaluated using generalized linear mixed models
  # - results are saved in "data/" folder

# Built with: R version 4.1.2

################################################################################

# Initializations and constants ------------------------------------------------
packages <- c("readr", "tidyverse", "reshape2", "DHARMa", "emmeans", "glmmTMB", 
              "broom", "sjPlot", "car", "rstudioapi") 

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))


# set contrasts
options(contrasts = rep ("contr.treatment", 2))

# set working directory to root folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd() #confirm directory is correct (top level of downloaded folder where script is located)

# ------------------------------------------------------------------------------

data <- read_csv("./data/steps2stability_data.csv", col_names = TRUE)

# designate factor variables
data$ID <-as.factor(data$ID)
data$group<-as.factor(data$group)
data$Sex <-as.factor(data$Sex)
data$Klgrade <-as.factor(data$Klgrade)
str(data)

# Convert from wide to long form
data_long <-melt(data)
data_long$variable <- as.factor(data_long$variable)

data_long$variable <- factor(
  data_long$variable,
  level = c(
    "peakacc_norm_HS",
    "peakacc_norm_TO" ,
    "peakgyr_norm_HS" ,
    "peakgyr_norm_TO",
    "stance_time" ,
    "foot_strike_angle" ,
    "fpa_classif"
  )
)

save(list = c("data", "data_long"), 
     file = "./data/stability_dfs.RData")

################################# MODELLING ####################################
# POISSON distribution with log link function ----------------------------------
mod1 <- lme4::glmer(value ~ 1 + group + variable + (1|ID), 
                    family = poisson (link = "log"), 
                    data = data_long)
# Diagnostics
mod1.simResid <- DHARMa::simulateResiduals(fittedModel = mod1, plot=F, 
                                             n = 1000, seed = 123)
plotQQunif(mod1.simResid) # KS test significant - indicates deviation from expected distrib (poisson)
plotResiduals(mod1.simResid) # points are well scattered but maybe some nonlinear trends (minor issue)
DHARMa::testDispersion(mod1.simResid) # indicates overdispersion is potential issue 
# dispersion = 2.497, p-value = 0.02


# NEGATIVE BINOMIAL distributions ----------------------------------------------
# Negative binomial -------------------
mod1.nb <- glmmTMB::glmmTMB(value ~ 1 + group + variable + (1|ID), 
                             family = "nbinom2", 
                             data = data_long)
# Diagnostics
mod1.nb.simResid <- DHARMa::simulateResiduals(fittedModel = mod1.nb, plot=F, 
                                                n = 1000, seed = 123)

plotQQunif(mod1.nb.simResid) # looks good
plotResiduals(mod1.nb.simResid) # Looks better than poison - nice scattering of points
DHARMa::testDispersion(mod1.nb) # looks good
# dispersion = 1.2299, p-value = 0.272

# Negative binomial with interaction effect ------------
mod1.nb.int <- glmmTMB::glmmTMB(value ~ 1 + group*variable + (1|ID), 
                                family = "nbinom2", 
                                data = data_long,
)

# Diagnostics
mod1.nb.int.simResid <- DHARMa::simulateResiduals(fittedModel = mod1.nb.int, plot=F, 
                                                    n = 1000, seed = 123)

plotQQunif(mod1.nb.int.simResid) 
plotResiduals(mod1.nb.int.simResid) 
DHARMa::testDispersion(mod1.nb.int) 
# dispersion = 1.1655, p-value = 0.456

# test significance of adding interaction
anova(mod1.nb, mod1.nb.int) # Interaction model is not significantly better - but it is analytically important so keeping.

# Recalculate with restricted maximum likelihood
mod1.nb.int <- glmmTMB::glmmTMB(value ~ 1 + group*variable + (1|ID), 
                                family = "nbinom2", 
                                data = data_long,
                                REML = TRUE
)

# Fixed main effects significance testing --------------------
# Wald Chi2 test
mod1.nb.int.FE <- car::Anova(mod1.nb.int, type = "III") # type III to match wald z statistic in model parameter tests
tab_df(mod1.nb.int.FE, show.rownames = TRUE, 
       digits = 3,
       file = "./data/fixed_effects_tbl.html"
       )

# Pairwise Comparisons --------------------------
# Interaction Model
emmeans(mod1.nb.int, ~ 1, type="response")
estMargMeans_int<-emmeans(mod1.nb.int, ~ group*variable, type = "response")
estMargMeans_var<-emmeans(mod1.nb.int, ~ variable, type = "response")

contrastsResults_groupBYvar <- as.data.frame(tidy(pairs(estMargMeans_int, adjust = "bonferroni", type = "response")))
contrastsResults_var <- as.data.frame(tidy(pairs(estMargMeans_var, adjust = "bonferroni", type = "response")))


# SUMMARY STATS ----------------------------------------------------------------

# function to calculate summary stats with various data splitting --------------
sumStat <- function(data_in){

  data.summarised.all <- data_in %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      lquart = quantile(value, type = 5, na.rm = TRUE)[2],
      uquart = quantile(value, type = 5, na.rm = TRUE)[4],
      min = min(value, na.rm = TRUE),
      max = max(value, na.rm = TRUE)
    )
  
  data.summarised.var <- data_in %>%
    group_by(variable) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      lquart = quantile(value, type = 5, na.rm = TRUE)[2],
      uquart = quantile(value, type = 5, na.rm = TRUE)[4],
      min = min(value, na.rm = TRUE),
      max = max(value, na.rm = TRUE)
    )
  
  data.summarised.grp <- data_in %>%
    group_by(group) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      lquart = quantile(value, type = 5, na.rm = TRUE)[2],
      uquart = quantile(value, type = 5, na.rm = TRUE)[4],
      min = min(value, na.rm = TRUE),
      max = max(value, na.rm = TRUE)
    )
  
  data.summarised.var.grp <- data_in %>%
    group_by(variable, group) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      lquart = quantile(value, type = 5, na.rm = TRUE)[2],
      uquart = quantile(value, type = 5, na.rm = TRUE)[4],
      min = min(value, na.rm = TRUE),
      max = max(value, na.rm = TRUE)
    )
 
 data_out <- list(data.summarised.all, data.summarised.var, 
                  data.summarised.grp, data.summarised.var.grp)
  return(data_out)
}
# ------------------------------------------------------------------------------

summarized_Data_List<- sumStat(data_long)

raw.sumStats <- summarized_Data_List[[1]]
raw.sumStats.var <- summarized_Data_List[[2]]
raw.sumStats.grp<-  summarized_Data_List[[3]]
raw.sumStats.var.grp <- summarized_Data_List[[4]]


tab_df(raw.sumStats.var.grp,
       digits = 0,
       file = "./data/stepcount_varBygrp.html")
tab_df(raw.sumStats.var,
       digits = 0,
       file = "./data/stepcount_var.html")
tab_df(raw.sumStats.grp,
       digits = 0,
       file = "./data/stepcount_grp.html")
tab_df(raw.sumStats,
       digits = 0,
       file = "./data/stepcount_ALL.html")


# PERCENT OF SIMILIARITY WHEN STABILITY OCCURS----------------------------------

percofStab <- read_csv("./data/percent_similar_data.csv", 
                 col_names = TRUE)
percofStab$ID <-as.factor(percofStab$ID)
percofStab$group<-as.factor(percofStab$group)
str(percofStab)
percofStab_long = melt(percofStab)

# Summary stats ---------------------------------------------------------------

summarized_Data_List.perc<- sumStat(percofStab_long)

raw.sumStats.perc <- summarized_Data_List[[1]]
raw.sumStats.var.perc <- summarized_Data_List[[2]]
raw.sumStats.grp.perc<-  summarized_Data_List[[3]]
raw.sumStats.var.grp.perc <- summarized_Data_List[[4]]

tab_df(raw.sumStats.var.grp.perc,
       digits = 3,
       file = "./data/StabilityStepCountsTable_varBygrp_percsim.html")
tab_df(raw.sumStats.var.perc,
       digits = 3,
       file = "./data/StabilityStepCountsTable_var_percsim.html")
tab_df(raw.sumStats.grp.perc,
       digits = 3,
       file = "./data/StabilityStepCountsTable_grp_percsim.html")
tab_df(raw.sumStats.perc,
       digits = 3,
       file = "./data/StabilityStepCountsTable_all_percsim.html")


# SAVE WORKSPACE ---------------------------------------------------------------

save.image("./data/analysis_results.RData")

# End