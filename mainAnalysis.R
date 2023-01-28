################################################################################
#                             Main Analysis Script
################################################################################

# Overview:
  # - Script acts on a dataset comprised on 30 participants over 7 gait variables
  #  and several demographic variables
  # - Data are cleaned and evaluated using generalized linear mixed models
  # - results are saved in "results" folder

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

data <- read.csv("steps2stability_data.csv")
ntotal_step<- read.csv("ntotal_steps.csv")

# designate factor variables
data$ID <- as.factor(data$ID)
data$group <- as.factor(data$group)
data$Sex <- as.factor(data$Sex)
data$Klgrade <- as.factor(data$Klgrade)
str(data)

# Replace censored data with participant-outcome specific total number of steps 
#   recorded plus 5.
uncensored_outcomes = matrix(NA, 30, length(vars))
vars = c("peakacc_norm_HS", "peakacc_norm_TO","peakgyr_norm_HS","peakgyr_norm_TO","stance_time", "foot_strike_angle", "fpa_classif")
for(v in 1:length(vars)){
  var = vars[v] # select variable of interest
  
  outcome_ = data[, c("ID", "group", "Age", "Sex", "BMI", "Klgrade", var)]
  colnames(outcome_)[colnames(outcome_) == var] = "steps"
  outcome_$status = as.numeric(!is.na(outcome_$steps)) # note that 1 = event, 0 = censored
  for(i in 1:nrow(outcome_)){
    if (outcome_$status[i] ==0){
      # if censored -> find index of nsteps for that participant and add it in +5 steps
      # Note '+5' was used given that the underlying consistency analysis was done using blocks of 5 steps. 
      # That is, our resolution was 5 steps, and all the steps to consistency point estimates are multiples of 5
      outcome_$steps[i] = ntotal_step[i,var] + 5
    }
  }
  
  uncensored_outcomes[,v] = outcome_$steps
  
}
uncensored_data <-cbind(as.data.frame(uncensored_outcomes), 
                                     data$ID, data$group,  data$Sex, data$Klgrade)
colnames(uncensored_data) <- c(vars,"ID","group","Sex","Klgrade")

# Convert from wide to long form
data_long <-melt(uncensored_data)
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

# write dataframes to file
dir.create(file.path(getwd(), "results"), showWarnings = FALSE) # if "/results/" exists ignore warning
save(list = c("uncensored_data", "data_long"), 
     file = "./results/stability_dfs.RData")

################################# MODELLING ####################################
# POISSON distribution with log link function ----------------------------------
mod1 <- lme4::glmer(value ~ 1 + group + variable + (1|ID), 
                    family = poisson (link = "log"), 
                    data = data_long)
# Diagnostics
mod1.simResid <- DHARMa::simulateResiduals(fittedModel = mod1, plot=F, 
                                             n = 1000, seed = 123)
plotQQunif(mod1.simResid) # KS test significant - indicates deviation from expected distrib (poisson)
plotResiduals(mod1.simResid) # nonlinearities
DHARMa::testDispersion(mod1.simResid) # indicates overdispersion is potential issue 

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
# Generally a better fit than Poisson - using neg binom

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
       file = "./results/fixed_effects_tbl.html"
       )

# Pairwise Comparisons --------------------------
# Interaction Model
emmeans(mod1.nb.int, ~ 1, type="response")
estMargMeans_int<-emmeans(mod1.nb.int, ~ group*variable, type = "response")
estMargMeans_var<-emmeans(mod1.nb.int, ~ variable, type = "response")

contrastsResults_groupBYvar <- as.data.frame(tidy(pairs(estMargMeans_int, adjust = "bonferroni", type = "response")))
contrastsResults_var <- as.data.frame(tidy(pairs(estMargMeans_var, adjust = "bonferroni", type = "response")))
sjPlot::tab_df(contrastsResults_var[,c("contrast", "ratio","std.error","df", "statistic","adj.p.value")], show.rownames = TRUE, 
               digits = 3,
               file = "./results/variable_pairwise_comps.html"
)

# nicer effect names
modelPredNames <- c("Group: KOA", "Outcome: AC-TO" , "Outcome: GYR-HS", 
                    "Outcome: GYR-TO",  "Outcome: ST",        
                    "Outcome: FSA", "Outcome: FPA",
                    "Interaction: KOA*AC-TO",
                    "Interaction: KOA*GYR-TO",
                    "Interaction: KOA*GYR-TO",
                    "Interaction: KOA*ST",
                    "Interaction: KOA*FSA",
                    "Interaction: KOA*FPA")
(stabilityModel.table <- mod1.nb.int %>%
    sjPlot::tab_model(
      show.intercept = FALSE, show.ci = 0.95, show.p = TRUE, show.stat = TRUE,
      file = './results/stability_glmer_table.html',
      pred.labels = modelPredNames,
      dv.labels = c(""),
      string.est = "Estimate", string.ci = "95% CI", string.p = "p-value",string.stat = "z score",
      digits = 2,
      digits.p = 3,
      emph.p = TRUE,
      CSS = list(
        css.body = paste('font-size:', 12, ';', sep = ""),
        css.firsttablerow = 'border-bottom:2.25px solid black; border-top:1px solid black;',
        css.thead = 'border-top: 1px; text-align:center; font-style: normal; font-weight: bold;',
        css.depvarhead = 'text-align: center; border-bottom: 2px solid; font-style: normal; font-weight: bold;',
        css.depvarheadnodv = 'font-weight: bold;',
        css.body = '+opacity: 0;',
        css.randomparts = '+border-bottom: 1px solid'
      )
    )) 


# SUMMARY STATS ----------------------------------------------------------------
median_table<- tbl_summary(uncensored_data[,c(vars,"group")], by=group,
) %>% add_overall() %>%
  modify_caption("Median steps to consistency point.") %>%
  as_gt() 

gt::gtsave(median_table, 
           "./results/median_table.html")

# PERCENT OF SIMILIARITY WHEN STABILITY OCCURS----------------------------------

percofStab <- read.csv("percent_similar_data.csv")
percofStab$ID <-as.factor(percofStab$ID)
percofStab$group<-as.factor(percofStab$group)
str(percofStab)
percofStab_long = melt(percofStab)

# Summary stats 

percent_sim_table<- tbl_summary(percofStab[,c(vars,"group")], by=group, missing="no",
                                statistic = all_continuous() ~ "{mean} ({sd})"
) %>% add_overall() %>%
  modify_caption("Percent of similar data at the consistency point.") %>%
  as_gt() 

gt::gtsave(percent_sim_table, 
           "./results/percent_sim_table.html")

# SAVE WORKSPACE ---------------------------------------------------------------

save.image("analysis_results.RData")

# End