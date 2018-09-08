

# PublicationBias ---------------------------------------------------------
funnel(mod, xlab = "Residual Value", back = "white") #funnel plot for qualitative evaluation publication bias
#dev.copy(tiff,'MAB_figures/MAB_male_funnel.tiff')##save for males
#dev.copy(tiff,'MAB_figures/MAB_female_funnel.tiff')##save for females
#dev.off()

#random effects model for evaluation of publication bias
modRanEff <- rma.uni(yi, vi,
                     #random = list(~1 | each, ~1 | exp),
                     mods = ~domain:hit2Grouped - 1,
                     method = "REML",
                     digits = 3,
                     data = dat)


fsn

#trim and fill
nomod <- rma.uni(yi, vi,
                 #random = list(~1 | each, ~1 | exp),
                 #mods = ~domain:hit2Grouped - 1,
                 method = "REML",
                 digits = 3,
                 data = dat)
trim <- trimfill(nomod) #0 papers >> what does it mean?
trim












# MetaForest: dataset preparation --------------------------------------------------------------
data$yi <- as.numeric(data$yi)
#data$exp <- as.factor(data$exp) #if I do not do this, the model does not run

#createDataset
data %>% 
  select(yi, year,
         species, speciesStrainGrouped, origin, sex, ageWeek, #animals
         
         model, mTimeLength, #models
         mCageGrouped, mLDGrouped, mControlGrouped,
         
         hit2Grouped, 
         
         domain, 
         testAuthorGrouped, 
         #varAuthorGrouped, >>removed because some levels have only a few values
         waterT, testLDGrouped, 
         freezingType, retentionGrouped, 
         
         seqGeneration, baseline, allocation, 
         #housing,  #remove because only 1 level
         blindExp, control, outAss,
         bias, blindRand, vi, exp) %>% 
  filter(testAuthorGrouped != "stepDownAvoidance") %>%
  droplevels() -> datForest


# MetaForest: Tuning ------------------------------------------------------------------
##Attention: Tuning assumes that convergence is reached

set.seed(3238) #set seeds to a random number


# Set up 10-fold grouped CV
fit_control <- trainControl(method = "cv", index = groupKFold(datForest$exp, k = 10))

# Set up a custom tuning grid for the three tuning parameters of MetaForest
rf_grid <- expand.grid(whichweights = c("random", "fixed", "unif"),
                       mtry = c(2, 4, 6),
                       min.node.size = c(2, 4, 6))

# Train the model
cv.mf.cluster <- train(y = datForest$yi, x = datForest[,-1],#from x remove yi and each
                       study = "exp", method = ModelInfo_mf(),
                       trControl = fit_control,
                       tuneGrid = rf_grid)
cv.mf.cluster

#cross validated R2 - .12, sd= .09
cv.mf.cluster$results[which(cv.mf.cluster$results$whichweights == "unif" & 
                              cv.mf.cluster$results$mtry == 6 & 
                              cv.mf.cluster$results$min.node.size == 2),] #details with the "best" model


cv.mf.cluster$finalModel #R2oob = .09

#convergence plot
plot(cv.mf.cluster$finalModel) 
#dev.copy(tiff,'MAB_figures/MAB_metaforest_convergencePlot.tiff')
#dev.off()

VarImpPlot(cv.mf.cluster$finalModel)
#dev.copy(tiff,'MAB_figures/MAB_metaforest_varImportance.tiff')
#dev.off()
##remove variables with negative variable importance?

# Metaforest plots --------------------------------------------------
variable <- c("year","species","speciesStrainGrouped","origin","sex","ageWeek",
              "model","mTimeLength","mCageGrouped","mLDGrouped","mControlGrouped",
              "hit2Grouped","domain","testAuthorGrouped","waterT","testLDGrouped",
              "freezingType","retentionGrouped","seqGeneration","baseline","allocation",
              "blindExp","control","outAss","bias","blindRand") 

a <- 1

for (i in c(1:length(variable))) {
  
  thisMod <- variable[i]
  WeightedScatter(datForest, yi = "yi", vi = "vi", vars = thisMod, summarize = TRUE)
  ggsave(paste0("MAB_figures/MAB_metaforest_WeightedScatter_", thisMod, ".tiff"))
  
  PartialDependence(cv.mf.cluster$finalModel, vars = thisMod) 
  ggsave(paste0("MAB_figures/MAB_metaforest_partDep_", thisMod, ".tiff"))
  
  a <- 1
  
}

