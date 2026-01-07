# Script created 03/12/2025

# This script explores whether the effect of snail pair (reflecting the fact that snails were paired in the latter part of the experiment, and were therefore non-independent) needed to be included in the final analysis, as detailed in Appendix S1: Section S5 of the associated manuscript

### Shell height --------------------------------------------------------------

sizedata <- read_csv("processed_data/sizedata.csv")

# Ensure all variables are coded correctly (important if loaded as processed data)
{sizedata$Snail_ID=as.factor(sizedata$Snail_ID)  # Individual snail identifier
  sizedata$Treat=as.factor(sizedata$Treat) # Thermal regime (cold or warm) 
  sizedata$Pop=as.factor(sizedata$Pop) # Population (NH or NC)
  sizedata$Egg_tag=as.factor(sizedata$Egg_tag) # Maternal line
  sizedata$Sex=as.factor(sizedata$Sex) # Male or Female
  sizedata$Measure=as.factor(sizedata$Measure) # Measurement time point
  sizedata$Is_Replacement=as.factor(sizedata$Is_Replacement) # Whether the snail was a reload added late to the experiment to compensate for early mortality (see methods)
  sizedata$Condo=as.factor(sizedata$Condo) # Condo
  
}



# Compared models where Condo is included vs where it is excluded

model=glmmTMB(m ~ Treat*Pop*Measure*Sex+ (1|Snail_ID), REML = TRUE, data=sizedata, family=gaussian(link="identity"))

model2=glmmTMB(m ~ Treat*Pop*Measure*Sex+ (1|Condo) + (1|Snail_ID), REML = TRUE, data=sizedata, family=gaussian(link="identity"))

anova(model,model2) # adding condo does nothing good to AIC




# Create separate version of dataframe where measurement timepoints pre-paring are excluded

sizedatapostpairing <- subset(sizedata,  !(Measure %in% c("1","2","3","4","5","6","7", "8")))


model=glmmTMB(m ~ Treat*Pop*Measure*Sex+ (1|Snail_ID), REML = TRUE, data=sizedatapostpairing, family=gaussian(link="identity"))

model2=glmmTMB(m ~ Treat*Pop*Measure*Sex+ (1|Condo) + (1|Snail_ID), REML = TRUE, data=sizedatapostpairing, family=gaussian(link="identity"))

anova(model,model2) # Condo is rejected based on model comparison



### Growth rate--------------------------------------------------------------------------------

growthdata <- read_csv("processed_data/growthdata.csv")

# Ensure all variables are coded correctly (important if loaded as processed data)
{growthdata$Snail_ID=as.factor(growthdata$Snail_ID)  # Individual snail identifier
  growthdata$Treat=as.factor(growthdata$Treat) # Thermal regime (cold or warm) 
  growthdata$Pop=as.factor(growthdata$Pop) # Population (NH or NC)
  growthdata$Egg_tag=as.factor(growthdata$Egg_tag) # Maternal line
  growthdata$Sex=as.factor(growthdata$Sex) # Male or Female
  growthdata$Measure=as.factor(growthdata$Measure) # Measurement timepoint
  growthdata$Is_Replacement=as.factor(growthdata$Is_Replacement) # Whether the snail was a reload added late to the experiment to compesate for early mortality (see methods)
  growthdata$Condo=as.factor(growthdata$Condo) # Condo
  }



# Compared models where Condo is included vs where it is excluded

model=glmmTMB(pgrowth~ Treat*Pop*Measure*Sex+ (1|Snail_ID), REML = TRUE, data=growthdata, family=gaussian(link="identity"),dispformula = ~ Measure*Treat*Pop,control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)))
model2=glmmTMB(pgrowth~ Treat*Pop*Measure*Sex+ (1|Condo) + (1|Snail_ID), REML = TRUE, data=growthdata, family=gaussian(link="identity"),dispformula = ~ Measure*Treat*Pop,control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)))

anova(model,model2)



# Create separate version of dataframe where measurement timepoints pre-paring are excluded

growthdatapostpairing <- subset(growthdata,  !(Measure %in% c("1","2","3","4","5","6","7", "8")))


model=glmmTMB(pgrowth~ Treat*Pop*Measure*Sex + (1|Snail_ID), REML = TRUE, data=growthdatapostpairing, family=gaussian(link="identity"),dispformula = ~ Measure*Treat*Pop,control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)))

model2=glmmTMB(pgrowth~ Treat*Pop*Measure*Sex+ (1|Condo) + (1|Snail_ID), REML = TRUE, data=growthdatapostpairing, family=gaussian(link="identity"),dispformula = ~ Measure*Treat*Pop,control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)))

anova(model,model2) # Condo is rejected based on model comparison



### Weight data--------------------------------------------------

weightdata <- read_csv("processed_data/weightdata.csv")

# Ensure all variables are coded correctly (important if loaded as processed data)
{
  weightdata$Snail_ID=as.factor(weightdata$Snail_ID)  # Individual snail identifier
  weightdata$Treat=as.factor(weightdata$Treat) # Thermal regime (cold or warm)
  weightdata$Pop=as.factor(weightdata$Pop) # Population (NH or NC)
  weightdata$Egg_tag=as.factor(weightdata$Egg_tag) # Maternal line
  weightdata$Sex=as.factor(weightdata$Sex) # Male or Female
  weightdata$Is_Replacement=as.factor(weightdata$Is_Replacement) # Whether the snail was a reload added late to the experiment to compesate for early mortality (see methods)
  weightdata$Month=as.factor(weightdata$Month) # Measurement timepoint
  weightdata$time=as.POSIXct(weightdata$Date,format="%Y-%m-%d")
  weightdata$Condo=as.factor(weightdata$Condo)
}


# Compare models with and without condo

model=glmmTMB(Tissue_weight~ Treat*Pop*Month*Sex+ (1|Snail_ID), REML = TRUE, data=weightdata, family=Gamma(link="log"))
model2=glmmTMB(Tissue_weight~ Treat*Pop*Month*Sex+ (1|Condo) + (1|Snail_ID), REML = TRUE, data=weightdata, family=Gamma(link="log"))

anova(model,model2)


# Create separate version of dataframe where measurement timepoints pre-paring are excluded

weightdatapostpairing <- subset(weightdata,  !(Month  %in% c("2022-10-01","2023-01-01","2023-02-01","2023-03-01")))


model=glmmTMB(Tissue_weight~ Treat*Pop*Month*Sex+ (1|Snail_ID), REML = TRUE, data=weightdatapostpairing, family=Gamma(link="log"))
model2=glmmTMB(Tissue_weight~ Treat*Pop*Month*Sex+ (1|Condo) + (1|Snail_ID), REML = TRUE, data=weightdatapostpairing, family=Gamma(link="log"))

anova(model,model2)

