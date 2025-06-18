

# Load required libraries------

{
library(plyr)
library(nlme)
library(lme4)
library(tidyr)
library(dplyr)
library(car)
library(readr)
library(lubridate)
library(emmeans)
library(multcomp)
library(multcompView)
library(glmmTMB)
library(data.table)
library(tidyverse)
library(sjPlot)
library(broom)
  
 
# set appropriate contrasts for type 3 SS ANOVA
options(contrasts = c("contr.sum", "contr.poly")) 
}


# Statistical outputs for all tests run in study, including summary statistics ------------------------------

# Data is passed from script 01. Data setup (if you run this, all data will be loaded into R already)

# alternatively, data can be loaded directly from the "processed data" folder using hashed read_csv commands


## Growth statistics ------------------------------------------------------------

### Shell height --------------------------------------------------------------

sizedata <- read_csv("processed_data/sizedata.csv")

# Ensure all variables are coded correctly (important if loaded as processed data)
{sizedata$Snail_ID=as.factor(sizedata$Snail_ID)  # Individual snail identifier
sizedata$Treat=as.factor(sizedata$Treat) # Thermal regime (cold or warm) 
sizedata$Pop=as.factor(sizedata$Pop) # Population (NH or NC)
sizedata$Egg_tag=as.factor(sizedata$Egg_tag) # Maternal line
sizedata$Sex=as.factor(sizedata$Sex) # Male or Female
sizedata$Measure=as.factor(sizedata$Measure) # Measurement timepoint
sizedata$Is_Replacement=as.factor(sizedata$Is_Replacement) # Whether the snail was a reload added late to the experiment to compesate for early mortality (see methods)
}

# Response variable "m" is height of the snail in mm

# Create glmmTMB model

model=glmmTMB(m ~ Treat*Pop*Measure+ (1|Snail_ID)+ (1|Sex)+ (1|Egg_tag), REML = TRUE, data=sizedata, family=gaussian(link="log"))

summary (model)

# Summary table
tab_df(tidy(Anova(model,type=3)), title = "Shell height", show.rownames = FALSE, show.type=FALSE,digits=2)
r2(model)

# Checks to see if other random effects structures were preferred (used in supplemental material)
{
model=glmmTMB(m ~ Treat*Pop*Measure+ (1|Snail_ID)+ (1|Sex)+ (1|Is_Replacement)+ (1|Egg_tag), REML = TRUE, data=sizedata, family=gaussian(link="log"))
model2=glmmTMB(m ~ Treat*Pop*Measure+ (1|Snail_ID)+ (1|Sex)+ (1|Is_Replacement), REML = TRUE, data=sizedata, family=gaussian(link="log"))
model3=glmmTMB(m ~ Treat*Pop*Measure+ (1|Snail_ID)+ (1|Sex)+ (1|Egg_tag), REML = TRUE, data=sizedata, family=gaussian(link="log"))
model4=glmmTMB(m ~ Treat*Pop*Measure+ (1|Snail_ID)+ (1|Is_Replacement)+ (1|Egg_tag), REML = TRUE, data=sizedata, family=gaussian(link="log"))
model5=glmmTMB(m ~ Treat*Pop*Measure+ (1|Snail_ID)+ (1|Sex), REML = TRUE, data=sizedata, family=gaussian(link="log"))
model6=glmmTMB(m ~ Treat*Pop*Measure+ (1|Snail_ID)+ (1|Is_Replacement), REML = TRUE, data=sizedata, family=gaussian(link="log"))
model7=glmmTMB(m ~ Treat*Pop*Measure+ (1|Snail_ID)+ (1|Egg_tag), REML = TRUE, data=sizedata, family=gaussian(link="log"))
model8=glmmTMB(m ~ Treat*Pop*Measure+ (1|Snail_ID), REML = TRUE, data=sizedata, family=gaussian(link="log"))
}

{anova_df=as.data.frame(anova(model,model2,model3,model4,model5,model6,model7,model8)) 
  anova_df$delta_aic = anova_df$AIC - min(anova_df$AI,na.rm=TRUE)
  anova_df=anova_df[,c(1,2,4,9)]
  anova_df <- anova_df[order(anova_df$delta_aic), ]
  tab_df(anova_df,show.rownames=TRUE)}




### Growth rate--------------------------------------------------------------------------------

#growthdata <- read_csv("processed_data/growthdata.csv")

# Ensure all variables are coded correctly (important if loaded as processed data)
{growthdata$Snail_ID=as.factor(growthdata$Snail_ID)  # Individual snail identifier
  growthdata$Treat=as.factor(growthdata$Treat) # Thermal regime (cold or warm) 
  growthdata$Pop=as.factor(growthdata$Pop) # Population (NH or NC)
  growthdata$Egg_tag=as.factor(growthdata$Egg_tag) # Maternal line
  growthdata$Sex=as.factor(growthdata$Sex) # Male or Female
  growthdata$Measure=as.factor(growthdata$Measure) # Measurement timepoint
  growthdata$Is_Replacement=as.factor(growthdata$Is_Replacement) # Whether the snail was a reload added late to the experiment to compesate for early mortality (see methods)

  }

# Response variable "pgrowth" refers to specific growth rate

hist(growthdata$pgrowth, breaks=30)

# Note there are a few values where growth is slightly negative - 
# potentially because of small measurement errors during periods where no growth was occurring/ handling damage.
# Therefore it seems reasonable to convert any negative growth values to zeros (because growth of shell length cant be negative)


growthdata$pgrowth2=growthdata$pgrowth
growthdata$pgrowth2 <- ifelse(growthdata$pgrowth2 < 0, 0, growthdata$pgrowth2)

hist(growthdata$pgrowth2, breaks=30)

# Create glmmTMB model
model=glmmTMB(pgrowth2 ~ Treat*Pop*Measure+ (1|Snail_ID)+ (1|Sex)+ (1|Is_Replacement), REML = TRUE, data=growthdata, family=tweedie(link="log"),dispformula=~Measure*Treat)
summary (model)

# Summary table
tab_df(tidy(Anova(model,type=3)), title = "Growth rate", show.rownames = FALSE, show.type=FALSE,digits=2)
r2(model)

# Checks to see if other random effects structures were preferred (used in supplemental material)

{model=glmmTMB(pgrowth2 ~ Treat*Pop*Measure+ (1|Snail_ID)+ (1|Sex)+ (1|Egg_tag) + (1|Is_Replacement), REML = TRUE, data=growthdata, , na.action = na.fail, family=tweedie(link="log"),dispformula=~Measure*Treat)
model2=glmmTMB(pgrowth2 ~ Treat*Pop*Measure+ (1|Snail_ID)+ (1|Sex)+ (1|Is_Replacement), REML = TRUE, data=growthdata, family=tweedie(link="log"),dispformula=~Measure*Treat)
model3=glmmTMB(pgrowth2 ~ Treat*Pop*Measure+ (1|Snail_ID)+ (1|Sex)+ (1|Egg_tag), REML = TRUE, data=growthdata, family=tweedie(link="log"),dispformula=~Measure*Treat)
model4=glmmTMB(pgrowth2 ~ Treat*Pop*Measure+ (1|Snail_ID)+ (1|Is_Replacement)+ (1|Egg_tag), REML = TRUE, data=growthdata, family=tweedie(link="log"),dispformula=~Measure*Treat)
model5=glmmTMB(pgrowth2 ~ Treat*Pop*Measure+ (1|Snail_ID)+ (1|Sex), REML = TRUE, data=growthdata, family=tweedie(link="log"),dispformula=~Measure*Treat)
model6=glmmTMB(pgrowth2 ~ Treat*Pop*Measure+ (1|Snail_ID)+ (1|Egg_tag), REML = TRUE, data=growthdata, family=tweedie(link="log"),dispformula=~Measure*Treat)
model7=glmmTMB(pgrowth2 ~ Treat*Pop*Measure+ (1|Snail_ID)+ (1|Is_Replacement), REML = TRUE, data=growthdata, family=tweedie(link="log"),dispformula=~Measure*Treat)
model8=glmmTMB(pgrowth2 ~ Treat*Pop*Measure+ (1|Snail_ID), REML = TRUE, data=growthdata, family=tweedie(link="log"),dispformula=~Measure*Treat)
}

{anova_df=as.data.frame(anova(model,model2,model3,model4,model5,model6,model7,model8)) 
  anova_df$delta_aic = anova_df$AIC - min(anova_df$AI,na.rm=TRUE)
  anova_df=anova_df[,c(1,2,4,9)]
  anova_df <- anova_df[order(anova_df$delta_aic), ]
  tab_df(anova_df,show.rownames=TRUE)}




### Weight data--------------------------------------------------

#weightdata <- read_csv("processed_data/weightdata.csv")

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
}


# Create glmmTMB model
model=glmmTMB(Tissue_weight~ Treat*Pop*Month+ (1|Snail_ID)+(1|Sex), REML = TRUE, data=weightdata, family=Gamma(link="log"))
summary (model)
Anova(model,type=3)

# Summary table
tab_df(tidy(Anova(model,type=3)), title = "Tissue weight", show.rownames = FALSE, show.type=FALSE,digits=2)

# Checks to see if other random effects structures were preferred (used in supplemental material)
{
model=glmmTMB(Tissue_weight~ Treat*Pop*Month+ (1|Snail_ID)+(1|Sex) + (1|Egg_tag)+(1|Is_Replacement), REML = TRUE, data=weightdata, family=Gamma(link="log"),control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)))
model2=glmmTMB(Tissue_weight~ Treat*Pop*Month+ (1|Snail_ID)+(1|Sex) + (1|Egg_tag), REML = TRUE, data=weightdata, family=Gamma(link="log"),control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)))
model3=glmmTMB(Tissue_weight~ Treat*Pop*Month+ (1|Snail_ID)+(1|Sex) +(1|Is_Replacement), REML = TRUE, data=weightdata, family=Gamma(link="log"),control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)))
model4=glmmTMB(Tissue_weight~ Treat*Pop*Month+ (1|Snail_ID)+(1|Egg_tag)+(1|Is_Replacement), REML = TRUE, data=weightdata, family=Gamma(link="log"),control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)))
model5=glmmTMB(Tissue_weight~ Treat*Pop*Month+ (1|Snail_ID)+(1|Is_Replacement), REML = TRUE, data=weightdata, family=Gamma(link="log"),control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)))
model6=glmmTMB(Tissue_weight~ Treat*Pop*Month+ (1|Snail_ID)+(1|Sex), REML = TRUE, data=weightdata, family=Gamma(link="log"),control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)))
model7=glmmTMB(Tissue_weight~ Treat*Pop*Month+ (1|Snail_ID)+(1|Egg_tag), REML = TRUE, data=weightdata, family=Gamma(link="log"),control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)))
model8=glmmTMB(Tissue_weight~ Treat*Pop*Month+ (1|Snail_ID), REML = TRUE, data=weightdata, family=Gamma(link="log"),control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)))
}

{anova_df=as.data.frame(anova(model,model2,model3,model4,model5,model6,model7,model8)) 
  anova_df$delta_aic = anova_df$AIC - min(anova_df$AI)
anova_df=anova_df[,c(1,2,4,9)]
anova_df <- anova_df[order(anova_df$delta_aic), ]
tab_df(anova_df,show.rownames=TRUE)}


## Feeding rate statistics (consumption data) --------------------------------------------

#consdata <- read_csv("processed_data/consdata.csv")

# Ensure all variables are coded correctly (important if loaded as processed data)
{consdata$Measure=as.factor(consdata$Measure) # Measurement timepoint
  consdata$Snail_ID=as.factor(consdata$Snail_ID)  # Individual snail identifier
  consdata$Treat=as.factor(consdata$Treat)  # Thermal regime (cold or warm)
  consdata$Pop=as.factor(consdata$Pop) # Population (NH or NC)
  consdata$Egg_tag=as.factor(consdata$Egg_tag) # Maternal line
  consdata$Sex=as.factor(consdata$Sex) # Male or Female
  consdata$Is_Replacement=as.factor(consdata$Is_Replacement) # Whether the snail was a reload added late to the experiment to compesate for early mortality (see methods)
  
}

# Response variable "c" refers to consumption in oysters / snail


boxplot(c~Pop*Measure,data=subset(consdata,consdata$Treat=="NH"))

boxplot(c~Pop*Measure,data=subset(consdata,consdata$Treat=="NC"))

# Zero consumption values occurs at several timepoints in cold (NH) regime,
# but never in warm regime 
# To avoid issues with total seperation of data at these timepoints, 
# we chose to analyse data separately for the two regimes,
# removing several timepoints for the cold regime only

### Analysis for Warm Regime ("NC") ----------------------------------------------

# Create glmmTMB model
modelNC=glmmTMB(c ~ Pop*Measure+ (1|Snail_ID), REML = FALSE, data=subset(consdata,consdata$Treat=="NC"),dispformula=~Measure, family=gaussian(link="identity"))
summary (modelNC)
Anova(modelNC,type=3)

# Summary table
tab_df(tidy(Anova(modelNC,type=3)), title = "Consumption (Warm regime)", show.rownames = FALSE, show.type=FALSE,digits=2)


# Checks to see if other random effects structures were preferred (used in supplemental material)

{
modelNC=glmmTMB(c ~ Pop*Measure+ (1|Snail_ID)+(1|Is_Replacement)+(1|Egg_tag), REML = FALSE, data=subset(consdata,consdata$Treat=="NC"),dispformula=~Measure, family=gaussian(link="identity"))
modelNC2=glmmTMB(c ~ Pop*Measure+ (1|Snail_ID)+(1|Egg_tag), REML = FALSE, data=subset(consdata,consdata$Treat=="NC"),dispformula=~Measure, family=gaussian(link="identity"))
modelNC3=glmmTMB(c ~ Pop*Measure+ (1|Snail_ID)+(1|Is_Replacement), REML = FALSE, data=subset(consdata,consdata$Treat=="NC"),dispformula=~Measure, family=gaussian(link="identity"))
modelNC4=glmmTMB(c ~ Pop*Measure+ (1|Snail_ID), REML = FALSE, data=subset(consdata,consdata$Treat=="NC"),dispformula=~Measure, family=gaussian(link="identity"))
}


{anova_df=as.data.frame(anova(modelNC,modelNC2,modelNC3,modelNC4)) 
  anova_df$delta_aic = anova_df$AIC - min(anova_df$AI)
  anova_df=anova_df[,c(1,2,4,9)]
  anova_df <- anova_df[order(anova_df$delta_aic), ]
  tab_df(anova_df,show.rownames=TRUE)}


### Analysis for Cold Regime ("NH") ----------------------------------------------

# Subset data, creating a version where timepoints where no feeding was observed in the cold regime are removed
consdataNH <- subset(consdata,  !(Measure %in% c("3","4","5", "6","7", "15")))

modelNH=glmmTMB(c ~ Pop*Measure+ (1|Snail_ID), REML = FALSE, data=subset(consdataNH,consdataNH$Treat=="NH"),dispformula=~Measure, family=gaussian(link="identity"))
summary (modelNH)
tab_df(tidy(Anova(modelNH,type=3)), title = "Consumption (Cold regime)", show.rownames = FALSE, show.type=FALSE,digits=2)

{modelNH=glmmTMB(c ~ Pop*Measure+ (1|Snail_ID)+(1|Is_Replacement)+(1|Egg_tag), REML = FALSE, data=subset(consdataNH,consdataNH$Treat=="NH"),dispformula=~Measure, family=gaussian(link="identity"))
  modelNH2=glmmTMB(c ~ Pop*Measure+ (1|Snail_ID)+(1|Egg_tag), REML = FALSE, data=subset(consdataNH,consdataNH$Treat=="NH"),dispformula=~Measure, family=gaussian(link="identity"))
  modelNH3=glmmTMB(c ~ Pop*Measure+ (1|Snail_ID)+(1|Is_Replacement), REML = FALSE, data=subset(consdataNH,consdataNH$Treat=="NH"),dispformula=~Measure, family=gaussian(link="identity"))
  modelNH4=glmmTMB(c ~ Pop*Measure+ (1|Snail_ID), REML = FALSE, data=subset(consdataNH,consdataNH$Treat=="NH"),dispformula=~Measure, family=gaussian(link="identity"))
}
{anova_df=as.data.frame(anova(modelNH,modelNH2,modelNH3,modelNH4)) 
  anova_df$delta_aic = anova_df$AIC - min(anova_df$AI)
  anova_df=anova_df[,c(1,2,4,9)]
  anova_df <- anova_df[order(anova_df$delta_aic), ]
  tab_df(anova_df,show.rownames=TRUE)}



## Reproductive traits----------

# Load processed data
#reproendstats <- read_csv("processed_data/reproendstats.csv")

# Ensure all variables are coded correctly (important if loaded as processed data)
{
reproendstats$Snail_ID=as.factor(reproendstats$Snail_ID) # Individual snail identifier
reproendstats$Sex=as.factor(reproendstats$Sex) # Male or Female
reproendstats$Treat=as.factor(reproendstats$Treat)  # Thermal regime (cold or warm)
reproendstats$Pop=as.factor(reproendstats$Pop) # Population (NH or NC)
reproendstats$Is_Replacement=as.factor(reproendstats$Is_Replacement) # Whether the snail was a reload added late to the experiment to compesate for early mortality (see methods)
reproendstats$Egg_tag=as.factor(reproendstats$Egg_tag) # Maternal line
}

### Size at first reproduction---------------------------------

# Create model
model=lm(firstrepro_size~ Treat*Pop, data=reproendstats)

# Summary table
tab_df(tidy(Anova(model,type=3)), title = "Size at first reproduction", show.rownames = FALSE, show.type=FALSE,digits=2)

# Checks to see if other random effects structures were preferred (used in supplemental material)
{
model=glmmTMB(firstrepro_size~ Treat*Pop + (1|Egg_tag)+(1|Is_Replacement), data=reproendstats, family=gaussian)
model2=glmmTMB(firstrepro_size~ Treat*Pop + (1|Is_Replacement), data=reproendstats, family=gaussian)
model3=glmmTMB(firstrepro_size~ Treat*Pop + (1|Egg_tag), data=reproendstats, family=gaussian)
model4=glmmTMB(firstrepro_size~ Treat*Pop, data=reproendstats, family=gaussian)}

{anova_df=as.data.frame(anova(model,model2,model3,model4)) 
  anova_df$delta_aic = anova_df$AIC - min(anova_df$AI,na.rm=TRUE)
  anova_df=anova_df[,c(1,2,4,9)]
  anova_df <- anova_df[order(anova_df$delta_aic), ]
  tab_df(anova_df,show.rownames=TRUE)}

### Age at first reproduction--------------------------

# Create model
model=lm(firstrepro_age~ Treat*Pop, data=reproendstats)

# Summary table
tab_df(tidy(Anova(model,type=3)), title = "Age at first reproduction", show.rownames = FALSE, show.type=FALSE,digits=2)

# Checks to see if other random effects structures were preferred (used in supplemental material)
{
model=glmmTMB(firstrepro_age~ Treat*Pop + (1|Egg_tag)+(1|Is_Replacement), data=reproendstats, family=gaussian)
model2=glmmTMB(firstrepro_age~ Treat*Pop + (1|Is_Replacement), data=reproendstats, family=gaussian)
model3=glmmTMB(firstrepro_age~ Treat*Pop + (1|Egg_tag), data=reproendstats, family=gaussian)
model4=glmmTMB(firstrepro_age~ Treat*Pop, data=reproendstats, family=gaussian)
}


{anova_df=as.data.frame(anova(model,model2,model3,model4)) 
  anova_df$delta_aic = anova_df$AIC - min(anova_df$AI,na.rm=TRUE)
  anova_df=anova_df[,c(1,2,4,9)]
  anova_df <- anova_df[order(anova_df$delta_aic), ]
  tab_df(anova_df,show.rownames=TRUE)}


### Total reproductive output-----------------------------------------

# Create model
model=glmmTMB(totalembryos ~ Treat*Pop + (1|Egg_tag), data=reproendstats, ziformula =~1, nbinom1(link = "log"))

# Summary table
tab_df(tidy(Anova(model,type=3)), title = "Total reproductive output", show.rownames = FALSE, show.type=FALSE,digits=2)


# Checks to see if other random effects structures were preferred (used in supplemental material)
{
model=glmmTMB(totalembryos ~ Treat*Pop + (1|Egg_tag)+(1|Is_Replacement), data=reproendstats, ziformula =~1, nbinom1(link = "log"))
model2=glmmTMB(totalembryos ~ Treat*Pop + (1|Egg_tag), data=reproendstats, ziformula =~1, nbinom1(link = "log"))
model3=glmmTMB(totalembryos ~ Treat*Pop + (1|Is_Replacement), data=reproendstats, ziformula =~1, nbinom1(link = "log"))
model4=glmmTMB(totalembryos ~ Treat*Pop, data=reproendstats, ziformula =~1, nbinom1(link = "log"))
}


{anova_df=as.data.frame(anova(model,model2,model3,model4)) 
  anova_df$delta_aic = anova_df$AIC - min(anova_df$AI,na.rm=TRUE)
  anova_df=anova_df[,c(1,2,4,9)]
  anova_df <- anova_df[order(anova_df$delta_aic), ]
  tab_df(anova_df,show.rownames=TRUE)}


### Hatching success-------------------------------------------------------

hatchingsuccess <- read_csv("processed_data/hatchingsuccess.csv")

hatchingsuccess$Treatment=as.factor(hatchingsuccess$Treatment) # Thermal regime (cold or warm)
hatchingsuccess$Population=as.factor(hatchingsuccess$Population) # Population (NH or NC)


model=glm(Prop_Hatchling_success ~ Treatment*Population, data=hatchingsuccess, family=quasibinomial(link = "identity"))

tab_df(tidy(Anova(model,type=3)), title = "Hatching success", show.rownames = FALSE, show.type=FALSE,digits=2)

# No interaction effect, therefore rerun without interaction term (equivelant to running a type-2 SS which is appropriate in the absence of significant interactions)

model=glm(Prop_Hatchling_success ~ Treatment+Population, data=hatchingsuccess, family=quasibinomial(link = "identity"))

tab_df(tidy(Anova(model,type=3)), title = "Hatching success", show.rownames = FALSE, show.type=FALSE,digits=2)






# Summary statistics referred to in text----------------------------------------


# Total number of embryos produced

data_summary(reproendstats, varname="totalembryos", groupnames=c("Pop","Treat"))


table(reproendstats$reproduced,reproendstats$Pop,reproendstats$Treat) 

# Size at first reproduction
data_summary(reproendstatsGeorgia, varname="firstrepro_size", groupnames=c("Pop","Treat"))

# Age at first reproduction
table(reproendstats$reproduced,reproendstats$Pop,reproendstats$Treat) 





