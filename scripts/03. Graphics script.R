# Script 01 needs to be run first to provide required data

# This script contains all the statistical tests from script 02, as well as ggplot code to make the figures.


# Load libraries and custom settings----
{
  library(plyr)
  library(nlme)
  library(lme4)
  library(tidyr)
  library(dplyr)
  library(car)
  library(lubridate)
  library(ggplot2)
  library(ggbreak)
  library(emmeans)
  library(multcomp)
  library(multcompView)
  library(ggrepel)
  library(Matrix)
  library(DHARMa)
  library(performance)
  library(glmmTMB)
  library(data.table)
  library(scales)
  library(tidyverse)
  library(cowplot)
  library(ggbreak)
  
  # set appropriate contrasts for type 3 SS ANOVA
  options(contrasts = c("contr.sum", "contr.poly")) 
  
   theme_mine <- function(base_size = 16) {
    # Starts with theme_grey and then modify some parts
    theme_bw(base_size = base_size) %+replace%
      theme(
        strip.background = element_blank(),
        #strip.text.x = element_text(size = 12),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 12),
        #axis.text.x = element_text(size=12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.ticks =  element_line(colour = "black"), 
        axis.title.x= element_text(size=14),
        axis.title.y= element_text(size=12,angle=90,vjust=4),
        panel.background = element_blank(), 
        #panel.border =element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        panel.spacing = unit(1.0, "lines"), 
        plot.background = element_blank(), 
        plot.margin = unit(c(0.5,  1, 0.5, 0.5), "lines"),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)
      )
  } #function for custom theme       
  
  #Set standard position dodge
  
  dodge=500000
  
  size=2
  
  #Set replacement labels for facets in figures 2 and 3
  regime_names <- c(
    `NC` = "Warm temperature regime",
    `NH` = "Cold temperature regime"
  )
  
  
  
  data_summary=function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
      c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE),
        se = sd(x[[col]]) / sqrt(length(x[[col]])), na.rm=TRUE)
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    return(data_sum)
  }
  #Function for data summary to create error bars and means on plots
  
  
  group.colors <-  data.table(NC = "#FF4E00", NH = "#16BDD8")
  
  
}



## Load temperature data-----------------------


# Lab temperature data
{labtemps <- read_csv("temperature_data/Lab_temperature_data.csv")
  labtemps$time <- as.POSIXct(labtemps$time, format = "%Y/%m/%d %H:%M")
  labtemps$pop = as.factor(labtemps$pop)}

# Field temperature data
{fieldtemps <- read_csv("processed_data/fieldtemps.csv")
fieldtemps$time <- as.POSIXct(fieldtemps$time, format = "%Y/%m/%d %H:%M")
fieldtemps$pop = as.factor(fieldtemps$pop)}

## Temperatures plot----

# Field temperatures indicated by geom_smooth

labtempplot=ggplot(labtemps, aes(y = temp, x = time)) +
  geom_point(aes(y = temp, x = time),    size = 0.2,    alpha = 0.1) +
  geom_smooth(aes(y = temp, x = time, colour = "Grey"),
              data = fieldtemps,
              se = FALSE,
              span = 0.1,) +  theme_mine() +
  theme(strip.text.x = element_text(size = 12)) +
  theme(legend.position = "none") +
  ylab("Temperature (C)") +
  xlab(NULL) +
  facet_grid(cols = vars(pop),labeller = as_labeller(regime_names))+
  scale_colour_manual(values = group.colors)+
    scale_x_datetime(
      breaks = as.POSIXct(seq.Date(
        from = as.Date("2022-08-01"),
        to = as.Date("2023-12-01"),
        by = "2 months"
      )),
      limits = as.POSIXct(c("2022-08-01 12:00", "2023-12-10 16:00")),
      labels = date_format("%b")
    )

labtempplot # Takes a while to load

## Size and Growth Rate figures----

### Shell size figure----

model=glmmTMB(m ~ Treat*Pop*Measure+ (1|Snail_ID)+ (1|Sex)+ (1|Egg_tag), REML = TRUE, data=sizedata, family=gaussian(link="log"))

{model_means=emmeans(model, list(pairwise ~ Pop|Treat|Measure), adjust = "tukey")
  model_means_cld <- cld(object = model_means,
                         adjust = "Tukey",
                         Letters = letters,
                         alpha = 0.05)
  model_means_cld$blank=replace(model_means_cld$.group,model_means_cld$.group=="  b","*")
  model_means_cld$blank2=replace(model_means_cld$blank,model_means_cld$blank==" a","")
  model_means_cld$blank=replace(model_means_cld$blank2,model_means_cld$blank2==" a ","")
  model_means_cld$blank2=replace(model_means_cld$blank,model_means_cld$Pop=="NH","")
  model_means_cld$asterisks=replace(model_means_cld$blank2,model_means_cld$blank=="*","*")
}


lengthplot={growthsummary <- data_summary(data=sizedata, varname="m", groupnames=c("Pop","Treat","Measure"))
timesummary=data.frame(aggregate(time~Measure,data=sizedata,median)) #Create summary of time data
growthsummary=merge(growthsummary,timesummary,by="Measure") #Recombine time data with summary growth data
model_means_cld2=merge(model_means_cld,timesummary,by="Measure") #Recombine time data with summary growth data
ggplot(data=growthsummary,aes(y=mean,x=time, linetype=Pop, shape=Pop))+
  geom_point(data=growthsummary, aes(y=mean,x=time,colour=Pop),position = position_dodge(width=dodge),size=size)+
  geom_text(data = model_means_cld2, aes(label=asterisks,x=time, y=exp(emmean)+1.5),size=4 )+
  geom_linerange(data = growthsummary, linetype="solid", aes(time, mean, ymin = mean - se, ymax = mean + se,colour=Pop),position = position_dodge(width=dodge),size=1)+
  geom_line(colour="black")+
  #geom_smooth(method="loess",se=FALSE, colour="black", span=0.4)+
  
  facet_grid(cols = vars(Treat))+
  theme_mine()+
  ylab("Shell height (mm)")+
  theme(legend.position="none")+
  scale_x_datetime(
    breaks = as.POSIXct(seq.Date(
      from = as.Date("2022-08-01"),
      to = as.Date("2023-12-01"),
      by = "2 months"
    )),
    limits = as.POSIXct(c("2022-08-01 12:00", "2023-12-10 16:00")),
    labels = date_format("%b")
  )+
  xlab(NULL)
}

lengthplot

### Specific growth rate figure ----

growthdata$pgrowth2=growthdata$pgrowth
growthdata$pgrowth2 <- ifelse(growthdata$pgrowth2 < 0, 0, growthdata$pgrowth2)
model=glmmTMB(pgrowth2 ~ Treat*Pop*Measure+ (1|Snail_ID)+ (1|Sex)+ (1|Is_Replacement), REML = TRUE, data=growthdata, family=tweedie(link="log"),dispformula=~Measure*Treat)

{model_means=emmeans(model, list(pairwise ~ Pop|Treat|Measure), adjust = "tukey")
  model_means_cld <- cld(object = model_means,
                         adjust = "Tukey",
                         Letters = letters,
                         alpha = 0.05)
  model_means_cld$blank=replace(model_means_cld$.group,model_means_cld$.group=="  b","*")
  model_means_cld$blank2=replace(model_means_cld$blank,model_means_cld$blank==" a","")
  model_means_cld$blank=replace(model_means_cld$blank2,model_means_cld$blank2==" a ","")
  model_means_cld$blank2=replace(model_means_cld$blank,model_means_cld$Pop=="NH","")
  model_means_cld$asterisks=replace(model_means_cld$blank2,model_means_cld$blank=="*","*")
}


pgrowthplot={growthsummary <- data_summary(growthdata, varname="pgrowth", groupnames=c("Pop","Treat","Measure"))
timesummary=data.frame(aggregate(time~Measure,data=growthdata,median)) #Create summary of time data
growthsummary=merge(growthsummary,timesummary,by="Measure") #Recombine time data with summary growth data
model_means_cld2=merge(model_means_cld,timesummary,by="Measure") #Recombine time data with summary growth data

ggplot(data=growthsummary,aes(y=mean,x=time, linetype=Pop, shape=Pop))+
  geom_point(data=growthsummary, aes(y=mean,x=time,colour=Pop),position = position_dodge(width=dodge),size=size)+
  geom_linerange(data = growthsummary, linetype="solid", aes(time, mean, ymin = mean - se, ymax = mean + se,colour=Pop),position = position_dodge(width=dodge),size=1)+
  geom_text(data = model_means_cld2, aes(label=asterisks,x=time, y=exp(emmean)+0.2),size=4 )+
  
  geom_line(colour="black")+
  #geom_smooth(method="loess",se=FALSE, colour="black", span=0.4)+
  facet_grid(cols = vars(Treat))+
  #scale_y_log10()+
  theme_mine()+
  theme(strip.text.x = element_text(size = 12)) +
  ylab("Specific growth rate
(% g day−1)")+
  theme(legend.position="none")+
  scale_colour_manual(values = group.colors)+
  scale_x_datetime(
    breaks = as.POSIXct(seq.Date(
      from = as.Date("2022-08-01"),
      to = as.Date("2023-12-01"),
      by = "2 months"
    )),
    limits = as.POSIXct(c("2022-08-01 12:00", "2023-12-10 16:00")),
    labels = date_format("%b")
  )+
  facet_grid(cols = vars(Treat),labeller = as_labeller(regime_names))+
  ggbreak::scale_y_break(c(1.3, 2.0), scales = c(0.4, 0.6),space=0.0)+ #Note - tried to add break, but doesnt display with cowplot
  xlab(NULL)}

pgrowthplot

### Weight data figure----
model=glmmTMB(Tissue_weight~ Treat*Pop*Month+ (1|Snail_ID)+(1|Sex), REML = TRUE, data=weightdata, family=Gamma(link="log"))


{model_means=emmeans(model, list(pairwise ~ Pop|Treat|Month), adjust = "tukey")
  model_means_cld <- cld(object = model_means,
                         adjust = "Tukey",
                         Letters = letters,
                         alpha = 0.05)
  model_means_cld$blank=replace(model_means_cld$.group,model_means_cld$.group=="  b","*")
  model_means_cld$blank2=replace(model_means_cld$blank,model_means_cld$blank==" a","")
  model_means_cld$blank=replace(model_means_cld$blank2,model_means_cld$blank2==" a ","")
  model_means_cld$blank2=replace(model_means_cld$blank,model_means_cld$Pop=="NH","")
  model_means_cld$asterisks=replace(model_means_cld$blank2,model_means_cld$blank=="*","*")
}

weightplot={
  weightsummary <- data_summary(weightdata, varname="Tissue_weight", groupnames=c("Pop","Treat","Month"))
  timesummary=data.frame(aggregate(time~Month,data=weightdata,median)) #Create summary of time data
  model_means_cld2=merge(model_means_cld,timesummary,by="Month") #Recombine time data with summary growth data
  
  weightsummary=merge(weightsummary,timesummary,by="Month") #Recombine time data with summary growth data
  ggplot(data=weightsummary,aes(y=mean,x=time, linetype=Pop, shape=Pop))+
    geom_point(data=weightsummary, aes(y=mean,x=time,colour=Pop),position = position_dodge(width=dodge),size=size)+
    geom_linerange(data = weightsummary, linetype="solid", aes(time, mean, ymin = mean - se, ymax = mean + se,colour=Pop),position = position_dodge(width=dodge),size=1)+
    geom_line(colour="black")+
    geom_text(data = model_means_cld2, aes(label=asterisks,x=time, y=exp(emmean+0.3)),size=4 )+
    #geom_smooth(method="loess",se=FALSE, colour="black", span=0.4)+
    facet_grid(cols = vars(Treat))+
    theme_mine()+
    ylab("Tissue weight (g)")+
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(size=12),) +
    scale_colour_manual(values = group.colors)+
    scale_x_datetime(
      breaks = as.POSIXct(seq.Date(
        from = as.Date("2022-08-01"),
        to = as.Date("2023-12-01"),
        by = "2 months"
      )),
      limits = as.POSIXct(c("2022-08-01 12:00", "2023-12-10 16:00")),
      labels = date_format("%b")
    )+
    xlab(NULL)
}

weightplot


## Consumption data figure ----

# Create separate data frame for NH data (see script 02 for details)

consdataNH <- subset(consdata,  !(Measure %in% c("3","4","5", "6","7", "15")))


modelNC=glmmTMB(c ~ Pop*Measure+ (1|Snail_ID), REML = FALSE, data=subset(consdata,consdata$Treat=="NC"),dispformula=~Measure, family=gaussian(link="identity"))

modelNH=glmmTMB(c ~ Pop*Measure+ (1|Snail_ID), REML = FALSE, data=subset(consdataNH,consdataNH$Treat=="NH"),dispformula=~Measure, family=gaussian(link="identity"))



{model_means=emmeans(modelNC, list(pairwise ~ Pop|Measure), adjust = "tukey")
  model_means_cld <- cld(object = model_means,
                         adjust = "Tukey",
                         Letters = letters,
                         alpha = 0.05)
  model_means_cld$blank=replace(model_means_cld$.group,model_means_cld$.group=="  b","*")
  model_means_cld$blank2=replace(model_means_cld$blank,model_means_cld$blank==" a","")
  model_means_cld$blank=replace(model_means_cld$blank2,model_means_cld$blank2==" a ","")
  model_means_cld$blank2=replace(model_means_cld$blank,model_means_cld$Pop=="NH","")
  model_means_cld$asterisks=replace(model_means_cld$blank2,model_means_cld$blank=="*","*")
  model_means_cldNC=model_means_cld
  model_means_cldNC$Treat="NC"
  
}

{model_means=emmeans(modelNH, list(pairwise ~ Pop|Measure), adjust = "tukey")
  model_means_cld <- cld(object = model_means,
                         adjust = "Tukey",
                         Letters = letters,
                         alpha = 0.05)
  model_means_cld$blank=replace(model_means_cld$.group,model_means_cld$.group=="  b","*")
  model_means_cld$blank2=replace(model_means_cld$blank,model_means_cld$blank==" a","")
  model_means_cld$blank=replace(model_means_cld$blank2,model_means_cld$blank2==" a ","")
  model_means_cld$blank2=replace(model_means_cld$blank,model_means_cld$Pop=="NH","")
  model_means_cld$asterisks=replace(model_means_cld$blank2,model_means_cld$blank=="*","*")
  model_means_cldNH=model_means_cld
  model_means_cldNH$Treat="NH"
  
}

model_means_cld=rbind(model_means_cldNH,model_means_cldNC)

consplot = { conssummary <- data_summary(data=consdata, varname="c", groupnames=c("Pop","Treat","Measure"))
timesummary=data.frame(aggregate(time~Measure,data=consdata,median)) #Create summary of time data
model_means_cld2=merge(model_means_cld,timesummary,by="Measure") #Recombine time data with summary growth data
conssummary=merge(conssummary,timesummary,by="Measure") #Recombine time data with summary growth data
consplot = ggplot(data=conssummary,aes(y=mean,x=time, linetype=Pop, shape=Pop))+
  geom_point(data=conssummary, aes(y=mean,x=time,colour=Pop),position = position_dodge(width=dodge),size=size)+
  geom_text(data = model_means_cld2, aes(label=asterisks,x=time, y=emmean+0.2),size=4 )+
  geom_linerange(data = conssummary, linetype="solid", aes(time, mean, ymin = mean - se, ymax = mean + se,colour=Pop),position = position_dodge(width=dodge),size=1)+
  geom_line(colour="black")+
  #geom_smooth(method="loess",se=FALSE, colour="black", span=0.4)+
  facet_grid(cols = vars(Treat))+
  theme_mine()+
  ylab("Feeding rate
(oysters day-1)")+
  theme(legend.position="none")+
  scale_colour_manual(values = group.colors)+
  scale_x_datetime(
    breaks = as.POSIXct(seq.Date(
      from = as.Date("2022-08-01"),
      to = as.Date("2023-12-01"),
      by = "2 months"
    )),
    limits = as.POSIXct(c("2022-08-01 12:00", "2023-12-10 16:00")),
    labels = date_format("%b")
  )+
  xlab(NULL)
}

consplot






## Reproductive output figures----
### Reproductive output throughout the experiment----

{

{reprosummary <- data_summary(fulldata, varname="embryos", groupnames=c("Pop","Treat","Measure"))
timesummary=data.frame(aggregate(time~Measure,data=fulldata,median)) #Create summary of time data
reprosummary=merge(reprosummary,timesummary,by="Measure") #Recombine time data with summary growth data
reproplot=ggplot(data=reprosummary,aes(y=mean,x=time, linetype=Pop, shape=Pop,fill=Pop))+
  geom_point(data=reprosummary, aes(y=mean,x=time,colour=Pop),position = position_dodge(width=dodge),size=size)+
  geom_linerange(data = reprosummary, linetype="solid", aes(time, mean, ymin = mean - se, ymax = mean + se,colour=Pop),position = position_dodge(width=dodge),size=1)+
  geom_line(colour="black")+
  facet_grid(cols = vars(Treat))+
  #scale_y_log10()+
  theme_mine()+
  coord_cartesian(ylim=c(0,350))+
  ylab("Embryos produced
       female-1")+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(size=12),plot.margin = unit(c(0.5,  1, 0.5, 0.5), "lines")) +
  scale_colour_manual(values = group.colors)+
  scale_x_datetime(
    breaks = as.POSIXct(seq.Date(
      from = as.Date("2022-08-01"),
      to = as.Date("2023-12-01"),
      by = "2 months"
    )),
    limits = as.POSIXct(c("2022-08-01 12:00", "2023-12-10 16:00")),
    labels = date_format("%b")
  )+
  xlab(NULL)
}






### Size at first reproduction plot ----


model=lm(firstrepro_size ~ Treat*Pop, data=reproendstats)

{model_means=emmeans(model, list(pairwise ~ Pop*Treat), adjust = "tukey")
  model_means_cld=cld(object = model_means,
                      adjust = "Tukey",
                      Letters = letters,
                      alpha = 0.05)
}

firstreprosizeplot=ggplot(data=reproendstats,aes(y=firstrepro_size,x=(Treat), shape=Pop))+
  geom_boxplot(aes(y=(firstrepro_size),x=(Treat),fill=Pop))+
  geom_text(data = model_means_cld, aes(label=.group, y=emmean),size=4 )+
  theme_mine()+
  theme(legend.position="none")+
  ylab("Size at first 
  reproduction")+
  scale_fill_manual(values = group.colors)+
  theme(axis.title.x=element_blank())


### Age at first reproduction plot----

model=lm((firstrepro_age) ~ Treat*Pop, data=reproendstats)

{
model_means=emmeans(model, list(pairwise ~ Pop*Treat), adjust = "tukey")
model_means_cld=cld(object = model_means,
                    adjust = "Tukey",
                    Letters = letters,
                    alpha = 0.05)
}


firstreproageplot=ggplot(data=reproendstats,aes(y=firstrepro_age,x=(Treat), shape=Pop))+
  geom_boxplot(aes(y=(firstrepro_age),x=(Treat),fill=Pop))+
  geom_text(data = model_means_cld, aes(label=.group, y=(emmean)),size=4 )+
  theme_mine()+
  theme(legend.position="none")+
  ylab("Time until first reproduction
       (Days)")+
  
  scale_fill_manual(values = group.colors)
  #theme(axis.title.x=element_blank())





### Total number of embryos produced plot----



model=glmmTMB(totalembryos ~ Treat*Pop + (1|Egg_tag), data=reproendstats, ziformula =~1, nbinom1(link = "log"))


{model_means=emmeans(model, list(pairwise ~ Pop*Treat), adjust = "tukey")
model_means_cld=cld(object = model_means,
                    adjust = "Tukey",
                    Letters = letters,
                    alpha = 0.05)
}

totalembryosplot=ggplot(data=reproendstats,aes(y=totalembryos,x=(Treat), shape=Pop))+
  geom_boxplot(aes(y=(totalembryos),x=(Treat),fill=Pop))+
  geom_text(data = model_means_cld, aes(label=.group, y=exp(emmean)),size=4 )+
  theme_mine()+
  theme(legend.position="none")+
  ylab("Total number of embryos produced")+
  scale_fill_manual(values = group.colors)+
  theme(axis.title.x=element_blank())



### Hatching success plot ----


model=glm(Prop_Hatchling_success ~ Treatment*Population, data=hatchingsuccess, family=quasibinomial(link = "identity"))

model_means=emmeans(model, list(pairwise ~ Population*Treatment), adjust = "tukey")
model_means_cld=cld(object = model_means,
                    adjust = "Tukey",
                    Letters = letters,
                    alpha = 0.05)

hatchingsuccesplot=ggplot(data=hatchingsuccess,aes(y=Prop_Hatchling_success,x=(Treatment), shape=Population))+
  geom_boxplot(aes(y=(Prop_Hatchling_success),x=(Treatment),fill=Population))+
  geom_text(data = model_means_cld, aes(label=.group, y=(1)),size=4 )+
  theme_mine()+
  theme(legend.position="none")+
  #theme(axis.title.x=element_blank())+
  scale_fill_manual(values = group.colors)+
  
  ylab("Hatching success
  (%)")

}
# Create Manuscript figures ----

## Create Figure legend----

{
legend <- get_plot_component(lengthplot + theme(legend.position="bottom", legend.box="horiztonal", legend.title = element_text(size=12),legend.text = element_text(size=10)) +
                               labs(
                                 colour = "Population", shape="Population",linetype="Population"
                               ), 'guide-box-bottom', return_all = TRUE)
}

## Figure 2 ----

fig2plot=cowplot::plot_grid(labtempplot,lengthplot,weightplot, ncol=1, align="v",hjust=-1, labels = c('A', 'B', 'C', 'D'))

cowplot::plot_grid(legend,fig2plot,  nrow=2, rel_heights=c(0.02,1))


## Figure 3 ----

fig3plot=cowplot::plot_grid(consplot,reproplot, ncol=1, align="v",hjust=-1, labels = c('B', 'C'))

fig3plot

fig3plot=cowplot::plot_grid(print(pgrowthplot),fig3plot,  nrow=2, rel_heights=c(1.2,2))


cowplot::plot_grid(legend,fig3plot,  nrow=2, rel_heights=c(0.02,1))

## Figure 4----

{
  legend <- get_plot_component(hatchingsuccesplot + theme(legend.position="bottom", legend.box="horiztonal", legend.title = element_text(size=12),legend.text = element_text(size=10)) +
                                 labs(
                                   colour = "Population", shape="Population",linetype="Population"
                                 ), 'guide-box-bottom', return_all = TRUE)
}



fig4plot=cowplot::plot_grid(totalembryosplot,firstreprosizeplot,firstreproageplot,hatchingsuccesplot, ncol=2, align="h", labels = c('A', 'B', 'C', 'D'))


cowplot::plot_grid(legend,fig4plot,  nrow=2, rel_heights=c(0.05,1))


