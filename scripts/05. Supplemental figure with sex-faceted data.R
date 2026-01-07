
# This script generates supplemental figure presenting length, weight and growth data faceted by sex. This is included in Appendix S2: Supplemental Figure S3.

# This script contains all the statistical tests from script 02, as well as ggplot code to make the figures.

# Script 01 needs to be run first to provide required data

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





## Size and Growth Rate figures----

### Shell size figure----

model=glmmTMB(m ~ Treat*Pop*Measure*Sex+ (1|Snail_ID), REML = TRUE, data=sizedata, family=gaussian(link="identity"))



{model_means=emmeans(model, list(pairwise ~ Pop|Treat|Measure|Sex), adjust = "tukey")
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


lengthplot2={growthsummary <- data_summary(data=sizedata, varname="m", groupnames=c("Pop","Treat","Measure","Sex"))
timesummary=data.frame(aggregate(time~Measure,data=sizedata,median)) #Create summary of time data
growthsummary=merge(growthsummary,timesummary,by="Measure") #Recombine time data with summary growth data
model_means_cld2=merge(model_means_cld,timesummary,by="Measure") #Recombine time data with summary growth data
ggplot(data=growthsummary,aes(y=mean,x=time, linetype=Pop, shape=Pop))+
  geom_point(data=growthsummary, aes(y=mean,x=time,colour=Pop),position = position_dodge(width=dodge),size=size)+
  geom_text(data = model_means_cld2, aes(label=asterisks,x=time, y=(emmean)+2),size=4 )+
  geom_linerange(data = growthsummary, linetype="solid", aes(time, mean, ymin = mean - se, ymax = mean + se,colour=Pop),position = position_dodge(width=dodge),size=1)+
  geom_line(colour="black")+
  #geom_smooth(method="loess",se=FALSE, colour="black", span=0.4)+
  facet_grid(cols = vars(Treat), rows= vars(Sex))+
  theme_mine()+
  theme(strip.text.x = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  ylab("Shell height (mm)")+
  theme(legend.position="none")+
  scale_x_datetime(
    breaks = as.POSIXct(seq.Date(
      from = as.Date("2022-08-01"),
      to = as.Date("2023-12-01"),
      by = "2 months"
    )),
    limits = as.POSIXct(c("2022-08-01 00:00", "2023-12-10 16:00")),
    labels = date_format("%b")
  )+
  xlab(NULL)
}

lengthplot2


### Specific growth rate figure ----

model=glmmTMB(pgrowth~ Treat*Pop*Measure*Sex+ (1|Snail_ID), REML = TRUE, data=growthdata, family=gaussian(link="identity"),dispformula = ~ Measure*Treat*Pop,control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)))

{model_means=emmeans(model, list(pairwise ~ Pop|Treat|Measure|Sex), adjust = "tukey")
  model_means_cld <- cld(object = model_means,
                         adjust = "Tukey",
                         Letters = letters,
                         alpha = 0.05)
  model_means_cld$blank=replace(model_means_cld$.group,model_means_cld$.group=="  b","*")
  model_means_cld$blank2=replace(model_means_cld$blank,model_means_cld$blank==" a","")
  model_means_cld$blank=replace(model_means_cld$blank2,model_means_cld$blank2==" a ","")
  model_means_cld$blank2=replace(model_means_cld$blank,model_means_cld$Pop=="NH","")
  model_means_cld$asterisks=replace(model_means_cld$blank2,model_means_cld$blank=="*","*")
} #Convoluted script to make asterisks only cld

# Plot for pgrowth includes broken y-axis - will tweak in illustrator for final publication version

pgrowthplot2={growthsummary <- data_summary(growthdata, varname="pgrowth", groupnames=c("Pop","Treat","Measure","Sex"))
timesummary=data.frame(aggregate(time~Measure,data=growthdata,median)) #Create summary of time data
growthsummary=merge(growthsummary,timesummary,by="Measure") #Recombine time data with summary growth data
model_means_cld2=merge(model_means_cld,timesummary,by="Measure") #Recombine time data with summary growth data

ggplot(data=growthsummary,aes(y=mean,x=time, linetype=Pop, shape=Pop))+
  geom_point(data=growthsummary, aes(y=mean,x=time,colour=Pop),position = position_dodge(width=dodge),size=size)+
  geom_linerange(data = growthsummary, linetype="solid", aes(time, mean, ymin = mean - se, ymax = mean + se,colour=Pop),position = position_dodge(width=dodge),size=1)+
  geom_text(data = model_means_cld2, aes(label=asterisks,x=time, y=(emmean)+0.2),size=4 )+
  
  geom_line(colour="black")+
  #geom_smooth(method="loess",se=FALSE, colour="black", span=0.4)+
  #scale_y_log10()+
  theme_mine()+
  theme(axis.text.x = element_text(size = 12)) +
  ylab("Specific growth rate (% g day−1)")+
  theme(legend.position="none")+
  scale_colour_manual(values = group.colors)+
  scale_x_datetime(
    breaks = as.POSIXct(seq.Date(
      from = as.Date("2022-08-01"),
      to = as.Date("2023-12-01"),
      by = "2 months"
    )),
    limits = as.POSIXct(c("2022-08-01 00:00", "2023-12-10 16:00")),
    labels = date_format("%b")
  ) +
  facet_grid(rows = vars(Sex), cols = vars(Treat))+
  scale_y_break(c(1.3, 2.0), scales = c(0.4, 0.6),space=0.0)+ #Note - tried to add break, but doesnt display with cowplot
  xlab(NULL)}


pgrowthplot2

### Weight data figure----
model=glmmTMB(Tissue_weight~ Treat*Pop*Month*Sex+ (1|Snail_ID), REML = TRUE, data=weightdata, family=Gamma(link="log"))


{model_means=emmeans(model, list(pairwise ~ Pop|Treat|Month|Sex), adjust = "tukey")
  model_means_cld <- cld(object = model_means,
                         adjust = "Tukey",
                         Letters = letters,
                         alpha = 0.05)
  model_means_cld$blank=replace(model_means_cld$.group,model_means_cld$.group=="  b","*")
  model_means_cld$blank2=replace(model_means_cld$blank,model_means_cld$blank==" a","")
  model_means_cld$blank=replace(model_means_cld$blank2,model_means_cld$blank2==" a ","")
  model_means_cld$blank2=replace(model_means_cld$blank,model_means_cld$Pop=="NH","")
  model_means_cld$asterisks=replace(model_means_cld$blank2,model_means_cld$blank=="*","*")
} #Convoluted script to make asterisks only cld

weightplot2={
  weightsummary <- data_summary(weightdata, varname="Tissue_weight", groupnames=c("Pop","Treat","Month","Sex"))
  timesummary=data.frame(aggregate(time~Month,data=weightdata,median)) #Create summary of time data
  model_means_cld2=merge(model_means_cld,timesummary,by="Month") #Recombine time data with summary growth data
  
  weightsummary=merge(weightsummary,timesummary,by="Month") #Recombine time data with summary growth data
  ggplot(data=weightsummary,aes(y=mean,x=time, linetype=Pop, shape=Pop))+
    geom_point(data=weightsummary, aes(y=mean,x=time,colour=Pop),position = position_dodge(width=dodge),size=size)+
    geom_linerange(data = weightsummary, linetype="solid", aes(time, mean, ymin = mean - se, ymax = mean + se,colour=Pop),position = position_dodge(width=dodge),size=1)+
    geom_line(colour="black")+
    geom_text(data = model_means_cld2, aes(label=asterisks,x=time, y=exp(emmean+0.3)),size=4 )+
    #geom_smooth(method="loess",se=FALSE, colour="black", span=0.4)+
    facet_grid(cols = vars(Treat),rows=vars(Sex))+
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
      limits = as.POSIXct(c("2022-08-01 00:00", "2023-12-10 16:00")),
      labels = date_format("%b")
    )+
    xlab(NULL)
}  

weightplot2


{
  legend <- get_plot_component(lengthplot2 +
                                 scale_colour_manual(values = group.colors, labels = c("Low latitude", "High latitude")) +
                                 scale_shape_manual(values = c(16, 17), labels = c("Low latitude", "High latitude")) +
                                 scale_linetype_manual(values = c("solid", "dashed"), labels = c("Low latitude", "High latitude")) +
                                 theme(
                                   legend.position = "bottom",
                                   legend.box = "horizontal",
                                   legend.title = element_text(size = 12),
                                   legend.text = element_text(size = 12)
                                 ) +
                                 labs(colour = "Population", shape = "Population", linetype = "Population"), 'guide-box-bottom', return_all = TRUE)
}



sexplot=cowplot::plot_grid(lengthplot2,weightplot2, ncol=1, align="v",hjust=-1, labels = c('A', 'B'))

sexplot=cowplot::plot_grid(sexplot,print(pgrowthplot2),  nrow=2, rel_heights=c(2,1.4))


cowplot::plot_grid(legend,sexplot,  nrow=2, rel_heights=c(0.02,1))

