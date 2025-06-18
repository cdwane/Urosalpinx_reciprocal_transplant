# This script reproduces plots contained in supplement 2
# including the Goergia population which was excluded from our main statistical analysis

# Script 01 needs to be run first to load the required data


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
  
  
  group.colors <-  data.table(NC = "#FF4E00", NH = "#16BDD8", GA = "#F4C811")
  
  
}



## Size and Growth Rate figures----


### Shell size figure----

lengthplot={growthsummary <- data_summary(data=sizedatageorgia, varname="m", groupnames=c("Pop","Treat","Measure"))
timesummary=data.frame(aggregate(time~Measure,data=sizedatageorgia,median)) #Create summary of time data
growthsummary=merge(growthsummary,timesummary,by="Measure") #Recombine time data with summary growth data
ggplot(data=growthsummary,aes(y=mean,x=time, linetype=Pop, shape=Pop))+
  geom_point(data=growthsummary, aes(y=mean,x=time,colour=Pop),position = position_dodge(width=dodge),size=size)+
  geom_linerange(data = growthsummary, linetype="solid", aes(time, mean, ymin = mean - se, ymax = mean + se,colour=Pop),position = position_dodge(width=dodge),size=1)+
  geom_line(colour="black")+
  scale_colour_manual(values = group.colors)+
  facet_grid(cols = vars(Treat))+
  theme_mine()+
  theme(strip.text.x = element_text(size = 12)) +
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
  facet_grid(cols = vars(Treat),labeller = as_labeller(regime_names))+
  xlab(NULL)
}

lengthplot

### Specific growth rate figure ----


pgrowthplot={growthsummary <- data_summary(growthdatageorgia, varname="pgrowth", groupnames=c("Pop","Treat","Measure"))
timesummary=data.frame(aggregate(time~Measure,data=growthdatageorgia,median)) #Create summary of time data
growthsummary=merge(growthsummary,timesummary,by="Measure") #Recombine time data with summary growth data
ggplot(data=growthsummary,aes(y=mean,x=time, linetype=Pop, shape=Pop))+
  geom_point(data=growthsummary, aes(y=mean,x=time,colour=Pop),position = position_dodge(width=dodge),size=size)+
  geom_linerange(data = growthsummary, linetype="solid", aes(time, mean, ymin = mean - se, ymax = mean + se,colour=Pop),position = position_dodge(width=dodge),size=1)+
  geom_line(colour="black")+
  facet_grid(cols = vars(Treat))+
  theme_mine()+
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
  ggbreak::scale_y_break(c(1.3, 2.0), scales = c(0.4, 0.6),space=0.0)+ #Note - tried to add break, but doesnt display with cowplot
  xlab(NULL)}

pgrowthplot

## Consumption data figure ----



consplot = {conssummary <- data_summary(data=consdataGeorgia, varname="c", groupnames=c("Pop","Treat","Month"))
timesummary=data.frame(aggregate(time~Month,data=consdataGeorgia,median)) #Create summary of time data
conssummary=merge(conssummary,timesummary,by="Month") #Recombine time data with summary growth data
ggplot(data=conssummary,aes(y=mean,x=time, linetype=Pop, shape=Pop))+
  geom_point(data=conssummary, aes(y=mean,x=time,colour=Pop),position = position_dodge(width=dodge),size=size)+
  geom_linerange(data = conssummary, linetype="solid", aes(time, mean, ymin = mean - se, ymax = mean + se,colour=Pop),position = position_dodge(width=dodge),size=1)+
  geom_line(colour="black")+
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


# How many females reproduced by the end of the experiment?


table(reproendstatsGeorgia$reproduced,reproendstatsGeorgia$Pop,reproendstatsGeorgia$Treat) 

#All Georgia females in NC treatment reproduced, none in NH treatment reproduced (but sample sizes only 11 and 7 respectively)



{reprosummary <- data_summary(fulldatageorgia, varname="embryos", groupnames=c("Pop","Treat","Measure"))
timesummary=data.frame(aggregate(time~Measure,data=fulldatageorgia,median)) #Create summary of time data
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

reproplot




### Size at first reproduction plot ----



reproendstatsGeorgia=subset(reproendstatsGeorgia,!is.na(reproendstatsGeorgia$firstrepro_size))


firstreprosizeplot=ggplot(data=reproendstatsGeorgia,aes(y=firstrepro_size,x=(Treat), shape=Pop))+
  geom_boxplot(aes(y=(firstrepro_size),x=(Treat),fill=Pop))+
  theme_mine()+
  theme(legend.position="none")+
  ylab("Size at first 
  reproduction")+
  scale_fill_manual(values = group.colors)+
  theme(axis.title.x=element_blank())


### Age at first reproduction plot----


firstreproageplot=ggplot(data=reproendstatsGeorgia,aes(y=firstrepro_age,x=(Treat), shape=Pop))+
  geom_boxplot(aes(y=(firstrepro_age),x=(Treat),fill=Pop))+
  theme_mine()+
  theme(legend.position="none")+
  ylab("Time until first reproduction
       (Days)")+
  scale_fill_manual(values = group.colors)+
  theme(axis.title.x=element_blank())





### Total number of embryos produced plot----

totalembryosplot=ggplot(data=reproendstatsGeorgia,aes(y=totalembryos,x=(Treat), shape=Pop))+
  geom_boxplot(aes(y=(totalembryos),x=(Treat),fill=Pop))+
  theme_mine()+
  theme(legend.position="none")+
  ylab("Total number of embryos produced")+
  scale_fill_manual(values = group.colors)+
  theme(axis.title.x=element_blank())


### Hatching success plot ----


hatchingsuccesplot=ggplot(data=hatchingsuccessGeorgia,aes(y=Prop_Hatchling_success,x=(Treatment), shape=Population))+
  geom_boxplot(aes(y=(Prop_Hatchling_success),x=(Treatment),fill=Population))+
  theme_mine()+
  theme(legend.position="none")+
  theme(axis.title.x=element_blank())+
  scale_fill_manual(values = group.colors)+
  
  ylab("Hatching success
  (%)")


# Create Manuscript figures ----

## Create Figure legend----

{
  legend <- get_plot_component(lengthplot + theme(legend.position="bottom", legend.box="horiztonal", legend.title = element_text(size=12),legend.text = element_text(size=10)) +
                                 labs(
                                   colour = "Population", shape="Population",linetype="Population"
                                 ), 'guide-box-bottom', return_all = TRUE)
}


## Figure S1 ----

{
figS1plot=cowplot::plot_grid(consplot,reproplot, ncol=1, align="v", labels = c('C','D', 'A'))
figS1plot=cowplot::plot_grid(print(pgrowthplot),fig3plot,  nrow=2, rel_heights=c(1.1,2))
figS1plot=cowplot::plot_grid(print(lengthplot),fig3plot,  align="v", nrow=2, rel_heights=c(1.1,3))
}
cowplot::plot_grid(legend,figS1plot,  nrow=2, rel_heights=c(0.02,1))

# Need to manually adjust alignments and add labels

## Figure S2----

{
  legend <- get_plot_component(hatchingsuccesplot + theme(legend.position="bottom", legend.box="horiztonal", legend.title = element_text(size=12),legend.text = element_text(size=10)) +
                                 labs(
                                   colour = "Population", shape="Population",linetype="Population"
                                 ), 'guide-box-bottom', return_all = TRUE)
}



figS2plot=cowplot::plot_grid(totalembryosplot,firstreprosizeplot,firstreproageplot,hatchingsuccesplot, ncol=2, align="hv", labels = c('A', 'B', 'C', 'D'))


cowplot::plot_grid(legend,figS2plot,  nrow=2, rel_heights=c(0.05,1))


