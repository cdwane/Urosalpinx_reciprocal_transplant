
#Load and prepare data for analysis-----


{ #Run here to start entire script


  
  # Or Run here to load each section seperately
## Load required libraries-----
{
{library(plyr)
 
library(dplyr)
 library(tidyr)
  library(readr)
  library(readxl)
  library(lubridate)
  
  # Create data_summary function 
  data_summary=function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
      c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE),
        se = sd(x[[col]]) / sqrt(length(x[[col]])), na.rm=TRUE)
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    return(data_sum)}
  
  
}
  
}
  
 

##  Load and prepare shell height / growth rate data for analysis --------------------------------
  { 
### Load raw data---------------------------------------------------------------------------------
{
  RTdata <- read_excel("raw_data/RTdata.xlsx")
}

### Create metadata file -------------------------------------------------------------------------
{
  RT_metadata<-as.data.frame(cbind(RTdata$Snail_ID,RTdata$Hatch_tag,RTdata$Alive,RTdata$Is_Replacement,RTdata$Treat,RTdata$Pop,RTdata$Egg_tag,RTdata$Condo,RTdata$Old_Condo,RTdata$Was_reloaded,RTdata$Death_date,RTdata$Sex,RTdata$transplanted.between.experiments.))
  colnames(RT_metadata)<-c("Snail_ID","Hatch_tag","Alive","Is_Replacement","Treat","Pop",
                           "Egg_tag","Condo","Old_Condo","Was_reloaded","Death_date","Sex","Is.Transplant")
  RT_metadata$Death_date=RTdata$Death_date
  RT_metadata=subset(RT_metadata,RT_metadata$Egg_tag!="NA") #Removes field individuals from metadata
  }

### Grab all measurement data and respective timepoints, remove a couple we don't use, and concatenate ------------------------------------------------------    
    
  {
  {
  RT_measurements <- cbind(RTdata$Snail_ID,RTdata[ , grepl("_mm", colnames(RTdata))])  #Grab all measurement data from the file
  
  RT_measurements=RT_measurements[ -c(8,12) ] # Remove data measured on 4/17/23 and 1/25/23 - we have two sets of measurements for both these months so removing this data makes gaps between measurement timepoints more consistent
  
  
  colnames(RT_measurements)<-c("Snail_ID","m.1","m.2","m.3","m.4","m.5","m.6","m.7","m.8","m.9","m.10","m.11","m.12","m.13","m.14","m.15","m.16","m.17")

  RT_measurements <- RT_measurements %>% mutate(across(c('m.1':'m.17'), as.numeric))
  
  RT_timepoints <- cbind(RTdata$Snail_ID,RTdata$Date_IN,RTdata[ , grepl("_date", colnames(RTdata))]) #Grab all timepoints for respective measurement data from the file
  
  RT_timepoints=RT_timepoints[ -c(10,14) ] # Remove data measured on 4/17/23 and 1/25/23 - we have two sets of measurements for both these months so removing this data makes gaps between measurement timepoints more consistent
  
  colnames(RT_timepoints)<-c("Snail_ID","time.1","reload_date","death_date","time.2","time.3","time.4","time.5","time.6","time.7","time.8","time.9","time.10","time.11","time.12","time.13","time.14","time.15","time.16","time.17")
  RT_timepoints = subset(RT_timepoints, select = -c(reload_date,death_date) )
  
  RT_timepoints <- RT_timepoints %>% mutate(across(c('time.1':'time.17'), as.POSIXct,format="%m/%d/%Y"))
  growthdata=merge(RT_measurements, RT_timepoints, by= "Snail_ID")  
}

### Clean up data by removing measurements likely to be erroneous---------------------------------------------------  
  
  #Exclude a few values likely to be errors - these appear on the original (handwritten) datasheets and are not typos but probably represent calliper user error
  
  growthdata[167,12] <- NA # NCNH721_314_40 Measure 11 length drops by 5mm vs previous recording but then goes back down the next month
  growthdata[266,10] <- NA # NHNC442_300_1 Measure 9 length goes up by 4mm relative to previous reading but then goes back down the next month
  growthdata[342,15] <- NA  # NHNH418_317_12 Measure 14 length goes up by 3mm relative to previous reading but then goes back down the next month

  growthdata[85,10] <- NA  # NCNC469_301_14 Measure 9 length goes down by 2mm relative to previous reading but then goes back down the next month

  growthdata[292,5] <- NA    # NHNC476_304_21, measurement 4 is a substantial outlier in growth analysis

### Calculate/ specific growth rate between each time point---------------------------------------------------------

#Following Lonthair et al. 2024

  {
    {growthdata$pgrowth.2= 100*(exp((log (growthdata$m.2)-log(growthdata$m.1))/(as.numeric(growthdata$time.2-growthdata$time.1)))-1)
    growthdata$pgrowth.3= 100*(exp((log (growthdata$m.3)-log(growthdata$m.2))/(as.numeric(growthdata$time.3-growthdata$time.2)))-1)
    growthdata$pgrowth.4= 100*(exp((log (growthdata$m.4)-log(growthdata$m.3))/(as.numeric(growthdata$time.4-growthdata$time.3)))-1)
    growthdata$pgrowth.5= 100*(exp((log (growthdata$m.5)-log(growthdata$m.4))/(as.numeric(growthdata$time.5-growthdata$time.4)))-1)
    growthdata$pgrowth.6= 100*(exp((log (growthdata$m.6)-log(growthdata$m.5))/(as.numeric(growthdata$time.6-growthdata$time.5)))-1)
    growthdata$pgrowth.7= 100*(exp((log (growthdata$m.7)-log(growthdata$m.6))/(as.numeric(growthdata$time.7-growthdata$time.6)))-1)
    growthdata$pgrowth.8= 100*(exp((log (growthdata$m.8)-log(growthdata$m.7))/(as.numeric(growthdata$time.8-growthdata$time.7)))-1)
    growthdata$pgrowth.9= 100*(exp((log (growthdata$m.9)-log(growthdata$m.8))/(as.numeric(growthdata$time.9-growthdata$time.8)))-1)
    growthdata$pgrowth.10= 100*(exp((log (growthdata$m.10)-log(growthdata$m.9))/(as.numeric(growthdata$time.10-growthdata$time.9)))-1)
    growthdata$pgrowth.11= 100*(exp((log (growthdata$m.11)-log(growthdata$m.10))/(as.numeric(growthdata$time.11-growthdata$time.10)))-1)
    growthdata$pgrowth.12= 100*(exp((log (growthdata$m.12)-log(growthdata$m.11))/(as.numeric(growthdata$time.12-growthdata$time.11)))-1)
    growthdata$pgrowth.13= 100*(exp((log (growthdata$m.13)-log(growthdata$m.12))/(as.numeric(growthdata$time.13-growthdata$time.12)))-1)
    growthdata$pgrowth.14= 100*(exp((log (growthdata$m.14)-log(growthdata$m.13))/(as.numeric(growthdata$time.14-growthdata$time.13)))-1)
    growthdata$pgrowth.15= 100*(exp((log (growthdata$m.15)-log(growthdata$m.14))/(as.numeric(growthdata$time.15-growthdata$time.14)))-1)
    growthdata$pgrowth.16= 100*(exp((log (growthdata$m.16)-log(growthdata$m.15))/(as.numeric(growthdata$time.16-growthdata$time.15)))-1)
    growthdata$pgrowth.17= 100*(exp((log (growthdata$m.17)-log(growthdata$m.16))/(as.numeric(growthdata$time.17-growthdata$time.16)))-1)
    }
} 






 ### Convert to long format and restore metadata----------------------------------------------------------------------------------------------

{growthdata=pivot_longer(growthdata, cols = c("m.1":"time.17","pgrowth.2":"pgrowth.17"), names_to = c(".value", "Measure"), names_sep = "\\.")

growthdata=merge(growthdata, RT_metadata, by= "Snail_ID") 
growthdata=subset(growthdata,growthdata$m!="NA") # Remove any rows where there is no growth measurement
}


#Ensure all columns are formatted correctly

{growthdata$Measure=as.numeric(growthdata$Measure)
  growthdata$Snail_ID=as.factor(growthdata$Snail_ID)  
  growthdata$Treat=as.factor(growthdata$Treat)
  growthdata$Pop=as.factor(growthdata$Pop)
  growthdata$Egg_tag=as.factor(growthdata$Egg_tag)
  growthdata$Alive=as.factor(growthdata$Alive)
  growthdata$Sex=as.factor(growthdata$Sex)
  growthdata$Measure=as.factor(growthdata$Measure)
  growthdata$Is_Replacement=as.factor(growthdata$Is_Replacement)}

  growthdataallindividuals=growthdata # Create backup version that includes all individuals including those tht died during the experiment
  
 growthdata=subset(growthdata,growthdata$Alive=="Y") # remove any individuals that died during the experiment
 
  {growthdata$Sex=as.character(growthdata$Sex) # Correct the number of Sex factors
    growthdata$Sex=as.factor(growthdata$Sex)}

 
 
 
### Create separate datafamea for size data, and for growth data (because I need to remove the NA values for the growth data as it can't be calculated for the first timepoint)----
  
 
{
 sizedatageorgia=growthdata  # Saves a version with Georgia (for use in supplemental material)
 sizedata=subset(growthdata,growthdata$Pop!="GA") #remove Georgia
 {sizedata$Pop=as.character(sizedata$Pop) #Correct number of pop factors (Because Georgia was removed form this version)
   sizedata$Pop=as.factor(sizedata$Pop)}
 {sizedata$Snail_ID=as.character(sizedata$Snail_ID) #Correct number of snails (Because Georgia was removed form this version)
   sizedata$Snail_ID=as.factor(sizedata$Snail_ID)}
}
 
 
 {
   growthdata$Measure=as.numeric(growthdata$Measure)
   growthdata=subset(growthdata,!is.na(growthdata$pgrowth)) #Removes all NA values for pgrowth (where one of the two values used to calculate pgrowth was missing)
   growthdata$Measure=as.factor(growthdata$Measure)
   
   growthdatageorgia=growthdata # Saves a version with Georgia (for use in supplemental material)

   growthdata=subset(growthdata,growthdata$Pop!="GA") #remove Georgia
   {growthdata$Pop=as.character(growthdata$Pop) #Correct number of pop factors (Because Georgia was removed form this version)
     growthdata$Pop=as.factor(growthdata$Pop)}
   {growthdata$Snail_ID=as.character(growthdata$Snail_ID) #Correct number of snails (Because Georgia was removed form this version)
     growthdata$Snail_ID=as.factor(growthdata$Snail_ID)}
 }


  }  }
  
   # Note: If correctly loaded, growthdata should contain 1787 rows, growthdatageorgia should contain 2172 rows
 # sizedata should contain 1918 rows, sizedatageorgia should contain 2331 rows
  
  #write_csv(sizedata, "processed_data/sizedata.csv", na = "")
  #write_csv(growthdata, "processed_data/growthdata.csv", na = "")
  
  #write_csv(sizedatageorgia, "processed_data/Georgia/sizedatageorgia.csv", na = "")
  #write_csv(growthdatageorgia, "processed_data/Georgia/growthdatageorgia.csv", na = "")
  
## Load and prepare consumption data for analysis ----------------------------------------------------------------------------------------
  
### Load Raw data again -----------------------------------------------------------------------------------------------------------------
  
{    
  RTdata <- read_excel("raw_data/RTdata.xlsx")
}
  
    
    ### Grab all consumption data and respective timepoints ------------------------------------------------------    

    {
  RT_consumption <- cbind(RTdata$Snail_ID,RTdata[ , grepl("_cons", colnames(RTdata))]) 
  colnames(RT_consumption)<-c("Snail_ID","c.1","c.2","c.3","c.4","c.5","c.6","c.7","c.8","c.9","c.10","c.11","c.12","c.13","c.14","c.15")
  RT_consumption <- RT_consumption %>% mutate(across(c('c.1':'c.15'), as.numeric))
  
  RT_dates <- cbind(RTdata$Snail_ID,RTdata[ , grepl("_loaded", colnames(RTdata))]) 
    colnames(RT_dates)<-c("Snail_ID","time.1","time.2","time.3","time.4","time.5","time.6","time.7","time.8","time.9","time.10","time.11","time.12","time.13","time.14","time.15")
    RT_dates <- RT_dates %>% mutate(across(c('time.1':'time.15'), as.POSIXct,format="%m/%d/%Y"))
  consdata=merge(RT_dates, RT_consumption, by= "Snail_ID")  
    }
    
  ### Create metadata ------------------------------------------------------    
  
  {
    RT_metadata<-as.data.frame(cbind(RTdata$Snail_ID,RTdata$Hatch_tag,RTdata$Alive,RTdata$Is_Replacement,RTdata$Treat,RTdata$Pop,RTdata$Egg_tag,RTdata$Condo,RTdata$Old_Condo,RTdata$Was_reloaded,RTdata$Death_date,RTdata$Sex,RTdata$transplanted.between.experiments.))
    colnames(RT_metadata)<-c("Snail_ID","Hatch_tag","Alive","Is_Replacement","Treat","Pop",
                             "Egg_tag","Condo","Old_Condo","Was_reloaded","Death_date","Sex","Is.Transplant")
    RT_metadata$Death_date=RTdata$Death_date
    RT_metadata=subset(RT_metadata,RT_metadata$Egg_tag!="NA") #Removes field individuals from metadata
  }
  
  
  
  consdata=merge(consdata, RT_metadata, by= "Snail_ID")
  
  
  ### Remove measurements which were taken from single individuals, past the point at which snails were paired---------------------------------------------------  

  { consdata[151,25:26] <- NA # Individual 768A_320_17
    consdata[226,25:26] <- NA # NHNH419_317_94
    consdata[157,27:31] <- NA # 443_300_2
    consdata[222, 26] <- NA # 416_317_10
    consdata[237, 24:26] <- NA # 433_310_93
    consdata[139,27:29] <- NA # NHGA762B_320_6
    consdata[157,27:31] <- NA # NHNC443_300_2
    consdata[51, 28:31] <- NA # NCNC478_304_16
  }
  
  # Note that consumption measurments from snails which were paired with field males were not recorded on the raw datasheet (since we cant assume that field males wouldnt have different feeding rates / patterns to lab snails)
  
  
  
  ### Convert to long format and restore metadata----------------------------------------------------------------------------------------------

    {
      consdata=pivot_longer(consdata, cols = c("time.1":"time.15","c.1":"c.15"), names_to = c(".value", "Measure"), names_sep = "\\.")
      consdata=subset(consdata,consdata$c!="NA")
      
      consdata$time=as.POSIXct(consdata$time,format="%m/%d/%Y")
      
      consdata$Month <- as.factor(consdata$Measure)
      
      levels(consdata$Month) <- list('09/22'="1", '10/22'="2",'11/22'="3",'01/23'="4", '02/23'="5",'03/23'="6","04/23"="7","05/23"="8","06/23"="9","07/23"="10","08/23"="11","09/23"="12","10/23"="13","11/23"="14","12/23"="15")
    }
    
    {consdata$Measure=as.numeric(consdata$Measure)
      consdata$Snail_ID=as.factor(consdata$Snail_ID)  
      consdata$Treat=as.factor(consdata$Treat)
      consdata$Pop=as.factor(consdata$Pop)
      consdata$Egg_tag=as.factor(consdata$Egg_tag)
      consdata$Alive=as.factor(consdata$Alive)
      consdata$Sex=as.factor(consdata$Sex)
    }
    
    {
      consdata$Measure=as.factor(consdata$Measure)
      consdata=subset(consdata,consdata$Alive=="Y") #remove dead individuals
      
    }
    
    ### Halve values for consumption after measurement 6 ------------------

  #(to account for the fact that snails were all paired past this point (and data from single snails were removed earlier)
    
  #table(consdata$Measure)
 
    
    consdata$Measure=as.numeric(consdata$Measure)
    consdata$c[consdata$Measure > 6] <- consdata$c[consdata$Measure > 6] / 2
    consdata$Measure=as.factor(consdata$Measure) #convert Measure back to factor

    {consdataGeorgia=consdata #Create version including Georgia for supplemental material
      consdata=subset(consdata,consdata$Pop!="GA")}
    
    {consdata$Pop=as.character(consdata$Pop) #Correct number of pop factors (Because Georgia was removed form this version)
      consdata$Pop=as.factor(consdata$Pop)}
    
  
  #If loaded correctly, consdata should contain 1154 rows, consdataGeorgia should contain 1361 rows
    
    #write_csv(consdata, "processed_data/consdata.csv",na = "")
    #write_csv(consdataGeorgia, "processed_data/Georgia/consdataGeorgia.csv",na = "")
    

## Load Weight data -------------------------------------------------------------------------------------

# Note: this script creates two versions of the dataset: one with only data from individuals included in the orginal subset of 48 individuals (described in the paper) excluding those which died, and the second including individuals which were added as extra replacement animals during the course of the experiment    
# Only "weightdata" is used in the paper, preliminary analysis determined results do not differ significantly between the two datasets, so for simplicity, data from the replacement snails was excluded
    
{
combined_RT1_weightdata <- read.csv("~/github/RT1-data-analysis/Rawdata/combined_RT1_weightdata.csv")
weightdata=combined_RT1_weightdata

{ weightdata$Snail_ID=as.factor(weightdata$Snail_ID)  
  weightdata$Treat=as.factor(weightdata$Treat)
  weightdata$Pop=as.factor(weightdata$Pop)
  weightdata$Egg_tag=as.factor(weightdata$Egg_tag)
  weightdata$Alive=as.factor(weightdata$Alive)
  weightdata$Sex=as.factor(weightdata$Sex)
  weightdata$Month=as.factor(weightdata$Month)
  weightdata$time=as.POSIXct(weightdata$Date,format="%Y-%m-%d")}


#February 2023 data- air weights for snails 418 and 739 appear to have been swapped by mistake in raw datasheet. 


# Fix swapped air weight values for indiviudals 739 and 418, 

weightdata[530,16] <- 0.1809

weightdata[72,16] <- 0.0839

# and modify tissue weights to reflect swapped data

weightdata[530,17] <- 0.1221


weightdata[72,17] <- 0.0519



    weightdata=subset(weightdata,weightdata$Alive=="Y") #remove dead individuals
    
    weightdata$Snail_ID=as.character(weightdata$Snail_ID)
    
    weightdata$Snail_ID=as.factor(weightdata$Snail_ID) #Restore correct number of individuals
    
    weightdata$Sex=as.character(weightdata$Sex)
    
    weightdata$Sex=as.factor(weightdata$Sex) #Restore correct number of Sexes (because some individuals that died were immature)

    
    # Step to make sure factor levels are specified correctly for Month.    
    weightdata <- weightdata %>%
      mutate(Month = as.Date(paste0("01-", Month), format = "%d-%b-%y")) %>%
      arrange(Month) 
  weightdata$Month <- as.Date(weightdata$Month, "%b-%y")    
    
    weightdata$Month <- factor(weightdata$Month, levels = sort(unique(weightdata$Month)))
    
  weightdatawithreplacements=weightdata # Retain a copy  which includes individuals that were used to replace dead individuals within the subset of snails used for weight data
  
  counts=table(weightdata$Snail_ID)
  weightdata <- weightdata[weightdata$Snail_ID %in% names(counts[counts >= 12]), ]
  
  # This step excludes any individuals which were not used at all of the 13 timepoints (i.e. excludes replacements, and the three individuals which were transferred between temperature treatments midway through the experiment)
  # However, it keeps three individuals where data was not collected at a single timepoint

  weightdata$Snail_ID=as.character(weightdata$Snail_ID)
  
  weightdata$Snail_ID=as.factor(weightdata$Snail_ID) #Restore correct number of individuals
  

    }
  
  # If loaded correctly, weightdatawithreplacements should contain 586 rows.
    # Weightdata should contain 515 rows
 
    #write_csv(weightdata, "processed_data/weightdata.csv", na = "")
    
    #write_csv(weightdatawithreplacements, "processed_data/weightdatawithreplacements.csv", na = "")
    

## Load Repro data----



{
RT_repro <- read_excel("raw_data/RT1_reprodata.xlsx",sheet = "RawReproOutput_ALR")

RT_repro$Num_embryos[is.na(RT_repro$Num_embryos)] <- 0 #Replace some NA values in Num_embryos with zeroes to avoid errors in counting total embryos later on



#Create looping script for each individual

individuallist <- unique(growthdatageorgia$Snail_ID)

columns=c("Measure","embryos","clutches","capsules","time","Snail_ID")


fulldata=data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(fulldata)=columns
RT_repro$Snail_ID=as.factor(RT_repro$Snail_ID)

for(i in individuallist[1:length(individuallist)]) {
  
  repro_ind <- data.frame(Measure=as.numeric())
  embryos=list()
  clutches=list()
  capsules=list()
  library(plyr)
  snail=i
  data_ind=subset(growthdata, growthdata$Snail_ID==snail)
  
  timelist <- data.frame(Measure=data_ind$Measure, time=data_ind$time)
  
  timelist=structure(list(Measure = structure(c(8L, 9L, 1L, 11L, 12L, 13L, 
10L, 15L, 16L, 2L, 3L, 4L, 5L, 6L, 7L, 14L, 17L), levels = c("1", 
"2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", 
"14", "15", "16", "17"), class = "factor"), time = structure(c(1678680000, 
1681099200, 1659412800, 1686024000, 1688356800, 1691985600, 1683518400, 
1696910400, 1699246800, 1661745600, 1664769600, 1665979200, 1669870800, 
1673413200, 1676264400, 1694664000, 1701666000), class = c("POSIXct", 
  "POSIXt"), tzone = "")), class = "data.frame", row.names = c(NA, 
-17L))
  
  measurelist=list(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)
  
  
  for(i in measurelist[1:length(measurelist)]) {
    
    timepoint=as.numeric(i)
    time1=subset(timelist,timelist$Measure==timepoint-1)
    time2=subset(timelist,timelist$Measure==timepoint)
    daterange <- interval(time1$time[1], time2$time[1])
    repro=subset(RT_repro, RT_repro$Snail_ID==snail)
    repro$Date_End=as.POSIXct(repro$Date_End,format="%m/%d/%Y")
    reprorange<-repro[which(repro$Date_End %within% daterange),]
    embryos[[timepoint]]=sum(reprorange$Num_embryos)
    clutches[[timepoint]]=length(unique(reprorange$Clutch_num))
    reprorange = reprorange[ -c(14:15) ]
    reprorange=unique(reprorange)
    capsules[[timepoint]]=sum(reprorange$Total_num_capsules)
  }
  
  embryos <- ldply (embryos, data.frame)
  colnames(embryos)[1] <- "embryos"
  
  clutches <- ldply (clutches, data.frame)
  colnames(clutches)[1] <- "clutches"
  
  capsules <- ldply (capsules, data.frame)
  colnames(capsules)[1] <- "capsules"
  
  
  
  embryos$Measure=rownames(embryos)
  clutches$Measure=rownames(clutches)
  capsules$Measure=rownames(capsules)
  embryos=merge(embryos,clutches,by="Measure", all = T)
  embryos=merge(embryos,capsules,by="Measure", all = T)
  embryos=merge(embryos,timelist,by="Measure", all = T)
  embryos$Snail_ID=snail
  

  
  
  fulldata=rbind(fulldata,embryos)
  
  
}

{fulldata=merge(fulldata,RT_metadata,by="Snail_ID", all = T)
fulldata=unique(fulldata)
fulldata=subset(fulldata,fulldata$Sex=="female")
fulldata=subset(fulldata,fulldata$Alive=="Y")
}

{
fulldata$Pop=as.factor(fulldata$Pop)
fulldata$Treat=as.factor(fulldata$Treat)
fulldata2=fulldata
fulldatageorgia=fulldata # Retain version of fulldata that includes Georgia, for use in creating supplemental figures
fulldata=subset(fulldata,fulldata$Pop!="GA") 
fulldata$Pop=as.character(fulldata$Pop)

fulldata$Pop=as.factor(fulldata$Pop)
}



#Determine age at first reproduction and size at first reproduction


{individuallist <- unique(RT_repro$Snail_ID)
firstrepro=list()
sizeatfirstrepro=list()
totalembryos <- data.frame(totalembryos=as.numeric(),totalclutches=as.numeric(),totalcapsules=as.numeric())
}

for(i in individuallist[1:length(individuallist)]) {
  repro=subset(RT_repro, RT_repro$Snail_ID==i)
  firstclutch=subset(repro,repro$Clutch_num==1) #selects the first clutch produced by said individual
  firstrepro[i]=((as.character(firstclutch$Date_End[1]))) #Selects the date of the first clutch for that individual and adds it to the list
  firstreproind=as.data.frame(firstrepro[i])
  firstreproind[,1]=as.POSIXct(firstreproind[,1],format="%Y-%m-%d")
  firstclutch=subset(growthdatageorgia,growthdatageorgia$Snail_ID==i) #selects all clutches produced by said individual from the full dataset
  closestfirstrepro= firstclutch$time[which.min(abs(firstclutch$time - firstreproind[1,]))] #selects the nearest size measurement to the timepoint at which first reproduction occured
  sizeatfirstrepro[i] <- firstclutch$m[which(firstclutch$time == closestfirstrepro)] #Adds size measurement to the list
}



individuallist <- unique(fulldatageorgia$Snail_ID) 
for(i in individuallist[1:length(individuallist)]) {
  #i="NCNH407_311_90"
  firstclutch=subset(fulldatageorgia,fulldatageorgia$Snail_ID==i) #selects all clutches produced by said individual from the full dataset
  totalembryos[i,1]=sum(firstclutch$embryos)
  totalembryos[i,2]=sum(firstclutch$clutches)
  totalembryos[i,3]=sum(firstclutch$capsules)
}

#firstrepro contains age at which the first clutch was laid
#sizeatfirstrepro contains the size at the nearest timepoint to the date of first reproduction

{
firstrepro <- ldply (firstrepro, data.frame)
colnames(firstrepro)[1] <- "Snail_ID"
sizeatfirstrepro <- ldply (sizeatfirstrepro, data.frame)
colnames(sizeatfirstrepro)[1] <- "Snail_ID"
reproendstats=merge(firstrepro,sizeatfirstrepro,by="Snail_ID", all = T)
colnames(reproendstats)[2] <- "firstrepro_date"
colnames(reproendstats)[3] <- "firstrepro_size"

totalembryos$Snail_ID=rownames(totalembryos)
reproendstats=merge(reproendstats,totalembryos,by="Snail_ID", all = T)
reproendstats=merge(reproendstats,RT_metadata,by="Snail_ID", all = F)

# Adda calculation of age at first reproduction (measured from the start of the experiment, NOT the hatching of individual snails)
reproendstats$firstrepro_age=as.numeric(difftime(reproendstats$firstrepro_date,"2022-08-02",units="days"))


}

#Add column indicating if the individual reproduced or not during the course of the experiment

reproendstats$reproduced="Y"
reproendstats$reproduced[is.na(reproendstats$firstrepro_size)] ="N"


{
  reproendstatsGeorgia=reproendstats #version with Georgia data
  reproendstats=subset(reproendstats,reproendstats$Pop!="GA")
  reproendstats$Pop=as.factor(reproendstats$Pop)
  reproendstats$Treat=as.factor(reproendstats$Treat)
  
}

# If loaded correctly, fulldata should contain 1071 rows, fulldatageorgia should contain 1377 rows
# reproendstats should contain 63 rows, reproendstatsGeorgia should contain 81 rows

#Note, from what I can see, no snails that reproduced died during the experiment

#write_csv(fulldata, "processed_data/fulldata.csv", na = "")

#write_csv(fulldatageorgia, "processed_data/Georgia/fulldatageorgia.csv", na = "")

#write_csv(reproendstats, "processed_data/reproendstats.csv", na = "")
#write_csv(reproendstatsGeorgia, "processed_data/Georgia/reproendstatsGeorgia.csv", na = "")


}


## Load Hatching success data----


{hatchingsuccess= read_excel("raw_data/RT1_reprodata.xlsx",sheet = "RawHatchingSuccess_ALR")
hatchingsuccess$Treatment=as.factor(hatchingsuccess$Treatment)
hatchingsuccess$Population=as.factor(hatchingsuccess$Population)


hatchingsuccess$Prop_Hatchling_success[hatchingsuccess$Prop_Hatchling_success>1] =1 #Remove values over 100% for hatching success because of miscounted embryos

hatchingsuccess$Prop_Hatchling_success=hatchingsuccess$Prop_Hatchling_success

hatchingsuccessGeorgia= hatchingsuccess 
hatchingsuccess=subset(hatchingsuccess,hatchingsuccess$Population!="GA")
hatchingsuccess$Population=as.character(hatchingsuccess$Population)
hatchingsuccess$Population=as.factor(hatchingsuccess$Population)
}
  
  #If loaded correctly, hatchingsuccess should contain 199 rows, hatchingsuccessGeorgia should contain 262 rows

    write_csv(hatchingsuccess, "processed_data/hatchingsuccess.csv", na = "")
    #write_csv(fulldata, "processed_data/Georgia/hatchingsuccessGeorgia.csv", na = "")
    
    
}


# Load and prepare field temperature data------------------------------------------ 
{
## Load Lab Temperature data-------------------------------------------------

{labtemps <- read_csv("temperature_data/Lab_temperature_data.csv")
labtemps$time <- as.POSIXct(labtemps$time, format = "%Y/%m/%d %H:%M")
labtemps$pop = as.factor(labtemps$pop)}

## Load field temperature data-----------------------------------------------

###NORTH CAROLINA ----------------------------------------------------------
{dataNC <-
  read_csv("temperature_data/NDBC_BFTN7_2019-2021_merged.csv")
dataNC$wtmp_degc <- na_if(dataNC$wtmp_degc, 999)
dataNC2 <- na.omit(dataNC)
dataNC2$day_avg <- paste(dataNC2$month, dataNC2$day, sep = "-")
dataNC2$day_avg <- factor(dataNC2$day_avg)
daily_mean_NC <- tapply(dataNC2$wtmp_degc, dataNC2$day_avg, mean)
daily_mean_NC <- data.frame(daily_mean_NC)
daily_mean_NC$date <- row.names(daily_mean_NC)
daily_mean_NC$temp = daily_mean_NC[, 1]

daily_mean_NC2 = daily_mean_NC # Duplicate data across both years
daily_mean_NC$year = c("2023")
daily_mean_NC2$year = c("2022")

field_NC = rbind(daily_mean_NC, daily_mean_NC2)
field_NC = field_NC[, 2:4]
field_NC <- field_NC %>% unite("date", c(1, 3), sep = "-")
field_NC$date <- as.POSIXct(field_NC$date, tz = "", "%m-%d-%Y")
field_NC$pop = c("NC")}

#Resulting file is called field_NC

###NEW HAMPSHIRE-----------------------------------------------------------------

#### Data from Great Bay buoy ---------------------------------------------------

{dataNH <-
  read_csv("temperature_data/grbgbwq2019-2021_merged.csv")
dataNH <- na.omit(dataNH)
dataNH$rdate <-
  strptime(as.character(dataNH$DateTimeStamp), "%m/%d/%Y%H:%M")

dataNH$md <-
  format(as.Date(dataNH$rdate, format = "%d/%m/%Y"), "%m")
dataNH$day <-
  format(as.Date(dataNH$rdate, format = "%d/%m/%Y"), "%d")
dataNH$moday <- factor(paste(dataNH$md, dataNH$day, sep = "-"))
daily_mean_NH <- data.frame(tapply(dataNH$Temp, dataNH$moday, mean))
daily_mean_NH$date <- row.names(daily_mean_NH)
daily_mean_NH$temp = daily_mean_NH[, 1]


# Duplicate data across both years
daily_mean_NH2 = daily_mean_NH 
daily_mean_NH$year = c("2023")
daily_mean_NH2$year = c("2022")
field_NH = rbind(daily_mean_NH, daily_mean_NH2)
field_NH = field_NH[, 2:4]
field_NH <- field_NH %>% unite("date", c(1, 3), sep = "-")
field_NH$date <- as.POSIXct(field_NH$date, tz = "", "%m-%d-%Y")
field_NH$pop = c("NH")}

#### Data from CMLN3 buoy (UNH lab) --------------------------------------------
{data_NH2 <-
  read_csv("temperature_data/CMLN3_2019-2021.csv")
  data_NH2$WTMP_C <- na_if(data_NH2$WTMP_C, 999)
  data_NH3 <- data_NH2 %>%
    na.omit() %>%
    mutate (date = paste(MM, DD, sep = "-")) %>%
    mutate (date = factor(date))
  daily_mean_NH3 <- tapply(data_NH3$WTMP_C, data_NH3$date, mean)
  daily_mean_NH3 <- data.frame(daily_mean_NH3)
  daily_mean_NH3$date <- row.names(daily_mean_NH3)
  daily_mean_NH4 = daily_mean_NH3
  daily_mean_NH3$year = c("2023")
  daily_mean_NH4$year = c("2022")
  field_NH2 = rbind(daily_mean_NH3, daily_mean_NH4)
  field_NH2$temp = field_NH2[, 1]
  field_NH2 = field_NH2[, 2:4]
  field_NH2 <- field_NH2 %>% unite("date", c(1, 2), sep = "-")
  field_NH2$date <- as.POSIXct(field_NH2$date, tz = "", "%m-%d-%Y")
  field_NH2$pop = c("NH")}

#### Merge Great Bay and UNH data-----------------------------------------------
{
field_NH_merged <- merge(field_NH, field_NH2, by = "date", all.y = TRUE, suffixes = c(".greatbay", ".CMLN3"))

# Replace values from great bay with those from offshore sensor where values are missing (because the Great Bay sensor is pulled out in the winter when freezing becomes a risk)


field_NH_merged <- within(field_NH_merged, {
  temp <- ifelse(is.na(temp.greatbay), temp.CMLN3, temp.greatbay)
  pop <- ifelse(is.na(pop.greatbay), pop.CMLN3, pop.greatbay)
})

field_NH_merged <- field_NH_merged[, c("date", "temp", "pop")]
}

#field_NH-merged contains combined NH data


### Combine Field data into single file ------

{fieldtemps = rbind(field_NH_merged, field_NC)
fieldtemps$time = fieldtemps$date

fieldtemps$pop = as.factor(fieldtemps$pop)
}


#write_csv(fieldtemps, "processed_data/fieldtemps.csv")

}

## Descriptive statistics on mortality


table(RT_metadata$Alive,RT_metadata$Pop,RT_metadata$Treat) 

metadata_subset=subset(RT_metadata, RT_metadata$Was_reloaded!="Y")

table(metadata_subset$Alive,metadata_subset$Pop,metadata_subset$Treat) 

table(metadata_subset$Death_date,metadata_subset$Pop,metadata_subset$Treat) 

table(RT_metadata$Death_date,RT_metadata$Pop,RT_metadata$Treat) 


