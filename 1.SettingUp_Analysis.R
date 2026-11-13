### 1. Setiing up the analysis
### Authors: Jorge Mora - Bianca Campganari
### After filtering and cleaning the data, the objective of this code is to set up the analysis

#### Setting working directory
setwd("/Users/jorgemora/Documents/JorgeMora/PhD/Classes/Fall_2026/EEOB_215/FinalProject_EEOB215/Exported_Data/output") #Working directory for JM

#### Loading libraries
library(tidyverse); library(vegan); library(glmmTMB); library(DHARMa)

#### Importing the data 

sensors <- read.csv("02_sensors/All_Sites.csv")

sensors <- sensors |> 
  mutate(date = as_date(datetime),
         hour = hour(datetime))

sensors_summary <- sensors |> # A situation here is the we are assuming that all the sensors have the same days of data for example
  group_by(Site) |> 
  summarise(
    mean_temp_in  = mean(temp_in, na.rm = TRUE),
    mean_temp_out = mean(temp_out, na.rm = TRUE),
    mean_rh_in    = mean(rh_in, na.rm = TRUE),
    mean_rh_out   = mean(rh_out, na.rm = TRUE),
    max_temp_in = max(temp_in, na.rm = TRUE),
    min_temp_in = min(temp_in, na.rm = TRUE),
    max_rh_in = max(rh_in, na.rm = TRUE),
    min_rh_in = min(rh_in, na.rm = TRUE)
  )

?na.rm=TRUE
  