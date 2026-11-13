### 0. Importing and cleaning data
### Authors: Jorge Mora - Bianca Campagnari
### The objective of this code is to import data and to clean data.
### The data is a results of a metabarcoding process for Osmia lignaria brood cells collected from three sample sites in Northem CA

### Setting working directory
setwd("/Users/jorgemora/Documents/JorgeMora/PhD/Classes/Fall_2026/EEOB_215/FinalProject_EEOB215/Exported_Data") #Working directory for JM

# Creating new folder "output" where we will save the results of the code
dir.create("output/01_cleandata",  recursive = TRUE, showWarnings = FALSE) 
dir.create("output/02_sensors", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures",   recursive = TRUE, showWarnings = FALSE)


## Loading Libraries
library(tidyverse)

## Importing the data

# Metadata
sample_metadata <- read_tsv("sample-metadata.tsv") |> 
  rename("SampleID" = 1) |> 
  filter(!str_detect(Site, "Thailand")) # Here we are removing etadata from samples that are not par od the analysis

#Saving the data
write_csv(sample_metadata, 'output/01_cleandata/metadata.csv')

# Plants
plants_filtered <- read_tsv("plants_filtered_table.tsv", skip = 1) |>
  rename("OTU_ID" = 1) |> 
  select(!matches("^BLANK|^EMPTY"))
  
plants_taxonomy <- read_tsv("plants_taxonomy.tsv") |> 
  rename("OTU_ID" = 1) |> 
  as.data.frame("plants_filtered & plants_taxonomy")

plants <- plants_filtered |> 
  left_join(plants_taxonomy, by = "OTU_ID") |> 
  filter(!str_detect(Taxon, "Unassigned|Theobroma")) |> # Here Im excluding theobroma because the sequencing is telling us that the bees foraged on Cacao in California. After a BLAST search    
  select(!matches("^TS")) |> #Removind samples that are not part of the analysis
  as.data.frame()

plants_per_site <- plants |> 
  pivot_longer( #We will use this to turn our wide data in long data
  cols = -c(OTU_ID, Taxon, Confidence), # We will keep these three columns
  names_to = "Sample_ID",      # We will take the colnames and convert them in variables within the column "Sample_ID"
  values_to = "Count") |>       # We will take the read of each ASV and we will assing it as the values, in this case the counts for each sample_ID
  mutate(
    Site = case_when( # We will create a Site column base on a conditional
      str_detect(Sample_ID, "JR") ~ "Jims Ranch",
      str_detect(Sample_ID, "SP") ~ "Steve Porch",
      str_detect(Sample_ID, "SG") ~ "Steve Garage",
      str_detect(Sample_ID, "VP") ~ "Vernal Pools"),
    Round = str_extract(Sample_ID, "R\\d+") 
  ) |>
  summarise(
    Total_Count = sum(Count, na.rm = TRUE),
    .by = c(Site, Round, Taxon) # .by is equivalent to group_by(), but they are not the same
  )

#### With plants we have a situacion, bunchs of read are identified as cacao

# top_ids <- plants |>
#   filter(str_detect(Taxon, "Theobroma")) |>
#   mutate(total = rowSums(across(where(is.numeric) & !Confidence))) |>
#   slice_max(total, n = 10) |>
#   pull(OTU_ID)

# Saving the data for plants
write_csv(plants, 'output/01_cleandata/plants.csv')
write_csv(plants_per_site, 'output/01_cleandata/plants_per_site.csv')


# For bacteria (16s)
bacteria_filtered <- read_tsv("16s_filtered_table.tsv", skip = 1) |>
  rename("OTU_ID" = 1) |> 
  select(!matches("^BLANK|^EMPTY"))

bacteria_taxonomy <- read_tsv("16s_taxonomy.tsv") |> 
  rename("OTU_ID" = 1) |> 
  as.data.frame("bacteria_filtered & bacteria_taxonomy")

bacteria <- bacteria_filtered |> 
  left_join(bacteria_taxonomy, by = "OTU_ID") |> 
  filter(!str_detect(Taxon, "Chloroplast|Mitochondria|Unassigned")) |> #Here we are filtering the taxon column to remove mitochondria, chloroplast and unassigned taxons
  select(!matches("^TS")) |> #Removind samples that are not part of the analysis
  as.data.frame()

bacteria_per_site <- bacteria |> 
  pivot_longer( #We will use this to turn our wide data in long data
    cols = -c(OTU_ID, Taxon, Confidence), # We will keep these three columns
    names_to = "Sample_ID",      # We will take the colnames and convert them in variables within the column "Sample_ID"
    values_to = "Count") |>       # We will take the read of each ASV and we will assing it as the values, in this case the counts for each sample_ID
  mutate(
    Site = case_when( # We will create a Site column base on a conditional
      str_detect(Sample_ID, "JR") ~ "Jims Ranch",
      str_detect(Sample_ID, "SP") ~ "Steve Porch",
      str_detect(Sample_ID, "SG") ~ "Steve Garage",
      str_detect(Sample_ID, "VP") ~ "Vernal Pools"),
    Round = str_extract(Sample_ID, "R\\d+") 
  ) |>
  summarise(
    Total_Count = sum(Count, na.rm = TRUE),
    .by = c(Site, Round, Taxon)
  )

# Saving the data for bacteria
write_csv(bacteria, 'output/01_cleandata/bacteria.csv')
write_csv(bacteria_per_site, 'output/01_cleandata/bacteria_per_site.csv')

#Fungal
fungi_filtered <- read_tsv("fungi_filtered_table.tsv", skip = 1) |>
  rename("OTU_ID" = 1) |> 
  select(!matches("^BLANK|^EMPTY"))

fungi_taxonomy <- read_tsv("fungi_taxonomy.tsv") |> 
  rename("OTU_ID" = 1) |> 
  as.data.frame("fungi_filtered & fungi_taxonomy")

fungi <- fungi_filtered |> 
  left_join(fungi_taxonomy, by = "OTU_ID") |> 
  filter(!str_detect(Taxon, "Unassigned")) |> 
  select(!matches("^TS")) |> #Removind samples that are not part of the analysis
  as.data.frame()

fungi_per_site <- fungi |> 
  pivot_longer( #We will use this to turn our wide data in long data
    cols = -c(OTU_ID, Taxon, Confidence), # We will keep these three columns
    names_to = "Sample_ID",      # We will take the colnames and convert them in variables within the column "Sample_ID"
    values_to = "Count") |>       # We will take the read of each ASV and we will assing it as the values, in this case the counts for each sample_ID
  mutate(
    Site = case_when( # We will create a Site column base on a conditional
      str_detect(Sample_ID, "JR") ~ "Jims Ranch",
      str_detect(Sample_ID, "SP") ~ "Steve Porch",
      str_detect(Sample_ID, "SG") ~ "Steve Garage",
      str_detect(Sample_ID, "VP") ~ "Vernal Pools"),
    Round = str_extract(Sample_ID, "R\\d+") 
  ) |>
  summarise(
    Total_Count = sum(Count, na.rm = TRUE),
    .by = c(Site, Round, Taxon)
  )

# Saving the data for fungi
write_csv(fungi, 'output/01_cleandata/fungi.csv')
write_csv(fungi_per_site, 'output/01_cleandata/fungi_per_site.csv')


#### Importing the data from the sensors data
### We have a folder, with folders per date, each one containing a .csv file with the sensors data (e.g. Atmostpheric pressure, Tempretarure in and out of the nest, Relative
### humidity in and out of the nest) 
### So, we need to create a code wich turn all those folder in a single dataframe with the data


## Helios nest -Vernal Pools

?list.files

files_helios <- list.files(    # list.files() looks inside a folder and returns the names of what it finds, this function just get the addresses of the files
  path = "osmiaCAM_envdata/envLogger_helios",  # Directory from where we will get the files
  pattern = "\\.csv$", # keep only names that end in .csv
  recursive = TRUE, # This is a logical arguments, which is use to tell the fuction to check the subfolder 
  full.names = TRUE # a logical value. If TRUE, the directory path is prepended to the file names to give a relative file path. If FALSE, the file names (rather than paths) are returned.
)

print(files_helios) # So we create a list of path for the data that we want to read

#Now, we have to use a function that read that list of file and turn it into a df

VernalPools <- read_csv(files_helios, 
                   col_names = c("datetime", "pressure", "temp_out", "rh_out", "temp_in", "rh_in"),
                   id = 'file_origin_path' ) |> 
  mutate(Site = "Vernal Pools",
         rh_out_flag = rh_out > 100,  # Some values for the relative humidity were over 100% which is physically imposible
         rh_out = pmin(rh_out, 100) )

# Saving the sensor data for Vernal Pools
write_csv(VernalPools, 'output/02_sensors/Vernal_Pools.csv')

## Jim's Ranch
files_JimsRanch <- list.files(
  path = "osmiaCAM_envdata/envLogger_JimsRanch",
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

JimsRanch <- read_csv(files_JimsRanch,
                      col_names = c("datetime", "pressure", "temp_out", "rh_out", "temp_in", "rh_in"),
                      id = "file_origin_path") |> 
  mutate(Site = "Jim Ranch",
         rh_out_flag = rh_out > 100,  # Some values for the relative humidity were over 100% which is physically imposible
         rh_out = pmin(rh_out, 100))

head(JimsRanch)

# Saving the sensor data for Vernal Pools
write_csv(JimsRanch, 'output/02_sensors/Jim_Ranch.csv')

## Osmia 3 - Steve Porche
files_osmia3 <- list.files(
  path = "osmiaCAM_envdata/envLogger_osmia3",
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

StevePorche <- read_csv(files_osmia3,
                      col_names = c("datetime", "pressure", "temp_out", "rh_out", "temp_in", "rh_in"),
                      id = "file_origin_path") |> 
  mutate(Site = "Steve Porche",
         rh_out_flag = rh_out > 100,  # Some values for the relative humidity were over 100% which is physically imposible
         rh_out = pmin(rh_out, 100) )

# Saving the sensor data for Steve Porche
write_csv(StevePorche, 'output/02_sensors/Steve_Porche.csv')

## Osmia 4 - Steve Garage
files_osmia4 <- list.files(
  path = "osmiaCAM_envdata/envLogger_osmia4",
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

SteveGarage <- read_csv(files_osmia4,
                      col_names = c("datetime", "pressure", "temp_out", "rh_out", "temp_in", "rh_in"),
                      id = "file_origin_path") |> 
  mutate(Site = "Steve Garage",
         rh_out_flag = rh_out > 100,  # Some values for the relative humidity were over 100% which is physically imposible
         rh_out = pmin(rh_out, 100)   # This happened because of condensation on the sensor. So we will keep the max RH at 100% but the other column will tell us which columne was an outlier
  )
# Saving the sensor data for Vernal Pools
write_csv(SteveGarage, 'output/02_sensors/Steve_Garage.csv')


#### Binding all the sensor data together
all_sites <- bind_rows(
  "Steve Garage" = SteveGarage,
  "Steve Porche"  = StevePorche,
  "Vernal Pools" = VernalPools,
  "Jims Ranch"   = JimsRanch,
  .id = "Site"
)

str(all_sites)

# Saving the sensor data for all the sites together
write_csv(all_sites, 'output/02_sensors/All_Sites.csv')
