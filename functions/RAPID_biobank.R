################################################################################
## The RAPID Biobank Sample Database
## May 2021
## Joseph.Shaw@gosh.nhs.uk
## These scripts allow the Excel spreadsheet of the RAPID 
## biobank to be read into R as a tibble.
################################################################################

#############################################################
# 1 - Loading the biobank database
#############################################################

# Load necessary packages
library(tidyverse)
library(readxl)
library(lubridate)

biobank_filepath <- "I:/Genetics/RAPID/Biobank Library/AA current biobank/"

# Use filename pattern recognition with a wildcard for the date 
# to select the current biobank.
biobank_current <- list.files(path= biobank_filepath, 
                              pattern="^RAPID_Biobank_database_.*.xlsx$")

# Loads the most recent version of the RAPID biobank.
RAPID_biobank <- read_excel(paste0(biobank_filepath,biobank_current),
                            sheet = "database",
                            # Have to specify the data types for each column 
                            # as not all cells have data
                            # and read_excel cannot guess correctly.
                            col_types = c(
                              #r_number
                              "numeric",
                              #study_id
                              "text",
                              #site
                              "text",
                              #date_of_blood_sample
                              "date",
                              #time_of_blood sample
                              "date",
                              #maternal_DOB
                              "date",
                              #expected_delivery_date
                              "date",
                              #gestation_weeks
                              "numeric",
                              #gestation_days
                              "numeric",
                              #date_of_first_spin
                              "date",
                              #time_first_spin
                              "date",
                              #first_interval
                              "text",
                              #date_of_2nd_spin
                              "date",
                              #time_of_2nd_spin
                              "date",
                              #second_ Interval
                              "text",
                              #vacutainer
                              "text",
                              #original_plasma_vol
                              "numeric",
                              #current_volume
                              "numeric",
                              #removed_for_plasma
                              "text",
                              #date_removed
                              "date",
                              #blood_Pellets
                              "numeric",
                              #maternal_episode_number
                              "text",
                              #maternal_episode_removed_for
                              "text",
                              #maternal_episode_date_removed
                              "date",
                              #paternal_episode_number
                              "text",
                              #paternal_episode_removed_for
                              "text",
                              #paternal_episode_date_removed
                              "date",
                              #fetal_episode_number
                              "text",
                              #fetal_episode_removed for_excluding_PAGE
                              "text",
                              #fetal_episode_date_removed
                              "date",
                              #indication
                              "text",
                              #risk
                              "text",
                              #confirmed_diagnosis
                              "text",
                              #reference
                              "text",
                              #mutation_genetic_info_fetus
                              "text",
                              #maternal_mutation
                              "text",
                              #paternal_mutation
                              "text",
                              #confirmed_fetal_sex
                              "text",
                              #ethnicity
                              "text",
                              #maternal_weight_kg
                              "numeric",
                              #maternal_height_cm
                              "numeric",
                              #diabetic
                              "text",
                              #smoking
                              "text",
                              #before_after_invasive
                              "text",
                              #partner_date_of_birth
                              "date",
                              #additional_comments
                              "text",
                              #report_acquired
                              "text",
                              #Recruited_to_EACH
                              "text",
                              #study_name
                              "skip",
                              #Recruited_to_PAGE
                              "skip",
                              #Recruited_to_BOOSTB4
                              "skip",
                              #Recruited_to_Qiagen
                              "text",
                              #EU_Consent
                              "text",
                              #Long_Term_follow_up
                              "text",
                              #DNA_storage_at_Sanger
                              "text",
                              #Checked_By_GOSH
                              "text",
                              #Fetal_Sex_Result_by_Lab_qPCR_dPCR
                              "text",
                              #Fetal_sex_double_checked_locally
                              "text",
                              #action_plasma
                              "text",
                              #tubes_plasma_original
                              "numeric",
                              #tubes_plasma_current
                              "numeric",
                              #Plasma_location_tray_number
                              "numeric",
                              #Plasma_location_Y
                              "numeric",
                              #Plasma_location_X
                              "numeric",
                              #action_blood
                              "text",
                              #tubes_blood_original
                              "numeric",
                              #tubes_blood_current
                              "numeric",
                              #Blood_location_tray_number
                              "numeric",
                              #Blood_location_box
                              "numeric",
                              #Reorganisation_notes
                              "text")) %>%
  # Convert gestation into numeric form
  mutate(gestation_days_as_weeks = gestation_days/7,
         Gestation_total_weeks = gestation_weeks + gestation_days_as_weeks,
         
         # Print the gestation in the way that fetal medicine units write it
         # I.e. "12+5" as 12 weeks and 5 days
         gestation_character = paste0(as.character(gestation_weeks), 
                                      "+", as.character(gestation_days)),
         Partner_sample_available = ifelse(paste0("Partner",study_id) %in% 
                                             study_id, "Yes", "No"),
         year = year(date_of_blood_sample),
         
         # Add sample type to allow easy exclusion of partner 
         # samples when searching.
         sample_type = ifelse(startsWith(study_id, "Partner"), 
                              "Partner", "Pregnant woman"))

#############################################################
# 2 - Searching the biobank database
#############################################################

# Finding relevant samples in the RAPID biobank Excel spreadsheet is tricky
# due to the use of free-type fields for sample indication/diagnosis etc and 
# lack of standardisation.

# There are several free-type columns that could have relevant information: 
# indication, mutation_genetic_info_fetus, risk, 
# confirmed_diagnosis, maternal_mutation, 
# paternal_mutation, additional_comments.

# The search_biobank function allows you to input various strings 
# that you want to search for and then uses grep to find matches in 
# the relevant columns, before filtering and and returning the results table. 
# Search terms can be input as a character vector.
# grep is set to ignore.case = TRUE, so search terms can be in any case.

# Example: search_biobank(c("sma ", "spinal", "atrophy", "smn1"))
# Example: search_biobank(c("CF", "cftr", "delta508", "cystic fibrosis"))

search_biobank <- function(search_terms) {
  search_hits <- unique(grep(paste(search_terms,collapse="|"), 
                 # Columns to look in
                 c(RAPID_biobank$indication, 
                   RAPID_biobank$mutation_genetic_info_fetus,
                   RAPID_biobank$risk, 
                   RAPID_biobank$confirmed_diagnosis,
                   RAPID_biobank$additional_comments, 
                   RAPID_biobank$maternal_mutation,
                   RAPID_biobank$paternal_mutation), 
                 value=TRUE, ignore.case = TRUE))
  
  search_results <- RAPID_biobank %>%
    filter(indication %in% search_hits 
           | mutation_genetic_info_fetus %in% search_hits
           | risk %in% search_hits
           | confirmed_diagnosis %in% search_hits
           | additional_comments %in% search_hits
           | maternal_mutation %in% search_hits
           | paternal_mutation %in% search_hits) %>%
    # Only show samples which were not discarded in the reorganisation process
    filter(r_number > 10444) %>%
    select(r_number, study_id, site, date_of_blood_sample, 
           maternal_DOB, indication, 
           mutation_genetic_info_fetus,confirmed_diagnosis, 
           additional_comments, maternal_mutation,
           paternal_mutation, tubes_plasma_current)
  
  return(search_results)
}

#############################################################
# 3 - Checking sample locations
#############################################################

# Plot the locations of the plasma tubes

ggplot(RAPID_biobank, 
       aes(x = Plasma_location_X, y = Plasma_location_Y))+
  geom_point(size = 3, alpha = 0.4)+
  facet_wrap(~Plasma_location_tray_number)+
  theme_bw()+
  scale_y_continuous(minor_breaks = seq(0.5 , 40.5, 1), 
                     breaks = seq(0, 40, 10))+
  scale_x_continuous(minor_breaks = seq(0.5 , 10.5, 1), 
                     breaks = seq(0, 10, 5))+
 theme(panel.grid.major = element_line(colour="white", size=0.5))

#############################################################
# 4 - Checking sample tube counts
#############################################################

# Blood and plasma tube counts are manually updated whenever tubes
# are removed, so that we have an accurate picture of how many tubes
# we have for each sample.

samples_removed <- read_excel(paste0(biobank_filepath,biobank_current),
                            sheet = "samples_removed") %>%
  # Freezer reorganisation was completed on 09/09/2020
  filter(removal_date > "2020-09-09") %>%
  # Extract R number from "epic_name" and convert to numeric
  mutate(r_number = as.numeric(
    readr::parse_number(epic_name, trim_ws = TRUE)))

plasma_removed <- samples_removed %>%
  filter(sample_type == "Plasma") %>%
  select(r_number, tubes_removed) %>%
  # Summarise the number of plasma tubes removed by r number
  group_by(r_number) %>% 
  summarise(plasma_tubes_removed = sum(tubes_removed),
            .groups="drop")

blood_removed <- samples_removed %>%
  filter(sample_type == "Blood") %>%
  select(r_number, tubes_removed) %>%
  group_by(r_number) %>% 
  summarise(blood_tubes_removed = sum(tubes_removed),
            .groups="drop")

# Check plasma tube values have been correctly updated

plasma_check_table <- plasma_removed %>% left_join (
  RAPID_biobank %>%
    select(r_number, tubes_plasma_original, tubes_plasma_current,
           action_plasma, Plasma_location_tray_number,
           Plasma_location_Y, Plasma_location_X,
           Reorganisation_notes),
  by = "r_number") %>%
  filter(action_plasma != "Sample Not Found") %>%
  mutate(plasma_check = ifelse(tubes_plasma_current == 
                                 (tubes_plasma_original - plasma_tubes_removed),
                               "TRUE",
                               "FALSE")) %>%
  select(r_number, tubes_plasma_original,
         plasma_tubes_removed, tubes_plasma_current, plasma_check,
         Plasma_location_tray_number,
         Plasma_location_Y, Plasma_location_X,
         Reorganisation_notes) %>%
  arrange(Plasma_location_tray_number)

# Check blood tube values have been correctly updated

blood_check_table <- blood_removed %>% left_join (
  RAPID_biobank %>%
    select(r_number, tubes_blood_original, tubes_blood_current,
           action_blood),
  by = "r_number") %>%
  filter(action_blood != "Sample Not Found") %>%
  mutate(blood_check = ifelse(tubes_blood_current == 
                                 (tubes_blood_original - blood_tubes_removed),
                               "TRUE",
                               "FALSE"))

#############################################################
# 5 - Tube tray locations
#############################################################

# Each tray is a a grid with 10 spaces across and 40 long.
# The coordinate order must go tray, y, x (for readability)

#################################
# Get all tray locations
#################################

# Create a dataframe with coordinates for every tray

tray_coordinates <- data.frame(tray = c(),
                           x_coordinate = c(),
                           y_coordinate = c())

# Create a data-frame of every position in the biobank trays
# Do for trays up to 54

for (i in 1:54) {
  for (j in 1:10) {
    for (k in 1:40) {
      temp_output <- data.frame(tray = c(i),
                                x_coordinate = c(j),
                                y_coordinate = c(k))
      tray_coordinates <- rbind(tray_coordinates, temp_output)
      rm(temp_output)
    }
  }
}

#################################
# Get all tube locations
#################################

# New strategy - make individual dataframes of long data for each r_number,
# 1 row per tube. Then rbind all those tables together.

make_coord_long <- function(single_sample) {
  
  new_data <- single_sample
  
  for(i in seq(2,(single_sample$tubes_plasma_current))) {
    
    new_row <- c(single_sample$r_number, 
                 single_sample$tubes_plasma_current, 
                 i, 
                 single_sample$tray, 
                 single_sample$position_y, 
                 single_sample$position_x + (i-1))
    
    new_data <- rbind(new_data, new_row)
    
    rm(new_row)
  }
  return(new_data)
}

# All relevant plasma samples
RAPID_plasma <- RAPID_biobank %>%
  #Remove samples which got separated during the reorganisation
  filter(is.na(Reorganisation_notes)) %>%
  select(r_number, tubes_plasma_current, Plasma_location_tray_number,
         Plasma_location_Y, Plasma_location_X) %>%
  # Remove rows with NAs
  drop_na() %>%
  # Remove rows with no plasma left
  filter(tubes_plasma_current != 0) %>%
  filter(Plasma_location_tray_number < 55) %>%
  mutate(tube = 1) %>%
  dplyr::rename(
    tray = Plasma_location_tray_number,
    position_y = Plasma_location_Y,
    position_x = Plasma_location_X) %>%
  select(r_number, tubes_plasma_current, tube,
         tray, position_y, position_x)

# For every row of RAPID_plasma, perform make_coord_long and then rbind together

total_coordinates <- data.frame()

for (i in 1:nrow(RAPID_plasma)) {
  
  next_sample <- make_coord_long(RAPID_plasma[i,])
  
  total_coordinates <- rbind(total_coordinates, next_sample)
  
  rm(next_sample)
  
}

# Now we need to apply the dimensions of the trays

total_test <- total_coordinates %>%
  # Apply the x dimension
  mutate(position_x_new = ifelse(position_x > 10,
                                 position_x -10,
                                 position_x),
         
         position_y_new = ifelse(position_x > 10,
                                 position_y+1,
                                 position_y),

        # Apply the y dimension
        position_y_newnew = ifelse(position_y_new > 40,
                                   position_y_new -40,
                                   position_y_new),
        
        tray_new = ifelse(position_y_new > 40,
                          tray+1,
                          tray),
        
        # Concatenate the coordinates in a single field
        coordinate_string = paste(tray_new, position_y_newnew, position_x_new,
                                  sep = "_"))


## Plot a graph of this

ggplot(total_test %>%
         filter(tray_new %in% c(14)), 
       aes(x = position_x_new, y = position_y_newnew, 
           label =  r_number))+
  geom_point(size = 5, pch = 21)+
  geom_text(size = 1) +
  facet_wrap(~tray_new)+
  theme_bw()+
  theme(panel.grid.major = element_line(colour="white", size=0.5)) +
  scale_y_continuous(minor_breaks = seq(0.5 , 40.5, 1), 
                     breaks = seq(0, 40, 10))+
  scale_x_continuous(minor_breaks = seq(0.5 , 10.5, 1), 
                     breaks = seq(0, 10, 5))

# Check for tubes in the same coordinates
total_test[duplicated(total_test$coordinate_string),]
