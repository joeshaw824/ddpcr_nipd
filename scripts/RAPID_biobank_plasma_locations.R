################################################################################
## RAPID Biobank Functions and Locations
## August 2021
## Joseph.Shaw@gosh.nhs.uk
## This script is for searching for clinical information and plasma
## samples in the RAPID Biobank database.
################################################################################

source("functions/RAPID_biobank.R")

#############################################################
# Searching the biobank database
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
# Checking sample tube counts
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
# Assume 58 trays 

for (i in 1:58) {
  for (j in 1:10) {
    for (k in 1:40) {
      temp_output <- data.frame(tray = c(i),
                                position_y  = c(k),
                                position_x = c(j))
      tray_coordinates <- rbind(tray_coordinates, temp_output)
      rm(temp_output)
    }
  }
}

#################################
# Get all tube locations
#################################

# All relevant plasma samples
RAPID_plasma <- RAPID_biobank %>%
  select(r_number, tubes_plasma_current, Plasma_location_tray_number,
         Plasma_location_Y, Plasma_location_X) %>%
  # Remove rows with NAs
  drop_na() %>%
  # Remove rows with no plasma left
  filter(tubes_plasma_current != 0) %>%
  mutate(tube = 1) %>%
  dplyr::rename(
    tray = Plasma_location_tray_number,
    position_y = Plasma_location_Y,
    position_x = Plasma_location_X) %>%
  select(r_number, tubes_plasma_current, tube,
         tray, position_y, position_x)

# New strategy - make individual dataframes of long data for each r_number,
# 1 row per tube. Then rbind all those tables together.

make_coord_long <- function(single_sample) {
  
  new_data <- single_sample
  
  for(i in seq(2,(single_sample$tubes_plasma_current))) {
    
    new_row <- c(
      # R number
      single_sample$r_number, 
      # Tubes of plasma
      single_sample$tubes_plasma_current,
      # Tube number
      i, 
      # Tray
      single_sample$tray, 
      # Position y
      single_sample$position_y, 
      # Position x
      single_sample$position_x + (i-1))
    
    new_data <- rbind(new_data, new_row)
    
    rm(new_row)
  }
  return(new_data)
}

# Identify the samples which only have one tube of plasma

single_tube_cases <- RAPID_plasma %>%
  filter(tubes_plasma_current == 1)

mutliple_tube_cases <- RAPID_plasma %>%
  filter(tubes_plasma_current > 1)

# For every row of RAPID_plasma with more than 1 tube of 
# plasma, perform make_coord_long and then rbind together

multiple_tube_coordinates <- data.frame()

for (i in 1:nrow(mutliple_tube_cases)) {
  
  next_sample <- make_coord_long(mutliple_tube_cases[i,])
  
  multiple_tube_coordinates <- rbind(multiple_tube_coordinates, next_sample)
  
  rm(next_sample)
  
}

# Bind to single tube cases
single_and_multiple_locations <- rbind(multiple_tube_coordinates, 
                                       single_tube_cases) %>%
  select(r_number, tray, position_y, position_x)

# Add on the extra tubes from the reorganisation that were stored in
# separate places

extra_tube_locations <- read.csv("data/extra_tube_locations.csv")

all_tube_locations <- rbind(single_and_multiple_locations, 
                            extra_tube_locations %>%
                              select(r_number, tray, 
                                     position_y, position_x))

#################################
# Apply tray dimensions
#################################

all_tubes_dimensions <- all_tube_locations %>%
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
                           tray)) %>%
  select(r_number, tray_new, 
         position_y_newnew, position_x_new) %>%
  dplyr::rename(
    tray = tray_new,
    position_y = position_y_newnew,
    position_x = position_x_new) %>%
  arrange(tray, position_y, position_x) %>%
  # Concatenate the coordinates in a single field
  mutate(coordinate_string = paste(tray, position_y, position_x,
                                   sep = "_"),
         status = "occupied") %>%
  select(tray, position_y, position_x, coordinate_string,
         status, r_number)


#################################
# Combine with all possible locations
#################################

empty_coordinates <- tray_coordinates %>%
  mutate(coordinate_string = paste(tray, position_y, position_x,
                                   sep = "_")) %>%
  filter(!coordinate_string %in% all_tubes_dimensions$coordinate_string) %>%
  mutate(status = "empty",
         r_number = NA)

biobank_plasma_locations <- rbind(all_tubes_dimensions,
                                  empty_coordinates) %>%
  arrange(tray, position_y, position_x)

write.csv(biobank_plasma_locations, 
          "analysis_outputs/biobank_plasma_locations.csv",
          row.names = FALSE)

#################################
# Checking against actual locations
#################################

live_locations <- read_excel(paste0(biobank_filepath, 
                                    "biobank_plasma_locations.xlsx"))

#################################
# Plot graphs
#################################

plasma_locations <- ggplot(live_locations %>%
                             filter(tray %in% (1) & 
                                      r_number != 0), 
                           aes(x = position_x, y = position_y, 
                               label = r_number))+
  geom_point(size = 4, pch = 1)+
  geom_text(size = 1) +
  facet_wrap(~tray)+
  theme_bw()+
  theme(panel.grid.major = element_line(colour="white", size=0.5)) +
  scale_y_continuous(minor_breaks = seq(0.5 , 40.5, 1), 
                     breaks = seq(0, 40, 5))+
  scale_x_continuous(minor_breaks = seq(0.5 , 10.5, 1), 
                     breaks = seq(0, 10, 5)) +
  theme(legend.position = "none") +
  labs(x = "", y = "", title = "RAPID Biobank Plasma Locations")



ggsave(filename = plasma_locations,
       device = "jpeg",
       plot = plasma_locations,
       path = "plots/")


live_locations_check <- live_locations %>%
  mutate(coordinate_string = paste(tray, position_y, position_x,
                                   sep = "_"))

live_locations_check[duplicated(live_locations_check$coordinate_string),]


ggplot(live_locations %>%
         filter(tray %in% c(50:58)), 
       aes(x = status, y = ))+
  geom_bar() +
  facet_wrap(~tray)+
  theme_bw()+
  theme(panel.grid.major = element_line(colour="white", size=0.5))

