################################################################################
## The RAPID Biobank Sample Database
## May 2021
## Joseph.Shaw@gosh.nhs.uk
## This script allow the Excel spreadsheet of the RAPID 
## biobank to be read into R as a tibble.
################################################################################

#############################################################
# Loading the biobank database
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
                              "Partner", "Pregnant woman"),
         
         # Add times
         
         time_of_blood_sample = format(time_of_blood_sample, "%I:%M"),
         date_of_blood_sample = as.POSIXct(date_of_blood_sample, format = "%Y-%m-%d", 
                           tz = ""),
         # Convert sampling time to date-time
         sampling_date_time = as.POSIXct(paste(date_of_blood_sample, 
                                               time_of_blood_sample), 
                                         format="%Y-%m-%d %H:%M",
                                         tz = "UTC"),
         
         time_first_spin = format(time_first_spin, "%I:%M"),
         date_of_first_spin = as.POSIXct(date_of_first_spin, format = "%Y-%m-%d", 
                                         tz = ""),
         
         first_spin_date_time = as.POSIXct(paste(date_of_first_spin, 
                                                 time_first_spin), 
                                           format="%Y-%m-%d %H:%M",
                                           tz = "UTC"),
         
         time_to_first_spin = round(difftime(first_spin_date_time, sampling_date_time,
                                             units = "hours")),
         
         date_of_2nd_spin = as.POSIXct(date_of_2nd_spin, format = "%Y-%m-%d", 
                                       tz = "") ,
         time_of_2nd_spin = format(time_of_2nd_spin, "%I:%M"),
         
         second_spin_date_time = as.POSIXct(paste(date_of_2nd_spin, 
                                                  time_of_2nd_spin), 
                                           format="%Y-%m-%d %H:%M",
                                           tz = "UTC"),
         
         # Assume samples were placed in storage after second spin
         time_to_storage = round(difftime(second_spin_date_time, sampling_date_time,
                                          units = "days")))

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



