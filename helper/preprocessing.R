library(tidyverse)

# READ_CDI
# Function that detects file type and calls corresponding readr function
read_cdi <- function(filename, save_unprocessed = F, verbose = F){
  if(verbose){print(filename)}
  filetype <- str_split(filename, "\\.")[[1]][2]
  # Define col_types for reading, except for standardized-score-CDI (to handle errors)
  col_types <- cols_only(
    labid = col_character(),
    subid = col_character(),
    `CDI-form` = col_character(),
    `CDI-agerange` = col_character(),
    `CDI-agedays` = col_integer(),
    vocab_nwords = col_integer(),
    `standardized-score-CDI` = col_guess(),
    `CDI-error` = col_character(),
    Notes = col_character()
  )
  # Read file
  df <- switch(filetype,
               csv = read_csv(paste0("data/cdi_raw/", filename),
                              col_types = col_types,
                              na = c("", "NA", "N/A")),
               tsv = read_tsv(paste0("data/cdi_raw/", filename),
                              col_types = col_types,
                              na = c("", "NA", "N/A"))) %>%
    as_tibble(.name_repair = "universal") # Replaces "-" by "." in column names
  # Save unprocessed data into global file if needed (typically for the first run)
  if(save_unprocessed){
    unprocessed <- "data/unprocessed_cdi_merged.csv"
    if(!file.exists(unprocessed)){
      ## Initialise unprocessed merge file (types don't matter)
      tibble(labid = character(),
             subid = character(),
             CDI.form = character(),
             CDI.agerange = character(),
             CDI.agedays = character(),
             vocab_nwords = character(),
             standardized.score.CDI = character(),
             CDI.error = character(),
             Notes = character(),
             .rows = 0) %>%
        write_csv(unprocessed)
    }
    ## Append data
    write_csv(df, unprocessed, append = T)
  }
  # Check standardized-score-CDI for input error
  if(!is_double(df$standardized.score.CDI)){
    df$standardized.score.CDI <- NA
  }
  return(df)
}

# READ_LABS
# Read lab information from the sign-up sheet
read_labs <- function(filename){
  # Define col_types for columns to import
  col_types <- cols_only(
    `lab ID code` = col_character(),
    `Language Community` = col_character()
  )
  # Read data, rename variables, recode language
  df <- read_csv(filename, col_types = col_types) %>%
    rename(labid = `lab ID code`,
           language = `Language Community`) %>%
    mutate(language_zone = case_when(str_detect(language, "British") ~ "British",
                                     str_detect(language, "American|Canadian") ~ "NAE",
                                     T ~ "Other")) %>%
    mutate_all(as_factor)
  return(df)
}
