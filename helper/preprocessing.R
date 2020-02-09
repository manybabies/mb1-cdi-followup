library(tidyverse)

# READ_CDI
# Function that detects file type and calls corresponding readr function
read_cdi <- function(filename){
  filetype <- str_split(filename, "\\.")[[1]][2]
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
  df <- switch(filetype,
               csv = read_csv(paste0("data/cdi/", filename), col_types = col_types),
               tsv = read_tsv(paste0("data/cdi/", filename), col_types = col_types))
  if(is_character(df$`standardized-score-CDI`)){
    df$`standardized-score-CDI` = NA
  }
  return(df)
}
