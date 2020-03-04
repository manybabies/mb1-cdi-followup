library(tidyverse)
library(lme4)
library(furrr)
plan(multiprocess, workers = 4)

# LR_TESTS.LMER
# Perform likelihood ratio tests by sequential term deletion and anova model comparison
lr_tests.lmer <- function(model){
  # Extract fixed and random effects =============
  fe <- formula(model, fixed.only = T) %>%
    terms.formula(keep.order = T) %>%
    attributes() %>%
    pluck("term.labels") %>%
    rev() # Put the first effect to remove first
  re <- formula(model, random.only = T) %>%
    terms.formula(keep.order = T) %>%
    attributes() %>%
    pluck("term.labels") %>%              # Gets each group of random effects
    str_split(" [+|] ", simplify = T) %>% # Get separate effects (we only have one group)
    {.[2:(length(.)-1)]}                  # Remove Intercept (first) and grouping variable (last)
  # Define model update function =================
  update.term_delete <- function(n_terms){
    ## Get fixed effects to remove
    remove.fe <- fe[1:n_terms]
    ## Add any corresponding random effects
    remove.re <- cross2(re, remove.fe) %>%
      transpose() %>%
      pmap_chr(~ if_else(setequal(str_split(.x, ":", simplify = T),
                                  str_split(.y, ":", simplify = T)),
                         .x, NA_character_)) %>%
      discard(is.na)
    ## Combine into formula update
    update.fe.formula <- remove.fe %>%
      prepend("~ .") %>%
      str_c(collapse = " - ")
    remove.re.formula <- re %>%
      prepend("- (1") %>%
      str_c(collapse = " + ") %>%
      str_c(" | labid)")
    add.re.formula <- setdiff(re, remove.re) %>%
      prepend("+ (1") %>%
      str_c(collapse = " + ") %>%
      str_c(" | labid)")
    remove.formula <- paste(update.fe.formula, remove.re.formula, add.re.formula) %>%
      formula()
    return(remove.formula)
  }
  # Run and compare models =======================
  lr <- 1:length(fe) %>%
    map(~ update(model, update.term_delete(.x)))
  return(exec(anova, !!!lr))
}
