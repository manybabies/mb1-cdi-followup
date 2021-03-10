library(tidyverse)
library(lme4)
library(brms); library(rstan); library(bridgesampling)
library(future); library(future.apply); library(furrr)
plan(multisession, workers = 8)
options("future" = T)

# LR_TESTS.LMER
# Perform likelihood ratio tests by sequential term deletion and anova model comparison
nested_models <- function(model, re_strat, verbose = F){
  # Extract formula (necessary for brmsfit) ======
  f <- formula(model)
  if(is.list(f)){
    f <- f$formula
  }
  f %<>% deparse() %>%
    str_c(collapse = "") %>%
    str_split("~|\\(|\\)", simplify = T) %>%
    str_trim() %>%
    discard(~ str_length(.x) == 0) %>%
    {.[-1]}
  fe.raw <- f[1]
  re.raw <- f[-1]
  # Extract fixed and random effects =============
  # TODO: extract terms directly by hand, following the output of terms.formula if possible
  fe <- fe.raw %>%
    str_split("\\+", simplify = T) %>%
    str_trim() %>%
    discard(~ str_length(.x) == 0) %>%
    rev() # Put the first effect to remove first
  re <- re.raw %>%
    pluck(1) %>%                          # Treat the first group of random effects only (for now)
    str_split(" [+|] ", simplify = T) %>% # Get separate effects
    str_trim() %>%
    discard(~ str_length(.x) == 0) %>%
    {.[-c(1, length(.))]}                 # Remove Intercept (first) and grouping variable (last)
  # Define model update function =================
  update.term_delete <- function(n_terms, fe, re){
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
  # Run models ===================================
  if(verbose){
    lr.verbose <- 1:length(fe) %>%
      map(quietly(~ update(model, update.term_delete(.x, fe, re)))) %>%
      prepend(model)
    return(lr.verbose)
  }
  nested_m <- 
  if(is.brmsfit(model)){
    nested_m <- 1:length(fe) %>%
      map(~ update(model, update.term_delete(.x, fe, re), refresh = 0)) %>%
      prepend(list(model))
  }else{
    nested_m <- 1:length(fe) %>%
      map(~ update(model, update.term_delete(.x, fe, re))) %>%
      prepend(model)
  }
  m_names <- c(fe, "Intercept") # if removing re together with fe
  return(list("models" = nested_m, "names" = m_names))
}

lr_tests.merMod <- function(model, re_strat, verbose = F){
  # Get (updated) nested models to compate, in descending order, and model names for summary
  nested <- nested_models(model, re_strat, verbose)
  nested_m <- nested[["models"]]
  m_names <- nested[["names"]]
  # Return model comparison with useful names
  return(exec(anova, !!!nested_m, model.names = m_names))
}

lr_tests.brmsfit <- function(model, re_strat, verbose = F){
  # Get (updated) nested models to compate, in descending order, and model names for summary
  nested <- nested_models(model, re_strat, verbose)
  nested_m <- nested[["models"]]
  m_names <- nested[["names"]] %>%
    discard(identical, "Intercept") # Remove "Intercept" from list
  # Bridge sample posteriors
  nested_bridges <- nested_m %>%
    future_map(bridge_sampler, silent = T)
  # Compare posteriors to compute Bayes factors
  bfs <- future_map2_dbl(nested_bridges[1:length(nested_bridges)-1],
                         nested_bridges[2:length(nested_bridges)],
                         ~ bayes_factor(.x, .y)$bf)
  # Create tibble to return Bayes factors and model names
  return(tibble(m_name = m_names, bf = bfs))
}
