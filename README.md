# mb1-cdi-followup
ManyBabies 1 Longitudinal CDI Followup

This repository contains the data, code and output for the ManyBabies 1 CDI Follow-up project. This project followed infants from the original ManyBabies1 project (which examined experimentally infants' preference for Infant-Directed Speech) and measured their vocabulary at 18 and/or 24 months using parental report checklists (CDIs).

https://manybabies.github.io/MB1/

To review the analysis, follow the rmd/html files in numerical order:

01_read_and_merge: Merging of the original data from MB1 with the CDI score data.
02_exclusion_and_new_variables: Implementing exclusion criteria and creation of new derivative variables.
03a_confirmatory_analysis: Pre-registered correlations and linear models.
03b_confirmatory_analysis: Pre-registered Bayes analysis. Note: this file requires additional software installation to run.
04_compute_percentiles_mb_cdi: Creation of new daily percentile scores for all participants across language groups.
05_main_exploratory: Main "exploratory" analysis using computed daily percentile scores for all language groups.
06_other-exploratory_analyses: Additional analyses not preregistered.
