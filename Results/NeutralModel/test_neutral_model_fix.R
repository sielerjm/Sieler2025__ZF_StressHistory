# Test script for the fixed neutral model function
# This script tests the modified neutral model function with the problematic data

# Set seed for reproducibility
set.seed(42)

# Load required libraries
library(minpack.lm)
library(Hmisc)
library(stats4)
library(phyloseq)
library(microViz)
library(tidyverse)

# Source the fixed neutral model function
source(here::here("Code", "Analysis", "NeutralModel", "neutralmodel_examplecode.R"))

# Load the phyloseq data
ps.tmp <- readRDS("/Users/michaelsieler/Dropbox/Mac (2)/Documents/Sharpton_Lab/Projects_Repository/Rules_of_Life/major-experiment-2023/Data/Robjects/pseq_uncleaned_05052025.rds")

# Clean the phyloseq object (same as in your original script)
ps.cleaned <-
    ps.tmp %>%
        ## Update Metadata
        ps_rename(Time = Timepoint) %>%
        microViz::ps_mutate(
            Treatment = case_when(
                Antibiotics == 0 & Temperature == 0 & Pathogen == 0 ~ "A- T- P-",
                Antibiotics == 0 & Temperature == 0 & Pathogen == 1 ~ "A- T- P+",
                Antibiotics == 1 & Temperature == 0 & Pathogen == 0 ~ "A+ T- P-",
                Antibiotics == 1 & Temperature == 0 & Pathogen == 1 ~ "A+ T- P+",
                Antibiotics == 0 & Temperature == 1 & Pathogen == 0 ~ "A- T+ P-",
                Antibiotics == 0 & Temperature == 1 & Pathogen == 1 ~ "A- T+ P+",
                Antibiotics == 1 & Temperature == 1 & Pathogen == 0 ~ "A+ T+ P-",
                Antibiotics == 1 & Temperature == 1 & Pathogen == 1 ~ "A+ T+ P+",
                TRUE ~ "Unknown"
            ), .after = "Pathogen"
        ) %>%
        microViz::ps_mutate(Sample = fecal.sample.number, .before = 1) %>%
        microViz::ps_mutate(Sample = gsub("^f", "", Sample)) %>%
        microViz::ps_filter(Treatment != "Unknown") %>%
        microViz::ps_mutate(
            History = case_when(
                Antibiotics + Temperature == 0 ~ 0,
                Antibiotics + Temperature == 1 ~ 1,
                Antibiotics + Temperature == 2 ~ 2,
            ), .after = "Treatment"
        ) %>%
        
        ## Additional metadata updates, factorizing metadata
        microViz::ps_mutate(
        # Create treatment code
            treatment_code = case_when(
              Antibiotics == 0 & Temperature == 0 & Pathogen == 0 ~ "Aneg_Tneg_Pneg",
              Antibiotics == 0 & Temperature == 0 & Pathogen == 1 ~ "Aneg_Tneg_Ppos",
              Antibiotics == 1 & Temperature == 0 & Pathogen == 0 ~ "Apos_Tneg_Pneg",
              Antibiotics == 1 & Temperature == 0 & Pathogen == 1 ~ "Apos_Tneg_Ppos",
              Antibiotics == 0 & Temperature == 1 & Pathogen == 0 ~ "Aneg_Tpos_Pneg",
              Antibiotics == 0 & Temperature == 1 & Pathogen == 1 ~ "Aneg_Tpos_Ppos",
              Antibiotics == 1 & Temperature == 1 & Pathogen == 0 ~ "Apos_Tpos_Pneg",
              Antibiotics == 1 & Temperature == 1 & Pathogen == 1 ~ "Apos_Tpos_Ppos"
            ),
            # Create treatment group factor
            treatment_group = case_when(
              Antibiotics == 0 & Temperature == 0 & Pathogen == 1 ~ "Parasite",
              Antibiotics == 1 & Temperature == 0 & Pathogen == 0 ~ "Antibiotics",
              Antibiotics == 1 & Temperature == 0 & Pathogen == 1 ~ "Antibiotics_Parasite",
              Antibiotics == 0 & Temperature == 1 & Pathogen == 0 ~ "Temperature",
              Antibiotics == 0 & Temperature == 1 & Pathogen == 1 ~ "Temperature_Parasite",
              Antibiotics == 1 & Temperature == 1 & Pathogen == 0 ~ "Antibiotics_Temperature",
              Antibiotics == 1 & Temperature == 1 & Pathogen == 1 ~ "Antibiotics_Temperature_Parasite",
              TRUE ~ "Control"
            ),
            # Convert to factor with appropriate levels
            treatment_group = factor(treatment_group, 
                                   levels = c("Control", "Parasite", 
                                              "Antibiotics", "Antibiotics_Parasite",
                                              "Temperature", "Temperature_Parasite",
                                            "Antibiotics_Temperature", "Antibiotics_Temperature_Parasite")
                                   ),
            treatment_code = factor(treatment_code, levels = treatment_order),
            # Create time point factor
            time_point = factor(Time, levels = c(0, 14, 18, 25, 29, 60)),
            # Create pathogen status factor
            pathogen_status = factor(ifelse(Pathogen == 1, "Exposed", "Unexposed"),
                                   levels = c("Unexposed", "Exposed")),
            # Create sex factor
            sex = factor(Sex, levels = c("M", "F"))
            )  %>%
    microViz::ps_mutate(Treatment = factor(Treatment, levels = treatment_order)) %>%
      microViz::ps_mutate(Exp_Type = case_when(
          Treatment %in% c("A- T- P-", "A- T- P+")  ~ "No prior stressor(s)",
          Treatment %in% c("A+ T- P-", "A+ T- P+")  ~ "Antibiotics",
          Treatment %in% c("A- T+ P-", "A- T+ P+") ~ "Temperature",
          Treatment %in% c("A+ T+ P-", "A+ T+ P+") ~ "Combined",
      )) %>%
      microViz::ps_mutate(Exp_Type = factor(Exp_Type, levels = c("No prior stressor(s)", "Antibiotics", "Temperature", "Combined"))) %>%
  # Fix names for taxonomic ranks not identified
  microViz::tax_fix(suffix_rank = "current", anon_unique = T, unknown = NA) %>% 
  # Filter for any samples that contain more than 5000 reads
  microViz::ps_filter(sample_sums(.) > 5000) %>%
  # Any taxa not found in at least 3 samples are removed
  microViz::tax_filter(min_prevalence = 3, undetected = 0) %>%
  # Remove any unwanted reads
  microViz::tax_select(c("Mitochondria", "Chloroplast", "Eukaryota"), deselect = TRUE) %>%
  microViz::tax_select(c("Bacteria, Phylum"), deselect = TRUE) 

# Filter for Day 60 data
ps.day60 <- ps.cleaned %>%
  microViz::ps_filter(Time == 60) 

# Extract OTU table for Day 60
otu_table_day60 <- as.data.frame(phyloseq::otu_table(ps.day60))
otu_table_day60 <- t(otu_table_day60)  # Transpose so samples are rows and taxa are columns

# Clean the OTU table
otu_table_clean <- otu_table_day60[rowSums(otu_table_day60) > 0, colSums(otu_table_day60) > 0]

cat("=== TESTING FIXED NEUTRAL MODEL FUNCTION ===\n")
cat("OTU table dimensions:", dim(otu_table_clean), "\n")
cat("Sample sums range:", range(rowSums(otu_table_clean)), "\n")

# Test the fixed neutral model function
cat("\n--- Testing with stats=TRUE ---\n")
tryCatch({
  neutral_stats <- sncm.fit(spp = otu_table_clean, stats = TRUE)
  cat("SUCCESS: Neutral model fitting completed!\n")
  cat("Migration rate (m):", neutral_stats$m, "\n")
  cat("R-squared:", neutral_stats$Rsqr, "\n")
  cat("AIC:", neutral_stats$AIC, "\n")
  print(neutral_stats)
}, error = function(e) {
  cat("ERROR in neutral model fitting:", e$message, "\n")
})

# Test with stats=FALSE to get predictions
cat("\n--- Testing with stats=FALSE ---\n")
tryCatch({
  tax_table_day60 <- as.data.frame(phyloseq::tax_table(ps.day60))
  neutral_predictions <- sncm.fit(spp = otu_table_clean, stats = FALSE, taxon = tax_table_day60[colnames(otu_table_clean), ])
  cat("SUCCESS: Neutral model predictions completed!\n")
  cat("Predictions dimensions:", dim(neutral_predictions), "\n")
  cat("First few predictions:\n")
  print(head(neutral_predictions, 5))
}, error = function(e) {
  cat("ERROR in neutral model predictions:", e$message, "\n")
})

cat("\n=== TEST COMPLETE ===\n")
cat("If both tests passed, the fixed neutral model function is working correctly.\n") 