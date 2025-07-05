# Sitka spruce tree genetic analysis
This repository contains the code for an analysis of diameter at breast height (dbh21) in Sitka spruce trees using pedigree-based linear mixed models. The analysis focuses on estimating heritability, breeding values, and treatment effects across three generations of trees.

## Analysis Overview
Model: Individual animal model using ASReml with pedigree-based random effects

Fixed Effects: Root treatment (control, treatment 2, treatment 3)

Random Effects: Additive genetic effect, permanent environment
Heritability Estimate: ~0.71 (SE = 0.055)

Accuracy of EBVs: Calculated for key individuals and across generations

Model Diagnostics: Residual plots, Q-Q plots, Shapiro-Wilk tests

## Key Results
Treatment Effects: Small but statistically significant differences between treatments

EBVs: Individuals with no phenotypic data but shared parentage had identical EBVs

Accuracy: Highest in generation 1, lowest in generation 2 due to missing phenotypes

## Sitka spruce diameter measurements by generation and by treatment
![image](https://github.com/user-attachments/assets/8139885d-3890-4131-bc76-213b829327af)

## Summary of accuracies of EBV estimation across treatments and generations. 
![image](https://github.com/user-attachments/assets/75d5302e-5502-4c0f-9838-8a25d7236308)
