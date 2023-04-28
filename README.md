# Droplet Digital PCR for Non-Invasive Prenatal Diagnosis

Shaw, J., Scotchman, E., Paternoster, B., Ramos, M., Nesbitt, S., Sheppard, S., Snowsill, T., Chitty, L. S. & Chandler, N. 2023. Non-invasive fetal genotyping for maternal alleles with droplet digital PCR: A comparative study of analytical approaches. Prenat Diagn.\
PMID: [36760169](https://pubmed.ncbi.nlm.nih.gov/36760169/)\
DOI: [10.1002/pd.6333](https://doi.org/10.1002/pd.6333)

These are analysis scripts for developing non-invasive prenatal diagnosis for maternal variants using ddPCR. ddPCR data was exported as csv files from Quantasoft (v1.7.4) and stored on the W drive at Great Ormond Street Hospital.

## Testing Structure

Each cfDNA test for a pathogenic variant comprises of two ddPCR assays: an assay to determine the variant fraction (vf_assay) and an assay to determine the fetal fraction (ff_assay).

The targets contained within each assay fall into 4 categories. The vf_assay detects the “variant” allele (i.e. “mutant”) and the “reference” allele (i.e. “wild type” or “normal”). The ff_assay detects a “maternal” allele, which is homozygous in the maternal cfDNA and heterozygous in the fetal cfDNA, and a “paternal” allele, which is not present in the maternal cfDNA and is heterozygous/hemizygous in the fetal cfDNA. 

## Bayesian Analysis
The Bayesian MCMC analysis section is compiled from R scripts and stan models kindly provided by Tristan Snowsill (University of Exeter), and run with cmdstan version 2.26.1.

## Naming Convention

New variables are named in a consistent format: target_category_qualifier.  

**Target examples:** variant, reference, paternal, maternal, fetal, difference, major_allele, minor_allele, ff_assay, vf_assay and total.  

**Category examples:** molecules, fraction, percent, positives.  

**Qualifier examples:** max, min.  

**Variable example:** fetal_percent_max.

## Sample Identification

Cell free DNA samples are named by the 4 or 5 digit research number (R number) assigned in the RAPID (Reliable, Accurate Prenatal, non-Invasive Diagnosis)
Biobank (example: 30216). Parental genomic DNA samples are named by an episode number from either the Omni (example: 18G07200) or Epic (example: 19RG-220G0191) laboratory information management systems.

## Disclaimer
This software is intended for research purposes only, and is not validated for clinical use. No guarantee is provided for third party usage.

## Licence
This project is licensed under the terms of the MIT license.
