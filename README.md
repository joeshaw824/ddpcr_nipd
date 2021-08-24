# ddpcr_nipd
Analysis scripts for developing non-invasive prenatal diagnosis of maternal variants using ddPCR. Raw ddPCR data is exported as csv files from Quantasoft and stored on the W drive at Great Ormond Street Hospital.

Testing Structure

Each cfDNA test for a pathogenic variant comprises of two ddPCR assays: an assay to determine the variant fraction (vf_assay) and an assay to determine the fetal fraction (ff_assay).

The targets contained within each assay fall into 4 categories. The vf_assay detects the “variant” allele (i.e. “mutant”) and the “reference” allele (i.e. “wild type” or “normal”). The ff_assay detects a “maternal” allele, which is homozygous in the maternal cfDNA and heterozygous in the fetal cfDNA, and a “paternal” allele, which is not present in the maternal cfDNA and is heterozygous/hemizygous in the fetal cfDNA. 

Naming Convention

New variables are named in a consistent format: target_category_qualifier.
Target examples: variant, reference, paternal, maternal, fetal, difference, major_allele, minor_allele, ff_assay, vf_assay and total.
Category examples: molecules, fraction, percent, positives.
Qualifier examples: max, min.
Variable example: fetal_percent_max.

Sample Identification

Cell free DNA samples are named by the research number (R number) assigned in the RAPID Biobank. Parental genomic DNA samples are named by an episode number from either Omni (example: 18G07200) or Epic (example: 19RG-220G0191).

I have aimed to keep a consistent coding style as per tidyverse guidance:
https://style.tidyverse.org/

