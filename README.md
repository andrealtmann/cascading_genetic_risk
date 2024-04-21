# Cascading Genetic Risk

This repositoroy contains the scripts to conduct the analysis described in /Towards Cascading Genetic Risk in Alzheimer's Disease/ (https://www.medrxiv.org/content/10.1101/2023.12.16.23300062v1). The required data can be obtained from the ADNI (https://adni.loni.usc.edu/data-samples/access-data/).

The R scripts contained in this repository:
* ggcoxzph_fixed.R - implements a fixed version of the ggcoxzph function based on  https://stats.stackexchange.com/questions/560975/how-to-interpret-schoenfield-residual-plot
* prepare_csf_clean.R - prepased a 'clean' version of the CSF biomaker measures using the UPENNBIOMK_MASTER_FINAL file from ADNI
* analysis_clean.R - read biomaker data (CSF and PET data) and prepare longitudinal data
* analysis_conversion1.R - analysis for A-T- to A+T- conversion
* analysis_conversion2.R - analysis for A+T- to A+T+ conversion
* screen_cutoffs_clean.R - wrapper script to screen cutoffs for tau PET and pTau

In order to conduct the analyis, run:
1. prepare_csf_clean.R - prepares the cleaned CSF data table
2. analysis_clean.R    - prepare longitudinal data
3. analysis_conversion1.R - Survival Model for A-T- to A+T- conversion
4. analysis_conversion1.R - Survival Model for A+T- to A+T+ conversion
5. screen_cutoffs_clean.R - Screen tau cutoffs 

Options (analysis_clean.R):
* to select the used biomakrers (CSF, PET, or both [default] )
* to add genetic population structure as a covariate (default: off)
* use GWAS summary stats by Bellenguez et al (2022) instead of Kunkle et al. (2019)
* select the PRS cutoff (default: 1e-08)

Options (analysis_conversion1.R):
* conduct a sensitivity analysis for conversion criteria (A+ instead of A+T-)
* Compute C-Index (default: off)
* Produce Figures (default: off)

Options (analysis_conversion2.R):
* conduct a sensitivity analysis for conversion criteria (T+ instead of A+T+)
* Compute C-Index (default: off)
* Produce Figures (default: off)

PRS scores and genetic resuls have been generated as discussed in the manuscript and are available upon request:
* ADNI1_3_genetic_summary_CEU80_rel.csv
* ADNI1-3_bellenguez_avgPRS_mod2.all.score

Following R libraries are required:
* data.table
* ADNIMERGE
* survival
* survminer
* cutpointr
* rms
