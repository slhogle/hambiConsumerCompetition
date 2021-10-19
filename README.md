# Effects of phenotypic variation on consumer coexistence and prey community structure

[Preprint available from bioRxiv](https://www.biorxiv.org/content/10.1101/2021.06.09.447767)

[Full text available from Ecology Letters]()

Supporting data and code is provided under the MIT License.

# Clone the repository
```{bash}
git clone https://github.com/slhogle/hambiConsumerCompetition.git
```

# Directory structure

1. `/sh` contains scripts from running analysis on the [puhti compute cluster](https://docs.csc.fi/computing/systems-puhti/)
2. `/r` contains R scripts
3. `/data_raw` contains unprocessed data
4. `/data` contains data that has been processed in some way for from `/dataRaw` for downstream use
5. `/figs` contains figures generated from R scripts
6. `/tables` contains summary tables generated from R scripts
7. `/rendered` contains any markdown or HTML that gets rendered in the process of running the scripts.

# 16S amplicon processing

## Download
Install the [NCBI SRA Toolkit](https://github.com/ncbi/sra-tools). You can download a prebuilt binary [from here.](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit)

Access SRA data following the [instructions here.](https://github.com/ncbi/sra-tools/wiki/HowTo:-Access-SRA-Data)

You will need to [setup your configurations](https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration), but afterwards you could basically do:

```{bash}
prefetch SRR14323812
fasterq-dump SRR14323812
```

You will need to do this for all the SRA accessions associated with BioProject: [PRJNA725120](https://www.ncbi.nlm.nih.gov/bioproject/725120). For example:

```{bash}
cut -f1 data/sra_records.tsv | while read ID; do
  prefetch $ID
  fasterq-dump $ID
done
```

## Quality control and mapping
Run these steps on an HPC cluster

1. [`amplicon_quality_control.sh`](sh/amplicon_quality_control.sh) -- Trim and filter reads
2. [`amplicon_mapping.sh`](sh/amplicon_mapping.sh) -- map QC'ed reads to the database in `data_raw/16S_amplicon/mapping.tar.gz`

# Analysis steps

```
R version 4.1.1 (2021-08-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Pop!_OS 21.04
```

Project dependencies are managed with [`pak`](https://pak.r-lib.org/) and [`renv`](https://rstudio.github.io/renv/index.html). This allows one to restore the
correct package versions for this project to run with the system and R version
listed above.

To setup packages to the versions used in this analysis, simply run:

```{r}
renv::restore()
```

from the main project directory.


Go through these steps in order to reproduce the analysis in the paper. 

## 1. Format data
1. [`rpkm2tab.R`](R/rpkm2tab.R) -- format bbmap output to tables
2. [`format_raw_data.R`](R/format_raw_data.R) -- format raw density data for downstream use

## 2. Process amplicon
2. [`correct_bias.R`](R/correct_bias.R) -- applying method from [this paper](https://elifesciences.org/articles/46923)
3. [`normalize_counts.R`](R/normalize_counts.R) -- normalizing sequencing counts for later 

## 3. Model consumer/prey densities
1. [`gam_ciliate_density.R`](R/gam_ciliate_density.R) -- fit GAM used in Table S2 and Fig 2
2. [`gam_OD600_density.R`](R/gam_OD600_density.R) -- fit GAM used in Table S3 and Fig 2
3. [`gam_nematode_density.R`](R/gam_nematode_density.R) -- fit GAM used in Table S4 and Fig 2
4. [`competitive_LV.R`](R/competitive_LV.R) -- parameterize LV competitive model
5. [`fig2.R`](R/fig2.R) -- Reproduces Fig. 2 from the main text

## 4. Prey community composition
1. [`figS1a.R`](R/figS1a.R) -- Reproduces Fig. S1A from supplementary
2. [`community_dissimilarity_plot.R`](R/community_dissimilarity_plot.R) -- Reproduces Fig. S1B from supplementary
3. [`community_dissimilarity_regression.R`](R/community_dissimilarity_regression.R) -- Reproduces Table S5
4. [`shannon_diversity.R`](R/shannon_diversity.R) -- Run divnet estimate of Shannon diversity and save result
5. [`shannon_change_point_regression.R`](R/shannon_change_point_regression.R) -- Perform multiple change point regression of the shannon diversity estimate. Reproduces Fig. S2.
6. [`shannon_sort_eq_regression_bayes.R`](R/shannon_sort_eq_regression_bayes.R) -- Performs regression using MCMC from Stan on the Shannon estimates in the two experimental phases. Reproduces Fig. S3 and Tables S6 and S7
  a. [`shannon_sort_eq_regression_freq.R`](R/shannon_sort_eq_regression_freq.R) -- Performs frequentist regression. No Stan required
7. [`ordination.R`](R/ordination.R) -- Ordination and jump lengths. Generate Fig. 3C and Table S8
8. [`fig3.R`](R/fig3.R) -- Generate final Fig 3

## 5. Bacteria traits
1. [`traits.R`](R/traits.R) -- Load and format trait data for use with JSDM models

## 6. Joint species distribution modeling
1. [`JSDM_format_sort.R`](R/JSDM_format_sort.R) -- formats sorting phase models. Models need to be uploaded to cluster afterward
2. [`JSDM_format_qeq.R`](R/JSDM_format_qeq.R) -- formats equilibrium phase models. Models need to be uploaded to cluster afterward
2. [`puhti_hmsc_mcmc.sh`](sh/puhti_hmsc_mcmc.sh) -- executes [`puhti_hmsc_mcmc.R`](R/puhti_hmsc_mcmc.R) and fits models using MCMC
3. [`puhti_hmsc_crossvaliation.sh`](sh/puhti_hmsc_crossvaliation.sh) executes [`puhti_hmsc_crossvaliation.R`](R/puhti_hmsc_crossvaliation.R) and performs 5-fold cross validation to determine predictive performance of the full models
4. [`JSDM_convergence.R`](R/JSDM_convergence.R) -- Various stats to demonstrate chains have converged
5. [`JSDM_inspect.R`](R/JSDM_inspect.R) -- Produces Table S9 and S10 and main text Fig. 4

## 7. Consumer feeding efficiency
1. [`fig5.R`](R/fig5.R) -- generates Figure 5 from main text
2. [`consumer_feeding_stats.R`](R/consumer_feeding_stats.R) -- reproduces Table S11 and S12