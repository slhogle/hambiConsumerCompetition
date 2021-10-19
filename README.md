# "Effects of phenotypic variation on consumer coexistence and prey community structure"

[Preprint available from bioRxiv](https://www.biorxiv.org/content/10.1101/2021.06.09.447767)

[Full text available from Ecology Letters]()

Supporting data and code is provided under the MIT License.

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

Go through these steps in order to reproduce the analysis in the paper. 

## 1. Format data
1. [`rpkm2tab.R`](r/rpkm2tab.R) -- format bbmap output to tables
2. [`format_raw_data.R`](r/format_raw_data.R) -- format raw density data for downstream use

## 2. Process amplicon
2. [`correct_bias.R`](r/correct_bias.R) -- applying method from [this paper](https://elifesciences.org/articles/46923)
3. [`normalize_counts.R`](r/normalize_counts.R) -- normalizing sequencing counts for later 

## 3. Model consumer/prey densities
1. [`gam_ciliate_density.R`](r/gam_ciliate_density.R) -- fit GAM used in Table S2 and Fig 2
2. [`gam_OD600_density.R`](r/gam_OD600_density.R) -- fit GAM used in Table S3 and Fig 2
3. [`gam_nematode_density.R`](r/gam_nematode_density.R) -- fit GAM used in Table S4 and Fig 2
4. [`competitive_LV.R`](r/competitive_LV.R) -- parameterize LV competitive model
5. [`fig2.R`](r/fig2.R) -- Reproduces Fig. 2 from the main text

## 4. Prey community composition
1. [`figS1a.R`](r/figS1a.R) -- Reproduces Fig. S1A from supplementary
2. [`community_dissimilarity_plot.R`](r/community_dissimilarity_plot.R) -- Reproduces Fig. S1B from supplementary
3. [`community_dissimilarity_regression.R`](r/community_dissimilarity_regression.R) -- Reproduces Table S5
4. [`shannon_diversity.R`](r/shannon_diversity.R) -- Run divnet estimate of Shannon diversity and save result
5. [`shannon_change_point_regression.R`](r/shannon_change_point_regression.R) -- Perform multiple change point regression of the shannon diversity estimate. Reproduces Fig. S2.
6. [`shannon_sort_eq_regression_bayes.R`](r/shannon_sort_eq_regression_bayes.R) -- Performs regression using MCMC from Stan on the Shannon estimates in the two experimental phases. Reproduces Fig. S3 and Tables S6 and S7
  a. [`shannon_sort_eq_regression_freq.R`](r/shannon_sort_eq_regression_freq.R) -- Performs frequentist regression. No Stan required
7. [`ordination.R`](r/ordination.R) -- Ordination and jump lengths. Generate Fig. 3C and Table S8
8. [`fig3.R`](r/fig3.R) -- Generate final Fig 3

## 5. Bacteria traits
1. ['traits.R'](r/traits.R) -- Load and format trait data for use with JSDM models

## 6. Joint species distribution modeling
1. JSDM_format_sort.R
2. JSDM_format_qeq.R
2. puhti_hmsc_mcmc.sh executes puhti_hmsc_mcmc.R
3. puhti_hmsc_crossvaliation.sh executes puhti_hmsc_crossvaliation.R
4. JSDM_convergence.R
5. JSDM_inspect.R

## 7. Consumer feeding efficiency