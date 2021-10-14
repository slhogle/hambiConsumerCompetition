# Supporting data and code for the paper:

"Within-species trait diversity reorganizes consumer competitive hierarchy and prey community structure in an experimental microbial ecosystem"

[Preprint available from bioRxiv](https://www.biorxiv.org/content/10.1101/2021.06.09.447767v1)

Data and code here is provided under the MIT License. Feel free to use or remix as you see fit.

# PROCESSING 16S AMPLICON DATA

## DOWNLOAD
Install the [NCBI SRA Toolkit](https://github.com/ncbi/sra-tools). You can download a prebuilt binary [from here.](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit)

Access SRA data following the [instructions here.](https://github.com/ncbi/sra-tools/wiki/HowTo:-Access-SRA-Data)

You will need to [setup your configurations](https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration), but afterwards you could basically do:

```{bash}
prefetch SRR14323812
fasterq-dump SRR14323812
```

You will need to do this for all the SRA accessions associated with BioProject: [PRJNA725120](https://www.ncbi.nlm.nih.gov/bioproject/725120). For example:

```{bash}
cut -f1 00_SRA/sra_records.tsv | while read ID; do
  prefetch $ID
  fasterq-dump $ID
done
```

This will take some time...

## QUALITY CONTROL 16S AMPLICON DATA
### [MultiQC](https://multiqc.info/) sequence quality reports:
1. [Before processing/quality control](01_amplicon/01_sequencing_qc/qcreports_preprocess)
2. [After processing/quality control](01_amplicon/01_sequencing_qc/qcreports_postprocess)
3. [Mapping statistics](01_amplicon/02_mapping/qcreports_mapping)

### SCRIPTS
1. [QC](01_amplicon/01_sequencing_qc/bin)
2. [Mapping with bbmap](01_amplicon/02_mapping/bin)

## CALIBRATION 16S AMPLICON COUNT TABLES
Applying method from [this paper](https://elifesciences.org/articles/46923)

[Rendered workflow](01_amplicon/03_analysis/bin/01_bias_correction.md)

## PREY COUNTS RLOG NORMALIZATION
Using rlog normalization method from [this paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)

[Rendered workflow](01_amplicon/03_analysis/bin/02_count_normalization.md)

## CONSUMER DENSITY GENERALIZED ADDITIVE MODELS
1. [Ciliate model fits/checks](02_consumer_bacteria_biomass/bin/01_ciliate_density_gam.md)
2. [Nematode model fits/checks](02_consumer_bacteria_biomass/bin/02_worm_density_gam.md)
3. [Bacterial model fits/checks](02_consumer_bacteria_biomass/bin/03_OD600_density_gam.md)
4. [Plots and estimated marginal means](02_consumer_bacteria_biomass/bin/04_consumer_density_results.md)
	- Reproduces Fig. 2 from the Main text
5. [Fitting density data to Lotka-Volterra](02_consumer_bacteria_biomass/bin/04_consumer_density_results.md)
	- Reproduces Fig. 2 from the Main text

## BACTERIAL COMMUNITY COMPOSITION
1. [Area plots](03_bacteria_community_composition/bin/01_figS1_areaplot.md)
	- Reproduces area plot of bacterial relative abundances from Fig. S1 from the supplementary.
2. [Beta diversity](03_bacteria_community_composition/bin/02_community_dissimilarity.md)
	- Reproduces the community dissimiarlity analysis used to determine the reproducibility of the biological replicates including Fig. 3A from the main text and the beta regression.
3. [Shannon Diversity](03_bacteria_community_composition/bin/03_diversity_divnet.md)
	- Calculates Shannon diversity using the [DivNet package.](https://github.com/adw96/DivNet)
4. [Change point regression](03_bacteria_community_composition/bin/04_change_point_regression.md)
	- Multiple change point regression of community Shannon diversity using the [MCP package.](https://github.com/lindeloev/mcp). Figure 3C from the main text.
5. [Regressions on Shannon diversity](03_bacteria_community_composition/bin/05_figS2_covariate_effects_diversity.md)
	- Regressions for effect of consumer treatment on ratees of Shannon diversity decline. Reproduces Fig. S2 from the supplementary.
6. [Putting it all together](03_bacteria_community_composition/bin/06_Fig3.md)
	- Assembles components of main text Figure 3. 

## LONGITUDINAL GAUSSIAN PROCESS REGRESSION
Using the really nice [lgpr package.](https://jtimonen.github.io/lgpr-usage/index.html)

1. [Prepare data for LGPR](04_bacteria_species_lgpr/bin/01_lgpr_pre.md)
2. [Run data on the cluster](04_bacteria_species_lgpr/bin/02_lgpr_puhti.R)
	- This is not markdown. Also see `02_lgpr_puhti.sh`
3. [Make the plots](04_bacteria_species_lgpr/bin/03_lgpr_post.md)
	- Actual model fits with the posterior samplings are multiple GB and too big to put here. 
	- Reproduces Fig. 4 from the main text and Figs. S3, S4, and S5 from the supplementary. Most of the work is done with functions from `plotfuns.R` and `processfuns.R`

## TRAITS
Exploring trait data for use with HMSC

1. [Reproduce Fig. S6](05_bacteria_traits/bin/01_format.md)

## JOINT SPECIES DISTRIBUTION MODELING WITH HMSC
Using the [`HMSC-r` package](https://github.com/hmsc-r/HMSC)

1. [lmer test](06_hmsc/bin/00_lme4_sort_qeq_testing.md)
	- HMSC is pretty slow. Under the hood, HMSC is a bunch of linear mixed effects models that share information across species. You can prevent some headache by fitting individual models using `lmer` from the `lme4`. It just gives you a more decent idea about what you can/can't model well using the HMSC approach.
2. [Format models: sorting phase](06_hmsc/bin/01_format_sort.md)
3. [Format models: equilibrium phase](06_hmsc/bin/01_format_qeq.md)
4. Run on cluster
	- HMSC needs quite a bit of time to run so that MCMC chains converge on a reasonable posterior. Also since we also do 5-fold cross validation this is fitting the models five different times. This is best done on a compute cluster and not on your laptop. 
	- [Model fits](06_hmsc/bin/02_mcmc_puhti.R)
	- [Predictive performance by cross validation](06_hmsc/bin/02_crossvalidation_puhti.R)
5. [Check HMSC convergence](06_hmsc/bin/03_HMSC_convergence.md)
	- Important to do conergence checks so we are sure model converged properly and posterior estimates are stable.
6. [Interpret HMSC results](06_hmsc/bin/04_HMSC_inspect.md)
	- Code for generating Fig. 5 from the main text and supplementary tables S14 and S15.

## PREY CLEARANCE
1. [Fig. 6 Main text](07_consumer_feeding_efficiency/bin/01_traits_plots.md)
2. [Prey clearance statistics](07_consumer_feeding_efficiency/bin/02_traits_stats.md)
