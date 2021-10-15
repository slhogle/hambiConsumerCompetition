
# Supporting data and code for the paper:

"Effects of phenotypic variation on consumer coexistence and prey community structure"

[Preprint available from bioRxiv](https://www.biorxiv.org/content/10.1101/2021.06.09.447767v1)

[Full text available from Ecology Letters]()

Data and code is provided under the MIT License.

# Directory structure

1. `/sh` contains scripts from running analysis on the [puhti compute cluster](https://docs.csc.fi/computing/systems-puhti/)
2. `/R` contains R scripts
3. `/dataRaw` contains unprocessed data
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
cut -f1 data/sraRecords.tsv | while read ID; do
  prefetch $ID
  fasterq-dump $ID
done
```

## Quality control and mapping
Run these steps on an HPC cluster

1. [`ampliconQualityControl.sh`](sh/ampliconQualityControl.sh) -- Trim and filter reads
2. [`ampliconMapping.sh`](sh/ampliconMapping.sh) -- map QC'ed reads to the database in 

# Analysis steps
Go through these steps in order to reproduce the analysis in the paper. 

Note you must untar `rawData/16SAmplicon/mapping.tar.gz` first, then `rawData/16SAmplicon/mapping/bbmapRPKM.tar.gz`

## 1. Process amplicon
1. [`rpkm2tab.R`](R/rpkm2tab.R) -- format bbmap output to tables
2. [`correctBias.R`](R/correctBias.R) -- applying method from [this paper](https://elifesciences.org/articles/46923)
3. [`normalizeCounts.R`](R/normalizeCounts.R) -- normalizing sequencing counts for later 

## 2. Model consumer/prey densities
1. [`ciliateDensityGam.R`](R/ciliateDensityGam.R) -- fit GAM used in Table S2 and Fig 2
2. [`OD600DensityGam.R`](R/OD600DensityGam.R) -- fit GAM used in Table S3 and Fig 2
3. [`nematodeDensityGam.R`](R/nematodeDensityGam.R) -- fit GAM used in Table S4 and Fig 2
4. [`competitiveLV.R`](R/competitiveLV.R) -- parameterize LV competitive model
4. [`Fig2.R`](R/Fig2.R) -- Reproduces Fig. 2 from the Main text

## 3. Prey community composition