# mmcomp
Snakemake workflow for yield-normalized (synchronized) comparative genome-centric metagenomics

## Workflow description
* The workflow is designed to take in metagenomic sequencing datasets, subsample the reads to specified depths and perform microbial genome recovery with [mmlong2](https://github.com/Serka-M/mmlong2)
* Afterwards, additional analysis is performed to investigate possible reasons behind differences in genome recovery efficiency, which includes micro-diversity analysis, assessment of non-prokaryotic DNA, and community composition analysis
* Multiple bioinformatics tools are used in the analysis to correct for possible tool-related biases

## Tentative usage
* Download the repo and update the 'config/config.yaml' file to point to the correct Conda environments, read datasets, and databases
