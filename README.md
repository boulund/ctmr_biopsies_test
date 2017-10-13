# CTMR biopsy metagenome pipeline
Metagenomics analysis pipeline under development, not for general use.

## Download the pipeline
Clone the repository

```
git clone https://github.com/boulund/ctmr_biopsies_test.git
```

## Run the pipeline on CTMR NAS
1. Setup the configuration of all paths to database files etc in `conf/ctmr_nas.config`.
2. Load the MetaPhlAn2 conda environment: `source activate metaphlan2`.
3. Run the pipeline: `nextflow run /path/to/main.nf -profile ctmr_nas --input_reads '/path/to/reads*_{1,2}.fq.gz'`.
4. The pipeline writes output files to `biopsy_pipeline_output`.
