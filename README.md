# Transcription factor ChIPseq data analysis pipeline

## Previous Version
https://github.com/RaimenChan/ChIPseq_pipeline

## This version has many changes:
1. The previous version required two replicates as input and selected peaks that appeared in both replicates for further analysis. This version analyzes each sample individually, without considering biological replicates.
2. Removed chipr analysis, correlation analysis, etc. These analyses require biological replicates.
3. Added MEME-ChIP.
4. Added `sample_target_gene_sorted.xlsx`, which can be used directly for publication.
5. New web output that summarizes all analysis results.

## Usage
Install the required packages.  
After modifying the parameters in `config.yaml` and the Snakemake file, run the following command:
```
snakemake -s TF_ChIPseq_snakemake_single_end_without_repeat_V5.py --cores 10 --keep-going --latency-wait 10
```
