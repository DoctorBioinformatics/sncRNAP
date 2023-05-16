
# sncRNAP

sncRNAs are known to be involved in post-transcriptional regulation of gene expression. Over the past few years, a group of sncRNA identification tools have been developed but none has shown the capacity to fully profile sncRNA and determine those that are differentially expressed in control (C) vs treated (T) samples. Therefore, a tool that profiles sncRNAs and identifies differentially expressed ones in group comparisons is needed. To profile sncRNAs and to determine those that are differentially expressed in C vs T samples, we developed sncRNAP, a Nextflow pipeline for the profiling and identification of differentially abundant sncRNAs from small RNAseq data. sncRNAP can be used for group comparisons such as C and T groups. The pipeline successfully identifies sncRNAs that are upregulated in C vs T samples (log2fc > 2, padj value < 0.005). Additionally, the pipeline profiles sncRNAs in each sample. In this manner, sncRNAP shows the length distribution, sncRNA counts per sample, overlapping sncRNA in each sample, sncRNA expression per sample and per sncRNA class,  volcano plots for each sncRNA class and heatmap for the top differentially expressed sncRNAs. Lastly, sncRNAP reports quality scores as well as the fasta sequence for the top identified candidates in the experiment.

![Presentation](https://user-images.githubusercontent.com/70538424/222744579-90cbe79e-3019-4ee1-81e5-69663298d648.png)

## Authors

- [Hesham Gibriel, PhD](https://github.com/hesham123457887)

- [Sharada Baindoor](https://github.com/SBaindoor)



## Installation

Run the below code to download sncRNAP:

```bash
git clone https://github.com/hesham123457887/sncRNAP.git
chmod 747 sncRNAP/bin/**
```

Run the below code to create a conda environment and install the required tools for the pipeline:

```bash
conda env create -f sncRNAP/environment.yml
conda activate sncRNAP
```
## Usage/Examples
Run the following to conduct a group comparison analysis:
```javascript
nextflow run sncRNAP \
--genome human \
--paired_samples FALSE \
--input_dir data/ \
--output_dir ./Results \
--layout layout/Layout.csv 
```

Run the following to conduct a paired group comparison analysis:
```javascript
nextflow run sncRNAP \
--genome human \
--paired_samples TRUE \
--input_dir data/ \
--output_dir ./Results \
--layout layout/Layout_paired_samples.csv 
```

## License

This project is licensed under the MIT License.
