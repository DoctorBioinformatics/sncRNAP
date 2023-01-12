
# sncRNAP

small non-coding RNA identification pipelines often lack the capacity to determine those that are differentially expressed in control vs treated samples.               
We developed sncRNAP, a Nextflow  pipeline for the detection of novel as well as known differentially expressed ncRNA and miRNAs from sRNAseq datasets.

![Fig2A](https://user-images.githubusercontent.com/82383235/211609047-df6e06f8-3fc5-43f0-b8e4-af0984fcd173.png)

## Authors

- [Hesham Gibriel, PhD](https://github.com/hesham123457887)

- [Sharada Baindoor](https://github.com/@sharadabaindoor1995)



## Installation

Run the below code to create a conda environment and install the required tools for the pipeline.

```bash
git clone https://github.com/hesham123457887/sncRNAP.git
conda env create -f sncRNAP/environment.yml
conda activate sncRNAP
```
## Usage/Examples
Run the following to conduct a group comparison analysis:
```javascript
nextflow run sncRNAP -profile conda --input_dir '*fq.gz' 
--output_dir ./Results --genome human 
--min_length 15 --trim_galore_max_length 50 
--mature "https://mirbase.org/ftp/CURRENT/mature.fa.gz" 
--hairpin "https://mirbase.org/ftp/CURRENT/hairpin.fa.gz" 
--layout ./layout.csv 

```


## License

This project is licensed under the MIT License.
