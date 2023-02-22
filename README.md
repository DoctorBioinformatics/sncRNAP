
# sncRNAP

small non-coding RNA identification pipelines often lack the capacity to determine those that are differentially expressed in control vs treated samples.               
We developed sncRNAP, a Nextflow  pipeline for the detection of novel as well as known differentially expressed ncRNA and miRNAs from sRNAseq datasets.

![Slide1](https://user-images.githubusercontent.com/70538424/220607066-d3e01d96-23cb-445f-8247-d5f09e7da5bf.jpg)

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
nextflow run sncRNAP \
--genome human \
--paired_samples false \
--input_dir input/ \
--output_dir ./Results \
--layout ./layout.csv 
```


## License

This project is licensed under the MIT License.
