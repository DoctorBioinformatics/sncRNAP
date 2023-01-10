
![Logo](https://dev-to-uploads.s3.amazonaws.com/uploads/articles/th5xamgrr6se0x5ro4g6.png)


# sncRNAP

small non-coding RNA identification pipelines often lack the capacity to determine those that are differentially expressed in control vs treated samples.               
We developed sncRNAP, a Nextflow  pipeline for the detection of novel as well as known differentially expressed ncRNA and miRNAs from sRNAseq datasets.

<img width="1440" alt="Screenshot 2023-01-10 at 18 31 56" src="https://user-images.githubusercontent.com/82383235/211608293-c101d063-08f2-4d97-a5c2-ba7eba179ec2.png">




## Authors

- [Hesham Gibriel, PhD](https://github.com/hesham123457887)

- [Sharada Baindoor](https://github.com/@sharadabaindoor1995)



## Installation

Run the below code to create a conda environment and install the required tools for the pipeline.

```bash
git clone https://github.com/hesham123457887/sncRNAP.git
conda env create -f sncRNAP/environment.yml
conda activate sncRNAP_env
```
## Usage/Examples
Run the following to conduct a group comparison analysis:
```javascript
nextflow run sncRNAP -profile conda --input '*fq.gz' 
--outdir ./Results --genome GRCm38 
--min_length 15 --trim_galore_max_length 50 
--mature "https://mirbase.org/ftp/CURRENT/mature.fa.gz" 
--hairpin "https://mirbase.org/ftp/CURRENT/hairpin.fa.gz" 
--layout ./layout.csv 

```


## License

This project is licensed under the MIT License.
