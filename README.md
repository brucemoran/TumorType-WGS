# TumorType-WGS
Classifying tumor types based on Whole Genome Sequencing (WGS) data

## Overview
This work is published in [Jiao et al. 2020](https://www.nature.com/articles/s41467-019-13825-8) as part of the [PCAWG](https://dcc.icgc.org/pcawg). Here the [original repo](https://github.com/ICGC-TCGA-PanCancer/TumorType-WGS/commit/b79d090321f04840ae564d0117d315b3e0df0ff1) attached to the paper has been converted to a working tool that can be run by users. NB that only the `DNN-Model` was used and the `RF-Model` has been removed.

## Usage

### Bash wrapper
Clone the repo and run with your files. Output goes to the VCF directory.
```
git clone https://github.com/brucemoran/TumorType-WGS && cd TumorType-WGS/
bash predict_from_vcf.sh \
  -v </path/to/your/vcf.vcf \
  -f </path/to/fasta.fa \
  -s <sample_name>
```

### About Input to `predict_cancer.py`
"The input file should be a csv file of dimensionality Nx3047, where N = number of samples, and 3047 are the features `mutation_distribution` and `mutation_types`."

This is parsed from VCF using `vcf2input.py`. Please use the same fasta used to align the data, version mismatches will be caught and cause the run to fail.
