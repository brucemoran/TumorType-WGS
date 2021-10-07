# TumorType-WGS
Classifying tumor types based on Whole Genome Sequencing (WGS) data

## Overview
This work is published in (Jiao et al. 2020)[https://www.nature.com/articles/s41467-019-13825-8] as part of the (PCAWG)[https://dcc.icgc.org/pcawg]
Here the original repo has been converted to a working tool that can be run in Docker. NB that only the `DNN-Model` was used and the `RF-Model` has been removed.


## Usage

### Bash wrapper
```
git clone https://github.com/brucemoran/TumorType-WGS && cd TumorType-WGS/
bash predict_from_vcf.sh \
  -v </path/to/your/vcf.vcf \
  -f </path/to/fasta.fa \
  -s <sample_name>
```

### Input
From the
"The input file should be a csv file of dimensionality Nx3047, where N = number of samples, and 3047 are the features ([mutation_distribution, mutation_types])."

This can be parsed from a VCF using `vcf2input.py`. To run:
```
/TumorType-WGS/DNN-Model/
```

## Running the DNN Model
N.B. I have no interest in RF Model so excluded it from this fork

```
To make a prediction run:

python predict_cancer.py --input_file <input.csv> --output <dir>
```
