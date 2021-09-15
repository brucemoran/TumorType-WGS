# TumorType-WGS
Classifying tumor types based on Whole Genome Sequencing (WGS) data

I have had some communication with PI and authors. I have also created a Dockerfile and docker is on the dockerhub: <dh>

## Docker
```
git clone https://github.com/brucemoran/TumorType-WGS
cd TumorType-WGS
docker build -t tumortype-wgs .
docker run -t tumortype-wgs
docker start tumortype-wgs
```

## Input
"The input file should be a csv file of dimensionality Nx3047, where N = number of samples, and 3047 are the features ([mutation_distribution, mutation_types])."

I am waiting for the script to parse these features from a 4.3 VCF.

## Running the DNN Model
N.B. I have no interest in RF Model so excluded it from this fork

```
To make a prediction run:

python predict_cancer.py --input_file <input.csv> --output <dir>
```
