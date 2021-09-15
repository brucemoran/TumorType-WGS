#! bash

mkdir -p /TumorType-WGS/DNN-Model/test/test_output
python3 /TumorType-WGS/DNN-Model/predict_cancer.py \
  --input /TumorType-WGS/DNN-Model/test/input_example.csv \
  --output "/TumorType-WGS/DNN-Model/test/test_output"
