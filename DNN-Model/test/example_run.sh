#! bash

perl pcawg_example_join.pl
mkdir -p /TumorType-WGS/DNN-Model/test/test_output
python3 /TumorType-WGS/DNN-Model/predict_cancer.py \
  --input /TumorType-WGS/DNN-Model/test/example_input.csv \
  --output "/TumorType-WGS/DNN-Model/test/test_output"
