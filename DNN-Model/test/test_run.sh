#! bash

perl pcawg_test_join.pl
mkdir -p /TumorType-WGS/DNN-Model/test/test_output
python3 /TumorType-WGS/DNN-Model/predict_cancer.py \
  --input_csv /TumorType-WGS/DNN-Model/test/test_input.csv \
  --output_dir "/TumorType-WGS/DNN-Model/test/test_output" \
  --sample_name "PCAWG_test_output"
