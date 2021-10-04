#! bash

##wrapper to run vcf2input.py, predict_cancer.py inside docker

if [[ $(docker --version | wc -l) != 1 ]]; then
  echo "Docker doesn't seem to be running, exiting..."
else

  DOCKER="brucemoran/tumortype-wgs"
  FOUND=$(docker images | grep ${DOCKER} | wc -l)

  if [[ $FOUND < 1 ]]; then
    docker pull ${DOCKER}
  fi

    ##input flags
    while getopts "v:f:s:" flag; do
      case "${flag}" in
          v) VCF=${OPTARG};;
          f) fasta=${OPTARG};;
          s) sample_name=${OPTARG};;
      esac
    done

  ##need to mount dirs of input and output
  MNL="--mount type=bind,source="
  MNO=",target=/mnt\n"
  MNT="$(echo -e "$MNL"$(dirname $VCF)"$MNO\n$MNL"$(dirname $fasta)"$MNO\n" | sort | uniq)"
  CMD="python3 /TumorType-WGS/DNN-Model/vcf2input.py \
        --vcf /mnt/$(basename "${VCF}") \
        --fasta /mnt/$(basename "${fasta}") \
        --sample_name ${sample_name} \
        --output_dir /mnt; \
       python3 /TumorType-WGS/DNN-Model/predict_cancer.py \
        --input_csv /mnt/${sample_name}.predict_cancer_input.csv \
        --output_dir /mnt"
  echo -e "Command to be run:\n$CMD"
  docker run ${MNT} ${DOCKER} bash -c "${CMD}"

fi
