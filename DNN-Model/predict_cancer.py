from __future__ import print_function
import tensorflow as tf
from tensorflow import keras
import random as rn
import numpy as np
import os

os.environ['PYTHONHASHSEED'] = '0'
# Setting the seed for numpy-generated random numbers
np.random.seed(45)

# Setting the graph-level random seed.
tf.random.set_seed(1337)

rn.seed(73)

from keras import backend as K

session_conf = tf.compat.v1.ConfigProto(
      intra_op_parallelism_threads=1,
      inter_op_parallelism_threads=1)

#Force Tensorflow to use a single thread
sess = tf.compat.v1.Session(graph=tf.compat.v1.get_default_graph(), config=session_conf)

K.set_session(sess)
import math
import pandas as pd
import argparse


import keras
from keras import backend as K
from keras.models import Sequential
from keras.layers import InputLayer, Input
from keras.layers import Dropout
from keras.layers import Dense
from keras.callbacks import TensorBoard
from tensorflow.keras.optimizers import Adam
from keras.models import load_model
from keras import regularizers


if __name__== "__main__":

    parser = argparse.ArgumentParser(description='Predict cancer type from mutation topology and mutation types', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input_csv', help='Input (use \'vcf2input.py\')', required = True)
    parser.add_argument('--output_dir', help='Path to directory in which to write output')
    args = parser.parse_args()

    if args.output_dir == None:
        args.output_dir = "./"

    d = pd.read_csv('/TumorType-WGS/DNN-Model/fold1_d.csv')
    factor = d.Factor
    cancer = d.Cancer
    factor_dict = dict(zip(factor, cancer))
    ##cat keras model.a* files and use that as model
    # filenames = ['model/ensemble_model.keras.aa', 'model/ensemble_model.keras.ab', 'model/ensemble_model.keras.ac', 'model/ensemble_model.keras.ad']
    # with open('ensemble_model.keras', 'w') as outfile:
    #     for fname in filenames:
    #         with open(fname) as infile:
    #             for line in infile:
    #                 outfile.write(line)
    os.system("cat /TumorType-WGS/DNN-Model/model/ensemble_model.keras.* >> /TumorType-WGS/DNN-Model/model/ensemble_model.keras")
    model = load_model('/TumorType-WGS/DNN-Model/model/ensemble_model.keras')
    input_file = args.input_csv
    output_dir = args.output_dir
    output_name = input_file.split("/")[-1].split(".")[0]
    data = pd.read_csv(input_file, index_col = [0])
    x_input = data.values
    if x_input.shape[-1] != 3047:
    	print('Input file not of correct dimensionality. See README for properly formatted file')
    predictions = model.predict(x_input)
    class_predictions = np.argmax(predictions, axis = 1)
    ### make predictions dataframe
    predictions_df = pd.DataFrame(data = predictions, index = data.index, columns = d.Cancer)
    predictions_df.to_csv(output_dir + "/" + output_name + ".cancer_prediction_probability.tsv", sep = '\t')
    cancer_classes = [factor_dict[i] for i in class_predictions]
    class_df = pd.DataFrame({'cancer_prediction':cancer_classes}, index = data.index)
    class_df.to_csv(output_dir + "/" + output_name + ".cancer_predictions.tsv", sep = '\t')
