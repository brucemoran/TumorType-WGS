#! python

import argparse
import time
import allel
import re
import concurrent.futures
import os
import sys
import pandas as pd
import numpy as np

from Bio import SeqIO
from multiprocessing import Pool

##help
parser = argparse.ArgumentParser(
    description='''Parse VCF into format for use with predict_cancer.py script ''')
parser.add_argument('--fasta', help='FASTA file for genome used to align VCF')
parser.add_argument('--vcf', help='VCF file from which to run predict_cancer.py')
parser.add_argument('--sampleID', help='String naming sample in VCF to tag output')
args=parser.parse_args()

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base)
                                 for base in reversed(seq))
    return(reverse_complement)


# def readAutoChrFASTA():
#     x = open(r'/TumorType-WGS/DNN-Model/chromosome.txt')
#     seq_list = {}
#     for i in range(23):
#         line = x.readline()
#         if(i < 9):
#             seq_list[line[0:4]] = line[5:]
#         else:
#             seq_list[line[0:5]] = line[6:]
#     x.close
#     return(seq_list)
#
#
##add to allow standard fasta format to be parsed into seq_list format
##replaces above readAutoChrFASTA()
def fasta_reader(fastafile):
    seq_list = {}
    for index, record in enumerate(SeqIO.parse(fastafile, "fasta")):
        seq_list[record.name] = record.seq
    return(seq_list)


def vcf2df(filename, seq_list):
    # read VCF files
    # filter variants in chromosomes in fasta
    # count only SNVs
    df = allel.vcf_to_dataframe(filename, fields='*', alt_number=2)
    print("--- %s Total Variants in the VCF file ---" % len(df))
    chr_list = dict.keys(seq_list)
    if df['CHROM'].iloc[0].find("chr") < 0:
        df['CHROM'] = "chr" + df['CHROM']
    vcf_df = df[(df['is_snp'] == True) & (df['CHROM'].isin(chr_list) == True)]
    print("--- %s SNVs in the VCF file ---" % len(vcf_df))
    return(vcf_df)


def vcf2bins(tmp1, sample_name):
    # set up data frame for SNV counts
    # load header
    filename = '/TumorType-WGS/DNN-Model/.1Mb.header.gz'
    binhd_df = pd.read_csv(filename, compression='gzip', header=None)
    # initial dataframe
    test = tmp1.CHROM + '.' + \
        tmp1.POS.apply(lambda x: int(round(float(x) / 1000000))).astype(str)
    bindf_tmp = pd.DataFrame(
        {"bins": pd.Series(pd.Categorical(test, categories=binhd_df.iloc[:, 0]))})
    bins_df = bindf_tmp.groupby('bins').size().reset_index(name=sample_name)
    # cap max bin count as 50
    tmp = bins_df.set_index("bins").to_numpy()
    tmp = np.clip(tmp, 0, 50)
    ndf = pd.DataFrame({'bins': bins_df.bins, sample_name: tmp[:, 0]})
    return(ndf)


def df2mut(tmp1, seq_list, sample_name):
    # Mutation types
    # Setup Dataframe
    # load header
    mut_df_header = pd.read_csv(
        '/TumorType-WGS/DNN-Model/Mut-Type-Header.csv')
    # initial dataframe
    changes = []
    ref_vcf = []
    ref_grch = []
    tmp1 = tmp1.reset_index()
    tmp1 = tmp1.to_numpy()
    for i in range(len(tmp1)):
        ref = tmp1[i, 4]
        ch = tmp1[i, 1]
        st = tmp1[i, 2]
        alt = tmp1[i, 5]
        ref_vcf.append(ref)
        ref_base = seq_list[ch][st - 1].upper()
        ref_grch.append(ref_base)
        ref_cont = seq_list[ch][st - 2:st + 1].upper()
        if(not ref == ref_base):
            print("Something wrong with the Genome Version, reference bases from VCF file doesn't match with records on hg19 genome\n")
            print(ch)
            print(st)
            print(ref)
            print(ref_base)
            print(ref_cont)
        if (re.search('[GT]', ref)):
            ref = reverse_complement(ref)
            ref_cont = reverse_complement(ref_cont)
            alt = reverse_complement(alt)
        change_sgl = ref + ".." + alt
        change_di1 = ref_cont[0:2] + ".." + ref_cont[0] + alt
        change_di2 = ref_cont[1:3] + ".." + alt + ref_cont[2]
        change_tri = ref_cont + ".." + ref_cont[0] + alt + ref_cont[2]
        changes.append(change_sgl)
        changes.append(change_di1)
        changes.append(change_di2)
        changes.append(change_tri)
    mutdf_tmp = pd.DataFrame({"bins": pd.Series(
        pd.Categorical(changes, categories=mut_df_header.iloc[:, 0]))})
    mut_df = mutdf_tmp.groupby('bins').size().reset_index(name=sample_name)
    if sum(mut_df[sample_name]) > 0:
        tmp_sgl = mut_df[sample_name].iloc[0:6] / \
            sum(mut_df[sample_name].iloc[0:6])
        tmp_di = mut_df[sample_name].iloc[6:54] / \
            sum(mut_df[sample_name].iloc[6:54])
        tmp_tri = mut_df[sample_name].iloc[54:150] / \
            sum(mut_df[sample_name].iloc[54:150])
        mut_df[sample_name].iloc[0:6] = tmp_sgl
        mut_df[sample_name].iloc[6:54] = tmp_di
        mut_df[sample_name].iloc[54:150] = tmp_tri
    return(mut_df)


if __name__ == '__main__':
    # load FASTA Files:
    start_time1 = time.time(args.fasta)
    seq_list = fasta_reader(args.fasta)
    print("--- %s seconds to load %s FASTA Files---" %
          ((time.time() - start_time1), len(seq_list)))

    # load VCF file:
    vcf_df = vcf2df(args.vcf)
    sample_name = args.sampleID

    # convert VCF to Bin counts
    start_time = time.time()
    bins_df = vcf2bins(vcf_df, sample_name)
    print(bins_df.shape)
    print("--- %s seconds to make SNV Counts Data Frame ---" %
          (time.time() - start_time))

    # convert VCF to Mutation Types
    start_time = time.time()
    mut_df = df2mut(vcf_df, seq_list, sample_name)

    print("--- %s seconds to make Mutation Type Data Frame ---" %
          (time.time() - start_time))

    cb_df = pd.concat([bins_df, mut_df])

    cb_df = cb_df.set_index('bins')

    cb_df = cb_df.transpose()
    cb_df.to_csv(args.sampleID + '.predict_cancer_input.csv', header=True)

    print("total --- %s seconds ---" % (time.time() - start_time1))
