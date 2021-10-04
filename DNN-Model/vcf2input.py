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
parser.add_argument('--fasta', help='FASTA file for genome used to align VCF', required = True)
parser.add_argument('--vcf', help='VCF file from which to run predict_cancer.py', required = True)
parser.add_argument('--output_dir', help='Path to directory in which to write output')
parser.add_argument('--sample_name', help='Sample naming to tag output')
args=parser.parse_args()

if args.output_dir == None:
    args.output_dir = "./"

if args.sample_name == None:
    args.sample_name = "TumorType_DNN"

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
        unw = ["X", "Y", "MT", "_", "GL", "KI"]
        if not any([ss in record.name for ss in unw]):
            seq_list[record.name] = record.seq
    return(seq_list)


def vcf2df(filename, seq_list):
    # read VCF files
    # filter variants in chromosomes in fasta
    # count only SNVs
    df = allel.vcf_to_dataframe(filename, fields='*', alt_number=2)
    print("--- %s Total Variants in the VCF file ---" % len(df))
    ##test for chr on seq_list, then add if not found in VCF or vice-versa
    seq_list_chr = "TRUE"
    try:
        seq_list["chr1"]
    except KeyError:
        seq_list_chr = "FALSE"
    chr_list = list(dict.keys(seq_list))
    if df['CHROM'].iloc[0].find("chr") < 0:
        if seq_list_chr == "FALSE":
            print("Fasta and VCF have no \'chr\' prefix...")
        else:
            print("Fasta has \'chr\' prefix, but VCF doesn't, fixing...")
            df['CHROM'] = "chr" + df['CHROM']
    else:
        if seq_list_chr == "FALSE":
            print("Fasta has no \'chr\' prefix, but VCF does, fixing...")
            df['CHROM'] = df['CHROM'].str.replace("chr", "")
    vcf_df = df[(df['is_snp'] == True) & (df['CHROM'].isin(chr_list) == True)]
    print("--- %s SNVs in the VCF file ---" % len(vcf_df))
    return(vcf_df)


def vcf2bins(tmp1, seq_list, sample_name):
    # set up data frame for SNV counts
    # load header
    ##set up bins from fasta seq_list set
    dkeys = list(dict.keys(seq_list))
    nochr_dkeys = [k.replace("chr", "") for k in dkeys]
    chr_bases = {}
    total_bases = 0
    seq_str = list(map(str, sorted(list(map(int, nochr_dkeys)))))
    if dkeys[0].find("chr") == 0:
        seq_str = ["chr" + k for k in seq_str]
    for i in seq_str:
        total_bases += len(seq_list[i])
        chr_bases[i] = len(seq_list[i])
    bin_bases = round(total_bases/2896.5)
    bases_df = pd.Series(chr_bases).apply(lambda x: int(round(float(x) / bin_bases)))
    binhd_tmp = {}
    a = 0
    for i in seq_str:
        for j in range(0,bases_df[i]):
            binhd_tmp[a] = str(i) + "." + str(j)
            a += 1
    binhd_df = pd.DataFrame(pd.Series(binhd_tmp))
    # initial dataframe
    test = tmp1.CHROM + '.' + \
        tmp1.POS.apply(lambda x: int(round(float(x) / bin_bases))).astype(str)
    bindf_tmp = pd.DataFrame(
        {"bins": pd.Series(pd.Categorical(test, categories=binhd_df.iloc[:, 0]))})
    bins_df = bindf_tmp.groupby('bins').size().reset_index(name = sample_name)
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
    tmp2 = tmp1.reset_index()
    tmp3 = tmp2.to_numpy()
    for i in range(len(tmp3)):
        ref = tmp3[i, 4]
        ch = tmp3[i, 1]
        st = tmp3[i, 2]
        alt = tmp3[i, 5]
        ref_vcf.append(ref)
        ref_base = seq_list[ch][st - 1].upper()
        ref_grch.append(ref_base)
        ref_cont = seq_list[ch][st - 2:st + 1].upper()
        if(not ref == ref_base):
            print("Mismatch in VCF and fasta genome versions, please ensure they are the same\n")
            print(ch, ":", st, "_", ref, "-", ref_base, "..", ref_cont)
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
        pd.Categorical(changes, categories=mut_df_header.iloc[:, 0]), dtype=object)})
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
    print("--- Loading FASTA File ---")
    start_time1 = time.time()
    seq_list = fasta_reader(args.fasta)
    print("--- %s seconds to load %s FASTA File ---" %
      ((time.time() - start_time1), len(seq_list)))

    # load VCF file:
    vcf_df = vcf2df(args.vcf, seq_list)
    sample_name = args.sample_name

    # convert VCF to Bin counts
    start_time = time.time()
    bins_df = vcf2bins(vcf_df, seq_list, sample_name)
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
    cb_df.to_csv(args.output_dir + "/" + sample_name + '.predict_cancer_input.csv', header=True)

    print("total --- %s seconds ---" % (time.time() - start_time1))
