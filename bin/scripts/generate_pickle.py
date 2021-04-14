"""
@author Jasmin
"""

import pickle
import os
import pandas as pd
from collections import defaultdict


def parse(data_path):
    """
    This function creates dictionaries for the bed and bigwig files of the provided data. The dictionaries are stored in
    pickle files. For the files of ChIP-seq and ATAC-seq a separate pickle file is created for each biosource.
    """

    print('-----Generate pickle files-----')

    # read linking_table to get all available genomes, biosources and tfs
    lt = pd.read_csv(os.path.join(data_path, 'linking_table.csv'), sep=';',
                     usecols=['genome', 'epigenetic_mark', 'biosource'])
    genomes = set(lt.values[:, 0])
    lt_tfs = set(lt.values[:, 2][lt.values[:, 2] != ('dnasei' or 'dna accessibility')])
    lt_biosources = set(x for x in lt.values[:, 1])

    # go through every folder for genomes in the linking_table
    for genome in genomes:

        print('-Processing genome', genome)

        # create folder structure for pickle files if it does not exist
        if not os.path.exists(os.path.join(data_path, 'pickledata')):
            os.mkdir(os.path.join(data_path, 'pickledata'))
        if not os.path.exists(os.path.join(data_path, 'pickledata', genome)):
            os.mkdir(os.path.join(data_path, 'pickledata', genome))
        if not os.path.exists(os.path.join(data_path, 'pickledata', genome, 'chip-seq')):
            os.mkdir(os.path.join(data_path, 'pickledata', genome, 'chip-seq'))
        if not os.path.exists(os.path.join(data_path, 'pickledata', genome, 'atac-seq')):
            os.mkdir(os.path.join(data_path, 'pickledata', genome, 'atac-seq'))

        # list all biosource folders for one genome
        biosources = [x for x in os.listdir(os.path.join(data_path, genome)) if x in lt_biosources]

        for biosource in biosources:
            print('--Processing biosource', biosource)
            # dictionary for bed data of one biosource
            try:
                bs_bed_dict = pickle.load(open(os.path.join(data_path, 'pickledata', genome, 'chip-seq', biosource + ".pickle"), "rb"))
            except FileNotFoundError:
                bs_bed_dict = {}

            # list all transcription factor folders for one biosource
            tfs = [x for x in os.listdir(os.path.join(data_path, genome, biosource, 'chip-seq')) if x in lt_tfs]

            for tf in tfs:
                print('---Processing transcription factor', tf)
                # dictionary for bed data of one tf
                if bs_bed_dict and tf in bs_bed_dict:
                    tf_bed_dict = bs_bed_dict[tf]
                else:
                    tf_bed_dict = {}

                # list all files for chip data of one tf
                files = [x for x in os.listdir(os.path.join(data_path, genome, biosource, 'chip-seq', tf)) if x.lower().endswith('.bed')]

                # test for .bed ending of the file
                # bed files are read in with the function read_bed and the data is saved in the dictionary tf_bed_dict
                # the key for tf_bed_dict is the path of the associated bigwig-file
                count=1

                for bed_f in files:
                    print('Processing file {} of {}'.format(count, len(files)))
                    count += 1
                    bw = os.path.join(data_path, genome, biosource, 'chip-seq', tf, bed_f).replace('.bed','.bw')
                    if os.path.exists(bw) and bw not in tf_bed_dict:
                        tf_bed_dict[bw] = read_bed(os.path.join(data_path, genome, biosource, 'chip-seq', tf, bed_f))

                # the bed data of the tf is stored in the chip dictionary for the biosource with tf as key
                if tf_bed_dict:
                    bs_bed_dict[tf] = tf_bed_dict

            # a pickle file is created that contains the chip data of one biosource
            # it is named after the biosource
            if bs_bed_dict:
                with open(os.path.join(data_path, 'pickledata', genome, 'chip-seq', biosource + '.pickle'), 'wb') as handle:
                    pickle.dump(bs_bed_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

            try:
                atac = pickle.load(open(os.path.join(data_path, 'pickledata', genome, 'atac-seq', biosource + ".pickle"), "rb"))
            except FileNotFoundError:
                atac={}

            # list all files for atac data of one biosource
            # save the path to the bigwig file with the greatest size in dictionary atac; key is chromosome
            # create pickle file containing the path of the atac bigwig named after the biosource
            atac_chr_dict = defaultdict(dict)
            for f in os.listdir(os.path.join(data_path, genome, biosource, 'atac-seq')):
                if f.lower().endswith(('.bigwig', '.bigWig', '.bw')):
                    chr = f.split('.')[-2]
                    if chr not in atac:
                        atac_chr_dict[chr][f] = os.stat(os.path.join(data_path, genome, biosource, 'atac-seq',f)).st_size

            for c in atac_chr_dict:
                atac[c] = os.path.join(data_path, genome, biosource, 'atac-seq',max(atac_chr_dict[c], key=lambda key: atac_chr_dict[c][key]))
            
            with open(os.path.join(data_path, 'pickledata', genome, 'atac-seq', biosource + '.pickle'),
                      'wb') as handle:
                pickle.dump(atac, handle, protocol=pickle.HIGHEST_PROTOCOL)

def read_bed(file):
    """
    This function reads in a bed-file and returns the contained information in a dictionary.
    :param file: The path of the bed-file
    :return: chromosome is a dictionary with the chromosome as key and [start, stop, score, peak] as value
    """
    chromosome = defaultdict(list)

    f = pd.read_csv(file, sep='\t', usecols=lambda c: c in {'seqnames', 'start', 'end', 'PEAK'}, index_col=False)

    for index, row in f.iterrows():
        start = int(row['start'])
        end = int(row['end'])

        try:
            peak = row['PEAK']
        except KeyError:
            peak = -1

        if peak == -1:
            peak = int((end - start) / 2)

        # tests if a chromosome is already a key in the dictionary
        # if yes, the values are appended to the existing values of the key
        # if not, the key is added with the values
        chromosome[row['seqnames']].append([start, end, peak])

    return chromosome
