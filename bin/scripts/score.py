# -*- coding: utf-8 -*-
"""
Spyder Editor
Dies ist eine temporäre Skriptdatei.
"""

import pickle
import pyBigWig
import os
import logging

def findarea(width, genom, biosource_ls, tf_ls, chr_list, outpath):
    """
    This function creates a dictionary containing the mean for both CHIP-seq and ATAC-seq scores over a specified area. 
    This area needs to be specified at first, with the help of the CHIP-seq pickle and the value "w".
    """

    logging.info('starting calculation of score')
    print('------Calculating score------')

    # path to pickledata
    picklepath = os.path.abspath(
        os.path.join(outpath, 'data', 'pickledata'))

    calculateddict = {}

    # go through beddict for each biosource, then each tf, then each chromosom, then every binding
    # get Peak and Area from beddict and calculate the scores
    for biosource in biosource_ls:
        print("Analyzing biosource: ", biosource)
        # load dictionarys contaning paths to chip and atac bigwig files
        try:
            atacdict = pickle.load(open(os.path.join(picklepath, genom, 'atac-seq', biosource + ".pickle"), "rb"))
        except FileNotFoundError:
            atacdict = None
            logging.warning('There is no ATAC/DNase-seq data for biosource ' + biosource)
            print('-There is no ATAC/DNase-seq data for biosource ' + biosource)

        try:
            chipdict = pickle.load(open(os.path.join(picklepath, genom, 'chip-seq', biosource + ".pickle"), "rb"))
        except FileNotFoundError:
            chipdict = None
            logging.warning('There is no ChIP-seq data for biosource ' + biosource)
            print('-There is no ChIP-seq data for biosource ' + biosource)

        if atacdict and chipdict:
            # generate key for biosource if it does not exist
            if biosource not in calculateddict:
                calculateddict[biosource] = {}

            for tf in chipdict:
                print("-Analyzing transcription factor: ", tf)
                # test if tf was requested by the user
                if tf in tf_ls:

                    # generate key for tf if it does not exist
                    if tf not in calculateddict[biosource]:
                        calculateddict[biosource][tf] = {}

                    count = 1
                    for file in chipdict[tf]:

                        print('Calculating file {} of {}'.format(count, len(chipdict[tf])))
                        count += 1

                        try:
                            # open chip bigwig for tf
                            chip = pyBigWig.open(file)

                            for chromosom in chipdict[tf][file]:

                                # test if chromosome was requested by user
                                if chromosom in chr_list:

                                    # generate key for chromosome if it does not exist
                                    if chromosom not in calculateddict[biosource][tf]:
                                        calculateddict[biosource][tf][chromosom] = []

                                    # open atac bigwig
                                    atac = pyBigWig.open(atacdict[chromosom])

                                    for binding in chipdict[tf][file][chromosom]:

                                        start = binding[0]
                                        peak = binding[2]

                                        # calculate the area to be analyzed
                                        peaklocation = start + peak
                                        peaklocationstart = peaklocation - width
                                        peaklocationend = peaklocation + width

                                        # call scores between start and end from atac and chip using pyBigWig
                                        if chromosom in chip.chroms() and chromosom in atac.chroms():
                                            chip_score = chip.intervals(chromosom, peaklocationstart, peaklocationend)
                                            atac_score = atac.intervals(chromosom, peaklocationstart, peaklocationend)

                                            # calculate mean of chip and atac scores and save them into a dict
                                            calculationls = [peaklocationstart, peaklocationend]
                                            for i in (chip_score, atac_score):
                                                calculationls.append(
                                                    calculate_mean(i, peaklocationstart, peaklocationend))
                                            if calculationls not in calculateddict[biosource][tf][chromosom]:
                                                calculateddict[biosource][tf][chromosom].append(calculationls)
                        except RuntimeError:
                            logging.warning('Unable to open file ' + str(file))
                            print('Unable to open file ' + file)

                    print('Finished analysis of ', tf)

                    # remove key if the value is empty
                    if calculateddict[biosource][tf]:
                        pass
                    else:
                        del calculateddict[biosource][tf]

            # remove key if the value is empty
            if calculateddict[biosource]:
                pass
            else:
                del calculateddict[biosource]
            print("Finished analysis of", biosource)

    logging.info('finished calculation of score')
    return calculateddict


def calculate_mean(i, peaklocationstart, peaklocationend):
    """
    This function is for calculating the mean over the specified area and returns this mean.
    """
    length = peaklocationend - peaklocationstart
    mean = 0
    if i:
        for interval in i:
            if interval[0] < peaklocationstart and interval[1] > peaklocationend:
                interval_length = peaklocationend - peaklocationstart
            else:
                if interval[1] > peaklocationend:
                    interval_length = peaklocationend - interval[0]
                elif interval[0] < peaklocationstart:
                    interval_length = interval[1] - peaklocationstart
                else:
                    interval_length = interval[1] - interval[0]
            mean += interval_length * interval[2]
    mean = mean / length
    return mean
