"""
Methods for normalizing ATAC-seq and ChIP-seq signal values.

Goes through the following steps:
- Read in linkage table .csv file as data frame
- Log-scale all new files and save paths to log-files that have already
log-scaled
- min-max-scale all files (to a range of 0-1) with global min and max values
derived from the log-scaled values of all files in current analysis run

To-do:
- Bin size problem:
    - Narrow and broad peaks -> narrow peaks only has summit? One value for
      multiple pos? Problem that needs solving


Use as follows:

import normalize_signal_values

normalize_signal_values.normalize_all("example_linkage_table_path.csv")


by Kristina Müller (kmlr81)
"""

import math
import os
import pandas as pd
import numpy
import pyBigWig


def normalize_all(linkage_table_path):
    """
    Method normalizes all files through log scaling first and then
    min-max scaling to a uniform range between 0 and 1.

    :param linkage_table_path: String with path to linkage table .csv
           file containing the files that are part of the current analysis run.
    """
    linkage_table = pd.read_csv(linkage_table_path, sep=';')
    file_paths = list(linkage_table["file_path"])
    column_names= list(linkage_table["format"])
    log_file_paths = []
    min_value = math.inf
    max_value = -math.inf

    # Check all files to see if they have been log-scaled before or not.
    # If they have, then add path of existing .ln file to log_file_paths,
    # else log-scale and then add new path
    for i in range(0, len(file_paths)):
        log_file_path = None

        if os.path.exists(file_paths[i]):
            try:
                log_file_path = log_scale_file(file_paths[i], column_names[i])
            except FileNotFoundError:
                # Add logging
                excluded_files.append(i)
            except RuntimeError as err:
                # Add logging
                excluded_files.append(i)
        else:
            excluded_files.append(i)
            raise FileNotFoundError(
                "The file {} does not exist or the file path is "
                "incorrect.".format(file_paths[i]))

        if log_file_path is not None:
            log_file_paths.append(log_file_path)

    # Find global min and global max values
    for log_path in log_file_paths:
        min_value, max_value = get_min_max(log_path, min_val=min_value,
                                           max_val=max_value)

    # Min-max-scale all files
    for j in range(0, len(file_paths)):
        min_max_scale_file(file_paths[j], log_file_paths[j], min_value,
                           max_value, column_names=column_names[j])


def log_scale_file(file_path, column_names=None):
    """
    Method takes a file and writes a new file from it, which is a copy of the
    old one but with log scaled signal values in case the old file has bigWig
    format or a one column file containing only log signal values in case the
    old file has bed or bedGraph format. The new file is named after the
    old one with an additional ending .ln.

    :param file_path: String of path to file to be read in
    :param column_names: List of Strings with column names in file. Is set to
           None, if not given
    :return: log_file_path: String containing path to log-scaled file
    """
    log_file_path = file_path + ".ln"

    if is_big_wig(file_path):
        bw = pyBigWig.open(file_path)
        chrs = bw.chroms()
        header = list(chrs.items())
        bw_log = pyBigWig.open(log_file_path, 'w')
        bw_log.addHeader(header)

        for chrom in chrs.keys():
            intervals = bw.intervals(chrom)
            chromosomes = [chrom] * len(intervals)
            starts = [interval[0] for interval in intervals]
            ends = [interval[1] for interval in intervals]
            signal_values = [interval[2] for interval in intervals]
            log_values = numpy.log(signal_values)
            bw_log.addEntries(chromosomes, starts, ends=ends, values=log_values)

        bw.close()
        bw_log.close()

    else:
        idx = get_value_index(column_names)

        #check if file has a header that would interfere with using
        # numpy.loadtxt() and skip first row in file if needed
        with open(file_path, 'r') as file:
            first_line = file.readline()
        skip_rows = 0 if is_float(first_line.split('\t')[idx]) else 1

        signal_values = numpy.loadtxt(file_path, usecols=[idx],
                                      skiprows=skip_rows)
        log_values = numpy.log(signal_values)
        log_values = [str(value) + "\n" for value in log_values]

        with open(log_file_path, 'w') as log_file:
            log_file.writelines(log_values)

    return log_file_path


def get_min_max(log_file_path, min_val=math.inf, max_val=-math.inf):
    """
    Method finds min value and max value in a .log file.

    :param log_file_path: String with path to file with log values
    :param min_val: Start value for min. Is set to math.Inf if not given
    :param max_val: Start value for max. Is set to -math.Inf if not given
    :return: min and max as float
    """
    if is_big_wig(log_file_path):
        log_bw = pyBigWig.open(log_file_path)
        header = log_bw.header()
        tmp_min = header['minVal']
        tmp_max = header['maxVal']
        log_bw.close()

    else:
        signal_values = numpy.loadtxt(log_file_path, usecols=[0])
        tmp_min = min(signal_values)
        tmp_max = max(signal_values)

    min_val = tmp_min if tmp_min < min_val else min_val
    max_val = tmp_max if tmp_max > max_val else max_val

    return min_val, max_val


def min_max_scale_file(file_path, log_file_path, min_val,
                       max_val, column_names=None):
    """
    Method min-max scales values in file to a range between 0 and 1.

    :param file_path: String with path to file to be scaled
    :param log_file_path: String with path to file with log-scaled values of
           original file
    :param max_val: Global max value
    :param min_val: Global min value
    :param column_names: List of strings with names of columns in file,
           only needs to be given, if file is not bigWig format
    """
    tmp_file_path = file_path + ".tmp"

    if is_big_wig(file_path) and is_big_wig(log_file_path):
        bw = pyBigWig.open(file_path)
        log_bw = pyBigWig.open(log_file_path)
        chrs = bw.chroms()
        header = list(chrs.items())
        min_max_bw = pyBigWig.open(tmp_file_path, 'w')
        min_max_bw.addHeader(header)

        for chrom in chrs.keys():
            intervals = bw.intervals(chrom)
            log_intervals = log_bw.intervals(chrom)
            chromosomes = [chrom] * len(intervals)
            starts = [interval[0] for interval in intervals]
            ends = [interval[1] for interval in intervals]
            min_max_values = [(interval[2] - min_val) / (max_val - min_val) for
                              interval in log_intervals]
            min_max_bw.addEntries(chromosomes, starts, ends=ends,
                              values=min_max_values)

        bw.close()
        log_bw.close()
        min_max_bw.close()

    else:
        log_values = numpy.loadtxt(log_file_path, usecols=[0])
        min_max_values = [(x - min_val) / (max_val - min_val) for x in
                          log_values]
        idx = get_value_index(column_names)
        cnt = 0

        # Check if file has head and if so copy first line from current file
        # to tmp file and then start printing files with values
        with open(file_path, 'r') as file:
            first_line = file.readline()
        skip_rows = False if is_float(first_line.split('\t')[idx]) else True

        with open(file_path, 'r') as file, open(tmp_file_path, 'w') as tmp_file:
            if skip_rows:
                tmp_file.write(file.readline())

            for line in file:
                line_split = line.strip().split("\t")
                line_split[idx] = str(min_max_values[cnt])
                cnt += 1
                line_split.append("\n")
                line_write = "\t".join(line_split)
                tmp_file.write(line_write)

    os.rename(tmp_file_path, file_path)


def get_value_index(column_names):
    """
    Method gets column index of signal value in bed or bedGraph file.

    :param column_names: List of strings with names of columns in file
    :return: index: Integer representing index of column with signal value in
             file
    """
    if "SIGNAL_VALUE" in column_names:
        index = column_names.index("SIGNAL_VALUE")
    else:
        index = column_names.index("VALUE")

    return index


def is_big_wig(file_path):
    """
    Method checks if a file has bigWig format or not.

    :param file_path: String representing path to file
    :return: True if file has bigWig format, False if not
    """
    try:
        file_ext = os.path.splitext(file_path)[1].lower()
        if file_ext == '.bw' or file_ext == '.bigwig':
            bw = pyBigWig.open(file_path)
            is_bw = bw.isBigWig()
            bw.close()
            return is_bw
        else:
            return False
    except RuntimeError:
        return False
    # option for logging here


def is_float(value):
    """
    Checks if a string value can be converted to float.

    :param value: String to be checked
    :return: True, if value can be converted to float, False, if not
    """
    try:
        float(value)
        return True
    except ValueError:
        return False
