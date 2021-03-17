"""
Methods for merging reverse/forward reads in ATAC-seq data.


Goes through the following workflow:

- Reads in linkage_table.csv and identifies all ATAC-seq files with
  forward/reverse reads that need merging
- Groups them into pairs to be merged
- Checks if file format is bigWig, if not converts files to bigWig
- Merges files
- Checks if bedGraph is an allowed format since that is the output format of
  merging tool, if not files are converted to bigWig
- Deletes all files that are no longer needed
- Adds entries for merged files to linkage_table.csv


Use as follows:

    import merge_reads

    merge_reads.merge_all(linkage_table_path, chrom_sizes_paths,
    conversion_tool_path, merge_tool_path, allowed_file_formats)


by Kristina Müller (kmlr81)
"""

import os
import pandas as pd
import logging


def merge_all(linkage_table_path, chrom_sizes_paths, allowed_file_formats,
              conversion_tool="bedGraphToBigWig", merge_tool="bigWigMerge"):
    """
    Method merges all forward/reverse ATAC-seq files after converting them to
    bigWig format if necessary, checks if merged files need to be converted
    to bigWig and does so if true, then deletes old files nad adds entries for
    new files to the linkage table .csv file.

    Default values for conversion_tool and merge_tool are the tools that are 
    installed earlier in the pipe. If different tools are to be used, 
    they have to be installed manually.
    
    :param linkage_table_path: String containing path to linking_table.csv
    :param chrom_sizes_paths: Array of Strings containing paths to all
           chrom.sizes files
    :param allowed_file_formats: List of Strings with allowed file formats
    :param conversion_tool: String with path to bedGraphToBigWig tool
    :param merge_tool: String with path to bigWigMerge tool
    """

    linkage_frame = read_linkage_table(linkage_table_path)
    pairs = find_pairs(linkage_frame)

    pairs_merge = []
    converted_idxs = []
    merged_files = []
    old_files = []
    rows = []

    # If file format of forward/reverse reads is bedGraph, then convert to
    # bigWig for merging tool and make new pairs_merge List with only bigWig
    # file paths for merging
    for idx, pair in enumerate(pairs):
        tmp_pair = []
        for i in range(0, len(pairs[idx])):
            file_ext = os.path.splitext(pair[i])[1].lower()
            if ".bedgraph" == file_ext:
                genome = linkage_frame.loc[linkage_frame["file_path"] ==
                                           pair[i]]["genome"]
                # Find the right chrom.sizes file for genome
                chrom_sizes = [el for el in chrom_sizes_paths if
                               genome.iloc[0] in el]
                bw_file_path = convert_bedgraph_to_bigwig(pair[i],
                                                          chrom_sizes[0],
                                                          conversion_tool)
                tmp_pair.append(bw_file_path)
                converted_idxs.append([idx, i])
            else:
                tmp_pair.append(pair[i])
        pairs_merge.append(tmp_pair)

    # Merge all file pairs with tool and save paths of merged files
    for pair in pairs_merge:
        merged_files.append(merge_pair(pair[0], pair[1], merge_tool))

    # add option to convert merged files to bigWig format
    allowed_file_formats = [file_format.lower() for file_format in
                            allowed_file_formats]
    if "bedgraph" not in allowed_file_formats and ("bigwig" in
                                                   allowed_file_formats or
                                                   "bw" in
                                                   allowed_file_formats):
        tmp_paths = []
        for j in range(0, len(merged_files)):
            genome = linkage_frame.loc[linkage_frame["file_path"] ==
                                       pairs[j][0]]["genome"]
            chrom_sizes = [el for el in chrom_sizes_paths if
                           genome.iloc[0] in el]
            bw_file_path = convert_bedgraph_to_bigwig(merged_files[j],
                                                      chrom_sizes[0],
                                                      conversion_tool)
            tmp_paths.append(bw_file_path)
            old_files.append(merged_files[j])
        merged_files = tmp_paths

    # Check if there are any bedGraph merged files that were converted to
    # bigWig and delete them since they are no longer needed
    if len(old_files) > 0:
        for file in old_files:
            delete_file(file)

    # Delete all forward/reverse files in bedGraph format now that they have
    # been converted to bigWig
    for pair in pairs:
        for k in range(0, len(pair)):
            delete_file(pair[k])

    # Delete all forward/reverse files in bigWig format now that they have
    # been merged
    for idx_pair in converted_idxs:
        delete_file(pairs_merge[idx_pair[0]][idx_pair[1]])

    # Make rows to append to linkage table .csv file for merged files
    for m in range(0, len(merged_files)):
        row = pd.DataFrame(linkage_frame.loc[linkage_frame['file_path'] ==
                                             pairs[m][0]])
        row['file_path'] = merged_files[m]
        row['filename'] = os.path.basename(merged_files[m])
        rows.append(row)

    # Concat rows to a tmp data frame and then append that to the linkage
    # table .csv
    if len(rows) > 0:
        tmp_frame = rows[0]

        for n in range(1, len(rows)):
            tmp_frame = pd.concat([tmp_frame, rows[n]])

        tmp_frame.to_csv(linkage_table_path, sep=';', index=False, header=False,
                         mode='a')


def read_linkage_table(linkage_table_path):
    """
    Method reads in .csv linkage table file via the given file path and
    returns a data frame with all information for forward/reverse read ATAC-seq
    files.

    :return: linkage_frame: A data frame containing all information for
             forward/reverse reads files in need of merging
    """
    full_lf = pd.read_csv(linkage_table_path, sep=';')
    linkage_frame = full_lf.loc[(full_lf['technique'] == 'atac-seq') &
                                ((full_lf['filename'].str.contains('forward')) |
                                 (full_lf['filename'].str.contains('reverse')))]
    return linkage_frame


def find_pairs(linkage_frame):
    """
    Method pairs files that need to be merged.

    :param linkage_frame: Data frame containing all important info for files
                          to be merged
    :return: paris: A list of two dimensional lists containing the file paths to
                    the two files that need to be merged with each other
    """

    pairs = []

    for i in range(0, len(linkage_frame["file_path"]) - 1):
        filename = os.path.basename(linkage_frame["file_path"][i])
        file_id = filename.split(".")[0]

        for j in range(i + 1, len(linkage_frame["file_path"])):
            if linkage_frame["genome"][i] == linkage_frame["genome"][j] and \
                    linkage_frame["biosource"][i] == \
                    linkage_frame["biosource"][j] and file_id in \
                    linkage_frame["file_path"][j]:
                if 'chromosome' in linkage_frame.columns:
                    if linkage_frame["chromosome"][i] == \
                            linkage_frame["chromosome"][j]:
                        pairs.append([linkage_frame["file_path"][i],
                                      linkage_frame["file_path"][j]])
                else:
                    pairs.append([linkage_frame["file_path"][i],
                                  linkage_frame["file_path"][j]])

    return pairs


def convert_bedgraph_to_bigwig(bg_file_path, chrom_sizes_path,
                               conversion_tool_path):
    """
    Method converts a bedGraph file into a bigWig file.

    :param bg_file_path: String containing path to bedGraph file
    :param chrom_sizes_path: String containing path to chrom.sizes file
    :param conversion_tool_path: String containing path to conversion tool
    :return: Path to new bigWig file
    """

    bw_file_path = bg_file_path.rsplit(".", maxsplit=1)[0] + ".bw"

    if os.path.exists(bg_file_path):
        command = conversion_tool_path + " \"" + bg_file_path + "\" \"" + \
                  chrom_sizes_path + "\" \"" + bw_file_path + "\""

        os.system(command)
    # else for logging

    return bw_file_path


def merge_pair(bw_file_path_1, bw_file_path_2, merge_tool_path):
    """
    Method merges two bigWig forward/reverse read files into one bedGraph file

    :param bw_file_path_1: String containing path to first file
    :param bw_file_path_2: String containing path to second file
    :param merge_tool_path: String containing path to merging tool
    :return: Path to new merged file
    """

    merged_file_path = bw_file_path_1.split("_")[0] + "_merged.bedGraph"

    if os.path.exists(bw_file_path_1) and os.path.exists(bw_file_path_2):
        command = merge_tool_path + " \"" + bw_file_path_1 + "\" \"" + \
                  bw_file_path_2 + "\" \"" + merged_file_path + "\""

        os.system(command)
    # else für logging

    return merged_file_path


def delete_file(file_path):
    """
    Method deletes a file.

    :param file_path: Path to file in that needs to be deleted
    """
    if os.path.exists(file_path):
        os.remove(file_path)
    # add logging with else and error message
