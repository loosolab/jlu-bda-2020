"""
Wrapper for the data preparation phase.

Goes through the following steps:
- call generate_linking_table.r to generate a linking table
- call download_deepblue_data.r to download the data
- create a list of all files that need to be validated and then call
  convert_files.sh to validate and convert them
- call merge_reads.py to merge forward reverse reads
- call sort_files.sh to sort files into the final folderstructure
- call generate_pickle.py to generate picklefiles
- create a list of files that need to be normalized and then call
  normalize_signal_values.py to normalize the data

Use as follows:

import generate_data
data = DataConfig(genome, chromosome, biosource, tf, path, csvname,
              datatype, localfiles, redoverification, offline, logfile)
data.pull_data()

by Jonathan Schaefer(JonnyCodewalker)
"""

import subprocess
import os
import logging
from csv import reader, writer
from scripts.merge_reads import merge_all
from scripts.generate_pickle import parse
from scripts.normalize_signal_values import normalize_all
from scripts.setup_logging import setup


class DataConfig:
    """Contains configuration data for the pull_data function.
    - genome: list of genomes
    - chromosome: list of chromosomes
    - biosource: list of biosources
    - epigenetic_mark: list of TFs
    - output_path: path for where to store data
    - csv_name: what to name the linking table
    - datatype: what filetype to convert to
    - localfiles: path to where localfiles are if used
    - redoverification: boolean to check if force verification is done
    - offline: boolean to check if offline mode is used
    - logilfe: path to the logfile used in the run
      """

    def __init__(self, genome, chromosome, biosource, epigenetic_mark,
                 output_path, csv_name, datatype, localfiles, redoverification, offline, logfile):
        """
        Inizialize the datastructure and log the parameters used.
        :param genome: list of genomes
        :param chromosome: list of chromosomes
        :param biosource: list of biosources
        :param epigenetic_mark: list of TFs
        :param output_path: path for where to store data
        :param csv_name: what to name the linking table
        :param datatype: what filetype to convert to
        :param localfiles: path to where localfiles are if used
        :param redoverification: boolean to check if force verification is done
        :param offline: boolean to check if offline mode is used
        :param logilfe: path to the logfile used in the run
        """
        self.genome = genome
        self.chromosome = chromosome
        self.biosource = biosource
        self.epigenetic_mark = epigenetic_mark
        self.binpath = os.path.abspath(os.path.dirname(__file__))
        self.outpath = output_path
        self.csvname = csv_name
        self.type = datatype
        self.localfiles = localfiles
        self.redoverification = redoverification
        self.offline = offline
        self.logfile = logfile

        logging.info("Genomes: " + '; '.join(self.genome))
        logging.info("Biosources: " + '; '.join(self.biosource))
        logging.info("Epigenetic marks: " + '; '.join(self.epigenetic_mark))
        logging.info("chromosomes: " + '; '.join(self.chromosome))

    def pull_data(self):
        """
        Recommended way to use this wrapper. Calls all needed functions.
        Depending on offlinemode, redoverification or if need data is used.
        """
        newdata = False
        if not self.offline:
            self.generate_csv()
            newdata = self.download_data()
        if newdata or self.redoverification:
            self.validate_convert_files()
            self.merge_forward_reverse()
            self.sort_files()
            self.generate_dictionaries()
        else:
            print(
                "No new data was downloaded, skipping validation, merging and sorting.")
        self.normalize()

    def generate_csv(self):
        """
        Ǵenerate a .csv containing all files fitting the parameters.

        Calls csv.r and handles the return value. csv.r uses genomes,
        biosoruces, chromosomes and epigenetic marks to generete a csv of data
        to be downloaded.

        This function also downloads the needed chrom.sizes files.

        """
        logging.info("starting csv creation")
        tool = os.path.join(self.binpath, "scripts",
                            "generate_linking_table.r")
        path = os.path.join(self.outpath, "data", "download")
        args = [tool, "-d", path, "-o", self.csvname, "-g"]
        args.extend(self.genome)
        args.append("-b")
        args.extend(self.biosource)
        args.append("-c")
        args.extend(self.chromosome)
        args.append("-m")
        args.extend(self.epigenetic_mark)
        args.append("--skip_verification")
        args.append("--append")
        args.append("--generate_chrom_sizes")
        rc = subprocess.call(args)
        if rc != 0 or not os.path.isfile(os.path.join(path, self.csvname)):
            logging.error("error generating .csv")
            raise Exception("csv.r could not create CSV")
        logging.info("finished csv creation")

    def download_data(self):
        """
        Download files from generated csv

        Calls export_from_csv.r and handles the return value.
        export_from_csv.r takes a csv generated by csv.r and downlaods the
        listed files into a given directory.

        If the returnvalue is 2 no new data was downloaded.
        """
        logging.info("starting data download")
        tool = os.path.join(self.binpath, "scripts",
                            "download_deepblue_data.r")
        csv = os.path.join(self.outpath, "data", "download", self.csvname)
        outdir = os.path.join(self.outpath, "data", "download")
        args = [tool, "-i", csv, "-o", outdir]
        args.append("-g")
        args.extend(self.genome)
        args.append("-b")
        args.extend(self.biosource)
        args.append("-c")
        args.extend(self.chromosome)
        args.append("-m")
        args.extend(self.epigenetic_mark)
        # check if localfiles are used and add parameter if they are
        if self.localfiles is not None:
            args.append("-l")
            args.append(self.localfiles)

        rc = subprocess.call(args)
        if rc == 2:
            logging.info("new new data was downloaded")
            logging.info("finished data download")
            return False
        if rc != 0:
            logging.error("export_from_csv.r could not download data")
            raise Exception("export_from_csv.r could not download data")
        logging.info("fXinished data download")
        return True

    def validate_convert_files(self):
        """
        Validates filetypes and converts them to requested filetype if needed.

        Before doing this the function goes through the linking table generated
        by generate_linking_table.r to create a list of all files needed to be
        validated for the given parameters

        Calls convert_files.sh and handles return value.
        convert_files.sh takes :
        -filetype: filetype to convert to
        -source_path: path to the downloaded data (./data/downlaod/)
        -out_path:path to put the  validated files (./data/temp/)
        -chrom_path: path to the chrom.sizes files
        -csv_name: name of the linking table file
        -logfile: path to the logfile
        """
        logging.info("starting file validation")
        tool = os.path.join(self.binpath, "scripts", "convert_files.sh")
        indir = os.path.join(self.outpath, "data", "download")
        outdir = os.path.join(self.outpath, "data")

        # create normalization.csv
        validation_csv = os.path.join(indir, "validation.csv")
        with open(validation_csv, 'w', newline='') as csvfile:
            outcsv = writer(csvfile, delimiter=';')
            found_biosources = []
            found_genomes = []
            found_chromosomes = []
            with open(os.path.join(indir, self.csvname), 'r') as csv:
                csv_reader = reader(csv, delimiter=';')
                header = next(csv_reader)
                outcsv.writerow(header)
                genome = header.index("genome")
                biosource = header.index("biosource")
                chromosome = header.index("chromosome")
                tf = header.index("epigenetic_mark")
                technique = header.index("technique")
                # find every chipseq that ifts the params and note
                # which combinations have been found.
                for row in csv_reader:
                    if (row[genome] in self.genome and
                            row[biosource] in self.biosource and
                            row[tf] in self.epigenetic_mark and
                            row[chromosome] in self.chromosome):

                        if row[genome] not in found_genomes:
                            found_genomes.append(row[genome])
                        if row[biosource] not in found_biosources:
                            found_biosources.append(row[biosource])
                        if row[chromosome] not in found_chromosomes:
                            found_chromosomes.append(row[chromosome])
                        outcsv.writerow(row)
                # go through the file again and check if atac/dnase
                # files fitting the parameters are found
                # reset the csv reader before doing so
                csv.seek(0)
                for row in csv_reader:
                    if row[technique] in ["atac-seq", "dnase-seq"]:
                        if (row[genome] in self.genome and
                                row[biosource] in self.biosource and
                                row[chromosome] in self.chromosome):
                            outcsv.writerow(row)

        # call convert_files.sh to actually validate and convert the files
        rc = subprocess.call(
            ["bash", tool, self.type, indir, outdir, self.csvname, self.logfile])
        if rc != 0:
            logging.error("convert_files.sh could not convert files")
            raise Exception("convert_files.sh could not convert files")

        # Keep backup of the last validation file for debugging purposes
        # os.rename does not overwrite on Windows meaning a check is needed.
        old_val = os.path.join(indir, "validaiton.csv.old")
        try:
            os.rename(validation_csv, old_val)
        except OSError:
            os.remove(old_val)
            os.rename(validation_csv, old_val)

        logging.info("finished file validation")

    def merge_forward_reverse(self):
        """
        merge forward/reverse read files into a single .bw
        creates a list of all the chrom.sizes files needed for the function
        Calls scripts.merge_files and handles return value

        """
        logging.info("starting forward reverse merging")
        bigwigMerge = subprocess.check_output(["which", "bigWigMerge"])
        bedgraphtobigwig = subprocess.check_output(
            ["which", "bedGraphToBigWig"])
        csvpath = os.path.join(self.outpath, "data", "temp", self.csvname)
        chromsizes = []
        for genome in self.genome:
            chromsizes.append(os.path.join(self.outpath, "data",
                                           "chromsizes", genome + ".chrom.sizes"))
        merge_all(csvpath, chromsizes, ["bigwig"],
                  bedgraphtobigwig, bigwigMerge)
        logging.info("finished forward reverse merging")

    def sort_files(self):
        """
        merge forward/reverse read files into a single .bw

        Calls scripts.merge_files and handles return value

        """
        logging.info("starting Filesorting")
        tool = os.path.join(self.binpath, "scripts", "sort_files.sh")
        outdir = os.path.join(self.outpath, "data")
        indir = os.path.join(outdir, "temp")
        csv = os.path.join(indir, self.csvname)

        rc = subprocess.call(
            ["bash", tool, indir, outdir, csv, self.csvname])
        if rc != 0:
            logging.error("sort_files.sh could not sort files")
            raise Exception("sort_files.sh could not sort files")
        logging.info("finished Filesorting")

    def normalize(self):
        """
        Normalize files to allow proper analysis
        Creates a normalize.csv and then calls scripts.normalize and handles
        the return value
        """
        logging.info("Creating normalization.csv")

        # create normalization.csv
        data = os.path.join(self.outpath, "data")
        norm_csv = os.path.join(data, "temp", "normalization.csv")
        set_header = False
        with open(norm_csv, 'w', newline='') as csvfile:
            outcsv = writer(csvfile, delimiter=';')
            found_biosources = []
            found_genomes = []
            found_chromosomes = []
            # first run goes through all files and appends ever chipseq it find
            for file in os.listdir(data):
                if file.endswith(".csv"):
                    with open(os.path.join(data, file), 'r') as csv:
                        csv_reader = reader(csv, delimiter=';')
                        header = next(csv_reader)
                        if not set_header:
                            outcsv.writerow(header)
                            set_header = True
                        genome = header.index("genome")
                        biosource = header.index("biosource")
                        chromosome = header.index("chromosome")
                        filename = header.index("filename")
                        tf = header.index("epigenetic_mark")
                        technique = header.index("technique")
                        # find every chipseq that ifts the params and note
                        # which combinations have been found.
                        for row in csv_reader:
                            if (row[genome] in self.genome and
                                    row[biosource] in self.biosource and
                                    row[tf] in self.epigenetic_mark and
                                    row[chromosome] in self.chromosome and
                                    row[filename].endswith(".bw")):

                                if row[genome] not in found_genomes:
                                    found_genomes.append(row[genome])
                                if row[biosource] not in found_biosources:
                                    found_biosources.append(row[biosource])
                                if row[chromosome] not in found_chromosomes:
                                    found_chromosomes.append(row[chromosome])
                                outcsv.writerow(row)
                        # go through the file again and check if atac/dnase
                        # files fitting are found
                        # reset the csv reader before doing so
                        csv.seek(0)
                        for row in csv_reader:
                            if row[technique] in ["atac-seq", "dnase-seq"]:
                                if (row[genome] in self.genome and
                                        row[biosource] in self.biosource and
                                        row[chromosome] in self.chromosome and
                                        row[filename].endswith(".bw")):
                                    outcsv.writerow(row)

        logging.info("starting Normalization")
        normalize_all(norm_csv)
        logging.info("finished Normalization")
        old_norm = os.path.join(data, "temp", "normalization.csv.old")

        # Keep backup of the last normalization file for debugging purposes
        # os.rename does not overwrite on Windows meaning a check is needed.
        try:
            os.rename(norm_csv, old_norm)
        except OSError:
            os.remove(old_norm)
            os.rename(norm_csv, old_norm)

        logging.info(
            f"cleaned up normalization, the run can be found at {old_norm}")

    def generate_dictionaries(self):
        """
        Generate pickle files for the downloaded data.
        Calls generate_pickle.py with the path the data is stored in.
        """
        path = os.path.join(self.outpath, "data")
        parse(path)
