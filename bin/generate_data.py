import subprocess
import os
import logging
from csv import reader, writer
from datetime import datetime
from scripts.merge_reads import merge_all
from scripts.generate_pickle import parse
from scripts.normalize_signal_values import normalize_all
# //TODO: angular insall in visualisation folder


class DataConfig:
    """Contains configuration data for the pull_data function.
    - chromsizes = path to folder of chromsizes
    - type = type of data (signal/peak)
      """

    def __init__(self, genome, chromosome, biosource, epigenetic_mark,
                 output_path, csv_name, datatype, localfiles, redoverification, offline):
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
        self.logfile = self.setup()

        logging.info(self.genome + self.biosource + self.epigenetic_mark)

    def setup(self):
        """
        Sets up logging and Ensures proper filestructure is given.
        """
        time = datetime.now().strftime("%d_%m_%Y_%H_%M_%S")
        temp = os.path.join(self.outpath, "data", "temp")
        result = os.path.join(self.outpath, "results")
        logs = os.path.join(self.outpath, "logs")
        download = os.path.join(self.outpath, "data", "download")
        chromsizes = os.path.join(self.outpath,
                                  "data", "chromsizes")
        if not os.path.exists(download):
            os.makedirs(download)
        if not os.path.exists(temp):
            os.makedirs(temp)
        if not os.path.exists(result):
            os.makedirs(result)
        if not os.path.exists(logs):
            os.makedirs(logs)
        if not os.path.exists(chromsizes):
            os.makedirs(chromsizes)

        logname = time + "_generate_data.log"
        logfile = os.path.join(logs, logname)
        logging.basicConfig(filename=logfile, level=logging.INFO)
        return logfile

    def pull_data(self):
        """ Recommended way to use this wrapper. Calls all needed functions.
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
        """Ǵenerate a .csv containing all files fitting the parameters.

        Calls csv.r and handles the return value. csv.r uses genomes,
        biosoruces and epigenetic marks to generete a csv of data to be
        downloaded.
        """
        tool = os.path.join(self.binpath, "scripts", "csv.r")
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

    def download_data(self):
        """Download files from .csv

        Calls export_from_csv.r and handles the return value.
        export_from_csv.r takes a csv generated by csv.r and downlaods the listed files into a given directory.

        """
        tool = os.path.join(self.binpath, "scripts", "export_from_csv.r")
        csv = os.path.join(self.outpath, "data", "download", self.csvname)
        outdir = os.path.join(self.outpath, "data", "download")
        if self.localfiles is not None:
            rc = subprocess.call(
                [tool, "-i", csv, "-o", outdir, "-l", self.localfiles])
        else:
            rc = subprocess.call(
                [tool, "-i", csv, "-o", outdir])
        if rc == 2:
            logging.info("new new data was downloaded")
            return False
        if rc != 0:
            logging.error("export_from_csv.r could not download data")
            raise Exception("export_from_csv.r could not download data")
        return True

    def validate_convert_files(self):
        """ Validates filetypes and converts them to requested filetype if needed

        Calls convert_files.sh and handles return value
        convert_files.sh takes a fileformat and a path with the
        convertable files
        """
        # TODO: outpath muss whitespace escapen/ mit "" übergeben werden
        tool = os.path.join(self.binpath, "scripts", "convert_files.sh")
        indir = os.path.join(self.outpath, "data", "download")
        outdir = os.path.join(self.outpath, "data")

        rc = subprocess.call(
            ["bash", tool, "bigwig", indir, outdir, self.csvname])
        if rc != 0:
            logging.error("convert_files.sh could not convert files")
            raise Exception("convert_files.sh could not convert files")

    def merge_forward_reverse(self):
        """ merge forward/reverse read files into a single .bw

        Calls scripts.merge_files and handles return value

        """
        bigwigMerge = subprocess.check_output(["which", "bigWigMerge"])
        bedgraphtobigwig = subprocess.check_output(
            ["which", "bedGraphToBigWig"])
        csvpath = os.path.join(self.outpath, "data", "temp", self.csvname)
        chromsizes = []
        for genome in self.genome:
            chromsizes.append(os.path.join(self.outpath, "data",
                                           "chromsizes", genome + ".chrom.sizes"))
        # TODO: self.chromsizes replace with array that contains paths to direct files
        merge_all(csvpath, chromsizes, ["bigwig"],
                  bedgraphtobigwig, bigwigMerge)

    def sort_files(self):
        """ merge forward/reverse read files into a single .bw

        Calls scripts.merge_files and handles return value

        """
        tool = os.path.join(self.binpath, "scripts", "sort_files.sh")
        outdir = os.path.join(self.outpath, "data")
        indir = os.path.join(outdir, "temp")
        csv = os.path.join(indir, self.csvname)

        rc = subprocess.call(
            ["bash", tool, indir, outdir, csv, self.csvname])
        if rc != 0:
            logging.error("sort_files.sh could not sort files")
            raise Exception("sort_files.sh could not sort files")

    def normalize(self):
        """ Normalize files to allow proper analysis
        Creates a normalize.csv and then
        Calls scripts.normalize and handles return value

        """
        data = os.path.join(self.outpath, "data")
        norm_csv = os.path.join(data, "temp", "normalization.csv")
        set_header = False
        with open(norm_csv, 'w', newline='') as csvfile:
            outcsv = writer(csvfile, delimiter=';')
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
                        for row in csv_reader:
                            if (row[genome] in self.genome and
                                row[biosource] in self.biosource and
                                    row[chromosome] in self.chromosome and
                                    row[filename].endswith(".bw")):
                                outcsv.writerow(row)
        normalize_all(norm_csv)

    def generate_dictionaries(self):
        """Generate pickle files for the downloaded data """
        path = os.path.join(self.outpath, "data")
        parse(path)
