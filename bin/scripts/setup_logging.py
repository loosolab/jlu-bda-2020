from datetime import datetime
import logging
import os


def setup(outpath):
    """
    Sets up logging and Ensures proper filestructure is given.
    """
    time = datetime.now().strftime("%d_%m_%Y_%H_%M_%S")
    temp = os.path.join(outpath, "data", "temp")
    result = os.path.join(outpath, "results")
    logs = os.path.join(outpath, "logs")
    download = os.path.join(outpath, "data", "download")
    chromsizes = os.path.join(outpath,
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

    logname = time + "_tfanalyzer.log"
    logfile = os.path.join(logs, logname)
    logging.basicConfig(filename=logfile, level=logging.INFO)
    return logfile
