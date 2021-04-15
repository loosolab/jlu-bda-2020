# TF Analyzer

## Introduction

Transcription factors are well-known for their important role in regulating gene expression. Modern methods such as **ChIP-seq** (Chromatin Immunoprecipitation sequencing) allow for high-resolution detection of their DNA binding sites. Histone modification is another identifier of regulatory activity on the genome, and techniques such as **DNase-seq** (DNase I hypersensitive sites sequencing) and **ATAC-seq** (Assay for Transposase-Accessible Chromatin using sequencing) provide genome-wide information on DNA accessibility. Combining the observed binding sites of known transcription factors with the degree of DNA accessibility within the same regions of the same sample may reveal new insights on their function as transcriptional activators, repressors, or both.

To this end, the TF Analyzer utilizes the [Deepblue Epigenomic Data Server](https://deepblue.mpi-inf.mpg.de/), which currently stores close to 40,000 ChIP-seq and roughly 3,800 DNase-seq experiments, to link corresponding DNA binding site and accessiblity experiments by their shared biosources (i.e. cell types). The downloaded files are normalized by means of both logarithmic and min/max scaling and individual scores are calculated for each transcription factor. The results are displayed as charts through a user-friendly web interface (see [Example case](#example-case)).

The latest working build can be cloned from the `dev` branch.

For further information on how this tool works, please check out our [wiki](https://github.com/loosolab/jlu-bda-2020/wiki/).
If you encounter any errors during a run, please refer to the "[Problems and Solutions](https://github.com/loosolab/jlu-bda-2020/wiki/Problems-and-Solutions)" section of the wiki before opening an Issue.

**Workflow**

 1. Link and download epigenomic data matching user-requested biosources, TFs, genomes and chromosomes
 2. Convert, normalize and sort downloaded data if necessary
 3. Calculate relations between ChIP and ATAC/DNase data for each TF and save results as CSV files
 4. Visualize data through a self-hosted web application

## Requirements

 - Python 3 and pip
 - [Conda](https://docs.conda.io/projects/conda/en/latest/index.html)

## Setup
First, we need to create a conda environment with the following dependencies:

**environment.yml**

```
channels:
  - bioconda
  - conda-forge
dependencies:
  - r-base=4.0.3
  - r-data.table=1.14.0
  - bioconductor-deepbluer=1.16.0
  - ucsc-bedgraphtobigwig=377
  - ucsc-bigwigtobedgraph=377
  - ucsc-bigwigmerge=377
  - pandas=1.2.3
  - nodejs=15.12.0
  - flask=1.1.2
  - flask-restful=0.3.8
  - flask-cors=3.0.8
  - bash=5.0.018
  - pybigwig=0.3.18
  - matplotlib=3.3.4
  - kneed=0.7.0
  - sklearn2pmml=0.69.0
  - natsort=7.1.1
  - tabulate=0.8.9
  - util-linux=2.36
  - psutil=5.8.0
  - regex=2021.4.4
```

```bash
$ conda env create -f environment.yml
```

The newly created environment must then be activated before running the pipeline:

```bash
$ conda activate tf_analyzer
```

Finally, from the repository root, navigate to `visualization/` and execute the following commands:

```bash
$ npm install -g @angular/cli

$ npm install --save-dev @angular-devkit/build-angular
```

**Note:** To update an existing environment, run

```bash
$ conda env update -f environment.yml
```

## How to use
The pipeline can be initiated by running
```bash
$ python bin/tf_analyzer.py
```

with these arguments:

`-g, --genome` genome (default: hg19)

`-b, --biosource` list of biosources (default: all)

`-t, --tf` list of transcription factors (default: all)

`-c, --chromosome` list of chromosomes (default: all chromosomes of the selected genome)

`-w, --width` single integer to determine the range of analysis on the chromosomes (peak +- w, default: 50)

`-cs, --component_size` single integer determining the component size for the analysis (will be calculated if not given)

`-o, --output_path` output directory for all data (default: `./`)

`--offline` runs the program in offline mode (will not query or download data from Deepblue)

`--redo_file_validation` downloaded files will be validated and sorted again (even if no new files were downloaded)

`--check_local_files` external directory containing Deepblue data (any files in this directory that would have to be downloaded will be copied from here into the output directory)

The following (optional) arguments will not initiate the pipeline but display information gathered from already existing results:

`--visualize` visualize existing results

`--list_chromosomes` print a list of all genomes with their associated chromosomes

`--list_downloaded_data` print a table containing all downloaded genomes, biosources, tfs and chromosomes

"all" refers to all data that is currently available on the Deepblue server.

Use the help argument (`-h` or `--help`) to display a more detailed list of available arguments.

**Note:** Please note that all arguments must be written in *lower case* letters. Multiple arguments (where applicable) must be divided by spaces (e.g. `-b "kidney" "huvec cell" "regulatory t cell"`).

The results can be found in the `plots` folder inside the output directory. To view the plots in detail, a web application is available at `http://localhost:4200/`. If the application is launched through a virtual machine, it can be accessed locally via SSH:

```bash
$ ssh -L 4200:localhost:4200 -L 5000:localhost:5000 user@server
```

## Example case

```bash
$ python bin/tf_analyzer.py -g hg19 -b kidney -t ar -c chr20
```

This command will download and analyze all data for transcription factor "AR", a known activator, sampled from chromosome 20 of the biosource "kidney".

**Results**

The web application running on `localhost:4200` lets you select the transcription factors from your last run that you want to look at. Three figures are currently available for each TF:

![ar_20_Contour](https://user-images.githubusercontent.com/26332337/114243219-7f96ea00-998c-11eb-9de7-26adf0dc3483.png)

3d Contour Plot: Displays the mean ATAC/DNase scores on the x-axis, the ChIP scores on the y-axis and the density of these values on the z-axis.

![ar_20_Scatter](https://user-images.githubusercontent.com/26332337/114243247-8aea1580-998c-11eb-80a6-bf1589a25f46.png)

Density Scatter Plot: Shows the mean ATAC/DNase scores on the x-axis and the ChIP scores on the y-axis. The histograms on both axes are showing the distribution of respective ATAC/ChIP scores. This plot also contains a heatmap, which shows the densities of the values used. (Zooming in may be required to see the details.)

The third output is a table containing the weights of the individual components of the analysis.

The results shown above are highly characteristic of a transcriptional activator due to a high degree of DNA accessibility and a high detection rate of DNA binding sites within the same genomic regions.

## License
This project is licensed under the [MIT license](https://github.com/loosolab/jlu-bda-2020/blob/main/LICENSE).
