# TF Analyzer

This tool allows for the analysis and classification of transcription factors by calculating the relations between TF binding sites and chromatin accessibility, using the [Deepblue Epigenomic Data Server](https://deepblue.mpi-inf.mpg.de/) to link epigenomic ChIP-seq and ATAC/DNAse-seq data by biosource (i.e. cell type).

The latest working build can be cloned from the `dev` branch.

**Workflow**

 1. Link and download epigenomic data matching user-requested biosources, TFs, genomes and chromosomes
 2. Convert, normalize and sort downloaded data if necessary
 3. Calculate relations between ChIP and ATAC data for each TF and save results as CSV files
 4. Visualize data
 
## Requirements

 - Python 3 and pip
 - [Conda](https://docs.conda.io/projects/conda/en/latest/index.html)
 - [npm](https://www.npmjs.com/)
 
## Setup
1. Run `conda env create -f environment.yml` to create a conda environment with the following modules:
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
```

2. Run `conda activate tf_analyzer` to activate the environment. This needs to stay activated in order for the tool to work.

3. From the repository root, navigate to `visualization/` and execute the following command: `npm install -g @angular/cli`.

## How to use
The pipeline can be initiated by running `python bin/tf_analyzer.py` with these arguments:

`-g, --genome` genome (default: hg19)

`-b, --biosource` list of biosources (default: all)

`-t, --tf` list of transcription factors (default: all)

`-c, --chromosome` list of chromosomes (default: all chromosomes of the selected genome)

`-w, --width` single integer to determine the range of analysis (peak +- w, default: 50)

`-r, --redo_analysis` (optional) existing results will be executed again and overwritten

`-cs, --component_size` (optional) single integer determining the component size for the analysis (will be calculated if not given)

`-o, --output_path` output directory for all data (default: `./`)

The following (optional) arguments will not initiate the pipeline but display information gathered from already existing results:

`--visualize` visualize existing results

`--list_chromosomes` print a list of all genomes with their associated chromosomes

`--list_downloaded_data` print a table containing all downloaded genomes, biosources, tfs and chromosomes

"all" refers to all data that is currently available on the Deepblue server.

Use the help argument (`-h` or `--help`) to display a more detailed list of available arguments.

The results (plot images and CSV files) can be found in the `results` folder inside the output directory. To view the plots in detail, a web application is available at `http://localhost:4200/`.
 
## Example case
`python bin/tf_analyzer.py -g hg19 -b kidney -t AR -c chr1`

This command will download and analyze all data for transcription factor "AR", a known activator, sampled from chromosome 1 of the biosource "kidney".

**Results**

For each transcription factor analyzed in the process, a set of three figures will be exported:

<img src="docs/img/rela1.png" width="500">

Density Scatter Plot: This scatter plot contains the means of the ATAC scores on the x-axis and the ChIP scores on the y-axis. The histograms on both axes are showing the distribution of respective ATAC/ChIP scores. This plot also contains a heatmap, which shows the densities of the values used. (Zooming in may be required to see the details.)

<img src="docs/img/rela4.png" width="500">

3d Contour Plot: This 3D contour plot shows the means of the ATAC scores on the x-axis, the ChIP scores on the y-axis and the density of these values on the z-axis.

The third output is a table containing the weights of the individual components of the analysis.

If any problems occur in the visualization, please check the web console of your browser.

## Known errors

### "Cannot find module..." when launching the web server
```
An unhandled exception occured: Cannot find module '@angular-devkit/build-angular/package.json'
Require stack:
[...]
```
To fix this, run `npm install --save-dev @angular-devkit/build-angular` inside the `visualization/` folder. Afterwards, run `ng serve` (in the same folder) to start the server.

## License
(to be added)
