#!/usr/bin/env Rscript
required_packages <- c("BiocManager","data.table","argparse")
not_installed <- setdiff(required_packages,rownames(installed.packages()))
if(length(not_installed) > 0) install.packages(not_installed,repos="http://cran.us.r-project.org")
if(!"DeepBlueR" %in% installed.packages()) BiocManager::install("DeepBlueR")
library(argparse)

parser <- ArgumentParser()

parser$add_argument("-g", "--genome", type="character", default="hg19",
                    help="Genome to search Deepblue database with [default: \"%(default)s\"]")
parser$add_argument("-c", "--chromosomes", nargs="+", type="character", default=NULL,
                    help="(List of) chromosomes to include (requires chr prefix) [default: all]")
parser$add_argument("-b", "--biosources", nargs="+", type="character", default=NULL,
                    help="(List of) biosources to include [default: all] (Refer to: https://deepblue.mpi-inf.mpg.de/)")
parser$add_argument("-t", "--type", nargs="+", type="character", default="peaks", choices=c("peaks","signal"),
                    help="Experiment file types allowed for CHiP-Seq data [default: \"%(default)s\"]")
parser$add_argument("-a", "--atactype", nargs="+", type="character", default="signal", choices=c("peaks","signal"),
                    help="Experiment file types allowed for ATAC/DNAse-Seq data [default: \"%(default)s\"]")
parser$add_argument("-m", "--marks", nargs="+", type="character", default=NULL,
                    help="(List of) epigenetic marks (i.e. transcription factors) to include [default: all] (Refer to: https://deepblue.mpi-inf.mpg.de/)")
parser$add_argument("-d", "--directory", type="character", default=".",
                    help="Output directory for CSV file. If an output filename (-o) containing \"/\" characters is provided, the filename will be used instead. [default: \"%(default)s\"]")
parser$add_argument("-o", "--output", type="character", default="linking_table",
                    help="Output file name without extension [default: \"%(default)s\"]")
parser$add_argument("--append", action="store_true",
                    help="Results are compared and, if needed, appended to an existing CSV file with the same filename (see -o).")
parser$add_argument("--skip_verification", action="store_true",
                    help="Input genomes, chromosomes, biosources and epigenetic marks will not be checked against the Deepblue database. Missing arguments will be treated as NULL.")
parser$add_argument("--generate_chrom_sizes", dest="chromsizes", action="store_true",
                    help="Generates a tab-separated chrom.sizes txt file for the requested genome.")

args <- parser$parse_args()

library(data.table)

suppressPackageStartupMessages(
  library(DeepBlueR)
)

# create_linking_table(args)
#
# The main function of this script takes the command-line arguments (if
# provided by the user) and generates a CSV with corresponding data from the
# Deepblue Epigenomic Server. The ChIP and ATAC/DNAse-seq data are linked by
# matching their biosources. The general workflow is like this:
# 0. If requested, check the command-line arguments against Deepblue database
# 1. Query ATAC/DNAse-seq experiments corresponding to genome and biosources
# 2. Query ChIP-seq experiments corresponding to remaining biosources and
#    transcription factors
# 3. Filter ATAC/DNAse data by remaining biosources from step 2 and build a list
#    of data.tables for each experiment and chromosome.
# Finally: Export the data.table as CSV.

create_linking_table <- function(genome,chrs,filter_biosources,chip_type,atac_type,chip_marks,outdir,outfile,append,skip,chromsizes) {
  
  # new_row(metadata)
  #
  # Accepts a list containing metadata of one experiment, which is the standard
  # output of deepblue_info(experiment_id). In the final CSV file, the first few
  # columns are fixed, so these are extracted explicitly with some added info
  # (sample info & extra metadata). Since all downloaded files are split into
  # chromosomes, this function uses expand.grid to create data.tables of the
  # same file for all desired chromosomes.
  # Example: If the user requested chr1, chr2 and chr3, then each experiment
  # found in the Deepblue database will add three rows of metadata to the list
  # of data.tables which new_row() returns. Only the `filename` and `chromosome`
  # columns will differ between the three rows.
  
  new_row <- function(metadata) {
    
    name <- metadata$name
    message(name)
    
    return_list <- apply(expand.grid(name,chrs),1,function(x) {
      
      dots <- strsplit(x[1],".",fixed=TRUE)[[1]]
      mylen <- length(dots)
      if(mylen > 1) {
        tmp <- dots[mylen]
        dots[mylen] <- x[2]
        dots[mylen+1] <- tmp
        filename <- paste(dots,collapse=".")
      } else {
        filename <- paste(x,collapse=".")
      }
      
      e <- metadata$`_id`
      sample_info <- metadata$sample_info
      extra <- metadata$extra_metadata
      meta <- data.table(
        experiment_id=e,
        genome=metadata$genome,
        biosource=tolower(metadata$sample_info$biosource_name),
        technique=tolower(metadata$technique),
        epigenetic_mark=tolower(metadata$epigenetic_mark),
        chromosome=x[2],
        filename,
        data_type=metadata$data_type,
        format=metadata$format,
        sample_id=metadata$sample_id,
        project=metadata$project,
        total_size=metadata$upload_info$total_size
        # as.data.table(sample_info),
        # as.data.table(extra)
      )
      return(meta)
      
    })
    
    return(return_list)
  }
  
  # verify_filters(values in query, values to check against, type of values)
  #
  # Automates the verification process of user input values (i.e. biosources,
  # TFs) which are checked against the data Deepblue provides. Both input_values
  # and all_values are character vectors while "type" is a string which is used
  # in case of an error to specify which type of data it originated from.
  # All values are converted to lower case to simplify the comparison. If all
  # user input values are invalid, the script will stop the execution.
  #
  # If the input values are NULL, this is likely due to missing arguments in
  # argparse, which is interpreted as if the user is requesting all available
  # data.
  
  verify_filters <- function(input_values,all_values,type="values") {
    message(paste("verifying",type,"..."))
    if(is.null(input_values)) {
      return(all_values)
    }
    input_values <- tolower(input_values)
    all_values <- tolower(all_values)
    removed_values <- input_values[!input_values %in% all_values]
    n_removed <- length(removed_values)
    if(n_removed > 0) warning(paste("dropped",n_removed,type,":",paste(removed_values,collapse=", ")))
    filter_values <- input_values[!input_values %in% removed_values]
    if(length(filter_values) == 0) {
      stop(paste("No",type,"given by user could be found in Deepblue database"))
    } else {
      return(filter_values)
    }
  }
  
  all_chroms <- suppressMessages(deepblue_chromosomes(genome = genome))
  
  if(!skip) {
    all_genomes <- suppressMessages(tolower(deepblue_list_genomes()$name))
    if(!genome %in% all_genomes) {
      stop(paste("No valid genomes provided by user. Available genomes:",paste(all_genomes,collapse=", ")))
    }
    
    chrs <- verify_filters(chrs,all_chroms$id,"chromosomes")
    
    all_biosources <- suppressMessages(deepblue_list_biosources()$name)
    filter_biosources <- verify_filters(filter_biosources,all_biosources,"biosources")
    
    tf_marks <- suppressMessages(deepblue_list_epigenetic_marks(extra_metadata = list(category="Transcription Factor Binding Sites"))$name)
    chip_marks <- verify_filters(chip_marks,tf_marks,"epigenetic marks")
  } else {
    # If --skip_verification is selected and...
    #    -g is missing       : defaults to hg19
    #    -b is missing (NULL): all biosources are used
    #    -m is missing (NULL): all epigenetic marks (including non-TFs) are used
    #    -c is missing (NULL): chrs stays NULL and must be assigned separately
    #                          (see below)
    if(is.null(chrs)) {
      chrs <- all_chroms$id
    }
  }
  
  if(chromsizes) {
    # Unrelated to the CSV generation process, this section creates a tab-
    # separated text file containing the sizes of all chromosomes of the current
    # genome as "{genome_name}.chrom.sizes".
    sizefile <- paste(genome,"chrom","sizes",sep=".")
    size_output <- paste(outdir,sizefile,sep="/")
    if(!file.exists(size_output)) {
      write.table(all_chroms,file=size_output,col.names=F,row.names=F,quote=F,sep="\t")
      message(paste(sizefile,"written to",normalizePath(outdir)))
    }
  }
    
  # 1st Step: Collect biosources for ATAC
  
  atac <- suppressMessages(deepblue_list_experiments(genome=genome, technique="ATAC-Seq", biosource=filter_biosources, type=atac_type))
  dnase <- suppressMessages(deepblue_list_experiments(genome=genome, technique="DNAse-Seq", biosource=filter_biosources, type=atac_type))
  
  # deepblue_list_experiments(...) returns "\n" if no experiments match the
  # arguments, a chr list of 2 if one experiment is found and otherwise a
  # data.frame of all results. is.list() treats data.frames and data.tables as
  # lists as well.
  
  if(is.list(atac)) atac <- atac$id else atac <- NULL
  if(is.list(dnase)) dnase <- dnase$id else dnase <- NULL
  
  if(is.null(atac) && is.null(dnase)) {
    stop(paste(genome,"No ATAC/DNAse-seq data available for given arguments",sep=": "))
  }
  access_experiments <- c(atac,dnase)
  message(paste("fetched",length(access_experiments),"ATAC/DNAse-seq experiments"))
  
  atac_metadata <- lapply(access_experiments,function(x){ suppressMessages(deepblue_info(x)) })
  
  # The biosources are extracted individually to create a filter for the next
  # step, since we only need ChIP experiments matching ATAC/DNAse biosources.
  atac_biosources <- unique(vapply(atac_metadata,function(x) { tolower(x$sample_info$biosource_name) },character(1L)))
  
  # 2nd Step: Collect ChiP experiments and add to csv list
  
  chips <- suppressMessages(deepblue_list_experiments(genome=genome, technique="ChIP-Seq", biosource=atac_biosources, type=chip_type, epigenetic_mark=chip_marks))

  if(!is.list(chips)) { 
  
      stop(paste(genome,"No ChIP-seq data available for given arguments",sep=": "))
    
  } else {
    
    chips <- chips$id
    message(paste("fetched",length(chips),"ChIP-seq experiments"))
    
    chip_metadata <- lapply(chips,function(c) {
      metadata <- suppressMessages(deepblue_info(c))
    })
    
    chip_biosources <- unique(tolower(vapply(chip_metadata,function(x){ x$sample_info$biosource_name },character(1L))))
    chip_list <- lapply(chip_metadata,new_row)
    
    # before: List of lists of data.tables
    chip_list <- unlist(chip_list,recursive=FALSE)
    # after:  List of data.tables
    
    # 3rd Step: Add ATAC experiments to csv if biosource has available ChiP data
    # Likewise, we only need ATAC/DNAse data that matches biosources retained
    # by our ChIP-seq experiments.
    
    atac_metadata <- atac_metadata[lapply(atac_metadata,function(x){ tolower(x$sample_info$biosource_name) }) %in% chip_biosources]
    atac_list <- unlist(lapply(atac_metadata,new_row),recursive=FALSE)
    
    csv_data <- rbindlist(c(chip_list,atac_list),fill=TRUE)
    
  }

  # If "output" argument contains "/", ignore "directory" argument
  # Otherwise, -o and -i arguments are pasted into the full path
  
  if(grepl("/",outfile)) {
    splitfile <- strsplit(outfile,"/")[[1]]
    outdir <- paste(splitfile[-length(splitfile)],collapse="/")
    output_file <- outfile
  } else {
    if(length(outdir) > 0) {
      outdir <- paste(strsplit(outdir,"/")[[1]],collapse="/") # remove terminating "/" characters
      output_file <- paste(paste(outdir,outfile,sep="/"))
    } else {
      output_file <- outfile
    }
  }
  output_file <- paste(sub("\\.csv$","",output_file),"csv",sep=".")
  
  if(is.null(csv_data)) {
    
    stop("No data available for CSV")
    
  } else {
    
    # check whether CSV with given filename already exists - if yes, add new rows
    if(append && file.exists(output_file)) {
      
      old_csv <- fread(file=output_file,header=TRUE,sep=";",colClasses=c("character"))
      old_filenames <- old_csv$filename
      csv_data.unique <- csv_data[!csv_data$filename %in% old_filenames]
      
      if(nrow(csv_data.unique) > 0) {
        
        csv_data <- rbind(old_csv,csv_data.unique,fill=TRUE)
        fwrite(csv_data,file=output_file,na="",sep=";")
        message(paste(nrow(csv_data.unique),"lines added to",normalizePath(output_file)))
        
      } else {
        
        message(paste("No new data was added to",normalizePath(output_file)))
        
      }
      
    } else {
      
      if(!dir.exists(outdir)) {
        dir.create(outdir,recursive = TRUE)
      }
      
      fwrite(csv_data,file=output_file,na="",sep=";")
      message(paste(nrow(csv_data),"lines written to",normalizePath(output_file)))
      
    }
  }
  
}

create_linking_table(args$genome,args$chromosomes,args$biosources,args$type,args$atactype,args$marks,args$directory,args$output,args$append,args$skip_verification,args$chromsizes)
