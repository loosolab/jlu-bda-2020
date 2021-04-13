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
  
  # new_row(metadata, output filename, chromosomes)
  #
  # Accepts a list containing metadata of one experiment, which is the standard
  # output of deepblue_info(experiment_id). Since all downloaded files are split into
  # chromosomes, this function uses expand.grid to create data.tables of the
  # same file for all desired chromosomes.
  # Example: If the user requested chr1, chr2 and chr3, then each experiment
  # found in the Deepblue database will add three rows of metadata to the list
  # of data.tables which new_row() returns. Only the `filename` and `chromosome`
  # columns will differ between the three rows.
  #
  # For each experiment (and chromosome), the number of regions available on
  # Deepblue is checked beforehand to avoid including non-existent data in the
  # CSV.
  # The resulting list of data.tables is rbound into a single data.table and
  # appended to the CSV file.
  
  new_row <- function(metadata,output_file,chrs) {
    
    name <- metadata$name
    message(name)
    
    return_list <- apply(expand.grid(name,chrs),1,function(x) {
      
      filename <- unname(x[1])
      chr <- unname(x[2])
      
      query_id <- suppressMessages(deepblue_select_experiments(experiment_name=filename, chromosome=chr))
      request_id <- suppressMessages(deepblue_count_regions(query_id))
      region_len <- as.integer(suppressMessages(deepblue_download_request_data(request_id, do_not_cache=TRUE)))
      
      if(region_len > 0) {
        
        dots <- strsplit(filename,".",fixed=TRUE)[[1]]
        mylen <- length(dots)
        if(mylen > 1) {
          tmp <- dots[mylen]
          dots[mylen] <- chr
          dots[mylen+1] <- tmp
          filename <- paste(dots,collapse=".")
        } else {
          filename <- paste(x,collapse=".")
        }
        
        meta <- data.table(
          experiment_id=metadata$`_id`,
          genome=metadata$genome,
          biosource=tolower(metadata$sample_info$biosource_name),
          technique=tolower(metadata$technique),
          epigenetic_mark=tolower(metadata$epigenetic_mark),
          chromosome=chr,
          filename,
          data_type=metadata$data_type,
          format=metadata$format,
          sample_id=metadata$sample_id,
          project=metadata$project,
          regions=region_len
        )
        return(meta)
        
      }
      
    })
    
    fwrite(rbindlist(return_list),file=output_file,na="",sep=";",append = TRUE)
    
    return(length(return_list[lengths(return_list) != 0]))
    
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
    if(is.null(input_values)) {
      return(all_values)
    }
    message(paste("verifying",type,"..."))
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
    all_genomes <- suppressMessages(deepblue_list_genomes()$name)
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
  
  if(grepl("/",outfile)) {
    # This section handles the output directory (-d) and filename (-o) args.
    # If -o (outfile) contains "/", ignore -d (outdir) argument.
    splitfile <- strsplit(outfile,"/")[[1]]
    outdir <- paste(splitfile[-length(splitfile)],collapse="/")
    output_file <- outfile
  } else {
    # Otherwise, -d and -o arguments are pasted into the full path.
    if(length(outdir) > 0) {
      outdir <- sub("/$","",outdir) # remove terminating "/" characters
      output_file <- paste(paste(outdir,outfile,sep="/"))
    } else {
      output_file <- outfile
    }
  }
  output_file <- paste(sub("\\.csv$","",output_file),"csv",sep=".")

  # --append: If a CSV file with the -o filename already exists, only new
  #           results will be appended.
  # With append set to TRUE, the id and chromosome columns of the existing CSV
  # are imported to check later whether a combination of a certain experiment
  # and chromosome is already contained in the old CSV.
  # If append is set to FALSE, the file will be removed. This is because
  # new_row() always appends its lines regardless. If file.create() was used
  # instead of .remove, fwrite(..., append=TRUE) would not write a header.
  
  if(file.exists(output_file)) {
    if(append) {
      old_data <- fread(file=output_file,header=TRUE,sep=";",select=c("experiment_id","chromosome"))
    } else {
      file.remove(output_file)
    }
  } else {
    append <- FALSE
    if(!dir.exists(outdir)) {
      dir.create(outdir,recursive = TRUE)
    }
  }

  # 1st Step: Collect biosources for ATAC
  
  atac_experiments <- suppressMessages(deepblue_list_experiments(genome=genome, technique=c("ATAC-seq","DNAse-seq"), biosource=filter_biosources, type=atac_type))
  
  # deepblue_list_experiments(...) returns "\n" if no experiments match the
  # arguments, a chr list of 2 if one experiment is found and otherwise a
  # data.frame of all results. is.list() treats data.frames and data.tables as
  # lists as well.
  
  if(is.list(atac_experiments)) {
    atac_experiments <- atac_experiments$id
  } else {
    stop(paste(genome,"No ATAC/DNAse-seq data available for given arguments",sep=": "))
  }
  message(paste("fetching",length(atac_experiments),"ATAC/DNAse-seq experiments ..."))
  
  atac_metadata <- lapply(atac_experiments,function(x){
    metadata <- suppressMessages(deepblue_info(x))
    return(metadata[!names(metadata) %in% c("extra_metadata","upload_info","columns")])
  })
  
  # The biosources are extracted individually to create a filter for the next
  # step, since we only need ChIP experiments matching ATAC/DNAse biosources.
  atac_biosources <- unique(vapply(atac_metadata,function(x) { tolower(x$sample_info$biosource_name) },character(1L)))
  
  # 2nd Step: Collect ChiP experiments and fetch metadata
  # Normally, it would suffice to call deepblue_list_experiments() once with all
  # biosources, but it may return experiments with misspelled biosources (e.g.
  # "t-47d" when "t47d" was requested), so the function is called once for each
  # biosource and the ChIP biosources are later overwritten by their ATAC
  # equivalents.
  
  chip_experiments <- rbindlist(
    lapply(atac_biosources,function(bsource) {
      experiments <- suppressMessages(deepblue_list_experiments(genome=genome, technique="ChIP-Seq", biosource=bsource, type=chip_type, epigenetic_mark=chip_marks))
      if(is.list(experiments)) {
        return(expand.grid(id=experiments$id, bsource=bsource))
      }
    })
  )

  if(nrow(chip_experiments) == 0) { 
  
    stop(paste(genome,"No ChIP-seq data available for given arguments",sep=": "))
    
  } else {
    
    message(paste("fetching",nrow(chip_experiments),"ChIP-seq experiments ..."))
    
    chip_metadata <- suppressWarnings(apply(chip_experiments,1,function(c) {
      metadata <- suppressMessages(deepblue_info(c[1]))
      metadata <- metadata[!names(metadata) %in% c("extra_metadata","upload_info","columns")]
      metadata$sample_info$biosource_name <- as.vector(c[2])
      return(metadata)
    }))
    
    chip_biosources <- unique(tolower(vapply(chip_metadata,function(x){ x$sample_info$biosource_name },character(1L))))

    # 3rd Step: Fetch ATAC metadata if biosource has available ChiP data
    # Likewise, we only need ATAC/DNAse data that matches biosources retained
    # by our ChIP-seq experiments.
    
    atac_metadata <- atac_metadata[lapply(atac_metadata,function(x){ tolower(x$sample_info$biosource_name) }) %in% chip_biosources]
    message(paste("kept",length(atac_metadata),"ATAC/DNAse-seq experiments"))
    
    if(length(atac_metadata) == 0) {
      stop(paste(genome,"No ATAC/ChIP-seq experiments could be linked"))
    }
    
    all_metadata <- c(chip_metadata,atac_metadata)
    line_count <- 0
    
    if(append) {
      for(m in all_metadata) {
        filtered_chrs <- chrs[!chrs %in% old_data[experiment_id==m$`_id`]$chromosome]
        if(length(filtered_chrs) > 0) {
          line_count <- line_count + new_row(m,output_file,filtered_chrs)
        }
      }
      message(paste(line_count,"lines added to",normalizePath(output_file)))
    } else {
      for(m in all_metadata) {
        line_count <- line_count + new_row(m,output_file,chrs)
      }
      message(paste(line_count,"lines written to",normalizePath(output_file)))
    }
    
  }
  
}

create_linking_table(args$genome,args$chromosomes,args$biosources,args$type,args$atactype,args$marks,args$directory,args$output,args$append,args$skip_verification,args$chromsizes)
