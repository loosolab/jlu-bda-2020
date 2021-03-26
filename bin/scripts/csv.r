#!/usr/bin/env Rscript
required_packages <- c("BiocManager","data.table","argparse")
not_installed <- setdiff(required_packages,rownames(installed.packages()))
if(length(not_installed) > 0) install.packages(not_installed,repos="http://cran.us.r-project.org")
if(!"DeepBlueR" %in% installed.packages()) BiocManager::install("DeepBlueR")
library(argparse)

parser <- ArgumentParser()

parser$add_argument("-g", "--genome", type="character", default="hg19",
                    help="Genome to search Deepblue database with [default: \"%(default)s\"]")
parser$add_argument("-c", "--chromosomes", nargs="+", type="character", default="all",
                    help="(List of) chromosomes to include (requires chr prefix) [default: all]")
parser$add_argument("-b", "--biosources", nargs="+", type="character", default="all",
                    help="(List of) biosources to include [default: all] (Refer to: https://deepblue.mpi-inf.mpg.de/)")
parser$add_argument("-t", "--type", nargs="+", type="character", default="peaks", choices=c("peaks","signal"),
                    help="Experiment file types allowed for CHiP-Seq data [default: \"%(default)s\"]")
parser$add_argument("-a", "--atactype", nargs="+", type="character", default="signal", choices=c("peaks","signal"),
                    help="Experiment file types allowed for ATAC/DNAse-Seq data [default: \"%(default)s\"]")
parser$add_argument("-m", "--marks", nargs="+", type="character", default="all",
                    help="(List of) epigenetic marks (i.e. transcription factors) to include [default: all] (Refer to: https://deepblue.mpi-inf.mpg.de/)")
parser$add_argument("-d", "--directory", type="character", default=".",
                    help="Output directory for CSV file. If an output filename (-o) containing \"/\" characters is provided, the filename will be used instead. [default: \"%(default)s\"]")
parser$add_argument("-o", "--output", type="character", default="linking_table",
                    help="Output file name without extension [default: \"%(default)s\"]")

args <- parser$parse_args()

library(data.table)

suppressPackageStartupMessages(
  library(DeepBlueR)
)

create_linking_table <- function(genome,chroms,filter_biosources,chip_type,atac_type,chip_marks,outdir,outfile) {
  
  # ! doc incomplete
  # new_row(metadata)
  # Expects a vector of filenames (strings) and a vector of chromosomes (strings)
  # Returns a vector of filenames with the given chromosomes inserted at the end but (if available) before the last extension
  # Example:
  # input:    filenames=c("MS034301.CM.signal_reverse.bedgraph","ENCFF001VOQ"), chroms=c("chr4","chr5")
  # output:   "MS034301.CM.signal_reverse.chr4.bedgraph,"MS034301.CM.signal_reverse.chr5.bedgraph",
  #           "ENCFF001VOQ.chr4","ENCFF001VOQ.chr5"

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
        total_size=metadata$upload_info$total_size,
        as.data.table(sample_info),
        as.data.table(extra)
      )
      return(meta)
      
    })
    
    return(return_list)
  }
  
  # ! doc missing
  
  verify_filters <- function(input_values,all_values,type="values") {
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
  
  all_genomes <- suppressMessages(deepblue_list_genomes()$name)
  if(!genome %in% all_genomes) {
    stop(paste("No valid genomes provided by user. Available genomes:",paste(all_genomes,collapse=", ")))
  }
  
  all_chroms <- suppressMessages(deepblue_chromosomes(genome = genome))
  if(chroms[1] != "all") {
    chrs <- verify_filters(chroms,all_chroms$id,"chromosomes")
  } else {
    chrs <- all_chroms$id
  }
  
  if(filter_biosources[1] != "all") {
    all_biosources <- suppressMessages(deepblue_list_biosources()$name)
    filter_biosources <- verify_filters(filter_biosources,all_biosources,"biosources")
  } else {
    filter_biosources <- NULL
  }
  
  tf_marks <- suppressMessages(deepblue_list_epigenetic_marks(extra_metadata = list(category="Transcription Factor Binding Sites"))$name)
  
  if(chip_marks[1] != "all") {
    chip_marks <- verify_filters(chip_marks,tf_marks,"epigenetic marks")
  } else {
    chip_marks <- tf_marks
  }
  
  # ATAC-Seq: signals (mark: DNA Accessibility)
  # DNAse-Seq: signals and peaks (marks: DNA Accessibility, DNAseI)
  
  # 1st Step: Collect biosources for ATAC
  
  atac <- suppressMessages(deepblue_list_experiments(genome=genome, technique="ATAC-Seq", biosource=filter_biosources, type=atac_type))
  dnase <- suppressMessages(deepblue_list_experiments(genome=genome, technique="DNAse-Seq", biosource=filter_biosources, type=atac_type))
  
  if(is.data.table(atac)) atac <- atac$id else atac <- NULL
  if(is.data.table(dnase)) dnase <- dnase$id else dnase <- NULL
  
  access_experiments <- c(atac,dnase)
  
  if(length(access_experiments) == 0) {
    
    warning(paste(genome,"No ATAC-seq data available for given arguments",sep=": "))
    
  } else {
    
    atac_metadata <- lapply(access_experiments,function(x){ suppressMessages(deepblue_info(x)) })
    atac_biosources <- unique(vapply(atac_metadata,function(x) { tolower(x$sample_info$biosource_name) },character(1L)))
    
    # 2nd Step: Collect ChiP experiments and add to csv list if biosource has available ATAC data
    
    chips <- suppressMessages(deepblue_list_experiments(genome=genome, technique="ChiP-Seq", biosource=atac_biosources, type=chip_type, epigenetic_mark=chip_marks))
    
    if(!is.data.table(chips)) { 
    
        warning(paste(genome,"No CHIP-seq data available for given arguments",sep=": "))
      
    } else {
      
      chips <- chips$id  
      chip_metadata <- lapply(chips,function(c) {
        metadata <- suppressMessages(deepblue_info(c))
      })
      
      chip_biosources <- unique(tolower(vapply(chip_metadata,function(x){ x$sample_info$biosource_name },character(1L))))
      
      chip_list <- lapply(chip_metadata,new_row)
      
      # before: List of lists of data.tables
      chip_list <- unlist(chip_list,recursive=FALSE)
      # after:  List of data.tables
      
      # 3rd Step: Add ATAC experiments to csv if biosource has available ChiP data
      
      atac_metadata <- atac_metadata[lapply(atac_metadata,function(x){ tolower(x$sample_info$biosource_name) }) %in% chip_biosources]
      
      atac_list <- unlist(lapply(atac_metadata,new_row),recursive=FALSE)
      
      csv_data <- rbindlist(c(chip_list,atac_list),fill=TRUE)
      
    }
  }
  
  # If "output" argument contains "/", ignore "directory" argument
  
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
  
  sizefile <- paste(genome,"chrom","sizes",sep=".")
  size_output <- paste(outdir,sizefile,sep="/")
  if(!file.exists(size_output)) {
    write.table(all_chroms,file=size_output,col.names=F,row.names=F,quote=F,sep="\t")
    message(paste(sizefile,"written to",normalizePath(outdir)))
  }
  
  if(is.null(csv_data)) {
    
    stop("No data available for CSV")
    
  } else {
    
    # check whether CSV with given filename already exists - if yes, add new rows
    if(file.exists(output_file)) {
      
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

create_linking_table(args$genome,args$chromosomes,args$biosources,args$type,args$atactype,args$marks,args$directory,args$output)
