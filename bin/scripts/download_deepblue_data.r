#!/usr/bin/env Rscript
required_packages <- c("BiocManager","data.table","argparse")
not_installed <- setdiff(required_packages,rownames(installed.packages()))
if(length(not_installed) > 0) install.packages(not_installed,repos="http://cran.us.r-project.org")
if(!"DeepBlueR" %in% installed.packages()) BiocManager::install("DeepBlueR")
library(argparse)
library(data.table)

parser <- ArgumentParser()

parser$add_argument("-i", "--input", type="character", default="linking_table.csv",
                    help="Input CSV file name with extension [default: \"%(default)s\"]")
parser$add_argument("-o", "--output", type="character", default=".",
                    help="Output directory name [default: \"%(default)s\"]")
parser$add_argument("-g", "--genome", type="character", default=NULL,
                    help="Selected genome to include [default: \"%(default)s\"]")
parser$add_argument("-c", "--chromosomes", nargs="+", type="character", default=NULL,
                    help="(List of) selected chromosomes to include (requires chr prefix) [default: all]")
parser$add_argument("-b", "--biosources", nargs="+", type="character", default=NULL,
                    help="(List of) selected biosources to include [default: all]")
parser$add_argument("-m", "--marks", nargs="+", type="character", default=NULL,
                    help="(List of) selected transcription factors to include [default: all]")
parser$add_argument("--chunk_size", dest="chunks", type="integer", default=10000000L,
                    help="Chunk size for ATAC/DNAse data [default: %(default)s]")
parser$add_argument("-l", "--check_local_files", dest="local_files", type="character", default=NULL,
                    help="External directory where local files are stored. Files in queue will be copied to the output directory and not downloaded. [default: NULL]")

args <- parser$parse_args()

suppressPackageStartupMessages(
  library(DeepBlueR)
)

# export_from_csv(csv filename, output directory, chunk size)
#
# Handles the entire download process, starting with the creation of a file
# queue which it iterates through. For ATAC files, each experiment is divided
# into chunks to prevent memory errors. If files that have not been downloaded
# remain, the function will call download_regions() for each file (and chunk) to
# request the corresponding regions from the Deepblue server.

export_from_csv <- function(csv_file,out_dir,genom,bios,chroms,marks,chunk_size,local_dir) {
  
  failed_warning <- function(filename,request_id,msg) {
    warning(paste(
      filename,": download with request id",request_id,"failed:",msg
    ))
  }
  
  # download_regions(experiment id, filename, chromosome, format, chunk start, output directory)
  #
  # Attempts to download a single file from the Deepblue database. The expected
  # number of regions is requested beforehand to check for mismatches later. If
  # 0 is returned, the function assumes that no data is available for this range
  # (or chromosome) and skips to the next file. The download is considered
  # successful only if the region count matches the downloaded data and the file
  # has been exported to its expected location. In this case TRUE is returned
  # and export_to_csv() will export the metadata of the current experiment.
  
  download_regions <- function(id,filename,chr,format,chunk,out_dir) {
    
    regions_length <- 0
    message(filename)
    
    query_id <- suppressMessages(deepblue_select_experiments(experiment_name = id, chromosome = chr, start = chunk, end = chunk + chunk_size))
    message(paste("query id:",query_id))
    
    # for ATAC/DNAse-seq files (chunk > 0), filter out regions with the value 0
    
    if(chunk > 0) {
      query_id <- suppressMessages(deepblue_filter_regions(query_id = query_id, field = "VALUE", operation = "!=", value = "0", type = "number"))
      message(paste("filtered query id:",query_id))
    }
    
    count_request <- suppressMessages(deepblue_count_regions(query_id))
    expected_regions <- as.integer(suppressMessages(deepblue_download_request_data(count_request,do_not_cache=TRUE)))
    
    if(expected_regions == 0) {
      
      # this is not considered an error
      warning(paste("No regions available for",filename,"(skipping) query_id:",query_id))
      message("skipping")
      
      # output_file <- paste(out_dir,"/",filename,".txt",sep="")
      # return(file.create(output_file)) # returns TRUE if file could be created
      return(TRUE)
      
    }
    
    request_id <- suppressMessages(deepblue_get_regions(query_id = query_id, output_format = format))
    message(paste("request id:",request_id))
    
    req_status <- suppressMessages(deepblue_info(request_id)$state)
    
    if(req_status != "removed") {
      
      regions <- try(suppressMessages(deepblue_download_request_data(request_id,do_not_cache=TRUE)))
      
      if(class(regions) != "GRanges") {
        
        if(class(regions) == "try-error") {
          
          err_msg <- regions[1]
          
        } else {
          
          err_msg <- paste("received",class(regions),"class object when calling deepblue_download_request_data")
          
        }
        
        failed_warning(filename,request_id,err_msg)
        return(FALSE)
        
      } else {
        
        regions_length <- length(regions)
        
        if(regions_length > 0) {
          
          if(regions_length != expected_regions) {
            
            err_msg <- paste("regions mismatch:",expected_regions,"regions expected,",regions_length,"downloaded")
            failed_warning(filename,request_id,err_msg)
            return(FALSE)
            
          } else {
            
            message(paste("regions:",regions_length))
            output_file <- paste(out_dir,"/",filename,".txt",sep="")
            
            deepblue_export_tab(regions,target.directory = out_dir,file.name = filename)
            
            if(!file.exists(output_file)) {
              
              err_msg <- paste("file was downloaded but not exported. regions had length",regions_length)
              failed_warning(filename,request_id,err_msg)
              return(FALSE)
              
            } else {
              
              # the download was successful
              message("done")
              return(TRUE)
              
            }
            
          }
          
        } else {
          
          failed_warning(filename,request_id,"returned empty file, expected",expected_regions,"regions")
          return(FALSE)
          
        }
        
      }
      
    } else {

      failed_warning(filename,request_id,paste("request status was:",req_status))
      return(FALSE)
      
    }
    
  }
  
  if(!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  
  # This section creates a queue of files that need to be downloaded.
  
  csv_data <- fread(
    file=csv_file,
    header=TRUE,
    sep=";"
  )
  
  # If genome, biosource, chromosome or tf arguments are provided, filter the
  # linking table by the respective values.
  
  if(is.null(genom)) {
    genom <- unique(csv_data$genome)[1]
  }
  if(is.null(bios)) {
    bios <- unique(csv_data$biosource)
  }
  if(is.null(marks)) {
    marks <- unique(csv_data$epigenetic_mark)
  }
  if(is.null(chroms)) {
    chroms <- unique(csv_data$chromosome)
  }
  
  data <- csv_data[genome %in% genom & biosource %in% bios & chromosome %in% chroms & epigenetic_mark %in% c(marks,"dnasei","dna accessibility")]

  # The script checks whether there are files that have already been downloaded
  # in the output folder. Is this the case then they are removed from the queue.
  # An existing [filename].meta.txt file is considered a flag for a successfully
  # downloaded file.
  
  message("creating queue ...")
  
  meta_files <- dir(path=out_dir, pattern="meta.txt")
  downloaded_files <- gsub(".meta.txt","",meta_files)
  queued_files <- data[!data$filename %in% downloaded_files]
  files_copied <- FALSE
  
  # --check_local_files <directory>
  # Meant to reduce download queue if files already exist inside an external
  # directory. Checks for *.meta.txt files corresponding to the filenames in
  # queue and copies them to the (-o)uput directory together with the actual
  # files if they exist.
  
  if(nrow(queued_files) > 0 && !is.null(local_dir)) {
    if(dir.exists(local_dir)) {
      local_dir <- sub("/$","",local_dir)
      ext_files <- dir(path=local_dir,pattern=".meta.txt")
      ext_files.compare <- gsub(".meta.txt","",ext_files)
      ext_files.inqueue <- ext_files.compare[ext_files.compare %in% queued_files$filename]
      message(paste("found",length(ext_files.inqueue),"experiment(s) in",local_dir))
      if(length(ext_files.inqueue) > 0) {
        for(experiment in ext_files.inqueue) {
          allfiles <- dir(path=local_dir,pattern=paste("^",experiment,sep=""))
          for(single_file in allfiles) {
            output_file <- file.path(local_dir,single_file)
            message(paste("copying",output_file))
            file.copy(output_file,out_dir)
          }
          files_copied <- TRUE
        }
        queued_files <- queued_files[!filename %in% ext_files.compare]
      }
    } else {
      warning(paste(local_dir,": directory not found"))
    }
  }
  
  n_files <- nrow(queued_files)
  
  if(n_files > 0) {
    
    # get chrom sizes for all genomes in the CSV
    genomes <- unique(queued_files$genome)
    chrom_sizes <- vector("list")
    for(genome in genomes) {
      chroms <- suppressMessages(deepblue_chromosomes(genome))
      chroms$id <- tolower(chroms$id)
      chrom_sizes[[genome]] <- chroms
    }
    
    failed_rows <- queued_files
    
    for(i in 1:n_files) {
      
      progress_msg <- paste("downloading experiment",i,"of",n_files)
      
      row <- queued_files[i,]
      filename <- row$filename
      chr <- tolower(sub(".*\\.(chr.*?)(\\..*|$)","\\1",filename))
      id <- row$experiment_id
      
      if(row$technique == "chip-seq") {
        
        message(progress_msg)
        no_errors <- download_regions(id,filename,chr,row$format,0,out_dir)
        
      } else {
        
        genom <- row$genome
        this_chrom_size <- chrom_sizes[[genom]][id == chr]$name
        chunks <- seq(1,this_chrom_size,by=chunk_size)
        n_chunks <- length(chunks)
        
        for(j in 1:n_chunks) {
          
          message(paste(progress_msg,"- chunk",j,"of",n_chunks))
          chunked_filename <- paste(filename,"chunk",chunks[j],sep="_")
          no_errors <- download_regions(id,chunked_filename,chr,row$format,chunks[j],out_dir)
          if(!no_errors) break
          
        }
      }
      
      if(no_errors) {
        failed_rows[i,] <- NA
        suppressMessages(
          deepblue_export_meta_data(id,target.directory=out_dir,file.name=filename)
        )
      }
      
    }
    
    failed_files <- failed_rows[!is.na(failed_rows$filename)]$filename
    n_failed <- length(failed_files)
    
    if(n_failed > 0) {
      
      # Terminate execution
      stop(paste("Download of",n_failed,"file(s) could not be completed due to previous errors:",paste(failed_files,collapse=", ")))
      
    }
    
  } else {
    
    # All files in queue have already been downloaded (or copied)
    message("No new files to download.")
    # exit with returncode 2 for generate_data.py error handling
    # if not otherwise specified, this will skip the validation and sorting of 
    # the files. if at least one file was copied from a (-l)ocal directory, exit
    # normally since validation may be required.
    if(!files_copied) {
      q(save="no",status=2)
    }
    
  }
  
}

export_from_csv(args$input,args$output,args$genome,args$biosources,args$chromosomes,args$marks,as.integer(args$chunks),args$local_files)