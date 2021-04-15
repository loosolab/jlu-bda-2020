#!/bin/env bash
#===============================================================================
#
#  FILE:  convert.sh
#
#  USAGE:  convert.sh [filetype to convert to] [path to files]
#					  [genome.chrom.size folder path] [name of csv]
#					  [logfile path]
#
#  DESCRIPTION:  Validate and convert  files in src_path according to the
#				parameters. At the same time ensure proper file-naming.
#  $1 = filetype to convert to
#  $2 = path to the input data (./data/download)
#  $3 = path to the output data(./data/)
#  $4 = name of the linking table
#  $5 = path to the logfile
#
#  NOTES:  ---
#  AUTHOR:  Jonathan Schaefer (JonnyCodewalker)
#===============================================================================


#==== Function =================================================================
#  Name: validate_file
#  Description: Validates that the file extension of a file fits the content of
#  a file by comparing the filename to the format provided in the .csv
#  If the file is a .bedfile the header of the file is looked at to ensure the
#  right column is cut when producing a .bedgraph.
#  sets the new_filename to properly fit the filetype.
#
#  $1 = filename of the file to check
#  $2 = format of the filecontent
#
#  Returns 2 if the fileformat is unusable or not recognized and 0 otherwise
#===============================================================================
validate_filetype () {
	local file_extension=${1##*.}
	local filetype
	case "$2" in
	"CHROMOSOME,START,END,VALUE")
		filetype="bedgraph"
		;;
	"CHROMOSOME,START,END,NAME,SCORE,STRAND,SIGNAL_VALUE,P_VALUE,Q_VALUE,PEAK" \
	| "CHROMOSOME,START,END,NAME,SCORE,STRAND,SIGNAL_VALUE,P_VALUE,Q_VALUE" \
	| "CHROMOSOME,START,END,NAME,SCORE"\
	| "CHROMOSOME,START,END,NAME")
		filetype="bed"
		# calculate the index of the value column to properly cut into bedgraph
		# since many files do not actually have the stated format internally
		header=$(head -n 1 "$source_path/$1")
		header=($header) #intentionally split at tab into array
		cut_value=0
		for i in "${!header[@]}"; do
		   if [[ "${header[$i]}" = "SIGNAL_VALUE" ]]; then
		   	   cut_value=$((i+1))
		   fi
		done
		;;
	*)
		echo "INFO:root:unrecognized file format, $1 can not be used" >> "$logfile"
		return 2
		;;
	esac
	if [ "$filetype" != "$file_extension" ]; then #filename gets new ending
		new_filename=$(basename "$1")
		new_filename="${new_filename}.$filetype"
	else
		new_filename=$1 # filename stays the same
	fi

	return 0
}

#==== Function =================================================================
#  Name: convert_file
#  Description: Converts a file into the requested format (atm only bigwig).
#  converts .bed to .bedgraph and then calls the appropriate ucsc tool
#  $1 = file to convert
#  §2 = filetype to convert to
#  $3 = genome of the filecontent, needed for chrom.sizes
#===============================================================================
convert_file() {
	local file_extension=${2##*.}
	local file_name=${2%.*}
	if [ -e "$out_path/$file_name.bw" ]; then
		return 0
	fi
	if [ "$3" == "bigwig" ] || [ "$3" == "bw" ]; then
		if  [ "$file_extension" == "bed" ]; then
			# convert file into a .bedgraph
			cut --fields 1-3,$cut_value "$1" > "$out_path/$file_name.bedgraph"
			file_extension="bedgraph"
		fi
		if [ "$file_extension" == "bedgraph" ]; then
			newfile="$out_path/$file_name.$file_extension"
			# fix for some bedgraphs not having proper format
			if [ "$(head -1 "$newfile" | tr '\t' '\n' | wc -l)" == 6 ]; then
				cut --fields 1-3,6 "$newfile" > "$6/tempfile"
				mv "$6/tempfile" "$newfile"
			fi
			# convert file into a .bw
			new_filename="$file_name.bw"
			tail -n +2 "$newfile" > "$6/tempfile"
			mv "$6/tempfile" "$newfile"
			bedGraphToBigWig "$newfile" \
				"$5/$4.chrom.sizes" "$out_path/$new_filename" &>/dev/null
			if [ $? == 255 ]; then # errors when overlapping regions are present
				echo "INFO:root:overlap found in $2" >> "$logfile"
				bedRemoveOverlap "$newfile" "$6/tempfile"
				mv "$6/tempfile" "$newfile"
				bedGraphToBigWig "$newfile" \
				"$5/$4.chrom.sizes" "$out_path/$new_filename"
			fi
		else
				echo "ERROR:root:unexpected file: $2 $file_extension" >> "$logfile"
		fi
	fi
}


#==== Function =================================================================
#  Name: merge_chunks
#  Description: merges file chunks into one file. This is needed because
#  large files downloaded from Deepblue are split into chunks.
#  $1 = input folder
#===============================================================================
merge_chunks() {
	folder=$1
	outfiles=()
	for file in $(ls -v "$folder"); do #bad, but I have not found a better way
		file="$folder/$file"
		temp=$(basename "$folder/$file")
		filename=${temp/_chunk*/}
		if [[ $file == *"chunk"* ]]; then
			if [[ $outfile != $folder/$filename ]]; then #create outfile
				outfile="$folder/$filename"
				echo "INFO:root:chunks merged into: $outfile" >> "$logfile"
				head -n1 "$file" > "$outfile"
			fi
			outfiles+=("$outfile")
			awk 'FNR>1' "$file" >> "$outfile"
			rm "$file" # remove the chunkfile after it has been merged
		fi
	done
	merged_files=("$(echo "${outfiles[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ')")
	for file in ${merged_files[*]}; do
		uniq "$file" > "tempfile"
		mv "tempfile" "$file"
	done
}

echo "---- File validation ----"
echo "Some files can take a long time to validate"
filetype=$1
source_path=$2
out_path=$3/temp
chrom_path=$3/chromsizes
csv_name=$4
logfile=$5

# create new linking table
new_link=$out_path/$csv_name
touch "$new_link"
export new_filename=""

# move chrom.sizes to proper folder
mv "$source_path"/*.chrom.sizes "$chrom_path" &>/dev/null
# Strip .txt ending of downloaded files
rename '.txt' '' "$source_path"/*.txt &>/dev/null
rename '.meta' '.meta.txt' "$source_path"/*.meta &>/dev/null

#Merge chunks
echo "merging chunks"
merge_chunks "$source_path"

echo "validating files"
headers=$(head -n 1 "$source_path/$csv_name")
echo "${headers:0:87}file_path;${headers:87}" > "$new_link" #add header
count=0
#===============================================================================
# Goes through all lines of the .csv and validates the file before attempting
# to convert it to the proper filetype.
# if the file exisst it is validated and skipped if bad format
# otherwise it is converted into a .bigwig and the entries are added into
# the created linking table
#===============================================================================
while IFS=";" read -r experiment_id	genome	biosource	technique	\
	epigenetic_mark chromosome	filename	data_type	format remaining
do
	if [ ! -e "$source_path/$filename" ]; then
		continue
	fi

	source_file="$source_path/$filename"
	validate_filetype "$filename" "$format"
	if [ $? == 2 ]; then # invalid file, does not contain needed information
		continue
	fi

	cp "$source_file" "$out_path/$new_filename"
	source_file="$out_path/$new_filename"
	convert_file "$source_file" "$new_filename" "$filetype" "$genome" "$chrom_path" "$out_path"
	filepath="$out_path/$new_filename"
	# add .bw file
	echo "$experiment_id;$genome;$biosource;$technique;$epigenetic_mark;\
$chromosome;$new_filename;$data_type;$filepath;$format;$remaining"\
	>> "$new_link"
	# add .bed file
	new_filename="${new_filename%.*}.bed"
	echo "$experiment_id;$genome;$biosource;$technique;$epigenetic_mark;\
$chromosome;$new_filename;$data_type;$filepath;$format;$remaining"\
	>> "$new_link"

	count=$((count+1))
	if ! (( count % 5 )) ; then
		echo "validated $count files"
	fi
done < <(tail --lines +2 "$source_path/validation.csv")
