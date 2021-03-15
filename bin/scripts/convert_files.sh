#!/bin/env bash
# TODO: wenn file schon richtig, konvertierung skippen
#===============================================================================
#
#  FILE:  convert.sh
#
#  USAGE:  convert.sh [filetype to convert to] [path to files]
#					  [genome.chrom.size folder path] [name of csv]
#
#  DESCRIPTION:  Validate and convert  files in src_path according to the
#				parameters. At the same time ensure proper file-naming.
#
#  NOTES:  ---
#  AUTHOR:  Jonathan Schäfer
#===============================================================================


#==== Function =================================================================
#  Name: validate_file
#  Description: Validates that the file extension of a file fits the content of
#  a file by comparing the filename to the format provided in the .csv
#  $1 = filename of the file to check
#  $2 = format of the filecontent
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
		;;
	*)
		echo "unrecognized file format"
		new_filename=$1 #filename stays the same
		return 2
		;;
	esac
	if [ "$filetype" != "$file_extension" ]; then #filename gets new ending
		new_filename=$(basename "$1")
		new_filename="${new_filename%.*}.$filetype"
		return 1
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
	if [ "$3" == "bigwig" ] || [ "$3" == "bw" ]; then
		if  [ "$file_extension" == "bed" ]; then
			cut --fields 1-3,7 "$1" > "$out_path/$file_name.bedgraph"
			file_extension="bedgraph"
		fi
		if [ "$file_extension" == "bedgraph" ]; then
			newfile="$out_path/$file_name.$file_extension"
			tail -n +2 "$newfile" > "$6/tempfile"
			mv "$6/tempfile" "$newfile"
			bedGraphToBigWig "$newfile" \
				"$5/$4.chrom.sizes" "$out_path/$file_name.bw"
			new_filename="$file_name.bw"
		else
				echo "unexpected file: $2 $file_extension" # TODO: proper error handling
		fi
	fi
}


#==== Function =================================================================
#  Name: merge_chunks
#  Description: merges Atac-seq file chunks into one file
#  $1 = input folder
#===============================================================================
merge_chunks() {
	folder=$1
	for file in $(ls -v "$folder"); do
		file=$folder/$file
		temp=$(basename "$folder/$file")
		filename=${temp/_chunk*/}
		if [[ $file == *"chunk"* ]]; then
			if [[ $outfile != $folder/$filename ]]; then
				outfile=$folder/$filename
				# TODO: echo "outfile: $outfile" Logging
				head -n1 "$file" > "$outfile"
			fi
			awk 'FNR>1' "$file" >> "$outfile"
		fi
	done
}

filetype=$1
source_path=$2
out_path=$3
chrom_path=$4
csv_name=$5

new_link=$out_path/$csv_name
touch "$new_link"
export new_filename=""

# Strip .txt ending of downloaded files
#rename s/'.txt'// $source_path/*.txt # TODO: error when doublequoting source_path/*

#Merge atac-seq chunks
#merge_chunks "$source_path"

headers=$(head -n 1 "$source_path/$csv_name")
echo "${headers:0:87}file_path;${headers:87}" > "$new_link"
#===============================================================================
# Goes through all lines of the .csv and validates the file before attempting
# to convert it to the proper filetype
#===============================================================================
while IFS=";" read -r experiment_id	genome	biosource	technique	\
	epigenetic_mark chromosome	filename	data_type	format remaining
do
	if [ ! -e "$source_path/$filename" ]; then
		continue
	fi

	source_file="$source_path/$filename"

	validate_filetype "$filename" "$format"
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
done < <(tail --lines +2 "$source_path/$csv_name")
