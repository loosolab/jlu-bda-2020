#!/bin/env bash

#===============================================================================
#
#  FILE:  sort.sh
#
#  USAGE:  sort.sh [data directory path] [destination directory path] [csv_path]
# 				   [csv name]
#
#  DESCRIPTION:  Sort files in data into a folderstructure in the destination
#				according to the parameters in the csv table, At the same time
#				adds the new filepath and cleans up the .csv. Files that do not
# 				exist anymore are deleted out of the .csv, this should only
#				affect the merging of forward/reverse reads, since those files
#				are only needed in the merged form.
#
#  $1 = path to the directory where the data is
#  $2 = path to the directory where the files are sorted into
#  $3 = path to the csv with the entries to sort
#  $4 = name of the linking table to write the entries to

#  NOTES:  ---
#  AUTHOR:  Jonathan Schaefer (JonnyCodewalker)
#===============================================================================
echo "sorting files"

source_path=$1
dest_path=$2
csv_path=$3
csv_name=$4
new_link=$dest_path/$csv_name

headers=$(head -n1 "$csv_path")
if [ ! -f "$new_link" ]; then
	echo "$headers" > "$new_link"
fi

#===============================================================================
# Goes through all lines of the .csv and sorts the files according to the
# sequencing technique, skipping non existing files to clean up the linkingtable
#===============================================================================

while IFS=";" read -r experiment_id	genome	biosource	technique	\
	epigenetic_mark	chromosome filename	data_type file_path	remaining
do
	if [ ! -e "$source_path/$filename" ]; then
		continue
	fi
	# lowercase comparison to ensure a match
	check=$(echo "$technique" | awk '{print tolower($0)}')
	if [ "$check" == "atac-seq" ]; then
		new_path="$dest_path/$genome/$biosource/$technique"
	elif [ "$check" == "chip-seq" ]; then
		new_path="$dest_path/$genome/$biosource/$technique/$epigenetic_mark"
	elif [ "$check" == "dnase-seq" ]; then
		new_path="$dest_path/$genome/$biosource/atac-seq"
	else
		new_path="$dest_path/$genome/$biosource/$technique/$epigenetic_mark"
	fi

	mkdir -p "$new_path"
	sourcefile="$source_path/$filename"
	file_path="$new_path/$filename"
	mv "$sourcefile" "$file_path"
	line="$experiment_id;$genome;$biosource;$technique;$epigenetic_mark;\
$chromosome;$filename;$data_type;$file_path;$remaining"

	# ensure that there are no double entries in the linking table
	if ! grep -Fxq "$line" "$new_link"; then
		echo "$line" >> "$new_link"
	fi
done < <(tail -n +2 "$csv_path")
