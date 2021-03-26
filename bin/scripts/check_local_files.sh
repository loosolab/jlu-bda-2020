#!/bin/env bash

# ==============================================================================================
#
# FILE: check_local_files.sh
#
# USAGE: check_local_files.sh [path to local files] [download path] [csv filename]
#
# DESCRIPTION: Goes through all filenames inside the csv and checks whether the files can be
#	found in the external folder (first argument), then copies that file into the designated
#	download folder (second argument) so export_from_csv.r does not need to download these
#	files again.
#	Note: meta.txt files are required for all files in order to be copied. export_from_csv.r
#	downloads a meta.txt file after each successful download. Therefore, this script is
#	meant to be used with files that have been downloaded with export_from_csv.r.
#
# AUTHOR: Patrick Groos https://github.com/fctsfrmspc
#
# ==============================================================================================

if [ ! $# -eq 3 ]; then
	echo "check_local_files.sh: expected 3 arguments, got $#"
	exit 1
fi

local_path=$1
dest_path=$2
csv_file=$3
copied=0
i=0

for file in $(cut -d ";" -f 7 $csv_file | sed -n '1d;p')
do
	#echo $file
	if [ ! -e "$dest_path/$file.meta.txt" ]; then
		if [ -e "$local_path/$file.meta.txt" ]; then
			if [ -e "$local_path/$file.txt" ]; then
				cp "$local_path/$file.txt" "$dest_path/$file.txt"
			fi
			cp "$local_path/$file.meta.txt" "$dest_path/$file.meta.txt"
			copied=$((copied+1))
		fi
	fi
	i=$((i+1))
done

echo "$copied of $i files were copied to `pwd`/$local_path"
