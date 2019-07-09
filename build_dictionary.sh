#!/usr/bin/env bash 

#####USES STANDARD IN#######
#build_dictionary.sh File.fasta
#Working as of 2017-01-19

newfile="$(basename $1 .fasta)"

mkdir $newfile
cd $newfile
mv ../$1 .

#Build dictionary
bwa index -p $newfile \
-a is $1

Dictionary="$newfile.dict"

if [ -f  "$Dictionary" ]; 
then 
	echo "$newfile.dict already created"
else
	echo "file does not exist"
	java -jar /local/cluster/picard-tools-2.0.1/dist/picard.jar CreateSequenceDictionary \
	R=$1 \
	O=$newfile.dict
fi

#Samtools: index ref seq
samtools faidx $1