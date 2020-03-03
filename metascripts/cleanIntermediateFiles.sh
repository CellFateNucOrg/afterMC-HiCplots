#! /bin/bash

all=$1

echo "deleting fastq_files"
rm -rf fasta_files
echo "deleting split_files"
rm -rf split_files
echo "deleting fastqFiles pass, fail and log files"
rm -rf fastqFiles/pass
rm -rf fastqFiles/fail
rm -rf fastqFiles/guppy_basecaller_log*
echo "deleting bcFastq fail files"
rm -rf bcFastq/fail
echo "deleting slurm output files"
rm -rf slurm*


if [[ "$all" == "all" ]]
then
	echo "deleting bcFastq"
	rm -rf bcFastq
fi
