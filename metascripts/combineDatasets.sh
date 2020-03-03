#! /bin/bash
# combining 2 or more different datasets using already basecalled data.
# example: ./combineDatasets.sh combine_HIC_3_6_7 20181013_HIC3PRQS/bcFastq/pass/noBC.fastq.gz 20190501_HIC6_7/bcFastq/pass/noBC.fastq.gz

if [ $# -lt 3 ]
then
	echo "Usage: $0 newDatasetName relative/path/to/file1 relative/path/to/file2 [ relative/path/to/file3 ] "
	exit 1
fi


NEWNAME=$1
echo "creating $NEWNAME dataset"
shift
ALLFILES=$@
#echo $ALLFILES

#workDir=`echo $PWD`
mkdir -p ${NEWNAME}/raw_files

cat ${ALLFILES[@]}  > ${NEWNAME}/raw_files/raw_${NEWNAME}_noBC_pass.fastq.gz


