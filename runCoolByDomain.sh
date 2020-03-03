#! /bin/bash
## Allocate resources
#SBATCH --time=0-01:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=2
##SBATCH --array=1
#SBATCH --partition=all

## job name
#SBATCH --job-name="domains"

source ${CONDA_ACTIVATE} MC-HiC-env
echo "conda activate is: " ${CONDA_ACTIVATE}

workDIR=./
coolFile=${workDIR}/mcool/combine_TEVneg_HIC_8_12.8_bw5kb.50000.cool
domains=( left center right )
binSize=50000

baseFileName=`basename ${coolFile}`
baseFileName=${baseFileName%.${binSize}.cool}

chrDomainDIR=${workDIR}/chrDomain
mkdir -p chrDomain

for domain in ${domains[@]}
do
  ranges=(`cut -f 1 ${workDIR}/${domain}.chrom.ranges`)
  pseudoChr=(`cut -f 2 ${workDIR}/${domain}.chrom.ranges`)
  roundedStart=(`cut -f 3 ${workDIR}/${domain}.chrom.ranges`)

  let maxIndex=${#ranges[@]}-1
  
  if [ -f ${chrDomainDIR}/${baseFileName}_${domain}Arm.bg2 ]
  then 
    rm ${chrDomainDIR}/${baseFileName}_${domain}Arm.bg2
  fi
  touch ${chrDomainDIR}/${baseFileName}_${domain}Arm.bg2

  for i in $(seq 0 $maxIndex)
  do

    echo "${pseudoChr[$i]} starting"
    cooler dump -r ${ranges[$i]} --join ${coolFile} --one-based-starts -o  ${chrDomainDIR}/${baseFileName}_${pseudoChr[$i]}_tmp.bg2

    echo "${pseudoChr[$i]} dumped"

    awk -F'\t' 'BEGIN {OFS=FS} {$1='${pseudoChr[$i]}'} {$4='${pseudoChr[$i]}'} {print}'  ${chrDomainDIR}/${baseFileName}_${pseudoChr[$i]}_tmp.bg2  >  ${chrDomainDIR}/${baseFileName}_${pseudoChr[$i]}.bg2

    awk -F'\t' -v chr=${pseudoChr[$i]} -v start=${roundedStart[$i]} 'BEGIN {OFS=FS} {$1=chr} {$2=$2-start} {$3=$3-start} {$4=chr} {$5=$5-start} {$6=$6-start} {print}'  ${chrDomainDIR}/${baseFileName}_${pseudoChr[$i]}_tmp.bg2  >  ${chrDomainDIR}/${baseFileName}_${pseudoChr[$i]}.bg2


    echo "${pseudoChr[$i]} finished"
    cat  ${chrDomainDIR}/${baseFileName}_${pseudoChr[$i]}.bg2 >>  ${chrDomainDIR}/${baseFileName}_${domain}Arm.bg2

  done
done

#rm  ${chrDomainDIR}/${baseFileName}_chr*_*.bg2

cooler load -f bg2 ${workDIR}/${domain}.chrom.sizes:50000 ${chrDomainDIR}/${baseFileName}_${domain}Arm.bg2 ${chrDomainDIR}/${baseFileName}_${domain}Arm.cool









coolFiles=(`ls hic_mats/*5000.cool`)
COOLER_RESOLUTIONS=5000,10000,20000,50000,100000,200000,500000


mkdir -p mcool

for coolFile in ${coolFiles[@]} 
do
	fileName=`basename $coolFile`
	mv $coolFile mcool/$fileName
	#mcoolFile=mcool/${fileName%cool}mcool
	cooler zoomify 	-n 4 -c 10000000 -r $COOLER_RESOLUTIONS mcool/$fileName
done

