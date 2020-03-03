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

workDIR=..
scriptDIR=.
coolFile=${workDIR}/mcool/combine_HIC_3_6_7_bw5kb.50000.cool
domains=( left center right )
binSize=50000

baseFileName=`basename ${coolFile}`
#fix name
baseFileName=${baseFileName%5kb.${binSize}.cool}50kb
echo "base file name: " $baseFileName

chrDomainDIR=${workDIR}/chrDomain
mkdir -p ${chrDomainDIR}

for domain in ${domains[@]}
do
  ranges=(`cut -f 1 ${scriptDIR}/${domain}.chrom.ranges`)
  pseudoChr=(`cut -f 2 ${scriptDIR}/${domain}.chrom.ranges`)
  roundedStart=(`cut -f 3 ${scriptDIR}/${domain}.chrom.ranges`)

  let maxIndex=${#ranges[@]}-1
  
  if [ -f ${chrDomainDIR}/${baseFileName}_${domain}Arm.bg2 ]
  then 
    rm ${chrDomainDIR}/${baseFileName}_${domain}Arm.bg2
  fi
  touch ${chrDomainDIR}/${baseFileName}_${domain}Arm.bg2

  for i in $(seq 0 $maxIndex)
  do

    echo "${pseudoChr[$i]} ${domain} starting"
    cooler dump -r ${ranges[$i]} --join ${coolFile} --one-based-starts -o  ${chrDomainDIR}/${baseFileName}_${pseudoChr[$i]}_${domain}_tmp.bg2

    echo "${pseudoChr[$i]} ${domain} dumped"

    awk -F'\t' -v chr=${pseudoChr[$i]} -v start=${roundedStart[$i]} 'BEGIN {OFS=FS} {$1=chr} {$2=$2-start} {$3=$3-start} {$4=chr} {$5=$5-start} {$6=$6-start} {print}'  ${chrDomainDIR}/${baseFileName}_${pseudoChr[$i]}_${domain}_tmp.bg2  >  ${chrDomainDIR}/${baseFileName}_${pseudoChr[$i]}_${domain}.bg2


    echo "${pseudoChr[$i]} ${domain} finished"
    cat  ${chrDomainDIR}/${baseFileName}_${pseudoChr[$i]}_${domain}.bg2 >>  ${chrDomainDIR}/${baseFileName}_${domain}Arm.bg2

  done

  rm  ${chrDomainDIR}/${baseFileName}_chr*_${domain}*.bg2

  cooler load -f bg2 ${scriptDIR}/${domain}.chrom.sizes:50000 ${chrDomainDIR}/${baseFileName}_${domain}Arm.bg2 ${chrDomainDIR}/${baseFileName}_${domain}Arm.cool

done


coolFiles=(`ls ${chrDomainDIR}/${baseFileName}_*Arm.cool`)
chrNames=( chrI chrII chrIII chrIV chrV chrX )

  ################
  # do decay plots
  #################
  # plotted per chr
  hicPlotDistVsCounts --matrices ${coolFiles[@]} --plotFile ${workDIR}/decayPlots/DistVsCounts_${baseFileName}_domains_perChr.pdf --perchr --skipDiagonal --maxdepth 20000000 --plotsize 6 5 --labels "center" "left" "right"

  # plotted all chr together
  hicPlotDistVsCounts --matrices ${coolFiles[@]} --plotFile ${workDIR}/decayPlots/DistVsCounts_${baseFileName}_domains.pdf --skipDiagonal --maxdepth 20000000 --plotsize 6 5 --labels "center" "left" "right"

  # autosomes only
  hicPlotDistVsCounts --matrices ${coolFiles[@]} --plotFile ${workDIR}/decayPlots/DistVsCounts_${baseFileName}_domains_autosomes.pdf --skipDiagonal --maxdepth 20000000 --chromosomeExclude chrX --plotsize 6 5 --labels "center" "left" "right"

  # Xchr only
  hicPlotDistVsCounts --matrices ${coolFiles[@]} --plotFile ${workDIR}/decayPlots/DistVsCounts_${baseFileName}_domains_Xchr.pdf --skipDiagonal --maxdepth 20000000 --chromosomeExclude chrI chrII chrIII chrIV chrV --plotsize 6 5 --labels "center" "left" "right"
