#! /bin/bash
## Allocate resources
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=4
#SBATCH --array=1
#SBATCH --partition=all

## job name
#SBATCH --job-name="hicExp"

source ${CONDA_ACTIVATE} MC-HiC-env

#"$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && cd ../../ && pwd )"
workDIR=${HOME}
scriptDIR=.
decayDIR=${workDIR}/decayPlots_compare

chromosomes=( 1 2 3 4 5 6 )
chrNames=( chrI chrII chrIII chrIV chrV chrX )
#resolutions=( 5000 10000 20000 50000 100000 200000 500000 )

h5Files=( ${workDIR}/combine_TEVneg_HIC_8_12.8/h5_files/combine_TEVneg_HIC_8_12.8_bw5kb.50000.h5 ${workDIR}/combine_TEVpos_HIC_10_12.5_13/h5_files/combine_TEVpos_HIC_10_12.5_13_bw5kb.50000.h5 )

mkdir -p ${decayDIR}

echo ${h5Files[@]}

resolution=50000

################
# do decay plots
#################
# plotted per chr
hicPlotDistVsCounts --matrices ${h5Files[@]} --plotFile ${decayDIR}/DistVsCounts_TEVposTEVneg_perChr.pdf --perchr --skipDiagonal --maxdepth 20000000 --plotsize 6 5 --labels "TEVcs-" "TEVcs+"

# plotted all chr together
hicPlotDistVsCounts --matrices ${h5Files[@]} --plotFile ${decayDIR}/DistVsCounts_TEVposTEVneg.pdf --skipDiagonal --maxdepth 20000000 --plotsize 6 5 --labels "TEVcs-" "TEVcs+"

# autosomes only
hicPlotDistVsCounts --matrices ${h5Files[@]} --plotFile ${decayDIR}/DistVsCounts_TEVposTEVneg_autosomes.pdf --skipDiagonal --maxdepth 20000000 --chromosomeExclude 6 --plotsize 6 5 --labels "TEVcs-" "TEVcs+"

# Xchr only
hicPlotDistVsCounts --matrices ${h5Files[@]} --plotFile ${decayDIR}/DistVsCounts_TEVposTEVneg_Xchr.pdf --skipDiagonal --maxdepth 20000000 --chromosomeExclude 1 2 3 4 5 --plotsize 6 5 --labels "TEVcs-" "TEVcs+"


domains=( left center right )
for domain in ${domains[@]}
do
  
  ################
  # do decay plots by domain
  #################
  
  coolFiles=( /home/mdas/combine_TEVneg_HIC_8_12.8/chrDomain/combine_TEVneg_HIC_8_12.8_bw50kb_${domain}Arm.cool /home/mdas/combine_TEVpos_HIC_10_12.5_13/chrDomain/combine_TEVpos_HIC_10_12.5_13_bw50kb_${domain}Arm.cool )   

  # plotted per chr
  hicPlotDistVsCounts --matrices ${coolFiles[@]} --plotFile ${decayDIR}/DistVsCounts_TEVposTEVneg_${domain}Arm_perChr.pdf --perchr --skipDiagonal --maxdepth 20000000 --plotsize 6 5 --labels "TEVcs-" "TEVcs+"

  # plotted all chr together
  hicPlotDistVsCounts --matrices ${coolFiles[@]} --plotFile ${decayDIR}/DistVsCounts_TEVposTEVneg_${domain}Arm.pdf --skipDiagonal --maxdepth 20000000 --plotsize 6 5 --labels "TEVcs-" "TEVcs+"

  # autosomes only
  hicPlotDistVsCounts --matrices ${coolFiles[@]} --plotFile ${decayDIR}/DistVsCounts_TEVposTEVneg_${domain}Arm_autosomes.pdf --skipDiagonal --maxdepth 20000000 --chromosomeExclude chrX --plotsize 6 5 --labels "TEVcs-" "TEVcs+"

  # Xchr only
  hicPlotDistVsCounts --matrices ${coolFiles[@]} --plotFile ${decayDIR}/DistVsCounts_TEVposTEVneg_${domain}Arm_Xchr.pdf --skipDiagonal --maxdepth 20000000 --chromosomeExclude chrI chrII chrIII chrIV chrV --plotsize 6 5 --labels "TEVcs-" "TEVcs+"

done

#  distance=2000000
#  prettyDist=`expr $distance / 1000000`
#  prettyDist=`echo thresh${prettyDist}Mb`
#  hicPlotSVL -m $h5File -pfn decayPlots/shortVlong_${prettyDist}_${baseName}.pdf -o decayPlots/shortVlong_pval_${prettyDist}_${baseName}.tsv -od decayPlots/shortVlong_data_${prettyDist}_${baseName}.tsv -d $distance -t $SLURM_CPUS_PER_TASK
#  
#  
#  
#  ################
#  # calculate eigenVectors
#  #################
#  echo "finding eigen vectors..."
#  hicPCA --matrix $h5File --outputFileName eigenVectors/${baseName}_pca1.bw eigenVectors/${baseName}_pca2.bw
#  
#  # convert eigen vectors to ce11 chr names
#  echo "converting to eigenVectors to ce11..."
#  ./convertBigWigChroms.py conversion.txt eigenVectors/${baseName}_pca1.bw eigenVectors/${baseName}_pca1_ce11.bw
#  ./convertBigWigChroms.py conversion.txt eigenVectors/${baseName}_pca2.bw eigenVectors/${baseName}_pca2_ce11.bw
#  
#  
#  ################
#  # detect loops
#  #################
#  for chr in ${chromosomes[@]}
#  do
#  	let j=${chr}-1  
#	echo "detecting loops..."
#  	hicDetectLoops --matrix $h5File --outFileName loops/${baseName}_loops_${chrNames[$j]}.bedgraph --chromosomes $chr
#  	
#	echo "plotting loops..."
#	hicPlotMatrix -m $h5File -o loops/${baseName}_loops_${chrNames[$j]}.png --log1p --region $chr --loops loops/${baseName}_loops_${chrNames[$j]}.bedgraph
#	
#  done
#
#done

#hicTransform --matrix ${h5File} -o eigenVectors/${baseName}_obsexp.h5

#hicCompartmentsPolarization -m eigenVectors/${baseName}_obsexp.h5 --pca eigenVectors/${baseName}_pca1.bw -o eigenVectors/${baseName}_compartments

