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

chromosomes=( 1 2 3 4 5 6 )
chrNames=( chrI chrII chrIII chrIV chrV chrX )
resolutions=( 5000 10000 20000 50000 100000 200000 500000 )

h5Files=( combine_TEVneg_HIC_8_12.8/h5_files/combine_TEVneg_HIC_8_12.8_bw5kb.50000.h5 combine_TEVpos_HIC_10_12.5_13/h5_files/combine_TEVpos_HIC_10_12.5_13_bw5kb.50000.h5 )

echo ${h5Files[@]}

resolution=50000
#for resolution in ${resolutions[@]}
#do
#  ##resolution=5000
#  ##prettyRes=`expr ${resolution} / 1000`
#  ##prettyRes=`echo res${prettyRes}kb`
#  coolFiles=(`ls mcool/*.${resolution}.cool | grep -v "_seq" | grep -v "_z"`)
#  
#  ##coolFile="/Users/semple/IdeaProjects/hicExplorerPlots/hic_mats/HiCmat_20191026_HIC13_noBC_pass_bw10kb_mf50_ice.10000.cool"
#  ##baseName="20191026_HIC13_noBC_pass_bw10kb_mf50_ice"
#  
#  mkdir -p h5_files
#  mkdir -p decayPlots
#  mkdir -p eigenVectors
#  mkdir -p loops
#  
#  echo $SLURM_ARRAY_TASK_ID "is slurm task id"
#  let i=$SLURM_ARRAY_TASK_ID-1
#  coolFile=${coolFiles[$i]}
#  
#  echo $coolFile
#  
#  ###############
#  # convert cool to h5 (not necessary?)
#  ###############
#  coolFileName=`basename ${coolFile}`
#  baseName=${coolFileName%.${resolution}.cool}
#  baseName=${baseName#HiCmat_}
#  baseName=`sed -nr 's/(.*)_noBC_pass(.*)_mf50(.*)/\1\2\3/p' <<< "$baseName"`
#  echo $baseName " is basename"
#  
#  coolPath=`dirname ${coolFile}`
#  
#  hicInfo -m $coolFile #--outFileName
#  
#  h5FileName=${coolFileName%cool}h5
#  h5Path=${coolPath%mcool}h5_files
#  
#  h5File=${h5Path}/${h5FileName}
#  
#  hicConvertFormat -m $coolFile --inputFormat cool --outputFormat h5 -o $h5File
  
  ################
  # do decay plots
  #################
  # plotted per chr
  hicPlotDistVsCounts --matrices ${h5Files[@]} --plotFile decayPlots/DistVsCounts_TEVposTEVneg_perChr.pdf --perchr --skipDiagonal --maxdepth 20000000 --plotsize 6 5 --labels "TEVcs-" "TEVcs+"

  # plotted all chr together
  hicPlotDistVsCounts --matrices ${h5Files[@]} --plotFile decayPlots/DistVsCounts_TEVposTEVneg.pdf --skipDiagonal --maxdepth 20000000 --plotsize 6 5 --labels "TEVcs-" "TEVcs+"

  # autosomes only
  hicPlotDistVsCounts --matrices ${h5Files[@]} --plotFile decayPlots/DistVsCounts_TEVposTEVneg_autosomes.pdf --skipDiagonal --maxdepth 20000000 --chromosomeExclude 6 --plotsize 6 5 --labels "TEVcs-" "TEVcs+"

 # Xchr only
    hicPlotDistVsCounts --matrices ${h5Files[@]} --plotFile decayPlots/DistVsCounts_TEVposTEVneg_Xchr.pdf --skipDiagonal --maxdepth 20000000 --chromosomeExclude 1 2 3 4 5 --plotsize 6 5 --labels "TEVcs-" "TEVcs+"

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

