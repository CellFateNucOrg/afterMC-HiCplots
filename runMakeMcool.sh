#! /bin/bash
## Allocate resources
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=1
#SBATCH --partition=all

## job name
#SBATCH --job-name="Mcool"

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

