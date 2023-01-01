#!/bin/bash
#SBATCH -J s3p1RC3k
#SBATCH -o s3p1RC3k_%j.txt
#SBATCH -e errs3p1RC3k_%j.err
#SBATCH -p skx-dev
##SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
##SBATCH -c 4
#SBATCH --export=ALL
#SBATCH --time=0:10:00
##SBATCH -A TG-CCR110036
##SBATCH -A TG-EAR160023
#SBATCH -A TG-EAR170019
#SBATCH --mail-user=shijia1019@gmail.com
#SBATCH --mail-type=all

export OMP_NUM_THREADS=2
export MV2_ENABLE_AFFINITY=0

cd /work/04149/tg834482/stampede2/nm0/intel18/PNM_newP2/
source SetEnv
cd bin
ibrun ./plmvcg_stampede2.out

