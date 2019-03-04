#!/bin/bash
#SBATCH -J s1p1C3k
#SBATCH -o s1p1C3k_%j.txt
#SBATCH -e errs1p1C3k_%j.err
#SBATCH -p skx-dev
##SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
##SBATCH -c 4
#SBATCH --export=ALL
#SBATCH --time=2:00:00
#SBATCH -A TG-EAR170019

export OMP_NUM_THREADS=2
export MV2_ENABLE_AFFINITY=0

cd /work/04149/tg834482/stampede2/release0/NormalModes/
source SetEnv
cd demos/
#ibrun valgrind --error-limit=no ./plmvcg_stampede2.out 
ibrun ./../bin/plmvcg_stampede2.out 
