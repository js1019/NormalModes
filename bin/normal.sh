#!/bin/bash
#SBATCH -J s3E80k
#SBATCH -o s3E80k_%j.txt
#SBATCH -e errs3E80k_%j.err
#SBATCH -p skx-normal
##SBATCH --ntasks=32
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=48
##SBATCH -c 4
#SBATCH --export=ALL
#SBATCH --time=24:00:00
##SBATCH -A TG-CCR110036
##SBATCH -A TG-EAR160023
#SBATCH -A TG-EAR170019
#SBATCH --mail-user=shijia1019@gmail.com
#SBATCH --mail-type=all

export OMP_NUM_THREADS=2
export MV2_ENABLE_AFFINITY=0

cd /work/04149/tg834482/stampede2/NormalModes/IOnew/PNMsG_0.0
source setEnv
cd trGbinE0
ibrun ./plmvcg_Stampede2_itl18.out

