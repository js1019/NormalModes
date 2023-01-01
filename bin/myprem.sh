#!/bin/bash

#PBS -q gmig
#PBS -l nodes=2:ppn=16
#PBS -N prem1M
#PBS -l walltime=1:00:00
#PBS -V
#PBS -n
#PBS -M shijia1019@gmail.com
#PBS -m ea


##export OMP_NUM_THREADS=1
##export MV2_ENABLE_AFFINITY=0

cd /home/shi126/NormalModes/CNT_evsl/PNM_3.0
source setEnv
cd bin/
mpirun ./plmvcg_conte.out
