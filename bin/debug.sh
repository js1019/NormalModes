#!/bin/bash

#PBS -q debug
#PBS -l nodes=1:ppn=16
#PBS -N Debug
#PBS -l walltime=0:30:00
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
