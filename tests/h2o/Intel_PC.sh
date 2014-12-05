#!/bin/bash
#BSUB -J BSUB-H2O
#BSUB -o BSUB-H2O-%J
#BSUB -n 2
#BSUB -R "span[ptile=2]"
#BSUB -x
# How to submit LSF job: bsub < script.sh
module load intel impi mkl
module list
set -x
date
hostname
cat /etc/*-release	# /etc/system-release

NTCHEM_DIR=${HOME}/ntchem/ntchem-rimp2
MODEL=h2o

DATA_DIR=${NTCHEM_DIR}/tests/${MODEL}
cd ${DATA_DIR}; if [ $? != 0 ] ; then echo '@@@ Directory error @@@'; exit; fi
eval `grep Name INPUT | tr [,] [\ ]`

NPROCS=2
export OMP_NUM_THREADS=8
export I_MPI_DEBUG=1
time mpirun -np $NPROCS ./rimp2.exe | tee ${Name}.out
ls -go
pwd
python ./check.py

