#!/bin/bash
#BSUB -J BSUB-TAXOL
#BSUB -o BSUB-TAXOL-%J
#BSUB -n 32
#BSUB -R "span[ptile=2]"
#BSUB -x
# How to submit LSF job: bsub < script.sh
module load intel impi mkl
module list
set -x
date

NTCHEM_DIR=${HOME}/ntchem/ntchem-rimp2
MODEL=taxol

DATA_DIR=${NTCHEM_DIR}/tests/${MODEL}
cd ${DATA_DIR}; if [ $? != 0 ] ; then echo '@@@ Directory error @@@'; exit; fi
eval `grep Name INPUT | tr [,] [\ ]`
NPROCS=32
export OMP_NUM_THREADS=8
export I_MPI_DEBUG=3
time mpirun -np $NPROCS ./rimp2.exe | tee ${Name}.out
ls -go
pwd
python ./check.py

