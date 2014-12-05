#!/bin/bash
#BSUB -J BSUB-C6H6
#BSUB -o BSUB-C6H6-%J
#BSUB -n 2
#BSUB -R "span[ptile=2]"
#BSUB -x
# How to submit LSF job: bsub < script.sh
module load intel impi mkl
module list
set -x
date
hostname
NTCHEM_DIR=${HOME}/ntchem/ntchem-rimp2
MOLECULE=c6h6
DATA_DIR=${NTCHEM_DIR}/tests/${MOLECULE}
cd ${DATA_DIR}; if [ $? != 0 ] ; then echo '@@@ Directory error @@@'; exit; fi
pwd
NPROCS=1
export OMP_NUM_THREADS=8
export I_MPI_DEBUG=1
time mpirun -np $NPROCS ./rimp2.exe | tee ${MOLECULE}_rimp2.out
ls -go
pwd
python ./check.py

