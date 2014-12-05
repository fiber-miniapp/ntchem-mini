#!/bin/bash

export PATH=/usr/apps/openmpi/1.4.2/intel/bin:$PATH
export LD_LIBRARY_PATH=/usr/apps/openmpi/1.4.2/intel/lib:$LD_LIBRARY_PATH

cd ${PBS_O_WORKDIR}

mpirun -n 32 -hostfile $PBS_NODEFILE ./Tsubame2-mpi_go_opn.sh
