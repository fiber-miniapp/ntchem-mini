#!/bin/bash -x
#
#PJM -o "taxol_rimp2_cc-pvdz.out"
#PJM -e "taxol_rimp2_cc-pvdz.err"
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "node=32"
#PJM --rsc-list "elapse=0:10:00"
#PJM --mpi "shape=32"
#PJM --mpi "use-rankdir"
#PJM --stg-transfiles all
#PJM --stgin "rank=* ./rimp2.exe %r:./"
#PJM --stgin "rank=* ./taxol_rimp2_cc-pvdz.inp %r:./INPUT"
#PJM --stgin "rank=* ./taxol_rimp2_cc-pvdz.Basis %r:./"
#PJM --stgin "rank=* ./taxol_rimp2_cc-pvdz.Basis_RIC %r:./"
#PJM --stgin "rank=* ./taxol_rimp2_cc-pvdz.Geom %r:./"
#PJM --stgin "rank=* ./taxol_rimp2_cc-pvdz.MO %r:./"
#PJM --stgin "rank=* ./taxol_rimp2_cc-pvdz.OrbEne %r:./"
#PJM --stgin "rank=* ./taxol_rimp2_cc-pvdz.SCF_Info %r:./"
#PJM --stgin "rank=* ./taxol_rimp2_cc-pvdz.TotEne %r:./"
#PJM -s
#

. /work/system/Env_base

export FLIB_FASTOMP=FALSE
export FLIB_CNTL_BARRIER_ERR=FALSE

mpiexec ./rimp2.exe
