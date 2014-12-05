#!/bin/bash -x
#
#PJM -o "h2o_10_rimp2.out"
#PJM -e "h2o_10_rimp2.err"
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "node=2"
#PJM --rsc-list "elapse=0:10:00"
#PJM --mpi "shape=2"
#PJM --mpi "use-rankdir"
#PJM --stg-transfiles all
#PJM --stgin "rank=* ./rimp2.exe %r:./"
#PJM --stgin "rank=* ./h2o_10_rimp2.inp %r:./INPUT"
#PJM --stgin "rank=* ./h2o_10_rimp2.Basis %r:./"
#PJM --stgin "rank=* ./h2o_10_rimp2.Basis_RIC %r:./"
#PJM --stgin "rank=* ./h2o_10_rimp2.Geom %r:./"
#PJM --stgin "rank=* ./h2o_10_rimp2.MO %r:./"
#PJM --stgin "rank=* ./h2o_10_rimp2.OrbEne %r:./"
#PJM --stgin "rank=* ./h2o_10_rimp2.SCF_Info %r:./"
#PJM --stgin "rank=* ./h2o_10_rimp2.TotEne %r:./"
#PJM -s
#

. /work/system/Env_base

export FLIB_FASTOMP=FALSE
export FLIB_CNTL_BARRIER_ERR=FALSE

mpiexec ./rimp2.exe
