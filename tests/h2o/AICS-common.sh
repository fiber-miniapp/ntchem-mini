#!/bin/sh
#
#$ -S /bin/sh
#$ -cwd
#$ -pe mpi 2
#$ -o h2o_10_rimp2.out
#$ -e h2o_10_rimp2.err

export RLIMIT_MEMLOCK=65536
ulimit -l 65536
. /opt/intel/ics/2011.05.23/ictvars.sh
export I_MPI_HYDRA_BOOTSTRAP=rsh
export I_MPI_HYDRA_BOOTSTRAP_EXEC=/usr/bin/rsh

nprocs=2
MOL=h2o_10_rimp2
curdir=`pwd`
bindir=$curdir
wrkdir=/ssd/katouda/ntchem

export OMP_NUM_THREADS=12

mpiexec.hydra -hostfile $TMPDIR/machines -n $nprocs -perhost 1 cp $curdir/$MOL.inp $wrkdir/INPUT
mpiexec.hydra -hostfile $TMPDIR/machines -n $nprocs -perhost 1 cp $curdir/$MOL.Geom $wrkdir/$MOL.Geom
mpiexec.hydra -hostfile $TMPDIR/machines -n $nprocs -perhost 1 cp $curdir/$MOL.Basis $wrkdir/$MOL.Basis
mpiexec.hydra -hostfile $TMPDIR/machines -n $nprocs -perhost 1 cp $curdir/$MOL.Basis_RIC $wrkdir/$MOL.Basis_RIC
mpiexec.hydra -hostfile $TMPDIR/machines -n $nprocs -perhost 1 cp $curdir/$MOL.SCF_Info $wrkdir/$MOL.SCF_Info
mpiexec.hydra -hostfile $TMPDIR/machines -n $nprocs -perhost 1 cp $curdir/$MOL.TotEne $wrkdir/$MOL.TotEne
mpiexec.hydra -hostfile $TMPDIR/machines -n $nprocs -perhost 1 cp $curdir/$MOL.OrbEne $wrkdir/$MOL.OrbEne
mpiexec.hydra -hostfile $TMPDIR/machines -n $nprocs -perhost 1 cp $curdir/$MOL.MO $wrkdir/$MOL.MO
cd $wrkdir
mpiexec.hydra -hostfile $TMPDIR/machines -n $nprocs -perhost 1 $bindir/rimp2.exe
mpiexec.hydra -hostfile $TMPDIR/machines -n $nprocs -perhost 1 rm $MOL.* INPUT
