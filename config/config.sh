#
export NTQC_TOP=/home/ra000004/mikami/ntchem/ntchem-mini-1.1
export HOSTTYPE=linux64_mpif90_omp_intel_proto
export LAPACK=
export BLAS=
export ATLAS=-mkl
export SCRATCH=/home/mikami/scr/ntchem
export PARALLEL=mpiomp
#
export TARGET=LINUX64
unset USE_MPI

# if you want to use MPICH, you can set the environmental variables as
# follos (see ./GA/README)
# 
# export MPI_USE=yes
# export MPI_INCLUDE=/usr/include
# export MPI_LIB=/usr/lib
# export LIBMPI=-lmpi

export LARGE_FILES=yes

