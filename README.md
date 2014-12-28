
NTChem-mini
=============

* version: 1.1 (based on NTChem/RI-MP2 dated 2014/10/31)
* update: 2014/12/28 (bug fix rimp2_rmp2energy_incore_v_mpiomp.F90)
* contact: miniapp@riken.jp

About NTChem and NTChem-mini
--------------------------------

NTChem is a high-performance software package for the molecular electronic
structure calculation for general purpose on the K computer.
It is a comprehensive new software of ab initio quantum chemistry made in AICS
from scratch.
NTChem contains not only standard quantum chemistry approaches but our own
original approaches. NTChem is expected to be a useful tool in various
computational studies for large and complicated molecular systems.

For more details, see NTChem official webpage
<http://labs.aics.riken.jp/nakajimat_top/ntchem_e.html> .

NTChem-mini includes a subset of NTChem known as NTChem/RI-MP2 that includes
the efficient scheme to account for the electron correlations based on the
resolution of identity second-order Møller–Plesset perturbation (RI-MP2)
method that enables the computation of the large molecular systems such as
nano-molecules and biological molecules.

NTChem-mini is a packaged version of NTChem/RI-MP2, including a suite of
test data, installation write up, and testing procedures.

For detail explanation of NTChem/RI-MP2, refer to the paper: M. Katouda
and T. Nakajima, "MPI/OpenMP Hybrid Parallel Algorithm of Resolution of
Identity Second-Order Møller–Plesset Perturbation Calculation for
Massively Parallel Multicore Supercomputers",
J. Chem. Theory Comput. 9, 5373 (2013).

Contact point for NTChem/RI-MP2: Dr. Michio Katouda <mkatouda@riken.jp>



Installation
------------

The installation steps of NTChem-mini on Linux platform are explained below.
Confirm first that the prerequisite software is available on the 
installing platform.

+ Prerequisite software:
  - Fortran compiler
  - MPI library,
  - BLAS library
  - LAPACK library

####step 1. Obtain the distribution package.

Obtain the NTChem-mini package file "ntchem-mini-1.1.tar.gz" from the Fiber
distribution site.
+ http://fiber-miniapp.github.io/

On an appropriate directory, extract the contents using tar command.

    $ tar -zxf ntchem-mini-1.1.tar.gz
    $ cd ./ntchem-mini
    $ TOP_DIR=${PWD}

The variable TOP_DIR is used to point the extracted top directory
hereafter .
This readme file (README.md) is included in the ${TOP_DIR} directory.


####step 2.  Choose config_mine for the installing platform.

In the ${TOP_DIR}/platforms subdirectory, there are several
config_mine.${TYPE} files where TYPE value represents the type of
the compiler or the system type of the installing platform.
These are small text script files which set the minimum path variables
for compilation.
Choose one of these config_mine.${TYPE} files that matches the installing
platform, and copy it as config_mine.
If the installing platform does not match any of the provided TYPE,
then create a new config_mine file, based on one of the included
config_mine.* files.

    $ ls -go platforms/config_mine.*
    -rwxr-xr-x 1  162 Nov  7 00:00 platforms/config_mine.K
    -rwxr-xr-x 1  140 Nov  7 00:00 platforms/config_mine.RICC
    -rwxr-xr-x 1  124 Nov  7 00:00 platforms/config_mine.Tsubame2
    -rwxr-xr-x 1  125 Nov  7 00:00 platforms/config_mine.c0
    -rwxr-xr-x 1  246 Nov  7 00:00 platforms/config_mine.gfortran
    -rwxr-xr-x 1  148 Nov  7 00:00 platforms/config_mine.intel
    -rwxr-xr-x 1  156 Nov  7 00:00 platforms/config_mine.pgi
    $ TYPE=intel
    $ cp config_mine.${TYPE} config_mine

Note that NTChem requires BLAS and LAPACK math libraries.
Typically, each TYPE has its own recommended BLAS and LAPACK libraries.
For example, K computer is best coupled with its native BLAS/LAPACK,
Intel with MKL, PGI with ACML, Gfortran(GNU) with ATLAS, etc.
Note that ATLAS library included in the common Linux distribution
may not perform well, and that freshly compiled ATLAS may perform
better on the target platform.  See comments in config_mine.gfortran.


####step 3.  Compile and build the executable program

In the ${TOP_DIR} directory, choose or edit the config_mine file as above,
execute config_mine script and make command.

    $ make clean
    $ ./config_mine
    $ make
    $ ls -go bin/rimp2.exe 
    -rwxr-xr-x 1 1556942 Oct 10 10:26 bin/rimp2.exe
    $ 

Compiling time of the source program should be short, i.e. within a minute.
After the successful compilation, the executable file "rimp2.ex" should be
created in the same directory.

#### Installation Example: Installing on K computer

Below is an example of compiling on K computer.

    $ make clean
    $ cp platforms/config_mine.K config_mine
    $ ./config_mine
    $ make
    $ ls -go rimp2.ex
	-rwxr-xr-x 1 1140539 Sept 22 13:35 2014 rimp2.ex


#### Installation Example: Intel platform with modulefiles

Below is an example of compiling with Intel compiler, Intel MPI and Intel MKL
with module environment.

    $ module load intel impi mkl
    $ cp platforms/config_mine.intel config_mine
    $ make clean
    $ ./config_mine
    $ make
    $ ls -go bin/rimp2.exe 
    -rwxr-xr-x 1 1556942 Oct 10 10:26 bin/rimp2.exe

#### Installation Example: Intel platform without modulefiles

If the installing platform does not support modulefiles,
then manually set the environment variables for using the compiler software.
Intel provides the initialization files for each of the components, and it
is usually simple enough to load the intel environment variables by reading
those initialization file.

Below is an example of compiling with Intel compiler, Intel MPI and Intel MKL
using the initialization files for sh/bash.
Note that the path variables
INTEL_DIR, I_MPI_ROOT and MKLROOT must be changed according to the
actual path on the installing platform.

    $ export INTEL_DIR=/usr/local/intel/composer_xe_2013
    $ source ${INTEL_DIR}/bin/compilervars.sh intel64
    $ export I_MPI_ROOT=/usr/local/intel/impi/4.1.0.024
    $ export I_MPI_F90=ifort
    $ export I_MPI_F77=ifort
    $ export I_MPI_CC=icc
    $ export I_MPI_CXX=icpc
    $ source ${I_MPI_ROOT}/bin64/mpivars.sh
    $ export MKLROOT=/usr/local/intel/composer_xe_2013/mkl
    $ source ${MKLROOT}/bin/mklvars.sh intel64
    $ make clean
    $ cp platforms/config_mine.intel config_mine
    $ ./config_mine
    $ make
    $ ls -go bin/rimp2.exe 
    -rwxr-xr-x 1 1556942 Oct 10 10:26 bin/rimp2.exe


Running the small test jobs
--------------------
In the ${TOP_DIR}/tests directory , there are several subdirectories
containing the test data along with the job script.

    $ cd ${TOP_DIR}/tests
	$ ls -go
	drwxr-xr-x 2    4096 Nov  5 17:53 c6h6
	drwxr-xr-x 2    4096 Nov  5 17:53 h2o
	drwxr-xr-x 2    4096 Nov  5 17:53 taxol

Each subdirectory contains a group of input files.
It also contains the job script examples for several platforms.
+ input files:
  - ${DataName}.Basis
  - ${DataName}.Basis_RIC
  - ${DataName}.Geom
  - ${DataName}.MO
  - ${DataName}.OrbEne
  - ${DataName}.SCF_Info
  - ${DataName}.TotEne
  - ${DataName}.inp
+ job scripts:
  - ${PLATFORM}.sh

The ${DataName}.inp is also aliases as "INPUT" file.
The string value of ${DataName} is given as the namelist variable of &prefix,
such as Name='taxol_rimp2_cc-pvdz' in "INPUT" file.
The string value of ${PLATFORM} represents the type of running platform
or maybe the name of job manager.

The parallel program model used in NTChem-mini is MPI and OpenMP.

Below is the list of subdirectories showing their data type, the computing
time and the required memory size on a reference system.

|directory| data type           |elapse time|memory size|# of nodes|
|---------|---------------------|-----------|-----------|----------|
| c6h6    | benzene RI-MP2/SVP  | 3 seconds |     7 MB  |         1|
| h2o     | H2O 10mer RI-MP2/cc-pVTZ| 8 seconds|  7 MB  |         1|
| taxol   | taxol RI-MP2/cc-pVDZ| 10 minutes|    30 MB  |        32|

#### Example: Running the test job on general Linux platform.
Below is an example of running NTChem-mini interactively using the taxol
test data on general Linux platform which supports the module environment.

    $ module load intel impi mkl
    $ cd ${TOP_DIR}/tests/taxol
    $ NPROCS=8
    $ export OMP_NUM_THREADS=2
    $ eval `grep Name INPUT | tr [,] [\ ]`
    $ time mpirun -np $NPROCS ./rimp2.exe | tee ${Name}.out

If the module environment is not supported on the plarform, then the
environment variables can be set explicitly as shown in the
"Installation Example:" section above.

For batch jobs, the job script file with corresponding directives
should be written.
Example job script files for some platforms can be found under the
${TOP_DIR}/tests/* subdirectories.

    $ cd ${TOP_DIR}/tests/taxol
    $ ls -go *sh
    -rw-r--r-- 1 1394 Sep 16 17:31 AICS-common.sh
    -rw-r--r-- 1  862 Sep 16 17:31 Kcomputer.sh
    -rw-r--r-- 1  627 Sep 16 17:31 RICC_mpc.sh
    -rwxr-xr-x 1  168 Sep 16 17:31 Tsubame2.sh



#### Example: Running the test job on K computer
Below in an example job script for K computer with necessary file staging
directives.

	pjsub Kcomputer.sh


#### Example: Running the test job on Intel platform with LSF

Below is an example of a batch job script with LSF directives 

    $ cat job.sh
    #BSUB -J BSUB-TAXOL
    #BSUB -o BSUB-TAXOL-%J
    #BSUB -n 32
    #BSUB -R "span[ptile=2]"
    #BSUB -x
    module load intel impi mkl
    set -x
    date
    MODEL=taxol
    DATA_DIR=${HOME}/ntchem/ntchem-rimp2/tests/${MODEL}
    cd ${DATA_DIR}; if [ $? != 0 ] ; then echo '@@@ Directory error @@@'; exit; fi
    eval `grep Name INPUT | tr [,] [\ ]`
    NPROCS=32
    export OMP_NUM_THREADS=8
    time mpirun -np $NPROCS ./rimp2.exe | tee ${Name}.out
    $
    $ bsub < job.sh
    $


#### Verifying the computed results

The data directory contains a small Python script file "check.py" to verify 
the computed results recorded in  ${Name}.out file.
Run this script file and confirm the output message.

    $ python ./check.py
    The result is OK



Running the large size performance tests 
-----------------------------

The default test data included in the above NTChem-mini package are limited
to small molecular sets in order to make the package compact.
For running the performance tests, it is recommended to use larger molecular
data.
A separate file containing the molecular data suitable for such performance
test has been made available for downloading from the following site.

http://hpci-aplfs.aics.riken.jp/fiber/ntchem-mini/ntchem-data-perf.tar.gz

#### Setting up the performance test data

After the previous Installation steps are completed, download the performance
data file ntchem-data-perf.tar.gz from the above URL.
Save it on the ${TOP_DIR} directory and extract the performance data files.
There should be new data file directories under tests_perf/ .

    $ cd ${TOP_DIR}
	$ tar -zxf ntchem-data-perf.tar.gz
	$ ls tests_perf/
	c54h18_cc-pvdz  c54h18x2_cc-pvdz  c60atc60h28
	c54h18_cc-pvtz  c54h18x2_cc-pvtz

+ directory name: data type
  - c54h18_cc-pvdz: nanographene C54H18 RI-MP2/cc-pVDZ
  - c54h18x2_cc-pvdz: nanographene dimer (C54H18)2 RI-MP2/cc-pVDZ
  - c54h18_cc-pvtz: nanographene C54H18 RI-MP2/cc-pVTZ
  - c54h18x2_cc-pvtz: nanographene dimer (C54H18)2 RI-MP2/cc-pVTZ
  - c60atc60h28: buckycatcher C60@C60H28 RI-MP2/cc-pVTZ
+ buckycatcher data _c60atc60h28_ is most appropriate for the
performance tests.

|directory   | data type           |elapse time|memory use|# of nodes|
|------------|---------------------|-----------|-----------|----------|
| c60atc60h28| buckycatcher C60@C60H28 RI-MP2/cc-pVTZ| 2 hours| 10 GB| 32|


#### Running the performance tests
Like the default test data included in the NTChem-mini package,
there is a set of job script files in the performance test data directory
covering the typical platforms.
Choose the one that matches the testing plarform.
Edit it, if necessary, in the same manner as was done for the
default test data.
Execute it as a batch job or as a background job.
The c60atc60h28 job takes about 2 hours using 32 nodes on K computer.


##### strong scaling tests

The strong scaling tests for NTChem-mini in general can be done by just
changing the number of MPI processes per job.
For example, a simple script such as below will define the strong scaling runs
for 8/16/32 processes.

    $ for NPROCS in 8 16 32
    $ do
    $ time -p mpirun -np $NPROCS ./rimp2.exe
    $ done

In this case, each execution of rimp2.exe reads the same "INPUT" file.

##### weak scaling tests

It is not straight forward to define the test data for evaluating the
weak scaling performance, since the computational workload of NTChem-mini
increases at N^5 order accroding to the number of atoms N.
In order to prepare for the weak scaling performance test, specific testing
configuration should be defined and the corresponding input data files
must be generated. Such configuration can be defined through
the discussion with the authors of NTChem.


License term
-----------------------------

The License term of NTChem-mini is provided in the BSD 2-Clause License.
Please refer to "LICENSE" file included in the NTChem-mini package.

Copyright (c) 2012-2014, Computational Molecular Science Research Team,
RIKEN Advanced Institute for Computational Science

