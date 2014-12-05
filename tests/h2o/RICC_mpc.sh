#!/bin/sh

#----- qsub option -----#
#MJS: -mpc
#MJS: -proc 2
#MJS: -thread 8
#MJS: -time 0:15:00
#MJS: -mem 9600mb
#MJS: -cwd
#MJS: -comp intel
#MJS: -o h2o_10_rimp2.out
#MJS: -e h2o_10_rimp2.err

#------- FTL command -------#

#<BEFORE>
# rimp2.exe, INPUT, h2o_10_rimp2.Basis, h2o_10_rimp2.Basis_RIC, h2o_10_rimp2.Geom, 
# h2o_10_rimp2.MO, h2o_10_rimp2.OrbEne, h2o_10_rimp2.SCF_Info, h2o_10_rimp2.TotEne
#</BEFORE>

#----- Program execution -----#

export KMP_LIBRARY=turnaround
export KMP_BLOCKTIME=infinite
export KMP_AFFINITY=compact

mpirun ./rimp2.exe
