#!/bin/sh

#----- qsub option -----#
#MJS: -mpc
#MJS: -proc 2
#MJS: -thread 8
#MJS: -time 0:10:00
#MJS: -mem 9600mb
#MJS: -cwd
#MJS: -comp intel
#MJS: -o c6h6_rimp2.out
#MJS: -e c6h6_rimp2.err

#------- FTL command -------#

#<BEFORE>
# rimp2.exe, INPUT, c6h6_rimp2.Basis, c6h6_rimp2.Basis_RIC, c6h6_rimp2.Geom, 
# c6h6_rimp2.MO, c6h6_rimp2.OrbEne, c6h6_rimp2.SCF_Info, c6h6_rimp2.TotEne
#</BEFORE>

#----- Program execution -----#

export KMP_LIBRARY=turnaround
export KMP_BLOCKTIME=infinite
export KMP_AFFINITY=compact

mpirun ./rimp2.exe
