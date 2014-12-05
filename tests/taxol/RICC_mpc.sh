#!/bin/sh

#----- qsub option -----#
#MJS: -mpc
#MJS: -proc 32
#MJS: -thread 8
#MJS: -time 0:10:00
#MJS: -mem 9600mb
#MJS: -cwd
#MJS: -comp intel
#MJS: -o taxol_rimp2_cc-pvdz.out
#MJS: -e taxol_rimp2_cc-pvdz.err

#------- FTL command -------#

#<BEFORE>
# rimp2.exe, INPUT, taxol_rimp2_cc-pvdz.Basis,
# taxol_rimp2_cc-pvdz.Basis_RIC, taxol_rimp2_cc-pvdz.Geom, 
# taxol_rimp2_cc-pvdz.MO, taxol_rimp2_cc-pvdz.OrbEne,
# taxol_rimp2_cc-pvdz.SCF_Info, taxol_rimp2_cc-pvdz.TotEne
#</BEFORE>

#----- Program execution -----#

export KMP_LIBRARY=turnaround
export KMP_BLOCKTIME=infinite
export KMP_AFFINITY=compact

mpirun ./rimp2.exe
