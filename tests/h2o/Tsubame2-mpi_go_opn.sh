#!/bin/bash

export PATH=/usr/apps/openmpi/1.4.2/intel/bin:$PATH
export LD_LIBRARY_PATH=/usr/apps/openmpi/1.4.2/intel/lib:$LD_LIBRARY_PATH

MOL=h2o_10_rimp2
curdir=`pwd`

cd /scr

if [ ! -d $OMPI_COMM_WORLD_RANK ] ; then
    mkdir $OMPI_COMM_WORLD_RANK
fi

cd $OMPI_COMM_WORLD_RANK

cp $curdir/rimp2.exe ./rimp2.exe
cp $curdir/$MOL.inp ./INPUT
cp $curdir/$MOL.Geom ./$MOL.Geom
cp $curdir/$MOL.Basis ./$MOL.Basis
cp $curdir/$MOL.Basis_RIC ./$MOL.Basis_RIC
cp $curdir/$MOL.SCF_Info ./$MOL.SCF_Info
cp $curdir/$MOL.TotEne ./$MOL.TotEne
cp $curdir/$MOL.OrbEne ./$MOL.OrbEne
cp $curdir/$MOL.MO ./$MOL.MO

./rimp2.exe

cd ../
rm -rf $OMPI_COMM_WORLD_RANK
