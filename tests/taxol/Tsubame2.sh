#!/bin/bash

mol=taxol_rimp2_cc-pvdz
out=$mol.out
err=$mol.err

t2sub -q S -l select=32:mpiprocs=1:ncpus=12:mem=16gb -l place=scatter -o $out -e $err ./Tsubame2-que.sh
