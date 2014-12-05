#!/bin/bash

mol=c6h6_rimp2
out=$mol.out
err=$mol.err

t2sub -q S -l select=2:mpiprocs=1:ncpus=12:mem=1gb -l place=scatter -o $out -e $err ./Tsubame2-que.sh
