# Makefile for parallel_mod
#                     30 Aug 2012 Modified by Yusuke Ootani
#
include ../../config/makeconfig

modules = \
mpi_module.o \

all: modcomp

modcomp : $(modules)

clean :
	$(RM) *.o *.mod

veryclean :
	$(RM) *.o *.mod 

.SUFFIXES:.f90
.f90.o: 
	$(F90C) $(F90FLAGS) $(OMPFLAGS) -c $<
