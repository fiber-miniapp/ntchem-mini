# Makefile for util_lib
#                     03 Sep 2007 Written  by Takahito Nakajima
#                     28 Sep 2012 Modified by Michio Katouda
#
include ../../config/makeconfig

objects_common = \
util_linout.o \
util_matout.o \
util_transchar.o \
util_tr1to2.o \
util_lowtrmzero.o \

objects_mpi = \
util_abortmpi.o \
util_initmpi.o \
util_finmpi.o \

objec77 = \

objects = $(objects_common) $(objects_mpi)

all: objcomp

objcomp : $(objects) $(objec77)
	$(RM) libutil.a
	$(AR) $(ARFLAGS) libutil.a $(objects) $(objec77)
	$(RANLIB) libutil.a

clean :
	$(RM) *.o *.mod *.a

veryclean :
	$(RM) *.o *.mod *.a

.SUFFIXES:.f90
.f90.o: 
	$(F90C) $(F90FLAGS) $(OMPFLAGS) -I$(NTQC_TOP)/src/parallel_mod -c $<

.SUFFIXES:.F90
.F90.o:
	$(F90C) $(F90FLAGS) $(OMPFLAGS) -I$(NTQC_TOP)/src/parallel_mod -c $<

.SUFFIXES:.f
.f.o: 
	$(F77C) $(F77FLAGS) -c $<
