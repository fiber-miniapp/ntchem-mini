#
#                          ____ Module Makefile ____
#
include ../config/makeconfig

all   :  mod            lib            exe
#        -------------- -------------- -------------- -------------- -------------- 
mod   :  eparallel_mod
#        -------------- -------------- -------------- -------------- -------------- 
lib   :  eutil_lib
#        -------------- -------------- -------------- -------------- -------------- 
exe   :  emp2
#        -------------- -------------- -------------- -------------- -------------- 
clean :  clmod          cllib          clexe
#        -------------- -------------- -------------- -------------- -------------- 
clmod :  cparallel_mod
#        -------------- -------------- -------------- -------------- -------------- 
cllib :  cutil_lib
#        -------------- -------------- -------------- -------------- -------------- 
clexe :	 cmp2
#        -------------- -------------- -------------- -------------- -------------- 

eparallel_mod : 
	cd parallel_mod; $(MAKE); cd ..
cparallel_mod : 
	cd parallel_mod; $(MAKE) clean; cd ..

eutil_lib : 
	cd util_lib; $(MAKE); cd ..
cutil_lib : 
	cd util_lib; $(MAKE) clean; cd ..

emp2 : 
	cd mp2; $(MAKE); cd ..
cmp2 : 
	cd mp2; $(MAKE) clean; cd ..
