# 2007/06/11 Takahito Nakajima

#include ./config/makeconfig

all: time_start check_compiler $(BIN) $(LIB) make_src copy_bin time_end

time_start:
	@$(RM) .timestamp_start
	@printf "  Start: " > .timestamp_start
	@date >> .timestamp_start

time_end:
	@$(RM) .timestamp_end
	@printf " Finish: " > .timestamp_end
	@date >> .timestamp_end
	@echo 
	@echo " Make report"
	@echo " -----------"
	@cat .timestamp_start
	@cat .timestamp_end
	@echo 

make_src:
	@echo
	@echo '+++++ Make +++++'
	@echo
	(cd ./src; $(MAKE))

copy_bin:
	find ./tests -name "*" -type d -exec $(CP) $(BIN)/*.exe {} \;

$(BIN):
	$(MKDIR) -p $(BIN)

$(LIB):
	$(MKDIR) -p $(LIB)

$(INCLUDE):
	$(MKDIR) -p $(INCLUDE)

clean:
	(cd ./src; $(MAKE) clean; cd ../)
	find ./tests -name "*.exe" -exec $(RM) {} \;

check_compiler:
	@printf '\n Checking paths for compilers.\n which %s : ' $(CC); which $(CC); \
	 printf ' which %s : ' $(CXX) ; which $(CXX); \
	 printf ' which %s : ' $(F77C) ; which $(F77C); \
	 printf ' which %s : ' $(F90C) ; which $(F90C)

distclean : veryclean

veryclean : 
#	(cd ./src; $(MAKE) distclean)
	find . \( -name '*.o' -o \
                  -name '*.a' -o \
                  -name '*.pyc' -o \
                  -name '*.mod' -o \
                  -name '*.exe' -o \
                  -name '*.out' -o \
                  -name '*.err' -o \
                  -name 'core' \) -exec $(RM) {} \; ;\
	$(RM) -f .config.sed config.history;\
	$(RM) -f .timestamp_start .timestamp_end ;\
	$(RM) -rf bin ;\
	$(RM) -rf lib ;\
	$(RM) -f config/config.sh config/makeconfig
#

test : 
	cd ./tests/c6h6 ;\
	rm c6h6_rimp2.out ;\
	./rimp2.exe > c6h6_rimp2.out ;\
	python ./check.py ;\
	rm c6h6_rimp2.out
