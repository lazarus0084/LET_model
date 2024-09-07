#Makefile for ILLITHERM-FV

PROGRAM = iitpave2
SOURCE = iitpave2.cc readinputs.cc LET.cc XipWip.cc besselroots.cc MATLAB_functions.cc lookup_gauss.cc arb_func.cc LET_response.cc arb_func_interp.cc
OBJECTS = iitpave2.o readinputs.o LET.o   XipWip.o besselroots.o MATLAB_functions.o lookup_gauss.o arb_func.o LET_response.o arb_func_interp.o
ARCHIVE = iitpave2.tar

CC = g++
CFLAGS = -g -lm -I.

DEPS = iitpave2.h
OPTIONS = -std=c++11

%.o: %.cpp $(DEPS)
	$(CC) $(OPTIONS) -c -o $@ $< $(CFLAGS)

$(PROGRAM):	$(OBJECTS)
	$(CC) $(OPTIONS) -o $(PROGRAM).out $(OBJECTS) $(CFLAGS)
	
clean: 
	rm -f $(OBJECTS) $(PROGRAM).out *.tar *.plt
	
rebuild:
	clean $(PROGRAM)
	
archive:
	@echo Creating code archive tar-file $(ARCHIVE) ...
	tar cf $(ARCHIVE) $(SOURCE) Makefile
	@ls -l $(ARCHIVE)

help:
	@echo Try:
	@echo make $(PROGRAM) .... to build the program named $(PROGRAM)
	@echo make clean .... to clean up, removing object files and program $(PROGRAM)
	@echo make archive .... to make an archive tar file you can transfer or submit