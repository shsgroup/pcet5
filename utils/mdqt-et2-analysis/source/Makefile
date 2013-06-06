#######################################################################
#
# Alexander V. Soudackov, Penn State, July 7, 2012
#
#######################################################################

# Set the F90 Compiler

F90 = ifort
F90FLAGS = -c -O3 -xHOST -ipo -zero -static -prec-div -fstack-security-check -fstack-protector -fp-model strict
#F90FLAGS = -c -g -check

# Set the F90 Loader
LD = $(F90) -o
LDFLAGS = -O3 -xHOST -ipo -zero -static -prec-div -fstack-security-check -fstack-protector -fp-model strict
#LDFLAGS = -L/opt/intel/Compiler/11.1/075/mkl/lib/em64t -Wl,--start-group /opt/intel/Compiler/11.1/075/mkl/lib/em64t/libmkl_intel_lp64.a /opt/intel/Compiler/11.1/075/mkl/lib/em64t/libmkl_core.a /opt/intel/Compiler/11.1/075/mkl/lib/em64t/libmkl_sequential.a -Wl,--end-group -lpthread -lm

#LDFLAGS = -g -check

SHELL=/bin/sh

.SUFFIXES: .f90 .mod .o

OBJS = analyze_et2_trajectories.o module_string_utilities.o module_sorting.o module_timers.o real_precision.o lpkbls.o odrpack95.o module_marcus.o

BINDIR = ../bin

EXE = $(BINDIR)/analyze_et2_trajectories.bin

all:	analyze_et2_trajectories

analyze_et2_trajectories: $(OBJS)
	$(LD) $(EXE) $(OBJS) $(LDFLAGS)

clean:
	-rm -f $(OBJS) *.mod

.f90.o: 
	$(F90) $(F90FLAGS) $*.f90

.f.o: 
	$(F90) $(F90FLAGS) $*.f

analyze_et2_trajectories.o: module_string_utilities.o module_sorting.o real_precision.o lpkbls.o odrpack95.o module_timers.o module_marcus.o
