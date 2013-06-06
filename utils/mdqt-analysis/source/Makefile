#######################################################################
#
# Alexander V. Soudackov, Penn State, March 15, 2011
#
#######################################################################

# Set the F90 Compiler

F90 = ifort
F90FLAGS = -c -O3 -xHOST -ipo -zero -static -prec-div -fstack-security-check -fstack-protector -fp-model strict
#F90FLAGS = -c -g -check

# Set the F90 Loader
LD = $(F90) -o
LDFLAGS = -fast -ftz -zero
#LDFLAGS = -g -check

SHELL=/bin/sh

.SUFFIXES: .f90 .mod .o

OBJS = analyze_trajectories.o module_string_utilities.o module_sorting.o module_timers.o

BINDIR = ../bin

EXE = $(BINDIR)/analyze_trajectories.bin

all:	analyze_trajectories

analyze_trajectories: $(OBJS)
	$(LD) $(EXE) $(OBJS) $(LDFLAGS)

clean:
	-rm -f $(OBJS) *.mod

.f90.o: 
	$(F90) $(F90FLAGS) $*.f90

analyze_trajectories.o: module_string_utilities.o module_sorting.o module_timers.o
