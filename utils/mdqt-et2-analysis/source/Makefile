#######################################################################
#
# Alexander V. Soudackov, Penn State, July 7, 2012
#
#######################################################################

# Set the F90 Compiler

F90 = ifort
F90FLAGS = -c -fast -ftz -zero
#F90FLAGS = -c -g -check

# Set the F90 Loader
LD = $(F90) -o
LDFLAGS = -fast -ftz -zero
#LDFLAGS = -g -check

SHELL=/bin/sh

.SUFFIXES: .f90 .mod .o

OBJS = analyze_et2_trajectories.o module_string_utilities.o module_sorting.o module_timers.o module_marcus.o

BINDIR = ../bin

EXE = $(BINDIR)/analyze_et2_trajectories.bin

all:	analyze_et2_trajectories

analyze_et2_trajectories: $(OBJS)
	$(LD) $(EXE) $(OBJS) $(LDFLAGS)

clean:
	-rm -f $(OBJS) *.mod

.f90.o: 
	$(F90) $(F90FLAGS) $*.f90

analyze_et2_trajectories.o: module_string_utilities.o module_sorting.o module_timers.o module_marcus.o
