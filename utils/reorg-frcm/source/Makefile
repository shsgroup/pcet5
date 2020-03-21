#############################################################
#
# Makefile for the modular version of PCET 5.x
#
# Alexander Soudackov, The Pennsylvania State University
# July 1, 2010
#
#############################################################

.SUFFIXES: .o .mod .f90 .f

include Makefile.machine

VERSION = 5.5

MODDIR  = ./mod

ifeq ($(FC),gfortran)
   MODDIR = .
else
   MODDIR  = ./mod
endif

#ifeq ($(FC),openf90)
#   MODDIR = .
#endif

MODOPT  = $(MODFLAG)$(MODDIR)
OBJDIR  = ./obj
EXEDIR  = ../bin
PROGRAM = reorg_frcm_$(VERSION)_$(ARCH)
EXE     = $(EXEDIR)/$(PROGRAM).x
LIB     = lib$(PROGRAM).a

# rules

$(OBJDIR)/%.o : %.f90
	$(FC) -c $(F90FLAGS) $(NOBACKSLASH) -o $(OBJDIR)/$*.o $(MODOPT) $*.f90

$(OBJDIR)/banner.o : banner.f90
	$(FC) -c $(F90FLAGS) $(BACKSLASH) -o $(OBJDIR)/banner.o $(MODOPT) banner.f90

$(OBJDIR)/%.o : %.f
	$(FC) -c $(F77FLAGS) $(NOBACKSLASH) -o $(OBJDIR)/$*.o $(MODOPT) $*.f

$(MODDIR)/%.mod : $(OBJDIR)/module_%.o
	@if [ ! -f $@ ]; then rm $< ; $(MAKE) $< ; fi


SRCOBJ = $(OBJDIR)/reorg_frcm.o  \
         $(OBJDIR)/deinitmat.o   \
         $(OBJDIR)/ereorg.o      \
         $(OBJDIR)/banner.o      \
         $(OBJDIR)/primat.o      \
         $(OBJDIR)/setjob.o

SRCSRC = reorg_frcm.f90  \
         banner.f90      \
         deinitmat.f90   \
         ereorg.f90      \
         primat.f        \
         setjob.f90


MODOBJ  = $(OBJDIR)/module_control.o     \
          $(OBJDIR)/module_cst.o         \
          $(OBJDIR)/module_elcm.o        \
          $(OBJDIR)/module_elmnts.o      \
          $(OBJDIR)/module_frcm.o        \
          $(OBJDIR)/module_geometry.o    \
          $(OBJDIR)/module_geosol.o      \
          $(OBJDIR)/module_parsol.o      \
          $(OBJDIR)/module_keys.o        \
          $(OBJDIR)/module_pardim.o      \
          $(OBJDIR)/module_strings.o     \
          $(OBJDIR)/module_timers.o


MODMOD  = $(MODDIR)/control.mod     \
          $(MODDIR)/cst.mod         \
          $(MODDIR)/elcm.mod        \
          $(MODDIR)/elmnts.mod      \
          $(MODDIR)/frcm.mod        \
          $(MODDIR)/geometry.mod    \
          $(MODDIR)/geosol.mod      \
          $(MODDIR)/parsol.mod      \
          $(MODDIR)/keys.mod        \
          $(MODDIR)/pardim.mod      \
          $(MODDIR)/strings.mod     \
          $(MODDIR)/timers.mod

MODSRC  = module_control.f90     \
          module_cst.f90         \
          module_elcm.f90        \
          module_elmnts.f90      \
          module_frcm.f          \
          module_geometry.f90    \
          module_geosol.f90      \
          module_parsol.f90      \
          module_keys.f90        \
          module_pardim.f90      \
          module_strings.f90     \
          module_timers.f90

note:
	@echo
	@echo '=== MAKEFILE for reorg_frcm 5.x ==='
	@echo 'Options:'
	@echo 'make build: build '
	@echo 'make all: build and link '
	@echo 'make clean: remove object and module files '
	@echo 'Have fun !!!'
	@echo


build:  $(MODMOD)
	@echo " "
	@echo "===> building object files "
	@echo " "
	$(MAKE) FC=$(FC) $(SRCOBJ) $(MODOBJ)

all:    build
	@echo " "
	@echo "===> linking" $(PROGRAM)
	@echo " "
	$(FC) -o $(EXE) $(MODOPT) $(SRCOBJ) $(MODOBJ) $(LDFLAGS)
	@echo "---------------------------------------------------------------------------------------"
	@echo "Executable: " $(EXE)
	@echo "---------------------------------------------------------------------------------------"
	@size $(EXE)
	@echo "---------------------------------------------------------------------------------------"

clean:
	rm -f *.o $(OBJDIR)/*.o* *.mod $(MODDIR)/*.mod *.stb *.dll *.a *.mod

#############################
# source files dependencies
#############################

include Makefile.dependencies
