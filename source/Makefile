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

VERSION = 5.3

###MODDIR  = ./mod

ifeq ($(FC),gfortran)
   MODDIR = .
else
   MODDIR  = ./mod
endif

ifeq ($(FC),openf90)
   MODDIR = .
endif

MODOPT  = $(MODFLAG)$(MODDIR)
OBJDIR  = ./obj
EXEDIR  = ../bin
PROGRAM = pcet_$(VERSION)_$(ARCH)
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


SRCOBJ = $(OBJDIR)/pcet.o        \
         $(OBJDIR)/deinitmat.o   \
         $(OBJDIR)/dynamics3.o   \
         $(OBJDIR)/dynamicset2.o \
         $(OBJDIR)/et2.o         \
         $(OBJDIR)/getmin.o      \
         $(OBJDIR)/initmat.o     \
         $(OBJDIR)/lbfgsb.o      \
         $(OBJDIR)/rdiag.o      \
         $(OBJDIR)/minima2.o     \
         $(OBJDIR)/minima3.o     \
         $(OBJDIR)/path2.o       \
         $(OBJDIR)/path3.o       \
         $(OBJDIR)/primat.o      \
         $(OBJDIR)/prngas.o      \
         $(OBJDIR)/prnsol.o      \
         $(OBJDIR)/quadpack.o    \
         $(OBJDIR)/rate2.o       \
         $(OBJDIR)/rate3.o       \
         $(OBJDIR)/rateb.o       \
         $(OBJDIR)/setjob.o      \
         $(OBJDIR)/setmat.o      \
         $(OBJDIR)/surface2.o    \
         $(OBJDIR)/surface3.o    \
         $(OBJDIR)/ugas.o        \
         $(OBJDIR)/usol.o        \
         $(OBJDIR)/wavef2.o      \
         $(OBJDIR)/wavef3.o      \
         $(OBJDIR)/weight2.o     \
         $(OBJDIR)/weight3.o     \
         $(OBJDIR)/wfn2prt.o     \
         $(OBJDIR)/wfn3prt.o     \
         $(OBJDIR)/zcoup.o

SRCSRC = pcet.f90        \
         banner.f90      \
         deinitmat.f90   \
         dynamics3.f90   \
         dynamicset2.f90 \
         et2.f90         \
         getmin.f90      \
         initmat.f90     \
         lbfgsb.f        \
         rdiag.f         \
         minima2.f90     \
         minima3.f90     \
         path2.f90       \
         path3.f90       \
         primat.f        \
         prngas.f90      \
         prnsol.f90      \
         quadpack.f      \
         rate2.f90       \
         rate3.f90       \
         rateb.f90       \
         setjob.f90      \
         setmat.f90      \
         surface2.f90    \
         surface3.f90    \
         ugas.f90        \
         usol.f90        \
         wavef2.f90      \
         wavef3.f90      \
         weight2.f90     \
         weight3.f90     \
         wfn2prt.f90     \
         wfn3prt.f90     \
         zcoup.f90


MODOBJ  = $(OBJDIR)/module_control.o     \
          $(OBJDIR)/module_control_dynamics.o \
          $(OBJDIR)/module_data_et2.o \
          $(OBJDIR)/module_cst.o         \
          $(OBJDIR)/module_dielectric_functions.o \
          $(OBJDIR)/module_eispack.o     \
          $(OBJDIR)/module_lapack_wrappers.o     \
          $(OBJDIR)/module_elcm.o        \
          $(OBJDIR)/module_elmnts.o      \
          $(OBJDIR)/module_fesmin_2d.o   \
          $(OBJDIR)/module_fesmin_3d.o   \
          $(OBJDIR)/module_feszz_2d.o    \
          $(OBJDIR)/module_feszz_3d.o    \
          $(OBJDIR)/module_frcm.o        \
          $(OBJDIR)/module_gasmat.o      \
          $(OBJDIR)/module_geogas.o      \
          $(OBJDIR)/module_geometry.o    \
          $(OBJDIR)/module_geosol.o      \
          $(OBJDIR)/module_keys.o        \
          $(OBJDIR)/module_laser.o       \
          $(OBJDIR)/module_monte_carlo.o \
          $(OBJDIR)/module_msevb_water.o \
          $(OBJDIR)/module_pardim.o      \
          $(OBJDIR)/module_parsol.o      \
          $(OBJDIR)/module_potential.o   \
          $(OBJDIR)/module_propagators_3d.o \
          $(OBJDIR)/module_propagators_et2.o \
          $(OBJDIR)/module_rk_parameters.o \
          $(OBJDIR)/module_quantum.o     \
          $(OBJDIR)/module_random_generators.o \
          $(OBJDIR)/module_rate_flux.o   \
          $(OBJDIR)/module_solmat.o      \
          $(OBJDIR)/module_sorting.o      \
          $(OBJDIR)/module_ssplib.o      \
          $(OBJDIR)/module_strings.o     \
          $(OBJDIR)/module_timers.o      \
          $(OBJDIR)/module_turbomole_wrappers.o


MODMOD  = $(MODDIR)/control.mod     \
          $(MODDIR)/control_dynamics.mod \
          $(MODDIR)/cst.mod         \
          $(MODDIR)/dielectric_functions.mod \
          $(MODDIR)/eispack.mod     \
          $(MODDIR)/data_et2.mod     \
          $(MODDIR)/lapack_wrappers.mod     \
          $(MODDIR)/elcm.mod        \
          $(MODDIR)/elmnts.mod      \
          $(MODDIR)/fesmin_2d.mod   \
          $(MODDIR)/fesmin_3d.mod   \
          $(MODDIR)/feszz_2d.mod    \
          $(MODDIR)/feszz_3d.mod    \
          $(MODDIR)/frcm.mod        \
          $(MODDIR)/gasmat.mod      \
          $(MODDIR)/geogas.mod      \
          $(MODDIR)/geometry.mod    \
          $(MODDIR)/geosol.mod      \
          $(MODDIR)/keys.mod        \
          $(MODDIR)/laser.mod       \
          $(MODDIR)/monte_carlo.mod \
          $(MODDIR)/msevb_water.mod \
          $(MODDIR)/pardim.mod      \
          $(MODDIR)/parsol.mod      \
          $(MODDIR)/potential.mod   \
          $(MODDIR)/propagators_3d.mod \
          $(MODDIR)/propagators_et2.mod \
          $(MODDIR)/rk_parameters.mod \
          $(MODDIR)/quantum.mod     \
          $(MODDIR)/random_generators.mod \
          $(MODDIR)/rate_flux.mod   \
          $(MODDIR)/solmat.mod      \
          $(MODDIR)/ssplib.mod      \
          $(MODDIR)/strings.mod     \
          $(MODDIR)/timers.mod      \
          $(MODDIR)/turbomole_wrappers.mod

MODSRC  = module_control.f90     \
          module_control_dynamics.f90 \
          module_cst.f90         \
          module_dielectric_functions.f90 \
          module_eispack.f90     \
          module_data_et2.f90     \
          module_lapack_wrappers.f90     \
          module_elcm.f90        \
          module_elmnts.f90      \
          module_fesmin_2d.f90   \
          module_fesmin_3d.f90   \
          module_feszz_2d.f90    \
          module_feszz_3d.f90    \
          module_frcm.f          \
          module_gasmat.f90      \
          module_geogas.f90      \
          module_geometry.f90    \
          module_geosol.f90      \
          module_keys.f90        \
          module_laser.f90       \
          module_monte_carlo.f90 \
          module_msevb_water.f90 \
          module_pardim.f90      \
          module_parsol.f90      \
          module_potential.f90   \
          module_propagators_3d.f90 \
          module_propagators_et2.f90 \
          module_rk_parameters.f90 \
          module_quantum.f90     \
          module_random_generators.f90 \
          module_rate_flux.f90   \
          module_solmat.f90      \
          module_sorting.f90      \
          module_ssplib.f90      \
          module_strings.f90     \
          module_timers.f90      \
          module_turbomole_wrappers.f90

note:
	@echo
	@echo '=== MAKEFILE for PCET 5.x ==='
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
	$(MAKE) FC=$(FC) $(SRCOBJ) $(OBJDIR)/banner.o $(MODOBJ)

all:    build
	@echo " "
	@echo "===> linking" $(PROGRAM)
	@echo " "
	$(FC) -o $(EXE) $(MODOPT) $(SRCOBJ) $(OBJDIR)/banner.o $(MODOBJ) $(LDFLAGS)
	@echo "---------------------------------------------------------------------------------------"
	@echo "Executable: " $(EXE)
	@echo "---------------------------------------------------------------------------------------"
	@size $(EXE)
	@echo "---------------------------------------------------------------------------------------"

clean:
	rm -f *.o $(OBJDIR)/*.o* *.mod $(MODDIR)/*.mod *.stb *.dll *.a

#############################
# source files dependencies
#############################

include Makefile.dep
