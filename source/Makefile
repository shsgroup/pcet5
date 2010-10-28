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

VERSION = 5.1

MODDIR  = mod
MODOPT  = -module $(MODDIR)
OBJDIR  = obj
EXEDIR  = ../bin
PROGRAM = pcet_$(VERSION)_$(ARCH)
EXE     = $(EXEDIR)/$(PROGRAM).x
LIB     = lib$(PROGRAM).a

# rules

$(OBJDIR)/%.o : %.f90
	$(FC) -c $(FCFLAGS) -o $(OBJDIR)/$*.o $(MODOPT) $*.f90

$(OBJDIR)/%.o : %.f
	$(FC) -c $(FCFLAGS) -o $(OBJDIR)/$*.o $(MODOPT) $*.f

$(MODDIR)/%.mod : $(OBJDIR)/module_%.o
	@if [ ! -f $@ ]; then rm $< ; $(MAKE) $< ; fi


SRCOBJ = $(OBJDIR)/pcet.o      \
         $(OBJDIR)/bubbli.o    \
         $(OBJDIR)/banner.o    \
         $(OBJDIR)/deinitmat.o \
         $(OBJDIR)/et2.o       \
         $(OBJDIR)/getmin.o    \
         $(OBJDIR)/initmat.o   \
         $(OBJDIR)/minima2.o   \
         $(OBJDIR)/minima3.o   \
         $(OBJDIR)/path2.o     \
         $(OBJDIR)/path3.o     \
         $(OBJDIR)/primat.o    \
         $(OBJDIR)/prngas.o    \
         $(OBJDIR)/prnsol.o    \
         $(OBJDIR)/quadpack.o  \
         $(OBJDIR)/rate2.o     \
         $(OBJDIR)/rate3.o     \
         $(OBJDIR)/rateb.o     \
         $(OBJDIR)/second.o    \
         $(OBJDIR)/setjob.o    \
         $(OBJDIR)/setmat.o    \
         $(OBJDIR)/surface2.o  \
         $(OBJDIR)/surface3.o  \
         $(OBJDIR)/ugas.o      \
         $(OBJDIR)/usol.o      \
         $(OBJDIR)/wavef2.o    \
         $(OBJDIR)/wavef3.o    \
         $(OBJDIR)/weight2.o   \
         $(OBJDIR)/weight3.o   \
         $(OBJDIR)/wfn2prt.o   \
         $(OBJDIR)/wfn3prt.o   \
         $(OBJDIR)/zcoup.o

SRCSRC = pcet.f90      \
         banner.f90    \
         bubbli.f90    \
         deinitmat.f90 \
         et2.f90       \
         getmin.f90    \
         initmat.f90   \
         minima2.f90   \
         minima3.f90   \
         path2.f90     \
         path3.f90     \
         primat.f      \
         prngas.f90    \
         prnsol.f90    \
         quadpack.f    \
         rate2.f90     \
         rate3.f90     \
         rateb.f90     \
         second.f90    \
         setjob.f90    \
         setmat.f90    \
         surface2.f90  \
         surface3.f90  \
         ugas.f90      \
         usol.f90      \
         wavef2.f90    \
         wavef3.f90    \
         weight2.f90   \
         weight3.f90   \
         wfn2prt.f90   \
         wfn3prt.f90   \
         zcoup.f90


MODOBJ  = $(OBJDIR)/module_control.o     \
          $(OBJDIR)/module_cst.o         \
          $(OBJDIR)/module_eispack.o     \
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
          $(OBJDIR)/module_monte_carlo.o \
          $(OBJDIR)/module_msevb_water.o \
          $(OBJDIR)/module_pardim.o      \
          $(OBJDIR)/module_parsol.o      \
          $(OBJDIR)/module_potential.o   \
          $(OBJDIR)/module_quantum.o     \
          $(OBJDIR)/module_rate_flux.o   \
          $(OBJDIR)/module_solmat.o      \
          $(OBJDIR)/module_ssplib.o      \
          $(OBJDIR)/module_strings.o


MODMOD  = $(MODDIR)/control.mod     \
          $(MODDIR)/cst.mod         \
          $(MODDIR)/eispack.mod     \
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
          $(MODDIR)/monte_carlo.mod \
          $(MODDIR)/msevb_water.mod \
          $(MODDIR)/pardim.mod      \
          $(MODDIR)/parsol.mod      \
          $(MODDIR)/potential.mod   \
          $(MODDIR)/quantum.mod     \
          $(MODDIR)/rate_flux.mod   \
          $(MODDIR)/solmat.mod      \
          $(MODDIR)/ssplib.mod      \
          $(MODDIR)/strings.mod

MODSRC  = module_control.f90     \
          module_cst.f90         \
          module_eispack.f90     \
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
          module_monte_carlo.f90 \
          module_msevb_water.f90 \
          module_pardim.f90      \
          module_parsol.f90      \
          module_potential.f90   \
          module_quantum.f90     \
          module_rate_flux.f90   \
          module_solmat.f90      \
          module_ssplib.f90      \
          module_strings.f90

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
	@echo "*** building object files ***"
	$(MAKE) FC=$(FC) $(SRCOBJ) $(MODOBJ)

all:    build
	@echo "*** linking $(PROGRAM) ***"
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $(EXE) $(MODOPT) $(SRCOBJ) $(MODOBJ)
	@echo "---------------------------------------------------------------------------------------"
	@size $(EXE)
	@echo "---------------------------------------------------------------------------------------"

clean:
	rm -f *.o $(OBJDIR)/*.o *.mod $(MODDIR)/*.mod *.stb *.dll *.a

#############################
# source files dependencies
#############################

$(OBJDIR)/bubbli.o : bubbli.f90 

$(OBJDIR)/banner.o : banner.f90 

$(OBJDIR)/deinitmat.o : deinitmat.f90 \
$(OBJDIR)/module_quantum.o \
$(OBJDIR)/module_geosol.o \
$(OBJDIR)/module_geogas.o \
$(OBJDIR)/module_solmat.o \
$(OBJDIR)/module_gasmat.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/et2.o : et2.f90 \
$(OBJDIR)/module_msevb_water.o \
$(OBJDIR)/module_potential.o \
$(OBJDIR)/module_elcm.o \
$(OBJDIR)/module_frcm.o \
$(OBJDIR)/module_geosol.o \
$(OBJDIR)/module_geogas.o \
$(OBJDIR)/module_control.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_strings.o \
$(OBJDIR)/module_keys.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/getmin.o : getmin.f90 

$(OBJDIR)/initmat.o : initmat.f90 \
$(OBJDIR)/module_quantum.o \
$(OBJDIR)/module_solmat.o \
$(OBJDIR)/module_gasmat.o \
$(OBJDIR)/module_control.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/minima2.o : minima2.f90 \
$(OBJDIR)/module_fesmin_2d.o \
$(OBJDIR)/module_control.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_keys.o \
$(OBJDIR)/module_strings.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/minima3.o : minima3.f90 \
$(OBJDIR)/module_fesmin_3d.o \
$(OBJDIR)/module_quantum.o \
$(OBJDIR)/module_control.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_keys.o \
$(OBJDIR)/module_strings.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/path2.o : path2.f90 \
$(OBJDIR)/module_feszz_2d.o \
$(OBJDIR)/module_quantum.o \
$(OBJDIR)/module_control.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_keys.o \
$(OBJDIR)/module_strings.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/path3.o : path3.f90 \
$(OBJDIR)/module_feszz_3d.o \
$(OBJDIR)/module_quantum.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_keys.o \
$(OBJDIR)/module_strings.o \
$(OBJDIR)/module_control.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/pcet.o : pcet.f90 \
$(OBJDIR)/module_pardim.o \
$(OBJDIR)/module_geometry.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_strings.o \
$(OBJDIR)/module_keys.o 

$(OBJDIR)/primat.o : primat.f 

$(OBJDIR)/prngas.o : prngas.f90 \
$(OBJDIR)/module_quantum.o \
$(OBJDIR)/module_strings.o \
$(OBJDIR)/module_keys.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/prnsol.o : prnsol.f90 \
$(OBJDIR)/module_quantum.o \
$(OBJDIR)/module_strings.o \
$(OBJDIR)/module_keys.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/rate2.o : rate2.f90 \
$(OBJDIR)/module_feszz_2d.o \
$(OBJDIR)/module_fesmin_2d.o \
$(OBJDIR)/module_solmat.o \
$(OBJDIR)/module_gasmat.o \
$(OBJDIR)/module_quantum.o \
$(OBJDIR)/module_control.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_strings.o \
$(OBJDIR)/module_keys.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/rate3.o : rate3.f90 \
$(OBJDIR)/module_feszz_3d.o \
$(OBJDIR)/module_fesmin_3d.o \
$(OBJDIR)/module_geogas.o \
$(OBJDIR)/module_solmat.o \
$(OBJDIR)/module_gasmat.o \
$(OBJDIR)/module_quantum.o \
$(OBJDIR)/module_control.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_strings.o \
$(OBJDIR)/module_keys.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/rateb.o : rateb.f90 \
$(OBJDIR)/module_feszz_3d.o \
$(OBJDIR)/module_fesmin_3d.o \
$(OBJDIR)/module_geogas.o \
$(OBJDIR)/module_solmat.o \
$(OBJDIR)/module_gasmat.o \
$(OBJDIR)/module_quantum.o \
$(OBJDIR)/module_control.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_strings.o \
$(OBJDIR)/module_keys.o \
$(OBJDIR)/module_rate_flux.o \
$(OBJDIR)/module_pardim.o 


$(OBJDIR)/second.o : second.f90 

$(OBJDIR)/setjob.o : setjob.f90 \
$(OBJDIR)/module_msevb_water.o \
$(OBJDIR)/module_frcm.o \
$(OBJDIR)/module_potential.o \
$(OBJDIR)/module_geometry.o \
$(OBJDIR)/module_geosol.o \
$(OBJDIR)/module_geogas.o \
$(OBJDIR)/module_parsol.o \
$(OBJDIR)/module_quantum.o \
$(OBJDIR)/module_control.o \
$(OBJDIR)/module_strings.o \
$(OBJDIR)/module_keys.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/setmat.o : setmat.f90 \
$(OBJDIR)/module_potential.o \
$(OBJDIR)/module_msevb_water.o \
$(OBJDIR)/module_elcm.o \
$(OBJDIR)/module_frcm.o \
$(OBJDIR)/module_quantum.o \
$(OBJDIR)/module_solmat.o \
$(OBJDIR)/module_gasmat.o \
$(OBJDIR)/module_geosol.o \
$(OBJDIR)/module_geogas.o \
$(OBJDIR)/module_control.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_strings.o \
$(OBJDIR)/module_keys.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/surface2.o : surface2.f90 \
$(OBJDIR)/module_feszz_2d.o \
$(OBJDIR)/module_quantum.o \
$(OBJDIR)/module_control.o \
$(OBJDIR)/module_solmat.o \
$(OBJDIR)/module_strings.o \
$(OBJDIR)/module_keys.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/surface3.o : surface3.f90 \
$(OBJDIR)/module_feszz_3d.o \
$(OBJDIR)/module_quantum.o \
$(OBJDIR)/module_control.o \
$(OBJDIR)/module_solmat.o \
$(OBJDIR)/module_strings.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_keys.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/ugas.o : ugas.f90 \
$(OBJDIR)/module_eispack.o \
$(OBJDIR)/module_solmat.o \
$(OBJDIR)/module_gasmat.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/usol.o : usol.f90 \
$(OBJDIR)/module_eispack.o \
$(OBJDIR)/module_solmat.o \
$(OBJDIR)/module_gasmat.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/wavef2.o : wavef2.f90 \
$(OBJDIR)/module_strings.o \
$(OBJDIR)/module_keys.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/wavef3.o : wavef3.f90 \
$(OBJDIR)/module_quantum.o \
$(OBJDIR)/module_strings.o \
$(OBJDIR)/module_keys.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/weight2.o : weight2.f90 \
$(OBJDIR)/module_feszz_2d.o \
$(OBJDIR)/module_quantum.o \
$(OBJDIR)/module_control.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/weight3.o : weight3.f90 \
$(OBJDIR)/module_feszz_3d.o \
$(OBJDIR)/module_quantum.o \
$(OBJDIR)/module_control.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/wfn2prt.o : wfn2prt.f90 \
$(OBJDIR)/module_feszz_2d.o \
$(OBJDIR)/module_feszz_3d.o \
$(OBJDIR)/module_quantum.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_control.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/wfn3prt.o : wfn3prt.f90 \
$(OBJDIR)/module_feszz_3d.o \
$(OBJDIR)/module_quantum.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/zcoup.o : zcoup.f90 


#############################
# module files dependencies
#############################

$(OBJDIR)/module_control.o : module_control.f90 

$(OBJDIR)/module_cst.o : module_cst.f90 

$(OBJDIR)/module_ssplib.o : module_ssplib.f90 

$(OBJDIR)/module_monte_carlo.o : module_monte_carlo.f90 \
$(OBJDIR)/module_cst.o

$(OBJDIR)/module_eispack.o : module_eispack.f90

$(OBJDIR)/module_elcm.o : module_elcm.f90 \
$(OBJDIR)/module_geosol.o \
$(OBJDIR)/module_parsol.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/module_elmnts.o : module_elmnts.f90 

$(OBJDIR)/module_fesmin_2d.o : module_fesmin_2d.f90 \
$(OBJDIR)/module_feszz_2d.o \
$(OBJDIR)/module_quantum.o \
$(OBJDIR)/module_control.o \
$(OBJDIR)/module_eispack.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/module_fesmin_3d.o : module_fesmin_3d.f90 \
$(OBJDIR)/module_feszz_3d.o \
$(OBJDIR)/module_quantum.o \
$(OBJDIR)/module_control.o \
$(OBJDIR)/module_eispack.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/module_feszz_2d.o : module_feszz_2d.f90 \
$(OBJDIR)/module_eispack.o \
$(OBJDIR)/module_quantum.o \
$(OBJDIR)/module_solmat.o \
$(OBJDIR)/module_gasmat.o \
$(OBJDIR)/module_control.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/module_feszz_3d.o : module_feszz_3d.f90 \
$(OBJDIR)/module_eispack.o \
$(OBJDIR)/module_quantum.o \
$(OBJDIR)/module_solmat.o \
$(OBJDIR)/module_gasmat.o \
$(OBJDIR)/module_control.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/module_frcm.o : module_frcm.f \
$(OBJDIR)/module_strings.o \
$(OBJDIR)/module_elmnts.o \
$(OBJDIR)/module_parsol.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/module_gasmat.o : module_gasmat.f90 

$(OBJDIR)/module_geogas.o : module_geogas.f90 \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/module_geometry.o : module_geometry.f90 \
$(OBJDIR)/banner.o \
$(OBJDIR)/module_control.o \
$(OBJDIR)/module_elmnts.o \
$(OBJDIR)/module_keys.o \
$(OBJDIR)/module_strings.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/module_geosol.o : module_geosol.f90 \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/module_keys.o : module_keys.f90 

$(OBJDIR)/module_msevb_water.o : module_msevb_water.f90 \
$(OBJDIR)/module_geogas.o \
$(OBJDIR)/module_elmnts.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/module_pardim.o : module_pardim.f90 

$(OBJDIR)/module_parsol.o : module_parsol.f90 

$(OBJDIR)/module_potential.o : module_potential.f90 \
$(OBJDIR)/module_geogas.o \
$(OBJDIR)/module_quantum.o \
$(OBJDIR)/module_strings.o \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_pardim.o 

$(OBJDIR)/module_quantum.o : module_quantum.f90 \
$(OBJDIR)/module_eispack.o \
$(OBJDIR)/module_control.o \
$(OBJDIR)/module_cst.o 

$(OBJDIR)/module_rate_flux.o : module_rate_flux.f90 \
$(OBJDIR)/module_cst.o \
$(OBJDIR)/module_monte_carlo.o \
$(OBJDIR)/module_ssplib.o 


$(OBJDIR)/module_solmat.o : module_solmat.f90 

$(OBJDIR)/module_strings.o : module_strings.f90 

