#   13 source files
#   16 dependencies

CONSTANTS : constants.f90                                              # M
FLUX : flux.f90                                                        # M
MONTE_CARLO : monte_carlo.f90                                          # M
RATE_PARS : rate_pars.f90                                              # M
SSPLIB : ssplib.f90                                                    # M

constants.f90 :                                                        # S
flux.f90 : CONSTANTS RATE_PARS                                         # S
monte_carlo.f90 : CONSTANTS RATE_PARS FLUX                             # S
quadpack.f :                                                           # S
rate_flux.f90 : CONSTANTS RATE_PARS FLUX MONTE_CARLO SSPLIB            # S
rate_pars.f90 : CONSTANTS                                              # S
second.f90 :                                                           # S
ssplib.f90 :                                                           # S
