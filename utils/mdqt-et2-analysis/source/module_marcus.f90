!!c*******************************************************************
!!c
!!c   MODULE:      marcus
!!c   ACTION:      contains data and subroutines
!!c                for calculation of Marcus nonadiabatic rate
!!c
!!c*******************************************************************

MODULE marcus

   implicit none
   private

   !-------------
   !! parameters

   integer, parameter :: QUAD=selected_real_kind(33,4931)

   real(kind=8), parameter :: pi=3.141592653d0
   real(kind=8), parameter :: au2ps=2.4189d-5
   real(kind=8), parameter :: au2kcal=627.5095d0
   real(kind=8), parameter :: kcal2au=1.d0/au2kcal
   real(kind=8), parameter :: kb=3.16683d-6         ! Boltzmann constant [a.u./K]

   real(kind=8) :: V             ! electronic coupling in kcal/mol
   real(kind=8) :: lambda        ! reorganization energy in kcal/mol
   real(kind=8) :: dG            ! reaction free energy (bias) in kcal/mol
   real(kind=8) :: beta          ! 1/kT in (kcal/mol)^(-1)
   real(kind=8) :: k_marcus      ! nonadiabatic rate constant in ps^(-1)

   !--------------------------------------------------------------------
   !! declaration of access rights for module subroutines and functions

   public :: k_marcus
   public :: set_marcus_parameters
   public :: calculate_marcus_rate_constant
   public :: marcus_diabatic_populations
   public :: fit_rate_constant

CONTAINS

   subroutine set_marcus_parameters(v_, lambda_, dg_, temp_)
      real(kind=8), intent(in) :: v_
      real(kind=8), intent(in) :: lambda_
      real(kind=8), intent(in) :: dg_
      real(kind=8), intent(in) :: temp_
      lambda = lambda_
      dg = dg_
      v = v_
      beta=1.d0/(kb*temp_*au2kcal)
   end subroutine set_marcus_parameters

   subroutine calculate_marcus_rate_constant()
      real(kind=8) :: prefactor, dgact
      prefactor = (2.d0*pi/au2ps)*(V*V*kcal2au*kcal2au)/sqrt(4.d0*pi*kcal2au*kcal2au*lambda/beta)
      dgact = (lambda+dG)*(lambda+dG)/(4.d0*lambda)
      k_marcus = prefactor*exp(-beta*dgact)
   end subroutine calculate_marcus_rate_constant

   function marcus_diabatic_populations(t_, k_) result(p)
      real(kind=8), intent(in) :: t_    ! time in ps
      real(kind=8), intent(in) :: k_    ! rate constant in ps^-1
      real(kind=8), dimension(2) :: p   ! diabatic populations
      p(1) = (1.d0 + exp(-beta*dG)*exp(-k_*(1.d0+exp(beta*dG))*t_))/(1.d0+exp(-beta*dG))
      p(2) = 1.d0 - p(1)
   end function marcus_diabatic_populations

   subroutine fit_rate_constant(n_,t_,pop_,kfit_,rcorr_)

      integer, intent(in) :: n_
      real(kind=8), dimension(n_), intent(in) :: t_
      real(kind=8), dimension(n_), intent(in) :: pop_

      real(kind=8), intent(out) :: kfit_
      real(kind=8), intent(out) :: rcorr_

      integer :: i
      real(kind=QUAD) :: popt, sx, sy, sxy, sxx, syy, bb
      real(kind=QUAD), dimension(n_) :: poptlog

      bb = exp(beta*dG) + 1.0_QUAD

      sxy = 0.0_QUAD
      sxx = 0.0_QUAD

      do i=1,n_
         popt = (exp(beta*dG) + 1.0_QUAD)*pop_(i) - exp(beta*dG)
         if (popt.gt.0.0_QUAD) then
            poptlog(i) = log(popt)/bb
         else
            write(*,*) "**** negative argument in LOG"
            poptlog(i) = 0.0_QUAD
         endif
         sx = sx + t_(i)
         sy = sy + poptlog(i)
         sxy = sxy + t_(i)*poptlog(i)
         sxx = sxx + t_(i)*t_(i)
         syy = syy + poptlog(i)*poptlog(i)
      enddo

      kfit_ = -sxy/sxx
      rcorr_ = (n_*sxy - sx*sy)/sqrt((n_*sxx - sx*sx)*(n_*syy - sy*sy))

   end subroutine fit_rate_constant

END MODULE marcus
