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
   real(kind=8), parameter :: ps2au=1.d0/au2ps
   real(kind=8), parameter :: au2kcal=627.5095d0
   real(kind=8), parameter :: kcal2au=1.d0/au2kcal
   real(kind=8), parameter :: kb=3.16683d-6         ! Boltzmann constant [a.u./K]

   real(kind=8) :: V              ! electronic coupling in kcal/mol
   real(kind=8) :: lambda         ! reorganization energy in kcal/mol
   real(kind=8) :: dG             ! reaction free energy (bias) in kcal/mol
   real(kind=8) :: taul           ! longest longitudinal relaxation time in ps
   real(kind=8) :: beta           ! 1/kT in (kcal/mol)^(-1)
   real(kind=8) :: k_marcus       ! Marcus nonadiabatic rate constant in ps^(-1)
   real(kind=8) :: k_rips_jortner ! Rips-Jortner nonadiabatic rate constant in ps^(-1)
   real(kind=8) :: k_zusman       ! Zusman nonadiabatic rate constant in ps^(-1)

   !--------------------------------------------------------------------
   !! declaration of access rights for module subroutines and functions

   public :: k_marcus, k_rips_jortner, k_zusman
   public :: set_marcus_parameters
   public :: calculate_marcus_rate_constant
   public :: calculate_rips_jortner_rate_constant
   public :: calculate_zusman_rate_constant
   public :: exp_diabatic_populations
   public :: fit_rate_constant

CONTAINS

   subroutine set_marcus_parameters(v_, lambda_, dg_, eps0_, eps8_, tau2_, temp_)
      real(kind=8), intent(in) :: v_
      real(kind=8), intent(in) :: lambda_
      real(kind=8), intent(in) :: dg_
      real(kind=8), intent(in) :: eps0_
      real(kind=8), intent(in) :: eps8_
      real(kind=8), intent(in) :: tau2_
      real(kind=8), intent(in) :: temp_
      lambda = lambda_
      dg = dg_
      v = v_
      beta=1.d0/(kb*temp_*au2kcal)
      taul = eps8_*tau2_/eps0_
   end subroutine set_marcus_parameters

   subroutine calculate_marcus_rate_constant()
      real(kind=8) :: prefactor, dgact
      prefactor = (2.d0*pi/au2ps)*(V*V*kcal2au*kcal2au)/sqrt(4.d0*pi*kcal2au*kcal2au*lambda/beta)
      dgact = (lambda+dG)*(lambda+dG)/(4.d0*lambda)
      k_marcus = prefactor*exp(-beta*dgact)
   end subroutine calculate_marcus_rate_constant


   subroutine calculate_rips_jortner_rate_constant()
      real(kind=8) :: prefactor1, prefactor2, dgact
      prefactor1 = (2.d0*pi/au2ps)*(V*V*kcal2au*kcal2au)/sqrt(4.d0*pi*kcal2au*kcal2au*lambda/beta)
      prefactor2 = 1.d0 + 4.d0*pi*(V*V*kcal2au*kcal2au)*taul*ps2au/(lambda*kcal2au)
      dgact = (lambda+dG)*(lambda+dG)/(4.d0*lambda)
      k_rips_jortner = prefactor1*prefactor2*exp(-beta*dgact)
   end subroutine calculate_rips_jortner_rate_constant


   subroutine calculate_zusman_rate_constant()
      real(kind=8) :: prefactor1, prefactor2, dgact
      prefactor1 = (2.d0*pi/au2ps)*(V*V*kcal2au*kcal2au)/sqrt(4.d0*pi*kcal2au*kcal2au*lambda/beta)
      prefactor2 = 1.d0 + 4.d0*pi*(V*V*kcal2au*kcal2au)*taul*ps2au/(lambda*kcal2au)*(1.d0 - dG*dG/(lambda*lambda))
      dgact = (lambda+dG)*(lambda+dG)/(4.d0*lambda)
      k_zusman = prefactor1*prefactor2*exp(-beta*dgact)
   end subroutine calculate_zusman_rate_constant


   function exp_diabatic_populations(t_, k_, a_) result(p)
      real(kind=8), intent(in) :: t_    ! time in ps
      real(kind=8), intent(in) :: k_    ! rate constant in ps^-1
      real(kind=8), intent(in) :: a_    ! fit defect
      real(kind=8), dimension(2) :: p   ! diabatic populations
      real(kind=8) :: ebg, defect
      ebg = exp(beta*dG)
      defect = exp(a_+a_*ebg)
      p(1) = (ebg + defect*exp(-k_*(1.d0+ebg)*t_))/(1.d0+ebg)
      p(2) = 1.d0 - p(1)
   end function exp_diabatic_populations


   subroutine fit_rate_constant(n_,t_,pop_,kfit_,rcorr_,defect_)

      integer, intent(in) :: n_
      real(kind=8), dimension(n_), intent(in) :: t_
      real(kind=8), dimension(n_), intent(in) :: pop_

      real(kind=8), intent(out) :: kfit_
      real(kind=8), intent(out) :: rcorr_
      real(kind=8), intent(out) :: defect_

      integer :: i, n
      real(kind=QUAD) :: popt, poptlog, ebg, bb
      real(kind=QUAD) :: sx, sy, sxy, sxx, syy, xav, yav, ssxx, ssyy, ssxy
      real(kind=QUAD), dimension(n_) :: pop

      ebg = exp(beta*dG)
      bb = ebg + 1.0_QUAD

      !-- adjust the population dynamics data so that pop(0) = 1
      do i=1,n_
         pop(i) = pop_(i)
      enddo
      !pop(1) = 1.0_QUAD

      sx  = 0.0_QUAD
      sy  = 0.0_QUAD
      sxy = 0.0_QUAD
      sxx = 0.0_QUAD
      syy = 0.0_QUAD

      n = 0

      open(22,file="poptlog.dat")

      do i=1,n_

         popt = bb*pop(i) - ebg
         if (popt.gt.0.0_QUAD) then
            poptlog = log(popt)/bb
            n = n + 1
         else
            cycle
         endif

         write(22,'(3g20.10)') t_(i), popt, poptlog

         sx = sx + t_(i)
         sy = sy + poptlog
         sxy = sxy + t_(i)*poptlog
         sxx = sxx + t_(i)*t_(i)
         syy = syy + poptlog*poptlog

      enddo

      xav = sx/n
      yav = sy/n

      close(22)

      ssxx = sxx - n*xav*xav
      ssyy = syy - n*yav*yav
      ssxy = sxy - n*xav*yav

      kfit_ = -ssxy/ssxx
      rcorr_ = ssxy*ssxy/ssxx/ssyy

      defect_ = yav + kfit_*xav

   end subroutine fit_rate_constant

END MODULE marcus
