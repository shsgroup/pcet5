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
   real(kind=8) :: kappa_ad       ! Rips-Jortner adiabaticity parameter, Eq.(3.35)
   real(kind=8) :: k_zusman       ! Zusman nonadiabatic rate constant in ps^(-1)
   real(kind=8) :: k_equil        ! equilibrium constant
   real(kind=8) :: n_equil_r      ! equilibrium reactant population
   real(kind=8) :: n_equil_p      ! equilibrium product population

   !--------------------------------------------------------------------
   !! declaration of access rights for module subroutines and functions

   public :: k_marcus, kappa_ad, k_rips_jortner, k_zusman, k_equil, n_equil_r, n_equil_p
   public :: set_marcus_parameters
   public :: calculate_equilibrium_quantities
   public :: calculate_marcus_rate_constant
   public :: calculate_rips_jortner_rate_constant
   public :: calculate_zusman_rate_constant
   public :: exp_diabatic_populations
   public :: fit_rate_constant_log1
   public :: fit_rate_constant_log
   public :: fit_rate_constant_odr

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
      kappa_ad = 4.d0*pi*(v*v*kcal2au*kcal2au)*taul*ps2au/(lambda*kcal2au)
   end subroutine set_marcus_parameters

   subroutine calculate_equilibrium_quantities()
      k_equil = exp(-beta*dg)
      n_equil_r = 1.d0/(1.d0 + k_equil)
      n_equil_p = k_equil/(1.d0 + k_equil)
   end subroutine calculate_equilibrium_quantities


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
      k_rips_jortner = (prefactor1/prefactor2)*exp(-beta*dgact)
   end subroutine calculate_rips_jortner_rate_constant


   subroutine calculate_zusman_rate_constant()

      real(kind=8) :: prefactor1, prefactor2, factorz, dgact, dG_sq, lambda_sq

      dG_sq = dG*dG
      lambda_sq = lambda*lambda

      if (dG_sq.gt.lambda_sq) then
         factorz = 1.d0 - (2.d0*lambda + dG)*(2.d0*lambda + dG)/lambda_sq
      elseif (dG_sq.lt.lambda_sq) then
         factorz = 1.d0 - dG_sq/lambda_sq
      endif

      if (abs(dG_sq-lambda_sq).gt.1.d-8) then
         prefactor1 = (2.d0*pi/au2ps)*(V*V*kcal2au*kcal2au)/sqrt(4.d0*pi*kcal2au*kcal2au*lambda/beta)
         prefactor2 = 1.d0 + 4.d0*pi*(V*V*kcal2au*kcal2au)*taul*ps2au/(lambda*kcal2au)/factorz
      else
         prefactor1 = (2.d0*pi/au2ps)*(V*V*kcal2au*kcal2au)
         prefactor2 = sqrt(4.d0*pi*kcal2au*kcal2au*lambda/beta) + 2.d0*pi*(V*V*kcal2au*kcal2au)*taul*ps2au
      endif

      dgact = (lambda+dG)*(lambda+dG)/(4.d0*lambda)
      k_zusman = (prefactor1/prefactor2)*exp(-beta*dgact)

   end subroutine calculate_zusman_rate_constant


   function exp_diabatic_populations(t_, k_, p0_, defect_) result(p)
      real(kind=8), intent(in)   :: t_       ! time in ps
      real(kind=8), intent(in)   :: k_       ! rate constant in ps^-1
      real(kind=8), intent(in)   :: p0_      ! initial population of the reactant state
      real(kind=8), intent(in)   :: defect_  ! fit defect
      real(kind=8), dimension(2) :: p        ! diabatic populations
      real(kind=8) :: k_eq, pop1_eq, pop2_eq, a

      k_eq = exp(-beta*dG)
      pop1_eq = 1.d0/(1.d0 + k_eq)
      pop2_eq = 1.d0 - pop1_eq

      a = exp(log(p0_-pop1_eq) - defect_)

      p(1) = pop1_eq + a*exp(-k_*t_/pop2_eq)
      p(2) = 1.d0 - p(1)

   end function exp_diabatic_populations


   subroutine fit_rate_constant_log1(n_,t_,pop_,kfit_,rcorr_,shift0_,defect_)

      integer, intent(in) :: n_
      real(kind=8), dimension(n_), intent(in) :: t_
      real(kind=8), dimension(n_), intent(in) :: pop_

      real(kind=8), intent(out) :: kfit_
      real(kind=8), intent(out) :: rcorr_
      real(kind=8), intent(out) :: shift0_, defect_

      integer :: i, n
      real(kind=8) :: popt, poptlog, ebg, bb
      real(kind=8) :: sx, sy, sxy, sxx, syy, xav, yav, ssxx, ssyy, ssxy
      real(kind=8), dimension(n_) :: pop

      ebg = exp(beta*dG)
      bb = ebg + 1.d0

      do i=1,n_
         pop(i) = pop_(i)
      enddo
      shift0_ = 1.d0 - pop(1)

      sx  = 0.d0
      sy  = 0.d0
      sxy = 0.d0
      sxx = 0.d0
      syy = 0.d0

      n = 0

      open(22,file="poptlog.dat")

      do i=1,n_

         popt = bb*pop(i) - ebg
         if (popt.gt.0.d0) then
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

   end subroutine fit_rate_constant_log1


   subroutine fit_rate_constant_log(n_,t_,pop_,kfit_,rcorr_,shift0_,defect_)

      integer, intent(in) :: n_
      real(kind=8), dimension(n_), intent(in) :: t_
      real(kind=8), dimension(n_), intent(in) :: pop_

      real(kind=8), intent(out) :: kfit_
      real(kind=8), intent(out) :: rcorr_
      real(kind=8), intent(out) :: shift0_, defect_

      integer :: i, n
      real(kind=8) :: popt, poptlog, pop1_0, pop1_eq, pop2_eq, k_eq
      real(kind=8) :: sx, sy, sxy, sxx, syy, xav, yav, ssxx, ssyy, ssxy
      real(kind=8), dimension(n_) :: dpop

      k_eq = exp(-beta*dG)
      pop1_eq = 1.d0/(1.d0 + k_eq)
      pop2_eq = 1.d0 - pop1_eq

      pop1_0 = pop_(1)
      shift0_ = 1.d0 - pop1_0

      do i=1,n_
         dpop(i) = pop_(i) - pop1_eq
      enddo

      sx  = 0.d0
      sy  = 0.d0
      sxy = 0.d0
      sxx = 0.d0
      syy = 0.d0

      n = 0

      open(22,file="poptlog.dat")

      do i=1,n_

         popt = dpop(i)
         if (popt.gt.0.d0) then
            poptlog = log(popt)
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

      kfit_ = -pop2_eq*ssxy/ssxx
      rcorr_ = ssxy*ssxy/ssxx/ssyy

      defect_ = log(pop1_0-pop1_eq) - (yav + kfit_*xav/pop2_eq)

   end subroutine fit_rate_constant_log


   subroutine fit_rate_constant_odr(n_,t_,pop_,kfit_,error_,shift0_)

      use odrpack95

      integer, intent(in) :: n_
      real(kind=8), dimension(n_), intent(in) :: t_
      real(kind=8), dimension(n_), intent(in) :: pop_

      real(kind=8), intent(out) :: kfit_
      real(kind=8), intent(out) :: error_, shift0_

      integer :: i
      real(kind=8) :: pop1_0, k_eq, pop1_eq, pop2_eq
      real(kind=8), dimension(n_,1) :: x
      real(kind=8), dimension(n_,1) :: y
      real(kind=8), dimension(2)    :: par
      real(kind=8), dimension(:,:), pointer :: delta

      allocate (delta(n_,1))

      k_eq = exp(-beta*dG)
      pop1_eq = 1.d0/(1.d0 + k_eq)
      pop2_eq = 1.d0 - pop1_eq

      pop1_0 = pop_(1)
      shift0_ = 1.d0 - pop1_0

      !-- population data
      do i=1,n_
         x(i,1) = t_(i)
         y(i,1) = pop_(i) - pop1_eq
         delta(i,1) = 0.d0
      enddo

      !-- initial values for parameters
      par(1) = pop1_0 - pop1_eq
      par(2) = -k_marcus/pop2_eq

      call odr(FCN=exp_model,&
      &        N=n_,&
      &        M=1,&
      &        NP=2,&
      &        NQ=1,&
      &        BETA=par,&
      &        Y=y,&
      &        X=x,&
      &        DELTA=delta,&
      &        JOB=00020,&
      &        MAXIT=10000,&
      &        IPRINT=0001,&
      &        LUNERR=-1,&
      &        LUNRPT=-1)

      kfit_ = -par(2)*pop2_eq

      error_ = 0.d0
      do i=1,n_
         error_ = error_ + delta(i,1)*delta(i,1)
      enddo
      error_ = sqrt(error_)/n_

      shift0_ = pop1_0 - pop1_eq - par(1)

      deallocate (delta)

   end subroutine fit_rate_constant_odr



   SUBROUTINE EXP_MODEL(N,M,NP,NQ,LDN,LDM,LDNP,PAR,XPLUSD,IFIXB,IFIXX,LDFIX,IDEVAL,F,FJACB,FJACD,ISTOP)

      use real_precision

      !-- INPUT ARGUMENTS
      !   (WHICH MUST NOT BE CHANGED BY THIS ROUTINE)

      INTEGER, intent(in)       :: IDEVAL,LDFIX,LDM,LDN,LDNP,M,N,NP,NQ
      INTEGER, intent(in)       :: IFIXB(NP),IFIXX(LDFIX,M)
      real(kind=R8), intent(in) :: PAR(NP),XPLUSD(LDN,M)

      !-- OUTPUT ARGUMENTS
      INTEGER :: ISTOP
      real(kind=R8) :: F(LDN,NQ), FJACB(LDN,LDNP,NQ), FJACD(LDN,LDM,NQ)

      !-- local variables
      integer :: i

      if (par(2).ge.0.d0) then
         istop = 1
         return
      else
         istop = 0
      endif

      !-- COMPUTE FUNCTION
      IF (MOD(IDEVAL,10).GE.1) THEN
         do i=1,N
            F(i,1) = PAR(1)*exp(PAR(2)*XPLUSD(i,1))
         enddo
      END IF

      !-- COMPUTE DERIVATIVES WITH RESPECT TO BETA
      IF (MOD(IDEVAL/10,10).GE.1) THEN
         !-- compute FJACB(I,K,L), I=1,...,N, K=1,...,NP, & L=1,...,NQ >
         do i=1,N
            FJACB(i,1,1) = exp(PAR(2)*XPLUSD(i,1))
            FJACB(i,2,1) = PAR(1)*XPLUSD(i,1)*exp(PAR(2)*XPLUSD(i,1))
         enddo
      END IF

      !-- COMPUTE DERIVATIVES WITH RESPECT TO DELTA
      IF (MOD(IDEVAL/100,10).GE.1) THEN
         !-- compute FJACD(I,J,L), I=1,...,N, J=1,...,M, & L=1,...,NQ >
         do i=1,N
            FJACD(i,1,1) = PAR(1)*PAR(2)*exp(PAR(2)*XPLUSD(i,1))
         enddo
      END IF
      RETURN
      
   END SUBROUTINE exp_model

END MODULE marcus
