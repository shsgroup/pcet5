module potential

!----------------------------------------------------------------
!  Gas-phase potentials
!----------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:36 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
!
!----------------------------------------------------------------

   use pardim
   use cst
   use strings
   use quantum
   use geogas

   implicit none
   private

   !-----------------------------------------------------------------------------
   ! constant potential for a fixed solute (ET model)
   ! (parameters for ...
   !-----------------------------------------------------------------------------
   type, private :: constant
      real(8) :: v12
      real(8) :: v34
      real(8) :: v13
      real(8) :: v24
      real(8) :: v14
      real(8) :: v23
      real(8), dimension(4) :: bias  !-energy bias parameters for diabatic states
   end type constant

   !-----------------------------------------------------------------------------
   ! one-dimensional harmonic potential for general D--H--A systems (DA is fixed)
   !-----------------------------------------------------------------------------
   type, private :: harm_1d
      real(8) :: r0dh
      real(8) :: r0ah
      real(8) :: omega_dh
      real(8) :: omega_ah
      real(8) :: vpt_1
      real(8) :: vpt_2
      real(8) :: vet_a
      real(8) :: vet_b
      real(8) :: vept_1
      real(8) :: vept_2
      real(8), dimension(4) :: corr  ! energy bias parameters for diabatic states
   end type harm_1d

   !-----------------------------------------------------------------------------
   ! one-dimensional PCET potential for Nocera-like N-H-O systems
   !-----------------------------------------------------------------------------
   type, private :: mm5_1d
      real(8) :: dnh
      real(8) :: doh
      real(8) :: bnh
      real(8) :: boh
      real(8) :: rnh
      real(8) :: roh
      real(8) :: drnh
      real(8) :: droh
      real(8) :: brnh
      real(8) :: broh
      real(8) :: croh
      real(8) :: crnh
      real(8) :: ksi
      real(8) :: vpt0_1
      real(8) :: gapt_1
      real(8) :: qpt0_1
      real(8) :: vpt0_2
      real(8) :: gapt_2
      real(8) :: qpt0_2
      real(8) :: vet0_a
      real(8) :: gaet_a
      real(8) :: qet0_a
      real(8) :: vet0_b
      real(8) :: gaet_b
      real(8) :: qet0_b
      real(8) :: VETPT0_1
      real(8) :: GAETPT_1
      real(8) :: QETPT0_1
      real(8) :: VETPT0_2
      real(8) :: GAETPT_2
      real(8) :: QETPT0_2
      real(8), dimension(4) :: corr  ! energy bias parameters for diabatic states
   end type mm5_1d

   !-----------------------------------------------------------------------------
   ! one-dimensional PCET potential for O-H-O systems with a significant
   ! difference of the O-H frequencies for reduced and oxydized states
   ! (designed for Tyr-O-H...O-A systems)
   !-----------------------------------------------------------------------------
   type, private :: mm5gen_1d

      !-- dissociation energies
      real(8) :: doh1a
      real(8) :: doh1b
      real(8) :: doh2a
      real(8) :: doh2b

      !-- beta parameters
      real(8) :: boh1a
      real(8) :: boh1b
      real(8) :: boh2a
      real(8) :: boh2b

      !-- equilibrium distances
      real(8) :: roh1a
      real(8) :: roh1b
      real(8) :: roh2a
      real(8) :: roh2b

      !-- repulsion parameters
      real(8) :: droh1a
      real(8) :: droh1b
      real(8) :: droh2a
      real(8) :: droh2b
      real(8) :: broh1a
      real(8) :: broh1b
      real(8) :: broh2a
      real(8) :: broh2b
      real(8) :: croh1a
      real(8) :: croh1b
      real(8) :: croh2a
      real(8) :: croh2b

      !-- Coulomb interaction smoothing parameter
      real(8) :: ksi

      !-- off-diagonal coupling parameters
      real(8) :: vpt0_1
      real(8) :: gapt_1
      real(8) :: qpt0_1
      real(8) :: vpt0_2
      real(8) :: gapt_2
      real(8) :: qpt0_2
      real(8) :: vet0_a
      real(8) :: gaet_a
      real(8) :: qet0_a
      real(8) :: vet0_b
      real(8) :: gaet_b
      real(8) :: qet0_b
      real(8) :: vetpt0_1
      real(8) :: gaetpt_1
      real(8) :: qetpt0_1
      real(8) :: vetpt0_2
      real(8) :: gaetpt_2
      real(8) :: qetpt0_2

      !-- energy bias parameters for diabatic states
      real(8), dimension(4) :: corr

   end type mm5gen_1d

   !-----------------------------------------------------------------------------
   ! two-dimensional LEPS-like potential for Nocera-like N-H-O PCET systems
   !-----------------------------------------------------------------------------
   type, private :: leps_2D
      real(8) :: dah     ! AH dissociation energy
      real(8) :: dbh     ! BH dissociation energy
      real(8) :: dab     ! AB dissociation energy
      real(8) :: bah     ! AH beta (Morse force constant)
      real(8) :: bbh     ! BH beta (Morse force constant)
      real(8) :: bab     ! AB beta (Morse force constant)
      real(8) :: rah     ! AH equilibrium distance
      real(8) :: rbh     ! BH equilibrium distance
      real(8) :: rab     ! AB equilibrium distance
      real(8) :: r0Dd1
      real(8) :: r0Dd2
      real(8) :: r0Aa1
      real(8) :: r0Aa2
      real(8) :: kah     ! AH Sato parameter
      real(8) :: kbh     ! BH Sato parameter
      real(8) :: kab     ! AB Sato parameter
      real(8) :: sigma   ! PT overlap integral
      real(8) :: VPT0
      real(8) :: ALPH
      real(8) :: VET0
      real(8) :: BETA
      real(8) :: kd1
      real(8) :: kd2
      real(8) :: ka1
      real(8) :: ka2
      real(8), dimension(4) :: corr  ! energy bias parameters for diabatic states
   end type leps_2D

   !--------------------------------------------------------------------------------
   ! two-dimensional hybrid (harmonic-LEPS) potential for DE--A-H-B--AE PCET systems
   !--------------------------------------------------------------------------------
   type, public :: hybrid_2D
      real(8) :: dah                 ! AH dissociation energy
      real(8) :: dbh                 ! BH dissociation energy
      real(8) :: omega_ah            ! A-H frequency
      real(8) :: omega_bh            ! B-H frequency
      real(8) :: omega_ab            ! A-B frequency
      real(8) :: rah                 ! A-H equilibrium distance
      real(8) :: rbh                 ! B-H equilibrium distance
      real(8) :: rab1                ! A-B equilibrium distance for ET reactants
      real(8) :: rab2                ! A-B equilibrium distance for ET products
      real(8) :: kah                 ! AH Sato parameter
      real(8) :: kbh                 ! BH Sato parameter
      real(8) :: sigma               ! PT overlap integral
      real(8) :: vpt0                ! PT off-diagonal coupling constant
      real(8) :: alph                ! PT off-diagonal coupling exponent
      real(8) :: vet0                ! ET off-diagonal coupling constant
      real(8) :: beta                ! ET off-diagonal coupling exponent
      real(8) :: vept0               ! EPT off-diagonal coupling constant
      real(8), dimension(4) :: corr  ! energy bias parameters for diabatic states
   end type hybrid_2D

   !-----------------------------------------------------------------------------
   type(constant), private, target, save :: constpar = constant(&
                                     & zero,zero,zero,zero,zero,zero,&
                                     & (/zero,zero,zero,zero/))

   type(mm5_1D),    private, target, save :: mm5par  = mm5_1D(&
                                     & zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,&
                                     & zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,&
                                     & zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,&
                                     & (/zero,zero,zero,zero/))

   type(mm5gen_1D), private, target, save :: mm5genpar = mm5gen_1D(&
                                     & zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,&
                                     & zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,&
                                     & zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,&
                                     & zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,&
                                     & zero,zero,zero,&
                                     & (/zero,zero,zero,zero/))

   type(leps_2D),   private, target, save :: lepspar = leps_2D(&
                                     & zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,&
                                     & zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,&
				     & zero,zero,zero,zero,zero,&
                                     & (/zero,zero,zero,zero/))

   type(hybrid_2D), private, target, save :: hybpar = hybrid_2D(&
                                     & zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,&
                                     & zero,zero,zero,zero,zero,zero,zero,&
                                     & (/zero,zero,zero,zero/))

   type(harm_1d), private, target, save :: harmpar = harm_1d(&
                                     & zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,&
                                     & (/zero,zero,zero,zero/))


   logical :: leadsp
   logical, public, save  :: coulomb=.true.
   integer, dimension(60) :: istart
   real(8), dimension(10) :: vpari

   character(120) :: line
   character( 80) :: string
   character(  8) :: par_symbol

   !----------------------------------------------------------------
   public  :: set_potential
   public  :: h0mat_mm5, dh0mat_mm5, d2h0mat_mm5
   public  :: h0mat_mm5gen, dh0mat_mm5gen, d2h0mat_mm5gen
   public  :: h0mat_leps0, h0mat_leps2
   public  :: h0mat_hyb
   public  :: h0mat_harm, dh0mat_harm, d2h0mat_harm
   public  :: h0mat_constant
   private :: sethyb, setleps, setmm5, setmm5gen, setharm, setconst

!----------------------------------------------------------------
contains
!----------------------------------------------------------------

   subroutine set_potential(iunit,id)

      implicit none
      integer, intent(in) :: iunit  ! parameter file unit
      integer, intent(in) :: id     ! id of the potential (unique)

      type(constant)  :: constpar_
      type(harm_1d)   :: harmpar_
      type(mm5_1d)    :: mm5par_
      type(mm5gen_1d) :: mm5genpar_
      type(leps_2d)   :: lepspar_
      type(hybrid_2d) :: hybpar_

      logical :: l_exist, l_opened

      inquire(unit=iunit,exist=l_exist,opened=l_opened)

      if (.not.l_exist) then
         write(*,*) "Error in set_potential: parameter file does not exist..."
         stop
      endif

      if (.not.l_opened) then
         write(*,*) "Error in set_potential: parameter file is not opened..."
         stop
      endif

      select case(id)

         case(0)
            call setconst(iunit,constpar_)
            constpar = constpar_

         case(1)
            call setmm5(iunit,mm5par_)
            mm5par = mm5par_

         case(3:4)
            call setleps(iunit,lepspar_)
            lepspar = lepspar_

         case(5)
            call sethyb(iunit,hybpar_)
            hybpar = hybpar_

         case(6)
            call setmm5gen(iunit,mm5genpar_)
            mm5genpar = mm5genpar_

         case(7)
            call setharm(iunit,harmpar_)
            harmpar = harmpar_

      end select

   end subroutine set_potential

   !======================================================================
   ! Reads parameters of the constant potential
   !======================================================================
   subroutine setconst(ifile,pars)

      implicit none

      integer, intent(in) :: ifile
      type(constant), intent(out) :: pars

      integer, parameter :: npars=7
      character(8),  dimension(npars) :: pa
      character(15), dimension(npars) :: dime

      integer :: io_status, i, nvalue, numpar
      real(8), dimension(10) :: vpari
      real(8)  :: vpar

      pa = (/ 'CORR    ',&
              'V12     ',&
              'V34     ',&
              'V13     ',&
              'V24     ',&
              'V14     ',&
              'V23     '   /)

      dime = (/ 'kcal/mol       ',&
                'kcal/mol       ',&
                'kcal/mol       ',&
                'kcal/mol       ',&
                'kcal/mol       ',&
                'kcal/mol       ',&
                'kcal/mol       '   /)

      !-- PT coupling parameters VPT
      pars%v12 = 50.d0
      pars%v34 = 50.d0

      !-- ET coupling parameters VET
      pars%v13 = 1.d0
      pars%v24 = 1.d0

      !-- EPT mixed coupling parameters VEPT
      pars%v14 = zero
      pars%v23 = zero

      ! Constant correction terms
      pars%bias(1) = zero
      pars%bias(2) = zero
      pars%bias(3) = zero
      pars%bias(4) = zero
      
      !**************************************************
      !***   Read parameters from the external file   ***
      !**************************************************

      write(6,'(/)')
      write(6,'(1x,''Parameters used in calculation:'')')
      write(6,'(1x,59(''=''))')
      write(6,'(5x,''parameter'',5x,''dimension'',5x,''value(s)'')')
      write(6,'(1x,59(''-''))')

      ! read from the parameter file (ifile)
      do

         read(ifile,'(a)',iostat=io_status) line

         if (io_status.lt.0) exit
         if (io_status.gt.0) then
            write(*,'(/5x,''msg from setharm: error in reading external file'')')
            stop
         endif

         if (line.eq.' ') cycle

         if (line(1:1).eq.'#'.or.&
             line(1:1).eq.'*'.or.&
             line(1:1).eq.'C'.or.&
             line(1:1).eq.'!'.or.&
             line(1:1).eq.'c')      cycle

         ! Clean the input data
         do i = 1,len(line)
            if (line(i:i).eq.tab.or.line(i:i).eq.comma) line(i:i)=space
         enddo

         ! Initialize ISTART to interpret blanks as zero's
         do i=1,60
            istart(i)=80
         enddo

         ! Find initial digit of all numbers,
         ! check for leading spaces followed
         ! by a character and store in ISTART

         leadsp=.true.
         nvalue=0
         do i=1,len(line)
            if (leadsp.and.line(i:i).ne.space) then
               nvalue = nvalue + 1
               istart(nvalue) = i
            endif
            leadsp=(line(i:i).eq.space)
         enddo

         ! Parameter symbol is read

         string = line(istart(1):istart(2)-1)
         par_symbol = string(1:8)

         ! Check for error in parameter symbol

         numpar = 0
         do i=1,npars
            if (par_symbol.eq.pa(i)) then
               numpar=i
               exit
            endif
         enddo

         ! All O.K.

         if (numpar.eq.1) then
            do i=1,nvalue-1
               vpari(i) = reada(line,istart(1+i))
            enddo
         else
            vpar = reada(line,istart(2))
         endif

         ! Begin to set the parameter PA(NUMPAR)

         select case(numpar)

            case(0)
               write(6,'(//5x,''SETHARM: unrecognized parameter name:  <'',a,''>'')') par_symbol
               stop

            case(1)
                do i=1,4
                   pars%bias(i) = vpari(i)
                enddo

            case(2)
                pars%v12 = vpar

            case(3)
                pars%v34 = vpar
            
            case(4)
                pars%v13 = vpar

            case(5)
                pars%v24 = vpar

            case(6)
                pars%v14 = vpar
            
            case(7)
                pars%v23 = vpar

         end select

         if (numpar.eq.1) then
            do i=1,nvalue-1
               write(6,'(5x,a8,2x,a15,3x,f12.5)') par_symbol,dime(numpar),vpari(i)
            enddo
         else
            write(6,'(5x,a8,2x,a15,3x,f12.5)') par_symbol,dime(numpar),vpar
         endif

      enddo

      write(6,'(1x,59(''='')/)')

      return

   end subroutine setconst

   !======================================================================
   ! Reads parameters of the model harmonic potential
   !======================================================================
   subroutine setharm(ifile,pars)

      implicit none

      integer, intent(in) :: ifile
      type(harm_1d), intent(out) :: pars

      integer, parameter :: npars=11
      character(8),  dimension(npars) :: pa
      character(15), dimension(npars) :: dime

      integer :: io_status, i, nvalue, numpar
      real(8), dimension(10) :: vpari
      real(8)  :: vpar

      pa = (/ 'R0DH    ',&
              'R0AH    ',&
              'OMEGA_DH',&
              'OMEGA_AH',&
              'CORR    ',&
              'VPT_1   ',&
              'VPT_2   ',&
              'VET_A   ',&
              'VET_B   ',&
              'VEPT_1  ',&
              'VEPT_2  '   /)

      dime = (/ 'Angstrom       ',&
                'Angstrom       ',&
                'cm^(-1)        ',&
                'cm^(-1)        ',&
                'kcal/mol       ',&
                'kcal/mol       ',&
                'kcal/mol       ',&
                'kcal/mol       ',&
                'kcal/mol       ',&
                'kcal/mol       ',&
                'kcal/mol       '   /)

      !-- equilibrium distances
      pars%r0dh = 1.d0
      pars%r0ah = 1.d0

      !-- frequencies
      pars%omega_dh = 3000.d0
      pars%omega_ah = 3000.d0

      !-- PT coupling parameters VPT
      pars%vpt_1 = 50.d0
      pars%vpt_2 = 50.d0

      !-- ET coupling parameters VET
      pars%vet_a = 1.d0
      pars%vet_b = 1.d0

      !-- EPT mixed coupling parameters VEPT
      pars%vept_1 = zero
      pars%vept_2 = zero

      ! Constant correction terms
      pars%corr(1) = zero
      pars%corr(2) = zero
      pars%corr(3) = zero
      pars%corr(4) = zero
      
      !**************************************************
      !***   Read parameters from the external file   ***
      !**************************************************

      write(6,'(/)')
      write(6,'(1x,''Parameters used in calculation:'')')
      write(6,'(1x,59(''=''))')
      write(6,'(5x,''parameter'',5x,''dimension'',5x,''value(s)'')')
      write(6,'(1x,59(''-''))')

      ! read from the parameter file (ifile)
      do

         read(ifile,'(a)',iostat=io_status) line

         if (io_status.lt.0) exit
         if (io_status.gt.0) then
            write(*,'(/5x,''msg from setharm: error in reading external file'')')
            stop
         endif

         if (line.eq.' ') cycle

         if (line(1:1).eq.'#'.or.&
             line(1:1).eq.'*'.or.&
             line(1:1).eq.'C'.or.&
             line(1:1).eq.'!'.or.&
             line(1:1).eq.'c')      cycle

         ! Clean the input data
         do i = 1,len(line)
            if (line(i:i).eq.tab.or.line(i:i).eq.comma) line(i:i)=space
         enddo

         ! Initialize ISTART to interpret blanks as zero's
         do i=1,60
            istart(i)=80
         enddo

         ! Find initial digit of all numbers,
         ! check for leading spaces followed
         ! by a character and store in ISTART

         leadsp=.true.
         nvalue=0
         do i=1,len(line)
            if (leadsp.and.line(i:i).ne.space) then
               nvalue = nvalue + 1
               istart(nvalue) = i
            endif
            leadsp=(line(i:i).eq.space)
         enddo

         ! Parameter symbol is read

         string = line(istart(1):istart(2)-1)
         par_symbol = string(1:8)

         ! Check for error in parameter symbol

         numpar = 0
         do i=1,npars
            if (par_symbol.eq.pa(i)) then
               numpar=i
               exit
            endif
         enddo

         ! All O.K.

         if (numpar.eq.5) then
            do i=1,nvalue-1
               vpari(i) = reada(line,istart(1+i))
            enddo
         else
            vpar = reada(line,istart(2))
         endif

         ! Begin to set the parameter PA(NUMPAR)

         select case(numpar)

            case(0)
               write(6,'(//5x,''SETHARM: unrecognized parameter name:  <'',a,''>'')') par_symbol
               stop

            case(1)
                pars%r0dh = vpar
            
            case(2)
                pars%r0ah = vpar
            
            case(3)
                pars%omega_dh = vpar

            case(4)
                pars%omega_ah = vpar
            
            case(5)
                do i=1,4
                   pars%corr(i) = vpari(i)
                enddo

            case(6)
                pars%vpt_1 = vpar

            case(7)
                pars%vpt_2 = vpar
            
            case(8)
                pars%vet_a = vpar

            case(9)
                pars%vet_b = vpar

            case(10)
                pars%vept_1 = vpar
            
            case(11)
                pars%vept_2 = vpar

         end select

         if (numpar.eq.5) then
            do i=1,nvalue-1
               write(6,'(5x,a8,2x,a15,3x,f12.5)') par_symbol,dime(numpar),vpari(i)
            enddo
         else
            write(6,'(5x,a8,2x,a15,3x,f12.5)') par_symbol,dime(numpar),vpar
         endif

      enddo

      write(6,'(1x,59(''='')/)')

      return

   end subroutine setharm


   !======================================================================
   ! Reads parameters of the model MM5 potential (Nocera-like NHO systems)
   !======================================================================
   subroutine setmm5(ifile,pars)

      implicit none

      integer, intent(in) :: ifile
      type(mm5_1D), intent(out) :: pars

      integer, parameter :: npars=32
      character(8),  dimension(npars) :: pa
      character(15), dimension(npars) :: dime

      integer :: io_status, i, nvalue, numpar
      real(8), dimension(10) :: vpari
      real(8)  :: vpar

      pa = (/ 'DNH     ',&
              'DOH     ',&
              'BNH     ',&
              'BOH     ',&
              'RNH     ',&
              'ROH     ',&
              'DRNH    ',&
              'DROH    ',&
              'BRNH    ',&
              'BROH    ',&
              'CRNH    ',&
              'CROH    ',&
              'KSI     ',&
              'CORR    ',&
              'VPT0_1  ',&
              'GAPT_1  ',&
              'QPT0_1  ',&
              'VPT0_2  ',&
              'GAPT_2  ',&
              'QPT0_2  ',&
              'VET0_A  ',&
              'GAET_A  ',&
              'QET0_A  ',&
              'VET0_B  ',&
              'GAET_B  ',&
              'QET0_B  ',&
              'VETPT0_1',&
              'GAETPT_1',&
              'QETPT0_1',&
              'VETPT0_2',&
              'GAETPT_2',&
              'QETPT0_2'   /)

      dime = (/ 'kcal/mole      ',&
                'kcal/mole      ',&
                'kcal/mole/A    ',&
                'kcal/mole/A    ',&
                'Angstrom       ',&
                'Angstrom       ',&
                'kcal/mole      ',&
                'kcal/mole      ',&
                'kcal/mole/A    ',&
                'kcal/mole/A    ',&
                'kcal/mole*A**9 ',&
                'kcal/mole*A**9 ',&
                '               ',&
                'kcal/mole      ',&
                'kcal/mole      ',&
                'kcal/mole/A**2 ',&
                'Angstrom       ',&
                'kcal/mole      ',&
                'kcal/mole/A**2 ',&
                'Angstrom       ',&
                'kcal/mole      ',&
                'kcal/mole/A**2 ',&
                'Angstrom       ',&
                'kcal/mole      ',&
                'kcal/mole/A**2 ',&
                'Angstrom       ',&
                'kcal/mole      ',&
                'kcal/mole/A**2 ',&
                'Angstrom       ',&
                'kcal/mole      ',&
                'kcal/mole/A**2 ',&
                'Angstrom       '   /)

      ! Morse parameters
      pars%dnh = 180.d0
      pars%doh = 180.d0
      pars%bnh = 1.55d0
      pars%boh = 1.55d0
      pars%rnh = 1.d0
      pars%roh = 1.d0

      ! Repulsion potential parameters
      pars%drnh = 800.d0
      pars%droh = 800.d0
      pars%brnh = 1.45d0
      pars%broh = 1.45d0
      pars%crnh = zero
      pars%croh = zero

      ! Coulomb smoothing parameter
      pars%ksi = 1.d0

      ! Proton coupling parameters VPT0*EXP(-GAPT*(Q-QPT0))
      pars%vpt0_1 = 15.d0
      pars%gapt_1 = zero
      pars%qpt0_1 = zero
      pars%vpt0_2 = 15.d0
      pars%gapt_2 = zero
      pars%qpt0_2 = zero

      ! Electron coupling parameters VET0*EXP(-GAET*(Q-QET0))
      pars%vet0_a = 5.d0
      pars%gaet_a = zero
      pars%qet0_a = zero
      pars%vet0_b = 5.d0
      pars%gaet_b = zero
      pars%qet0_b = zero

      ! Mixed coupling parameters VETPT0*EXP(-GAETPT*(Q-QETPT0))
      pars%vetpt0_1 = 5.d0
      pars%gaetpt_1 = zero
      pars%qetpt0_1 = zero
      pars%vetpt0_2 = 5.d0
      pars%gaetpt_2 = zero
      pars%qetpt0_2 = zero

      ! Constant correction terms
      pars%corr(1) = zero
      pars%corr(2) = zero
      pars%corr(3) = zero
      pars%corr(4) = zero
      
      !**************************************************
      !***   Read parameters from the external file   ***
      !**************************************************

      write(6,'(/)')
      write(6,'(1x,''Parameters used in calculation:'')')
      write(6,'(1x,59(''=''))')
      write(6,'(5x,''parameter'',5x,''dimension'',5x,''value(s)'')')
      write(6,'(1x,59(''-''))')

      ! read from the parameter file (ifile)
      do

         read(ifile,'(a)',iostat=io_status) line

         if (io_status.lt.0) exit
         if (io_status.gt.0) then
            write(*,'(/5x,''msg from setharm: error in reading external file'')')
            stop
         endif

         if (line.eq.' ') cycle

         if (line(1:1).eq.'#'.or.&
             line(1:1).eq.'*'.or.&
             line(1:1).eq.'C'.or.&
             line(1:1).eq.'c')      cycle

         ! Clean the input data
         do i = 1,len(line)
            if (line(i:i).eq.tab.or.line(i:i).eq.comma) line(i:i)=space
         enddo

         ! Initialize ISTART to interpret blanks as zero's
         do i=1,60
            istart(i)=80
         enddo

         ! Find initial digit of all numbers,
         ! check for leading spaces followed
         ! by a character and store in ISTART

         leadsp=.true.
         nvalue=0
         do i=1,len(line)
            if (leadsp.and.line(i:i).ne.space) then
               nvalue = nvalue + 1
               istart(nvalue) = i
            endif
            leadsp=(line(i:i).eq.space)
         enddo

         ! Parameter symbol is read

         string = line(istart(1):istart(2)-1)
         par_symbol = string(1:8)

         ! Check for error in parameter symbol

         numpar = 0
         do i=1,npars
            if (par_symbol.eq.pa(i)) then
               numpar=i
               exit
            endif
         enddo

         ! All O.K.

         if (numpar.eq.14) then
            do i=1,nvalue-1
               vpari(i) = reada(line,istart(1+i))
            enddo
         else
            vpar = reada(line,istart(2))
         endif

         ! Begin to set the parameter PA(NUMPAR)

         select case(numpar)

            case(0)
               write(6,'(//5x,''SETPAR: unrecognized parameter name:  <'',a,''>'')') par_symbol
               stop

            case(1)
                pars%dnh = vpar
            
            case(2)
                pars%doh = vpar
            
            case(3)
                pars%bnh = vpar

            case(4)
                pars%boh = vpar
            
            case(5)
                pars%rnh = vpar
            
            case(6)
                pars%roh = vpar
            
            case(7)
                pars%drnh = vpar

            case(8)
                pars%droh = vpar
            
            case(9)
                pars%brnh = vpar
            
            case(10)
                pars%broh = vpar
            
            case(11)
                pars%crnh = vpar

            case(12)
                pars%croh = vpar

            case(13)
                pars%ksi  = vpar

            case(14)
                do i=1,4
                   pars%corr(i) = vpari(i)
                enddo

            case(15)
                pars%vpt0_1 = vpar
            
            case(16)
                pars%gapt_1 = vpar
            
            case(17)
                pars%qpt0_1 = vpar
            
            case(18)
                pars%vpt0_2 = vpar
            
            case(19)
                pars%gapt_2 = vpar

            case(20)
                pars%qpt0_2 = vpar
            
            case(21)
                pars%vet0_a = vpar

            case(22)
                pars%gaet_a = vpar
            
            case(23)
                pars%qet0_a = vpar
            
            case(24)
                pars%vet0_b = vpar

            case(25)
                pars%gaet_b = vpar
            
            case(26)
                pars%qet0_b = vpar
            
            case(27)
                pars%vetpt0_1 = vpar
            
            case(28)
                pars%gaetpt_1 = vpar

            case(29)
                pars%qetpt0_1 = vpar

            case(30)
                pars%vetpt0_2 = vpar

            case(31)
                pars%gaetpt_2 = vpar

            case(32)
                pars%qetpt0_2 = vpar

         end select

         if (numpar.eq.14) then
            do i=1,nvalue-1
               write(6,'(5x,a8,2x,a15,3x,f12.5)') par_symbol,dime(numpar),vpari(i)
            enddo
         else
            write(6,'(5x,a8,2x,a15,3x,f12.5)') par_symbol,dime(numpar),vpar
         endif

      enddo

      write(6,'(1x,59(''='')/)')

      return

   end subroutine setmm5

   !======================================================================
   ! Reads parameters of the MMGEN potential (general O-H...O systems)
   !======================================================================
   subroutine setmm5gen(ifile,pars)

      implicit none

      integer, intent(in) :: ifile
      type(mm5gen_1D), intent(out) :: pars

      integer, parameter :: npars=44
      character(8),  dimension(npars) :: pa
      character(15), dimension(npars) :: dime

      integer :: io_status, i, nvalue, numpar
      real(8), dimension(10) :: vpari
      real(8)  :: vpar

      pa = (/ 'DOH1A   ',&
              'DOH1B   ',&
              'DOH2A   ',&
              'DOH2B   ',&
              'BOH1A   ',&
              'BOH1B   ',&
              'BOH2A   ',&
              'BOH2B   ',&
              'ROH1A   ',&
              'ROH1B   ',&
              'ROH2A   ',&
              'ROH2B   ',&
              'DROH1A  ',&
              'DROH1B  ',&
              'DROH2A  ',&
              'DROH2B  ',&
              'BROH1A  ',&
              'BROH1B  ',&
              'BROH2A  ',&
              'BROH2B  ',&
              'CROH1A  ',&
              'CROH1B  ',&
              'CROH2A  ',&
              'CROH2B  ',&
              'KSI     ',&
              'CORR    ',&
              'VPT0_1  ',&
              'GAPT_1  ',&
              'QPT0_1  ',&
              'VPT0_2  ',&
              'GAPT_2  ',&
              'QPT0_2  ',&
              'VET0_A  ',&
              'GAET_A  ',&
              'QET0_A  ',&
              'VET0_B  ',&
              'GAET_B  ',&
              'QET0_B  ',&
              'VETPT0_1',&
              'GAETPT_1',&
              'QETPT0_1',&
              'VETPT0_2',&
              'GAETPT_2',&
              'QETPT0_2'/)

      dime = (/ 'kcal/mole      ',&
                'kcal/mole      ',&
                'kcal/mole      ',&
                'kcal/mole      ',&
                'kcal/mole/A    ',&
                'kcal/mole/A    ',&
                'kcal/mole/A    ',&
                'kcal/mole/A    ',&
                'Angstrom       ',&
                'Angstrom       ',&
                'Angstrom       ',&
                'Angstrom       ',&
                'kcal/mole      ',&
                'kcal/mole      ',&
                'kcal/mole      ',&
                'kcal/mole      ',&
                'kcal/mole/A    ',&
                'kcal/mole/A    ',&
                'kcal/mole/A    ',&
                'kcal/mole/A    ',&
                'kcal/mole*A**9 ',&
                'kcal/mole*A**9 ',&
                'kcal/mole*A**9 ',&
                'kcal/mole*A**9 ',&
                '               ',&
                'kcal/mole      ',&
                'kcal/mole      ',&
                'kcal/mole/A**2 ',&
                'Angstrom       ',&
                'kcal/mole      ',&
                'kcal/mole/A**2 ',&
                'Angstrom       ',&
                'kcal/mole      ',&
                'kcal/mole/A**2 ',&
                'Angstrom       ',&
                'kcal/mole      ',&
                'kcal/mole/A**2 ',&
                'Angstrom       ',&
                'kcal/mole      ',&
                'kcal/mole/A**2 ',&
                'Angstrom       ',&
                'kcal/mole      ',&
                'kcal/mole/A**2 ',&
                'Angstrom       '/)

      ! Morse parameters
      pars%doh1a = 100.d0
      pars%doh1b = 100.d0
      pars%doh2a = 100.d0
      pars%doh2b = 100.d0
      pars%boh1a = 1.55d0
      pars%boh1b = 1.55d0
      pars%boh2a = 1.55d0
      pars%boh2b = 1.55d0
      pars%roh1a = 1.d0
      pars%roh1b = 1.d0
      pars%roh2a = 1.d0
      pars%roh2b = 1.d0

      ! Repulsion potential parameters
      pars%droh1a = 800.d0
      pars%droh1b = 800.d0
      pars%droh2a = 800.d0
      pars%droh2b = 800.d0
      pars%broh1a = 1.45d0
      pars%broh1b = 1.45d0
      pars%broh2a = 1.45d0
      pars%broh2b = 1.45d0
      pars%croh1a = zero
      pars%croh1b = zero
      pars%croh2a = zero
      pars%croh2b = zero

      ! Coulomb smoothing parameter
      pars%ksi = 1.d0

      ! Proton coupling parameters VPT0*EXP(-GAPT*(Q-QPT0))
      pars%vpt0_1 = 15.d0
      pars%gapt_1 = zero
      pars%qpt0_1 = zero
      pars%vpt0_2 = 15.d0
      pars%gapt_2 = zero
      pars%qpt0_2 = zero

      ! Electron coupling parameters VET0*EXP(-GAET*(Q-QET0))
      pars%vet0_a = 5.d0
      pars%gaet_a = zero
      pars%qet0_a = zero
      pars%vet0_b = 5.d0
      pars%gaet_b = zero
      pars%qet0_b = zero

      ! Mixed coupling parameters VETPT0*EXP(-GAETPT*(Q-QETPT0))
      pars%vetpt0_1 = 5.d0
      pars%gaetpt_1 = zero
      pars%qetpt0_1 = zero
      pars%vetpt0_2 = 5.d0
      pars%gaetpt_2 = zero
      pars%qetpt0_2 = zero

      ! Constant correction terms
      pars%corr(1) = zero
      pars%corr(2) = zero
      pars%corr(3) = zero
      pars%corr(4) = zero
      
      !**************************************************
      !***   Read parameters from the external file   ***
      !**************************************************

      write(6,'(/)')
      write(6,'(1x,''Parameters used in calculation:'')')
      write(6,'(1x,59(''=''))')
      write(6,'(5x,''parameter'',5x,''dimension'',5x,''value(s)'')')
      write(6,'(1x,59(''-''))')

      ! read from the parameter file (ifile)
      do

         read(ifile,'(a)',iostat=io_status) line

         if (io_status.lt.0) exit
         if (io_status.gt.0) then
            write(*,'(/5x,''msg from setharm: error in reading external file'')')
            stop
         endif

         if (line.eq.' ') cycle

         if (line(1:1).eq.'#'.or.&
             line(1:1).eq.'*'.or.&
             line(1:1).eq.'C'.or.&
             line(1:1).eq.'c')      cycle

         ! Clean the input data
         do i = 1,len(line)
            if (line(i:i).eq.tab.or.line(i:i).eq.comma) line(i:i)=space
         enddo

         ! Initialize ISTART to interpret blanks as zero's
         do i=1,60
            istart(i)=80
         enddo

         ! Find initial digit of all numbers,
         ! check for leading spaces followed
         ! by a character and store in ISTART

         leadsp=.true.
         nvalue=0
         do i=1,len(line)
            if (leadsp.and.line(i:i).ne.space) then
               nvalue = nvalue + 1
               istart(nvalue) = i
            endif
            leadsp=(line(i:i).eq.space)
         enddo

         ! Parameter symbol is read

         string = line(istart(1):istart(2)-1)
         par_symbol = string(1:8)

         ! Check for error in parameter symbol

         numpar = 0
         do i=1,npars
            if (par_symbol.eq.pa(i)) then
               numpar=i
               exit
            endif
         enddo

         ! All O.K.

         if (numpar.eq.26) then
            do i=1,nvalue-1
               vpari(i) = reada(line,istart(1+i))
            enddo
         else
            vpar = reada(line,istart(2))
         endif

         ! Begin to set the parameter PA(NUMPAR)

         select case(numpar)

            case(0)
               write(6,'(//5x,''SETPAR: unrecognized parameter name:  <'',a,''>'')') par_symbol
               stop

            case(1)
                pars%doh1a = vpar
            
            case(2)
                pars%doh1b = vpar

            case(3)
                pars%doh2a = vpar
            
            case(4)
                pars%doh2b = vpar
            
            case(5)
                pars%boh1a = vpar

            case(6)
                pars%boh1b = vpar
            
            case(7)
                pars%boh2a = vpar

            case(8)
                pars%boh2b = vpar
            
            case(9)
                pars%roh1a = vpar
            
            case(10)
                pars%roh1b = vpar
            
            case(11)
                pars%roh2a = vpar
            
            case(12)
                pars%roh2b = vpar
            
            case(13)
                pars%droh1a = vpar

            case(14)
                pars%droh1b = vpar
            
            case(15)
                pars%droh2a = vpar

            case(16)
                pars%droh2b = vpar
            
            case(17)
                pars%broh1a = vpar
            
            case(18)
                pars%broh1b = vpar
            
            case(19)
                pars%broh2a = vpar
            
            case(20)
                pars%broh2b = vpar
            
            case(21)
                pars%croh1a = vpar

            case(22)
                pars%croh1b = vpar

            case(23)
                pars%croh2a = vpar

            case(24)
                pars%croh2b = vpar

            case(25)
                pars%ksi  = vpar

            case(26)
                do i=1,4
                   pars%corr(i) = vpari(i)
                enddo

            case(27)
                pars%vpt0_1 = vpar
            
            case(28)
                pars%gapt_1 = vpar
            
            case(29)
                pars%qpt0_1 = vpar
            
            case(30)
                pars%vpt0_2 = vpar
            
            case(31)
                pars%gapt_2 = vpar

            case(32)
                pars%qpt0_2 = vpar
            
            case(33)
                pars%vet0_a = vpar

            case(34)
                pars%gaet_a = vpar
            
            case(35)
                pars%qet0_a = vpar
            
            case(36)
                pars%vet0_b = vpar

            case(37)
                pars%gaet_b = vpar
            
            case(38)
                pars%qet0_b = vpar
            
            case(39)
                pars%vetpt0_1 = vpar
            
            case(40)
                pars%gaetpt_1 = vpar

            case(41)
                pars%qetpt0_1 = vpar

            case(42)
                pars%vetpt0_2 = vpar

            case(43)
                pars%gaetpt_2 = vpar

            case(44)
                pars%qetpt0_2 = vpar

         end select

         if (numpar.eq.26) then
            do i=1,nvalue-1
               write(6,'(5x,a8,2x,a15,3x,f12.5)') par_symbol,dime(numpar),vpari(i)
            enddo
         else
            write(6,'(5x,a8,2x,a15,3x,f12.5)') par_symbol,dime(numpar),vpar
         endif

      enddo

      write(6,'(1x,59(''='')/)')

      return

   end subroutine setmm5gen

   !======================================================================
   ! Reads parameters of the model 2D LEPS potential
   !======================================================================
   subroutine setleps(ifile,pars)
   
      implicit none

      integer, intent(in) :: ifile
      type(leps_2d), intent(out) :: pars

      integer, parameter :: npars=30
      character(8),  dimension(npars) :: pa
      character(15), dimension(npars) :: dime

      integer :: io_status, i, nvalue, numpar
      real(8), dimension(10) :: vpari
      real(8)  :: vpar, fd1, fd1au, vpark, fd2, fd2au, fa1, fa1au, fa2, fa2au

      pa = (/ 'DAH     ',&
              'DBH     ',&
              'DAB     ',&
              'BAH     ',&
              'BBH     ',&
              'BAB     ',&
              'RAH     ',&
              'RBH     ',&
              'RAB     ',&
              'R0Dd1   ',&
              'R0Dd2   ',&
	      'R0Aa1   ',&
              'R0Aa2   ',&
	      'KAH     ',&
              'KBH     ',&
              'KAB     ',&
              'SIGMA   ',&
              'VPT0    ',&
              'ALPH    ',&
              'VET0    ',&
              'BETA    ',&
              'CORR    ',&
              'FD1     ',&
              'FD2     ',&
              'FA1     ',&
              'FA2     ',&
              'KD1     ',&
              'KD2     ',&
              'KA1     ',&
              'KA2     '   /)

      dime = (/  'kcal/mole      ',&
                 'kcal/mole      ',&
                 'kcal/mole      ',&
                 'kcal/mole/A    ',&
                 'kcal/mole/A    ',&
                 'kcal/mole/A    ',&
                 'Angstrom       ',&
                 'Angstrom       ',&
                 'Angstrom       ',&
                 'Angstrom       ',&
                 'Angstrom       ',&
                 'Angstrom       ',&
                 'Angstrom       ',&
                 '               ',&
                 '               ',&
                 '               ',&
                 '               ',&
                 'kcal/mole      ',&
                 '1/Angstrom     ',&
                 'kcal/mole      ',&
                 '1/Angstrom     ',&
                 'kcal/mole      ',&
                 '1/cm           ',&
                 '1/cm           ',&
                 '1/cm           ',&
                 '1/cm           ',&
                 'kcal/mole/A^2  ',&
                 'kcal/mole/A^2  ',&
                 'kcal/mole/A^2  ',&
                 'kcal/mole/A^2  '    /)

      ! Morse parameters
      pars%dah = 180.d0
      pars%dbh = 180.d0
      pars%dab = 180.d0
      pars%bah = 1.55d0
      pars%bbh = 1.55d0
      pars%bab = 1.55d0
      pars%rah = 1.d0
      pars%rbh = 1.d0
      pars%rab = 1.d0
      pars%r0dd1 = 1.d0
      pars%r0aa1 = 1.d0
      pars%r0dd2 = 1.d0
      pars%r0aa2 = 1.d0

      ! Sato Parameters
      pars%kah = 0.5d0
      pars%kbh = 0.5d0
      pars%kab = 0.5d0

      ! VPT parameters
      pars%vpt0 = 1.d0
      pars%alph = 1.d0

      ! Sigma, ET coupling
      pars%sigma = 0.5d0
      pars%vet0  = 0.1d0
      pars%beta  = zero

      ! Constant correction terms
      pars%corr(1) = zero
      pars%corr(2) = zero
      pars%corr(3) = zero
      pars%corr(4) = zero

      ! Force constants for ET-PT donors (acceptors)
      pars%kd1 = 100.d0
      pars%kd2 = 100.d0
      pars%ka1 = 100.d0
      pars%ka2 = 100.d0


      !**************************************************
      !***   Read parameters from the external file   ***
      !**************************************************

      write(6,'(/)')
      write(6,'(1x,''Parameters used in calculation:'')')
      write(6,'(1x,59(''=''))')
      write(6,'(5x,''parameter'',5x,''dimension'',5x,''value(s)'')')
      write(6,'(1x,59(''-''))')

      ! read from the parameter file (ifile)
      do

         read(ifile,'(a)',iostat=io_status) line

         if (io_status.lt.0) exit
         if (io_status.gt.0) then
            write(*,'(/5x,''msg from setleps: error in reading external file'')')
            stop
         endif

         if (line.eq.' ') cycle

         if (line(1:1).eq.'#'.or.&
             line(1:1).eq.'*'.or.&
             line(1:1).eq.'C'.or.&
             line(1:1).eq.'c')      cycle

         ! Clean the input data
         do i = 1,len(line)
            if (line(i:i).eq.tab.or.line(i:i).eq.comma) line(i:i)=space
         enddo

         ! Initialize ISTART to interpret blanks as zero's
         do i=1,60
            istart(i)=80
         enddo

         ! Find initial digit of all numbers,
         ! check for leading spaces followed
         ! by a character and store in ISTART

         leadsp=.true.
         nvalue=0
         do i=1,len(line)
            if (leadsp.and.line(i:i).ne.space) then
               nvalue = nvalue + 1
               istart(nvalue) = i
            endif
            leadsp=(line(i:i).eq.space)
         enddo

         ! Parameter symbol is read

         string = line(istart(1):istart(2)-1)
         par_symbol = string(1:8)

         ! Check for error in parameter symbol

         numpar = 0
         do i=1,npars
            if (par_symbol.eq.pa(i)) then
               numpar=i
               exit
            endif
         enddo

         ! All O.K.

         if (numpar.eq.22) then
            do i=1,nvalue-1
               vpari(i) = reada(line,istart(1+i))
            enddo
         else
            vpar = reada(line,istart(2))
         endif

         ! Begin to set the parameter PA(NUMPAR)

         select case(numpar)

            case(0)
               write(6,'(//5x,''SETPAR: unrecognized parameter name:  <'',a,''>'')') par_symbol
               stop

            case(1)
                pars%dah = vpar

            case(2)
                pars%dbh = vpar

            case(3)
                pars%dab = vpar

            case(4)
                pars%bah = vpar

            case(5)
                pars%bbh = vpar

            case(6)
                pars%bab = vpar

            case(7)
                pars%rah = vpar

            case(8)
                pars%rbh = vpar

            case(9)
                pars%rab = vpar

            case(10)
                pars%r0dd1 = vpar
            
	    case(11)
                pars%r0dd2 = vpar

            case(12)
                pars%r0aa1 = vpar
            
	    case(13)
                pars%r0aa2 = vpar

            case(14)
                pars%kah = vpar

            case(15)
                pars%kbh  = vpar

            case(16)
                pars%kab  = vpar

            case(17)
                pars%sigma = vpar

            case(18)
                pars%vpt0 = vpar

            case(19)
                pars%alph = vpar

            case(20)
                pars%vet0 = vpar

            case(21)
                pars%beta = vpar

            case(22)
                do i=1,4
                   pars%corr(i) = vpari(i)
                enddo

            case(23)
                fd1 = vpar
                fd1au = fd1*cm2ev*ev2au
                pars%kd1 = dm*fd1au*fd1au*au2cal/bohr2a/bohr2a
                vpark = pars%kd1

            case(24)
                fd2 = vpar
                fd2au = fd2*cm2ev*ev2au
                pars%kd2 = dm*fd2au*fd2au*au2cal/bohr2a/bohr2a
                vpark = pars%kd2

            case(25)
                fa1 = vpar
                fa1au = fa1*cm2ev*ev2au
                pars%ka1 = am*fa1au*fa1au*au2cal/bohr2a/bohr2a
                vpark = pars%ka1

            case(26)
                fa2 = vpar
                fa2au = fa2*cm2ev*ev2au
                pars%ka2 = am*fa2au*fa2au*au2cal/bohr2a/bohr2a
                vpark = pars%ka2

            case(27)
                pars%kd1 = vpar

            case(28)
                pars%kd2 = vpar

            case(29)
                pars%ka1 = vpar

            case(30)
                pars%ka2 = vpar

         end select

         if (numpar.eq.22) then

            do i=1,nvalue-1
               write(6,'(5x,a8,2x,a15,3x,f12.5)') par_symbol,dime(numpar),vpari(i)
            enddo

         elseif (numpar.ge.23.and.numpar.le.26) then

            write(6,'(5x,a8,2x,a15,3x,f12.5)') par_symbol,dime(numpar),vpar
            write(6,'(5x,a8,2x,a15,3x,f12.5)') 'K'//par_symbol(2:),'kcal/mol/A^2   ',vpark

         else

            write(6,'(5x,a8,2x,a15,3x,f12.5)') par_symbol,dime(numpar),vpar

         endif
      
      enddo

      write(6,'(1x,59(''='')/)')

      return

   end subroutine setleps

   !======================================================================
   !     Reads parameters of the model hybrid potential
   !======================================================================
   subroutine sethyb(ifile,pars)

      implicit none

      integer, intent(in) :: ifile
      type(hybrid_2d), intent(out) :: pars

      integer, parameter :: npars=18
      character(8),  dimension(npars) :: pa
      character(15), dimension(npars) :: dime

      integer :: io_status, i, nvalue, numpar
      real(8), dimension(10) :: vpari
      real(8)  :: vpar

      pa = (/ 'OMEGA_AH',&
              'OMEGA_BH',&
              'OMEGA_AB',&
              'DAH     ',&
              'DBH     ',&
              'RAH     ',&
              'RBH     ',&
              'RAB1    ',&
              'RAB2    ',&
              'KAH     ',&
              'KBH     ',&
              'SIGMA   ',&
              'VPT0    ',&
              'ALPH    ',&
              'VET0    ',&
              'BETA    ',&
              'VEPT0   ',&
              'CORR    '   /)

      dime = (/  '1/cm                ',&
                 '1/cm                ',&
                 '1/cm                ',&
                 'kcal/mole           ',&
                 'kcal/mole           ',&
                 'Angstroem           ',&
                 'Angstroem           ',&
                 'Angstroem           ',&
                 'Angstroem           ',&
                 '                    ',&
                 '                    ',&
                 '                    ',&
                 'kcal/mole           ',&
                 '1/Angstroem         ',&
                 'kcal/mole           ',&
                 '1/Angstroem         ',&
                 'kcal/mole           ',&
                 'kcal/mole           '    /)


      !==== gas phase hybrid potential parameters (kcal/mole, Angstroms)
      !     (default values)

      !- Frequencies
      pars%omega_ah = 3000.d0
      pars%omega_bh = 3000.d0
      pars%omega_ab =  300.d0

      !- bond dissociation energies
      pars%dah  = 100.d0
      pars%dbh  = 100.d0

      !- Equilibrium distances
      pars%rah  = 1.d0
      pars%rbh  = 1.d0
      pars%rab1 = 3.d0
      pars%rab2 = 3.d0

      !- scaling Sato factors and overlap integral
      pars%kah   = half
      pars%kbh   = half
      pars%sigma = half

      !- VPT coupling parameters
      pars%vpt0 = one
      pars%alph = one

      !- ET coupling parameters
      pars%vet0  = 0.1d0
      pars%beta  = zero

      !- EPT coupling parameters
      pars%vept0  = zero

      !- Constant correction terms (energy biases)
      pars%corr(1) = zero
      pars%corr(2) = zero
      pars%corr(3) = zero
      pars%corr(4) = zero


      !**************************************************
      !***   Read parameters from the external file   ***
      !**************************************************

      write(6,'(/)')
      write(6,'(1x,''Parameters of the hybrid potential:'')')
      write(6,'(1x,59(''=''))')
      write(6,'(5x,''parameter'',5x,''dimension'',5x,''value(s)'')')
      write(6,'(1x,59(''-''))')

      ! read from the parameter file (ifile)
      do

         read(ifile,'(a)',iostat=io_status) line

         if (io_status.lt.0) exit
         if (io_status.gt.0) then
            write(*,'(/5x,''msg from SETHYB: error in reading external file'')')
            stop
         endif

         if (line.eq.' ') cycle

         if (line(1:1).eq.'#'.or.&
             line(1:1).eq.'*'.or.&
             line(1:1).eq.'C'.or.&
             line(1:1).eq.'c')      cycle

         ! Clean the input data
         do i = 1,len(line)
            if (line(i:i).eq.tab.or.line(i:i).eq.comma) line(i:i)=space
         enddo

         ! Initialize ISTART to interpret blanks as zero's
         do i=1,60
            istart(i)=80
         enddo

         ! Find initial digit of all numbers,
         ! check for leading spaces followed
         ! by a character and store in ISTART

         leadsp=.true.
         nvalue=0
         do i=1,len(line)
            if (leadsp.and.line(i:i).ne.space) then
               nvalue = nvalue + 1
               istart(nvalue) = i
            endif
            leadsp=(line(i:i).eq.space)
         enddo

         ! Parameter symbol is read

         string = line(istart(1):istart(2)-1)
         par_symbol = string(1:8)

         ! Check for error in parameter symbol

         numpar = 0
         do i=1,npars
            if (par_symbol.eq.pa(i)) then
               numpar=i
               exit
            endif
         enddo

         ! All O.K.

         if (numpar.eq.18) then
            do i=1,nvalue-1
               vpari(i) = reada(line,istart(1+i))
            enddo
         else
            vpar = reada(line,istart(2))
         endif

         ! Begin to set the parameter PA(NUMPAR)

         select case(numpar)

            case(0)
               write(6,'(//5x,''SETPAR: unrecognized parameter name:  <'',a,''>'')') par_symbol
               stop

            case(1)
               pars%omega_ah = vpar

            case(2)
               pars%omega_bh = vpar

            case(3)
               pars%omega_ab = vpar

            case(4)
               pars%dah = vpar

            case(5)
               pars%dbh = vpar

            case(6)
               pars%rah = vpar

            case(7)
               pars%rbh = vpar

            case(8)
               pars%rab1 = vpar

            case(9)
               pars%rab2 = vpar

            case(10)
               pars%kah = vpar

            case(11)
               pars%kbh = vpar

            case(12)
               pars%sigma = vpar

            case(13)
               pars%vpt0 = vpar

            case(14)
               pars%alph = vpar

            case(15)
               pars%vet0 = vpar

            case(16)
               pars%beta = vpar

            case(17)
               pars%vept0 = vpar

            case(18)
               do i=1,4
                  pars%corr(i) = vpari(i)
               enddo

         end select

         if (numpar.eq.18) then
            do i=1,nvalue-1
               write(6,'(5x,a8,2x,a15,3x,f12.5)') par_symbol,dime(numpar),vpari(i)
            enddo
         else
            write(6,'(5x,a8,2x,a15,3x,f12.5)') par_symbol,dime(numpar),vpar
         endif

      enddo

      write(6,'(1x,59(''='')/)')

      return

   end subroutine sethyb


   !-- Gas-phase Hamiltonians --------------------------------------------------


   subroutine h0mat_constant(h0_)
   !======================================================================!
   ! Calculates the gas phase Hamiltonian matrix using
   ! constant potential for a fixed solute (kcal/mole)
   !======================================================================!
      implicit none

      real(8), intent(out), dimension(4,4) :: h0_

      integer :: pdgas, nhgas, pagas, i, ia, ja
      real(8)  :: dp, dp2, q, ua, ub, odh, oah, rdh, rah

      !-- local variables: parameters
      real(8) :: v12
      real(8) :: v34
      real(8) :: v13
      real(8) :: v24
      real(8) :: v14
      real(8) :: v23
      real(8), dimension(4) :: bias

      !-- pointer to the real set of parameters
      type(constant), pointer :: par

      par => constpar

      v12  = par%v12
      v34  = par%v34
      v13  = par%v13
      v24  = par%v24
      v14  = par%v14
      v23  = par%v23
      bias = par%bias

      h0_ = zero

      !==== DIAGONAL ELEMENTS

      h0_(1,1) = bias(1)
      h0_(2,2) = bias(2)
      h0_(3,3) = bias(3)
      h0_(4,4) = bias(4)

      !==== OFF-DIAGONAL ELEMENTS

      !-- PT couplings

      h0_(1,2) = v12
      h0_(2,1) = v12
      h0_(3,4) = v34
      h0_(4,3) = v34

      !-- ET couplings

      h0_(1,3) = v13
      h0_(3,1) = v13
      h0_(2,4) = v24
      h0_(4,2) = v24

      !-- Mixed EPT couplings

      h0_(1,4) = v14
      h0_(4,1) = v14
      h0_(2,3) = v23
      h0_(3,2) = v23

      return

   end subroutine h0mat_constant


   subroutine h0mat_harm(h0_)
   !======================================================================!
   ! Calculates the gas phase Hamiltonian matrix using
   ! harmonic potential for general D---H---A systems (kcal/mole)
   !======================================================================!
      implicit none

      real(8), intent(out), dimension(4,4) :: h0_

      integer :: pdgas, nhgas, pagas, i, ia, ja
      real(8)  :: dp, dp2, q, ua, ub, odh, oah, rdh, rah

      !-- local variables: parameters
      real(8) :: r0dh
      real(8) :: r0ah
      real(8) :: omega_dh
      real(8) :: omega_ah
      real(8) :: vpt_1
      real(8) :: vpt_2
      real(8) :: vet_a
      real(8) :: vet_b
      real(8) :: vept_1
      real(8) :: vept_2
      real(8), dimension(4) :: corr

      !-- pointer to the real set of parameters
      type(harm_1d), pointer :: par

      par => harmpar

      r0dh     = par%r0dh
      r0ah     = par%r0ah
      omega_dh = par%omega_dh
      omega_ah = par%omega_ah
      vpt_1    = par%vpt_1
      vpt_2    = par%vpt_2
      vet_a    = par%vet_a
      vet_b    = par%vet_b
      vept_1   = par%vept_1
      vept_2   = par%vept_2
      corr     = par%corr

      h0_ = zero
      pdgas = iptgas(1)
      nhgas = iptgas(2)
      pagas = iptgas(3)

      !-- DP - a distance between PT donor (PDGAS) and acceptor (PAGAS)
      !   Note that PT interface is oriented along the x-axis.
      !   Q - antisymmetrical proton coordinate

      dp  = dabs(xyzgas(1,pagas)-xyzgas(1,pdgas))
      dp2 = dp/2.d0
      q   = xyzgas(1,nhgas)

      rdh = (dp2 + q - r0dh)*a2bohr
      rah = (dp2 - q - r0ah)*a2bohr

      odh = omega_dh*cm2ev*ev2au
      oah = omega_ah*cm2ev*ev2au

      !==== DIAGONAL ELEMENTS

      !-- Harmonic potentials for proton transfer EVB states

      ua = half*pm*odh*odh*rdh*rdh*au2cal
      ub = half*pm*oah*oah*rah*rah*au2cal

      h0_(1,1) = ua + corr(1)
      h0_(2,2) = ub + corr(2)
      h0_(3,3) = ua + corr(3)
      h0_(4,4) = ub + corr(4)

      !==== OFF-DIAGONAL ELEMENTS

      !-- PT couplings

      h0_(1,2) = vpt_1
      h0_(2,1) = vpt_1
      h0_(3,4) = vpt_2
      h0_(4,3) = vpt_2

      !-- ET couplings

      h0_(1,3) = vet_a
      h0_(3,1) = vet_a
      h0_(2,4) = vet_b
      h0_(4,2) = vet_b

      !-- Mixed EPT couplings

      h0_(1,4) = vept_1
      h0_(4,1) = vept_1
      h0_(2,3) = vept_2
      h0_(3,2) = vept_2

      return

   end subroutine h0mat_harm


   subroutine dh0mat_harm(dh0_)
   !======================================================================!
   ! Calculates the derivative of the gas phase Hamiltonian matrix using
   ! harmonic potential for general D---H---A systems (kcal/mole)
   !======================================================================!
      implicit none

      real(8), intent(out), dimension(4,4) :: dh0_

      integer :: pdgas, nhgas, pagas, i, ia, ja
      real(8)  :: dp, dp2, q, dua, dub, odh, oah, rdh, rah

      !-- local variables: parameters
      real(8) :: r0dh
      real(8) :: r0ah
      real(8) :: omega_dh
      real(8) :: omega_ah
      real(8) :: vpt_1
      real(8) :: vpt_2
      real(8) :: vet_a
      real(8) :: vet_b
      real(8) :: vept_1
      real(8) :: vept_2
      real(8), dimension(4) :: corr

      !-- pointer to the real set of parameters
      type(harm_1d), pointer :: par

      par => harmpar

      r0dh     = par%r0dh
      r0ah     = par%r0ah
      omega_dh = par%omega_dh
      omega_ah = par%omega_ah
      vpt_1    = par%vpt_1
      vpt_2    = par%vpt_2
      vet_a    = par%vet_a
      vet_b    = par%vet_b
      vept_1   = par%vept_1
      vept_2   = par%vept_2
      corr     = par%corr

      dh0_ = zero
      pdgas = iptgas(1)
      nhgas = iptgas(2)
      pagas = iptgas(3)

      !-- DP - a distance between PT donor (PDGAS) and acceptor (PAGAS)
      !   Note that PT interface is oriented along the x-axis.
      !   Q - antisymmetrical proton coordinate

      dp  = dabs(xyzgas(1,pagas)-xyzgas(1,pdgas))
      dp2 = dp/2.d0
      q   = xyzgas(1,nhgas)

      rdh = (dp2 + q - r0dh)*a2bohr
      rah = (dp2 - q - r0ah)*a2bohr

      odh = omega_dh*cm2ev*ev2au
      oah = omega_ah*cm2ev*ev2au

      !==== derivatives of the DIAGONAL ELEMENTS

      !-- Harmonic potentials for proton transfer EVB states

      dua = pm*odh*odh*rdh
      dub = pm*oah*oah*rah

      dh0_(1,1) = dua
      dh0_(2,2) = dub
      dh0_(3,3) = dua
      dh0_(4,4) = dub

      !==== OFF-DIAGONAL ELEMENTS are zero (constant couplings)

      return

   end subroutine dh0mat_harm


   subroutine d2h0mat_harm(d2h0_)
   !======================================================================!
   ! Calculates the second derivative of the Hamiltonian matrix using
   ! harmonic potential for general D---H---A systems (kcal/mole)
   !======================================================================!
      implicit none

      real(8), intent(out), dimension(4,4) :: d2h0_

      integer :: pdgas, nhgas, pagas, i, ia, ja
      real(8)  :: dp, dp2, q, d2ua, d2ub, odh, oah, rdh, rah

      !-- local variables: parameters
      real(8) :: r0dh
      real(8) :: r0ah
      real(8) :: omega_dh
      real(8) :: omega_ah
      real(8) :: vpt_1
      real(8) :: vpt_2
      real(8) :: vet_a
      real(8) :: vet_b
      real(8) :: vept_1
      real(8) :: vept_2
      real(8), dimension(4) :: corr

      !-- pointer to the real set of parameters
      type(harm_1d), pointer :: par

      par => harmpar

      r0dh     = par%r0dh
      r0ah     = par%r0ah
      omega_dh = par%omega_dh
      omega_ah = par%omega_ah
      vpt_1    = par%vpt_1
      vpt_2    = par%vpt_2
      vet_a    = par%vet_a
      vet_b    = par%vet_b
      vept_1   = par%vept_1
      vept_2   = par%vept_2
      corr     = par%corr

      d2h0_ = zero
      pdgas = iptgas(1)
      nhgas = iptgas(2)
      pagas = iptgas(3)

      !-- DP - a distance between PT donor (PDGAS) and acceptor (PAGAS)
      !   Note that PT interface is oriented along the x-axis.
      !   Q - antisymmetrical proton coordinate

      dp  = dabs(xyzgas(1,pagas)-xyzgas(1,pdgas))
      dp2 = dp/2.d0
      q   = xyzgas(1,nhgas)

      rdh = (dp2 + q - r0dh)*a2bohr
      rah = (dp2 - q - r0ah)*a2bohr

      odh = omega_dh*cm2ev*ev2au
      oah = omega_ah*cm2ev*ev2au

      !==== second derivatives of the DIAGONAL ELEMENTS

      !-- Harmonic potentials for proton transfer EVB states

      d2ua = pm*odh*odh
      d2ub = pm*oah*oah

      d2h0_(1,1) = d2ua
      d2h0_(2,2) = d2ub
      d2h0_(3,3) = d2ua
      d2h0_(4,4) = d2ub

      !==== OFF-DIAGONAL ELEMENTS are zero (constant couplings)

      return

   end subroutine d2h0mat_harm


   subroutine h0mat_mm5(h0_)
   !======================================================================!
   ! Calculates the gas phase Hamiltonian matrix using
   ! MM5 potential for Nocera like systems (kcal/mole)
   !======================================================================!
      implicit none

      real(8), intent(out), dimension(4,4) :: h0_

      integer :: pdgas, nhgas, pagas, i, ia, ja
      real(8)  :: dp, dp2, q, vma, vmb, vc1a, vc1b, vc2a, vc2b, dxc, vcci
      real(8)  :: xoh, xnh, vra, vrb, vpt1, vpt2, veta, vetb, vetpt1, vetpt2
      real(8), dimension(4) :: vcc

      ! local variables: parameters
      real(8) :: dnh
      real(8) :: doh
      real(8) :: bnh
      real(8) :: boh
      real(8) :: rnh
      real(8) :: roh
      real(8) :: drnh
      real(8) :: droh
      real(8) :: brnh
      real(8) :: broh
      real(8) :: croh
      real(8) :: crnh
      real(8) :: ksi
      real(8) :: vpt0_1
      real(8) :: gapt_1
      real(8) :: qpt0_1
      real(8) :: vpt0_2
      real(8) :: gapt_2
      real(8) :: qpt0_2
      real(8) :: vet0_a
      real(8) :: gaet_a
      real(8) :: qet0_a
      real(8) :: vet0_b
      real(8) :: gaet_b
      real(8) :: qet0_b
      real(8) :: vetpt0_1
      real(8) :: gaetpt_1
      real(8) :: qetpt0_1
      real(8) :: vetpt0_2
      real(8) :: gaetpt_2
      real(8) :: qetpt0_2
      real(8), dimension(4) :: corr

      ! pointer to the real set of parameters
      type(mm5_1d), pointer :: par

      par => mm5par

      dnh      = par%dnh
      doh      = par%doh
      bnh      = par%bnh
      boh      = par%boh
      rnh      = par%rnh
      roh      = par%roh
      drnh     = par%drnh
      droh     = par%droh
      brnh     = par%brnh
      broh     = par%broh
      croh     = par%croh
      crnh     = par%crnh
      ksi      = par%ksi
      vpt0_1   = par%vpt0_1
      gapt_1   = par%gapt_1
      qpt0_1   = par%qpt0_1
      vpt0_2   = par%vpt0_2
      gapt_2   = par%gapt_2
      qpt0_2   = par%qpt0_2
      vet0_a   = par%vet0_a
      gaet_a   = par%gaet_a
      qet0_a   = par%qet0_a
      vet0_b   = par%vet0_b
      gaet_b   = par%gaet_b
      qet0_b   = par%qet0_b
      vetpt0_1 = par%vetpt0_1
      gaetpt_1 = par%gaetpt_1
      qetpt0_1 = par%qetpt0_1
      vetpt0_2 = par%vetpt0_2
      gaetpt_2 = par%gaetpt_2
      qetpt0_2 = par%qetpt0_2
      corr     = par%corr

      h0_ = zero
      pdgas = iptgas(1)
      nhgas = iptgas(2)
      pagas = iptgas(3)

      ! DP - a distance between PT donor (PDGAS) and acceptor (PAGAS)
      ! Note that PT interface is oriented along the x-axis.
      ! Q - antisymmetrical proton coordinate

      dp  = dabs(xyzgas(1,pagas)-xyzgas(1,pdgas))
      dp2 = dp/2.d0
      q   = xyzgas(1,nhgas)

      !==== DIAGONAL ELEMENTS

      ! Morse potentials for proton transfer EVB states

      vma = dnh*(1.d0-dexp(-bnh*(q+dp2-rnh)))**2.d0
      vmb = doh*(1.d0-dexp(-boh*(dp2-q-roh)))**2.d0

      ! Coulomb interaction terms

      vc1a = zero
      vc1b = zero
      vc2a = zero
      vc2b = zero
      do i=1,4
         vcc(i) = zero
      enddo

      if (coulomb) then

         do i=1,natgas
            if (i.ne.pdgas.and.i.ne.nhgas) then
               dxc = distance(nhgas,i,xyzgas)
               vc1a = vc1a + chrgas(1,nhgas)*chrgas(1,i)/dxc
            endif
         enddo
         vc1a = vc1a*e2*ev2cal

         do i=1,natgas
            if (i.ne.pagas.and.i.ne.nhgas) then
               dxc = distance(nhgas,i,xyzgas)
               vc1b = vc1b + chrgas(2,nhgas)*chrgas(2,i)/dxc
            endif
         enddo
         vc1b = vc1b*e2*ev2cal

         do i=1,natgas
            if (i.ne.pdgas.and.i.ne.nhgas) then
               dxc = distance(nhgas,i,xyzgas)
               vc2a = vc2a + chrgas(3,nhgas)*chrgas(3,i)/dxc
            endif
         enddo
         vc2a = vc2a*e2*ev2cal

         do i=1,5
            if (i.ne.pagas.and.i.ne.nhgas) then
               dxc = distance(nhgas,i,xyzgas)
               vc2b = vc2b + chrgas(4,nhgas)*chrgas(4,i)/dxc
            endif
         enddo
         vc2b = vc2b*e2*ev2cal

         ! Constant Coulomb interaction terms

         do i=1,nelst
            vcci = zero
            do ia=1,natgas-1
               do ja=ia+1,natgas
                  if (ia.ne.nhgas.and.ja.ne.nhgas) then
                     dxc = distance(ia,ja,xyzgas)
                     vcci = vcci + chrgas(i,ia)*chrgas(i,ja)/dxc
                  endif
               enddo
            enddo
            vcc(i) = vcci*e2*ev2cal
         enddo

      endif

      ! Repulsion terms

      xoh = dabs(dp2 - q)
      xnh = dabs(dp2 + q)
      vra = droh*dexp(-broh*xoh) + croh/xoh**9
      vrb = drnh*dexp(-brnh*xnh) + crnh/xnh**9

      h0_(1,1) = vma + vc1a + vra + vcc(1) + corr(1)
      h0_(2,2) = vmb + vc1b + vrb + vcc(2) + corr(2)
      h0_(3,3) = vma + vc2a + vra + vcc(3) + corr(3)
      h0_(4,4) = vmb + vc2b + vrb + vcc(4) + corr(4)

      !==== OFF-DIAGONAL ELEMENTS

      ! Proton coupling 1a-1b/2a-2b

      vpt1 = vpt0_1*dexp(-gapt_1*(q - qpt0_1)**2)
      vpt2 = vpt0_2*dexp(-gapt_2*(q - qpt0_2)**2)

      ! Electron coupling 1a-2a/1b-2b

      veta = vet0_a*dexp(-gaet_a*(q - qet0_a)**2)
      vetb = vet0_b*dexp(-gaet_b*(q - qet0_b)**2)

      ! Mixed coupling 1a-2b/1b-2a

      vetpt1 = vetpt0_1*dexp(-gaetpt_1*(q - qetpt0_1)**2)
      vetpt2 = vetpt0_2*dexp(-gaetpt_2*(q - qetpt0_2)**2)

      h0_(1,2) = vpt1
      h0_(2,1) = vpt1
      h0_(1,3) = veta
      h0_(3,1) = veta
      h0_(1,4) = vetpt1
      h0_(4,1) = vetpt1
      h0_(2,3) = vetpt2
      h0_(3,2) = vetpt2
      h0_(2,4) = vetb
      h0_(4,2) = vetb
      h0_(3,4) = vpt2
      h0_(4,3) = vpt2

      return

   end subroutine h0mat_mm5


   subroutine dh0mat_mm5(dh0_)
   !======================================================================!
   ! Computes the derivative of gas phase Hamiltonian (au/Bohr)
   ! at proton position Q
   !======================================================================!
      real(8), intent(out), dimension(4,4) :: dh0_

      ! local variables
      integer :: pdgas, nhgas, pagas, i, j
      real(8) :: dp, dp2, q, dvma, dvmb, dvc1a, dvc1b, dvc2a, dvc2b
      real(8) :: xxi, dxc, dxc3, xoh, xnh, dvra, dvrb
      real(8) :: vpt1, vpt2, dvpt1, dvpt2, veta, vetb, dveta, dvetb
      real(8) :: vetpt1, vetpt2, dvetpt1, dvetpt2

      ! local variables: parameters
      real(8) :: dnh
      real(8) :: doh
      real(8) :: bnh
      real(8) :: boh
      real(8) :: rnh
      real(8) :: roh
      real(8) :: drnh
      real(8) :: droh
      real(8) :: brnh
      real(8) :: broh
      real(8) :: croh
      real(8) :: crnh
      real(8) :: ksi
      real(8) :: vpt0_1
      real(8) :: gapt_1
      real(8) :: qpt0_1
      real(8) :: vpt0_2
      real(8) :: gapt_2
      real(8) :: qpt0_2
      real(8) :: vet0_a
      real(8) :: gaet_a
      real(8) :: qet0_a
      real(8) :: vet0_b
      real(8) :: gaet_b
      real(8) :: qet0_b
      real(8) :: vetpt0_1
      real(8) :: gaetpt_1
      real(8) :: qetpt0_1
      real(8) :: vetpt0_2
      real(8) :: gaetpt_2
      real(8) :: qetpt0_2
      real(8), dimension(4) :: corr

      ! pointer to the real set of parameters
      type(mm5_1D), pointer :: par

      par => mm5par

      dnh      = par%dnh
      doh      = par%doh
      bnh      = par%bnh
      boh      = par%boh
      rnh      = par%rnh
      roh      = par%roh
      drnh     = par%drnh
      droh     = par%droh
      brnh     = par%brnh
      broh     = par%broh
      croh     = par%croh
      crnh     = par%crnh
      ksi      = par%ksi
      vpt0_1   = par%vpt0_1
      gapt_1   = par%gapt_1
      qpt0_1   = par%qpt0_1
      vpt0_2   = par%vpt0_2
      gapt_2   = par%gapt_2
      qpt0_2   = par%qpt0_2
      vet0_a   = par%vet0_a
      gaet_a   = par%gaet_a
      qet0_a   = par%qet0_a
      vet0_b   = par%vet0_b
      gaet_b   = par%gaet_b
      qet0_b   = par%qet0_b
      vetpt0_1 = par%vetpt0_1
      gaetpt_1 = par%gaetpt_1
      qetpt0_1 = par%qetpt0_1
      vetpt0_2 = par%vetpt0_2
      gaetpt_2 = par%gaetpt_2
      qetpt0_2 = par%qetpt0_2
      corr     = par%corr

      dh0_ = zero
      pdgas = iptgas(1)
      nhgas = iptgas(2)
      pagas = iptgas(3)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! DP - a distance between PT donor (PDGAS) and acceptor (PAGAS)
      ! Note that PT interface is oriented along the x-axis.
      ! Q - antisymmetrical proton coordinate
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      dp  = dabs(xyzgas(1,pagas)-xyzgas(1,pdgas))
      dp2 = dp/2.d0
      q   = xyzgas(1,nhgas)

      !===> DERIVATIVES OF THE DIAGONAL ELEMENTS

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Derivatives of the Morse potentials
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      dvma =  2.d0*dnh*bnh*(1.d0-dexp(-bnh*(q+dp2-rnh)))*dexp(-bnh*(q+dp2-rnh))
      dvmb = -2.d0*doh*boh*(1.d0-dexp(-boh*(dp2-q-roh)))*dexp(-boh*(dp2-q-roh))

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Derivatives of the Coulomb interaction terms
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      dvc1a = zero
      dvc1b = zero
      dvc2a = zero
      dvc2b = zero

      if (coulomb) then

         do i=1,natgas
            if (i.ne.pdgas.and.i.ne.nhgas) then
               xxi = q - xyzgas(1,i)
               dxc = distance(nhgas,i,xyzgas)
               dxc3 = dxc**(3.d0)
               dvc1a = dvc1a - chrgas(1,nhgas)*chrgas(1,i)*xxi/dxc3
            endif
         enddo
         dvc1a = dvc1a*e2*ev2cal

         do i=1,natgas
            if (i.ne.pagas.and.i.ne.nhgas) then
               xxi = q - xyzgas(1,i)
               dxc = distance(nhgas,i,xyzgas)
               dxc3 = dxc**(3.d0)
               dvc1b = dvc1b - chrgas(2,nhgas)*chrgas(2,i)*xxi/dxc3
            endif
         enddo
         dvc1b = dvc1b*e2*ev2cal

         do i=1,natgas
            if (i.ne.pdgas.and.i.ne.nhgas) then
               xxi = q - xyzgas(1,i)
               dxc = distance(nhgas,i,xyzgas)
               dxc3 = dxc**(3.d0)
               dvc2a = dvc2a - chrgas(3,nhgas)*chrgas(3,i)*xxi/dxc3
            endif
         enddo
         dvc2a = dvc2a*e2*ev2cal

         do i=1,natgas
            if (i.ne.pagas.and.i.ne.nhgas) then
               xxi = q - xyzgas(1,i)
               dxc = distance(nhgas,i,xyzgas)
               dxc3 = dxc**(3.d0)
               dvc2b = dvc2b - chrgas(4,nhgas)*chrgas(4,i)*xxi/dxc3
            endif
         enddo
         dvc2b = dvc2b*e2*ev2cal

      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Derivatives of the repulsion terms (kcal/mole)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      xoh = dabs(dp2 - q)
      xnh = dabs(dp2 + q)
      dvra =  broh*droh*dexp(-broh*xoh) + 9.d0*croh/xoh**10
      dvrb = -brnh*drnh*dexp(-brnh*xnh) - 9.d0*crnh/xnh**10

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Diagonal elements (kcal/mole)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      dh0_(1,1) = dvma + dvc1a + dvra
      dh0_(2,2) = dvmb + dvc1b + dvrb
      dh0_(3,3) = dvma + dvc2a + dvra
      dh0_(4,4) = dvmb + dvc2b + dvrb

      !===> OFF-DIAGONAL ELEMENTS (kcal/mole)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Proton coupling 1a-1b/2a-2b
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      vpt1 = vpt0_1*dexp(-gapt_1*(q - qpt0_1)**2)
      vpt2 = vpt0_2*dexp(-gapt_2*(q - qpt0_2)**2)
      dvpt1 = -2.d0*gapt_1*(q - qpt0_1)*vpt1
      dvpt2 = -2.d0*gapt_2*(q - qpt0_2)*vpt2

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Electron coupling 1a-2a/1b-2b
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      veta = vet0_a*dexp(-gaet_a*(q - qet0_a)**2)
      vetb = vet0_b*dexp(-gaet_b*(q - qet0_b)**2)
      dveta = -2.d0*gaet_a*(q - qet0_a)*veta
      dvetb = -2.d0*gaet_b*(q - qet0_b)*vetb

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Mixed coupling 1a-2b/1b-2a
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      vetpt1 = vetpt0_1*dexp(-gaetpt_1*(q - qetpt0_1)**2)
      vetpt2 = vetpt0_2*dexp(-gaetpt_2*(q - qetpt0_2)**2)
      dvetpt1 = -2.d0*gaetpt_1*(q - qetpt0_1)*vetpt1
      dvetpt2 = -2.d0*gaetpt_2*(q - qetpt0_2)*vetpt2

      dh0_(1,2) = dvpt1
      dh0_(2,1) = dvpt1
      dh0_(1,3) = dveta
      dh0_(3,1) = dveta
      dh0_(1,4) = dvetpt1
      dh0_(4,1) = dvetpt1
      dh0_(2,3) = dvetpt2
      dh0_(3,2) = dvetpt2
      dh0_(2,4) = dvetb
      dh0_(4,2) = dvetb
      dh0_(3,4) = dvpt2
      dh0_(4,3) = dvpt2

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Convert units: kcal/mole/Angstroem --> au/Bohr
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do i=1,4
         do j=1,4
            dh0_(i,j) = dh0_(i,j)*cal2au/a2bohr
         enddo
      enddo

      return

   end subroutine dh0mat_mm5

   subroutine d2h0mat_mm5(d2h0_)
   !======================================================================C
   ! Computes the second derivative of the gas phase Hamiltonian
   ! (au/A**2) at proton position Q
   !======================================================================C
      implicit real(8) (a-h,o-z)
      implicit integer (i-n)

      real(8), intent(out), dimension(4,4) :: d2h0_

      integer :: pdgas,nhgas,pagas

      ! local variables: parameters
      real(8) :: dnh
      real(8) :: doh
      real(8) :: bnh
      real(8) :: boh
      real(8) :: rnh
      real(8) :: roh
      real(8) :: drnh
      real(8) :: droh
      real(8) :: brnh
      real(8) :: broh
      real(8) :: croh
      real(8) :: crnh
      real(8) :: ksi
      real(8) :: vpt0_1
      real(8) :: gapt_1
      real(8) :: qpt0_1
      real(8) :: vpt0_2
      real(8) :: gapt_2
      real(8) :: qpt0_2
      real(8) :: vet0_a
      real(8) :: gaet_a
      real(8) :: qet0_a
      real(8) :: vet0_b
      real(8) :: gaet_b
      real(8) :: qet0_b
      real(8) :: vetpt0_1
      real(8) :: gaetpt_1
      real(8) :: qetpt0_1
      real(8) :: vetpt0_2
      real(8) :: gaetpt_2
      real(8) :: qetpt0_2
      real(8), dimension(4) :: corr

      ! pointer to the real set of parameters
      type(mm5_1d), pointer :: par

      par => mm5par

      dnh      = par%dnh
      doh      = par%doh
      bnh      = par%bnh
      boh      = par%boh
      rnh      = par%rnh
      roh      = par%roh
      drnh     = par%drnh
      droh     = par%droh
      brnh     = par%brnh
      broh     = par%broh
      croh     = par%croh
      crnh     = par%crnh
      ksi      = par%ksi
      vpt0_1   = par%vpt0_1
      gapt_1   = par%gapt_1
      qpt0_1   = par%qpt0_1
      vpt0_2   = par%vpt0_2
      gapt_2   = par%gapt_2
      qpt0_2   = par%qpt0_2
      vet0_a   = par%vet0_a
      gaet_a   = par%gaet_a
      qet0_a   = par%qet0_a
      vet0_b   = par%vet0_b
      gaet_b   = par%gaet_b
      qet0_b   = par%qet0_b
      vetpt0_1 = par%vetpt0_1
      gaetpt_1 = par%gaetpt_1
      qetpt0_1 = par%qetpt0_1
      vetpt0_2 = par%vetpt0_2
      gaetpt_2 = par%gaetpt_2
      qetpt0_2 = par%qetpt0_2
      corr     = par%corr

      d2h0_ = zero
      pdgas = iptgas(1)
      nhgas = iptgas(2)
      pagas = iptgas(3)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! DP - a distance between PT donor (PDGAS) and acceptor (PAGAS)
      ! Note that PT interface is oriented along the x-axis.
      ! Q - antisymmetrical proton coordinate
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      dp  = dabs(xyzgas(1,pagas)-xyzgas(1,pdgas))
      dp2 = dp/2.d0
      q   = xyzgas(1,nhgas)

      !===> DERIVATIVES OF THE DIAGONAL ELEMENTS

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Derivatives of the Morse potentials
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      d2vma = -2.d0*dnh*bnh*bnh*(-2.d0+dexp(bnh*(q+dp2-rnh)))*dexp(-2.d0*bnh*(q+dp2-rnh))
      d2vmb = -2.d0*doh*boh*boh*(-2.d0+dexp(boh*(dp2-q-roh)))*dexp(-2.d0*boh*(dp2-q-roh))

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Derivatives of the Coulomb interaction terms
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      d2vc1a = zero
      d2vc1b = zero
      d2vc2a = zero
      d2vc2b = zero

      if (coulomb) then

         do i=1,natgas
            if (i.ne.pdgas.and.i.ne.nhgas) then
               xxi = q - xyzgas(1,i)
               xxi2 = xxi*xxi
               dxc = distance(nhgas,i,xyzgas)
               dxc2 = dxc*dxc
               dxc3 = dxc2*dxc
               dxc5 = dxc2*dxc3
               d2vc1a = d2vc1a - chrgas(1,nhgas)*chrgas(1,i)*(dxc2 - 3.d0*xxi2)/dxc5
            endif
         enddo
         d2vc1a = d2vc1a*e2*ev2cal

         do i=1,natgas
            if (i.ne.pagas.and.i.ne.nhgas) then
               xxi = q - xyzgas(1,i)
               xxi2 = xxi*xxi
               dxc = distance(nhgas,i,xyzgas)
               dxc2 = dxc*dxc
               dxc3 = dxc2*dxc
               dxc5 = dxc2*dxc3
               d2vc1b = d2vc1b - chrgas(2,nhgas)*chrgas(2,i)*(dxc2 - 3.d0*xxi2)/dxc5
            endif
         enddo
         d2vc1b = d2vc1b*e2*ev2cal

         do i=1,natgas
            if (i.ne.pdgas.and.i.ne.nhgas) then
               xxi = q - xyzgas(1,i)
               xxi2 = xxi*xxi
               dxc = distance(nhgas,i,xyzgas)
               dxc2 = dxc*dxc
               dxc3 = dxc2*dxc
               dxc5 = dxc2*dxc3
               d2vc2a = d2vc2a - chrgas(3,nhgas)*chrgas(3,i)*(dxc2 - 3.d0*xxi2)/dxc5
            endif
         enddo
         d2vc2a = d2vc2a*e2*ev2cal

         do i=1,natgas
            if (i.ne.pagas.and.i.ne.nhgas) then
               xxi = q - xyzgas(1,i)
               xxi2 = xxi*xxi
               dxc = distance(nhgas,i,xyzgas)
               dxc2 = dxc*dxc
               dxc3 = dxc2*dxc
               dxc5 = dxc2*dxc3
               d2vc2b = d2vc2b - chrgas(4,nhgas)*chrgas(4,i)*(dxc2 - 3.d0*xxi2)/dxc5
            endif
         enddo
         d2vc2b = d2vc2b*e2*ev2cal

      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Second derivatives of the repulsion terms
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      xoh = dabs(dp2 - q)
      xnh = dabs(dp2 + q)
      d2vra =  broh*broh*droh*dexp(-broh*xoh) + 90.d0*croh/xoh**11.d0
      d2vrb =  brnh*brnh*drnh*dexp(-brnh*xnh) + 90.d0*crnh/xnh**11.d0

      !====> Diagonal elements

      d2h0_(1,1) = d2vma + d2vc1a + d2vra
      d2h0_(2,2) = d2vmb + d2vc1b + d2vrb
      d2h0_(3,3) = d2vma + d2vc2a + d2vra
      d2h0_(4,4) = d2vmb + d2vc2b + d2vrb

      !====> Off-diagonal elements

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Proton coupling 1a-1b/2a-2b
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      vpt1 = vpt0_1*dexp(-gapt_1*(q - qpt0_1)*(q - qpt0_1))
      vpt2 = vpt0_2*dexp(-gapt_2*(q - qpt0_2)*(q - qpt0_2))
      d2vpt1 = 2.d0*gapt_1*vpt1*(-1.d0+2.d0*gapt_1*(q - qpt0_1)*(q - qpt0_1))
      d2vpt2 = 2.d0*gapt_2*vpt2*(-1.d0+2.d0*gapt_2*(q - qpt0_2)*(q - qpt0_2))

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Electron coupling 1a-2a/1b-2b
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      veta = vet0_a*dexp(-gaet_a*(q - qet0_a)*(q - qet0_a))
      vetb = vet0_b*dexp(-gaet_b*(q - qet0_b)*(q - qet0_b))
      d2veta = 2.d0*gaet_a*veta*(-1.d0+2.d0*gaet_a*(q - qet0_a)*(q - qet0_a))
      d2vetb = 2.d0*gaet_b*vetb*(-1.d0+2.d0*gaet_b*(q - qet0_b)*(q - qet0_b))

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Mixed coupling 1a-2b/1b-2a
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      vetpt1 = vetpt0_1*dexp(-gaetpt_1*(q - qetpt0_1)*(q - qetpt0_1))
      vetpt2 = vetpt0_2*dexp(-gaetpt_2*(q - qetpt0_2)*(q - qetpt0_2))
      d2vetpt1 = 2.d0*gaetpt_1*vetpt1*(-1.d0+2.d0*gaetpt_1*(q - qetpt0_1)*(q - qetpt0_1))
      d2vetpt2 = 2.d0*gaetpt_2*vetpt2*(-1.d0+2.d0*gaetpt_2*(q - qetpt0_2)*(q - qetpt0_2))

      d2h0_(1,2) = d2vpt1
      d2h0_(2,1) = d2vpt1
      d2h0_(1,3) = d2veta
      d2h0_(3,1) = d2veta
      d2h0_(1,4) = d2vetpt1
      d2h0_(4,1) = d2vetpt1
      d2h0_(2,3) = d2vetpt2
      d2h0_(3,2) = d2vetpt2
      d2h0_(2,4) = d2vetb
      d2h0_(4,2) = d2vetb
      d2h0_(3,4) = d2vpt2
      d2h0_(4,3) = d2vpt2

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Convert units: kcal/mole/Angstroem**2 --> au/Bohr**2
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do i=1,4
         do j=1,4
            d2h0_(i,j) = d2h0_(i,j)*cal2au/(a2bohr*a2bohr)
         enddo
      enddo

      return

   end subroutine d2h0mat_mm5

   subroutine h0mat_mm5gen(h0_)
   !======================================================================!
   ! Calculates the gas phase Hamiltonian matrix using
   ! MMGEN potential for general O-H...O systems (kcal/mole)
   !======================================================================!
      implicit none

      real(8), intent(out), dimension(4,4) :: h0_

      integer :: pdgas, nhgas, pagas, i, ia, ja
      real(8)  :: dp, dp2, q, vm1a, vm1b, vm2a, vm2b, vc1a, vc1b, vc2a, vc2b, dxc, vcci
      real(8)  :: xoha, xohb, vr1a, vr1b, vr2a, vr2b, vpt1, vpt2, veta, vetb, vetpt1, vetpt2
      real(8), dimension(4) :: vcc

      ! local variables: parameters
      real(8) :: doh1a
      real(8) :: doh1b
      real(8) :: doh2a
      real(8) :: doh2b
      real(8) :: boh1a
      real(8) :: boh1b
      real(8) :: boh2a
      real(8) :: boh2b
      real(8) :: roh1a
      real(8) :: roh1b
      real(8) :: roh2a
      real(8) :: roh2b
      real(8) :: droh1a
      real(8) :: droh1b
      real(8) :: droh2a
      real(8) :: droh2b
      real(8) :: broh1a
      real(8) :: broh1b
      real(8) :: broh2a
      real(8) :: broh2b
      real(8) :: croh1a
      real(8) :: croh1b
      real(8) :: croh2a
      real(8) :: croh2b
      real(8) :: ksi
      real(8) :: vpt0_1
      real(8) :: gapt_1
      real(8) :: qpt0_1
      real(8) :: vpt0_2
      real(8) :: gapt_2
      real(8) :: qpt0_2
      real(8) :: vet0_a
      real(8) :: gaet_a
      real(8) :: qet0_a
      real(8) :: vet0_b
      real(8) :: gaet_b
      real(8) :: qet0_b
      real(8) :: vetpt0_1
      real(8) :: gaetpt_1
      real(8) :: qetpt0_1
      real(8) :: vetpt0_2
      real(8) :: gaetpt_2
      real(8) :: qetpt0_2
      real(8), dimension(4) :: corr

      ! pointer to the real set of parameters
      type(mm5gen_1d), pointer :: par

      par => mm5genpar

      doh1a    = par%doh1a
      doh1b    = par%doh1b
      doh2a    = par%doh2a
      doh2b    = par%doh2b
      boh1a    = par%boh1a
      boh1b    = par%boh1b
      boh2a    = par%boh2a
      boh2b    = par%boh2b
      roh1a    = par%roh1a
      roh1b    = par%roh1b
      roh2a    = par%roh2a
      roh2b    = par%roh2b
      droh1a   = par%droh1a
      droh1b   = par%droh1b
      droh2a   = par%droh2a
      droh2b   = par%droh2b
      broh1a   = par%broh1a
      broh1b   = par%broh1b
      broh2a   = par%broh2a
      broh2b   = par%broh2b
      croh1a   = par%croh1a
      croh1b   = par%croh1b
      croh2a   = par%croh2a
      croh2b   = par%croh2b
      ksi      = par%ksi
      vpt0_1   = par%vpt0_1
      gapt_1   = par%gapt_1
      qpt0_1   = par%qpt0_1
      vpt0_2   = par%vpt0_2
      gapt_2   = par%gapt_2
      qpt0_2   = par%qpt0_2
      vet0_a   = par%vet0_a
      gaet_a   = par%gaet_a
      qet0_a   = par%qet0_a
      vet0_b   = par%vet0_b
      gaet_b   = par%gaet_b
      qet0_b   = par%qet0_b
      vetpt0_1 = par%vetpt0_1
      gaetpt_1 = par%gaetpt_1
      qetpt0_1 = par%qetpt0_1
      vetpt0_2 = par%vetpt0_2
      gaetpt_2 = par%gaetpt_2
      qetpt0_2 = par%qetpt0_2
      corr     = par%corr

      h0_ = zero
      pdgas = iptgas(1)
      nhgas = iptgas(2)
      pagas = iptgas(3)

      ! DP - a distance between PT donor (PDGAS) and acceptor (PAGAS)
      ! Note that PT interface is oriented along the x-axis.
      ! Q - antisymmetrical proton coordinate

      dp  = dabs(xyzgas(1,pagas)-xyzgas(1,pdgas))
      dp2 = dp/2.d0
      q   = xyzgas(1,nhgas)

      !==== DIAGONAL ELEMENTS

      ! Morse potentials for proton transfer EVB states

      vm1a = doh1a*(1.d0-dexp(-boh1a*(q+dp2-roh1a)))**2.d0
      vm1b = doh1b*(1.d0-dexp(-boh1b*(dp2-q-roh1b)))**2.d0
      vm2a = doh2a*(1.d0-dexp(-boh2a*(q+dp2-roh2a)))**2.d0
      vm2b = doh2b*(1.d0-dexp(-boh2b*(dp2-q-roh2b)))**2.d0

      ! Coulomb interaction terms

      vc1a = zero
      vc1b = zero
      vc2a = zero
      vc2b = zero
      do i=1,4
         vcc(i) = zero
      enddo

      if (coulomb) then

         do i=1,natgas
            if (i.ne.pdgas.and.i.ne.nhgas) then
               dxc = distance(nhgas,i,xyzgas)
               vc1a = vc1a + chrgas(1,nhgas)*chrgas(1,i)/dxc
            endif
         enddo
         vc1a = vc1a*e2*ev2cal

         do i=1,natgas
            if (i.ne.pagas.and.i.ne.nhgas) then
               dxc = distance(nhgas,i,xyzgas)
               vc1b = vc1b + chrgas(2,nhgas)*chrgas(2,i)/dxc
            endif
         enddo
         vc1b = vc1b*e2*ev2cal

         do i=1,natgas
            if (i.ne.pdgas.and.i.ne.nhgas) then
               dxc = distance(nhgas,i,xyzgas)
               vc2a = vc2a + chrgas(3,nhgas)*chrgas(3,i)/dxc
            endif
         enddo
         vc2a = vc2a*e2*ev2cal

         do i=1,5
            if (i.ne.pagas.and.i.ne.nhgas) then
               dxc = distance(nhgas,i,xyzgas)
               vc2b = vc2b + chrgas(4,nhgas)*chrgas(4,i)/dxc
            endif
         enddo
         vc2b = vc2b*e2*ev2cal

         ! Constant Coulomb interaction terms

         do i=1,nelst
            vcci = zero
            do ia=1,natgas-1
               do ja=ia+1,natgas
                  if (ia.ne.nhgas.and.ja.ne.nhgas) then
                     dxc = distance(ia,ja,xyzgas)
                     vcci = vcci + chrgas(i,ia)*chrgas(i,ja)/dxc
                  endif
               enddo
            enddo
            vcc(i) = vcci*e2*ev2cal
         enddo

      endif

      ! Repulsion terms

      xoha = dabs(dp2 - q)
      xohb = dabs(dp2 + q)
      vr1a = droh1a*dexp(-broh1a*xoha) + croh1a/xoha**9
      vr1b = droh1b*dexp(-broh1b*xohb) + croh1b/xohb**9
      vr2a = droh2a*dexp(-broh2a*xoha) + croh2a/xoha**9
      vr2b = droh2b*dexp(-broh2b*xohb) + croh2b/xohb**9

      h0_(1,1) = vm1a + vc1a + vr1a + vcc(1) + corr(1)
      h0_(2,2) = vm1b + vc1b + vr1b + vcc(2) + corr(2)
      h0_(3,3) = vm2a + vc2a + vr2a + vcc(3) + corr(3)
      h0_(4,4) = vm2b + vc2b + vr2b + vcc(4) + corr(4)

      !==== OFF-DIAGONAL ELEMENTS

      ! Proton coupling 1a-1b/2a-2b

      vpt1 = vpt0_1*dexp(-gapt_1*(q - qpt0_1)**2)
      vpt2 = vpt0_2*dexp(-gapt_2*(q - qpt0_2)**2)

      ! Electron coupling 1a-2a/1b-2b

      veta = vet0_a*dexp(-gaet_a*(q - qet0_a)**2)
      vetb = vet0_b*dexp(-gaet_b*(q - qet0_b)**2)

      ! Mixed coupling 1a-2b/1b-2a

      vetpt1 = vetpt0_1*dexp(-gaetpt_1*(q - qetpt0_1)**2)
      vetpt2 = vetpt0_2*dexp(-gaetpt_2*(q - qetpt0_2)**2)

      h0_(1,2) = vpt1
      h0_(2,1) = vpt1
      h0_(1,3) = veta
      h0_(3,1) = veta
      h0_(1,4) = vetpt1
      h0_(4,1) = vetpt1
      h0_(2,3) = vetpt2
      h0_(3,2) = vetpt2
      h0_(2,4) = vetb
      h0_(4,2) = vetb
      h0_(3,4) = vpt2
      h0_(4,3) = vpt2

      return

   end subroutine h0mat_mm5gen

   subroutine dh0mat_mm5gen(dh0_)
   !======================================================================!
   ! Computes the derivative of gas phase Hamiltonian (au/Bohr)
   ! at proton position Q
   !======================================================================!
      real(8), intent(out), dimension(4,4) :: dh0_

      ! local variables
      integer :: pdgas, nhgas, pagas, i, j
      real(8) :: dp, dp2, q, dvm1a, dvm1b, dvm2a, dvm2b, dvc1a, dvc1b, dvc2a, dvc2b
      real(8) :: xxi, dxc, dxc3, xoha, xohb, dvr1a, dvr1b, dvr2a, dvr2b
      real(8) :: vpt1, vpt2, dvpt1, dvpt2, veta, vetb, dveta, dvetb
      real(8) :: vetpt1, vetpt2, dvetpt1, dvetpt2

      ! local variables: parameters
      real(8) :: doh1a
      real(8) :: doh1b
      real(8) :: doh2a
      real(8) :: doh2b
      real(8) :: boh1a
      real(8) :: boh1b
      real(8) :: boh2a
      real(8) :: boh2b
      real(8) :: roh1a
      real(8) :: roh1b
      real(8) :: roh2a
      real(8) :: roh2b
      real(8) :: droh1a
      real(8) :: droh1b
      real(8) :: droh2a
      real(8) :: droh2b
      real(8) :: broh1a
      real(8) :: broh1b
      real(8) :: broh2a
      real(8) :: broh2b
      real(8) :: croh1a
      real(8) :: croh1b
      real(8) :: croh2a
      real(8) :: croh2b
      real(8) :: ksi
      real(8) :: vpt0_1
      real(8) :: gapt_1
      real(8) :: qpt0_1
      real(8) :: vpt0_2
      real(8) :: gapt_2
      real(8) :: qpt0_2
      real(8) :: vet0_a
      real(8) :: gaet_a
      real(8) :: qet0_a
      real(8) :: vet0_b
      real(8) :: gaet_b
      real(8) :: qet0_b
      real(8) :: vetpt0_1
      real(8) :: gaetpt_1
      real(8) :: qetpt0_1
      real(8) :: vetpt0_2
      real(8) :: gaetpt_2
      real(8) :: qetpt0_2
      real(8), dimension(4) :: corr

      ! pointer to the real set of parameters
      type(mm5gen_1d), pointer :: par

      par => mm5genpar

      doh1a    = par%doh1a
      doh1b    = par%doh1b
      doh2a    = par%doh2a
      doh2b    = par%doh2b
      boh1a    = par%boh1a
      boh1b    = par%boh1b
      boh2a    = par%boh2a
      boh2b    = par%boh2b
      roh1a    = par%roh1a
      roh1b    = par%roh1b
      roh2a    = par%roh2a
      roh2b    = par%roh2b
      droh1a   = par%droh1a
      droh1b   = par%droh1b
      droh2a   = par%droh2a
      droh2b   = par%droh2b
      broh1a   = par%broh1a
      broh1b   = par%broh1b
      broh2a   = par%broh2a
      broh2b   = par%broh2b
      croh1a   = par%croh1a
      croh1b   = par%croh1b
      croh2a   = par%croh2a
      croh2b   = par%croh2b
      ksi      = par%ksi
      vpt0_1   = par%vpt0_1
      gapt_1   = par%gapt_1
      qpt0_1   = par%qpt0_1
      vpt0_2   = par%vpt0_2
      gapt_2   = par%gapt_2
      qpt0_2   = par%qpt0_2
      vet0_a   = par%vet0_a
      gaet_a   = par%gaet_a
      qet0_a   = par%qet0_a
      vet0_b   = par%vet0_b
      gaet_b   = par%gaet_b
      qet0_b   = par%qet0_b
      vetpt0_1 = par%vetpt0_1
      gaetpt_1 = par%gaetpt_1
      qetpt0_1 = par%qetpt0_1
      vetpt0_2 = par%vetpt0_2
      gaetpt_2 = par%gaetpt_2
      qetpt0_2 = par%qetpt0_2
      corr     = par%corr

      dh0_ = zero
      pdgas = iptgas(1)
      nhgas = iptgas(2)
      pagas = iptgas(3)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! DP - a distance between PT donor (PDGAS) and acceptor (PAGAS)
      ! Note that PT interface is oriented along the x-axis.
      ! Q - antisymmetrical proton coordinate
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      dp  = dabs(xyzgas(1,pagas)-xyzgas(1,pdgas))
      dp2 = dp/2.d0
      q   = xyzgas(1,nhgas)

      !===> DERIVATIVES OF THE DIAGONAL ELEMENTS

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Derivatives of the Morse potentials
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      dvm1a =  2.d0*doh1a*boh1a*(1.d0-dexp(-boh1a*(q+dp2-roh1a)))*dexp(-boh1a*(q+dp2-roh1a))
      dvm1b = -2.d0*doh1b*boh1b*(1.d0-dexp(-boh1b*(dp2-q-roh1b)))*dexp(-boh1b*(dp2-q-roh1b))
      dvm1a =  2.d0*doh2a*boh2a*(1.d0-dexp(-boh2a*(q+dp2-roh2a)))*dexp(-boh2a*(q+dp2-roh2a))
      dvm1b = -2.d0*doh2b*boh2b*(1.d0-dexp(-boh2b*(dp2-q-roh2b)))*dexp(-boh2b*(dp2-q-roh2b))

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Derivatives of the Coulomb interaction terms
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      dvc1a = zero
      dvc1b = zero
      dvc2a = zero
      dvc2b = zero

      if (coulomb) then

         do i=1,natgas
            if (i.ne.pdgas.and.i.ne.nhgas) then
               xxi = q - xyzgas(1,i)
               dxc = distance(nhgas,i,xyzgas)
               dxc3 = dxc**(3.d0)
               dvc1a = dvc1a - chrgas(1,nhgas)*chrgas(1,i)*xxi/dxc3
            endif
         enddo
         dvc1a = dvc1a*e2*ev2cal

         do i=1,natgas
            if (i.ne.pagas.and.i.ne.nhgas) then
               xxi = q - xyzgas(1,i)
               dxc = distance(nhgas,i,xyzgas)
               dxc3 = dxc**(3.d0)
               dvc1b = dvc1b - chrgas(2,nhgas)*chrgas(2,i)*xxi/dxc3
            endif
         enddo
         dvc1b = dvc1b*e2*ev2cal

         do i=1,natgas
            if (i.ne.pdgas.and.i.ne.nhgas) then
               xxi = q - xyzgas(1,i)
               dxc = distance(nhgas,i,xyzgas)
               dxc3 = dxc**(3.d0)
               dvc2a = dvc2a - chrgas(3,nhgas)*chrgas(3,i)*xxi/dxc3
            endif
         enddo
         dvc2a = dvc2a*e2*ev2cal

         do i=1,natgas
            if (i.ne.pagas.and.i.ne.nhgas) then
               xxi = q - xyzgas(1,i)
               dxc = distance(nhgas,i,xyzgas)
               dxc3 = dxc**(3.d0)
               dvc2b = dvc2b - chrgas(4,nhgas)*chrgas(4,i)*xxi/dxc3
            endif
         enddo
         dvc2b = dvc2b*e2*ev2cal

      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Derivatives of the repulsion terms (kcal/mole)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      xoha = dabs(dp2 - q)
      xohb = dabs(dp2 + q)
      dvr1a =  broh1a*droh1a*dexp(-broh1a*xoha) + 9.d0*croh1a/xoha**10
      dvr1b = -broh1b*droh1b*dexp(-broh1b*xohb) - 9.d0*croh1b/xohb**10
      dvr2a =  broh2a*droh2a*dexp(-broh2a*xoha) + 9.d0*croh2a/xoha**10
      dvr2b = -broh2b*droh2b*dexp(-broh2b*xohb) - 9.d0*croh2b/xohb**10

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Diagonal elements (kcal/mole)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      dh0_(1,1) = dvm1a + dvc1a + dvr1a
      dh0_(2,2) = dvm1b + dvc1b + dvr1b
      dh0_(3,3) = dvm2a + dvc2a + dvr2a
      dh0_(4,4) = dvm2b + dvc2b + dvr2b

      !===> OFF-DIAGONAL ELEMENTS (kcal/mole)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Proton coupling 1a-1b/2a-2b
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      vpt1 = vpt0_1*dexp(-gapt_1*(q - qpt0_1)**2)
      vpt2 = vpt0_2*dexp(-gapt_2*(q - qpt0_2)**2)
      dvpt1 = -2.d0*gapt_1*(q - qpt0_1)*vpt1
      dvpt2 = -2.d0*gapt_2*(q - qpt0_2)*vpt2

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Electron coupling 1a-2a/1b-2b
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      veta = vet0_a*dexp(-gaet_a*(q - qet0_a)**2)
      vetb = vet0_b*dexp(-gaet_b*(q - qet0_b)**2)
      dveta = -2.d0*gaet_a*(q - qet0_a)*veta
      dvetb = -2.d0*gaet_b*(q - qet0_b)*vetb

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Mixed coupling 1a-2b/1b-2a
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      vetpt1 = vetpt0_1*dexp(-gaetpt_1*(q - qetpt0_1)**2)
      vetpt2 = vetpt0_2*dexp(-gaetpt_2*(q - qetpt0_2)**2)
      dvetpt1 = -2.d0*gaetpt_1*(q - qetpt0_1)*vetpt1
      dvetpt2 = -2.d0*gaetpt_2*(q - qetpt0_2)*vetpt2

      dh0_(1,2) = dvpt1
      dh0_(2,1) = dvpt1
      dh0_(1,3) = dveta
      dh0_(3,1) = dveta
      dh0_(1,4) = dvetpt1
      dh0_(4,1) = dvetpt1
      dh0_(2,3) = dvetpt2
      dh0_(3,2) = dvetpt2
      dh0_(2,4) = dvetb
      dh0_(4,2) = dvetb
      dh0_(3,4) = dvpt2
      dh0_(4,3) = dvpt2

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Convert units: kcal/mole/Angstroem --> au/Bohr
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do i=1,4
         do j=1,4
            dh0_(i,j) = dh0_(i,j)*cal2au/a2bohr
         enddo
      enddo

      return

   end subroutine dh0mat_mm5gen

   subroutine d2h0mat_mm5gen(d2h0_)
   !======================================================================C
   ! Computes the second derivative of the gas phase Hamiltonian
   ! (au/A**2) at proton position Q
   !======================================================================C
      implicit real(8) (a-h,o-z)
      implicit integer (i-n)

      real(8), intent(out), dimension(4,4) :: d2h0_

      integer :: pdgas,nhgas,pagas

      ! local variables: parameters
      real(8) :: doh1a
      real(8) :: doh1b
      real(8) :: doh2a
      real(8) :: doh2b
      real(8) :: boh1a
      real(8) :: boh1b
      real(8) :: boh2a
      real(8) :: boh2b
      real(8) :: roh1a
      real(8) :: roh1b
      real(8) :: roh2a
      real(8) :: roh2b
      real(8) :: droh1a
      real(8) :: droh1b
      real(8) :: droh2a
      real(8) :: droh2b
      real(8) :: broh1a
      real(8) :: broh1b
      real(8) :: broh2a
      real(8) :: broh2b
      real(8) :: croh1a
      real(8) :: croh1b
      real(8) :: croh2a
      real(8) :: croh2b
      real(8) :: ksi
      real(8) :: vpt0_1
      real(8) :: gapt_1
      real(8) :: qpt0_1
      real(8) :: vpt0_2
      real(8) :: gapt_2
      real(8) :: qpt0_2
      real(8) :: vet0_a
      real(8) :: gaet_a
      real(8) :: qet0_a
      real(8) :: vet0_b
      real(8) :: gaet_b
      real(8) :: qet0_b
      real(8) :: vetpt0_1
      real(8) :: gaetpt_1
      real(8) :: qetpt0_1
      real(8) :: vetpt0_2
      real(8) :: gaetpt_2
      real(8) :: qetpt0_2
      real(8), dimension(4) :: corr

      ! pointer to the real set of parameters
      type(mm5gen_1d), pointer :: par

      par => mm5genpar

      doh1a    = par%doh1a
      doh1b    = par%doh1b
      doh2a    = par%doh2a
      doh2b    = par%doh2b
      boh1a    = par%boh1a
      boh1b    = par%boh1b
      boh2a    = par%boh2a
      boh2b    = par%boh2b
      roh1a    = par%roh1a
      roh1b    = par%roh1b
      roh2a    = par%roh2a
      roh2b    = par%roh2b
      droh1a   = par%droh1a
      droh1b   = par%droh1b
      droh2a   = par%droh2a
      droh2b   = par%droh2b
      broh1a   = par%broh1a
      broh1b   = par%broh1b
      broh2a   = par%broh2a
      broh2b   = par%broh2b
      croh1a   = par%croh1a
      croh1b   = par%croh1b
      croh2a   = par%croh2a
      croh2b   = par%croh2b
      ksi      = par%ksi
      vpt0_1   = par%vpt0_1
      gapt_1   = par%gapt_1
      qpt0_1   = par%qpt0_1
      vpt0_2   = par%vpt0_2
      gapt_2   = par%gapt_2
      qpt0_2   = par%qpt0_2
      vet0_a   = par%vet0_a
      gaet_a   = par%gaet_a
      qet0_a   = par%qet0_a
      vet0_b   = par%vet0_b
      gaet_b   = par%gaet_b
      qet0_b   = par%qet0_b
      vetpt0_1 = par%vetpt0_1
      gaetpt_1 = par%gaetpt_1
      qetpt0_1 = par%qetpt0_1
      vetpt0_2 = par%vetpt0_2
      gaetpt_2 = par%gaetpt_2
      qetpt0_2 = par%qetpt0_2
      corr     = par%corr

      d2h0_ = zero
      pdgas = iptgas(1)
      nhgas = iptgas(2)
      pagas = iptgas(3)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! DP - a distance between PT donor (PDGAS) and acceptor (PAGAS)
      ! Note that PT interface is oriented along the x-axis.
      ! Q - antisymmetrical proton coordinate
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      dp  = dabs(xyzgas(1,pagas)-xyzgas(1,pdgas))
      dp2 = dp/2.d0
      q   = xyzgas(1,nhgas)

      !===> DERIVATIVES OF THE DIAGONAL ELEMENTS

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Derivatives of the Morse potentials
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      d2vm1a = -2.d0*doh1a*boh1a*boh1a*(-2.d0+dexp(boh1a*(q+dp2-roh1a)))*dexp(-2.d0*boh1a*(q+dp2-roh1a))
      d2vm1b = -2.d0*doh1b*boh1b*boh1b*(-2.d0+dexp(boh1b*(dp2-q-roh1b)))*dexp(-2.d0*boh1b*(dp2-q-roh1b))
      d2vm2a = -2.d0*doh2a*boh2a*boh2a*(-2.d0+dexp(boh2a*(q+dp2-roh2a)))*dexp(-2.d0*boh2a*(q+dp2-roh2a))
      d2vm2b = -2.d0*doh2b*boh2b*boh2b*(-2.d0+dexp(boh2b*(dp2-q-roh2b)))*dexp(-2.d0*boh2b*(dp2-q-roh2b))

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Derivatives of the Coulomb interaction terms
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      d2vc1a = zero
      d2vc1b = zero
      d2vc2a = zero
      d2vc2b = zero

      if (coulomb) then

         do i=1,natgas
            if (i.ne.pdgas.and.i.ne.nhgas) then
               xxi = q - xyzgas(1,i)
               xxi2 = xxi*xxi
               dxc = distance(nhgas,i,xyzgas)
               dxc2 = dxc*dxc
               dxc3 = dxc2*dxc
               dxc5 = dxc2*dxc3
               d2vc1a = d2vc1a - chrgas(1,nhgas)*chrgas(1,i)*(dxc2 - 3.d0*xxi2)/dxc5
            endif
         enddo
         d2vc1a = d2vc1a*e2*ev2cal

         do i=1,natgas
            if (i.ne.pagas.and.i.ne.nhgas) then
               xxi = q - xyzgas(1,i)
               xxi2 = xxi*xxi
               dxc = distance(nhgas,i,xyzgas)
               dxc2 = dxc*dxc
               dxc3 = dxc2*dxc
               dxc5 = dxc2*dxc3
               d2vc1b = d2vc1b - chrgas(2,nhgas)*chrgas(2,i)*(dxc2 - 3.d0*xxi2)/dxc5
            endif
         enddo
         d2vc1b = d2vc1b*e2*ev2cal

         do i=1,natgas
            if (i.ne.pdgas.and.i.ne.nhgas) then
               xxi = q - xyzgas(1,i)
               xxi2 = xxi*xxi
               dxc = distance(nhgas,i,xyzgas)
               dxc2 = dxc*dxc
               dxc3 = dxc2*dxc
               dxc5 = dxc2*dxc3
               d2vc2a = d2vc2a - chrgas(3,nhgas)*chrgas(3,i)*(dxc2 - 3.d0*xxi2)/dxc5
            endif
         enddo
         d2vc2a = d2vc2a*e2*ev2cal

         do i=1,natgas
            if (i.ne.pagas.and.i.ne.nhgas) then
               xxi = q - xyzgas(1,i)
               xxi2 = xxi*xxi
               dxc = distance(nhgas,i,xyzgas)
               dxc2 = dxc*dxc
               dxc3 = dxc2*dxc
               dxc5 = dxc2*dxc3
               d2vc2b = d2vc2b - chrgas(4,nhgas)*chrgas(4,i)*(dxc2 - 3.d0*xxi2)/dxc5
            endif
         enddo
         d2vc2b = d2vc2b*e2*ev2cal

      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Second derivatives of the repulsion terms
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      xoha = dabs(dp2 - q)
      xohb = dabs(dp2 + q)
      d2vr1a =  broh1a*broh1a*droh1a*dexp(-broh1a*xoha) + 90.d0*croh1a/xoha**11.d0
      d2vr1b =  broh1b*broh1b*droh1b*dexp(-broh1b*xohb) + 90.d0*croh1b/xohb**11.d0
      d2vr2a =  broh2a*broh2a*droh2a*dexp(-broh2a*xoha) + 90.d0*croh2a/xoha**11.d0
      d2vr2b =  broh2b*broh2b*droh2b*dexp(-broh2b*xohb) + 90.d0*croh2b/xohb**11.d0

      !====> Diagonal elements

      d2h0_(1,1) = d2vm1a + d2vc1a + d2vr1a
      d2h0_(2,2) = d2vm1b + d2vc1b + d2vr1b
      d2h0_(3,3) = d2vm2a + d2vc2a + d2vr2a
      d2h0_(4,4) = d2vm2b + d2vc2b + d2vr2b

      !====> Off-diagonal elements

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Proton coupling 1a-1b/2a-2b
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      vpt1 = vpt0_1*dexp(-gapt_1*(q - qpt0_1)*(q - qpt0_1))
      vpt2 = vpt0_2*dexp(-gapt_2*(q - qpt0_2)*(q - qpt0_2))
      d2vpt1 = 2.d0*gapt_1*vpt1*(-1.d0+2.d0*gapt_1*(q - qpt0_1)*(q - qpt0_1))
      d2vpt2 = 2.d0*gapt_2*vpt2*(-1.d0+2.d0*gapt_2*(q - qpt0_2)*(q - qpt0_2))

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Electron coupling 1a-2a/1b-2b
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      veta = vet0_a*dexp(-gaet_a*(q - qet0_a)*(q - qet0_a))
      vetb = vet0_b*dexp(-gaet_b*(q - qet0_b)*(q - qet0_b))
      d2veta = 2.d0*gaet_a*veta*(-1.d0+2.d0*gaet_a*(q - qet0_a)*(q - qet0_a))
      d2vetb = 2.d0*gaet_b*vetb*(-1.d0+2.d0*gaet_b*(q - qet0_b)*(q - qet0_b))

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Mixed coupling 1a-2b/1b-2a
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      vetpt1 = vetpt0_1*dexp(-gaetpt_1*(q - qetpt0_1)*(q - qetpt0_1))
      vetpt2 = vetpt0_2*dexp(-gaetpt_2*(q - qetpt0_2)*(q - qetpt0_2))
      d2vetpt1 = 2.d0*gaetpt_1*vetpt1*(-1.d0+2.d0*gaetpt_1*(q - qetpt0_1)*(q - qetpt0_1))
      d2vetpt2 = 2.d0*gaetpt_2*vetpt2*(-1.d0+2.d0*gaetpt_2*(q - qetpt0_2)*(q - qetpt0_2))

      d2h0_(1,2) = d2vpt1
      d2h0_(2,1) = d2vpt1
      d2h0_(1,3) = d2veta
      d2h0_(3,1) = d2veta
      d2h0_(1,4) = d2vetpt1
      d2h0_(4,1) = d2vetpt1
      d2h0_(2,3) = d2vetpt2
      d2h0_(3,2) = d2vetpt2
      d2h0_(2,4) = d2vetb
      d2h0_(4,2) = d2vetb
      d2h0_(3,4) = d2vpt2
      d2h0_(4,3) = d2vpt2

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Convert units: kcal/mole/Angstroem**2 --> au/Bohr**2
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do i=1,4
         do j=1,4
            d2h0_(i,j) = d2h0_(i,j)*cal2au/(a2bohr*a2bohr)
         enddo
      enddo

      return

   end subroutine d2h0mat_mm5gen

   subroutine h0mat_leps0(h0_,dh0_,d2h0_,dgh0_,dg2h0_)
   !======================================================================!
   ! Calculates the gas phase Hamiltonian matrix and
   ! its first and second derivatives using
   ! modified LEPS potential (kcal/mole)
   !----------------------------------------------------------------------!

      implicit real(8) (a-h,o-z)
      implicit integer (i-n)

      real(8), intent(out), dimension(4,4) :: h0_, dh0_, d2h0_, dgh0_, dg2h0_

      ! local variables and arrays
      real(8), dimension(4,4) :: h, dh, d2h, dgh, dg2h
      real(8) :: sl(4,4), e(4), v(4,4), w1(4), w2(4)
      integer :: pdgas, nhgas, pagas, edgas, eagas
      real(8) mt,jah,jbh,jab

      ! local variables for parameters
      real(8) :: dah
      real(8) :: dbh
      real(8) :: dab
      real(8) :: bah
      real(8) :: bbh
      real(8) :: bab
      real(8) :: rah
      real(8) :: rbh
      real(8) :: rab
      real(8) :: r0Dd1
      real(8) :: r0Dd2
      real(8) :: r0Aa1
      real(8) :: r0Aa2
      real(8) :: kah
      real(8) :: kbh
      real(8) :: kab
      real(8) :: sigma
      real(8) :: VPT0
      real(8) :: ALPH
      real(8) :: VET0
      real(8) :: BETA
      real(8) :: kd1
      real(8) :: kd2
      real(8) :: ka1
      real(8) :: ka2
      real(8), dimension(4) :: corr

      ! pointer to the real set of parameters
      type(leps_2D), pointer :: par

      par => lepspar

      dah   = par%dah
      dbh   = par%dbh
      dab   = par%dab
      bah   = par%bah
      bbh   = par%bbh
      bab   = par%bab
      rah   = par%rah
      rbh   = par%rbh
      rab   = par%rab
      r0Dd1  = par%r0Dd1
      r0Aa1  = par%r0Aa1
      r0Dd2  = par%r0Dd2
      r0Aa2  = par%r0Aa2
      kah   = par%kah
      kbh   = par%kbh
      kab   = par%kab
      sigma = par%sigma
      VPT0  = par%VPT0
      ALPH  = par%ALPH
      VET0  = par%VET0
      BETA  = par%BETA
      kd1   = par%kd1
      kd2   = par%kd2
      ka1   = par%ka1
      ka2   = par%ka2
      corr  = par%corr

      h0_    = zero
      h      = zero
      dh0_   = zero
      dh     = zero
      d2h0_  = zero
      d2h    = zero
      dgh0_  = zero
      dgh    = zero
      dg2h0_ = zero
      dg2h   = zero
      sl     = zero

      pdgas = iptgas(1)
      nhgas = iptgas(2)
      pagas = iptgas(3)
      edgas = ietgas(1)
      eagas = ietgas(2)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! DP - a distance between PT donor (PDGAS) and acceptor (PAGAS)
      ! Note that PT interface is oriented along the x-axis.
      ! Q - antisymmetrical proton coordinate
      ! MT - Total mass (mass of donor and acceptor)
      ! XP - Conventional PT coordinate
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      mt = dm + am
      dp = dabs(xyzgas(1,pagas)-xyzgas(1,pdgas))
      q   = xyzgas(1,nhgas)

      !=================================================================
      ! Proton interface coordinates including COM
      ! XD = -(ma/mt)*dp    donor coordinate
      ! XA =  (md/mt)*dp    acceptor coordinate
      !=================================================================
      xd = xyzgas(1,pdgas)
      xa = xyzgas(1,pagas)

      !=================================================================
      ! Electron transfer pair coordinates including COM
      ! XED = electron donor coordinate
      ! XEA = electron acceptor coordinate
      !=================================================================
      xed = xyzgas(1,edgas)
      xea = xyzgas(1,eagas)

      !=================================================================
      ! Distances between ET donor/acceptor and PT donor/acceptor
      ! RDd = distance between ET donor and PT donor
      ! RAa = distance between ET acceptor and PT acceptor
      !=================================================================
      RDd = dabs(xed - xd)
      RAa = dabs(xea - xa)

      !=================================================================
      ! Coventional PT coordinate
      ! XP =  PT coordinate
      !=================================================================
      xp = dabs(xd - q)

      !=================================================================
      ! E1AH(E1BH,E1AB) - singlet state energy
      ! E3AH(E3BH,E3AB) - triplet state energy
      ! QAH(QBH,QAB) - Coulomb integral
      ! JAH(JBH,JAB) - exchange integral
      !=================================================================
      e1ah=dah*(dexp(-2.d0*bah*(xp-rah))-2.d0*dexp(-bah*(xp-rah)))
      e3ah=dah*(dexp(-2.d0*bah*(xp-rah))+2.d0*dexp(-bah*(xp-rah)))/2.d0
      qah=(e1ah*(1.d0+kah)+e3ah*(1.d0-kah))/2.d0
      jah=(e1ah*(1.d0+kah)-e3ah*(1.d0-kah))/2.d0

      rb = dp - xp
      e1bh=dbh*(dexp(-2.d0*bbh*(rb-rbh))-2.d0*dexp(-bbh*(rb-rbh)))
      e3bh=dbh*(dexp(-2.d0*bbh*(rb-rbh))+2.d0*dexp(-bbh*(rb-rbh)))/2.d0
      qbh =(e1bh*(1.d0+kbh)+e3bh*(1.d0-kbh))/2.d0
      jbh =(e1bh*(1.d0+kbh)-e3bh*(1.d0-kbh))/2.d0

      e1ab=dab*(dexp(-2.d0*bab*(dp-rab))-2.d0*dexp(-bab*(dp-rab)))
      e3ab=dab*(dexp(-2.d0*bab*(dp-rab))+2.d0*dexp(-bab*(dp-rab)))/2.d0
      qab=(e1ab*(1.d0+kab)+e3ab*(1.d0-kab))/2.d0
      jab=(e1ab*(1.d0+kab)-e3ab*(1.d0-kab))/2.d0

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! SCALING OF THE COULOMB AND EXCHANGE INTEGRALS
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      qah = qah/(1.d0+kah)
      jah = jah/(1.d0+kah)
      qbh = qbh/(1.d0+kbh)
      jbh = jbh/(1.d0+kbh)
      qab = qab/(1.d0+kab)
      jab = jab/(1.d0+kab)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Coulomb interaction between proton and electron donor and acceptor
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      vc1 = zero
      vc2 = zero

      do i=1,natgas
         if (i.ne.pdgas.and.i.ne.pagas.and.i.ne.nhgas) then
            dxc = distance(nhgas,i,xyzgas)
            vc1 = vc1 + chrgas(1,nhgas)*chrgas(1,i)/dxc
         endif
      enddo
      vc1 = vc1*e2*ev2cal

      do i=1,natgas
         if (i.ne.pdgas.and.i.ne.pagas.and.i.ne.nhgas) then
            dxc = distance(nhgas,i,xyzgas)
            vc2 = vc2 + chrgas(3,nhgas)*chrgas(3,i)/dxc
         endif
      enddo
      vc2 = vc2*e2*ev2cal

      !===================================================================
      ! Harmonic potentials
      ! UD1 - Harm. Potential between the ET/PT donor at 1st ET state
      ! UD2 - Harm. Potential between the ET/PT donor at 2nd ET state
      ! UA1 - Harm. Potential between the ET/PT acceptor at 1st ET state
      ! UA2 - Harm. Potential between the ET/PT acceptor at 2nd ET state
      !===================================================================
      UD1 = 0.5d0*kd1*(RDd-R0Dd1)*(RDd-R0Dd1)
      UD2 = 0.5d0*kd2*(RDd-R0Dd2)*(RDd-R0Dd2)
      UA1 = 0.5d0*ka1*(RAa-R0Aa1)*(RAa-R0Aa1)
      UA2 = 0.5d0*ka2*(RAa-R0Aa2)*(RAa-R0Aa2)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Coulombic and harmonic potentials (Non proton interface)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      vc1h = vc1 + ud1 + ua1
      vc2h = vc2 + ud2 + ua2

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Coupling (VET)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      RDA = DP + RDd + RAa
      vet = vet0*dexp(-beta*0.5d0*rda)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Hamiltonian Matrix in the non-orthogonal basis
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      h(1,1) = qah + qbh + qab - (jab+jbh)/2.d0 + jah + vc1h
      h(2,2) = qah + qbh + qab - (jab+jah)/2.d0 + jbh + vc1h
      h(3,3) = qah + qbh + qab - (jab+jbh)/2.d0 + jah + vc2h
      h(4,4) = qah + qbh + qab - (jab+jah)/2.d0 + jbh + vc2h

      h(1,2) = (qah+qbh+qab)/2.d0 - jab + (jah+jbh)/2.d0 + vc1h*sigma
      h(2,1) = h(1,2)
      h(1,3) = vet
      h(3,1) = vet

      vetovl = vet*sigma

      h(1,4) = vetovl
      h(4,1) = vetovl
      h(2,3) = vetovl
      h(3,2) = vetovl
      h(2,4) = vet
      h(4,2) = vet
      h(3,4) = (qah+qbh+qab)/2.d0 - jab + (jah+jbh)/2.d0 + vc2h*sigma
      h(4,3) = h(3,4)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Calculation of Derivatives
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      !=============================================================
      ! First derivative with respect to proton coordinate
      !=============================================================

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Derivatives of Morse and AntiMorse functions
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      DE1AH=-2.d0*DAH*BAH*(DEXP(-2.d0*BAH*(XP-RAH))-DEXP(-BAH*(XP-RAH)))
      DE3AH=-DAH*BAH*(DEXP(-2.d0*BAH*(XP-RAH))+DEXP(-BAH*(XP-RAH)))
      DQAH=(DE1AH*(1.d0+KAH)+DE3AH*(1.d0-KAH))/2.d0
      DJAH=(DE1AH*(1.d0+KAH)-DE3AH*(1.d0-KAH))/2.d0

      DE1BH=2.d0*DBH*BBH*(DEXP(-2.d0*BBH*(RB-RBH))-DEXP(-BBH*(RB-RBH)))
      DE3BH=DBH*BBH*(DEXP(-2.d0*BBH*(RB-RBH))+DEXP(-BBH*(RB-RBH)))
      DQBH=(DE1BH*(1.d0+KBH)+DE3BH*(1.d0-KBH))/2.d0
      DJBH=(DE1BH*(1.d0+KBH)-DE3BH*(1.d0-KBH))/2.d0

      DE1AB = zero
      DE3AB = zero
      DQAB  = zero
      DJAB  = zero

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! SCALING OF THE COULOMB AND EXCHANGE INTEGRALS
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DQAH = DQAH/(1.d0+KAH)
      DJAH = DJAH/(1.d0+KAH)
      DQBH = DQBH/(1.d0+KBH)
      DJBH = DJBH/(1.d0+KBH)
      DQAB = DQAB/(1.d0+KAB)
      DJAB = DJAB/(1.d0+KAB)


      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Diagonal elements of Leps Terms of first derivative
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      Dh(1,1)= DQAH + DQBH + DQAB - (DJAB + DJBH)/2.d0 + DJAH
      Dh(2,2)= DQAH + DQBH + DQAB - (DJAB + DJAH)/2.d0 + DJBH
      Dh(3,3)= DQAH + DQBH + DQAB - (DJAB + DJBH)/2.d0 + DJAH
      Dh(4,4)= DQAH + DQBH + DQAB - (DJAB + DJAH)/2.d0 + DJBH

      !=================================================================
      ! Derivative of Coulombic terms
      !=================================================================

      DVC1 = zero
      DVC2 = zero

      DO I=1,NATGAS
         IF (I.NE.PDGAS.AND.I.NE.NHGAS.AND.I.NE.PAGAS) THEN
            XXI = Q - XYZGAS(1,I)
            DXC = DISTANCE(NHGAS,I,XYZGAS)
            DXC3 = DXC**(3.D0)
            DVC1 = DVC1 - CHRGAS(1,NHGAS)*CHRGAS(1,I)*XXI/DXC3
         ENDIF
      ENDDO
      DVC1 = DVC1*E2*EV2CAL

      DO I=1,NATGAS
         IF (I.NE.PDGAS.AND.I.NE.NHGAS) THEN
            XXI = Q - XYZGAS(1,I)
            DXC = DISTANCE(NHGAS,I,XYZGAS)
            DXC3 = DXC**(3.D0)
            DVC2 = DVC2 - CHRGAS(3,NHGAS)*CHRGAS(3,I)*XXI/DXC3
         ENDIF
      ENDDO
      DVC2 = DVC2*E2*EV2CAL

      !================================================================
      ! Diagonal derivatives (kcal/mol/A) including Coloumb Energy
      !================================================================

      DH(1,1) = Dh(1,1) + DVC1
      DH(2,2) = Dh(2,2) + DVC1
      DH(3,3) = Dh(3,3) + DVC2
      DH(4,4) = Dh(4,4) + DVC2

      !================================================================
      ! Off-Diagonal derivatives (kcal/mol/A)
      !================================================================

      DH(1,2)=(DQAH+DQBH+DQAB)/2.d0-DJAB+(DJAH+DJBH)/2.d0+DVC1*sigma
      DH(3,4)=(DQAH+DQBH+DQAB)/2.d0-DJAB+(DJAH+DJBH)/2.d0+DVC2*sigma

      DH(1,3) = zero
      DH(1,4) = zero
      DH(2,1) = DH(1,2)
      DH(2,3) = zero
      DH(2,4) = zero
      DH(3,1) = DH(3,1)
      DH(3,2) = DH(2,3)
      DH(4,1) = DH(1,4)
      DH(4,2) = DH(2,4)
      DH(4,3) = DH(3,4)

      !=============================================================
      ! Second Derivatives with respect to the proton coordinate
      !=============================================================

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Morse and Antimorse for second deriv.
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      D2E1AH=-2.d0*DAH*BAH*BAH*(-2.d0*DEXP(-2.d0*BAH*(XP-RAH)) + DEXP(-BAH*(XP-RAH)))
      D2E3AH=-DAH*BAH*BAH*(-2.d0*DEXP(-2.d0*BAH*(XP-RAH)) - DEXP(-BAH*(XP-RAH)))
      D2QAH=(D2E1AH*(1.d0+KAH)+D2E3AH*(1.d0-KAH))/2.d0
      D2JAH=(D2E1AH*(1.d0+KAH)-D2E3AH*(1.d0-KAH))/2.d0

      D2E1BH=2.d0*DBH*BBH*BBH*(2.D0*DEXP(-2.d0*BBH*(RB-RBH))-DEXP(-BBH*(RB-RBH)))
      D2E3BH=DBH*BBH*BBH*(2.D0*DEXP(-2.d0*BBH*(RB-RBH))+DEXP(-BBH*(RB-RBH)))
      D2QBH=(D2E1BH*(1.d0+KBH)+D2E3BH*(1.d0-KBH))/2.d0
      D2JBH=(D2E1BH*(1.d0+KBH)-D2E3BH*(1.d0-KBH))/2.d0

      D2E1AB = zero
      D2E3AB = zero
      D2QAB  = zero
      D2JAB  = zero

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! SCALING OF THE COULOMB AND EXCHANGE INTEGRALS
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      D2QAH = D2QAH/(1.d0+KAH)
      D2JAH = D2JAH/(1.d0+KAH)
      D2QBH = D2QBH/(1.d0+KBH)
      D2JBH = D2JBH/(1.d0+KBH)
      D2QAB = D2QAB/(1.d0+KAB)
      D2JAB = D2JAB/(1.d0+KAB)

      !===============================================================
      ! Diagonal second derivatives (kcal/mol/A^2) with Leps parameters
      !===============================================================
      D2H(1,1)= D2QAH+D2QBH+D2QAB-(D2JAB+D2JBH)/2.d0+D2JAH
      D2H(2,2)= D2QAH+D2QBH+D2QAB-(D2JAB+D2JAH)/2.d0+D2JBH
      D2H(3,3)= D2QAH+D2QBH+D2QAB-(D2JAB+D2JBH)/2.d0+D2JAH
      D2H(4,4)= D2QAH+D2QBH+D2QAB-(D2JAB+D2JAH)/2.d0+D2JBH

      !===============================================================
      ! 2-nd Derivative of Coulombic
      !===============================================================
      D2VC1 = zero
      D2VC2 = zero

      DO I=1,NATGAS
         IF (I.NE.PDGAS.AND.I.NE.NHGAS.AND.I.NE.PAGAS) THEN
            XXI = Q - XYZGAS(1,I)
            XXI2 = XXI*XXI
            DXC = DISTANCE(NHGAS,I,XYZGAS)
            DXC2 = DXC*DXC
            DXC3 = DXC2*DXC
            DXC5 = DXC2*DXC3
            D2VC1 = D2VC1 - CHRGAS(1,NHGAS)*CHRGAS(1,I)*(DXC2 - 3.D0*XXI2)/DXC5
         ENDIF
      ENDDO
      D2VC1 = D2VC1*E2*EV2CAL

      DO I=1,NATGAS
         IF (I.NE.PDGAS.AND.I.NE.NHGAS.AND.I.NE.PAGAS) THEN
            XXI = Q - XYZGAS(1,I)
            XXI2 = XXI*XXI
            DXC = DISTANCE(NHGAS,I,XYZGAS)
            DXC2 = DXC*DXC
            DXC3 = DXC2*DXC
            DXC5 = DXC2*DXC3
            D2VC2 = D2VC2 - CHRGAS(3,NHGAS)*CHRGAS(3,I)*(DXC2 - 3.D0*XXI2)/DXC5
         ENDIF
      ENDDO
      D2VC2 = D2VC2*E2*EV2CAL

      !================================================================
      ! Diagonal elements (kcal/mol/A^2) including Coloumb energy
      !================================================================
      D2H(1,1) = D2h(1,1) + D2VC1
      D2H(2,2) = D2h(2,2) + D2VC1
      D2H(3,3) = D2h(3,3) + D2VC2
      D2H(4,4) = D2h(4,4) + D2VC2

      !================================================================
      ! Off-Diagonal elements (kcal/mol/A^2)
      !================================================================
      D2H(1,2)=(D2QAH+D2QBH+D2QAB)/2.d0-D2JAB+(D2JAH+D2JBH)/2.d0+D2VC1*sigma
      D2H(3,4)=(D2QAH+D2QBH+D2QAB)/2.d0-D2JAB+(D2JAH+D2JBH)/2.d0+D2VC2*sigma

      D2H(1,3) = zero
      D2H(1,4) = zero
      D2H(2,1) = D2H(1,2)
      D2H(2,3) = zero
      D2H(2,4) = zero
      D2H(3,1) = D2H(3,1)
      D2H(3,2) = D2H(2,3)
      D2H(4,1) = D2H(1,4)
      D2H(4,2) = D2H(2,4)
      D2H(4,3) = D2H(3,4)

      !=============================================================
      ! First Derivative with respect to the gating coordinate
      !=============================================================
      DXP = AM/MT
      DGE1AH=-2.d0*DAH*DXP*BAH*(DEXP(-2.d0*BAH*(XP-RAH))-DEXP(-BAH*(XP-RAH)))
      DGE3AH=-DAH*DXP*BAH*(DEXP(-2.d0*BAH*(XP-RAH))+DEXP(-BAH*(XP-RAH)))
      DGQAH=(DGE1AH*(1.d0+KAH)+DGE3AH*(1.d0-KAH))/2.d0
      DGJAH=(DGE1AH*(1.d0+KAH)-DGE3AH*(1.d0-KAH))/2.d0

      DRB = 1.d0 - DXP
      DGE1BH=-2.d0*DBH*DRB*BBH*(DEXP(-2.d0*BBH*(RB-RBH))-DEXP(-BBH*(RB-RBH)))
      DGE3BH=-DBH*DRB*BBH*(DEXP(-2.d0*BBH*(RB-RBH))+DEXP(-BBH*(RB-RBH)))
      DGQBH=(DGE1BH*(1.d0+KBH)+DGE3BH*(1.d0-KBH))/2.d0
      DGJBH=(DGE1BH*(1.d0+KBH)-DGE3BH*(1.d0-KBH))/2.d0

      DGE1AB=-2.d0*DAB*BAB*(DEXP(-2.d0*BAB*(DP-RAB))-DEXP(-BAB*(DP-RAB)))
      DGE3AB=-DAB*BAB*(DEXP(-2.d0*BAB*(DP-RAB))+DEXP(-BAB*(DP-RAB)))
      DGQAB=(DGE1AB*(1.d0+KAB)+DGE3AB*(1.d0-KAB))/2.d0
      DGJAB=(DGE1AB*(1.d0+KAB)-DGE3AB*(1.d0-KAB))/2.d0

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! SCALING OF THE COULOMB AND EXCHANGE INTEGRALS
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DGQAH = DGQAH/(1.d0+KAH)
      DGJAH = DGJAH/(1.d0+KAH)
      DGQBH = DGQBH/(1.d0+KBH)
      DGJBH = DGJBH/(1.d0+KBH)
      DGQAB = DGQAB/(1.d0+KAB)
      DGJAB = DGJAB/(1.d0+KAB)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Diagonal elements derivatives with Leps parameters
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DGh(1,1)= DGQAH + DGQBH + DGQAB - (DGJAB + DGJBH)/2.d0 + DGJAH
      DGh(2,2)= DGQAH + DGQBH + DGQAB - (DGJAB + DGJAH)/2.d0 + DGJBH
      DGh(3,3)= DGQAH + DGQBH + DGQAB - (DGJAB + DGJBH)/2.d0 + DGJAH
      DGh(4,4)= DGQAH + DGQBH + DGQAB - (DGJAB + DGJAH)/2.d0 + DGJBH

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Derivatives of harmonic potentials
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DRDd = -AM/MT
      DRAa = -DM/MT

      DUD1 = DRDd * kd1 * (RDd-R0Dd1)
      DUD2 = DRDd * kd2 * (RDd-R0Dd2)
      DUA1 = DRAa * ka1 * (RAa-R0Aa1)
      DUA2 = DRAa * ka2 * (RAa-R0Aa2)

      !===========================================================
      ! Diagonal elements (kcal/mol/A) including Harmonics
      !===========================================================
      DGH(1,1) = DGH(1,1) + DUD1
      DGH(2,2) = DGH(2,2) + DUD2
      DGH(3,3) = DGH(3,3) + DUA1
      DGH(4,4) = DGH(4,4) + DUA2

      !===========================================================
      ! Off-Diagonal elements (kcal/mol/A)
      !===========================================================
      DGH(1,2) = (DQAH+DQBH+DQAB)/2.d0-DJAB+(DJAH+DJBH)/2.d0+DUD1*sigma
      DGH(3,4) = (DQAH+DQBH+DQAB)/2.d0-DJAB+(DJAH+DJBH)/2.d0+DUA1*sigma
      DGH(1,3) = zero
      DGH(1,4) = zero
      DGH(2,1) = DGH(1,2)
      DH(2,3)  = zero
      DH(2,4)  = zero
      DGH(3,1) = DGH(3,1)
      DGH(3,2) = DGH(2,3)
      DGH(4,1) = DGH(1,4)
      DGH(4,2) = DGH(2,4)
      DGH(4,3) = DGH(3,4)

      !=============================================================
      ! Second derivatives with respect to the gating coordinate
      !=============================================================
      D2GE1AH=-2.d0*DAH*DXP*DXP*BAH*BAH*(-2.d0*DEXP(-2.d0*BAH*(XP-RAH)) + DEXP(-BAH*(XP-RAH)))
      D2GE3AH=-DAH*DXP*DXP*BAH*BAH*(-2.d0*DEXP(-2.d0*BAH*(XP-RAH)) - DEXP(-BAH*(XP-RAH)))
      D2GQAH=(D2GE1AH*(1.d0+KAH)+D2GE3AH*(1.d0-KAH))/2.d0
      D2GJAH=(D2GE1AH*(1.d0+KAH)-D2GE3AH*(1.d0-KAH))/2.d0

      D2GE1BH=-2.d0*DBH*DRB*DRB*BBH*BBH*(-2.d0*DEXP(-2.d0*BBH*(RB-RBH)) + DEXP(-BBH*(RB-RBH)))
      D2GE3BH=-DBH*DRB*DRB*BBH*BBH*(-2.d0*DEXP(-2.d0*BBH*(RB-RBH)) - DEXP(-BBH*(RB-RBH)))
      D2GQBH=(D2GE1BH*(1.d0+KBH)+D2GE3BH*(1.d0-KBH))/2.d0
      D2GJBH=(D2GE1BH*(1.d0+KBH)-D2GE3BH*(1.d0-KBH))/2.d0

      D2GE1AB=-2.d0*DAB*BAB*BAB*(-2.d0*DEXP(-2.d0*BAB*(DP-RAB)) + DEXP(-BAB*(DP-RAB)))
      D2GE3AB=-DAB*BAB*BAB*(-2.d0*DEXP(-2.d0*BAB*(DP-RAB)) - DEXP(-BAB*(DP-RAB)))
      D2GQAB=(D2GE1AB*(1.d0+KAB)+D2GE3AB*(1.d0-KAB))/2.d0
      D2GJAB=(D2GE1AB*(1.d0+KAB)-D2GE3AB*(1.d0-KAB))/2.d0

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! SCALING OF THE COULOMB AND EXCHANGE INTEGRALS
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      D2GQAH = D2GQAH/(1.d0+KAH)
      D2GJAH = D2GJAH/(1.d0+KAH)
      D2GQBH = D2GQBH/(1.d0+KBH)
      D2GJBH = D2GJBH/(1.d0+KBH)
      D2GQAB = D2GQAB/(1.d0+KAB)
      D2GJAB = D2GJAB/(1.d0+KAB)

      !==========================================================
      ! Diagonal elements second derivatives with Leps parameters
      !==========================================================
      DG2H(1,1)= D2GQAH+D2GQBH+D2GQAB-(D2GJAB+D2JBH)/2.d0+D2GJAH
      DG2H(2,2)= D2GQAH+D2GQBH+D2GQAB-(D2GJAB+D2GJAH)/2.d0+D2GJBH
      DG2H(3,3)= D2GQAH+D2GQBH+D2GQAB-(D2GJAB+D2GJBH)/2.d0+D2GJAH
      DG2H(4,4)= D2GQAH+D2GQBH+D2GQAB-(D2GJAB+D2GJAH)/2.d0+D2GJBH

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! 2-nd Derivatives of harmonic potentials
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      D2UD1 = DRDd *DRDd* kd1
      D2UD2 = DRDd *DRDd* kd2
      D2UA1 = DRAa *DRAa* ka1
      D2UA2 = DRAa *DRAa* ka2

      !================================================================
      ! Diagonal second derivatives (kcal/mol/A^2) INCLUDING HARMONICS
      !================================================================
      D2H(1,1) = D2h(1,1) + D2UD1
      D2H(2,2) = D2h(2,2) + D2UD2
      D2H(3,3) = D2h(3,3) + D2UA1
      D2H(4,4) = D2h(4,4) + D2UA2

      !===============================================================
      ! Off-Diagonal elements (kcal/mol/A^2)
      !===============================================================
      DG2H(1,2)=(D2QAH+D2QBH+D2QAB)/2.d0-D2JAB+(D2JAH+D2JBH)/2.d0+D2UD1*sigma
      DG2H(3,4)=(D2QAH+D2QBH+D2QAB)/2.d0-D2JAB+(D2JAH+D2JBH)/2.d0+D2UA1*sigma

      DG2H(1,3) = zero
      DG2H(1,4) = zero
      DG2H(2,1) = DG2H(1,2)
      DG2H(2,3) = zero
      DG2H(2,4) = zero
      DG2H(3,1) = DG2H(3,1)
      DG2H(3,2) = DG2H(2,3)
      DG2H(4,1) = DG2H(1,4)
      DG2H(4,2) = DG2H(2,4)
      DG2H(4,3) = DG2H(3,4)

      !=============================================================
      ! Lowdin orthogonalization
      !=============================================================

      ! Overlap matrix (S)

      ! / 1  s  0  0 \
      ! | s  1  0  0 |
      ! | 0  0  1  s |
      ! \ 0  0  s  1 /

      ! S^(-1/2) matrix

      s1 = 1.d0/dsqrt(1.d0-sigma)
      s2 = 1.d0/dsqrt(1.d0+sigma)

      do i=1,4
         sl(i,i) = 0.5d0*(s1+s2)
      enddo

      sl(1,2) = 0.5d0*(-s1+s2)
      sl(2,1) = 0.5d0*(-s1+s2)
      sl(3,4) = 0.5d0*(-s1+s2)
      sl(4,3) = 0.5d0*(-s1+s2)


      ! Hamiltonian matrix in the orthogonalized basis

      do i=1,4
         do j=1,4
            h0ij = zero
            dh0ij = zero
            d2h0ij = zero
            dgh0ij = zero
            dg2h0ij = zero
            do k=1,4
               do l=1,4
                  h0ij    = h0ij + sl(i,k)*h(k,l)*sl(l,j)
                  dh0ij   = dh0ij + sl(i,k)*dh(k,l)*sl(l,j)
                  d2h0ij  = d2h0ij + sl(i,k)*d2h(k,l)*sl(l,j)
                  dgh0ij  = dgh0ij + sl(i,k)*dgh(k,l)*sl(l,j)
                  dg2h0ij = dg2h0ij + sl(i,k)*dg2h(k,l)*sl(l,j)
               enddo
            enddo
            h0_(i,j)    = h0ij
            dh0_(i,j)   = dh0ij
            d2h0_(i,j)  = d2h0ij
            dgh0_(i,j)  = dgh0ij
            dg2h0_(i,j) = dg2h0ij
         enddo
      enddo

      h0_(1,1) = h0_(1,1) + corr(1)
      h0_(2,2) = h0_(2,2) + corr(2)
      h0_(3,3) = h0_(3,3) + corr(3)
      h0_(4,4) = h0_(4,4) + corr(4)

      fd1 = 0.5d0*(1.d0 - DTanh(h0_(1,2)))
      fd2 = 0.5d0*(1.d0 - DTanh(h0_(3,4)))

      h0_(1,2) = VPT0*h0_(1,2)*fd1
      h0_(2,1) = h0_(1,2)
      h0_(3,4) = VPT0*h0_(3,4)*fd2
      h0_(4,3) = h0_(3,4)

      dh0_(1,2) = VPT0*dh0_(1,2)
      dh0_(2,1) = dh0_(1,2)
      dh0_(3,4) = VPT0*dh0_(3,4)
      dh0_(4,3) = dh0_(3,4)

      d2h0_(1,2) = VPT0*d2h0_(1,2)
      d2h0_(2,1) = d2h0_(1,2)
      d2h0_(3,4) = VPT0*d2h0_(3,4)
      d2h0_(4,3) = d2h0_(3,4)

      dgh0_(1,2) = VPT0*dgh0_(1,2)
      dgh0_(2,1) = dgh0_(1,2)
      dgh0_(3,4) = VPT0*dgh0_(3,4)
      dgh0_(4,3) = dgh0_(3,4)

      dg2h0_(1,2) = VPT0*dg2h0_(1,2)
      dg2h0_(2,1) = dg2h0_(1,2)
      dg2h0_(3,4) = VPT0*dg2h0_(3,4)
      dg2h0_(4,3) = dg2h0_(3,4)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! 1-st Derivative of off diagonal elements including damping
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Sech1 = 1.d0/DCosh(h0_(1,2))
      Sech2 = 1.d0/DCosh(h0_(3,4))

      dfd1=-0.5d0*(Sech1*Sech1)*Dh0_(1,2)
      dfd2=-0.5d0*(Sech2*Sech2)*Dh0_(3,4)

      DH0_(1,2)=fd1*dh0_(1,2) + dfd1*h0_(1,2)
      DH0_(3,4)=fd2*dh0_(3,4) + dfd2*h0_(3,4)

      DGH0_(1,2)=fd1*dgh0_(1,2) + dfd1*h0_(1,2)
      DGH0_(3,4)=fd2*dgh0_(3,4) + dfd2*h0_(3,4)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! 2-nd Deriv. of off diagonal elements including damping
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      dSech1 = -Sech1*dtanh(h0_(1,2))
      dSech2 = -Sech2*dtanh(h0_(3,4))
      d2f1 = -Sech1*dSech1*D2h0_(1,2)
      d2f2 = -Sech2*dSech2*D2h0_(3,4)

      D2H0_(1,2)=fd1*d2h0_(1,2) + dfd1*dh0_(1,2) + d2f1*h0_(1,2) + dfd1*d2h0_(1,2)
      D2H0_(3,4)=fd2*d2h0_(3,4) + dfd2*dh0_(3,4) + d2f2*h0_(3,4) + dfd2*d2h0_(3,4)

      DG2H0_(1,2)=fd1*dg2h0_(1,2) + dfd1*dgh0_(1,2) + d2f1*h0_(1,2) + dfd1*dg2h0_(1,2)
      DG2H0_(3,4)=fd2*dg2h0_(3,4) + dfd2*dgh0_(3,4) + d2f2*h0_(3,4) + dfd2*dg2h0_(3,4)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Convert units for derivatives: kcal/mole/Angstroem -> au/Bohr
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do i=1,4
         do j=1,4
            dh0_(i,j)   = dh0_(i,j)  *cal2au/a2bohr
            d2h0_(i,j)  = d2h0_(i,j) *cal2au/a2bohr
            dgh0_(i,j)  = dgh0_(i,j) *cal2au/a2bohr
            dg2h0_(i,j) = dg2h0_(i,j)*cal2au/a2bohr
         enddo
      enddo

      return

   end subroutine h0mat_leps0


   subroutine h0mat_leps2(h0_,dh0_,d2h0_,dgh0_,dg2h0_)
   !======================================================================!
   ! Calculates the gas phase Hamiltonian matrix and
   ! its first and second derivatives using
   ! modified LEPS potential (kcal/mole)
   !======================================================================!
      implicit real(8) (a-h,o-z)
      implicit integer (i-n)

      real(8), intent(out), dimension(4,4) :: h0_, dh0_, d2h0_, dgh0_, dg2h0_

      ! local variables and arrays
      real(8), dimension(4,4) :: h, dh, d2h, dgh, dg2h
      integer :: pdgas, nhgas, pagas, edgas, eagas
      real(8)  :: mt, jah, jbh, jab

      ! local variables for parameters
      real(8) :: dah
      real(8) :: dbh
      real(8) :: dab
      real(8) :: bah
      real(8) :: bbh
      real(8) :: bab
      real(8) :: rah
      real(8) :: rbh
      real(8) :: rab
      real(8) :: r0Dd1
      real(8) :: r0Aa1
      real(8) :: r0Dd2
      real(8) :: r0Aa2
      real(8) :: kah
      real(8) :: kbh
      real(8) :: kab
      real(8) :: sigma
      real(8) :: VPT0
      real(8) :: ALPH
      real(8) :: VET0
      real(8) :: BETA
      real(8) :: kd1
      real(8) :: kd2
      real(8) :: ka1
      real(8) :: ka2
      real(8), dimension(4) :: corr

      ! pointer to the real set of parameters
      type(leps_2D), pointer :: par

      par => lepspar

      dah   = par%dah
      dbh   = par%dbh
      dab   = par%dab
      bah   = par%bah
      bbh   = par%bbh
      bab   = par%bab
      rah   = par%rah
      rbh   = par%rbh
      rab   = par%rab
      r0Dd1  = par%r0Dd1
      r0Aa1  = par%r0Aa1
      r0Dd2  = par%r0Dd2
      r0Aa2  = par%r0Aa2
      kah   = par%kah
      kbh   = par%kbh
      kab   = par%kab
      sigma = par%sigma
      VPT0  = par%VPT0
      ALPH  = par%ALPH
      VET0  = par%VET0
      BETA  = par%BETA
      kd1   = par%kd1
      kd2   = par%kd2
      ka1   = par%ka1
      ka2   = par%ka2
      corr  = par%corr

      h0_    = zero
      h      = zero
      dh0_   = zero
      dh     = zero
      d2h0_  = zero
      d2h    = zero
      dgh0_  = zero
      dgh    = zero
      dg2h0_ = zero
      dg2h   = zero

      pdgas = iptgas(1)
      nhgas = iptgas(2)
      pagas = iptgas(3)
      edgas = ietgas(1)
      eagas = ietgas(2)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! DP - a distance between PT donor (PDGAS) and acceptor (PAGAS)
      ! Note that PT interface is oriented along the x-axis.
      ! Q - antisymmetrical proton coordinate
      ! MT - Total mass (mass of donor and acceptor)
      ! XP - Conventional PT coordinate
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      mt = dm + am
      dp = dabs(xyzgas(1,pagas)-xyzgas(1,pdgas))
      q   = xyzgas(1,nhgas)

      !=================================================================
      ! Proton interface coordinates including COM
      ! XD = -(ma/mt)*dp    donor coordinate
      ! XA =  (md/mt)*dp    acceptor coordinate
      !=================================================================
      xd = xyzgas(1,pdgas)
      xa = xyzgas(1,pagas)

      !=================================================================
      ! Electron transfer pair coordinates including COM
      ! XED = electron donor coordinate
      ! XEA = electron acceptor coordinate
      !=================================================================
      xed = xyzgas(1,edgas)
      xea = xyzgas(1,eagas)

      !=================================================================
      ! Distances between ET donor/acceptor and PT donor/acceptor
      ! RDd = distance between ET donor and PT donor
      ! RAa = distance between ET acceptor and PT acceptor
      !=================================================================
      RDd = xd  - xed
      RAa = xea - xa

      !=================================================================
      ! Coventional PT coordinate
      ! XP =  PT coordinate
      !=================================================================
      xp = q - xd

      !=================================================================
      ! E1AH(E1BH,E1AB) - singlet state energy
      ! E3AH(E3BH,E3AB) - triplet state energy
      ! QAH(QBH,QAB) - Coulomb integral
      ! JAH(JBH,JAB) - exchange integral
      !=================================================================
      e1ah=dah*(dexp(-2.d0*bah*(xp-rah))-2.d0*dexp(-bah*(xp-rah)))
      e3ah=dah*(dexp(-2.d0*bah*(xp-rah))+2.d0*dexp(-bah*(xp-rah)))/2.d0
      qah=(e1ah*(1.d0+kah)+e3ah*(1.d0-kah))/2.d0
      jah=(e1ah*(1.d0+kah)-e3ah*(1.d0-kah))/2.d0

      rb = dp - xp
      e1bh=dbh*(dexp(-2.d0*bbh*(rb-rbh))-2.d0*dexp(-bbh*(rb-rbh)))
      e3bh=dbh*(dexp(-2.d0*bbh*(rb-rbh))+2.d0*dexp(-bbh*(rb-rbh)))/2.d0
      qbh =(e1bh*(1.d0+kbh)+e3bh*(1.d0-kbh))/2.d0
      jbh =(e1bh*(1.d0+kbh)-e3bh*(1.d0-kbh))/2.d0

      e1ab=dab*(dexp(-2.d0*bab*(dp-rab))-2.d0*dexp(-bab*(dp-rab)))
      e3ab=dab*(dexp(-2.d0*bab*(dp-rab))+2.d0*dexp(-bab*(dp-rab)))/2.d0
      qab=(e1ab*(1.d0+kab)+e3ab*(1.d0-kab))/2.d0
      jab=(e1ab*(1.d0+kab)-e3ab*(1.d0-kab))/2.d0

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! SCALING OF THE COULOMB AND EXCHANGE INTEGRALS
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      qah = qah/(1.d0+kah)
      jah = jah/(1.d0+kah)
      qbh = qbh/(1.d0+kbh)
      jbh = jbh/(1.d0+kbh)
      qab = qab/(1.d0+kab)
      jab = jab/(1.d0+kab)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Coulomb interaction between proton and electron donor and acceptor
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      vc1 = zero
      vc2 = zero

      do i=1,natgas
         if (i.ne.pdgas.and.i.ne.pagas.and.i.ne.nhgas) then
            dxc = distance(nhgas,i,xyzgas)
            vc1 = vc1 + chrgas(1,nhgas)*chrgas(1,i)/dxc
         endif
      enddo
      vc1 = vc1*e2*ev2cal

      do i=1,natgas
         if (i.ne.pdgas.and.i.ne.pagas.and.i.ne.nhgas) then
            dxc = distance(nhgas,i,xyzgas)
            vc2 = vc2 + chrgas(3,nhgas)*chrgas(3,i)/dxc
         endif
      enddo
      vc2 = vc2*e2*ev2cal

      !===================================================================
      ! Harmonic potentials
      ! UD1 - Harm. Potential between the ET/PT donor at 1st ET state
      ! UD2 - Harm. Potential between the ET/PT donor at 2nd ET state
      ! UA1 - Harm. Potential between the ET/PT acceptor at 1st ET state
      ! UA2 - Harm. Potential between the ET/PT acceptor at 2nd ET state
      !===================================================================
      ud1 = 0.5d0*kd1*(RDd-R0Dd1)*(RDd-R0Dd1)
      ud2 = 0.5d0*kd2*(RDd-R0Dd2)*(RDd-R0Dd2)
      ua1 = 0.5d0*ka1*(RAa-R0Aa1)*(RAa-R0Aa1)
      ua2 = 0.5d0*ka2*(RAa-R0Aa2)*(RAa-R0Aa2)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Coulombic and harmonic potentials (Non proton interface)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !vc1h = vc1 + ud1 + ua1
      !vc2h = vc2 + ud2 + ua2

      vc1h = ud1 + ua1
      vc2h = ud2 + ua2

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Coupling (VET)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      rda = dabs(xea - xed)
      vet = vet0*dexp(-beta*0.5d0*rda)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Hamiltonian Matrix in the non-orthogonal basis
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      deltar = dp - rab

      h(1,1) = qah + qbh + qab - (jab+jbh)/2.d0 + jah + vc1h
      h(2,2) = qah + qbh + qab - (jab+jah)/2.d0 + jbh + vc1h
      h(3,3) = qah + qbh + qab - (jab+jbh)/2.d0 + jah + vc2h
      h(4,4) = qah + qbh + qab - (jab+jah)/2.d0 + jbh + vc2h

      h(1,2) = -vpt0*dexp(-alph*deltar)
      h(2,1) = h(1,2)
      h(1,3) = vet
      h(3,1) = vet

      vetovl = vet*sigma

      h(1,4) = vetovl
      h(4,1) = vetovl
      h(2,3) = vetovl
      h(3,2) = vetovl
      h(2,4) = vet
      h(4,2) = vet
      h(3,4) = -vpt0*dexp(-alph*deltar)
      h(4,3) = h(3,4)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Calculation of Derivatives
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      !=============================================================
      ! First derivative with respect to proton coordinate
      !=============================================================

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Derivatives of Morse and AntiMorse functions
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      de1ah=-2.d0*dah*bah*(dexp(-2.d0*bah*(xp-rah))-dexp(-bah*(xp-rah)))
      de3ah=-dah*bah*(dexp(-2.d0*bah*(xp-rah))+dexp(-bah*(xp-rah)))
      dqah=(de1ah*(1.d0+kah)+de3ah*(1.d0-kah))/2.d0
      djah=(de1ah*(1.d0+kah)-de3ah*(1.d0-kah))/2.d0

      de1bh=2.d0*dbh*bbh*(dexp(-2.d0*bbh*(rb-rbh))-dexp(-bbh*(rb-rbh)))
      de3bh=dbh*bbh*(dexp(-2.d0*bbh*(rb-rbh))+dexp(-bbh*(rb-rbh)))
      dqbh=(de1bh*(1.d0+kbh)+de3bh*(1.d0-kbh))/2.d0
      djbh=(de1bh*(1.d0+kbh)-de3bh*(1.d0-kbh))/2.d0

      de1ab = zero
      de3ab = zero
      dqab  = zero
      djab  = zero

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! SCALING OF THE COULOMB AND EXCHANGE INTEGRALS
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      dqah = dqah/(1.d0+kah)
      djah = djah/(1.d0+kah)
      dqbh = dqbh/(1.d0+kbh)
      djbh = djbh/(1.d0+kbh)
      dqab = dqab/(1.d0+kab)
      djab = djab/(1.d0+kab)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Diagonal elements of Leps Terms of first derivative
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      dh(1,1)= dqah + dqbh + dqab - (djab + djbh)/2.d0 + djah
      dh(2,2)= dqah + dqbh + dqab - (djab + djah)/2.d0 + djbh
      dh(3,3)= dqah + dqbh + dqab - (djab + djbh)/2.d0 + djah
      dh(4,4)= dqah + dqbh + dqab - (djab + djah)/2.d0 + djbh

      !=================================================================
      ! Derivative of Coulombic terms
      !=================================================================
      dvc1 = zero
      dvc2 = zero

      do i=1,natgas
         if (i.ne.pdgas.and.i.ne.nhgas.and.i.ne.pagas) then
            xxi = q - xyzgas(1,i)
            dxc = distance(nhgas,i,xyzgas)
            dxc3 = dxc**(3.d0)
            dvc1 = dvc1 - chrgas(1,nhgas)*chrgas(1,i)*xxi/dxc3
         endif
      enddo
      dvc1 = dvc1*e2*ev2cal

      do i=1,natgas
         if (i.ne.pdgas.and.i.ne.nhgas) then
            xxi = q - xyzgas(1,i)
            dxc = distance(nhgas,i,xyzgas)
            dxc3 = dxc**(3.d0)
            dvc2 = dvc2 - chrgas(3,nhgas)*chrgas(3,i)*xxi/dxc3
         endif
      enddo
      dvc2 = dvc2*e2*ev2cal

      !================================================================
      ! Diagonal derivatives (kcal/mol/A) including Coloumb Energy
      !================================================================
      dh(1,1) = dh(1,1) + dvc1
      dh(2,2) = dh(2,2) + dvc1
      dh(3,3) = dh(3,3) + dvc2
      dh(4,4) = dh(4,4) + dvc2

      !================================================================
      ! Off-Diagonal derivatives (kcal/mol/A)
      !================================================================
      dh(1,2) = zero
      dh(3,4) = zero
      dh(1,3) = zero
      dh(1,4) = zero
      dh(2,1) = dh(1,2)
      dh(2,3) = zero
      dh(2,4) = zero
      dh(3,1) = dh(3,1)
      dh(3,2) = dh(2,3)
      dh(4,1) = dh(1,4)
      dh(4,2) = dh(2,4)
      dh(4,3) = dh(3,4)

      !=============================================================
      ! Second Derivatives with respect to the proton coordinate
      !=============================================================

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Morse and Antimorse for second deriv.
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      d2e1ah=-2.d0*dah*bah*bah*(-2.d0*dexp(-2.d0*bah*(xp-rah)) + dexp(-bah*(xp-rah)))
      d2e3ah=-dah*bah*bah*(-2.d0*dexp(-2.d0*bah*(xp-rah)) - dexp(-bah*(xp-rah)))
      d2qah=(d2e1ah*(1.d0+kah)+d2e3ah*(1.d0-kah))/2.d0
      d2jah=(d2e1ah*(1.d0+kah)-d2e3ah*(1.d0-kah))/2.d0

      d2e1bh=2.d0*dbh*bbh*bbh*(2.d0*dexp(-2.d0*bbh*(rb-rbh))-dexp(-bbh*(rb-rbh)))
      d2e3bh=dbh*bbh*bbh*(2.d0*dexp(-2.d0*bbh*(rb-rbh))+dexp(-bbh*(rb-rbh)))
      d2qbh=(d2e1bh*(1.d0+kbh)+d2e3bh*(1.d0-kbh))/2.d0
      d2jbh=(d2e1bh*(1.d0+kbh)-d2e3bh*(1.d0-kbh))/2.d0

      d2e1ab = zero
      d2e3ab = zero
      d2qab  = zero
      d2jab  = zero

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! SCALING OF THE COULOMB AND EXCHANGE INTEGRALS
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      d2qah = d2qah/(1.d0+kah)
      d2jah = d2jah/(1.d0+kah)
      d2qbh = d2qbh/(1.d0+kbh)
      d2jbh = d2jbh/(1.d0+kbh)
      d2qab = d2qab/(1.d0+kab)
      d2jab = d2jab/(1.d0+kab)

      !===============================================================
      ! Diagonal second derivatives (kcal/mol/A^2) with Leps parameters
      !===============================================================
      d2h(1,1) = d2qah+d2qbh+d2qab-(d2jab+d2jbh)/2.d0+d2jah
      d2h(2,2) = d2qah+d2qbh+d2qab-(d2jab+d2jah)/2.d0+d2jbh
      d2h(3,3) = d2qah+d2qbh+d2qab-(d2jab+d2jbh)/2.d0+d2jah
      d2h(4,4) = d2qah+d2qbh+d2qab-(d2jab+d2jah)/2.d0+d2jbh

      !===============================================================
      ! 2-nd Derivative of Coulombic
      !===============================================================
      d2vc1 = zero
      d2vc2 = zero

      do i=1,natgas
         if (i.ne.pdgas.and.i.ne.nhgas.and.i.ne.pagas) then
            xxi = q - xyzgas(1,i)
            xxi2 = xxi*xxi
            dxc = distance(nhgas,i,xyzgas)
            dxc2 = dxc*dxc
            dxc3 = dxc2*dxc
            dxc5 = dxc2*dxc3
            d2vc1 = d2vc1 - chrgas(1,nhgas)*chrgas(1,i)*(dxc2 - 3.d0*xxi2)/dxc5
         endif
      enddo
      d2vc1 = d2vc1*e2*ev2cal

      do i=1,natgas
         if (i.ne.pdgas.and.i.ne.nhgas.and.i.ne.pagas) then
            xxi = q - xyzgas(1,i)
            xxi2 = xxi*xxi
            dxc = distance(nhgas,i,xyzgas)
            dxc2 = dxc*dxc
            dxc3 = dxc2*dxc
            dxc5 = dxc2*dxc3
            d2vc2 = d2vc2 - chrgas(3,nhgas)*chrgas(3,i)*(dxc2 - 3.d0*xxi2)/dxc5
         endif
      enddo
      d2vc2 = d2vc2*e2*ev2cal

      !================================================================
      ! Diagonal elements (kcal/mol/A^2) including Coloumb energy
      !================================================================
      d2h(1,1) = d2h(1,1) + d2vc1
      d2h(2,2) = d2h(2,2) + d2vc1
      d2h(3,3) = d2h(3,3) + d2vc2
      d2h(4,4) = d2h(4,4) + d2vc2

      !================================================================
      ! Off-Diagonal elements (kcal/mol/A^2)
      !================================================================
      d2h(1,2) = zero
      d2h(3,4) = zero
      d2h(1,3) = zero
      d2h(1,4) = zero
      d2h(2,1) = d2h(1,2)
      d2h(2,3) = zero
      d2h(2,4) = zero
      d2h(3,1) = d2h(3,1)
      d2h(3,2) = d2h(2,3)
      d2h(4,1) = d2h(1,4)
      d2h(4,2) = d2h(2,4)
      d2h(4,3) = d2h(3,4)

      !=============================================================
      ! first derivative with respect to the gating coordinate
      !=============================================================
      dxp = am/mt
      dge1ah=-2.d0*dah*dxp*bah*(dexp(-2.d0*bah*(xp-rah))-dexp(-bah*(xp-rah)))
      dge3ah=-dah*dxp*bah*(dexp(-2.d0*bah*(xp-rah))+dexp(-bah*(xp-rah)))
      dgqah=(dge1ah*(1.d0+kah)+dge3ah*(1.d0-kah))/2.d0
      dgjah=(dge1ah*(1.d0+kah)-dge3ah*(1.d0-kah))/2.d0

      drb = 1.d0 - dxp
      dge1bh=-2.d0*dbh*drb*bbh*(dexp(-2.d0*bbh*(rb-rbh))-dexp(-bbh*(rb-rbh)))
      dge3bh=-dbh*drb*bbh*(dexp(-2.d0*bbh*(rb-rbh))+dexp(-bbh*(rb-rbh)))
      dgqbh=(dge1bh*(1.d0+kbh)+dge3bh*(1.d0-kbh))/2.d0
      dgjbh=(dge1bh*(1.d0+kbh)-dge3bh*(1.d0-kbh))/2.d0

      dge1ab=-2.d0*dab*bab*(dexp(-2.d0*bab*(dp-rab))-dexp(-bab*(dp-rab)))
      dge3ab=-dab*bab*(dexp(-2.d0*bab*(dp-rab))+dexp(-bab*(dp-rab)))
      dgqab=(dge1ab*(1.d0+kab)+dge3ab*(1.d0-kab))/2.d0
      dgjab=(dge1ab*(1.d0+kab)-dge3ab*(1.d0-kab))/2.d0

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! scaling of the coulomb and exchange integrals
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      dgqah = dgqah/(1.d0+kah)
      dgjah = dgjah/(1.d0+kah)
      dgqbh = dgqbh/(1.d0+kbh)
      dgjbh = dgjbh/(1.d0+kbh)
      dgqab = dgqab/(1.d0+kab)
      dgjab = dgjab/(1.d0+kab)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! diagonal elements derivatives with leps parameters
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      dgh(1,1)= dgqah + dgqbh + dgqab - (dgjab + dgjbh)/2.d0 + dgjah
      dgh(2,2)= dgqah + dgqbh + dgqab - (dgjab + dgjah)/2.d0 + dgjbh
      dgh(3,3)= dgqah + dgqbh + dgqab - (dgjab + dgjbh)/2.d0 + dgjah
      dgh(4,4)= dgqah + dgqbh + dgqab - (dgjab + dgjah)/2.d0 + dgjbh

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Derivatives of harmonic potentials
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DRDd = -am/mt
      DRAa = -dm/mt

      dud1 = DRDd * kd1 * (RDd-R0Dd1)
      dud2 = DRDd * kd2 * (RDd-R0Dd2)
      dua1 = DRAa * ka1 * (RAa-R0Aa1)
      dua2 = DRAa * ka2 * (RAa-R0Aa2)

      !===========================================================
      ! Diagonal elements (kcal/mol/A) including Harmonics
      !===========================================================
      dgh(1,1) = dgh(1,1) + dud1
      dgh(2,2) = dgh(2,2) + dud2
      dgh(3,3) = dgh(3,3) + dua1
      dgh(4,4) = dgh(4,4) + dua2

      !===========================================================
      ! Off-Diagonal elements (kcal/mol/A)
      !===========================================================
      dgh(1,2) = -alph*h(1,2)
      dgh(3,4) = -alph*h(3,4)
      dgh(1,3) = zero
      dgh(1,4) = zero
      dgh(2,1) = dgh(1,2)
      dgh(2,3) = zero
      dgh(2,4) = zero
      dgh(3,1) = dgh(3,1)
      dgh(3,2) = dgh(2,3)
      dgh(4,1) = dgh(1,4)
      dgh(4,2) = dgh(2,4)
      dgh(4,3) = dgh(3,4)

      !=============================================================
      ! Second derivatives with respect to the gating coordinate
      !=============================================================
      d2ge1ah=-2.d0*dah*dxp*dxp*bah*bah*(-2.d0*dexp(-2.d0*bah*(xp-rah)) + dexp(-bah*(xp-rah)))
      d2ge3ah=-dah*dxp*dxp*bah*bah*(-2.d0*dexp(-2.d0*bah*(xp-rah)) - dexp(-bah*(xp-rah)))
      d2gqah=(d2ge1ah*(1.d0+kah)+d2ge3ah*(1.d0-kah))/2.d0
      d2gjah=(d2ge1ah*(1.d0+kah)-d2ge3ah*(1.d0-kah))/2.d0

      d2ge1bh=-2.d0*dbh*drb*drb*bbh*bbh*(-2.d0*dexp(-2.d0*bbh*(rb-rbh)) + dexp(-bbh*(rb-rbh)))
      d2ge3bh=-dbh*drb*drb*bbh*bbh*(-2.d0*dexp(-2.d0*bbh*(rb-rbh)) - dexp(-bbh*(rb-rbh)))
      d2gqbh=(d2ge1bh*(1.d0+kbh)+d2ge3bh*(1.d0-kbh))/2.d0
      d2gjbh=(d2ge1bh*(1.d0+kbh)-d2ge3bh*(1.d0-kbh))/2.d0

      d2ge1ab=-2.d0*dab*bab*bab*(-2.d0*dexp(-2.d0*bab*(dp-rab)) + dexp(-bab*(dp-rab)))
      d2ge3ab=-dab*bab*bab*(-2.d0*dexp(-2.d0*bab*(dp-rab)) - dexp(-bab*(dp-rab)))
      d2gqab=(d2ge1ab*(1.d0+kab)+d2ge3ab*(1.d0-kab))/2.d0
      d2gjab=(d2ge1ab*(1.d0+kab)-d2ge3ab*(1.d0-kab))/2.d0

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! SCALING OF THE COULOMB AND EXCHANGE INTEGRALS
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      d2gqah = d2gqah/(1.d0+kah)
      d2gjah = d2gjah/(1.d0+kah)
      d2gqbh = d2gqbh/(1.d0+kbh)
      d2gjbh = d2gjbh/(1.d0+kbh)
      d2gqab = d2gqab/(1.d0+kab)
      d2gjab = d2gjab/(1.d0+kab)

      !==========================================================
      ! Diagonal elements second derivatives with Leps parameters
      !==========================================================
      dg2h(1,1)= d2gqah+d2gqbh+d2gqab-(d2gjab+d2jbh)/2.d0+d2gjah
      dg2h(2,2)= d2gqah+d2gqbh+d2gqab-(d2gjab+d2gjah)/2.d0+d2gjbh
      dg2h(3,3)= d2gqah+d2gqbh+d2gqab-(d2gjab+d2gjbh)/2.d0+d2gjah
      dg2h(4,4)= d2gqah+d2gqbh+d2gqab-(d2gjab+d2gjah)/2.d0+d2gjbh

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! 2-nd derivatives of harmonic potentials
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      d2ud1 = drdd *drdd* kd1
      d2ud2 = drdd *drdd* kd2
      d2ua1 = draa *draa* ka1
      d2ua2 = draa *draa* ka2

      !================================================================
      ! Diagonal second derivatives (kcal/mol/A^2) INCLUDING HARMONICS
      !================================================================
      d2h(1,1) = d2h(1,1) + d2ud1
      d2h(2,2) = d2h(2,2) + d2ud2
      d2h(3,3) = d2h(3,3) + d2ua1
      d2h(4,4) = d2h(4,4) + d2ua2

      !===============================================================
      ! Off-Diagonal elements (kcal/mol/A^2)
      !===============================================================
      dg2h(1,2) = -alph*dgh(1,2)
      dg2h(3,4) = -alph*dgh(3,4)

      dg2h(1,3) = zero
      dg2h(1,4) = zero
      dg2h(2,1) = dg2h(1,2)
      dg2h(2,3) = zero
      dg2h(2,4) = zero
      dg2h(3,1) = dg2h(3,1)
      dg2h(3,2) = dg2h(2,3)
      dg2h(4,1) = dg2h(1,4)
      dg2h(4,2) = dg2h(2,4)
      dg2h(4,3) = dg2h(3,4)

      ! Hamiltonian matrix and its derivatives

      h0_ = h
      h0_(1,1) = h0_(1,1) + corr(1)
      h0_(2,2) = h0_(2,2) + corr(2)
      h0_(3,3) = h0_(3,3) + corr(3)
      h0_(4,4) = h0_(4,4) + corr(4)

      dh0_   = dh
      d2h0_  = d2h
      dgh0_  = dgh
      dg2h0_ = dg2h

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Convert units for derivatives: kcal/mole/Angstroem -> au/Bohr
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do i=1,4
         do j=1,4
            dh0_(i,j)   = dh0_(i,j)  *cal2au/a2bohr
            d2h0_(i,j)  = d2h0_(i,j) *cal2au/a2bohr
            dgh0_(i,j)  = dgh0_(i,j) *cal2au/a2bohr
            dg2h0_(i,j) = dg2h0_(i,j)*cal2au/a2bohr
         enddo
      enddo

      return
   
   end subroutine h0mat_leps2


   subroutine h0mat_hyb(h0_,dh0_,d2h0_,dgh0_,dg2h0_)
   !======================================================================C
   ! Calculates the gas phase Hamiltonian matrix and
   ! its first and second derivatives using
   ! hybrid diabatic potentials (kcal/mole)
   !======================================================================C
      implicit none
      real(8), intent(out), dimension(4,4) :: h0_, dh0_, d2h0_, dgh0_, dg2h0_

      integer :: pdgas, nhgas, pagas, edgas, eagas
      integer :: i, j, k, l
      real(8)  :: mt, mah, mbh, mab
      real(8)  :: dp, q, xd, xa, xed, xea, xp, rb, rda, vet
      real(8)  :: deltar1, deltar2, vpt1, vpt2, vept
      real(8)  :: cspeed, convm, cjoule, unit
      real(8)  :: fcah, fcbh, fcab

      real(8) :: bah, bbh, e1ah, e3ah, qah, jah
      real(8) :: e1bh, e3bh, qbh, jbh, uah, ubh, uab1, uab2, de1ah, de3ah, dqah
      real(8) :: djah, de1bh, de3bh, dqbh, djbh, d2e1ah, d2e3ah, d2qah, d2jah
      real(8) :: d2e1bh, d2e3bh, d2qbh, d2jbh, dxp, dge1ah, dge3ah, dgqah, dgjah
      real(8) :: drb, dge1bh, dge3bh, dgqbh, dgjbh, d2ge1ah, d2ge3ah, d2gqah, d2gjah
      real(8) :: d2ge1bh, d2ge3bh, d2gqbh, d2gjbh
      real(8) :: s1, s2, h0ij, dh0ij, d2h0ij, dgh0ij, dg2h0ij

      real(8), dimension(4,4) :: h, dh, d2h, dgh, dg2h, sl

      ! local variables for parameters
      real(8) :: omega_ah            ! A-H frequency
      real(8) :: omega_bh            ! B-H frequency
      real(8) :: omega_ab            ! A-B frequency
      real(8) :: dah                 ! A-H dissocation energy
      real(8) :: dbh                 ! B-H dissocation energy
      real(8) :: rah                 ! A-H equilibrium distance
      real(8) :: rbh                 ! B-H equilibrium distance
      real(8) :: rab1                ! A-B equilibrium distance for ET reactants
      real(8) :: rab2                ! A-B equilibrium distance for ET products
      real(8) :: kah                 ! A-H Sato parameter
      real(8) :: kbh                 ! B-H Sato parameter
      real(8) :: sigma               ! PT overlap integral
      real(8) :: vpt0                ! PT off-diagonal coupling constant
      real(8) :: alph                ! PT off-diagonal coupling exponent
      real(8) :: vet0                ! ET off-diagonal coupling constant
      real(8) :: beta                ! ET off-diagonal coupling exponent
      real(8) :: vept0               ! EPT off-diagonal coupling constant
      real(8), dimension(4) :: corr  ! energy bias parameters for diabatic states

      ! pointer to the real set of parameters
      type(hybrid_2D), pointer :: par

      par => hybpar

      omega_ah = par%omega_ah
      omega_bh = par%omega_bh
      omega_ab = par%omega_ab
      dah      = par%dah
      dbh      = par%dbh
      rah      = par%rah
      rbh      = par%rbh
      rab1     = par%rab1
      rab2     = par%rab2
      kah      = par%kah
      kbh      = par%kbh
      sigma    = par%sigma
      vpt0     = par%vpt0
      alph     = par%alph
      vet0     = par%vet0
      beta     = par%beta
      vept0    = par%vept0
      corr     = par%corr

      h0_    = zero
      dh0_   = zero
      d2h0_  = zero
      dgh0_  = zero
      dg2h0_ = zero

      pdgas = iptgas(1)
      nhgas = iptgas(2)
      pagas = iptgas(3)
      edgas = ietgas(1)
      eagas = ietgas(2)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! DP - a distance between PT donor (PDGAS) and acceptor (PAGAS)
      ! Note that PT interface is oriented along the x-axis.
      ! Q - antisymmetrical proton coordinate
      ! MT - Total mass (mass of donor and acceptor)
      ! XP - Conventional PT coordinate
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      mt = dm + am
      dp = dabs(xyzgas(1,pagas)-xyzgas(1,pdgas))
      q  = xyzgas(1,nhgas)

      !=================================================================
      ! Proton interface coordinates including COM
      ! XD = -(ma/mt)*dp    donor coordinate
      ! XA =  (md/mt)*dp    acceptor coordinate
      !=================================================================
      xd = xyzgas(1,pdgas)
      xa = xyzgas(1,pagas)

      !=================================================================
      ! Electron transfer pair coordinates
      ! XED = electron donor coordinate
      ! XEA = electron acceptor coordinate
      !=================================================================
      xed = xyzgas(1,edgas)
      xea = xyzgas(1,eagas)

      !=================================================================
      ! Coventional PT coordinate (A-H distance)
      ! XP =  PT coordinate
      !=================================================================
      xp = q - xd
      rb = dp - xp

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Couplings
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      rda = abs(xea - xed)
      vet = vet0*exp(-half*beta*rda)
      deltar1 = dp - rab1
      deltar2 = dp - rab2
      vpt1 = -vpt0*exp(-alph*deltar1)
      vpt2 = -vpt0*exp(-alph*deltar2)
      vept = vet*sigma

      ! unit = cm2ev*ev2au*cm2ev*ev2au*au2cal/bohr2a/bohr2a
      cspeed = 2.99792458d0
      convm = 0.910953d0
      cjoule = 4.3597482d0
      unit = 1.d-12*4.d0*pi*pi*convm*cspeed*cspeed*au2cal/cjoule

      ! beta Morse parameters for AH and BH bonds
      mah = dm*pm/(dm+pm)
      bah = a2bohr*ev2au*cm2ev*omega_ah*sqrt(half*mah*au2cal/dah)

      mbh = am*pm/(am+pm)
      bbh = a2bohr*ev2au*cm2ev*omega_bh*sqrt(half*mbh*au2cal/dbh)

      mab = am*dm/mt
      fcab = mab*omega_ab*omega_ab*unit

      !=================================================================
      ! e1ah(e1bh,e1ab) - singlet state energy
      ! e3ah(e3bh,e3ab) - triplet state energy
      ! qah(qbh,qab)    - Coulomb integral
      ! jah(jbh,jab)    - exchange integral
      !=================================================================
      e1ah = dah*(exp(-two*bah*(xp-rah))-two*exp(-bah*(xp-rah)))
      e3ah = half*dah*(exp(-two*bah*(xp-rah))+two*exp(-bah*(xp-rah)))
      qah = half*(e1ah*(one+kah)+e3ah*(one-kah))
      jah = half*(e1ah*(one+kah)-e3ah*(one-kah))

      rb = dp - xp
      e1bh = dbh*(exp(-two*bbh*(rb-rbh))-two*exp(-bbh*(rb-rbh)))
      e3bh = half*dbh*(exp(-two*bbh*(rb-rbh))+two*exp(-bbh*(rb-rbh)))
      qbh = half*(e1bh*(one+kbh)+e3bh*(one-kbh))
      jbh = half*(e1bh*(one+kbh)-e3bh*(one-kbh))

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! scaling of the coulomb and exchange integrals
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      qah = qah/(one+kah)
      jah = jah/(one+kah)
      qbh = qbh/(one+kbh)
      jbh = jbh/(one+kbh)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! AH and BH diabatic LEPS potentials
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      uah = qah + jah + qbh - half*jbh
      ubh = qbh + jbh + qah - half*jah

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Harmonic diabatic potentials for AB bond
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      uab1 = half*fcab*deltar1*deltar1
      uab2 = half*fcab*deltar2*deltar2

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Hamiltonian Matrix in the non-orthgonal LEPS basis
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      h(1,1) = uah + uab1 + corr(1)
      h(2,2) = ubh + uab1 + corr(2)
      h(3,3) = uah + uab2 + corr(3)
      h(4,4) = ubh + uab2 + corr(4)

      h(1,2) = vpt1
      h(2,1) = h(1,2)
      h(1,3) = vet
      h(3,1) = vet
      h(1,4) = vept
      h(4,1) = vept
      h(2,3) = vept
      h(3,2) = vept
      h(2,4) = vept
      h(4,2) = vept
      h(3,4) = vpt2
      h(4,3) = h(3,4)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Calculation of Derivatives
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      !=============================================================
      ! First derivative with respect to the proton coordinate
      !=============================================================

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Derivatives of Morse and AntiMorse functions
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      de1ah = -two*dah*bah*(exp(-two*bah*(xp-rah))-exp(-bah*(xp-rah)))
      de3ah = -dah*bah*(exp(-two*bah*(xp-rah))+exp(-bah*(xp-rah)))
      dqah = half*(de1ah*(one+kah)+de3ah*(one-kah))
      djah = half*(de1ah*(one+kah)-de3ah*(one-kah))

      de1bh = two*dbh*bbh*(exp(-two*bbh*(rb-rbh))-exp(-bbh*(rb-rbh)))
      de3bh = dbh*bbh*(exp(-two*bbh*(rb-rbh))+exp(-bbh*(rb-rbh)))
      dqbh = half*(de1bh*(one+kbh)+de3bh*(one-kbh))
      djbh = half*(de1bh*(one+kbh)-de3bh*(one-kbh))

      dqah = dqah/(one+kah)
      djah = djah/(one+kah)
      dqbh = dqbh/(one+kbh)
      djbh = djbh/(one+kbh)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Diagonal first derivatives: LEPS terms
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      dh(1,1) = dqah + dqbh - half*djbh + djah
      dh(2,2) = dqah + dqbh - half*djah + djbh
      dh(3,3) = dqah + dqbh - half*djbh + djah
      dh(4,4) = dqah + dqbh - half*djah + djbh

      !================================================================
      ! Off-Diagonal derivatives (kcal/mol/A)
      !================================================================
      dh(1,2) = zero
      dh(3,4) = zero
      dh(1,3) = zero
      dh(1,4) = zero
      dh(2,1) = dh(1,2)
      dh(2,3) = zero
      dh(2,4) = zero
      dh(3,1) = dh(3,1)
      dh(3,2) = dh(2,3)
      dh(4,1) = dh(1,4)
      dh(4,2) = dh(2,4)
      dh(4,3) = dh(3,4)

      !=============================================================
      ! Second Derivatives with respect to the proton coordinate
      !=============================================================

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Morse and Anti-Morse terms
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      d2e1ah = -two*dah*bah*bah*(-two*exp(-two*bah*(xp-rah))+exp(-bah*(xp-rah)))
      d2e3ah = -dah*bah*bah*(-two*exp(-two*bah*(xp-rah))-exp(-bah*(xp-rah)))
      d2qah = half*(d2e1ah*(one+kah)+d2e3ah*(one-kah))
      d2jah = half*(d2e1ah*(one+kah)-d2e3ah*(one-kah))

      d2e1bh = two*dbh*bbh*bbh*(two*exp(-two*bbh*(rb-rbh))-exp(-bbh*(rb-rbh)))
      d2e3bh = dbh*bbh*bbh*(two*exp(-two*bbh*(rb-rbh))+exp(-bbh*(rb-rbh)))
      d2qbh = half*(d2e1bh*(one+kbh)+d2e3bh*(one-kbh))
      d2jbh = half*(d2e1bh*(one+kbh)-d2e3bh*(one-kbh))

      d2qah = d2qah/(one+kah)
      d2jah = d2jah/(one+kah)
      d2qbh = d2qbh/(one+kbh)
      d2jbh = d2jbh/(one+kbh)

      !================================================================
      ! Off-Diagonal second derivatives (kcal/mol/A^2)
      !================================================================
      d2h(1,2) = zero
      d2h(3,4) = zero
      d2h(1,3) = zero
      d2h(1,4) = zero
      d2h(2,1) = d2h(1,2)
      d2h(2,3) = zero
      d2h(2,4) = zero
      d2h(3,1) = d2h(3,1)
      d2h(3,2) = d2h(2,3)
      d2h(4,1) = d2h(1,4)
      d2h(4,2) = d2h(2,4)
      d2h(4,3) = d2h(3,4)

      !=============================================================
      ! First Derivative with respect to the gating coordinate
      !=============================================================

      dxp = am/mt
      dge1ah = -two*dah*dxp*bah*(exp(-two*bah*(xp-rah))-exp(-bah*(xp-rah)))
      dge3ah = -dah*dxp*bah*(exp(-two*bah*(xp-rah))+exp(-bah*(xp-rah)))
      dgqah = half*(dge1ah*(one+kah)+dge3ah*(one-kah))
      dgjah = half*(dge1ah*(one+kah)-dge3ah*(one-kah))

      drb = one - dxp
      dge1bh = -two*dbh*drb*bbh*(exp(-two*bbh*(rb-rbh))-exp(-bbh*(rb-rbh)))
      dge3bh = -dbh*drb*bbh*(exp(-two*bbh*(rb-rbh))+exp(-bbh*(rb-rbh)))
      dgqbh = half*(dge1bh*(one+kbh)+dge3bh*(one-kbh))
      dgjbh = half*(dge1bh*(one+kbh)-dge3bh*(one-kbh))

      dgqah = dgqah/(one+kah)
      dgjah = dgjah/(one+kah)
      dgqbh = dgqbh/(one+kbh)
      dgjbh = dgjbh/(one+kbh)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! diagonal elements derivatives
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      dgh(1,1) = dgqah + dgqbh - half*dgjbh + dgjah + fcab*deltar1
      dgh(2,2) = dgqah + dgqbh - half*dgjah + dgjbh + fcab*deltar1
      dgh(3,3) = dgqah + dgqbh - half*dgjbh + dgjah + fcab*deltar2
      dgh(4,4) = dgqah + dgqbh - half*dgjah + dgjbh + fcab*deltar2

      !===========================================================
      ! Off-Diagonal elements (kcal/mol/A)
      !===========================================================
      dgh(1,2) = -alph*vpt1
      dgh(3,4) = -alph*vpt2
      dgh(1,3) = zero
      dgh(1,4) = zero
      dgh(2,1) = dgh(1,2)
      dgh(2,3) = zero
      dgh(2,4) = zero
      dgh(3,1) = dgh(3,1)
      dgh(3,2) = dgh(2,3)
      dgh(4,1) = dgh(1,4)
      dgh(4,2) = dgh(2,4)
      dgh(4,3) = dgh(3,4)

      !=============================================================
      ! Second derivatives with respect to the gating coordinate
      !=============================================================
      d2ge1ah = -two*dah*dxp*dxp*bah*bah*(-two*exp(-two*bah*(xp-rah))+exp(-bah*(xp-rah)))
      d2ge3ah = -dah*dxp*dxp*bah*bah*(-two*exp(-two*bah*(xp-rah))-exp(-bah*(xp-rah)))
      d2gqah = half*(d2ge1ah*(one+kah)+d2ge3ah*(one-kah))
      d2gjah = half*(d2ge1ah*(one+kah)-d2ge3ah*(one-kah))

      d2ge1bh = -two*dbh*drb*drb*bbh*bbh*(-two*exp(-two*bbh*(rb-rbh))+exp(-bbh*(rb-rbh)))
      d2ge3bh = -dbh*drb*drb*bbh*bbh*(-two*exp(-two*bbh*(rb-rbh))-exp(-bbh*(rb-rbh)))
      d2gqbh = half*(d2ge1bh*(one+kbh)+d2ge3bh*(one-kbh))
      d2gjbh = half*(d2ge1bh*(one+kbh)-d2ge3bh*(one-kbh))

      d2gqah = d2gqah/(one+kah)
      d2gjah = d2gjah/(one+kah)
      d2gqbh = d2gqbh/(one+kbh)
      d2gjbh = d2gjbh/(one+kbh)

      dg2h(1,1) = d2gqah + d2gqbh - half*d2gjbh + d2gjah + fcab
      dg2h(2,2) = d2gqah + d2gqbh - half*d2gjah + d2gjbh + fcab
      dg2h(3,3) = d2gqah + d2gqbh - half*d2gjbh + d2gjah + fcab
      dg2h(4,4) = d2gqah + d2gqbh - half*d2gjah + d2gjbh + fcab

      !===============================================================
      ! Off-Diagonal elements (kcal/mol/A^2)
      !===============================================================
      dg2h(1,2) = -alph*dgh(1,2)
      dg2h(3,4) = -alph*dgh(3,4)
      dg2h(1,3) = zero
      dg2h(1,4) = zero
      dg2h(2,1) = dg2h(1,2)
      dg2h(2,3) = zero
      dg2h(2,4) = zero
      dg2h(3,1) = dg2h(3,1)
      dg2h(3,2) = dg2h(2,3)
      dg2h(4,1) = dg2h(1,4)
      dg2h(4,2) = dg2h(2,4)
      dg2h(4,3) = dg2h(3,4)

      !=============================================================
      ! Lowdin orthogonalization
      !=============================================================

      ! Overlap matrix (S)

      ! / 1  s  0  0 \
      ! | s  1  0  0 |
      ! | 0  0  1  s |
      ! \ 0  0  s  1 /

      ! S^(-1/2) matrix
      sl = zero

      s1 = one/sqrt(one-sigma)
      s2 = one/sqrt(one+sigma)

      do i=1,4
         sl(i,i) = half*(s1+s2)
      enddo

      sl(1,2) = half*(s2-s1)
      sl(2,1) = half*(s2-s1)
      sl(3,4) = half*(s2-s1)
      sl(4,3) = half*(s2-s1)

      ! Hamiltonian matrix and its derivatives in the orthogonalized basis

      do i=1,4
         do j=1,4
            h0ij = zero
            dh0ij = zero
            d2h0ij = zero
            dgh0ij = zero
            dg2h0ij = zero
            do k=1,4
               do l=1,4
                  h0ij    = h0ij    + sl(i,k)*h(k,l)   *sl(l,j)
                  dh0ij   = dh0ij   + sl(i,k)*dh(k,l)  *sl(l,j)
                  d2h0ij  = d2h0ij  + sl(i,k)*d2h(k,l) *sl(l,j)
                  dgh0ij  = dgh0ij  + sl(i,k)*dgh(k,l) *sl(l,j)
                  dg2h0ij = dg2h0ij + sl(i,k)*dg2h(k,l)*sl(l,j)
               enddo
            enddo
            h0_(i,j)    = h0ij
            dh0_(i,j)   = dh0ij
            d2h0_(i,j)  = d2h0ij
            dgh0_(i,j)  = dgh0ij
            dg2h0_(i,j) = dg2h0ij
         enddo
      enddo

      !h0_(1,1) = h0_(1,1) + corr(1)
      !h0_(2,2) = h0_(2,2) + corr(2)
      !h0_(3,3) = h0_(3,3) + corr(3)
      !h0_(4,4) = h0_(4,4) + corr(4)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Convert units for derivatives: kcal/mole/Angstroem -> au/Bohr
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do i=1,4
         do j=1,4
            dh0_(i,j)   = dh0_(i,j)  *cal2au/a2bohr
            d2h0_(i,j)  = d2h0_(i,j) *cal2au/a2bohr
            dgh0_(i,j)  = dgh0_(i,j) *cal2au/a2bohr
            dg2h0_(i,j) = dg2h0_(i,j)*cal2au/a2bohr
         enddo
      enddo

      return

   end subroutine h0mat_hyb


   function distance(n1,n2,xyz) result(r)

      implicit none
      integer, intent(in) :: n1, n2
      real(8), intent(in), dimension(:,:) :: xyz
      real(8) :: r

      integer :: j
      real(8) :: dxyz

      r = zero
      do j=1,3
         dxyz = xyz(j,n1) - xyz(j,n2)
         r = r + dxyz*dxyz
      enddo
      
      r = dsqrt(r)

      return

   end function distance

end module potential
