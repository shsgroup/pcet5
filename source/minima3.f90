subroutine minima3
!====================================================================
!  Sets parameters for the routine which finds minima
!  on the spesified three-dimensional free energy surface
!  using canonical Newton-Raphson minimization algorithm
!  for solvent coordinates and a simple minimization on the grid
!  for the gating coordinate
!
!  KEYWORD: MIN3
!
!  OPTIONS:
!
!  ADIAB=N   - search on the N-th adiabatic free energy surface
!
!  DIAB2=N/M - search on the N-th ET diabatic surface
!              within M-th set. "M" must be either 1 (1a/1b)
!              or 2 (2a/2b).
!
!  DIAB4=N/M - search on the N-th ET/PT diabatic surface
!              within M-th set. "M" can be 1 (1a), 2 (1b)
!              3 (2a) or 4 (2b).
!
!  ZP0=<VALUE>/<SCALE> - initial approximation for ZP coordinate
!                        and a scaling factor. Scaling factor
!                        is intended to equalize the scales for
!                        both coordinates, sometimes it makes
!                        life easier for the minimization
!                        procedure)
!
!  ZE0=<VALUE>/<SCALE> - the same for ZE coordinate.
!
!  R0=<VALUE> - initial value for the gating coordinate.
!
!  ACC=<VALUE> - the desired accuracy of the minimization
!
!  SLIM=<VALUE> - the maximum length of the Newton-Raphson step
!
!  IPRINT=N - printing level in LBFGS minimization routine (default=-1)
!
!  MCORR=m  - number of corrections used in the limited memory matrix.
!             Values of m < 3  are not recommended, and large values
!             of m can result in excessive computing time. The range
!             3 <= m <= 20 is recommended (default=5)

!  FACTR=<value> - accuracy of LBFGS minimization (factor of machine precision)
!                  1.d+12 for low accuracy;
!                  1.d+7  for moderate accuracy; 
!                  1.d+1  for extremely high accuracy
!                  Default value is 1.d+6
!
!  PGTOL=<value> - projected gradient tolerance in LBFGS minimization (positive)
!                  The iteration will stop when max{|proj g_i|, i=1,...,n} <= pgtol
!                  where pg_i is the ith component of the projected gradient.
!
!  MAXIT=N  - maximum number of Newton-Raphson steps allowed
!
!  WEIGHTS - prints out the results of the wavefunction analysis
!            to the standard output file
!
!  WEIGHTS=<filename> - prints out the results of the wavefunction
!                       analysis to the external file
!
!  WAVEFUN=<filename> - dump the proton vibrational wavefunctions
!                       to the external file <filename>
!
!--------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2012-03-13 22:01:18 $
!  $Revision: 5.3 $
!  $Log: not supported by cvs2svn $
!  Revision 5.2  2010/10/28 21:29:36  souda
!  First (working and hopefully bug-free) source of PCET 5.x
!
!
!====================================================================

   use pardim
   use strings
   use keys
   use cst
   use control
   use quantum
   use fesmin_3d
   use solmat

   implicit none

   character(1024) :: options
   character( 40)  :: fname
   character(  5)  :: mode
   logical         :: adiab, diab2, diab4, weights, wavefun, gfix

   integer :: ierr, ikey, iadiab, istate, iset, idiab2, idiab21, islash
   integer :: idiab4, idiab41, izp0, ize0, ir0, izp01, ize01, ir01
   integer :: islim, iacc, imaxit, iweights, lenf, ispa
   integer :: nout, kgmin, iwavefun, i, j
   integer :: iprin, imcorr, ifactr, ipgtol

   real(8) :: zp0, zps, ze0, zes, r0
   real(8) :: zpmin, zemin, z1min, z2min, rmin, femin

   real(8), dimension(3)   :: grad
   real(8), dimension(3,3) :: hess

   real(8) :: dzpmin, dzemin, d2zpmin, d2zemin, d2zpzemin

   adiab = .false.
   diab2 = .false.
   diab4 = .false.
   gfix  = .false.

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Extract options
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ikey = index(keywrd,' MIN3(')

   if (ikey.eq.0) then
      write(*,'(/1x,"*** (in MINIMA3): You MUST specify options for MIN3 keyword ***"/)')
      stop
   else
      call getopt(keywrd,ikey+6,options)
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Define a surface type (ADIAB, DIAB2, DIAB4) and
   ! determine a particular state to search a minimum on
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (index(options,' ADIAB').ne.0) then

      mode = 'ADIAB'
      adiab=.true.

      iadiab = index(options,' ADIAB=')
      if (iadiab.ne.0) then
         istate = reada(options,iadiab+7)
      else
         istate = 1
      endif
      iset = 1

   elseif (index(options,' DIAB2').ne.0) then

      mode = 'DIAB2'
      diab2 = .true.

      idiab2 = index(options,' DIAB2=')

      if (idiab2.ne.0) then

         idiab21 = idiab2+7
         islash = index(options(idiab21:),'/')
         istate = reada(options(idiab21:idiab21+islash-2),1)
         iset = reada(options,idiab21+islash)

      else

         istate = 1
         iset = 1

      endif

   elseif (index(options,' DIAB4').ne.0) then

      mode = 'DIAB4'
      diab4=.true.

      idiab4 = index(options,' DIAB4=')

      if (idiab4.ne.0) then

         idiab41 = idiab4+7
         islash = index(options(idiab41:),'/')
         istate = reada(options(idiab41:idiab41+islash-2),1)
         iset = reada(options,idiab41+islash)

      else

         istate = 1
         iset = 1

      endif

   ELSE

      !~~~~~~~~~~~~~~~
      ! Default values
      !~~~~~~~~~~~~~~~
      mode = 'ADIAB'
      adiab=.true.
      istate = 1
      iset = 1

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Print out the information
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   write(6,'(/1x,"The minimum will be located on the ",1x,i2,a3,1x,a5," free energy surface",/,&
             &1x,"within the ",i1,a3," set"/)') istate,th(istate),mode,iset,th(iset)

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Extract initial values of solvent coordinates
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   izp0 = index(options,' ZP0=')
   ize0 = index(options,' ZE0=')
   ir0  = index(options,' R0=')

   if (izp0.ne.0.and.ize0.ne.0) then

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! initial value and scaling factor for zp coordinate
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      izp01 = izp0 + 5
      islash = index(options(izp01:),'/')
      zp0 = reada(options(izp01:izp01+islash-2),1)
      zps = reada(options,izp01+islash)
      if (iminim.eq.2) zps = 1.d0

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! initial value and scaling factor for ze coordinate
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ize01 = ize0 + 5
      islash = index(options(ize01:),'/')
      ze0 = reada(options(ize01:ize01+islash-2),1)
      zes = reada(options,ize01+islash)
      if (iminim.eq.2) zes = 1.d0

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! initial value for the gating coordinate
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ir01 = ir0 + 4
      r0 = reada(options(ir01:),1)
      if (ir0.eq.0.or.r0.eq.0.d0) gfix = .true.

      write(6,'(/1x,"initial values for solvent coordinates:")')
      write(6,'( 1x,"zp0=",f10.3," kcal/mol  ;  scaling factor: ",f10.3,/,&
                &1x,"ze0=",f10.3," kcal/mol  ;  scaling factor: ",f10.3,/)')&
                &zp0,zps,ze0,zes

      if (.not.gfix) then
         if (r0.lt.glist(1)*bohr2a.or.r0.gt.glist(npntsg)*bohr2a) then
            write(*,'(/1x,"*** (in minima3): r0 value is out of bounds ***"/)')
            stop
         endif
         write(6,'(/1x,"initial values for the gating coordinate:")')
         write(6,'( 1x,"R0=",f10.3," a "/)') r0
      else
         write(6,'(/1x,"Calculation will be done for a single gating distance (from the input file)")')
      endif

   else

      write(*,'(/1x,"*** (in minima3): ",/,&
      &"you must specify at least zp0=, ze0= options for the min3 keyword ***"/)')
      stop

   endif


   if (iminim.eq.1) then

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Maximum length of the Newton-Raphson step
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      islim = index(options,' SLIM=')
      if (islim.ne.0) then
         slim = reada(options,islim+6)
      else
         slim = 1.d-1
      endif

      write(6,'(/1x,"Maximum length of the newton-raphson step: ",g13.6)') slim

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Accuracy of the minimization
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      iacc = index(options,' ACC=')
      if (iacc.ne.0) then
         acc = reada(options,iacc+5)
      else
         acc = 1.d-6
      endif

      write(6,'(/1x,"Accuracy of the minimization: ",g13.6)') acc

   elseif (iminim.eq.2) then

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Printing level in LBFGS routine
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      iprin = index(options,' IPRINT=')
      if (iprin.ne.0) then
         iprint = reada(options,iprin+8)
      else
         iprint = -1
      endif
      write(6,'(/1x,"Printing level in LBFGS minimization: ",i10)') iprint

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Printing level in LBFGS routine
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      imcorr = index(options,' MCORR=')
      if (imcorr.ne.0) then
         mcorr = reada(options,imcorr+7)
      else
         mcorr = 5
      endif
      write(6,'(/1x,"Number of corrections used in the limited memory matrix (LBFGS): ",i10)') mcorr

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! LBFGS tolerance
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ifactr = index(options,' FACTR=')
      if (ifactr.ne.0) then
         factr = reada(options,ifactr+7)
      else
         factr = 1.d+6
      endif
      write(6,'(/1X,"LBFGS tolerance: ",G13.6)') factr

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! LBFGS projected gradient tolerance
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ipgtol = index(options,' PGTOL=')
      if (ipgtol.ne.0) then
         pgtol = reada(options,ipgtol+7)
      else
         pgtol = 1.d-6
      endif
      write(6,'(/1X,"LBFGS gradient tolerance: ",G13.6)') pgtol

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Maximum number of iterations
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   imaxit = index(options,' MAXIT=')
   if (imaxit.ne.0) then
      maxit = reada(options,imaxit+7)
   else
      maxit = 100
   endif

   write(6,'(/1x,"Maximum number of iterations: ",i10)') maxit

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Warning...
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   if (method.eq.2) then
      write(6,'(/1x,"***warning: using method=2 in minimum search is not effective",/,&
                &1x,"the calculation might take a very long time.",/,&
                &1x,"So, be prepared..."/)')
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Call minimization routine
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (gfix) then
      call fes3gmin(mode,iset,istate,1,zp0,zps,ze0,zes,&
                       ierr,zpmin,zemin,femin,&
                       dzpmin,dzemin,d2zpmin,d2zemin,d2zpzemin)
   else
      call fes3min(mode,iset,istate,zp0,zps,ze0,zes,r0,&
                  &ierr,zpmin,zemin,kgmin,rmin,femin,grad,hess)
   endif


   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Print out the results of the minimization
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   write(6,'(/1x,"RESULTS OF THE MINIMIZATION (kcal/mol):"/)')
   write(6,'( 1x,"Solvent coordinates:             ",2f13.6)') zpmin,zemin
   call zpze_to_z1z2(zpmin,zemin,z1min,z2min)
   write(6,'( 1x,"Transformed solvent coordinates: ",2f13.6)') z1min,z2min
   if (.not.gfix) write(6,'( 1x,"Gating grid point: ",  i3,f13.6)') kgmin,rmin
   write(6,'( 1x,"Free energy at the min: ",f13.6)') femin
   if (gfix) then
      write(6,'( 1x,"Gradient: ",             2g20.10)') dzpmin,dzemin
      write(6,'( 1x,"Hessian: ",/,           (2f13.6))') d2zpmin, d2zpzemin, d2zpzemin, d2zemin
   else
      write(6,'( 1x,"Gradient: ",             3g20.10)') grad
      write(6,'( 1x,"Hessian: ",/,           (3f13.6))') ((hess(i,j),j=1,3),i=1,3)
   endif
   write(6,'(/)')

   if (ierr.ne.0) then
      write(6,'(/1x,"Minimization failed: the further output is meaningless..."/)')
      return
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Some additional output
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   weights = index(options,' WEIGHTS').ne.0
   wavefun = index(options,' WAVEFUN').ne.0

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Print out the results of the wavefunction analysis
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (weights) then

      iweights = index(options,' WEIGHTS=')

      if (iweights.eq.0) then
         fname  = 'weights3'
         lenf = 8
      else
         ispa = index(options(iweights+9:),' ')
         fname = options(iweights+9:iweights+ispa+7)
         lenf = ispa - 1
         call locase(fname,lenf)
         write(6,'(/1x,"Results of the wavefunction analysis at the minimum are written",/,&
         &" to the external file <",a,">")') fname(1:lenf)
      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Calculate and print out the weights
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      nout = 1
      lenf = lenf + ljob + 1
      fname = job(1:ljob)//'/'//fname(1:lenf)
      open(1,file=fname(1:lenf),status='new')
      call weight3(nout,kgmin,zpmin,zemin,mode,iset,istate)
      close(nout)

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Print out the proton vibrational wavefunctions
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (wavefun) then

      iwavefun = index(options,' WAVEFUN=')

      if (iwavefun.eq.0) then

         fname  = 'vibfun'
         lenf = 6

      else

         ispa = index(options(iwavefun+9:),' ')
         fname = options(iwavefun+9:iwavefun+ispa+7)
         lenf = ispa - 1
         call locase(fname,lenf)

      endif

      fname  = job(1:ljob)//'/'//fname(1:lenf)//'.prot'
      lenf = lenf + 5 + ljob + 1
      nout  = 1
      open(1,file=fname(1:lenf),status='new')
      write(6,'(/1x,"vibrational wavefunctions at the minimum are written",/,&
                  &" to the external files <",a,">")') fname(1:lenf)//'.*'

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! print out the wavefunctions
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call wfn3prt(nout,zpmin,zemin,mode,iset,5)
      close(nout)

   endif

   return

end subroutine minima3

