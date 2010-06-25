subroutine minima2
!====================================================================
!  Sets parameters for the routine which finds minima
!  on the spesified two-dimensional free energy surface
!  using canonical Newton-Raphson minimization algorithm.
!
!  KEYWORD: MIN2
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
!  ACC=<VALUE> - the desired accuracy of the minimization
!
!  SLIM=<VALUE> - the maximum length of the Newton-Raphson step
!
!  MAXIT=N  - maximum number of Newton-Raphson steps allowed
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
!  souda
!  2010/06/25 20:02:35
!  4.1
!  Exp
!  minima2.f90,v 4.1 2010/06/25 20:02:35 souda Exp
!  minima2.f90,v
!  Revision 4.1  2010/06/25 20:02:35  souda
!  Release 4.1
!
!  Revision 1.5  2008/04/11 00:07:19  souda
!  length of string OPTIONS increased to 1024
!  to accomodate more options (not critical)
!
!  Revision 1.4  2007/03/12 23:08:03  souda
!  Modifications related to using LBFGS minimization method.
!
!  Revision 1.3  2004/05/15 03:32:45  souda
!  Added Borgis-Hynes rate routine
!
!  Revision 1.2  2004/01/15 20:09:53  souda
!  minor bugs
!
!  Revision 1.1.1.1  2004/01/13 19:32:44  souda
!  Initial PCET-4.0 Release
!
!
!====================================================================
   use pardim
   use strings
   use keys
   use cst
   use control
   use fesmin_2d

   implicit none
   character(1024) :: options
   character( 40)  :: fname, fnamep, fnameg
   character(  5)  :: mode
   logical         :: adiab, diab2, diab4, weights, wavefun

   integer :: ierr, ikey, iadiab, istate, iset, idiab2, idiab21
   integer :: islash, idiab4, idiab41, izp0, ize0, izp01, ize01
   integer :: islim, iacc, imaxit, iweights, ispa, lenf, nout
   integer :: iwavefun, noutg
   integer :: iprin, imcorr, ifactr, ipgtol
   
   real(8) :: zp0, zps, ze0, zes
   real(8) :: zpmin, zemin, femin, dzpmin, dzemin, d2zpmin, d2zemin, d2zpzemin

   adiab = .false.
   diab2 = .false.
   diab4 = .false.

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Extract options
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ikey = index(keywrd,' MIN2(')

   if (ikey.eq.0) then
      write(*,'(/1x,"*** (in minima2): you must specify options for min2 keyword ***"/)')
      stop 'try again (minima2)...'
   else
      call getopt(keywrd,ikey+6,options)
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Define a surface type (ADIAB, DIAB2, DIAB4) and
   ! determine a particular state to search a minimum on
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (index(options,' ADIAB').ne.0) then

      mode = 'ADIAB'
      adiab = .true.

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

         idiab21 = idiab2 + 7
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

   else

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
   write(6,'(/1X,"The minimum will be located on the ",&
             &1X,I2,A3,1X,A5," free energy surface",/,&
             &1X,"within the ",I1,A3," set"/)')&
             &istate,th(istate),mode,iset,th(iset)

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Extract initial values of solvent coordinates
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   izp0 = index(options,' ZP0=')
   ize0 = index(options,' ZE0=')

   if (izp0.ne.0.and.ize0.ne.0) then

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! initial value and scaling factor for zp coordinate
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      izp01 = izp0+5
      islash = index(options(izp01:),'/')
      zp0 = reada(options(izp01:izp01+islash-2),1)
      zps = reada(options,izp01+islash)
      if (iminim.eq.2) zps = 1.d0

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! initial value and scaling factor for ze coordinate
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ize01 = ize0+5
      islash = index(options(ize01:),'/')
      ze0 = reada(options(ize01:ize01+islash-2),1)
      zes = reada(options,ize01+islash)
      if (iminim.eq.2) zes = 1.d0

      write(6,'(/1x,"initial values for solvent coordinates:")')
      write(6,'( 1x,"zp0=",f8.3," kcal/mol  ;  scaling factor: ",f8.3,/,&
                &1x,"ze0=",f8.3," kcal/mol  ;  scaling factor: ",f8.3,/)')&
                &zp0,zps,ze0,zes

   else

      write(*,'(/1x,"*** (in minima2): ",/,&
      &"you must specify both zp0= and ze0= options for the min2 keyword ***"/)')
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

      write(6,'(/1x,"Maximum length of the Newton-Raphson step: ",g12.6)') slim

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Accuracy of the minimization
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IACC = INDEX(OPTIONS,' ACC=')
      IF (IACC.NE.0) THEN
         ACC = READA(OPTIONS,IACC+5)
      ELSE
         ACC = 1.D-6
      ENDIF
      WRITE(6,'(/1X,"Accuracy of the minimization: ",G12.6)') ACC

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
      write(6,'(/1X,"LBFGS tolerance: ",G12.6)') factr

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! LBFGS projected gradient tolerance
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ipgtol = index(options,' PGTOL=')
      if (ipgtol.ne.0) then
         pgtol = reada(options,ipgtol+7)
      else
         pgtol = 1.d-6
      endif
      write(6,'(/1X,"LBFGS gradient tolerance: ",G12.6)') pgtol

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

   call fes2min(mode,iset,istate,zp0,zps,ze0,zes,&
               &ierr,zpmin,zemin,femin,&
               &dzpmin,dzemin,d2zpmin,d2zemin,d2zpzemin)

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Print out the results of the minimization
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   write(6,'(/1x,"RESULTS OF THE MINIMIZATION (kcal/mol):"/)')
   write(6,'(1x,"Coordinates at the min: ",2f12.6)') zpmin,zemin
   write(6,'(1x,"Free energy at the min: ", f12.6)') femin
   write(6,'(1x,"Gradient: ",              2f12.6)') dzpmin,dzemin
   write(6,'(1x,"Hessian (11,22,12): ",    3f12.6)') d2zpmin,d2zemin,d2zpzemin
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
         nout = 6
      else
         ispa = index(options(iweights+9:),' ')
         fname = options(iweights+9:iweights+ispa+7)
         lenf = ispa - 1
         call locase(fname,lenf)
         fname = job(1:ljob)//'/'//fname(1:lenf)
         lenf = lenf + ljob + 1
         nout = 1
         open(1,file=fname(1:lenf),status='new')
         write(6,'(/1x,"Results of the wavefunction analysis at the minimum are written",/,&
         &" to the external file <",a,">")') fname(1:lenf)
      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Calculate and print out the weightings
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call weight2(nout,zpmin,zemin,mode,iset,istate)
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

      fnamep = job(1:ljob)//'/'//fname(1:lenf)//'.prot'
      fnameg = job(1:ljob)//'/'//fname(1:lenf)//'.gate'
      lenf = lenf + 5 + ljob + 1
      nout  = 1
      noutg = 2
      open(1,file=fnamep(1:lenf),status='new')
      open(2,file=fnameg(1:lenf),status='new')
      write(6,'(/1x,"vibrational wavefunctions at the minimum are written",/,&
      &" to the external files <",a,"> and <",a,">")') fnamep(1:lenf),fnameg(1:lenf)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Print out the wavefunctions
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call wfn2prt(nout,noutg,zpmin,zemin,mode,iset,10,5)
      close(nout)
      close(noutg)

   endif

   return

end subroutine minima2

