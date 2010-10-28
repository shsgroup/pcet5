subroutine wavef2

!=====================================================================
!
!   Calculates and writes out vibrational wavefunctions
!   and potential energy profiles at specified values
!   of solvent coordinates
!
!   KEYWORD: WAVEFUN2
!
!   OPTIONS:
!
!   ZP=<VALUE> - value of the solvent coordinate ZP
!   ZE=<VALUE> - value of the solvent coordinate ZE
!
!   ADIAB - adiabatic electron/proton vibrational free energy surfaces
!
!   DIAB2=<ISET> - ET diabatic free energy surfaces (two output files)
!
!   DIAB4=<ISET> - diabatic free energy surfaces (four output files)
!
!          ISET - set of states (1 for ADIAB, 1/2 for DIAB2)
!
!   NVIBST - number of proton vibrational states to print
!
!   NGATST - number of gating vibrational states to print
!
!   OUTPUT=<filename> - name of the external output file
!
!---------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:37 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
!
!=====================================================================

   use pardim
   use keys
   use strings

   implicit none
   character(1024) :: options
   character(  40) :: fname, fnamep, fnameg
   character(   5) :: mode

   integer :: ikey, iset, idiab2, idiab4, izp, ize
   integer :: invibst, nvibst, ingatst, ngatst
   integer :: ioutput, lenf, ispa, noutp, noutg
   
   real(8) :: zp, ze

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Extract options
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ikey = index(keywrd,' WAVEFUN2(')

   if (ikey.eq.0) then
      write(*,'(/1x,''*** (in WAVEF2): You MUST specify options for WAVEF keyword ***''/)')
      stop
   else
      call getopt(keywrd,ikey+10,options)
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Define a surface type (ADIAB, DIAB2, DIAB4) and
   ! determine a particular state to search a minimum for
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (index(options,' ADIAB').ne.0) then

      mode = 'ADIAB'
      iset = 1

   elseif (index(options,' DIAB2').ne.0) then

      mode = 'DIAB2'

      idiab2 = index(options,' DIAB2=')

      if (idiab2.ne.0) then
         iset = reada(options,idiab2+7)
      else
         iset = 1
      endif

   elseif (index(options,' DIAB4').ne.0) then

      mode = 'DIAB4'

      idiab4 = index(options,' DIAB4=')

      if (idiab4.ne.0) then
         iset = reada(options,idiab4+7)
      else
         iset = 1
      endif

   else

      !~~~~~~~~~~~~~~~
      ! Default values
      !~~~~~~~~~~~~~~~
      mode = 'ADIAB'
      iset = 1

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Print out the information
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    write(6,'(/1x,''The output of vibrational wavefunctions for'',1x,a5,'' states'',/,&
              &1x,''within the '',i1,a3,'' set''/)') mode,iset,th(iset)

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Extract values of solvent coordinates
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   izp = index(options,' ZP=')
   ize = index(options,' ZE=')

   if (izp.ne.0.and.ize.ne.0) then

      zp = reada(options,izp+4)
      ze = reada(options,ize+4)

      write(6,'(/1x,''Values for solvent coordinates:'')')
      write(6,'( 1x,''zp='',f8.3,'' kcal/mol'',1x,''ze='',f8.3,'' kcal/mol''/)') zp,ze

   else

      write(*,'(/1x,''*** (in WAVEF2): '',/,&
      &''You MUST specify both ZP= and ZE= options for the WAVEFUN keyword ***''/)')
      stop

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Number of vibrational wavefunctions to print out
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   invibst = index(options,' NVIBST=')
   if (invibst.ne.0) then
      nvibst = reada(options,invibst+8)
   else
      nvibst = 5
   endif

   write(6,'(/1x,''Number of vibrational wavefunctions to write out: '',i10)') nvibst

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Number of gating vibrational wavefunctions to print out
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ingatst = index(options,' NGATST=')
   if (ingatst.ne.0) then
      ngatst = reada(options,ingatst+8)
   else
      ngatst = 5
   endif

   write(6,'(/1x,''Number of gating wavefunctions to write out: '',i10)') ngatst

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Print out the proton vibrational wavefunctions
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ioutput = index(options,' OUTPUT=')

   if (ioutput.eq.0) then

      fname  = 'vibfun'
      lenf = 6

   else

      ispa = index(options(ioutput+8:),' ')
      fname = options(ioutput+8:ioutput+ispa+6)
      lenf = ispa - 1
      call locase(fname,lenf)

   endif

   fnamep = job(1:ljob)//'/'//fname(1:lenf)//'.prot'
   fnameg = job(1:ljob)//'/'//fname(1:lenf)//'.gate'
   lenf = lenf + 5 + ljob + 1
   noutp = 1
   noutg = 2
   open(noutp,file=fnamep(1:lenf),status='new')
   open(noutg,file=fnameg(1:lenf),status='new')
   write(6,'(/1x,''Vibrational wavefunctions are written'',/,&
               &'' to the external files <'',a,''>'')') fname(1:lenf)//'.*'

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Print out the wavefunctions
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   call wfn2prt(noutp,noutg,zp,ze,mode,iset,nvibst,ngatst)
   close(noutp)
   close(noutg)

   return

end subroutine wavef2

