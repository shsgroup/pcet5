subroutine setjob

!====================================================================C
!   Analyses specified keywords and sets options
!   and parameters for the current job
!---------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2012-03-13 22:06:35 $
!  $Revision: 5.8 $
!  $Log: not supported by cvs2svn $
!  Revision 5.7  2011/03/01 23:52:15  souda
!  insignificant change in output
!
!  Revision 5.6  2011/02/25 19:11:25  souda
!  Now using a separate set of dielectric constant for solvent dynamics.
!
!  Revision 5.5  2011/01/04 19:59:14  souda
!  added: now the custom delta and kappa FRCM parameters can be
!         specified even if the solvent name is explicitly given.
!
!  Revision 5.4  2010/12/15 21:24:56  souda
!  various fixes (non-critical)
!
!  Revision 5.3  2010/11/10 21:14:21  souda
!  Last addition/changes. IMPORTANT: Solvation keyword is back to SOLV( for compatibility reasons
!
!  Revision 5.2  2010/10/28 21:29:37  souda
!  First (working and hopefully bug-free) source of PCET 5.x
!
!
!====================================================================C

   use pardim
   use cst
   use keys
   use strings
   use control
   use parsol
   use geosol
   use geometry
   use frcm

   implicit none

   character(1024) :: options
   character(  40) :: fname
   logical :: ok
   logical :: quantum_proton

   integer :: ikey, ingrid, ingrids, npntsi, inprst, ilims, ialim, islash, itset
   integer :: irhmin, npntspow2, iarlim, igquant, ingast, icoulomb
   integer :: idist, ideltag, ivept, ipars, ispa, lenf
   integer :: igeo, ixyz, imajor, iminor, ikappa, idelta, igeom, itread
   integer :: natsol_, natgas_, nhpt, k, ieps0, ieps8, npntssol2

   real(8) :: x, ba, ba2

   !==================================================
   ! Set control parameters
   !==================================================


   !==================================================
   ! Total charge of the solute
   !==================================================

   ikey = index(keywrd,' CHARGE=')
   if (ikey.ne.0) then
      charge = reada(keywrd,ikey+8)
   else
      charge = 0.d0
   endif
   write(6,'(/1x,"Total charge of the solute: ",f13.6)') charge


   !===============================================================
   ! Solvent model
   !===============================================================
   ! FRCM    - Frequency Resolved Cavity Model (Basilevsky, Rostov)
   ! ELLIPSE - Simple electrostatic model with ellipsoidal cavity
   ! TSET    - Reorganization energy matrix is reconstructed from
   !           the values given in input (new in version 5.2)
   !===============================================================

   ikey = index(keywrd,' SOLV(')

   if (ikey.eq.0) then

      write(*,'(/1x,"*** (in SETJOB): You MUST specify the SOLV keyword with options ***"/)')
      stop

   else

      call getopt(keywrd,ikey+6,options)

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Read values of dielectric constants (EPS0,EPS8) at 298.15 K
   ! [Y. Marcus, Ion Solvation, Wiley, NY, 1985, p.137-138]
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (index(options,' H2O').ne.0.or.index(options,' WATER').ne.0) then

      eps0 = 78.39d0
      eps8 = 1.78d0

   elseif (index(options,' CH2CL2').ne.0.or.index(options,' DICHLOROETHANE').ne.0) then

      eps0 = 8.93d0
      eps8 = 2.02d0

   elseif (index(options,' MEOH').ne.0.or.index(options,' METHANOL').ne.0) then

      eps0 = 33.7d0
      eps8 = 5.6d0

   elseif (index(options,' ETOH').ne.0.or.index(options,' ETHANOL').ne.0) then

      eps0 = 24.6d0
      eps8 = 1.85d0

   elseif (index(options,' CH3CN').ne.0) then

      eps0 = 37.5d0
      eps8 = 1.80d0

   elseif (index(options,' DCE').ne.0) then

      eps0 = 10.4d0
      eps8 = 2.08d0

   elseif (index(options,' THF').ne.0) then

      eps0 = 7.58d0
      eps8 = 1.97d0

   elseif (index(options,' NBZ').ne.0) then

      eps0 = 34.82d0
      eps8 = 2.40d0

   elseif (index(options,' DMF').ne.0) then

      eps0 = 36.7d0
      eps8 = 2.04d0

   else

      ieps0 = index(options,' EPS0=')
      ieps8 = index(options,' EPS8=')

      if (ieps0.ne.0) then
         eps0 = reada(options,ieps0+6)
      else
         write(*,'(/1x,"*** (in SETJOB): You MUST specify EPS0 constant in SOLV ***"/)')
         stop
      endif

      if (ieps8.ne.0) then
         eps8 = reada(options,ieps8+6)
      else
         write(*,'(/1x,"*** (in SETJOB): You MUST specify EPS8 constant in SOLV ***"/)')
         stop
      endif

   endif

   !-- initialize inverse Pekar factor f_0
   f0 = four*pi*eps0*eps8/(eps0 - eps8)

   !-- symmetrization options for reorganization energy matrices

   if (index(options,' SYMT').ne.0) then
      symt=.true.
   else
      symt=.false.
   endif

   if (index(options,' SYMPT').ne.0) then
      sympt=.true.
   else
      sympt=.false.
   endif

   if (index(options,' SYMET').ne.0) then
      symet=.true.
   else
      symet=.false.
   endif

   if (index(options,' REDDENS').ne.0) then
      reddens=.true.
   else
      reddens=.false.
   endif

   if (index(options,' NOSYMD').ne.0) then
      nosymd=.true.
   else
      nosymd=.false.
   endif

   write(6,'(/1x,"Static solvent parameters:")')
   write(6,'(10x,"Static dielectric constant:  ",f13.6)') eps0
   write(6,'(10x,"Optical dielectric constant: ",f13.6)') eps8
   write(6,'(10x,"Inverse Pekar factor f0:     ",f13.6)') f0
   write(6,'(10x,"Total charge of the solute:  ",f13.6)') charge

   if (symt) then
      write(6,'("Solvation T matrices symmetrized for PCET in a symmetric system")')
   endif

   if (sympt) then
      write(6,'("Solvation T matrices symmetrized for PT in a symmetric system")')
   endif

   if (symet) then
      write(6,'("Solvation T matrices symmetrized for ET in a symmetric system")')
   endif

   if (reddens) then
      write(6,'("Reduced basis used for solvation matrices")')
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Solvation model
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (index(options,' TSET').ne.0) then

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Reconstruction of reorganization energy matrix from
      ! input values of partial reorganization energies
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(6,'(/1x,"Reorganization energy matrices will be reconstructed from partial reorganization energies specified in input")')
      isolv = 0


   elseif (index(options,' ELLIPSE').ne.0.or.index(options,' ELCM').ne.0) THEN

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Simple electrostatic model with ellipsoidal cavity
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(6,'(/1x,"Ellipsoidal model will be used for solvation calculations")')
      isolv = 1

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Read values for ellipsoid parameters
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! A - major semiaxis
      ! B - minor semiaxis
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      imajor = index(options,' A=')
      if (imajor.ne.0) then
         a = reada(options,imajor+3)
      else
         write(*,'(/1x,"*** (in SETJOB): You MUST specify the major axis of ellipsoid ***"/)')
         stop
      endif

      iminor = index(options,' B=')
      if (iminor.ne.0) then
         b = reada(options,iminor+3)
      else
         write(*,'(/1x,"*** (in SETJOB): You MUST specify the minor axis of ellipsoid ***"/)')
         stop
      endif

      ba  = b/a
      ba2 = ba*ba
      l0  = dsqrt(1.d0/(1.d0-ba2))
      r   = 2.d0*a/l0

      write(6,'(/1x,"Ellipsoidal cavity parameters:")')
      write(6,'(10x,"Major semiaxis:      ",f13.6)') a
      write(6,'(10x,"Minor semiaxis:      ",f13.6)') b
      write(6,'(10x,"Interfocal distance: ",f13.6)') r

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   elseif (index(options,' FRCM').ne.0) then

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Frequency Resolved Cavity Model (Basilevsky, Chudinov, Rostov)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(6,'(/1x,"FRCM model will be used for solvation calculations")')
      isolv = 2

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Read values for cavity parameters in FRCM
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! KAPPA - factor for VdW radii
      ! DELTA - the width of the layer between two cavities
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      !-- default values for different solvents

      if (index(options,' H2O').ne.0.or.index(options,' WATER').ne.0) then
         kappa = 0.9d0
         delta = 1.0d0
      elseif (index(options,' CH2CL2').ne.0) then
         kappa = 0.9d0
         delta = 1.2d0
      elseif (index(options,' MEOH').ne.0) then
         kappa = 0.9d0
         delta = 1.1d0
      elseif (index(options,' ETOH').ne.0) then
         kappa = 0.9d0
         delta = 1.4d0
      elseif (index(options,' CH3CN').ne.0) then
         kappa = 0.9d0
         delta = 1.8d0
      elseif (index(options,' DCE').ne.0) then
         kappa = 0.9d0
         delta = 2.3d0
      elseif (index(options,' THF').ne.0) then
         kappa = 0.9d0
         delta = 2.3d0
      elseif (index(options,' NBZ').ne.0) then
         kappa = 0.9d0
         delta = 3.7d0
      elseif (index(options,' DMF').ne.0) then
         kappa = 0.9d0
         delta = 2.7d0
      else
         kappa = 0.9d0
         delta = 1.0d0
      endif

      ikappa = index(options,' KAPPA=')
      idelta = index(options,' DELTA=')

      if (ikappa.ne.0) then
         kappa = reada(options,ikappa+7)
         write(*,'(/1x,"The custom value of KAPPA parameter was specified in the input"/)')
      endif

      if (idelta.ne.0) then
         delta = reada(options,idelta+7)
         write(*,'(/1x,"The custom value of DELTA parameter was specified in the input"/)')
      endif

      write(6,'(/1x,"FRCM cavity parameters:")')
      write(6,'(10x,"Factor for VdW radii (KAPPA):     ",f13.6)') kappa
      write(6,'(10x,"Width of the inner layer (DELTA): ",f13.6," Angstroms")') delta

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


   else

      write(*,'(/1x,"*** (in SETJOB): You MUST specify the solvation model (TSET/ELLIPSE/ELCM/FRCM) with SOLV ***"/)')
      stop

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Read external file with the geometry of the complex
   ! for solvation calculations and charges for EVB states
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   igeom = index(options,' GEOM=')

   if (igeom.ne.0) then

      !-- extract the name of the geometry input file
      ispa = index(options(igeom+6:),' ')
      fname = options(igeom+6:igeom+ispa+4)
      lenf = ispa - 1
      call locase(fname,lenf)
      write(6,'(/1x,"Geometry for solvation calculations from the file <",a,">")') fname(1:lenf)

      !-- copy the input file to the output directory
      call system("cp "//fname(1:lenf)//" "//job(1:ljob))

      !-- open the file
      open(1,file=fname(1:lenf),status='old')

      !-- scan input file for number of atoms
      call scan_geo(1,natsol_)

      !-- allocate the geometry/charges arrays
      call alloc_geosol(natsol_)

      !-- read the geometry
      rewind 1
      call geoin0(1,natsol_,natsol,labsol,xyzsol,chrsol)

      if (isolv.eq.2) then

         rewind 1
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! In case of FRCM solvation model:
         ! - read specific FRCM keywords
         ! - initialize internal COMMON-blocks used in FRCM
         !   calculations of the reorganization energy matrices
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         call frcminit(1,eps0,eps8,kappa,delta,charge,natsol,labsol,xyzsol,chrsol)

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! - use specific FRCM keywords
         ! - build spheres around the atoms
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! call getsfe(1,xyzsol,natsol)
         ! now it is called from within FRCMINIT

      endif
      close(1)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! For ellipsoidal model check whether all the point
      ! charges are on the straight line (x-axis)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (isolv.eq.1) then

         ok = straight(natsol,xyzsol)
         if (.not.ok) then
            write(*,'(/1x,"*** (in SETJOB): ",/,&
            &"The geometry for solvation calculations using an ellipsoidal model",/,&
            &" should be LINEAR !"/)')
            stop
         endif

      endif

   else

      write(*,'(/1x,"*** (in SETJOB): You MUST specify file with the geometry (GEOM=) in SOLV keyword ***"/)')
      stop

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Write out the cartesian coordinates to the external
   ! file if option XYZOUT is specified
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ixyz = index(options,' XYZOUT=')

   if (ixyz.ne.0.and.isolv.gt.0) then

      ispa = index(options(ixyz+8:),' ')
      fname = options(ixyz+8:ixyz+ispa+6)
      lenf = ispa - 1
      call locase(fname,lenf)
      write(6,'(/1x,"Geometry for solvation calculations is written to the file <",A,">",/,&
      &"in xyz format (use MOLDEN to view)")') fname(1:lenf)
      open(1,file=job(1:ljob)//'/'//fname(1:lenf),status='new')
      call xyzout(1,natsol,labsol,xyzsol)
      close(1)

   endif

   return

CONTAINS

!=======================================================================
!  Internal procedures
!=======================================================================
   
   subroutine scan_geo(iunit,n)

   ! scans geometry input file (iunit) and extracts
   ! number of atoms (n)

      implicit none
      integer, intent(in)  :: iunit
      integer, intent(out) :: n

      character(80) :: line
      integer :: io_status

      ! title records
      read(iunit,'(a)') line
      read(iunit,'(a)') line
      read(iunit,'(a)') line

      n = 0

      do
         read(iunit,'(a)',iostat=io_status) line
         if (io_status.ne.0) then
            write(*,*) "Error/eof while scanning geometry input file"
            close(iunit)
            stop
         endif
         if (line.eq. ' ') exit
         n = n + 1
      enddo

      if (n.le.0) then
         write(*,*) "No atoms in the geometry input file"
         close(iunit)
         stop
      endif

      return

   end subroutine scan_geo


   !=======================================================================
   logical function straight(nat,xyz)
   ! Checks whether all the atoms are on the straight line (x-axis)

      implicit none

      integer, intent(in) :: nat
      real*8,  intent(in), dimension(:,:) :: xyz

      integer :: i, j

      straight = .true.

      do i=1,nat
         do j=2,3
            if (dabs(xyz(j,i)).gt.1.d-7) then
               straight = .false.
               return
            endif
         enddo
      enddo

      return
   end function straight

   !=======================================================================
   logical function power2(n)
   ! Checks whether the given integer (N) is a power of two

      implicit none
      integer, intent(in) :: n

      integer, parameter :: itwo = 2
      integer :: nw

      nw = n
      if (nw.lt.2) then
         power2 = .false.
         return
      else
         power2 = .true.
      endif

      do while (nw.gt.1)
         if (mod(nw,itwo).ne.0) then
            power2 = .false.
            return
         else
            nw = nw/itwo
         endif
      enddo

      return
   end function power2

end subroutine setjob
