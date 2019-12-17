module geometry

!-----------------------------------------------------------------------
! Routines related to the geometry manipulations
!-----------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2012-03-13 22:05:20 $
!  $Revision: 5.5 $
!  $Log: not supported by cvs2svn $
!  Revision 5.4  2010/12/15 21:24:55  souda
!  various fixes (non-critical)
!
!  Revision 5.3  2010/11/24 22:37:40  souda
!  workaround for optional input of custom VdW radii for FRCM calculation
!
!  Revision 5.2  2010/10/28 21:29:36  souda
!  First (working and hopefully bug-free) source of PCET 5.x
!
!
!-----------------------------------------------------------------------

   use pardim
   use cst
   use strings
   use keys
   use elmnts
   use control

   !---------------------------------------------------------------------
   implicit none
   public

!---------------------------------------------------------------------
contains

   !---------------------------------------------------------------------
   subroutine read0(in)
   !---------------------------------------------------------------------
      implicit none
      integer, intent(in) :: in

      character( 1 ) :: keychar
      character(480) :: line
      character(481) :: keywr0
      character(2048) :: keywrd_tmp
      integer        :: length, icount, i, io_status

      length = 1
      keywrd_tmp(1:1) = space

      !-- read keywords section until either an empty line
      !   or end of file occurs

      !-- read lines with titles (length < 120)

      read(in,'(a)') title
      read(in,'(a)') title2

      !-- read a line with keywords (length < 480)
      !   empty line or end of file terminates the loop

      do
         read(in,'(a)',iostat=io_status) line
         if (line.eq.' '.or.&
            &io_status.lt.0.or.&
            &index(line,'END').ne.0.or.&
            &index(line,'end').ne.0.or.&
            &index(line,'End').ne.0.or.&
            &index(line,'EnD').ne.0) exit

         !-- remove double spaces and forbidden symbols

         icount = len(line)
         call rmdbsp(line,icount)

         if (line(1:1).ne.space) then
            keywr0 = space//line
            icount = icount + 1
         else
            keywr0 = line
         endif

         keywrd_tmp = keywrd_tmp(1:length)//keywr0(1:icount)
         length = length + icount

      enddo

      keywrd = keywrd_tmp(2:length)

      do i=1024,1,-1
         keychar = keywrd(i:i)
         if (index(allow1,keychar).ne.0.or.index(allow2,keychar).ne.0) then
            length = i
            exit
         endif
      enddo

      !-- convert keywrd string to uppercase

      call upcase(keywrd,length)

      return

   end subroutine read0


   !---------------------------------------------------------------------
   subroutine geoin0(in,nat_,nat,labels_,xyz,chr)
   !---------------------------------------------------------------------
   ! Reads geometry of the reaction complex
   ! and charges for EVB states.
   !---------------------------------------------------------------------
   !  Parameters:
   !
   !  IN  - channel number for the input file
   !
   !  NAT_ - number of atoms including dummy atoms (input)
   !
   !  NAT  - number of atoms excluding dummy atoms (output)
   !
   !  LABELS(NAT) - atomic numbers of atoms.
   !
   !  XYZ(3,NAT)  - cartesian coordinates of the reaction complex.
   !
   !  CHR(4,NAT)  - charges on atoms for the EVB states 1a,1b,2a,and 2b.
   !--------------------------------------------------------------------
   !  The input file is supposed to have the following format:
   !
   !  Card 1: Keyword string specifying the geometry format
   !          XYZ - cartesian coordinates;
   !          INT - internal coordinates (default);
   !          <...> - any FRCM keywords.
   !
   !  Card 2: Title (up to 80 characters)
   !
   !  Card 3: Comment (up to 80 characters)
   !
   !  Cards 4..N: geometry specification
   !
   !  Card N+1: empty line
   !
   !  Card N+2: optional line with RADIUS keyword (custom VdW radii)
   !
   !  Cards N+3..M: VdW radii for FRCM
   !
   !  Card M+1: empty line
   !
   !  Card M+2: Title for EVB charges (with CHARGES keyword)
   !
   !  Cards M+3..L: charges on atoms for EVB states 1a, 1b, 2a, and 2b
   !
   !  Card L+1: empty line
   !
   !--------------------------------------------------------------------
      implicit none

      integer, intent(in)    :: in
      integer, intent(in)    :: nat_
      integer, intent(inout) :: nat
      integer, intent(inout) :: labels_(nat_)
      real(8), intent(inout) :: xyz(3,nat_), chr(nelst,nat_)

      logical   :: xyzkey, intkey
      character(241) :: key
      character( 80) :: tit, kom
      character( 81) :: line
      integer, dimension(nat_)   :: nai, nbi, nci
      real(8), dimension(3,nat_) :: geo
      real(8), dimension(nat_)   :: qatom

      integer :: length, i, j, natom, nchr, io_status
      real(8) :: tach, dif

      key  = ' '
      line = ' '
      tit  = ' '
      kom  = ' '

      !-- Read a line with keywords (length < 80)
      read(in,'(a)') key(1:80)

      !-- Remove double spaces and forbidden symbols

      if (index(key(1:80),' +') .ne.0) then

         !-- read second keyword line
         read(in,'(a)')key(81:160)

         if (index(key(81:160),' +') .ne.0) then
            ! read third keyword line
            read(in,'(a)') key(161:240)
         endif

      endif

      line=' '

      if (key(1:1).ne.space) then
         key(1:241) = space//key(1:240)
      endif

      length = 241
      call rmdbsp(key,length)
      call upcase(key,length)

      xyzkey = index(key,' XYZ').ne.0
      intkey = index(key,' INT').ne.0

      if (xyzkey.and.intkey) then
         write(*,'(/1x,''*** input error (in geo): '',&
            &''you should specify either xyz or int keyword,''/,&
            &''but not both of them. It is confusing...''/)')
         stop 'try again...'
      endif

      if ((.not.xyzkey).and.(.not.intkey)) then
         write(*,'(/1x,''*** input error (in geo): '',&
            &''you must specify either xyz or int keyword''/)')
         stop 'try again...'
      endif

      !-- Read title and comments cards
      read(in,'(a)') tit
      read(in,'(a)') kom

      !-- Read in geometry specification
      call getcoo(in,labels_,geo,nai,nbi,nci,natom,xyzkey)

      !-- Check data

      do i=1,natom

         if (labels_(i) .le. 0 ) then
            write(6,'('' atomic number of '',i3,'' ?'')') labels_(i)
            stop 'sorry, you have to try again...'
         endif

         if (intkey) then
            if (  nai(i).ge.i.or. nbi(i).ge.i.or. nci(i).ge.i &
           .or. (nai(i).eq.nbi(i))   .and. i.gt.1 &
           .or. (nai(i).eq.nci(i).or.nbi(i).eq.nci(i))  .and. i.gt.2) then
               write(6,'('' atom number '',i3,'' is illdefined'')') i
               stop 'sorry, you have to try again...'
            endif
         endif

      enddo

      if (intkey) then

         !-- convert angles to radians
         do i=1,natom
            do j=2,3
               geo(j,i) = geo(j,i) * 0.01745329252d0
            enddo
         enddo

      endif

      !-- Print out Z-matrix in MOPAC format

      if (intkey) call geoput(geo,nai,nbi,nci,natom,labels_)

      !-- Convert internal coordinates to cartesian ones

      call convrt(geo,xyz,labels_,nai,nbi,nci,natom,xyzkey)
      nat = natom

      !-- skip optional sections read by FRCM routines

      do

         read(in,'(a)',iostat=io_status) line

         if (io_status.gt.0) then
            write(6,'(/" ERROR during reading geometry input file")')
            stop 'in the GEOIN routine...'
         endif

         if (index(line,"CHARGES").ne.0) then
            backspace in
            exit
         endif

         if (io_status.lt.0) then
            write(6,'(/" Unexpected end of the geometry input file: CHARGES section is missing.")')
            stop 'in the GEOIN routine...'
         endif

      enddo


      !-- Read in EVB charges on atoms
      !   (Note that dummy atoms are not in the list!!!)

      call getchr(in,nchr,chr)

      !-- Check whether a number of charges is consistent
      !   with a number of atoms

      if (nchr.ne.nat) then
         write(*,'(/1x,''*** input error (in geo): '',&
            &''number of charges is not equal to a number of atoms''/)')
         stop 'try again...'
      endif

      !-- Check whether a total charge for all EVB states
      !   is equal to the total charge specified in KEYWORDS

      do i=1,nelst

         tach = 0.d0
         do j=1,nchr
            qatom(j) = chr(i,j)
            tach = tach + qatom(j)
         enddo

         dif = dabs(tach - charge)

         if (dif.gt.1.0d-2) then

            write(*,'(/'' CAUTION! Sum of atomic charges ('',F10.3,'')'',&
            &'' in column '',I2,'' of the input file'',&
            &'' is different''/'' from the charge'',&
            &'' specified in keyword ('',F10.3,''). The program'',&
            &'' considers this''/'' difference (''&
            &,F10.3,'' as result of an error in the input file'')')&
            tach,i,charge,dif
            write(6,'('' program execution is terminated'')')
            write(6,'('' check charges in the input file'')')
            stop 'In GEO...'

         elseif (dif.gt.1.d-4) then

            !-- Normalization of charge will be performed

            write(6,'(/'' ATTENTION! Sum of atomic charges ('',F12.5,&
            &'') is slightly different''/&
            &'' from the charge specified in keyword ('',F12.5,'').'')')&
            tach,charge

            call normch(nchr,charge,qatom)

            do j=1,nchr
               chr(i,j) = qatom(j)
            enddo
            write(6,'(" a compensation of atomic charges is performed"&
            &," to gain consistency"/" in charges")')

         endif

      enddo

      return

   end subroutine geoin0


   !---------------------------------------------------------------------
   subroutine geoin(in,dm,am,nat_,nat,ipt,iet,labels_,xyz,chr)
   !---------------------------------------------------------------------
   ! Reads geometry of the reaction complex
   ! and charges for EVB states.
   !---------------------------------------------------------------------
   !  Parameters:
   !
   !  IN  - channel number for the input file
   !
   !  DM  - mass of the donor atom (group)
   !
   !  AM  - mass of the acceptor atom (group)
   !
   !  NAT_ - number of atoms including dummy atoms (input)
   !         number of atoms excluding dummy atoms (output)
   !
   !  IPT(3) - numbers of three atoms identifying the
   !           proton transfer interface; the first
   !           atom defines the proton donor, the second
   !           atom defines the transfering hydrogen, and
   !           the third atom defines the proton acceptor.
   !           Note that PT state "a" always corresponds
   !           to the hydrogen bonded to the proton donor.
   !
   !  IET(2) - numbers of two atoms identifying the
   !           electron transfer interface; the first
   !           atom defines the electron donor, and
   !           the second atom defines the electron acceptor.
   !           Note that ET state "1" always corresponds
   !           to the electron localized on the donor.
   !           (meaningful only for gas phase model geometry)
   !
   !  LABELS(NAT) - atomic numbers of atoms.
   !
   !  XYZ(3,NAT)  - cartesian coordinates of the reaction complex.
   !
   !  CHR(4,NAT)  - charges on atoms for the EVB states 1a,1b,2a,and 2b.
   !--------------------------------------------------------------------
   !  The input file is supposed to have the following format:
   !
   !  Card 1: Keyword string specifying the geometry format
   !          XYZ - cartesian coordinates;
   !          INT - internal coordinates (default);
   !          <...> - any FRCM keywords.
   !
   !  Card 2: Title (up to 80 characters)
   !
   !  Card 3: Comment (up to 80 characters)
   !
   !  Cards 4..N: geometry specification
   !
   !  Card N+1: empty line
   !
   !  Card N+2: optional line with RADIUS keyword (custom VdW radii)
   !
   !  Cards N+3..M: VdW radii for FRCM
   !
   !  Card M+1: empty line
   !
   !  Card M+2: Title for EVB charges (with CHARGES keyword)
   !
   !  Cards M+3..L: charges on atoms for EVB states 1a, 1b, 2a, and 2b
   !
   !  Card L+1: empty line
   !
   !  Card L+2: Title for PT interface (with PT keyword)
   !
   !  Card L+3: three integers identifying the atoms in the PT interface
   !
   !  Card L+4..: empty line
   !
   !  Card L+5: Title for ET interface (with ET keyword)
   !
   !  Card L+6: two integers identifying the ET donor and acceptor atoms
   !
   !--------------------------------------------------------------------
      implicit none

      integer, intent(in)    :: in
      real(8), intent(in)    :: dm, am
      integer, intent(in)    :: nat_
      integer, intent(inout) :: nat
      integer, intent(inout) :: ipt(3), iet(2)
      integer, intent(inout) :: labels_(nat_)
      real(8), intent(inout) :: xyz(3,nat_), chr(nelst,nat_)

      logical   :: xyzkey, intkey
      character(241) :: key
      character( 80) :: tit, kom
      character( 81) :: line
      integer, dimension(nat_)   :: nai, nbi, nci
      real(8), dimension(3,nat_) ::  geo
      real(8), dimension(nat_)   :: qatom

      integer :: length, i, j, natom, nchr, ndp, nda, io_status
      real(8) :: tach, dif, totm, xdp, ydp, zdp, xda, yda, zda
      real(8) :: xcntr, ycntr, zcntr

      key  = ' '
      line = ' '
      tit  = ' '
      kom  = ' '

      ! Read a line with keywords (length < 80)
      read(in,'(a)') key(1:80)

      ! Remove double spaces and forbidden symbols

      if (index(key(1:80),' +') .ne.0) then

         ! read second keyword line
         read(in,'(a)')key(81:160)

         if (index(key(81:160),' +') .ne.0) then
            ! read third keyword line
            read(in,'(a)') key(161:240)
         endif

      endif

      line=' '

      if (key(1:1).ne.space) then
         key(1:241) = space//key(1:240)
      endif

      length = 241
      call rmdbsp(key,length)
      call upcase(key,length)

      xyzkey = index(key,' XYZ').ne.0
      intkey = index(key,' INT').ne.0

      if (xyzkey.and.intkey) then
         write(*,'(/1x,''*** input error (in geo): '',&
            &''you should specify either xyz or int keyword,''/,&
            &''but not both of them. it is confusing...''/)')
         stop 'try again...'
      endif

      if ((.not.xyzkey).and.(.not.intkey)) then
         write(*,'(/1x,''*** input error (in geo): '',&
            &''you must specify either xyz or int keyword''/)')
         stop 'try again...'
      endif

      ! Read title and comments cards
      read(in,'(a)') tit
      read(in,'(a)') kom

      ! Read in geometry specification
      call getcoo(in,labels_,geo,nai,nbi,nci,natom,xyzkey)

      ! Check data

      do i=1,natom

         if (labels_(i) .le. 0 ) then
            write(6,'('' atomic number of '',i3,'' ?'')') labels_(i)
            stop 'sorry, you have to try again...'
         endif

         if (intkey) then
            if (  nai(i).ge.i.or. nbi(i).ge.i.or. nci(i).ge.i &
           .or. (nai(i).eq.nbi(i))   .and. i.gt.1 &
           .or. (nai(i).eq.nci(i).or.nbi(i).eq.nci(i))  .and. i.gt.2) then
               write(6,'('' atom number '',i3,'' is illdefined'')') i
               stop 'sorry, you have to try again...'
            endif
         endif

      enddo

      if (intkey) then

         ! convert angles to radians
         do i=1,natom
            do j=2,3
               geo(j,i) = geo(j,i) * 0.01745329252d0
            enddo
         enddo

      endif

      ! Print out Z-matrix in MOPAC format

      if (intkey) call geoput(geo,nai,nbi,nci,natom,labels_)

      ! Convert internal coordinates to cartesian ones

      call convrt(geo,xyz,labels_,nai,nbi,nci,natom,xyzkey)
      nat = natom

      !-- skip optional sections read by FRCM routines

      do

         read(in,'(a)',iostat=io_status) line

         if (io_status.gt.0) then
            write(6,'(/" ERROR during reading geometry input file")')
            stop 'in the GEOIN routine...'
         endif

         if (index(line,"CHARGES").ne.0) then
            backspace in
            exit
         endif

         if (io_status.lt.0) then
            write(6,'(/" Unexpected end of the geometry input file: CHARGES section is missing.")')
            stop 'in the GEOIN routine...'
         endif

      enddo


      ! Read in EVB charges on atoms
      ! (Note that dummy atoms are not in the list!!!)

      call getchr(in,nchr,chr)

      ! Check whether a number of charges is consistent
      ! with a number of atoms

      if (nchr.ne.nat) then
         write(*,'(/1x,''*** input error (in geo): '',&
            &''number of charges is not equal to a number of atoms''/)')
         stop 'try again...'
      endif

      ! Check whether a total charge for all EVB states
      ! is equal to the total charge specified in KEYWORDS

      do i=1,nelst

         tach=0.d0
         do j=1,nchr
            qatom(j) = chr(i,j)
            tach=tach + qatom(j)
         enddo

         dif = dabs(tach-charge)

         if (dif.gt.1.0d-2) then

            write(*,'(/'' CAUTION! Sum of atomic charges ('',F10.3,'')'',&
            &'' in column '',I2,'' of the input file'',&
            &'' is different''/'' from the charge'',&
            &'' specified in keyword ('',F10.3,''). The program'',&
            &'' considers this''/'' difference (''&
            &,F10.3,'' as result of an error in the input file'')')&
            tach,i,charge,dif
            write(6,'('' program execution is terminated'')')
            write(6,'('' check charges in the input file'')')
            stop 'In GEO...'

         elseif (dif.gt.1.d-4) then

            ! Normalization of charge will be performed

            write(6,'(/'' ATTENTION! Sum of atomic charges ('',F12.5,&
            &'') is slightly different''/&
            &'' from the charge specified in keyword ('',F12.5,'').'')')&
            tach,charge

            call normch(nchr,charge,qatom)

            do j=1,nchr
               chr(i,j)=qatom(j)
            enddo
            write(6,'(" a compensation of atomic charges is performed"&
            &," to gain consistency"/" in charges")')

         endif

      enddo

      ! Read in three integers identifying the atoms
      ! in the PT interface
      ! (Proton donor, Hydrogen, Proton Acceptor)

      read(in,'(a)') line
      read(in,*) (ipt(i),i=1,3)

      !-- check if indices are within atom range
      if (.not.check_range(nat,ipt)) then
         write(6,'(/" The indices of atoms in the PT interface are outside the range of atoms")')
         stop 'in the GEOIN routine...'
      endif

      ! Read in two integers identifying the atoms
      ! in the ET interface
      ! (Electron Donor, Electron Acceptor)

      read(in,'(a)') line
      read(in,*) (iet(i),i=1,2)

      !-- check if indices are within atom range
      if (.not.check_range(nat,iet)) then
         write(6,'(/" The indices of atoms in the ET interface are outside the range of atoms")')
         stop 'in the GEOIN routine...'
      endif

      ! Put the coordinat system origin at the center of mass
      ! of the proton transfer interface

      ndp = ipt(1)
      nda = ipt(3)

      totm = dm + am

      xdp = xyz(1,ndp)
      ydp = xyz(2,ndp)
      zdp = xyz(3,ndp)

      xda = xyz(1,nda)
      yda = xyz(2,nda)
      zda = xyz(3,nda)

      xcntr = (dm*xdp + am*xda)/totm
      ycntr = (dm*ydp + am*yda)/totm
      zcntr = (dm*zdp + am*zda)/totm

      do i=1,nat
         xyz(1,i) = xyz(1,i) - xcntr
         xyz(2,i) = xyz(2,i) - ycntr
         xyz(3,i) = xyz(3,i) - zcntr
      enddo

      !     Rotate the coordinate system to put atoms
      !     DP and DA of the PT interface on X axis.
      !
      !     --------+-----------------*--------+----------->
      !             DP              c.o.m.     DA          x-axis
      !

      call rotate(nat,ndp,nda,xyz)

      return

   end subroutine geoin

   !---------------------------------------------------------------------
   subroutine getcoo(iread,labels_,geo,na,nb,nc,natoms_,xyzkey)
   !---------------------------------------------------------------------
   !     GET reads in the geometry. the element is specified by it's
   !         chemical symbol, or, optionally, by it's atomic number.
   !
   !  on input IREAD  = channel number for read
   !           XYZKEY = flag for cartesian coordinates input
   !
   ! on output LABELS = atomic numbers of all atoms, including dummies.
   !           NATOMS = total number of atoms, including dummies
   !           GEO    = internal coordinates, in angstroms, and degrees.
   !           NA     = integer array of atoms (see data input)
   !           NB     = integer array of atoms (see data input)
   !           NC     = integer array of atoms (see data input)
   !
   !---------------------------------------------------------------------
      implicit none
      
      integer, intent(in)  :: iread
      logical, intent(in)  :: xyzkey
      integer, intent(out) :: natoms_
      real(8), intent(inout), dimension(:,:) :: geo
      integer, intent(inout), dimension(:)   :: labels_, na, nb, nc

      integer, dimension(40) :: istart
      logical :: leadsp, ele_found

      character(120) :: line, string
      character(  2) :: ele, elsymb
      character(  1), parameter :: nine='9', zero='0', lbra='(', rbra=')'

      integer :: numat_, io_status, j, k, i, nvalue, label

      NATOMS_ = 0
      NUMAT_ = 0

      ! read a line
      ! empty line or end of file terminates the loop

      do

         read(iread,'(a)',iostat=io_status) line

         if (line.eq.' '.or.io_status.lt.0) exit

         if (io_status.gt.0) then
            write(6,'( '' error during read at atom number '', i3 )') natoms_
            j = natoms_ - 1
            write(6,'('' data currently read in are '')')
            do k=1,j
               write(6,'(i4,2x,3(f12.5,4x),3(i2,1x))')&
               labels_(k),(geo(j,k),j=1,3),na(k),nb(k),nc(k)
            enddo
            stop 'at the getcoo routine...'
         endif

         if (natoms_.gt.maxatm) then
            write(*,'(/1x,''*** input error (in getcoo): '',/,&
            &''****  max. number of atoms allowed:'',i4)') maxatm
            stop 'sorry, i can not handle it...'
         endif

         natoms_ = natoms_ + 1

         ! clean the input data

         do i=1,120
            if (line(i:i).eq.tab.or.line(i:i).eq.comma) line(i:i)=space
         enddo

         ! initialize istart to interpret blanks as zero's

         do i=1,10
            istart(i)=80
         enddo

         ! find initial digit of all numbers, check for leading spaces
         ! followed by a character and store in ISTART

         leadsp=.true.
         nvalue=0
         do i=1,120
            if (leadsp.and.line(i:i).ne.space) then
              nvalue=nvalue+1
              istart(nvalue)=i
            endif
            leadsp=(line(i:i).eq.space)
         enddo

         ! establish the element's name and isotope,
         ! check for errors or e.o.data

         string = line(istart(1):istart(2)-1)

         if ( string(1:1) .ge. zero .and. string(1:1) .le. nine) then

            ! atomic number used:
            label = reada(string,1)
            if (label.eq.0) then
               natoms_ = natoms_ - 1
               exit
            endif
            
            if (label.lt.0.or.label.gt.108) then
               write(6,'(''  illegal atomic number'')')
               j = natoms_ - 1
               write(6,'('' data currently read in are '')')
               do k=1,j
                  write(6,'(i4,2x,3(f12.5,4x),3(i2,1x))')&
                  labels_(k),(geo(j,k),j=1,3),na(k),nb(k),nc(k)
               enddo
               stop 'at the getcoo routine...'
            endif

         else

            ! atomic symbol used: colon is not necessary
            ele = string(1:2)

            ! check for error in atomic symbol
            ele_found = .false.
            do i=1,108

               if (elsym(i)(3:3).eq.space) then
                  elsymb = elsym(i)(4:4)//space
               else
                  elsymb = elsym(i)(3:4)
               endif

               if (ele.eq.elsymb) then
                  label = i
                  ele_found = .true.
                  exit
               endif

            enddo

            if (.not.ele_found) then
               write(6,'(''  unrecognized element name:  <'',a,''>'')') ele
               j = natoms_ - 1
               write(6,'('' data currently read in are '')')
               do k=1,j
                  write(6,'(i4,2x,3(f12.5,4x),3(i2,1x))')&
                  labels_(k),(geo(j,k),j=1,3),na(k),nb(k),nc(k)
               enddo
               stop 'at the get routine...'
            endif

         endif

         ! All O.K.

         if (label.ne.99.and.label.ne.108) numat_ = numat_ + 1

         labels_(natoms_)   = label
         geo(1,natoms_)    = reada(line,istart(2))

         if (xyzkey) then

            ! cartesian coordinates

            geo(2,natoms_) = reada(line,istart(3))
            geo(3,natoms_) = reada(line,istart(4))
            na(natoms_)    = 0
            nb(natoms_)    = 0
            nc(natoms_)    = 0

         else

            ! internal coordinates

            geo(2,natoms_)  = reada(line,istart(4))
            geo(3,natoms_)  = reada(line,istart(6))
            na(natoms_)     = reada(line,istart(8))
            nb(natoms_)     = reada(line,istart(9))
            nc(natoms_)     = reada(line,istart(10))

         endif

      ! read next line
      enddo

      ! all data read in, clean up and return

      if (natoms_.ge.3) then
         if (na(3).eq.0) then
             nb(3) = 1
             na(3) = 2
         endif
      endif

      return

   end subroutine getcoo

   !---------------------------------------------------------------------
   subroutine getchr(iread,nchr,chr)
   !---------------------------------------------------------------------
   !  GETCHR reads in the charges for EVB states.
   !
   !  on input  IREAD  = input channel number
   !
   !  on output NCHR   = total number of point charges
   !            CHR    = charges for EVB states
   !---------------------------------------------------------------------
      implicit none
      integer, intent(in)  :: iread
      integer, intent(out) :: nchr
      real(8), intent(inout), dimension(:,:) :: chr

      integer, dimension(40) :: istart
      logical                :: leadsp
      character(120)         :: line, tit
      integer                :: io_status, i, nvalue

      NCHR = 0

      ! read title
      read(iread,'(a)',iostat=io_status) tit
      if (io_status.gt.0) then
         write(*,'('' error in getchr: title record '')')
         stop 'at the getchr routine...'
      endif

      ! read a line
      ! empty line or end of file terminates the loop
      do

         read(iread,'(a)',iostat=io_status) line

         if (line.eq.' '.or.io_status.lt.0) exit

         if (io_status.gt.0) then
            write(*,'('' error in getchr: charge number '',i3)') nchr
            stop 'at the getchr routine...'
         endif

         if (nchr.gt.maxatm) then
            write(*,'(/1x,''*** input error (in getchr): '',/,&
            &''****  max. number of charges allowed:'',i4)') maxatm
            stop 'sorry, i can not handle it...'
         endif

         nchr = nchr + 1

         ! clean the input data

         do i=1,120
            if (line(i:i).eq.tab.or.line(i:i).eq.comma) line(i:i) = space
         enddo

         ! initialize istart to interpret blanks as zero's

         do i=1,10
            istart(i) = 80
         enddo

         ! find initial digit of all numbers, check for leading spaces
         ! followed by a character and store in istart

         leadsp = .true.
         nvalue = 0
         do i=1,120
            if (leadsp.and.line(i:i).ne.space) then
              nvalue = nvalue + 1
              istart(nvalue) = i
            endif
            leadsp = (line(i:i).eq.space)
         enddo

         ! assign charges (1a, 1b, 2a, 2b)
         do i=1,nelst
            chr(i,nchr) = reada(line,istart(i))
         enddo

      ! read next line
      enddo

      return

   end subroutine getchr

   !---------------------------------------------------------------------
   subroutine normch(numat_,charge_,qatom)
   !---------------------------------------------------------------------
   ! Mertush's charge compensation.
   !---------------------------------------------------------------------
      implicit none
      integer, intent(in) :: numat_
      real(8), intent(in) :: charge_
      real(8), intent(inout), dimension(:) :: qatom

      integer :: i
      real(8) :: qq2, dc, qqn, qqp, zqn, zqp, aa, charge1

      qq2 = 0.0d0
      do i=1,numat_
         qq2 = qq2 + qatom(i)
      enddo

      write (6,'(/17x,''charge of solute in keyword = '',f12.5)') charge_
      write (6,'(14x,''sum of atomic charges on input = '',f12.5)') qq2

      dc = dabs(qq2-charge_)
      if (dabs(charge_).gt.1.d0) dc = dc/dabs(charge_)

      if (dabs(charge_).lt.1.0d-5) then

         qqn = 0.0d0
         do i=1,numat_
            if (qatom(i).lt.0.0d0) qqn = qqn + qatom(i)
         enddo

         qqp = qq2 - qqn

         if (dabs(qqn).gt.1.d-10.and.dabs(qqp).gt.1.d-10) then

            ! Mertush's charge compensation

            zqn = 1.d0 - 0.5d0*qq2/qqn
            zqp = 1.d0 - 0.5d0*qq2/qqp

            do i=1,numat_
               if (qatom(i).lt.0.d0) then
                   qatom(i) = qatom(i)*zqn
               else
                   qatom(i) = qatom(i)*zqp
               endif
            enddo

         endif

      else

         aa = charge_/qq2
         do i=1,numat_
            qatom(i) = qatom(i)*aa
         enddo

      endif

      charge1 = 0.0d0

      do i=1,numat_
          charge1 = charge1+qatom(i)
      enddo

      write (6,'(4x,''sum of atomic charges after compensation = '',f12.5)') charge1

      return

   end subroutine normch

   !---------------------------------------------------------------------
   subroutine convrt(geo,coord,labels_,na,nb,nc,natoms_,xyzkey)
   !---------------------------------------------------------------------
   !     CONVRT (GMETRY IN AMPAC-PACKAGE)
   !     Computes coordinates from bond-angles and lengths.
   !
   !**** It is adapted from the program written by M.J.S. Dewar.
   !
   !     Two separate options exist within CONVRT. These are:
   !
   ! (A) If XYZKEY=.TRUE. then GEO is assumed to be in cartesian
   !     rather than internal coordinates, and COORD is then set
   !     equal to GEO.
   !
   ! (C) Normal conversion from internal to cartesian coordinates is done
   !
   !  ON INPUT:
   !         GEO    = array of internal coordinates.
   !         NATOMS = number of atoms, including dummies.
   !         NA     = array of atom labels for bond lengths.
   !         NB     = array of atom labels for bond angles.
   !         NC     = array of atom labels for dihedrals.
   !
   !  ON OUTPUT:
   !         COORD    = array of cartesian coordinates
   !         NATOMS = number of atoms, excluding dummies.
   !
   !---------------------------------------------------------------------
      implicit none

      real(8), intent(in), dimension(:,:)    :: geo
      integer, intent(in), dimension(:)      :: na, nb, nc
      logical, intent(in)                    :: xyzkey
      integer, intent(inout), dimension(:)   :: labels_
      integer, intent(inout)                 :: natoms_
      real(8), intent(inout), dimension(:,:) :: coord

      integer :: i, j, mb, mc, ma , k
      real(8) :: ccos, cosa, xb, yb, zb, rbc, xa, ya, za
      real(8) :: xyb, xpa, ypa, xpb
      real(8) :: costh, sinth, sinph, cosph, xqa, zqa, yza
      real(8) :: coskh, sinkh, sina, sind, cosd, xd, yd, zd
      real(8) :: xpd, ypd, zpd, xqd, yqd, zqd, xrd

      ! Option (A)

      if (xyzkey) then
         do i=1,3
            do j=1,natoms_
               coord(i,j)=geo(i,j)
            enddo
         enddo
         return
      endif

      ! Option (C)

      coord(1,1) = 0.d0
      coord(2,1) = 0.d0
      coord(3,1) = 0.d0
      coord(1,2) = geo(1,2)
      coord(2,2) = 0.d0
      coord(3,2) = 0.d0

      if (natoms_.gt.2) then

         ccos = dcos(geo(2,3))

         if (na(3).eq.1) then
             coord(1,3) = coord(1,1) + geo(1,3)*ccos
         else
             coord(1,3) = coord(1,2) - geo(1,3)*ccos
         endif

         coord(2,3) = geo(1,3)*dsin(geo(2,3))
         coord(3,3) = 0.d0

         do i=4,natoms_

            cosa = dcos(geo(2,i))
            mb = nb(i)
            mc = na(i)
            xb = coord(1,mb) - coord(1,mc)
            yb = coord(2,mb) - coord(2,mc)
            zb = coord(3,mb) - coord(3,mc)
            rbc = 1.00/dsqrt(xb*xb+yb*yb+zb*zb)

            if (dabs(cosa).ge.0.999999d0) then

               ! atoms mc, mb, and (i) are collinear

               rbc = geo(1,i)*rbc*cosa
               coord(1,i) = coord(1,mc) + xb*rbc
               coord(2,i) = coord(2,mc) + yb*rbc
               coord(3,i) = coord(3,mc) + zb*rbc
               cycle

            endif

            ! the atoms are not collinear

            ma = nc(i)
            xa = coord(1,ma) - coord(1,mc)
            ya = coord(2,ma) - coord(2,mc)
            za = coord(3,ma) - coord(3,mc)

            ! rotate about the z-axis to make yb=0, and xb positive.  if xyb is
            ! too small, first rotate the y-axis by 90 degrees.

            xyb=dsqrt(xb*xb+yb*yb)
            k=-1

            if (xyb.le.0.1d0) then
               xpa = za
               za = -xa
               xa = xpa
               xpb = zb
               zb = -xb
               xb = xpb
               xyb = dsqrt(xb*xb+yb*yb)
               k = +1
            endif

            ! rotate about the y-axis to make zb vanish

            costh = xb/xyb
            sinth = yb/xyb
            xpa = xa*costh + ya*sinth
            ypa = ya*costh - xa*sinth
            sinph = zb*rbc
            cosph = dsqrt(dabs(1.d0-sinph*sinph))
            xqa = xpa*cosph + za*sinph
            zqa = za*cosph - xpa*sinph

            ! rotate about the x-axis to make za=0, and ya positive.

            yza = dsqrt(ypa*ypa + zqa*zqa)

            if (yza.lt.1.d-1 ) then

               if (yza.lt.1.d-5) goto 21

                  write(6,'(//20x,'' calculation abandoned at this point'')')
                  write(6,'(//10x,'' three atoms being used to define the'',/&
                  &10x,'' coordinates of a fourth atom, whose bond-angle is'',/&
                  &10x,'' not zero or 180 degreees, are in an almost straight'',/&
                  &10x,'' line.'')')
                  write(6,'(//20x,''the faulty atom is atom number'',i4)')i
                  call geoput(geo,na,nb,nc,natoms_,labels_)
                  write(6,'(//20x,''cartesian coordinates up to faulty atom'')')
                  write(6,'(//5x,''i'',12x,''x'',12x,''y'',12x,''z'')')
                  do j=1,i
                     write(6,'(i6,f16.5,2f13.5)') j,(coord(k,j),k=1,3)
                  enddo
                  write(6,'(//6x,'' atoms'',i3,'','',i3,'', and'',i3,&
                  &'' are within'',f11.4,'' Angstroms of a straight line'')')&
                  mc,mb,ma,yza
                  stop

            endif

            coskh = ypa/yza
            sinkh = zqa/yza

            goto 22

  21        continue

            ! angle too small to be important

            coskh = 1.d0
            sinkh = 0.d0

  22        continue

            ! coordinates :-   A=(xqa,yza,0),   B=(rbc,0,0),  C=(0,0,0)
            ! none are negative.
            ! the coordinates of i are evaluated in the new frame.

            sina =  dsin(geo(2,i))
            sind = -dsin(geo(3,i))
            cosd = dcos(geo(3,i))
            xd = geo(1,i)*cosa
            yd = geo(1,i)*sina*cosd
            zd = geo(1,i)*sina*sind

            ! transform the coordinates back to the origional system.

            ypd = yd*coskh  -  zd*sinkh
            zpd = zd*coskh  +  yd*sinkh
            xpd = xd*cosph  - zpd*sinph
            zqd = zpd*cosph +  xd*sinph
            xqd = xpd*costh - ypd*sinth
            yqd = ypd*costh + xpd*sinth

            if (k.ge.1) then
               xrd = -zqd
               zqd =  xqd
               xqd =  xrd
            endif

            coord(1,i) = xqd + coord(1,mc)
            coord(2,i) = yqd + coord(2,mc)
            coord(3,i) = zqd + coord(3,mc)

         enddo

      endif

      ! now remove the dummy atom coordinates, if any.

      j = 0
      do i=1,natoms_
         if (labels_(i).ne.99) then
            j = j + 1
            do k=1,3
               coord(k,j) = coord(k,i)
            enddo
            labels_(j) = labels_(i)
         endif
      enddo
      natoms_ = j

      return

   end subroutine convrt

   !---------------------------------------------------------------------
   subroutine rotate(nat,ndp,nda,xyz)
   !---------------------------------------------------------------------
   !   Rotates the coordinate system to put atoms
   !   NDP and NDA on the X axis.
   !
   !   --------x--------------+---------x-----------> x-axis
   !          NDP             0        NDA
   !
   !   (originally written by Helene Decornez, University of Notre Dame)
   !   (modified by Alexander Soudackov, University of Notre Dame)
   !---------------------------------------------------------------------
      implicit none

      integer, intent(in)                    :: nat, ndp, nda
      real(8), intent(inout), dimension(:,:) :: xyz

      integer :: i, ia, j
      real(8) :: zda, yda, xda, rr2, rr, rrp2, rrp
      real(8) :: cosphi, sinphi, costhe, sinthe, costhe1, sinthe1
      real(8), dimension(nat) :: xt, yt, zt, xt1, yt1, zt1
      logical :: ok

      ! Copy to storage array

      do ia=1,nat
         xt(ia) = xyz(1,ia)
         yt(ia) = xyz(2,ia)
         zt(ia) = xyz(3,ia)
      enddo

      ! Calculate rotation angles
      !
      ! PHI - the angle between the projection of 0-DA
      !     on xy plane and x-axis
      !
      ! THE - the angle between 0-DA and the projection
      !     of 0-DA on xy plane

      zda = zt(nda)
      yda = yt(nda)
      xda = xt(nda)
      rr2 = xda*xda + yda*yda + zda*zda
      rr = dsqrt(rr2)
      rrp2 = xda*xda + yda*yda
      rrp = dsqrt(rrp2)

      if (dabs(rrp).lt.1.d-12) then
         cosphi = 1.d0
         sinphi = 0.d0
      else
         cosphi = xda/rrp
         sinphi = yda/rrp
      endif
      costhe = zda/rr
      sinthe = rrp/rr

      ! First rotate on -PHI around z-axis

      do i=1,nat
         xt1(i) =  xt(i)*cosphi + yt(i)*sinphi
         yt1(i) = -xt(i)*sinphi + yt(i)*cosphi
         zt1(i) =  zt(i)
      enddo

      ! Then rotate on -PI/2+THE around y-axis

      costhe1 = sinthe
      sinthe1 = costhe
      do i=1,nat
         xyz(1,i) =  xt1(i)*costhe1 + zt1(i)*sinthe1
         xyz(2,i) =  yt1(i)
         xyz(3,i) = -xt1(i)*sinthe1 + zt1(i)*costhe1
      enddo

      ! Check whether the atom NDA is on a positive x-semiaxis
      ! If not then rotate 180 degrees about z-axis

      ok = xyz(1,nda).gt.0.d0
      if (.not.ok) then
         do ia=1,nat
            do j=1,2
               xyz(j,ia) = -xyz(j,ia)
            enddo
         enddo
      endif

      return

   end subroutine rotate

   !---------------------------------------------------------------------
   subroutine printb
   !---------------------------------------------------------------------
   ! Prints banner and related information for a current job
   !---------------------------------------------------------------------
      implicit none
      integer :: i, itit,itit2
      integer :: length
      integer :: lengt1, lengt2, lengt3, lengt4, lengt5
      integer :: lengt6, lengt7, lengt8, lengt9, lengt10

      !-- Determine the actual lengths of TITLE and TITLE2

      do i=len(title),1,-1
         if (title(i:i).ne.' ') then
            itit = i
            exit
         endif
      enddo

      do i=len(title),1,-1
         if (title2(i:i).ne.' ') then
            itit2 = i
            exit
         endif
      enddo

      !-- Print banner

      call banner(trim(program_version))

      write(6,'(3x,"Title of the job:  ",a)') title (1:itit)
      write(6,'(22x,a)')                      title2(1:itit2)
      write(6,'(3x,"Start date:        ",a)') trim(strdat)
      write(6,'(3x,"Job identificator: ",a)') trim(job)
      write(6,'(1x,79("-"))')

      !-- Print keywords

      do i=len(keywrd),1,-1
         if (keywrd(i:i).ne.' ') then
            length = i
            exit
         endif
      enddo

      if (length.le.80) then

         write(6,'(1x,a)') keywrd(2:length)

      elseif (length.le.160) then

         do i=80,1,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt1 = i
               exit
            endif
         enddo

         write(6,'(1x,a)') keywrd(2:lengt1)
         write(6,'(1x,a)') keywrd(lengt1+1:length)

      elseif (length.le.240) then

         do i=160,81,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt2 = i
               exit
            endif
         enddo

         do i=80,1,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt1 = i
               exit
            endif
         enddo

         write(6,'(1x,a)') keywrd(2:lengt1)
         write(6,'(1x,a)') keywrd(lengt1+1:lengt2)
         write(6,'(1x,a)') keywrd(lengt2+1:length)

      elseif (length.le.320) then

         do i=240,161,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt3 = i
               exit
            endif
         enddo

         do i=160,81,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt2 = i
               exit
            endif
         enddo

         do i=80,1,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt1 = i
               exit
            endif
         enddo

         write(6,'(1x,a)') keywrd(2:lengt1)
         write(6,'(1x,a)') keywrd(lengt1+1:lengt2)
         write(6,'(1x,a)') keywrd(lengt2+1:lengt3)
         write(6,'(1x,a)') keywrd(lengt3+1:length)

      elseif (length.le.400) then

         do i=320,241,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt4 = i
               exit
            endif
         enddo

         do i=240,161,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt3 = i
               exit
            endif
         enddo

         do i=160,81,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt2 = i
               exit
            endif
         enddo

         do i=80,1,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt1 = i
               exit
            endif
         enddo

         write(6,'(1x,a)') keywrd(2:lengt1)
         write(6,'(1x,a)') keywrd(lengt1+1:lengt2)
         write(6,'(1x,a)') keywrd(lengt2+1:lengt3)
         write(6,'(1x,a)') keywrd(lengt3+1:lengt4)
         write(6,'(1x,a)') keywrd(lengt4+1:length)

      elseif (length.le.480) then

         do i=400,321,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt5 = i
               exit
            endif
         enddo

         do i=320,241,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt4 = i
               exit
            endif
         enddo

         do i=240,161,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt3 = i
               exit
            endif
         enddo

         do i=160,81,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt2 = i
               exit
            endif
         enddo

         do i=80,1,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt1 = i
               exit
            endif
         enddo

         write(6,'(1x,a)') keywrd(2:lengt1)
         write(6,'(1x,a)') keywrd(lengt1+1:lengt2)
         write(6,'(1x,a)') keywrd(lengt2+1:lengt3)
         write(6,'(1x,a)') keywrd(lengt3+1:lengt4)
         write(6,'(1x,a)') keywrd(lengt4+1:lengt5)
         write(6,'(1x,a)') keywrd(lengt5+1:length)


      elseif (length.le.560) then

         do i=480,401,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt6 = i
               exit
            endif
         enddo

         do i=400,321,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt5 = i
               exit
            endif
         enddo

         do i=320,241,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt4 = i
               exit
            endif
         enddo

         do i=240,161,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt3 = i
               exit
            endif
         enddo

         do i=160,81,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt2 = i
               exit
            endif
         enddo

         do i=80,1,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt1 = i
               exit
            endif
         enddo

         write(6,'(1x,a)') keywrd(2:lengt1)
         write(6,'(1x,a)') keywrd(lengt1+1:lengt2)
         write(6,'(1x,a)') keywrd(lengt2+1:lengt3)
         write(6,'(1x,a)') keywrd(lengt3+1:lengt4)
         write(6,'(1x,a)') keywrd(lengt4+1:lengt5)
         write(6,'(1x,a)') keywrd(lengt5+1:lengt6)
         write(6,'(1x,a)') keywrd(lengt6+1:length)


      elseif (length.le.640) then

         do i=560,481,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt7 = i
               exit
            endif
         enddo

         do i=480,401,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt6 = i
               exit
            endif
         enddo

         do i=400,321,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt5 = i
               exit
            endif
         enddo

         do i=320,241,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt4 = i
               exit
            endif
         enddo

         do i=240,161,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt3 = i
               exit
            endif
         enddo

         do i=160,81,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt2 = i
               exit
            endif
         enddo

         do i=80,1,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt1 = i
               exit
            endif
         enddo

         write(6,'(1x,a)') keywrd(2:lengt1)
         write(6,'(1x,a)') keywrd(lengt1+1:lengt2)
         write(6,'(1x,a)') keywrd(lengt2+1:lengt3)
         write(6,'(1x,a)') keywrd(lengt3+1:lengt4)
         write(6,'(1x,a)') keywrd(lengt4+1:lengt5)
         write(6,'(1x,a)') keywrd(lengt5+1:lengt6)
         write(6,'(1x,a)') keywrd(lengt6+1:lengt7)
         write(6,'(1x,a)') keywrd(lengt7+1:length)


      elseif (length.le.720) then

         do i=640,561,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt8 = i
               exit
            endif
         enddo

         do i=560,481,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt7 = i
               exit
            endif
         enddo

         do i=480,401,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt6 = i
               exit
            endif
         enddo

         do i=400,321,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt5 = i
               exit
            endif
         enddo

         do i=320,241,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt4 = i
               exit
            endif
         enddo

         do i=240,161,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt3 = i
               exit
            endif
         enddo

         do i=160,81,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt2 = i
               exit
            endif
         enddo

         do i=80,1,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt1 = i
               exit
            endif
         enddo

         write(6,'(1x,a)') keywrd(2:lengt1)
         write(6,'(1x,a)') keywrd(lengt1+1:lengt2)
         write(6,'(1x,a)') keywrd(lengt2+1:lengt3)
         write(6,'(1x,a)') keywrd(lengt3+1:lengt4)
         write(6,'(1x,a)') keywrd(lengt4+1:lengt5)
         write(6,'(1x,a)') keywrd(lengt5+1:lengt6)
         write(6,'(1x,a)') keywrd(lengt6+1:lengt7)
         write(6,'(1x,a)') keywrd(lengt7+1:lengt8)
         write(6,'(1x,a)') keywrd(lengt8+1:length)


      elseif (length.le.800) then

         do i=720,641,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt9 = i
               exit
            endif
         enddo

         do i=640,561,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt8 = i
               exit
            endif
         enddo

         do i=560,481,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt7 = i
               exit
            endif
         enddo

         do i=480,401,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt6 = i
               exit
            endif
         enddo

         do i=400,321,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt5 = i
               exit
            endif
         enddo

         do i=320,241,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt4 = i
               exit
            endif
         enddo

         do i=240,161,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt3 = i
               exit
            endif
         enddo

         do i=160,81,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt2 = i
               exit
            endif
         enddo

         do i=80,1,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt1 = i
               exit
            endif
         enddo

         write(6,'(1x,a)') keywrd(2:lengt1)
         write(6,'(1x,a)') keywrd(lengt1+1:lengt2)
         write(6,'(1x,a)') keywrd(lengt2+1:lengt3)
         write(6,'(1x,a)') keywrd(lengt3+1:lengt4)
         write(6,'(1x,a)') keywrd(lengt4+1:lengt5)
         write(6,'(1x,a)') keywrd(lengt5+1:lengt6)
         write(6,'(1x,a)') keywrd(lengt6+1:lengt7)
         write(6,'(1x,a)') keywrd(lengt7+1:lengt8)
         write(6,'(1x,a)') keywrd(lengt8+1:lengt9)
         write(6,'(1x,a)') keywrd(lengt9+1:length)

      else

         do i=800,721,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt10 = i
               exit
            endif
         enddo

         do i=720,641,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt9 = i
               exit
            endif
         enddo

         do i=640,561,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt8 = i
               exit
            endif
         enddo

         do i=560,481,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt7 = i
               exit
            endif
         enddo

         do i=480,401,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt6 = i
               exit
            endif
         enddo

         do i=400,321,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt5 = i
               exit
            endif
         enddo

         do i=320,241,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt4 = i
               exit
            endif
         enddo

         do i=240,161,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt3 = i
               exit
            endif
         enddo

         do i=160,81,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt2 = i
               exit
            endif
         enddo

         do i=80,1,-1
            if (keywrd(i:i).eq.' '.or.keywrd(i:i).eq.',') then
               lengt1 = i
               exit
            endif
         enddo

         write(6,'(1x,a)') keywrd(2:lengt1)
         write(6,'(1x,a)') keywrd(lengt1+1:lengt2)
         write(6,'(1x,a)') keywrd(lengt2+1:lengt3)
         write(6,'(1x,a)') keywrd(lengt3+1:lengt4)
         write(6,'(1x,a)') keywrd(lengt4+1:lengt5)
         write(6,'(1x,a)') keywrd(lengt5+1:lengt6)
         write(6,'(1x,a)') keywrd(lengt6+1:lengt7)
         write(6,'(1x,a)') keywrd(lengt7+1:lengt8)
         write(6,'(1x,a)') keywrd(lengt8+1:lengt9)
         write(6,'(1x,a)') keywrd(lengt9+1:lengt10)
         write(6,'(1x,a)') keywrd(lengt10+1:length)

      endif

      write(6,'(1x,79(''-''))')

      return

   end subroutine printb

   !---------------------------------------------------------------------
   subroutine geoput(geo,na,nb,nc,natoms_,labels_)
   !---------------------------------------------------------------------
   ! GEOPUT prints the current geometry in MOPAC Z-matrix format
   !---------------------------------------------------------------------
      implicit none

      integer, intent(in)                 :: natoms_
      integer, intent(in), dimension(:)   :: labels_, na, nb, nc
      real(8), intent(in), dimension(:,:) :: geo

      real(8), parameter :: degree = 57.295779510d0
      integer :: i
      real(8) :: w, x

      write (6,20)
      20 FORMAT (/6X,'ATOM',4X,'CHEMICAL',3X,'BOND LENGTH',4X,'BOND ANGLE'&
      ,4X ,'TWIST ANGLE',/5X,'NUMBER',3X,'SYMBOL', 5X,'(angstroms)',5&
      X,'(degrees)',5X,'(degrees)',/6X,'(I)',20X,'NA:I',10X,'NB:NA:I',5&
      X,'NC:NB:NA:I',4X,'NA',2X,'NB',2X,'NC',/)

      write (6,'(8x,1h1,3x,a4)') elsym(labels_(1))

      if (labels_(2).ne.0)&
         write (6,'(8x,1h2,3x,a4,f16.5,36x,i2)') elsym(labels_(2)),geo(1,2),na(2)

      w = geo(2,3) * degree

      if (labels_(3).ne.0)&
        write (6,'(8x,1h3,3x,a4,f16.5,2x,f15.5,19x,2(i2,2x))')&
        elsym(labels_(3)),geo(1,3),w,na(3),nb(3)

      if (natoms_.lt.4) return

      do i=4,natoms_
         w = geo(2,i) * degree
         x = geo(3,i) * degree
         write (6,'(7x,i2 ,3x,a4,f16.5,2x,f15.5,2x,f12.5,2x,i5,2i4)')&
         i,elsym(labels_(i)),geo(1,i),w,x,na(i),nb(i),nc(i)
      enddo

      return

   end subroutine geoput

   !---------------------------------------------------------------------
   subroutine xyzout(out,nat,lab,xyz)
   !---------------------------------------------------------------------
   !  Writes out the cartesian coordinates to the external file
   !  in xyz format (can be viewed with MOLDEN)
   !---------------------------------------------------------------------
      implicit none

      integer, intent(in)                 :: out, nat
      integer, intent(in), dimension(:)   :: lab
      real(8), intent(in), dimension(:,:) :: xyz

      character(4) :: symb
      integer :: i, k

      write(out,'(i5)') nat
      write(out,'(1x)')

      do i=1,nat

         if (lab(i).eq.103) then
            symb = '  De'
         elseif (lab(i).eq.104) then
            symb = '  Ae'
         elseif (lab(i).eq.105) then
            symb = '  Dp'
         elseif (lab(i).eq.106) then
            symb = '  Ap'
         elseif (lab(i).eq.107) then
            symb = '  xx'
         else
            symb = elsym(lab(i))
         endif

         write(out,'(1x,a4,3x,3f15.6)') symb, (xyz(k,i),k=1,3)

      enddo

      return

   end subroutine xyzout

   !---------------------------------------------------------------------
   ! check if elements of iarr are within the range ( < n )
   !---------------------------------------------------------------------
   function check_range(n,iarr) result(flag)

      integer, intent(in) :: n
      integer, intent(in), dimension(:) :: iarr
      logical :: flag

      integer :: i, isize

      flag = .true.
      isize = size(iarr)

      do i=1,isize
         if (iarr(i).lt.1.or.iarr(i).gt.n) then
            flag = .false.
            exit
         endif
      enddo

   end function check_range


   !---------------------------------------------------------------------
end module geometry
