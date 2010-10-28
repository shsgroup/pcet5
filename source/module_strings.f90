module strings

!-------------------------------------------------------------------
!  Routines working with character strings
!-------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:36 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
!
!--------------------------------------------------------------------

   implicit none
   public

   character( 1), parameter :: space = " "
   character( 1), parameter :: comma = ","
   character( 1), parameter :: tab   = char(9)

   character(52), parameter :: allow1 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
   character(19), parameter :: allow2 = "1234567890()=,.:-_+"

!--------------------------------------------------------------------
contains

   real*8 function digit(string,istart)
   !=======================================================================
   ! Fortran function to convert numeric field to double precision
   ! number.  the string is assumed to be clean (no invalid digit
   ! or character combinations from istart to the first nonspace,
   ! nondigit, nonsign, and nondecimal point character).
   !=======================================================================
      implicit none

      character(*), intent(in) :: string
      integer,      intent(in) :: istart

      real*8  :: c1,c2,deciml
      logical :: sign
      integer :: i0,i9,ineg,ipos,idot,ispc
      integer :: i,j,l,idig,n

      ! define ascii values of numeric field characters
      i0   = ichar('0')
      i9   = ichar('9')
      ineg = ichar('-')
      ipos = ichar('+')
      idot = ichar('.')
      ispc = ichar(' ')

      c1 = 0.d0
      c2 = 0.d0
      sign = .true.
      l = len(string)

      ! determine the contribution to the number greater than one

      idig=0

      do i=istart,l

         n = ichar(string(i:i))

         if (n.ge.i0.and.n.le.i9) then

            idig = idig + 1
            c1 = c1*1.d1 + n - i0

         elseif (n.eq.ineg.or.n.eq.ipos.or.n.eq.ispc) then

            if (n.eq.ineg) sign=.false.

         elseif (n.eq.idot) then

            exit

         else

            digit = c1 + c2
            if (.not.sign) digit = -digit
            return

         endif

      enddo

      ! determine the contribution to the number less than than one

      deciml = 1.d0
      do j=i+1,l
         n = ichar(string(j:j))
         if (n.ge.i0.and.n.le.i9) then
            deciml = deciml/1.d1
            c2 = c2 + (n-i0)*deciml
         elseif (n.ne.ispc) then
            exit
         endif
      enddo

      ! put the pieces together

      digit = c1 + c2
      if (.not.sign) digit = -digit
      return

   end function digit

   subroutine getopt(keywrd,ibra,options)
   !======================================================================
   ! Extracts OPTIONS string from the string KEYWRD
   ! starting from position IBRA
   ! (Options should be in parentheses)
   !======================================================================
      implicit none

      character(*), intent(in) :: keywrd
      character(*), intent(inout) :: options
      integer, intent(in) ::  ibra

      integer :: iket, length

      iket = index(keywrd(ibra:),')')

      if (iket.eq.0) then
         write(*,'(/1x,"*** Input error (in GETOPT): ",/,&
        &"Options MUST be in parentheses ***"/)')
         stop 'Try again...'
      endif

      options = keywrd(ibra:ibra+iket-2)
      length = iket - 1

      ! Clean the OPTIONS string, remove all spaces
      ! and replace commas/tabs with spaces

      call clnopt(options,length)

      return

   end subroutine getopt

   subroutine locase(string,length)
   !======================================================================
   ! Converts character string STRING to lowercase
   !======================================================================
      implicit none

      character(*), intent(inout) :: string
      integer,      intent(in)    :: length

      integer :: ilowa, ilowz, icapa, icapz, i, iline

      ilowa = ichar('a')
      ilowz = ichar('z')
      icapa = ichar('A')
      icapz = ichar('Z')

      do i=1,length
         iline=ichar(string(i:i))
         if (iline.ge.icapa.and.iline.le.icapz) then
            string(i:i)=char(iline+ilowa-icapa)
         endif
      enddo

      return

   end subroutine locase

   integer function nblen(string)

      implicit none
      character(*), intent(in) :: string
      integer n

      n = 1
      do while(string(n:n).ne.space.and.string(n:n).ne.comma)
         n = n + 1
      enddo

      nblen = n - 1
      return

   end function nblen

   subroutine rmdbsp(line,length)
   !======================================================================
   ! Removes double spaces and forbidden symbols
   !======================================================================
      implicit none

      character(*), intent(inout) :: line
      integer,      intent(inout) :: length

      integer :: i, icount
      logical :: curr, next
      character(80) ::  string

      icount = 0
      
      do i=1,length
      
         curr = index(allow1,line(i:i)).ne.0.or.index(allow2,line(i:i)).ne.0
         if (i.lt.80) then
            next = index(allow1,line(i+1:i+1)).ne.0.or.index(allow2,line(i+1:i+1)).ne.0
         else
            next = .false.
         endif
         
         if (curr.or.next) then
            icount = icount + 1
            string(icount:icount) = line(i:i)
         endif
         
      enddo

      line = string(1:icount)
      length = icount

      return

   end subroutine rmdbsp

   subroutine rmsp(string,length)
   !======================================================================
   ! Removes all spaces from the character string STRING
   !======================================================================
      implicit none

      character(*), intent(inout) :: string
      integer,      intent(inout) :: length

      character(160) :: temp
      integer :: i, ilen

      ilen = 0
      do i=1,length
         if (string(i:i).ne.space) then
            ilen = ilen + 1
            temp(ilen:ilen) = string(i:i)
         endif
      enddo

      length = ilen
      string = temp(1:ilen)

      return

   end subroutine rmsp

   subroutine clnopt(options,length)
   !=======================================================================
   !     Clean the OPTIONS string, remove all spaces
   !     and replace commas/tabs with spaces
   !=======================================================================

      character(*), intent(inout) :: options
      integer, intent(inout) :: length

      integer :: i

      call rmsp(options,length)
      do i=1,length
         if (options(i:i).eq.tab.or.options(i:i).eq.comma) options(i:i)=space
      enddo
      options=space//options(1:length)//space

      return

   end subroutine clnopt

   character*3 function th(i)

      integer :: i
      integer :: lastdigit

      lastdigit = mod(i,10)

      if (lastdigit.eq.1) then
         th = '-st'
      elseif (lastdigit.eq.2) then
         th = '-nd'
      elseif (lastdigit.eq.3) then
         th = '-rd'
      else
         th = '-th'
      endif

   end function th

   subroutine upcase(string,length)
   !======================================================================
   ! Converts character string STRING to uppercase
   !======================================================================
      implicit none

      character(*), intent(inout) :: string
      integer,      intent(in)    :: length

      integer :: ilowa, ilowz, icapa
      integer :: i, iline

      ilowa = ichar('a')
      ilowz = ichar('z')
      icapa = ichar('A')

      do i=1,length
         iline = ichar(string(i:i))
         if (iline.ge.ilowa.and.iline.le.ilowz) then
            string(i:i) = char(iline+icapa-ilowa)
         endif
      enddo

      return

   end subroutine upcase

   real*8 function reada(string,istart)
   !======================================================================
   ! FORTRAN FUNCTION TO EXTRACT NUMBER FROM STRING (from MOPAC)
   !======================================================================
      implicit none

      character(*), intent(in) :: string
      integer,      intent(in) :: istart

      logical :: defalt, expnnt

      integer :: i0, i9, idot, ineg, ipos, icapd, icape, ismld, ismle
      integer :: l, i, j, iadd, n

      ! DEFINE ASCII VALUES OF NUMERIC FIELD CHARACTERS
      i0    = ichar('0')
      i9    = ichar('9')
      idot  = ichar('.')
      ineg  = ichar('-')
      ipos  = ichar('+')
      icapd = ichar('D')
      icape = ichar('E')
      ismld = ichar('d')
      ismle = ichar('e')

      l = len(string)

      ! find the start of the numeric field
      do i=istart,l

         iadd = 0
         n = ichar(string(i:i))

         ! signal start of numeric field if digit found
         if (n.ge.i0.and.n.le.i9) goto 20

         ! account for consecutive signs [- and(or) +]
         if (n.eq.ineg.or.n.eq.ipos) then
            iadd = iadd + 1
            if (i+iadd.gt.l) goto 50
            n = ichar(string(i+iadd:i+iadd))
            if (n.ge.i0.and.n.le.i9) goto 20
         endif

         ! account for consecutive decimal points (.)
         if (n.eq.idot) then
            iadd = iadd + 1
            if (i+iadd.gt.l) goto 50
            n = ichar(string(i+iadd:i+iadd))
            if (n.ge.i0.and.n.le.i9) goto 20
         endif
      
      enddo
      goto 50

      ! find the end of the numeric field

   20 expnnt = .false.

      do j=i+1,l

         iadd = 0
         n = ichar(string(j:j))

         ! continue search for end if digit found
         if (n.ge.i0.and.n.le.i9) cycle

         ! continue search for end if sign found and expnnt true
         if (n.eq.ineg.or.n.eq.ipos) then
            if (.not.expnnt) goto 40
            iadd = iadd + 1
            if (j+iadd.gt.l) goto 40
            n = ichar(string(j+iadd:j+iadd))
            if (n.ge.i0.and.n.le.i9) cycle
         endif

         if (n.eq.idot) then
            iadd = iadd + 1
            if (j+iadd.gt.l) goto 40
            n = ichar(string(j+iadd:j+iadd))
            if (n.ge.i0.and.n.le.i9) cycle
            if (n.eq.icape.or.n.eq.ismle.or.n.eq.icapd.or.n.eq.ismld) cycle
         endif

         if (n.eq.icape.or.n.eq.ismle.or.n.eq.icapd.or.n.eq.ismld) then
            if (expnnt) goto 40
            expnnt = .true.
            cycle
         endif

         goto 40

      enddo
      j = l + 1

   40 n = ichar(string(j-1:j-1))
      if (n.eq.icape.or.n.eq.ismle.or.n.eq.icapd.or.n.eq.ismld) j=j-1

      ! found the end of the numeric field (it runs 'i' thru 'j-1')

      n = 0
      n = n + index(string(i:j-1),'e')
      n = n + index(string(i:j-1),'E')
      n = n + index(string(i:j-1),'d')
      n = n + index(string(i:j-1),'D')

      if (n.eq.0) then
         reada = digit(string(i:j-1),1)
      else
         reada = digit(string(:i+n-2),i)*1.d1**digit(string(:j-1),i+n)
      endif
      defalt = .false.

      return

      ! default value returned because no numeric field found
   50 reada = 0.d0
      defalt = .true.
      return

   end function reada


end module strings
