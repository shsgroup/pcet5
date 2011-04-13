!!c*******************************************************************
!!c
!!c   MODULE:      strings_utilities
!!c   PKG VERSION: DLPROTEIN-2.1w
!!c   ACTION:      contains subroutines and functions working on
!!c                character strings and their conversions
!!c
!!c------------------------------------------------------------------
!!c
!!c   $Author: souda $
!!c   $Date: 2011-04-13 22:18:57 $
!!c   $Revision: 1.1 $
!!c   $Log: not supported by cvs2svn $
!!c   Revision 3.4  2008/03/07 23:22:52  souda
!!c   (1) pbc_mod module added (periodic boundar conditions)
!!c
!!c   (2) config2pdb now has the option of folding back into the box;
!!c       (it will ask the question)
!!c
!!c   UPDATE and ENJOY!
!!c
!!c   Revision 3.3  2007/11/01 22:57:47  souda
!!c   Some nasty bugs fixed (affects only GROMOS force field):
!!c
!!c   Revision 3.2  2007/06/09 00:05:40  souda
!!c   minor changes
!!c
!!c   Revision 3.1  2007/02/20 00:51:21  souda
!!c   MInor changes and fixes
!!c
!!c
!!c*******************************************************************

!!---------------------------------------------------------
!! This module combines all the routines dealing
!! with reading the data from the input files
!!---------------------------------------------------------
MODULE string_utilities

   IMPLICIT NONE
   private

   !-------------
   !! parameters

   character(len=1), parameter, dimension(0:9) :: n=(/"0","1","2","3","4","5","6","7","8","9"/)
   character(len=1), parameter :: dot=".", d="d", e="e"
   character(len=1), parameter :: space=" "

   !--------------------------------------------------------------------
   !! declaration of access rights for module subroutines and functions

   PUBLIC :: dblstr, intstr, readrec, strip, lowcase

CONTAINS

   !! strips the string from leading and trailing blanks
   subroutine strip(string)
      character(len=*), intent(inout) :: string
      integer :: i, len_string
      character(len=1) :: tab=achar(9)
      len_string = len_trim(string)
      do i=1,len_string
         if (string(1:1).eq.space.or.string(1:1).eq.tab) then
	    string = string(2:)
	 else
	    exit
	 endif
      enddo
      string = trim(string)
   end subroutine strip

   !! converts string to lower case
   subroutine lowcase(string)
      character(len=*), intent(inout) :: string
      character(len=1) :: letter
      integer :: i, iletter
      integer :: lowa = ichar("a"), lowz = ichar("z"),&
               & capa = ichar("A"), capz = ichar("Z")
      do i=1,len(string)
         iletter = ichar(string(i:i))
         if (iletter.ge.capa.and.iletter.le.capz) string(i:i) = char(iletter+lowa-capa)
      enddo
   end subroutine lowcase

   !!---------------------------------------------------------
   !! extracts first integer from a character string
   !!---------------------------------------------------------
   function intstr(word) result(ivalue)

      character(len=*) :: word
      character(len=1) :: ksn="+"
      integer :: ivalue

      logical :: flag, final
      integer :: isn, lst, j

      isn = 1
      ivalue = 0
      final=.false.

      do lst=1,len(word)

         flag=.false.

         do j=0,9
            if (n(j).eq.word(lst:lst)) then
               ivalue = 10*ivalue + j
               flag=.true.
            endif
         enddo

         if (flag.and.ksn.eq.'-') isn = -1
         ksn = word(lst:lst)

         if (flag) then
            final=.true.
         else if (final) then
            ivalue = isn*ivalue
            return
         endif

      enddo

      ivalue = isn*ivalue

   end function intstr

   !!---------------------------------------------------------
   !! extracts first real number from a character string
   !!---------------------------------------------------------
   function dblstr(word) result(rvalue)

      character(len=*) :: word
      real(kind=8) :: rvalue

      character(len=100) :: work
      character(len=1) :: ksn="+"
      logical :: flag, ldot, start
      real(kind=8) :: sn, ten, one

      integer :: iexp, idum, lst, word_length, i, j

      rvalue = 0.d0
      ten = 10.d0
      sn = 1.d0
      one = 1.d0
      iexp = 0
      idum = 0
      start = .false.
      ldot = .false.
      work = space

      word_length = len(word)

      do lst=1,word_length

         flag=.false.

         do j=0,9
            if (n(j).eq.word(lst:lst)) then
               rvalue = ten*rvalue+one*real(j)
               flag=.true.
               start=.true.
            endif
         enddo

         if (dot.eq.word(lst:lst)) then
            flag=.true.
            ten = 1.d0
            ldot =.true.
            start=.true.
         endif

         if (flag.and.ksn.eq.'-') sn=-1.d0

         ksn = word(lst:lst)
         if (ksn.eq."D") ksn="d"
         if (ksn.eq."E") ksn="e"

         if (ldot) one = one/10.d0

         if (start) then
            if (d.eq.ksn.or.e.eq.ksn) then
               do i = 1,word_length-lst
                  work(i:i) = word(i+lst:i+lst)
               enddo
               iexp = intstr(work)
               exit
            endif
            if (.not.flag) exit
         endif

      enddo

      rvalue = sn*rvalue*(10.d0**iexp)

   end function dblstr


   !!---------------------------------------------------------
   !! This subroutine reads a record and extracts the entries
   !! from it according to the code of the record
   !!---------------------------------------------------------
   subroutine readrec(record,codeline,carr,iarr,rarr)

      character(len=*), intent(in) :: record
      character(len=*), intent(in) :: codeline
      character(len=*), dimension(:), intent(inout) :: carr
      integer,          dimension(:), intent(inout) :: iarr
      real(kind=8),     dimension(:), intent(inout) :: rarr

      logical :: check=.true.
      character(len=500) :: record_tmp
      character(len=7) :: subroutine_name="readrec"
      integer :: i, k, i_carr, i_iarr, i_rarr, ispace

      carr = ''
      iarr = 0
      rarr = 0.d0

      i_carr = 0
      i_iarr = 0
      i_rarr = 0

      record_tmp = record
      call strip(record_tmp)
      record_tmp = trim(record_tmp)//space

      do i=1,len_trim(codeline)

         if (codeline(i:i).eq.'c') then

            !! character string
            ispace = scan(record_tmp,space)
            i_carr = i_carr + 1

            if (i_carr.gt.size(carr)) then
               write(*,*) trim(subroutine_name),": size of carr array is too small..."
	       stop
            endif

            !-- check if the first character in record_tmp is a digit

	    check = .true.
	    do k=0,9
	       if (record_tmp(1:1).eq.n(k)) then
	          check = .false.
		  exit
	       endif
	    enddo

            !if (.not.check) then
            !   write(*,*) trim(subroutine_name),"Character string starts with a digit (",n(k),"): something is wrong..."
	    !   stop
            !endif

            carr(i_carr) = record_tmp(1:ispace-1)
            record_tmp = adjustl(record_tmp(ispace+1:))

         elseif (codeline(i:i).eq.'i') then

            !! integer
            ispace = scan(record_tmp,space)
            i_iarr = i_iarr + 1
            if (i_iarr.gt.size(iarr)) then
               write(*,*) trim(subroutine_name),": size of iarr array is too small..."
	       stop
            endif
            iarr(i_iarr) = intstr(record_tmp(:ispace))
            record_tmp = adjustl(record_tmp(ispace+1:))

         elseif (codeline(i:i).eq.'r') then

            !! float
            ispace = scan(record_tmp,space)
            i_rarr = i_rarr + 1
            if (i_rarr.gt.size(rarr)) then
               write(*,*) trim(subroutine_name),": size of rarr array is too small..."
	       stop
            endif
            rarr(i_rarr) = dblstr(record_tmp(:ispace))
            record_tmp = adjustl(record_tmp(ispace+1:))

         else

            write(*,*) trim(subroutine_name),": Wrong code in the codeline..."
	    stop

         endif

      enddo

   end subroutine readrec

END MODULE string_utilities
