subroutine initmat
!===================================================================C
!  Initializes the matrices on the grid
!--------------------------------------------------------------------
!
!  souda
!  2010/06/25 20:02:35
!  4.1
!  Exp
!  initmat.f90,v 4.1 2010/06/25 20:02:35 souda Exp
!  initmat.f90,v
!  Revision 4.1  2010/06/25 20:02:35  souda
!  Release 4.1
!
!  Revision 1.1.1.1  2004/01/07 17:14:57  souda
!  Initial PCET-4.0 Release
!
!
!===================================================================C
   use pardim
   use control
   use gasmat
   use solmat
   use quantum

   implicit none
   logical :: deriv, derivg

   ! do we need to allocate the arrays for derivatives?

   deriv  = method.eq.2.and.(igas.eq.1.or.igas.eq.4.or.igas.eq.5).and.isolv.eq.1
   derivg = gquant.and.mgquant.eq.1.and.(igas.eq.4.or.igas.eq.5)

   ! allocate arrays for matrices on the grid
   ! by calling corresponding module procedures

   call alloc_gasmat(npnts,npntsg,deriv,derivg)
   call alloc_solmat(npnts,npntsg,deriv)
   call alloc_pquant(npnts,deriv)
   call alloc_gquant(npntsg,derivg)

   return

end subroutine initmat

