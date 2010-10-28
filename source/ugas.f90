subroutine ugas(k, kg, hd, ha12, ha, hdes, ha12es, haes)

!===================================================================C
!
!  Calculates the gas-phase and electronically solvated
!  energy profiles along the proton coordinate
!
!  June 14, 2000: Alexander Soudackov, University of Notre Dame
!
!--------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:37 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
!
!===================================================================C

   use pardim
   use gasmat
   use solmat
   use eispack

   implicit none
      
   integer, intent(in) :: k, kg
   real(8),  intent(out), dimension(nelst) :: hd, ha12, ha, hdes, ha12es, haes

   integer :: i, ierr

   real(8), dimension(nelst,nelst) :: h0k, h0k_tmp, cu
   real(8), dimension(nelst) :: work, work1

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Gas-phase 4x4 Hamiltonian (H0K)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   h0k = h0(:,:,k,kg)

   !do i=1,4
   !   do j=1,4
   !      h0k(i,j) = h0(i,j,k,kg)
   !   enddo
   !enddo
   !call primat(h0k,4,4,4,5,6,0,'Gas phase Hamiltonian')

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Gas-phase diabatic states
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   do i=1,nelst
      hd(i) = h0k(i,i)
   enddo

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Gas-phase ET-diabatic/PT-adiabatic states
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   cu = 0.d0
   h0k_tmp = h0k
   h0k_tmp(1,3) = 0.d0
   h0k_tmp(3,1) = 0.d0
   h0k_tmp(1,4) = 0.d0
   h0k_tmp(4,1) = 0.d0
   h0k_tmp(2,3) = 0.d0
   h0k_tmp(3,2) = 0.d0
   h0k_tmp(2,4) = 0.d0
   h0k_tmp(4,2) = 0.d0
   call rs(nelst,nelst,h0k_tmp,ha12,4,cu,work,work1,ierr)

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Gas-phase adiabatic states
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   cu = 0.d0
   call rs(nelst,nelst,h0k,ha,4,cu,work,work1,ierr)

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Electronic solvation of the gas-phase Hamiltonian
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !tinf44 = tinf(1,1,k) + tinf(2,2,k) + tinf(3,3,k) + &
   !&        2.d0*tinf(2,3,k) - 2.d0*tinf(1,2,k) - 2.d0*tinf(1,3,k)

   h0k(1,1) = h0k(1,1) - tinf(1,1,k,kg)/2.d0
   h0k(2,2) = h0k(2,2) - tinf(2,2,k,kg)/2.d0
   h0k(3,3) = h0k(3,3) - tinf(3,3,k,kg)/2.d0
   h0k(4,4) = h0k(4,4) - tinf(4,4,k,kg)/2.d0

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Gas-phase diabatic states with electronic solvation
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   do i=1,nelst
      hdes(i) = h0k(i,i)
   enddo

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Gas-phase ET-diabatic/PT-adiabatic states
   ! with electronic solvation
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   cu = 0.d0
   work = 0.d0
   work1 = 0.d0
   h0k_tmp = h0k
   h0k_tmp(1,3) = 0.d0
   h0k_tmp(3,1) = 0.d0
   h0k_tmp(1,4) = 0.d0
   h0k_tmp(4,1) = 0.d0
   h0k_tmp(2,3) = 0.d0
   h0k_tmp(3,2) = 0.d0
   h0k_tmp(2,4) = 0.d0
   h0k_tmp(4,2) = 0.d0
   call rs(nelst,nelst,h0k_tmp,ha12es,4,cu,work,work1,ierr)

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Gas-phase adiabatic states with electronic solvation
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   cu = 0.d0
   work = 0.d0
   work1 = 0.d0
   call rs(nelst,nelst,h0k,haes,4,cu,work,work1,ierr)

   return

end subroutine ugas
