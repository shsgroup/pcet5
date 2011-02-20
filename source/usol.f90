subroutine usol(mode,k,kg,zp,ze,u,cu)

!======================================================================!
!
!  Calculates the free energy U(i,q,ZP,ZE) at the grid point K
!
!  MODE  - calculation mode:
!        = 'ADIAB' - adiabatic energies and wavefunctions
!        = 'DIAB2' - ET diabatic energies and wavefunctions
!        = 'DIAB4' - diabatic energies and wavefunctions
!
!  K     - the proton position (GRID POINT);
!
!  KG    - gating distance (GRID POINT);
!
!  ZP,ZE - the (PT) and (ET) medium coordinates (kcal/mole);
!
!  SW    - self-energy of the inertial polarization (kcal/mole).
!
!  U     - the free energies (kcal/mole):
!        if MODE = 'ADIAB' then adiabatic energies (sorted)
!        if MODE = 'DIAB2' then ET diabatic energies
!                   U(1),U(2) - 1a/1b pair
!                   U(3),U(4) - 2a/2b pair
!        if MODE = 'DIAB4' then diabatic energies
!
!  CU    - the corresponding eigenvectors
!
!-----------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2011-02-20 00:58:11 $
!  $Revision: 5.3 $
!  $Log: not supported by cvs2svn $
!  Revision 5.2  2010/10/28 21:29:37  souda
!  First (working and hopefully bug-free) source of PCET 5.x
!
!
!======================================================================!

   use pardim
   use gasmat
   use solmat
   use eispack

   implicit none

   character(5), intent(in) :: mode
   integer,      intent(in) :: k, kg
   real*8,       intent(in) :: zp, ze
   real*8,       intent(out), dimension(nelst)       :: u
   real*8,       intent(out), dimension(nelst,nelst) :: cu

   integer :: i, j, ierr
   real*8  :: sw

   real*8, allocatable, dimension(:,:) :: cu1, cu2, h0k, hs, hs1, hs2
   real*8, allocatable, dimension(:)   :: w, w1, w2, work1, work2

   allocate (h0k(4,4),hs(4,4))

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Gas-phase 4x4 Hamiltonian (H0K)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   do i=1,4
      do j=1,4
         h0k(i,j) = h0(i,j,k,kg)
      enddo
   enddo

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Electronic solvation of the gas-phase Hamiltonian
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   !tinf44 = tinf(1,1,k) + tinf(2,2,k) + tinf(3,3,k)&
   !      &+ 2.d0*tinf(2,3,k) - 2.d0*tinf(1,2,k) - 2.d0*tinf(1,3,k)

   h0k(1,1) = h0k(1,1) - tinf(1,1,k,kg)/2.d0
   h0k(2,2) = h0k(2,2) - tinf(2,2,k,kg)/2.d0
   h0k(3,3) = h0k(3,3) - tinf(3,3,k,kg)/2.d0
   h0k(4,4) = h0k(4,4) - tinf(4,4,k,kg)/2.d0

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Solvated Hamiltonian (HS)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   hs = h0k
   hs(1,1) = 0.d0
   hs(2,2) = h0k(2,2) - h0k(1,1) + zp
   hs(3,3) = h0k(3,3) - h0k(1,1) + ze
   hs(4,4) = h0k(4,4) - h0k(1,1) + zp + ze

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Self-energy term
   ! (summation over linearly independent solvent variables)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   sw = selfen(k,kg,zp,ze)

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Energies and eigenvectors
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   u = 0.d0
   cu = 0.d0

   if (mode.eq.'DIAB4') then

      do i=1,4
         u(i) = hs(i,i) + sw + h0k(1,1) - tr(1,1,k,kg)/2.d0
         cu(i,i) = 1.d0
      enddo

   elseif (mode.eq.'DIAB2') then

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Diagonalization of the electronically diabatic
      ! blocks of the total solvated Hamiltonian
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      allocate (hs1(2,2),hs2(2,2))
      allocate (cu1(2,2),cu2(2,2))
      allocate (w1(2),w2(2))
      allocate (work1(2),work2(2))

      cu1 = 0.d0
      cu2 = 0.d0
      w1 = 0.d0
      w2 = 0.d0
      work1 = 0.d0
      work2 = 0.d0

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! block 1a-1b
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      hs1(1,1) = hs(1,1)
      hs1(2,2) = hs(2,2)
      hs1(1,2) = hs(1,2)
      hs1(2,1) = hs(2,1)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! block 2a-2b
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      hs2(1,1) = hs(3,3)
      hs2(2,2) = hs(4,4)
      hs2(1,2) = hs(3,4)
      hs2(2,1) = hs(4,3)

      call rs(2,2,hs1,w1,2,cu1,work1,work2,ierr)
      call rs(2,2,hs2,w2,2,cu2,work1,work2,ierr)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Total ET diabatic free energies (kcal/mole)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do i=1,2
         u(i) = w1(i)   + sw + h0k(1,1) - tr(1,1,k,kg)/2.d0
      enddo
      do i=3,4
         u(i) = w2(i-2) + sw + h0k(1,1) - tr(1,1,k,kg)/2.d0
      enddo

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! ET diabatic electronic wavefunctions (eigenvectors)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      cu(1:2,1:2) = cu1
      !do i=1,2
      !   do j=1,2
      !      cu(j,i) = cu1(j,i)
      !   enddo
      !enddo

      cu(1:2,3:4) = cu2
      !do i=3,4
      !   do j=1,2
      !      cu(j,i) = cu2(j,i-2)
      !   enddo
      !enddo

      deallocate (hs1,hs2)
      deallocate (cu1,cu2)
      deallocate (w1,w2)
      deallocate (work1,work2)

   elseif (mode.eq.'ADIAB') then

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Diagonalization of the total solvated Hamiltonian
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      allocate (w(4),work1(4),work2(4))

      w = 0.d0
      work1 = 0.d0
      work2 = 0.d0
      call rs(4,4,hs,w,4,cu,work1,work2,ierr)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Total adiabatic free energies (kcal/mole)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do i=1,4
         u(i) = w(i) + sw + h0k(1,1) - tr(1,1,k,kg)/2.d0
      enddo
      deallocate (w,work1,work2)

   else

      write(*,'(/1x,''*** (in USOL): '',/,&
      &''MODE is neither ADIAB, nor DIAB2, nor DIAB4 ... ***''/)')
      deallocate (h0k,hs)
      stop

   endif

   deallocate (h0k,hs)

end subroutine usol



















subroutine usol_diab(k,kg,u)

   !======================================================================!
   !  Calculates the diabatic electronic energies U(i,q) on the grid
   !  K     - the proton position (GRID POINT);
   !  KG    - gating distance (GRID POINT);
   !  U     - the diabatic electronic energies (kcal/mole):
   !======================================================================!

   use pardim
   use gasmat
   use solmat

   implicit none

   integer, intent(in) :: k, kg
   real*8,  intent(out), dimension(nelst) :: u

   integer :: i, j
   real*8, allocatable, dimension(:,:) :: h0k, hs, hs1

   allocate (h0k(4,4),hs(4,4))

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Gas-phase 4x4 Hamiltonian (H0K)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   do i=1,4
      do j=1,4
         h0k(i,j) = h0(i,j,k,kg)
      enddo
   enddo

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Electronic solvation of the gas-phase Hamiltonian
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   !tinf44 = tinf(1,1,k) + tinf(2,2,k) + tinf(3,3,k)&
   !      &+ 2.d0*tinf(2,3,k) - 2.d0*tinf(1,2,k) - 2.d0*tinf(1,3,k)

   h0k(1,1) = h0k(1,1) - tinf(1,1,k,kg)/2.d0
   h0k(2,2) = h0k(2,2) - tinf(2,2,k,kg)/2.d0
   h0k(3,3) = h0k(3,3) - tinf(3,3,k,kg)/2.d0
   h0k(4,4) = h0k(4,4) - tinf(4,4,k,kg)/2.d0

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Reduced Hamiltonian (HS)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   hs = h0k
   hs(1,1) = 0.d0
   hs(2,2) = h0k(2,2) - h0k(1,1)
   hs(3,3) = h0k(3,3) - h0k(1,1)
   hs(4,4) = h0k(4,4) - h0k(1,1)

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Energies
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   u = 0.d0
   do i=1,4
      u(i) = hs(i,i) + h0k(1,1) - tr(1,1,k,kg)/2.d0
   enddo

   deallocate (h0k,hs)

end subroutine usol_diab
