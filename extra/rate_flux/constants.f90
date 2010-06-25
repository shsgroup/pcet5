module constants

   implicit none
   public

   ! conversion factors and constants
   real(8), parameter :: pi = 3.141592654d0
   real(8), parameter :: au2kcal = 627.5095d0
   real(8), parameter :: hbarps = 4.7685388d-2/(pi*au2kcal)   ! Planck constant in au*ps
   real(8), parameter :: kb = 1.9872159d-3                    ! Boltzmann constant in kcal/(mol*K)
   real(8), parameter :: ev2cm = 8065.54477d0
   real(8), parameter :: au2ev = 27.2113834d0
   real(8), parameter :: bohr2a = 0.529177d0
   real(8), parameter :: dalton = 1822.888d0
   
   real(8), parameter :: zero  = 0.0d0,&
                       & half  = 0.5d0,&
		       & one   = 1.0d0,&
		       & two   = 2.0d0,&
		       & three = 3.0d0,&
		       & four  = 4.0d0,&
		       & five  = 5.0d0

end module constants
