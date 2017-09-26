subroutine GTOeval(info)
use MOL_info
use GRID_INFO
    implicit none
real(8),parameter   :: ans2bohr = 1.889725989
real(8),parameter   :: PI = 3.14159265359
real(8),parameter   :: covrad(26) = [0.38,0.32,1.34,0.9,0.82,0.77,0.75,&
                                     0.73,0.71,0.69,1.54,1.30,1.18,1.11,&
                                     1.06,1.02,0.99,0.97,1.96,1.74,1.44,&
                                     1.36,1.25,1.27,1.39,1.25 ]
!INCLUDE 'parameter.h'
    integer    :: i,j,k,l,n,m,start,info
    integer    :: a,b,c
    real(8)    :: val,Nor
    real(8)    :: val1(3)
    real(8)    :: COORE(3),COORR(3)
    real(8)    :: valt,rho
    real(8)    :: drv1(3)

   real(8),allocatable   ::  valRec(:,:)

! COORE(1) = 0.0
! COORE(2) = 0.5
! COORE(3) =-1.0

valt = 0
  do k = 1,size(Grids)
    allocate(Grids(k)%val0(nconts))
    allocate(Grids(k)%val1(nconts,3))
    l =1
    rho = 0
    do i = 1,natoms
       COORR = Grids(k)%coor-atoms(i)%coor*ans2bohr
       do j = 1,size(atoms(i)%shell)
          do a =atoms(i)%shell(j)%angMoment,0,-1
             do b = atoms(i)%shell(j)%angMoment-a,0,-1
                c = atoms(i)%shell(j)%angMoment-a-b
                  call GTOeval0(COORR,a,b,c,atoms(i)%shell(j)%nGauss,atoms(i)%shell(j)%contrCoeff,&
                  atoms(i)%shell(j)%exponents,val)
                  call GTOeval1(COORR,a,b,c,atoms(i)%shell(j)%nGauss,atoms(i)%shell(j)%contrCoeff,&
                  atoms(i)%shell(j)%exponents,val1)
                  Grids(k)%val0(l) = val
                  Grids(k)%val1(l,1) = val1(1)
                  Grids(k)%val1(l,2) = val1(2)
                  Grids(k)%val1(l,3) = val1(3)
                l = l +1
             enddo 
          enddo
       enddo
    enddo
  enddo

end subroutine




!/*
! * deriv 0: exp(-ar^2) x^n                     ! Done
! * deriv 1: exp(-ar^2)[nx^{n-1} - 2ax^{n+1}]   ! Done
! * deriv 2: exp(-ar^2)[n(n-1)x^{n-2} - 2a(2n+1)x^n + 4a^2x^{n+2}]  !TODO
! * deriv 3: exp(-ar^2)[n(n-1)(n-2)x^{n-3} - 2a3n^2x^{n-1} + 4a^2(3n+3)x^{n+1} - 8a^3x^{n+3}]
! * deriv 4: exp(-ar^2)[n(n-1)(n-2)(n-3)x^{n-4} - 2a(4n^3-6n^2+2)x^n{-2}
! *                     + 4a^2(6n^2+6n+3)x^n - 8a(4n+6)x^{n+2} + 16a^4x^{n+4}]
! */

subroutine GTOeval0(coor,a,b,c,nGauss,coeff,expo,val)
implicit none
real(8),parameter   :: ans2bohr = 1.889725989
real(8),parameter   :: PI = 3.14159265359
real(8),parameter   :: covrad(26) = [0.38,0.32,1.34,0.9,0.82,0.77,0.75,&
                                     0.73,0.71,0.69,1.54,1.30,1.18,1.11,&
                                     1.06,1.02,0.99,0.97,1.96,1.74,1.44,&
                                     1.36,1.25,1.27,1.39,1.25 ]
!INCLUDE 'parameter.h'
real(8)      :: coor(3)
integer      :: a,b,c,nGauss
real(8)      :: coeff(nGauss),expo(nGauss)
real(8)      :: Nor1,val,Nor2
real(8)      :: rs
integer      :: Lx1,Ly1,Lz1,Lx2,Ly2,Lz2,Lt
integer,external ::factorial,factorial2
integer      :: i,j,k,l

!====================
!  F(r) = N \sum Ci (x^{a}y^{b}z^{c}*exp(-expo*r2))
!  r2 = x^2 + y^2 + z^2 
!  
!
! * deriv 0: exp(-ar^2) x^n
!===================
Lt = a+b+c
Lx1 = factorial2(a)
Ly1 = factorial2(b)
Lz1 = factorial2(c)
val = 0
do i = 1,nGauss
   Nor2 = (2**(2*Lt+1.5)*expo(i)**(Lt+1.5))**0.5 * (Lx1 * Ly1 * Lz1 * PI **(1.5))**(-0.5)
   rs = coor(1)**2 + coor(2)**2 + coor(3)**2
   val =coor(1)**a * coor(2)**b * coor(3)**c * &
         DEXP(-expo(i)*rs) * coeff(i) * Nor2+ val
enddo
end subroutine

subroutine GTOeval1(coor,a,b,c,nGauss,coeff,expo,val)
implicit none
real(8),parameter   :: ans2bohr = 1.889725989
real(8),parameter   :: PI = 3.14159265359
real(8),parameter   :: covrad(26) = [0.38,0.32,1.34,0.9,0.82,0.77,0.75,&
                                     0.73,0.71,0.69,1.54,1.30,1.18,1.11,&
                                     1.06,1.02,0.99,0.97,1.96,1.74,1.44,&
                                     1.36,1.25,1.27,1.39,1.25 ]
!INCLUDE 'parameter.h'
real(8)      :: coor(3)
integer      :: a,b,c,nGauss
real(8)      :: coeff(nGauss),expo(nGauss)
real(8)      :: Nor1,Nor2
real(8)      :: val(3)
real(8)      :: rs
integer      :: Lx1,Ly1,Lz1,Lx2,Ly2,Lz2,Lt
integer,external ::factorial,factorial2
integer      :: i,j,k,l

!====================
!  F(r) = N \sum Ci (x^{a}y^{b}z^{c}*exp(-expo*r2))
!  r2 = x^2 + y^2 + z^2 
!  
! * deriv 1: exp(-ar^2)[nx^{n-1} - 2ax^{n+1}]
!
!===================
Lt = a+b+c
Lx1 = factorial2(a)
Ly1 = factorial2(b)
Lz1 = factorial2(c)

val = 0
Nor1=1.0
do i = 1,nGauss
   Nor2 = (2**(2*Lt+1.5)*expo(i)**(Lt+1.5))**0.5 * (Lx1 * Ly1 * Lz1 * PI **(1.5))**(-0.5)
   rs = coor(1)**2 + coor(2)**2 + coor(3)**2

   if (a .eq. 0) then
   val(1) = -Nor1*2*expo(i)*coor(1)**(a+1) *&
            coor(2)**b * coor(3)**c *DEXP(-expo(i)*rs) * coeff(i) * Nor2+ val(1)
   else
   val(1) = Nor1*( a*coor(1)**(a-1) - 2*expo(i)*coor(1)**(a+1) )*&
            coor(2)**b * coor(3)**c *DEXP(-expo(i)*rs) * coeff(i) * Nor2+ val(1)
   endif

   if (b .eq. 0) then
   val(2) = -Nor1*2*expo(i)*coor(2)**(b+1)*&
            coor(1)**a * coor(3)**c *DEXP(-expo(i)*rs) * coeff(i) * Nor2+ val(2)
       
   else
   val(2) = Nor1*( b*coor(2)**(b-1) - 2*expo(i)*coor(2)**(b+1) )*&
            coor(1)**a * coor(3)**c *DEXP(-expo(i)*rs) * coeff(i) * Nor2+ val(2)
   endif


   if (c .eq. 0) then
   val(3) = -Nor1*2*expo(i)*coor(3)**(c+1)*&
            coor(2)**b * coor(1)**a *DEXP(-expo(i)*rs) * coeff(i) * Nor2+ val(3)
   
   else
   val(3) = Nor1*( c*coor(3)**(c-1) - 2*expo(i)*coor(3)**(c+1) )*&
            coor(2)**b * coor(1)**a *DEXP(-expo(i)*rs) * coeff(i) * Nor2+ val(3)
   endif
enddo
end subroutine

