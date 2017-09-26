!================================================
! Do DFT calculation
subroutine  DFT_calc(n,Ex,Ec,Fxca,Fxcb,info)
!  Give DFT-energy and KS-Matrix
!    
!================================================

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
    integer    :: i,j,k,l,n,m,info
    integer    :: igrid
    real(8)    :: rho_a,rho_b
    real(8)    :: sigma_aa,sigma_ab,sigma_bb,sigma_t
    real(8)    :: Et 
    real(8)    :: Exc,ex,ec
    real(8)    :: ext,ect
    real(8)    :: rho_t,gamma_aa_t
    real(8)    :: drv_a(3),drv_b(3)
    real(8)    :: Fxca(n,n),Fxcb(n,n)
    real(8),allocatable :: Fxc_at(:,:),Fxc_bt(:,:)
    real(8)    :: d1rho,d1sig,d2rho,d2rhosig,d2sig

    real(8)    :: vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab, &
       v2rhoa2,v2rhob2,v2rhoab,&
       v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,&
       v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,&
       v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,&
       v2sigmaab2,v2sigmaabbb,v2sigmabb2
    rho_t = 0
    gamma_aa_t = 0
    Exc = 0
    Ex = 0
    Ec = 0
    allocate(Fxc_at(nconts,nconts),Fxc_bt(nconts,nconts))
    Fxc_at = 0
    Fxc_bt = 0
    Fxca = 0
    Fxcb = 0

    do igrid = 1,ngrids 
        rho_a=0
        rho_b=0
        sigma_t=0
        sigma_aa=0
        sigma_ab=0
        sigma_bb=0
        drv_a=0
        drv_b=0
        ! 1   Calculat rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb
        do j = 1,nconts
           do k =j,nconts
             ! Calculate /rho_a
             if (j .ne. k) then 
              rho_a = rho_a+Pa(j,k)*grids(igrid)%val0(j)*grids(igrid)%val0(k)*2
              ! Calculate /delta /rho_a
              drv_a(1) = drv_a(1)+2*Pa(j,k)*&
              (grids(igrid)%val0(j)*grids(igrid)%val1(k,1)+grids(igrid)%val1(j,1)*grids(igrid)%val0(k))
              drv_a(2) = drv_a(2)+2*Pa(j,k)*&
              (grids(igrid)%val0(j)*grids(igrid)%val1(k,2)+grids(igrid)%val1(j,2)*grids(igrid)%val0(k))
              drv_a(3) = drv_a(3)+2*Pa(j,k)*&
              (grids(igrid)%val0(j)*grids(igrid)%val1(k,3)+grids(igrid)%val1(j,3)*grids(igrid)%val0(k))
             else 
              rho_a = rho_a+Pa(j,k)*grids(igrid)%val0(j)*grids(igrid)%val0(k)
              drv_a(1) = drv_a(1)+Pa(j,k)*&
              (grids(igrid)%val0(j)*grids(igrid)%val1(k,1)+grids(igrid)%val1(j,1)*grids(igrid)%val0(k))
              drv_a(2) = drv_a(2)+Pa(j,k)*&
              (grids(igrid)%val0(j)*grids(igrid)%val1(k,2)+grids(igrid)%val1(j,2)*grids(igrid)%val0(k))
              drv_a(3) = drv_a(3)+Pa(j,k)*&
              (grids(igrid)%val0(j)*grids(igrid)%val1(k,3)+grids(igrid)%val1(j,3)*grids(igrid)%val0(k))
             endif
                
              if (Multi == -1) then  ! For Close shell
                 rho_b = rho_a
                 drv_b = drv_a
              else                  ! For Open shell 
                 ! Calculate /rho_b
                if (j.ne.k) then
                 rho_b = rho_b + Pb(j,k) *grids(igrid)%val0(j)*grids(igrid)%val0(k)*2
                 ! Calculate /delta /rho_b
                 drv_b(1) = drv_b(1)+2*Pb(j,k)*&
                 (grids(igrid)%val0(j)*grids(igrid)%val1(k,1)+grids(igrid)%val1(j,1)*grids(igrid)%val0(k))
                 drv_b(2) = drv_b(2)+2*Pb(j,k)*&
                 (grids(igrid)%val0(j)*grids(igrid)%val1(k,2)+grids(igrid)%val1(j,2)*grids(igrid)%val0(k))
                 drv_b(3) = drv_b(3)+2*Pb(j,k)*&
                 (grids(igrid)%val0(j)*grids(igrid)%val1(k,3)+grids(igrid)%val1(j,3)*grids(igrid)%val0(k))
                else
                 rho_b = rho_b + Pb(j,k) *grids(igrid)%val0(j)*grids(igrid)%val0(k)
                 ! Calculate /delta /rho_b
                 drv_b(1) = drv_b(1)+Pb(j,k)*&
                 (grids(igrid)%val0(j)*grids(igrid)%val1(k,1)+grids(igrid)%val1(j,1)*grids(igrid)%val0(k))
                 drv_b(2) = drv_b(2)+Pb(j,k)*&
                 (grids(igrid)%val0(j)*grids(igrid)%val1(k,2)+grids(igrid)%val1(j,2)*grids(igrid)%val0(k))
                 drv_b(3) = drv_b(3)+Pb(j,k)*&
                 (grids(igrid)%val0(j)*grids(igrid)%val1(k,3)+grids(igrid)%val1(j,3)*grids(igrid)%val0(k))
                endif
              endif
           enddo
        enddo
        if (Multi == -1) then ! For Close shell
            ! Calculate gamma
            sigma_t=(drv_a(1)+drv_b(1))**2 + &
                    (drv_a(2)+drv_b(2))**2 + &
                    (drv_a(3)+drv_b(3))**2
        else ! For Open shell
            ! Calculate gamma_aa
            sigma_aa=drv_a(1)**2 + drv_a(2)**2 + drv_a(3)**2 
            ! Calculate gamma_bb
            sigma_bb=drv_b(1)**2 + drv_b(2)**2 + drv_b(3)**2
            ! Calculate gamma_ab
            sigma_ab=drv_a(1)*drv_b(1)+drv_a(2)*drv_b(2)+drv_a(3)*drv_b(3)
        endif
        call XCinterface(n,Multi,1,rho_a,rho_b,sigma_aa,sigma_bb,sigma_ab,sigma_t,igrid,&
                          drv_a,drv_b,1,Ext,Ect,Fxc_at,Fxc_bt)
        Ex = Ext*grids(igrid)%weight+Ex
        Ec = Ect*grids(igrid)%weight+Ec
        Fxca = Fxca + Fxc_at*grids(igrid)%weight
        Fxcb = Fxcb + Fxc_bt*grids(igrid)%weight
    enddo
deallocate(Fxc_at,Fxc_bt)
end subroutine

!=========
subroutine XCinterface(n,imult,idrv,rho_a,rho_b,sigma_aa,sigma_bb,sigma_ab,sigma_t,igrd,&
                       drv_a,drv_b,ifunc,ex,ec,Fxca,Fxcb)
!  XClib interface
!=========
use GRID_INFO
integer :: imult,idrv !Multipcility,Deriv: 1: 1st-deriv. 2: 2rd-deriv
real(8) :: rho_a,rho_b,sigma_aa,sigma_bb,sigma_ab,sigma_t
real(8) :: weight,E
integer :: i,j,k

integer :: n,igrd
real(8) :: drv_a(3),drv_b(3)
real(8) :: Fxca(n,n),Fxcb(n,n)
real(8) :: ex,ec
real(8)    :: d1rho,d1sig,d2rho,d2rhosig,d2sig

real(8)    :: d1rhoc,d1sigc,d2rhoc,d2rhosigc,d2sigc
real(8)    :: d1rhox,d1sigx,d2rhox,d2rhosigx,d2sigx

real(8)    :: vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab, &
   v2rhoa2,v2rhob2,v2rhoab,&
   v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,&
   v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,&
   v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,&
   v2sigmaab2,v2sigmaabbb,v2sigmabb2

real(8)    :: vrhoac,vrhobc,vsigmaaac,vsigmabbc,vsigmaabc, &
   v2rhoa2c,v2rhob2c,v2rhoabc,&
   v2rhoasigmaaac,v2rhoasigmaabc,v2rhoasigmabbc,&
   v2rhobsigmabbc,v2rhobsigmaabc,v2rhobsigmaaac,&
   v2sigmaaa2c,v2sigmaaaabc,v2sigmaaabbc,&
   v2sigmaab2c,v2sigmaabbbc,v2sigmabb2c

real(8)    :: vrhoax,vrhobx,vsigmaaax,vsigmabbx,vsigmaabx, &
   v2rhoa2x,v2rhob2x,v2rhoabx,&
   v2rhoasigmaaax,v2rhoasigmaabx,v2rhoasigmabbx,&
   v2rhobsigmabbx,v2rhobsigmaabx,v2rhobsigmaaax,&
   v2sigmaaa2x,v2sigmaaaabx,v2sigmaaabbx,&
   v2sigmaab2x,v2sigmaabbbx,v2sigmabb2x


real(8) :: Temp_A,Temp_B(3),Temp_C(3)



!1. Call XC_fucntional 
!   get E_ex and E_corr
!
if (imult .eq. -1) then ! For Close shell
!    call rks_x_pbe(1,idrv,rho_a+rho_b,sigma_t,ex,d1rhox,d1sigx,d2rhox,d2rhosigx,d2sigx)
    call rks_c_pbe(1,idrv,rho_a+rho_b,sigma_t,ec,d1rhoc,d1sigc,d2rhoc,d2rhosigc,d2sigc)
    d1rho = d1rhox + d1rhoc
    d1sig = (d1sigx + d1sigc)/4.0
    d2rho = d2rhox + d2rhoc
    d2rhosig = d2rhosigx + d2rhosigc
    d2sig = d2sigc + d2sigx
else ! For Open shell
!      call uks_x_pbe                         &
!       (1,idrv,rho_a,rho_b,sigma_aa,sigma_bb,sigma_ab,  &
!       ex,vrhoax,vrhobx,vsigmaaax,vsigmabbx,vsigmaabx,          &
!       v2rhoa2x,v2rhob2x,v2rhoabx,                            &
!       v2rhoasigmaaax,v2rhoasigmaabx,v2rhoasigmabbx,          &
!       v2rhobsigmabbx,v2rhobsigmaabx,v2rhobsigmaaax,          &
!       v2sigmaaa2x,v2sigmaaaabx,v2sigmaaabbx,                 &
!       v2sigmaab2x,v2sigmaabbbx,v2sigmabb2x)
      call uks_c_pbe                         &
       (1,idrv,rho_a,rho_b,sigma_aa,sigma_bb,sigma_ab,  &
       ec,vrhoac,vrhobc,vsigmaaac,vsigmabbc,vsigmaabc,          &
       v2rhoa2c,v2rhob2c,v2rhoabc,                            &
       v2rhoasigmaaac,v2rhoasigmaabc,v2rhoasigmabbc,          &
       v2rhobsigmabbc,v2rhobsigmaabc,v2rhobsigmaaac,          &
       v2sigmaaa2c,v2sigmaaaabc,v2sigmaaabbc,                 &
       v2sigmaab2c,v2sigmaabbbc,v2sigmabb2c)

       vrhoa = vrhoax + vrhoac
       vrhob = vrhobx + vrhobc
       vsigmaaa = vsigmaaax + vsigmaaac
       vsigmabb = vsigmabbx + vsigmabbc
       vsigmaab = vsigmaabx + vsigmaabc
       v2rhoa2 =  v2rhoa2x + v2rhoa2c
       v2rhob2 = v2rhob2x + v2rhob2c
       v2rhoab = v2rhoabx + v2rhoabc
       v2rhoasigmaaa = v2rhoasigmaaax + v2rhoasigmaaac
       v2rhoasigmaab = v2rhoasigmaabx + v2rhoasigmaabc
       v2rhoasigmabb = v2rhoasigmabbx + v2rhoasigmabbc
       v2rhobsigmabb = v2rhobsigmabbx + v2rhobsigmabbc
       v2rhobsigmaab = v2rhobsigmaabx + v2rhobsigmaabc
       v2rhobsigmaaa = v2rhobsigmaaax + v2rhobsigmaaac
       v2sigmaaa2 = v2sigmaaa2x + v2sigmaaa2c
       v2sigmaaaab = v2sigmaaaabx + v2sigmaaaabc
       v2sigmaaabb = v2sigmaaabbx + v2sigmaaabbc
       v2sigmaab2 = v2sigmaab2x + v2sigmaab2c
       v2sigmaabbb = v2sigmaabbbx + v2sigmaabbbc
       v2sigmabb2 =v2sigmabb2x + v2sigmabb2c

endif  
Fxca = 0
Fxcb = 0

do i = 1,n
   do j = i,n
      Temp_C(1) = (grids(igrd)%val0(j)*grids(igrd)%val1(i,1)+&
                   grids(igrd)%val1(j,1)*grids(igrd)%val0(i))
 
      Temp_C(2) = (grids(igrd)%val0(j)*grids(igrd)%val1(i,2)+&
                   grids(igrd)%val1(j,2)*grids(igrd)%val0(i))

      Temp_C(3) = (grids(igrd)%val0(j)*grids(igrd)%val1(i,3)+&
                   grids(igrd)%val1(j,3)*grids(igrd)%val0(i))
      if (imult .eq. -1) then ! For Close shell
           Temp_A = d1rho*grids(igrd)%val0(i)*grids(igrd)%val0(j)
           Temp_B = 2*d1sig*(drv_a+drv_b)
           Fxca(i,j) = Temp_A + (Temp_B(1)*Temp_C(1) + Temp_B(2)*Temp_C(2) +Temp_B(3)*Temp_C(3)) 
           Fxcb(i,j) = Fxca(i,j)

      else ! For Open shell
           Temp_A = vrhoa*grids(igrd)%val0(i)*grids(igrd)%val0(j)
           Temp_B = 2*vsigmaaa*drv_a + vsigmaab*drv_b
           Fxca(i,j) = Temp_A + (Temp_B(1)*Temp_C(1) + Temp_B(2)*Temp_C(2) +Temp_B(3)*Temp_C(3)) 
          
           Temp_A = vrhob*grids(igrd)%val0(i)*grids(igrd)%val0(j)
           Temp_B = 2*vsigmabb*drv_b + vsigmaab*drv_a
           Fxcb(i,j) = Temp_A + (Temp_B(1)*Temp_C(1) + Temp_B(2)*Temp_C(2) +Temp_B(3)*Temp_C(3)) 
      endif
      Fxca(j,i) = Fxca(i,j)
      Fxcb(j,i) = Fxcb(i,j)
   enddo
enddo
end subroutine
