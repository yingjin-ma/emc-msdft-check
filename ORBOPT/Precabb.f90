!*   See Olsen,Yeager, Joergensen:                                      *
!*    "Optimization and characterization of an MCSCF state"             *
!*                                                                      *
!*     Eq. C12j

! =====================================================================
! yma: Actually it is Hess_ta,ta when t is virtual orbitals
!
! I add this part in order to match Molcas and let the LR code runs
! The reason for this part I need to check later.  2017.1.7
!

      Subroutine precabb(nsub,occ,cls,tot,norb,nocc,ncls,D,p,&
                         FI,FA,F,U,diag)

! 1) Actually all Hess_sa,tb is calculated, but only diagional part is used
! 2) I didn't used the packed form as Molcas

        integer::nsub,nocc,norb,ncls
        integer::occ(nsub),cls(nsub),tot(nsub)
        double precision::d(nocc,nocc)
        double precision::p(nocc,nocc,nocc,nocc)
        double precision::U(norb,norb,norb,norb)
        double precision::FI(norb,norb),FA(norb,norb),F(norb,norb)

        double precision::diag(norb-nocc,nocc-ncls)

!       only active
        double precision::da(nocc-ncls,nocc-ncls)
        double precision::pa(nocc-ncls,nocc-ncls,nocc-ncls,nocc-ncls)

        integer act(nsub),vir(nsub)
        double precision dv,dv1,dv2
        double precision,allocatable::TM1(:,:)
        double precision,allocatable::GM1(:,:,:,:)

        double precision,allocatable::GM2(:,:,:,:)
        double precision,allocatable::GM3(:,:,:,:)

        act=occ-cls
        vir=tot-occ 
        nact=nocc-ncls
        nvir=norb-nocc

        da=0.0d0
        pa=0.0d0
        diag=0.0d0

        allocate(GM2(nvir,nvir,nact,nact))
        allocate(GM3(nvir,nact,nact,nvir))
        GM2=0.0d0
        GM3=0.0d0

! act RDMs and partly active integrals

        i0=0
        do i=1,nsub
          j0=0
          do j=1,nsub
            da(1:act(i),1:act(j))=d(cls(i)+1:occ(i),cls(j)+1:occ(j))

            k0=0
            do k=1,nsub
              l0=0
              do l=1,nsub
                pa(1:act(i),1:act(j),1:act(k),1:act(l)) &
               =p (cls(i)+1:occ(i),cls(j)+1:occ(j),     &
                   cls(k)+1:occ(k),cls(l)+1:occ(l))

                ix=occ(i)+i0
                jx=occ(j)+j0
                kx=cls(k)+k0
                lx=cls(l)+l0                 

                GM2=U(ix+1:ix+vir(i),jx+1:jx+vir(j), &
                      kx+1:kx+act(k),lx+1:lx+act(l))

                GM3=U(ix+1:ix+vir(i),kx+1:kx+act(k), &
                      lx+1:lx+act(l),jx+1:jx+vir(j))

                l0=l0+tot(l)
              end do
              k0=k0+tot(k)
            end do
            j0=j0+tot(j)
          end do
          i0=i0+tot(i)
        end do
         
        ! GM2,GM3, to be check later 
        !call print_gat(nvir,nvir,nact,nact,GM2,127) 
        !call print_gat(nvir,nact,nact,nvir,GM3,128) 
  
!        write(6,*)" == Pa start == "
!        call print_gat(nact,nact,nact,nact,pa,6) 
!        write(6,*)" == Pa  done =="

        allocate(GM1(nvir,nact,nvir,nact))
        GM1=0.0d0          

        dv=0.0d0
! only for no symmetry/C1
        do i=1,nact           ! a
          do j=1,nact         ! b   
            do k=1,nvir        ! s  
              do l=1,nvir      ! t 
                dv=0.0d0
                do m=1,nact     ! c
                  do n=1,nact   ! d

                    dv =dv+2.0d0*pa(i,n,m,j)*GM2(k,l,m,n)

                    !write(*,*)i,n,m,j,"rdens1",pa(i,n,m,j)

                    !write(*,*)"stcd",GM2(k,l,m,n)

                    !write(*,*)"Part1",2.0d0*pa(i,n,m,j)*GM2(k,l,m,n)

                    !write(*,*)"dv",dv

                  end do
                end do
                GM1(k,i,l,j)=GM1(k,i,l,j) + dv
                !write(*,*)"GM1(k,i,l,j)",k,i,l,j,GM1(k,i,l,j)    
                !write(*,*)
              end do
            end do
          end do
        end do

        !write(*,*)"================= no-packed 1 ======================"

! only for no symmetry/C1
        do i=1,nact           ! a
          do j=1,nact         ! b   
            do k=1,nvir        ! s  
              do l=1,nvir      ! t 
                dv=0.0d0
                do m=1,nact     ! c
                  do n=1,nact   ! d

                    dv =dv+2.0d0*pa(i,m,n,j)*GM3(k,n,m,l)

                    !write(*,*)"i,m,n,j",i,m,n,j,"rdens1",pa(i,m,n,j)

                    !write(*,*)"k,n,m,l",k,n,m,l,"sdtc",GM3(k,n,m,l)

                    !write(*,*)"Part2",2.0d0*pa(i,m,n,j)*GM3(k,n,m,l)

                    !write(*,*)"dv",dv

                  end do
                end do
                GM1(k,i,l,j)=GM1(k,i,l,j) + dv
                !write(*,*)"GM1(k,i,l,j)",k,i,l,j,GM1(k,i,l,j)
                !write(*,*)
              end do
            end do
          end do
        end do

        !write(*,*)"================== no-packed 2 ====================="

! only for no symmetry/C1
        do i=1,nact           ! a
          do j=1,nact         ! b   
            do k=1,nvir        ! s  
              do l=1,nvir      ! t 
                dv=0.0d0
                do m=1,nact     ! c
                  do n=1,nact   ! d

                    dv =dv+2.0d0*pa(i,j,m,n)*GM3(k,n,m,l)

                    !write(*,*)"i,j,m,n",i,j,m,n,"rdens1",pa(i,j,m,n)

                    !write(*,*)"k,n,m,l",k,n,m,l,"sdtc",GM3(k,n,m,l)

                    !write(*,*)"Part2",2.0d0*pa(i,j,m,n)*GM3(k,n,m,l)

                    !write(*,*)"dv",dv

                  end do
                end do
                GM1(k,i,l,j)=GM1(k,i,l,j) + dv
                !write(*,*)"GM1(k,i,l,j)",k,i,l,j,GM1(k,i,l,j)
                !write(*,*)
              end do
            end do
          end do
        end do

        dv=0.0d0
! only for no symmetry/C1
        do i=1,nact           ! a
          do j=1,nact         ! b   
            do k=1,nvir        ! s  
              do l=1,nvir      ! t 

                dv=0.0d0

                ix=i+ncls
                jx=j+ncls
                kx=k+nocc
                lx=l+nocc

!                write(*,*)"da(i,j)",da(i,j)  
!                write(*,*)"FI(k,l)",FI(kx,lx),"F(i,j)",F(ix,jx)    
 
                dv = dv + 2.0d0*da(i,j)*FI(kx,lx)
                dv = dv - 2.0d0*F(jx,ix)
                GM1(k,i,l,j)= GM1(k,i,l,j) + dv

!                write(*,*)"GM1(k,i,l,j)",k,i,l,j,GM1(k,i,l,j)    
!                write(*,*)

              end do
            end do
          end do
        end do

        call print_gat(nvir,nact,nvir,nact,GM1,129) 

        do i=1,nvir
          do j=1,nact
            Diag(i,j)=GM1(i,j,i,j)  
          end do
        end do  

        deallocate(GM1)
        deallocate(GM2)
        deallocate(GM3)

      End Subroutine   

