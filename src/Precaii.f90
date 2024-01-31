!*   See Olsen,Yeager, Joergensen:                                      *
!*    "Optimization and characterization of an MCSCF state"             *
!*                                                                      *
!*     Eq. C.12a

! =====================================================================
! yma: Actually it is Hess_ia,ia when i is closed shell
!      and I don't know why Hess_ia,ia not match Hess_ai,ai 
!      in closed shell case??
!
! I add this part in order to match Molcas and let the LR code runs
! The reason for this part I need to check later.  2016.10.25
!

      Subroutine precaii(nsub,occ,cls,tot,norb,nocc,ncls,D,p,&
                         FI,FA,F,U,Diag)

! Actually all Hess_ai,bj is calculated, but only diagional part is used

        integer::nsub,nocc,norb,ncls
        integer::occ(nsub),cls(nsub),tot(nsub)
        double precision::diag(ncls,nocc-ncls)
        double precision::d(nocc,nocc)
        double precision::p(nocc,nocc,nocc,nocc)
        double precision::U(norb,norb,norb,norb)
        double precision::FI(norb,norb),FA(norb,norb),F(norb,norb)

!       only active
        double precision::da(nocc-ncls,nocc-ncls)
        double precision::pa(nocc-ncls,nocc-ncls,nocc-ncls,nocc-ncls)

        integer act(nsub)
        double precision dv,dv1,dv2
        double precision,allocatable::TM1(:,:)
        double precision,allocatable::GM1(:,:,:,:)

        double precision,allocatable::GM2(:,:,:,:)
        double precision,allocatable::GM3(:,:,:,:)

        act=occ-cls
        nact=nocc-ncls
 
        da=0.0d0
        pa=0.0d0      
        diag=0.0d0

        allocate(GM2(nact,nact,ncls,ncls))
        allocate(GM3(nact,ncls,ncls,nact))
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
 
                ix=cls(i)+i0
                jx=cls(j)+j0

                GM2=U(ix+1:ix+act(i),jx+1:jx+act(j), &
                      k0+1:k0+cls(k),l0+1:l0+cls(l))

                GM3=U(ix+1:ix+act(i),k0+1:k0+cls(k), &
                      l0+1:l0+cls(l),jx+1:jx+act(j))

                l0=l0+tot(l)
              end do
              k0=k0+tot(k)
            end do
            j0=j0+tot(j)
          end do
          i0=i0+tot(i)
        end do

        !GM2 : a,b,i,j
        !GM3 : a,i,j,b , actually  (a,i,b,j)

        allocate(GM1(nact,ncls,nact,ncls)) 
        GM1=0.0d0

        dv=0.0d0
! only for no symmetry/C1
        do i=1,ncls   
          do j=1,ncls 
            !do k=nact,nact        !a 
            do k=1,nact        !a 
              !do l=nact,nact      !b
              do l=1,nact      !b
                dv=0.0d0
                do m=1,nact    !c
                  do n=1,nact  !d

! Althought it match molcas, I still doublt the results
                    dv =dv+2.0d0*pa(l,n,k,m)*GM2(m,n,i,j)
                    dv =dv+1.0d0*pa(l,k,n,m)*GM3(m,i,j,n)

                    dv =dv+2.0d0*pa(l,m,k,n)*GM2(m,n,i,j)
                    dv =dv+1.0d0*pa(l,m,k,n)*GM3(m,i,j,n)

!                    write(*,*)"rdens1",pa(l,n,k,m),pa(l,m,k,n)
!                    write(*,*)"rdens2",pa(l,k,n,m),pa(l,m,k,n)
!                    write(*,*)"cdij",GM2(m,n,i,j)
!                    write(*,*)"cidj",GM3(m,i,j,n)

!                    write(*,*)"dv",dv

                  end do
                end do
                GM1(k,i,l,j)= GM1(k,i,l,j) + dv
              end do
            end do
          end do
        end do

!        write(*,*)"dv start",dv
!        stop

        do i=1,ncls
          do j=1,ncls
            !do k=nact,nact        !a 
            do k=1,nact        !a 
              !do l=nact,nact      !b
              do l=1,nact      !b
                dv=0.0d0
                do m=1,nact    !c

                    dv1=0.0d0
                    if(k.eq.m)then
                      dv1=1.0d0-da(k,m)
                    else
                      dv1=     -da(k,m)  
                    end if                  
!                    write(*,*)"dv1/rdens (ac)",dv1
                    dv1=dv1*2.0d0

!                    write(*,*)"cibj",GM3(m,i,j,l)
!                    write(*,*)"bicj",GM3(l,i,j,m)
!                    write(*,*)"bcij",GM2(l,m,i,j)

                    dv2=0.0d0
                    dv2=4.0d0*GM3(m,i,j,l)-GM3(l,i,j,m)-GM2(l,m,i,j) 
!                    write(*,*)"dv2 1/2",dv2
                    dv =dv+dv1*dv2

                    dv1=0.0d0
                    if(l.eq.m)then
                      dv1=1.0d0-da(l,m)
                    else
                      dv1=     -da(l,m)
                    end if
!                    write(*,*)"dv1/rdens (bc)",dv1
                    dv1=dv1*2.0d0

!                    write(*,*)"cjai",GM3(m,j,i,k)
!                    write(*,*)"ajci",GM3(k,j,i,m)
!                    write(*,*)"acij",GM2(k,m,i,j)

                    dv2=0.0d0
                    dv2=4.0d0*GM3(m,j,i,k)-GM3(k,j,i,m)-GM2(k,m,i,j)
!                    write(*,*)"dv2 2/2",dv2
                    dv =dv+dv1*dv2

!                    write(*,*)"dv",dv

                end do
                GM1(k,i,l,j)= GM1(k,i,l,j) + dv
              end do
            end do
          end do
        end do

!        write(*,*)"dv",dv

        do i=1,ncls
          do j=1,ncls
            !do k=nact,nact        !a 
            do k=1,nact        !a 
              !do l=nact,nact      !b
              do l=1,nact      !b

                    dv=0.0d0

                    kx=k+ncls
                    lx=l+ncls

                    dv1=0.0d0
                    if(k.eq.l)then
                      dv1=4.0d0
                    end if
                    dv2=0.0d0
                    dv2=FI(i,j)+FA(i,j)
!                    write(6,*)"dv -- 0",dv,"rdens",da(k,l)
                    dv =dv-dv1*dv2
                    dv =dv+2.0d0*da(k,l)*FI(i,j)  ! ?
!                    write(6,*)"dv -- 1",dv

                    dv1=4.0d0*FI(kx,lx)+4.0d0*FA(kx,lx)-2.0d0*F(lx,kx)
!                    write(6,*)"dv -- 2",dv
                    dv =dv+dv1

                    GM1(k,i,l,j)= GM1(k,i,l,j) + dv

              end do
            end do
          end do
        end do

!        write(*,*)"dv",dv
!        stop

        write(6,*)" Full Hess_ia,bi"
!        call print_gat(nact,ncls,nact,ncls,GM1,6)

        do i=1,ncls
          do j=1,nact
            diag(i,j)=GM1(j,i,j,i)
          end do
        end do

!        write(6,*)" Diag"
!        call print_mat(ncls,nact,diag,6)

        !deallocate(TM1)
        deallocate(GM1,GM2,GM3)

      End subroutine precaii

