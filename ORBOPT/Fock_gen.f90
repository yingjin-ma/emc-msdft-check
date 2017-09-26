      Subroutine fock_gen2(nsub,occ,tot,nact,norb,T,U,D,P,A)

! Use the same way as Molcas/<Molecular Electronic-Structure theory>
!   in order to check the SA-Hessian issue

        integer::nsub,nact,norb
        integer::occ(nsub)
        integer::tot(nsub)

        double precision::T(norb,norb)
        double precision::U(norb,norb,norb,norb)

        double precision::D(nact,nact)
        double precision::P(nact,nact,nact,nact)

        double precision::A(norb,norb),Qmat(norb,norb)

        double precision,allocatable::TM1(:,:),TM2(:,:),TM3(:,:)
        double precision,allocatable::TM4(:,:),TM5(:,:),TM6(:,:)
        double precision,allocatable::GM1(:,:,:,:),GM2(:,:,:,:)

        A=0.0d0

!        call fock_I() 
!        call fock_A()
        call Qmat_gen(nsub,occ,tot,nact,norb,U,P,Qmat)     

 
      End subroutine fock_gen2

      Subroutine Qmat_gen(nsub,occ,tot,nact,norb,U,P,Qmat)

! The Qmat    
        integer::nsub,nact,norb
        integer::occ(nsub)
        integer::tot(nsub)

        double precision::U(norb,norb,norb,norb)
        double precision::P(nact,nact,nact,nact)
        double precision::Qmat(norb,norb)

        double precision,allocatable::TM1(:,:),TM2(:,:),TM3(:,:)
        double precision,allocatable::TM4(:,:),TM5(:,:),TM6(:,:)
        double precision,allocatable::GM1(:,:,:,:),GM2(:,:,:,:)
    
        Qmat=0.0d0  

  


      End Subroutine Qmat_gen

      Subroutine fock_A(nsub,occ,cls,tot,nocc,norb,U,D,FA)

! The so-called   active fock matrix (only for c1)
        integer::nsub,norb
        integer::cls(nsub)
        integer::occ(nsub)
        integer::tot(nsub)
        double precision::D(nocc,nocc)
        double precision::U(norb,norb,norb,norb)

        double precision,allocatable::TM1(:,:),TM2(:,:)

        double precision::FA(norb,norb)        

        FA=0.0d0

!        call print_mat(nocc,nocc,D,6)

        allocate(TM1(norb,norb));TM1=0.0d0

        i0=0 
        do i=1,nsub         
          j0=0
          do j=1,nsub           
 
            do i1=cls(i)+1,occ(i)
              do j1=cls(j)+1,occ(j)
                ix=i1+i0
                jx=j1+j0 
!                write(6,*)"ix,jx",ix,jx,D(ix,jx)
!                call print_mat(norb,norb,U(:,:,ix,jx),6)
!                write(6,*)"J/K"
!                call print_mat(norb,norb,U(:,ix,jx,:),6)                 
                  
                TM1=TM1+D(ix,jx)*(1.0d0*U(:,:,ix,jx)-0.5d0*U(:,ix,jx,:))
                ! C.10 should be corrected : 2.0 -> 1.0 
                 
              end do  
            end do
 
            j0=j0+tot(j)
          end do
          i0=i0+tot(i)
        end do

        FA=TM1
        deallocate(TM1)  


      End Subroutine fock_A


      Subroutine fock_I(nsub,cls,tot,norb,T,U,FI)

! The so-called inactive fock matrix (only for C1)
        integer::nsub,norb
        integer::cls(nsub)
        integer::tot(nsub)
        double precision::T(norb,norb)
        double precision::U(norb,norb,norb,norb)

        double precision,allocatable::TM1(:,:),TM2(:,:)

        double precision::FI(norb,norb) 

        FI=0.0d0

        FI=T 
         
        allocate(TM1(norb,norb));TM1=0.0d0

        i0=0               
        do i=1,nsub
          do j=1,cls(i) 
            ix=j+i0
            TM1=TM1+2.0d0*U(:,:,ix,ix)-U(:,ix,ix,:)
          end do
          i0=i0+tot(i)
        end do 

        FI=FI+TM1   
 
        deallocate(TM1) 

!        call print_mat(norb,norb,F1,6)  

      End Subroutine fock_I 
 
! to Stefan : Eq.24 / Eq.31
!             Eq.31 can share with this formula when U+ is already operated on MO-integrals

      Subroutine fock_gen(nsub,occ,tot,nact,norb,T,U,D,P,A)

! A is the Fock matrix 

        integer::nsub,nact,norb
        integer::occ(nsub)
        integer::tot(nsub)
        
        double precision::T(norb,norb)
        double precision::U(norb,norb,norb,norb)

        double precision::D(nact,nact)
        double precision::P(nact,nact,nact,nact)
 
        double precision::A(norb,norb)

        double precision,allocatable::TM1(:,:),TM2(:,:),TM3(:,:)
        double precision,allocatable::TM4(:,:),TM5(:,:),TM6(:,:)
        double precision,allocatable::GM1(:,:,:,:),GM2(:,:,:,:)        
  
        A=0.0d0

        allocate(TM1(norb,nact));TM1=0.0d0
        allocate(TM2(norb,nact));TM2=0.0d0
        i0=0
        ia=0
        do i=1,nsub
          ir=tot(i)
          ij=occ(i)
          TM1(i0+1:i0+ir,ia+1:ia+ij) = T(i0+1:i0+ir,i0+1:i0+ij)
          i0=i0+tot(i)
          ia=ia+occ(i)
        end do
        call MXMG(norb,nact,nact,TM1,D,TM2,"NN")
        i0=0
        ia=0
        do i=1,nsub
          ir=tot(i)
          ij=occ(i)
          A(i0+1:i0+ir,i0+1:i0+ij)=TM2(i0+1:i0+ir,ia+1:ia+ij)
          i0=i0+tot(i)
          ia=ia+occ(i)
        end do
        deallocate(TM1)
        deallocate(TM2)         


!       Continue to construct the A matrix base on T/U integrals
        allocate(TM5(norb,nact)); TM5=0.0d0
        k0=0
        ka=0
        do k=1,nsub; do kk=1,occ(k)
          l0=0
          la=0
          do l=1,nsub; do ll=1,occ(l)

! J^{kl}*P^{lk} 
            allocate(TM1(norb,norb));TM1=0.0d0
            allocate(TM2(norb,nact));TM2=0.0d0
            allocate(TM3(nact,nact));TM3=0.0d0
            allocate(TM4(norb,nact));TM4=0.0d0
            TM1=     U(:,:,kk+k0,ll+l0)
            TM3=     P(:,kk+ka,ll+la,:)

            i0=0
            ia=0
            do i=1,nsub
              ir=tot(i)
              ij=occ(i)

              j0=0
              ja=0
              do j=1,nsub
                jr=tot(j)
                jj=occ(j)

              TM2(i0+1:i0+ir,ja+1:ja+jj)=TM1(i0+1:i0+ir,j0+1:j0+jj)

                j0=j0+tot(j)
                ja=ja+occ(j)
              end do

              i0=i0+tot(i)
              ia=ia+occ(i)
            end do

            call MXMG(norb,nact,nact,TM2,TM3,TM4,"NN")

            TM5=TM5+TM4

            deallocate(TM1)
            deallocate(TM2)
            deallocate(TM3)
            deallocate(TM4)

          end do
          l0=l0+tot(l)
          la=la+occ(l)
          end do
        end do
        k0=k0+tot(k)
        ka=ka+occ(k)
        end do

        allocate(TM1(norb,norb));TM1=0.0d0
        i0=0
        ia=0
        do i=1,nsub
          ir=tot(i)
          ij=occ(i)
          TM1(i0+1:i0+ir,i0+1:i0+ij)=TM5(i0+1:i0+ir,ia+1:ia+ij)
          i0=i0+tot(i)
          ia=ia+occ(i)
        end do

        A=A+TM1
        A=2.0d0*A
        deallocate(TM1)
        deallocate(TM5)

!        write(2,*)"The specific A-mat :"
!        call print_mat(norb,norb,A,2)  ! checked
!        call flush(2)

      End Subroutine fock_gen 

