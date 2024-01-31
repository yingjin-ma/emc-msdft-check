! to Stefan : eq.8 / eq.20 / eq.30
!             eq.20/30 can share with this formula if when orbital change (T or dR) is already integrated into MO-integrals 
               

      Subroutine energy_gen(nsub,occ,tot,nact,norb,T,U,D,P,energy)

! energy is the calculated energy basing on the inputed integrals (T,U) and rdms (D,P)

        integer::nsub,nact,norb
        integer::occ(nsub)
        integer::tot(nsub)

        double precision::T(norb,norb)
        double precision::U(norb,norb,norb,norb)

        double precision::D(nact,nact)
        double precision::P(nact,nact,nact,nact)

        double precision::energy,renergy

        double precision,allocatable::TM1(:,:),TM2(:,:),TM3(:,:)
        double precision,allocatable::TM4(:,:),TM5(:,:),TM6(:,:)
!        double precision,allocatable::GM1(:,:,:,:),GM2(:,:,:,:)
 
        double precision dv


            energy=0.0d0
            renergy=0.0d0 

            allocate(TM1(nact,nact));TM1=0.0d0
            allocate(TM2(nact,nact));TM2=0.0d0
            i0=0
            ia=0
            do i=1,nsub
              ij=occ(i)
              TM1(ia+1:ia+ij,ia+1:ia+ij)=T(i0+1:i0+ij,i0+1:i0+ij)
              i0=i0+tot(i)
              ia=ia+occ(i)
            end do
            call MXM(nact,TM1,D,TM2)
            call trace(nact,TM2,dv)
            energy=dv
            deallocate(TM1)
            deallocate(TM2)

            write(2,*)"The 1-body energy",dv
            call flush(2)

            k0=0
            ka=0
            do k=1,nsub; do kk=1,occ(k)
              l0=0
              la=0
              do l=1,nsub; do ll=1,occ(l)

! \sum_trace(J^{kl}*P^{lk})
                allocate(TM1(norb,norb));TM1=0.0d0
                allocate(TM2(nact,nact));TM2=0.0d0
                allocate(TM3(nact,nact));TM3=0.0d0
                allocate(TM4(nact,nact));TM4=0.0d0

                TM1=             U(:,:,kk+k0,ll+l0)  ! bk
                TM3=             P(:,kk+ka,ll+la,:)  ! bk

!                write(113,*)k,"k;",kk,"kk;",l,'l;',ll,"ll;"
!                write(113,*)"TM1"
!                call print_mat(norb,norb,TM1,113)
!                write(113,*)"TM3"
!                call print_mat(nact,nact,TM3,113)

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
                    TM2(ia+1:ia+ij,ja+1:ja+jj)=TM1(i0+1:i0+ij,j0+1:j0+jj)
                    j0=j0+tot(j)
                    ja=ja+occ(j)
                  end do
                  i0=i0+tot(i)
                  ia=ia+occ(i)
                end do
!                write(113,*)"TM2"
!                call print_mat(nact,nact,TM2,113)

                call MXM(nact,TM2,TM3,TM4)
!                write(113,*)"TM4"
!                call print_mat(nact,nact,TM4,113)

                call trace(nact,TM4,dv)
!              write(113,*)"k,kk,l,ll ",k,kk,l,ll," The 2-body energy",dv

                 energy= energy  + dv*0.5d0 
                renergy= renergy + dv

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
 
            write(2,*)"The 2-body contribution", renergy
            call flush(2)

      End Subroutine energy_gen


