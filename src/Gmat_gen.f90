! to Stefan : <r|G_ij|s> in Eq.25 / Eq.32
!             Eq.32 can share with this formula when U+&U is already operated on MO-integrals

      Subroutine Gmat_gen(nsub,occ,tot,nact,norb,T,U,D,P,Gmat,group)

        integer::nsub,nact,norb
        integer::occ(nsub)
        integer::tot(nsub)
        integer::group(8,8)

        double precision::T(norb,norb)
        double precision::U(norb,norb,norb,norb)

        double precision::D(nact,nact)
        double precision::P(nact,nact,nact,nact)

        double precision::Gmat(norb,norb,norb,norb)

        double precision,allocatable::TM1(:,:),TM2(:,:),TM3(:,:)
        double precision,allocatable::TM4(:,:),TM5(:,:),TM6(:,:)
        double precision,allocatable::GM1(:,:,:,:),GM2(:,:,:,:)

        double precision dv,dtmp 

        Gmat=0.0d0        

        allocate(GM1(norb,norb,norb,norb)) ! the k,ij,l | The i,j can be further reduced
        GM1=0.0d0
        ioffset=0
        ioffset1=0
        do i=1,nsub; do ii=1,occ(i)
          i0=ioffset
          i1=tot(i)
          i2=occ(i)
          ix=ioffset1

          joffset=0
          joffset1=0
          do j=1,nsub; do jj=1,occ(j)
            j0=joffset
            j1=tot(j)
            j2=occ(j)
            jx=joffset1

            GM1(:,i0+ii,j0+jj,:)=T*D(ix+ii,jx+jj)*2.0d0
          end do
          joffset=joffset+tot(j)
          joffset1=joffset1+occ(j)
          end do
        end do
        ioffset=ioffset+tot(i)
        ioffset1=ioffset1+occ(i)
        end do
        Gmat=GM1
!         write(*,*)"mat2%G || GM1" 
!         call print_GAT(norb,norb,norb,norb,GM1,315) 
        deallocate(GM1)


        allocate(GM1(norb,norb,norb,norb))
        GM1=0.0d0

        ioffset=0
        ioffset1=0
        do i=1,nsub; do ii=1,occ(i)
          i0=ioffset
          i1=tot(i)
          i2=occ(i)
          ix=ioffset1

          joffset=0
          joffset1=0
          do j=1,nsub; do jj=1,occ(j)
            j0=joffset
            j1=tot(j)
            j2=occ(j)
            jx=joffset1

            allocate(TM1(norb,norb));TM1=0.0d0

            koffset=0
            koffset1=0
            do k=1,nsub; do kk=1,occ(k)
              k0=koffset
              k1=tot(k)
              k2=occ(k)
              kx=koffset1

              loffset=0
              loffset1=0
              do l=1,nsub
! this is for JP
                if(group(i,j).eq.group(k,l))then
                  do ll=1,occ(l)
                    l0=loffset
                    l1=tot(l)
                    l2=occ(l)
                    lx=loffset1
                    dtmp=P(ix+ii,kx+kk,lx+ll,jx+jj)
                    TM1=TM1+U(:,:,k0+kk,l0+ll)*dtmp*2.0d0
                  end do
                end if

! this is for KQ
                if(group(i,k).eq.group(j,l))then
                  do ll=1,occ(l)
                    l0=loffset
                    l1=tot(l)
                    l2=occ(l)
                    lx=loffset1
                    dtmp=P(ix+ii,jx+jj,lx+ll,kx+kk)
                    TM1=TM1+U(:,k0+kk,l0+ll,:)*dtmp*4.0d0
                  end do
                end if
                loffset=loffset+tot(l)
                loffset1=loffset1+occ(l)
              end do
            end do
            koffset=koffset+tot(k)
            koffset1=koffset1+occ(k)
            end do
            GM1(:,i0+ii,j0+jj,:)=TM1
            deallocate(TM1)

          end do
          joffset=joffset+tot(j)
          joffset1=joffset1+occ(j)
          end do
        end do
        ioffset=ioffset+tot(i)
        ioffset1=ioffset1+occ(i)
        end do

!        call print_gat(norb,norb,norb,norb,GM1,317) 

        Gmat=Gmat+GM1
        deallocate(GM1)        


      End Subroutine Gmat_gen

