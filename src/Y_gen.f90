      Subroutine Dh_gen(nsub,occ,tot,nact,norb,T,D,dH,group)

        integer::nsub,nact,norb
        integer::occ(nsub)
        integer::tot(nsub)
        integer::group(8,8)
        double precision::T(norb,norb)
        double precision::D(nact,nact)
        double precision dH(norb,norb,norb,norb)
        double precision,allocatable::GM1(:,:,:,:),GM2(:,:,:,:)
        double precision dv,dtmp
        write(2,*)"Enter the Dh_gen"
        write(*,*) "nsub",occ(1)
        call flush(2)
        allocate(GM1(norb,norb,norb,norb)) ! the k,ij,l | The i,j can be further reduced
        GM1=0.0d0
        ioffset=0
        ioffset1=0
        call cpu_time(start_timeDhgen)
        do i=1,nsub; 
          do ii=1,occ(i)
            i0=ioffset
            ix=ioffset1
            joffset=0
            joffset1=0
            do j=1,nsub; 
              do jj=1,occ(j)
                j0=joffset
                jx=joffset1
                GM1(i0+ii,j0+jj,:,:)=T*D(ix+ii,jx+jj)*2.0d0
              end do
              joffset=joffset+tot(j)
              joffset1=joffset1+occ(j)
            end do
          end do
          ioffset=ioffset+tot(i)
          ioffset1=ioffset1+occ(i)
        end do
        call cpu_time(end_timeDhgen)
        ! call matgendh(occ,tot,GM1,T,D,nsub,nact,norb)
        elapsed_timeDhgen = end_timeDhgen - start_timeDhgen
        write(*, *) "Dhgen_time: ", elapsed_timeDhgen, " 秒"
        Dh=GM1
        deallocate(GM1)
      End Subroutine Dh_gen 


      Subroutine Y_gen(nsub,occ,tot,nact,norb,U,P,Y,group)

        ! Use the same way as Molcas/<Molecular Electronic-Structure theory>
        !   in order to check the SA-Hessian issue
        use omp_lib
        use date_time

        integer::nsub,nact,norb
        integer::occ(nsub)
        integer::tot(nsub)
        integer::group(8,8)

        double precision::U(norb,norb,norb,norb)
        double precision::P(nact,nact,nact,nact)

        double precision::Y(norb,norb,norb,norb)

        double precision,allocatable::TM1(:,:),TM2(:,:),TM3(:,:)
        double precision,allocatable::TM4(:,:),TM5(:,:),TM6(:,:)
        double precision,allocatable::GM1(:,:,:,:),GM2(:,:,:,:)

        double precision dv,dtmp
        call cpu_time(start_timeYgen)

        Y=0.0d0         
         
        allocate(GM1(norb,norb,norb,norb))
        GM1=0.0d0

        
        allocate(TM1(norb,norb));
        allocate(TM2(norb,norb));

        ! ioffset=0
        ! ioffset1=0
        ! do i=1,nsub;
        !   do ii=1,occ(i)
        !   i0=ioffset
        !   ix=ioffset1

        !   joffset=0
        !   joffset1=0
        !   do j=1,nsub; do jj=1,occ(j)
        !     j0=joffset
        !     jx=joffset1
        !     TM1=0.0d0
        !     TM2=0.0d0
        !     koffset=0
        !     koffset1=0
        !     do k=1,nsub; do kk=1,occ(k)
        !       k0=koffset
        !       kx=koffset1
        !       loffset=0
        !       loffset1=0
        !       do l=1,nsub
        !         if(group(i,j).eq.group(k,l))then
        !           do ll=1,occ(l)
        !             l0=loffset
        !             lx=loffset1
        !             dtmp=P(ix+ii,kx+kk,lx+ll,jx+jj) ! ij,kl
        !             TM1=TM1+U(:,:,k0+kk,l0+ll)*dtmp
        !           end do
        !         end if
        !         if(group(i,k).eq.group(j,l))then
        !           do ll=1,occ(l)
        !             l0=loffset
        !             lx=loffset1
        !             dtmp=P(ix+ii,jx+jj,lx+ll,kx+kk)& ! ik,jl 
        !                 +P(ix+ii,lx+ll,jx+jj,kx+kk)  ! ik,lj
        !             TM2=TM2+U(:,k0+kk,l0+ll,:)*dtmp
        !           end do
        !         end if
        !         loffset=loffset+tot(l)
        !         loffset1=loffset1+occ(l)
        !       end do
        !     end do
        !     koffset=koffset+tot(k)
        !     koffset1=koffset1+occ(k)
        !     end do
        !     GM1(i0+ii,:,j0+jj,:)=(TM1+TM2)*2.0d0
        !   end do
        !   joffset=joffset+tot(j)
        !   joffset1=joffset1+occ(j)
        !   end do
        ! end do
        ! ioffset=ioffset+tot(i)
        ! ioffset1=ioffset1+occ(i)
        ! end do

        walltime(50) = wtime()
        ioffset=0
        ioffset1=0
        do i=1,nsub
          joffset=0
          joffset1=0
          do j=1,nsub
            koffset=0
            koffset1=0
            do k=1,nsub
              loffset=0
              loffset1=0
              do l=1,nsub
                !$OMP PARALLEL DEFAULT(shared) PRIVATE(ii, jj, kk, ll, TM1, TM2, dtmp)
                !$OMP DO
                do ii=1,occ(i)
                  !$OMP PARALLEL DO
                  do jj=1,occ(j)
                    TM1=0.0d0
                    TM2=0.0d0
                    do kk=1,occ(k)
                      if(group(i,j).eq.group(k,l))then
                        do ll=1,occ(l)
                          dtmp=P(ioffset1+ii,koffset1+kk,loffset1+ll,joffset1+jj) ! ij,kl
                          TM1=TM1+U(:,:,koffset+kk,loffset+ll)*dtmp
                        end do
                      end if
                      if(group(i,k).eq.group(j,l))then
                        do ll=1,occ(l)
                          dtmp=P(ioffset1+ii,joffset1+jj,loffset1+ll,koffset1+kk)& ! ik,jl
                              +P(ioffset1+ii,loffset1+ll,joffset1+jj,koffset1+kk)  ! ik,lj
                          TM2=TM2+U(:,koffset+kk,loffset+ll,:)*dtmp
                        end do
                      end if
                    end do
                    GM1(ioffset+ii,:,joffset+jj,:)=(TM1+TM2)*2.0d0
                  end do
                end do
                !$OMP END PARALLEL
                loffset=loffset+tot(l)
                loffset1=loffset1+occ(l)
              end do
              koffset=koffset+tot(k)
              koffset1=koffset1+occ(k)
              
            end do
            joffset=joffset+tot(j)
            joffset1=joffset1+occ(j)
          end do
          ioffset=ioffset+tot(i)
          ioffset1=ioffset1+occ(i)
        end do

        walltime(51) = wtime()
        write(*,*)"*******************Y_time*****",walltime(51)-walltime(50)



        ! noffset=0
        ! noffset1=0
        ! do n=1,nsub
        !   moffset=0
        !   moffset1=0 
        !   do m=1,nsub
        !     kooffset=0
        !     kooffset1=0
        !     do ko=1,nsub
        !       joffset=0
        !       joffset1=0
        !       do j=1,nsub
        !         koffset=0
        !         koffset1=0
        !         do k=1,nsub
        !           loffset=0
        !           loffset1=0
        !           do l=1,nsub
        !             do n1=1,tot(n)
        !               do m1=1,occ(m)               
        !                 do ko1=1,tot(ko)
        !                   do j1=1,occ(j)
        !                     do k1=1,occ(k)
        !                       do l1=1,occ(l)
        !                         nr=n1+noffset1
        !                         mr=m1+moffset1
        !                         kor=ko1+kooffset1
        !                         jr=j1+joffset1
        !                         kr=k1+koffset1
        !                         lr=l1+loffset1

        !                         ni=n1+noffset
        !                         mi=m1+moffset
        !                         koi=ko1+kooffset
        !                         ji=j1+joffset
        !                         ki=k1+koffset
        !                         li=l1+loffset

        !                         GM1(mi,koi,ji,ni)=GM1(mi,koi,ji,ni)+P(mr,kr,lr,jr)*U(ni,koi,ki,li)&
        !                         +(P(mr,jr,lr,kr)+P(mr,lr,jr,kr))*U(ni,ki,li,koi)

        !                       end do
        !                     end do
        !                   end do
        !                 end do
        !               end do
        !             end do
        !             loffset=loffset+tot(l)
        !             loffset1=loffset1+occ(l)
        !           end do
        !           koffset=koffset+tot(k)
        !           koffset1=koffset1+occ(k)
        !         end do
        !         joffset=joffset+tot(j)
        !         joffset1=joffset1+occ(j)
        !       end do
        !       kooffset=kooffset+tot(ko)
        !       kooffset1=kooffset1+occ(ko)
        !     end do
        !     moffset=moffset+tot(m)
        !     moffset1=moffset1+occ(m)
        !   end do
        !   noffset=noffset+tot(n)
        !   noffset1=noffset1+occ(n)
        ! end do





        call cpu_time(end_timeYgen)
        ! call matgendh(occ,tot,GM1,T,D,nsub,nact,norb)
        elapsed_timeYgen = end_timeYgen - start_timeYgen
        ! write(*, *) "Ygen_time: ", elapsed_timeYgen, " 秒"
        Y=GM1
        deallocate(TM1)
        deallocate(GM1)
        !        call print_gat(norb,norb,norb,norb,Y,318) 

      End Subroutine Y_gen

