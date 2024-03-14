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
        allocate(TM1(norb,norb));
        allocate(TM2(norb,norb));
        allocate(TM3(norb,norb));

        ! write(*,*) U(:,:,1,1)

        ! GM1=0.0d0
        ! ioffset=0
        ! ioffset1=0
        
        ! do i=1,nsub; 
        !   do ii=1,occ(i)
        !   i0=ioffset
        !   i1=tot(i)
        !   i2=occ(i)
        !   ix=ioffset1

        !   joffset=0
        !   joffset1=0
        !   do j=1,nsub; do jj=1,occ(j)
        !     j0=joffset
        !     j1=tot(j)
        !     j2=occ(j)
        !     jx=joffset1
        !     TM1=0.0d0
        !     TM2=0.0d0
        !     TM3=0.0d0
        !     koffset=0
        !     koffset1=0
        !     do k=1,nsub; do kk=1,occ(k)
        !       k0=koffset
        !       k1=tot(k)
        !       k2=occ(k)
        !       kx=koffset1
        !       loffset=0
        !       loffset1=0
        !       do l=1,nsub
        !         if(group(i,j).eq.group(k,l))then
        !           do ll=1,occ(l)
        !             l0=loffset
        !             l1=tot(l)
        !             l2=occ(l)
        !             lx=loffset1
        !             dtmp=P(ix+ii,kx+kk,lx+ll,jx+jj) ! ij,kl
        !             TM1=TM1+U(:,:,k0+kk,l0+ll)*dtmp
        !             ! if(ii==1 .and. jj==1 .and. kk==1 .and.ll==1)then
        !             !   ! write(*,*) "U(:,:,k0+kk,l0+ll)" ,U(:,:,k0+kk,l0+ll)
        !             !   ! write(*,*) U(:,1,k0+kk,l0+ll)
        !             !   write(*,*) TM1(:,1)
        !             ! end if
        !           end do
        !         end if
        !         if(group(i,k).eq.group(j,l))then
        !           do ll=1,occ(l)
        !             l0=loffset
        !             l1=tot(l)
        !             l2=occ(l)
        !             lx=loffset1
        !             dtmp=P(ix+ii,jx+jj,lx+ll,kx+kk)& ! ik,jl 
        !                 +P(ix+ii,lx+ll,jx+jj,kx+kk)  ! ik,lj
        !             TM2=TM2+U(:,k0+kk,l0+ll,:)*dtmp
        !             ! if(ii==1 .and. jj==1 .and. kk==1 .and.ll==1)then
        !             !   ! write(*,*) "U(:,:,k0+kk,l0+ll)" ,U(:,:,k0+kk,l0+ll)
        !             !   ! write(*,*) U(:,1,k0+kk,l0+ll)
        !             !   write(*,*) TM2(:,1)
        !             ! end if
        !           end do
        !           if(ii==1 .and. jj==1 .and. kk==1)then
        !             write(*,*) "*********23"
        !           end if
        !         end if
        !         loffset=loffset+tot(l)
        !         loffset1=loffset1+occ(l)
        !       end do
        !     end do
        !     koffset=koffset+tot(k)
        !     koffset1=koffset1+occ(k)
        !     end do
        !     GM1(i0+ii,:,j0+jj,:)=(TM1+TM2)*2.0d0
        !     if(ii==1 .and. jj==1)then
        !       ! write(*,*) "U(:,:,k0+kk,l0+ll)" ,U(:,:,k0+kk,l0+ll)
        !       ! write(*,*) U(:,1,k0+kk,l0+ll)
        !       write(*,*) GM1(:,1,1,1)
        !     end if
            
        !   end do
        !   joffset=joffset+tot(j)
        !   joffset1=joffset1+occ(j)
        !   end do
        ! end do
        ! ioffset=ioffset+tot(i)
        ! ioffset1=ioffset1+occ(i)
        ! end do
        ! call cpu_time(end_timeYgen)
        call matgeny(occ,tot,GM1,T,D,nsub,nact,norb,group,P,U);
        elapsed_timeYgen = end_timeYgen - start_timeYgen
        write(*, *) "Ygen_time: ", elapsed_timeYgen, " 秒"
        Y=GM1
        deallocate(TM1)
        deallocate(TM2)
        deallocate(TM3)
        deallocate(GM1)

        !        call print_gat(norb,norb,norb,norb,Y,318) 

      End Subroutine Y_gen

