
      Subroutine PmultU&
                 (nsub,occ,tot,nact,norb,U,P,Y,group,&
                  isym,ith,jsym,jth,itype)

! Use the same way as Molcas (step-by-step basing on rotation-pairs)
!   in order to check the SA-Hessian issue

! temporary only for C1 symmetry

        integer::nsub,nact,norb
        integer::occ(nsub)
        integer::tot(nsub)
        integer::group(8,8)
        integer::isym,ith,jsym,jth,itype

        double precision::U(norb,norb,norb,norb)
        double precision::P(nact,nact,nact,nact)

        double precision::Y(norb,norb)

        double precision,allocatable::  LM1(:),  LM2(:)
        double precision,allocatable::TM1(:,:),TM2(:,:),TM3(:,:)
        double precision,allocatable::TM4(:,:),TM5(:,:),TM6(:,:)
        double precision,allocatable::GM1(:,:,:,:),GM2(:,:,:,:)
        double precision,allocatable::GM3(:,:,:,:)

        double precision dv,dtmp,fact

!        Y=0.0d0
!*                                                                      *
!************************************************************************
!*                                                                      *
        itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
!*                                                                      *
!************************************************************************
!*                                                                      *
        ndens1=  nact*(  nact+1)/2
        ndens2=ndens1*(ndens1+1)/2

        i0=0
        i1=0
        do i=0,isym-1
          i0=i0+tot(i)
          i1=i1+occ(i)
        end do
        iorb=ith+i0  ! The active orbital index
        iact=ith+i1  ! The active orbital index in RDMs

        j0=0
        j1=0
        do j=0,jsym-1
          j0=j0+tot(j)
          j1=j1+occ(j)
        end do
        jorb=jth+j0  ! The active orbital index
        jact=jth+j1  ! The active orbital index in RDMs


        if(itype.eq.1)then
 
          do i=1,occ(isym)
            do j=1,occ(isym)

!              write(2,*)"i",i,"j",j

              allocate(TM1(norb,norb));TM1=0.0d0
              TM1=U(:,:,i+i0,j+j0)     
!              write(2,*)"The square 2-integrals" 
!              call print_mat(norb,norb,TM1,2)
         
              allocate(TM2(nact,nact));TM2=0.0d0
              TM2=P(iact,:,:,jact)
!              write(2,*)"The square 2-RDMs if itype.eq.1" 
!              call print_mat(nact,nact,TM2,2)

              Y=Y+TM1*TM2(i,j)*2.0d0
              deallocate(TM2)
              deallocate(TM1)

            end do 
          end do 


        else if(itype.eq.2)then

! The j and k packing, inorder to match Molcas - I        
          allocate(GM1(nact,nact,nact,nact))
          allocate(GM2(nact,nact,nact,nact))
          allocate(GM3(nact,nact,nact,nact))
          GM1=0.0d0
          GM2=0.0d0
          GM3=0.0d0
          GM1=P

          do i=1,nact
            do j=1,nact
              do k=1,j
                do l=1,nact
                  if(.true.)then ! If molcas' j,k packing
                    if(j.eq.k)then
                      GM2(i,j,k,l)=GM1(i,j,k,l)
                    else
                      GM2(i,j,k,l)=GM1(i,j,k,l)+GM1(i,k,j,l)
                    end if
                  else            ! if no packing
                    GM2(i,j,k,l)=GM1(i,j,k,l)
                  end if
                end do
              end do
            end do
          end do

          allocate(LM1(ndens2))
          allocate(LM2(ndens2))
          LM1=0.0d0
          LM2=0.0d0

          ij=0
          ijkl=1
          do i=1,nact
            do j=1,i
              ij=ij+1
              kl=0
              do k=1,nact
                do l=1,k
                  kl=kl+1
                  if(ij.ge.kl)then
                    LM1(ijkl)=GM2(i,k,l,j)
!                    write(6,*)ijkl,"ijkl",LM1(ijkl),"LM1(ijkl)"
!                    call flush(6)
                    ijkl=ijkl+1
                  end if
                end do
              end do
            end do
          end do
        
          do i=1,ndenS2
            write(2,*)i,"LM1",LM1(i)
          end do 

          Do iB=1,nact
            Do jB=1,iB
              iDij=iTri(ib,jB)
              Do kB=1,ib
                Do lB=1,kB
                  iDkl=iTri(kB,lB)
                  fact=1.0d00
                  if(iDij.ge.iDkl .and. kB.eq.lB) fact=2.0d00
                  if(iDij.lt.iDkl .and. iB.eq.jB) fact=2.0d00
                  iijkl=itri(iDij,iDkl)
                  LM2(iijkl)=Fact*LM1(iijkl)
!                  write(6,*)iijkl,"iijkl",LM1(iijkl),  &
!                           "LM1(iijkl-1)",LM2(iijkl),"LM2(iijkl)"
                  GM3(iB,jB,kB,lB)=LM2(iijkl)
                  GM3(kB,lB,iB,jB)=LM2(iijkl)
                  GM3(jB,iB,lB,kB)=LM2(iijkl)
                  GM3(lB,kB,iB,jB)=LM2(iijkl)
                End Do
              End Do
            End Do
          End Do

!          do i=1,ndenS2
!            write(2,*)i,"LM2",LM2(i)
!          end do 
 
!          write(2,*)" The packed(RASSCF) rdm-2 tensor " 
!          call print_gat(nact,nact,nact,nact,GM2,2)

          if(iact.eq.jact)then
            do i=1,nact
              do j=i+1,nact
                do k=1,nact 
                  do l=1,k
                  ! GM3(l,i,k,j)=GM3(k,i,l,j)
                  ! if(k.lt.l)then
                  !GM3(k,i,l,j)=GM3(l,i,k,j)
                  GM3(k,j,l,i)=GM3(l,j,k,i) !origin
                  !GM3(l,i,k,j)=GM3(k,j,l,i)  !
                  !GM3(k,i,l,j)=0
                  !GM3(l,i,k,j)=0
                  ! end if 
                  end do
                end do
              end do 
            end do
          end if

          write(2,*)" The packed  (MCLR) rdm-2 tensor " 
          call print_gat(nact,nact,nact,nact,GM3,2)

          deallocate(GM1)
          deallocate(GM2)

!          stop

          do i=1,occ(isym)
            do j=1,occ(jsym)

              write(2,*)"i",i,"j",j
              write(2,*)"iact",iact,"jact",jact

              allocate(TM1(norb,norb));TM1=0.0d0
              allocate(TM3(norb,norb));TM3=0.0d0
              TM1=U(:,i+i0,j+j0,:)
 
              ! 
              do ii=1,norb
                do jj=1,norb
                  TM3(ii,jj)=TM1(jj,ii)
                end do
              end do  

              write(2,*)"The square 2-integrals"
              call print_mat(norb,norb,TM1,2)

              allocate(TM2(nact,nact));TM2=0.0d0
              !TM2=P(iact,iact,:,:)
              !TM2=GM3(iact,iact,:,:)
              TM2=GM3(:,iact,:,jact)
              write(2,*)"The square 2-RDMs if itype.eq.2"
              call print_mat(nact,nact,TM2,2)

              Y=Y+TM1*TM2(i,j)*2.0d0
!              if(i.ne.j)then
!                if(iact.ne.jact)then
!                  Y=Y+TM3*TM2(i,j)*2.0d0
!                end if
!              end if  
 
              deallocate(TM3)
              deallocate(TM2)
              deallocate(TM1)
 
              write(2,*)" --- The squared PmultU ---" 
              call print_mat(norb,norb,Y,2)

            end do
          end do

          deallocate(GM3) 

!          Do Ks=1,nsym
!            iOpt=1
!            JLB=1
!            JLBas=0
!            ijkl=nOrb(js)*nash(ks)
!            If (ijkl.ne.0) Then
!
!              jlBB=0
!              Do LB=nish(ks)+1,nB(KS)
!                kkc=nA(ks)+lb-nish(ks)
!               Do JB=nish(ks)+1,nB(KS)
!                 kkb=nA(ks)+jb-nish(ks)
!                 Call EXCH(js,ks,js,ks,jb,lb,Temp1,Scr)
!                 ipT=1
!                 If (LB.gt.nISH(ks).and.jb.gt.nish(ks)) Then
!                   rDens2=sign*4.0d0*Work(ipG2-1+
!    &              itri(itri(iib,kkc),itri(kkb,iib)))
!
!                   write(6,*)"LB,JB",LB,JB," rdens2 ", rdens2   ! yma
!                   call flush(6)
!
!                   do i=1,nO**2
!            write(6,*)i," The 2nd round/start, Temp1 (loop) ", Temp1(i)
!                   end do
!                   do i=1,nO     ! yma
!                     do j=1,nO
!                     write(6,*)i,j," The 2nd, Temp2 (loop) ",temp2(i,j)
!                     end do
!                   end do
!                   call flush(6)
!
!                   Call DaXpY_(nO**2,rDens2,Temp1(ipT),1,Temp2,1)

!                   do i=1,nO**2
!            write(6,*)i," The 2nd round/after, Temp1 (loop) ", Temp1(i)
!                   end do
!                   do i=1,nO     ! yma
!                     do j=1,nO
!                       write(6,*)i,j," The 2nd, Temp2 (loop) ",temp2(i,j)
!                     end do
!                   end do
!                   call flush(6)
!
!                 End If
!               End Do
!             End Do
!           End If
!         End Do

 
        else
          write(6,*)"Unknow PmultU type"
        end if 

!        write(2,*)"The squared PmultU " 
!        call print_mat(norb,norb,Y,2)
!        stop

      End Subroutine PmultU

