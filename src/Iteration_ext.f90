! The extended Iteration, which is origional from Iteration
! The "extended" means all the operators will be updated 
!      using the 1-index during the iterations
! 

      Subroutine Iteration_ext(icycle,pre_solved_internal)

        use global_control
        use matrix
 
        integer::icycle
        logical::pre_solved_internal
        double precision d1,d2,d3,d4,d5,R_threshold
        double precision iter_threshold
        double precision,allocatable::deltR(:,:)

        type::iter
          double precision R20,R21,Rsub,Rmax
          double precision,allocatable::Edn(:,:) 
          double precision,allocatable::Edn1(:,:) 
          double precision,allocatable::deltR(:,:) 
          double precision,allocatable::A_Ap(:,:) 
          double precision,allocatable::Rn(:,:) 
          double precision,allocatable::Rn1(:,:) 
          double precision,allocatable::Rnp(:,:) 
          double precision,allocatable::Dn(:,:) 
          double precision,allocatable::Dn1(:,:) 
        end type iter
        type(iter),allocatable::mat3(:)

        double precision,allocatable::Tsv(:,:),Usv(:,:,:,:)

        allocate(mat3(orb%nsub))

        do i=1,orb%nsub
          allocate(mat3(i)%Edn(orb%total(i),orb%total(i)))
          mat3(i)%Edn=0.0d0
          allocate(mat3(i)%Edn1(orb%total(i),orb%total(i)))
          mat3(i)%Edn1=0.0d0
          allocate(mat3(i)%deltR(orb%total(i),orb%total(i)))
          mat3(i)%deltR=0.0d0
          allocate(mat3(i)%A_Ap(orb%total(i),orb%total(i)))
          mat3(i)%A_Ap=0.0d0
          allocate(mat3(i)%Rn(orb%total(i),orb%total(i)))
          mat3(i)%Rn=0.0d0
          allocate(mat3(i)%Rn1(orb%total(i),orb%total(i)))
          mat3(i)%Rn1=0.0d0
          allocate(mat3(i)%Rnp(orb%total(i),orb%total(i)))
          mat3(i)%Rnp=0.0d0
          allocate(mat3(i)%Dn(orb%total(i),orb%total(i)))
          mat3(i)%Dn=0.0d0
          allocate(mat3(i)%Dn1(orb%total(i),orb%total(i)))
          mat3(i)%Dn1=0.0d0
        end do 

        !write(*,*)"++++++++++"

!       Sum|delta_R| threshold 
        if(icycle.eq.1)then
          iter_threshold=0.01d0 
        else if(icycle.eq.2)then
          if(threshold(1).gt.0.01d0)then
            iter_threshold=0.01d0
          else
            iter_threshold=0.001d0
          end if 
        else if(icycle.eq.3)then
          if(threshold(icycle-1).lt.threshold(icycle-2))then
            iter_threshold=threshold(icycle-1)/10.0d0
          else
            iter_threshold=threshold(icycle-1)
          end if
        else
          if(threshold(icycle-1).lt.threshold(icycle-3))then
            if(threshold(icycle-1).lt.threshold(icycle-2))then
              iter_threshold=threshold(icycle-1)/100.0d0
            else
              iter_threshold=threshold(icycle-1)
            end if
          else
            if(threshold(icycle-1).lt.threshold(icycle-2))then
              iter_threshold=threshold(icycle-1)/10.0d0
            else
              iter_threshold=threshold(icycle-1)
            end if
          end if
          if(iter_threshold.lt.1.0e-8)then
            iter_threshold=1.0e-8
          end if
        end if 

!       R_threshold
        if(norb.lt.30)then
          R_threshold=5.0d0
        else if(norb.lt.60)then
          R_threshold=10.0d0
        else if(norb.lt.120)then
          R_threshold=20.0d0
        else 
          R_threshold=30.0d0
        end if

!        iter_threshold=1.0e-5
!        R_threshold=40.d0

        allocate(deltR(norb,norb))
        deltR=0.0d0 
!        write(*,*)"In iteration",orb%closed(1)
!        call print_mat(norb,norb,mat1(1)%R)
!        call print_mat(norb,norb,mat1(1)%Hdiag)
        !
        ioffset=0
        do i=1,orb%nsub
          mat3(i)%Dn=mat1(i)%Hdiag
          mat3(i)%Rn=mat1(i)%R
          do j=1,orb%total(i)
            do k=1,orb%total(i) 
              if(j.ne.k)then 
                if(j.gt.orb%occ(i).and.k.gt.orb%occ(i))then
                else if(j.le.orb%closed(i).and.k.le.orb%closed(i))then
                else
                  if(icycle.gt.0)then
                    mat3(i)%A_Ap(j,k)=mat1(i)%A(j,k)-mat1(i)%A(k,j)
                    if(dabs(mat3(i)%A_Ap(j,k)).lt.1.0d-9)then
                      mat3(i)%A_Ap(j,k)=0.0d0
                    end if 
                  else  ! this should not happen, why it is here?  
                    mat3(i)%A_Ap(j,k)=0.0d0 
                  end if
!                  write(*,*)j,k,mat3(i)%A_Ap(j,k)
                  mat3(i)%Rn(j,k)=-1.0d0*mat3(i)%A_Ap(j,k)/mat3(i)%Dn(j,k)
                end if
              end if
              if(pre_solved_internal)then  ! partly inversed Hessian was used
                if(j.le.orb%occ(i).and.k.le.orb%occ(i))then
                  mat3(i)%Rn(j,k)=0.0d0
                  mat3(i)%A_Ap(j,k)=0.0d0 
                end if
              end if 
              mat2%R(j+ioffset,k+ioffset)=mat3(i)%Rn(j,k)
            end do
          end do
          ioffset=ioffset+orb%total(i)
!          call print_mat(orb%total(i),orb%total(i),mat3(i)%A_Ap,6)
!          call print_mat(orb%total(i),orb%total(i),mat3(i)%Rn,6)
          !mat3(i)%Rn=0.0d0 
        end do
 
        d3=0.0d0
        do i=1,orb%nsub
          mat3(i)%R20=0.0d0
          do j=1,orb%total(i)
            do k=1,orb%total(i)
              mat3(i)%R20=mat3(i)%R20+dabs(mat3(i)%Rn(j,k))
            end do   
          end do
          d3=d3+mat3(i)%R20 
        end do
        do while(dabs(d3).gt.sqrt(norb/1.0d0)) 
          d3=d3/2.0d0 
          do i=1,orb%nsub
            mat3(i)%Rn=mat3(i)%Rn/2.0d0
          end do 
        end do        
        ioffset=0
        do i=1,orb%nsub
          do j=1,orb%total(i)
            do k=1,orb%total(i) 
              mat2%R(j+ioffset,k+ioffset)=mat3(i)%Rn(j,k)
            end do
          end do
          ioffset=ioffset+orb%total(i)
        end do

!        write(*,*)" === mat2%R === "
!        call print_MAT(norb,norb,mat2%R) 
!        stop

        ! At first, save the currnt T and U  
        allocate(Tsv(norb,norb))
        allocate(Usv(norb,norb,norb,norb))
        Tsv=T
        Usv=U

        iloop=0
        do 
          iloop=iloop+1

          ! R-> T-> U
          d3=0.0d0
          do i=1,orb%nsub
            mat3(i)%R20=0.0d0
            do j=1,orb%total(i)
              do k=1,orb%total(i)
                mat3(i)%R20=mat3(i)%R20+dabs(mat3(i)%Rn(j,k))
              end do   
            end do
            d3=d3+mat3(i)%R20 
            !call CAL_T(orb%total(i),mat3(i)%Rn,mat3(i)%R20,150,mat1(i)%T) 
            !call T_TO_U(orb%total(i),mat1(i)%T,mat1(i)%U)
          end do

          call CAL_T(norb,mat2%R,d3,150,mat2%T)
          call T_TO_U(norb,mat2%T,mat2%U)

          if(mod(iloop,10).eq.5)then
          ! E'=UB-BU 
            !write(*,*)"=======",method 
            call transform_1index&
                 (orb%nsub,orb%act,orb%total,nact,norb,&
                  MAT2%U,T,U,T1dx,U1dx,orb%grouptable)
            call matrix_copy(norb,norb,T1dx,T)
            call tensor_copy(norb,norb,norb,norb,U1dx,U) 
            !call transform_INT(.false.)  ! Thus only generate new FCIDUMP and FCIDUMP_ACTIVE
            !call deallocate_preonescf()
            !call PreOneSCF(icycle)
            !call operator_update() !rewrite
            !write(6,*)"==  After new working equation  =====",icycle
            !call flush(6)
            !call Initial_value()
            call derivaties(0)
            !write(6,*)"==  After derivatie  =====",icycle
            !call flush(6)
          end if

          do i=1,orb%nsub
            mat3(i)%Edn=mat1(i)%Ederi
            !do j=1,orb%total(i)
             ! do k=1,orb%total(i)
                !write(*,*)mat1(i)%Ederi(j,k),mat3(i)%Rn(j,k)
              !end do
            !end do
          end do

          !stop
          !
          !
          if(iloop.gt.2)then
            do i=1,orb%nsub
              mat3(i)%Dn=0.0d0
              do j=1,orb%total(i)
                do k=1,orb%total(i)         
                  if(j.ne.k)then 
                    if(j.gt.orb%occ(i).and.k.gt.orb%occ(i))then
                    else if(j.le.orb%occ(i).and.k.le.orb%occ(i)&
                                    .and.act_omit.eqv..true. )then
!!!                    else if(j.le.orb%closed(i).and.k.le.orb%closed(i))then
                    else if(dabs(mat3(i)%A_Ap(j,k)).lt.1.0d-9)then !only this should work
                    else
                      d1=mat3(i)%Edn(j,k)-mat3(i)%Edn1(j,k)
                      d2=mat3(i)%Rn(j,k)-mat3(i)%Rn1(j,k)
                      mat3(i)%Dn(j,k)=d1/d2 
                      !write(*,*)j,k,mat3(i)%Edn(j,k),mat3(i)%Edn1(j,k)
                      !write(*,*)j,k,mat3(i)%Rn(j,k),mat3(i)%Rn1(j,k)
                      !write(*,*)"D",mat3(i)%Dn(j,k)
                      if(mat3(i)%Dn(j,k).lt.0.67d0*mat3(i)%Dn1(j,k))then
                        mat3(i)%Dn(j,k)=mat3(i)%Dn1(j,k)*0.67d0
                      else if(mat3(i)%Dn(j,k).gt.1.5d0*mat3(i)%Dn1(j,k))then
                        mat3(i)%Dn(j,k)=mat3(i)%Dn1(j,k)*1.5d0
                      else
                        mat3(i)%Dn(j,k)=mat3(i)%Dn1(j,k)       
                      end if
                      !if(mat3(i)%Dn(j,k).lt.0)then
                      !  mat3(i)%Dn(j,k)=-0.1d0*mat3(i)%Dn(j,k)
                      !end if
!                      write(708,*)i,j,k,"D",mat3(i)%Dn(j,k)
                    end if
                  end if
                end do  
              end do    
              !mat3(i)%Rn1=mat3(i)%Rn
            end do 
            !stop
          end if

          !write(*,*)"+++++++++++++++++++"
 
          ioffset=0
          do i=1,orb%nsub
            do j=1,orb%total(i)
              do k=1,orb%total(i)
                if(j.ne.k)then 
                  if(j.gt.orb%occ(i).and.k.gt.orb%occ(i))then
                  else if(j.le.orb%occ(i).and.k.le.orb%occ(i)&
                                    .and.act_omit.eqv..true. )then
!!!                  else if(j.le.orb%closed(i).and.k.le.orb%closed(i))then
                  else if(dabs(mat3(i)%A_Ap(j,k)).lt.1.0d-9)then !only this should work
                  else
                    mat3(i)%deltR(j,k)=-1.0d0*mat3(i)%Edn(j,k)/mat3(i)%Dn(j,k)
!                    write(705,*)i,j,k,mat3(i)%Edn(j,k),mat3(i)%Dn(j,k)
                  end if 
                end if
                deltR(j+ioffset,k+ioffset)=mat3(i)%deltR(j,k)
              end do
            end do
            ioffset=ioffset+orb%total(i)
          end do

!=====================================================
          

          ioffset=0
          d1=0.0d0
          do i=1,orb%nsub
            mat3(i)%Rsub=0.0d0
            do j=1,orb%total(i)
              do k=1,orb%total(i)
                mat3(i)%Rsub=mat3(i)%Rsub+dabs(deltR(j+ioffset,k+ioffset))
                if(deltR(j+ioffset,k+ioffset).gt.0.1d0)then
                  deltR(j+ioffset,k+ioffset)=0.1d0
                else if(deltR(j+ioffset,k+ioffset).lt.-0.1d0)then
                  deltR(j+ioffset,k+ioffset)=-0.1d0
                end if
                mat3(i)%deltR(j,k)=deltR(j+ioffset,k+ioffset)
              end do
            end do
            ioffset=ioffset+orb%total(i)
            d1=d1+mat3(i)%Rsub   
          end do

          !do i=1,orb%nsub
          !!  do while(mat3(i)%Rsub.gt.4.0d0)
          !    mat3(i)%Rsub=mat3(i)%Rsub/2.0d0
          !    mat3(i)%deltR=mat3(i)%deltR/2.0d0
          !  end do
          !end do

          ioffset=0
          do i=1,orb%nsub
            do j=1,orb%total(i)
              do k=1,orb%total(i)
                deltR(j+ioffset,k+ioffset)=mat3(i)%deltR(j,k)
                mat2%deltR(j+ioffset,k+ioffset)=mat3(i)%deltR(j,k)
              end do
            end do 
            ioffset=ioffset+orb%total(i)
          end do
         
!====================================================

          !  
          ioffset=0
          do i=1,orb%nsub         
            !mat3(i)%Rnp=mat3(i)%Rn+mat3(i)%deltR
            do j=1,orb%total(i)
              do k=1,orb%total(i)
                if(j.ne.k)then
                  if(j.gt.orb%occ(i).and.k.gt.orb%occ(i))then
!!!                  else if(j.le.orb%closed(i).and.k.le.orb%closed(i))then
                  else
                    if(mat3(i)%Dn(j,k).gt.0.0d0)then
                      mat3(i)%Rnp(j,k)=mat3(i)%Rn(j,k)+mat3(i)%deltR(j,k) 
                    else if(mat3(i)%Dn(j,k).lt.0.0d0)then
                      mat3(i)%Rnp(j,k)=mat3(i)%Rn(j,k)-mat3(i)%deltR(j,k)
                    else
                      mat3(i)%Rnp(j,k)=mat3(i)%Rn(j,k) 
                    end if
                  end if
                end if 
                !write(*,*)j,k,mat3(i)%Rn(j,k),mat3(i)%deltR(j,k)
                mat2%R(j+ioffset,k+ioffset)=mat3(i)%Rnp(j,k)
!                write(709,609)i,j,k,mat3(i)%Dn(j,k),mat3(i)%Rnp(j,k),mat3(i)%deltR(j,k)
              end do
            end do
            mat3(i)%Rn1=mat3(i)%Rn
            mat3(i)%Rn=mat3(i)%Rnp
            mat3(i)%Edn1=mat3(i)%Edn
            mat3(i)%Dn1=mat3(i)%Dn
            ioffset=ioffset+orb%total(i)
          end do
609       format(I3,I3,I3,f8.4,f8.4,f8.4)

          if(pre_solved_internal)then
            ioffset=0
            do i=1,orb%nsub         
              do j=1,orb%occ(i)
                do k=1,orb%occ(i)
                  mat2%R(j+ioffset,k+ioffset)=mat2%Rocc(j,k)
                end do 
              end do
              ioffset=ioffset+orb%total(i)
            end do
          end if

!          write(*,*)"========================" 
!          call print_mat(norb,norb,mat2%R,6) 

          ! R-> T-> U
          d4=0.0d0
          do i=1,orb%nsub
            mat3(i)%R21=0.0d0
            do j=1,orb%total(i)
              do k=1,orb%total(i)
                mat3(i)%R21=mat3(i)%R21+dabs(mat3(i)%Rn(j,k))
              end do   
            end do
            d4=d4+mat3(i)%R21
          end do

          !write(*,*)iloop,d3,d4,d3-d4

          if(dabs(d3-d4).lt.iter_threshold)then
            threshold(icycle)=dabs(d3-d4)
            Rabs(icycle)=dabs(d4)
            exit 
          end if
                      
          if(dabs(d4).gt.R_threshold)then
            threshold(icycle)=dabs(d3-d4)
            Rabs(icycle)=dabs(d4)
            exit
          end if    
          if(iloop.gt.150)then
            threshold(icycle)=dabs(d3-d4)
            Rabs(icycle)=dabs(d4)
            exit
          end if 

        end do

        T=Tsv
        U=Usv
        deallocate(Tsv,Usv)

        write(*,2435,advance='no')"    ",Rabs(icycle)  
        write(*,2436)"   ",method(icycle)  
2435    format(A4,f8.4)
2436    format(A3,A7)

        do i=1,orb%nsub
          mat1(i)%R=mat3(i)%Rn
          do j=1,orb%total(i)
            do k=1,orb%total(i)
              write(701,*)j,k,mat1(i)%R(j,k)
            end do
          end do   
        end do

!!        call GRAM_SCHMIDT(norb,mat2%U)

!        do i=1,norb
!          do j=1,norb
!            write(702,*)i,j,mat2%U(i,j)
!          end do
!        end do 

        deallocate(deltR)
        deallocate(mat3)

      end Subroutine Iteration_ext

