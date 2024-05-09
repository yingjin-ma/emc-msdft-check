Subroutine Solver_gmres(icycle)
        !定义
          use global_control
          use matrix
  
          integer::icycle
          double precision,allocatable::x_estimate(:)
          double precision,allocatable::xxx(:)
          double precision,allocatable::xD(:)
          real(4),allocatable::xxxs(:)
          
   
          integer,allocatable::porb(:)
          double precision d0,dv,redundant
          double precision,allocatable::D(:),X(:),Y(:),XP(:),YP(:) 
          double precision,allocatable::H(:,:),HP(:,:),HPA(:,:),HP2(:,:)
  
          double precision,allocatable::T1(:,:),T2(:,:),T3(:,:) 
  
          ! Used for the augmented Hessian part
          double precision lamda   ! damping parameter that control the level shift
  
          integer,allocatable::valid(:,:)
          double precision,allocatable::G(:)
          logical Lvalid
          real(4) , allocatable :: SHP(:,:)
          double precision :: sum_h
          real(4) , allocatable :: xS(:)
          
          
  
          ! Parameter initialization
          lamda=1.0d0 ! "1" means the standard AH 
          ! The larger lamda, the shorter step
  
          ! initialize of part of Hessian matrix that will be used in this solver
  
          ndim1=norb
  
          allocate(X(ndim1**2)); X=0.0d0
          allocate(Y(ndim1**2)); Y=0.0d0
          allocate(G(ndim1**2)); G=0.0d0
          allocate(H(ndim1**2,ndim1**2)); H=0.0d0
  
        !end
        !计算Hessian
          iHessian=1
          if(iHessian.eq.0)then 
            call Hessian(norb,mat2%A,mat2%G,mat2%H,mat2%Hdiag)
          else
            if(.NOT.ALLOCATED(mat2%Dh))then
              allocate(mat2%Dh(norb,norb,norb,norb)) ! use the same form as G
            else
              mat2%Dh=0.0d0 
            end if
            if(.NOT.ALLOCATED(mat2%Y))then
              allocate(mat2%Y(norb,norb,norb,norb)) ! use the same form as G
            else
              mat2%Y=0.0d0 
            end if
            ! Same way as in "Molecular electronic-structure theory" 
            write(2,*)"Dh_gen start"
            call flush(2)
            write(2,*)"orb%nsub" ,  orb%nsub
            call flush(2)
            write(2,*)"orb%act"  ,  orb%act
            call flush(2)
            write(2,*)"orb%total",  orb%total
            call flush(2)
            write(2,*)"nact"     ,     nact
            call flush(2)
            write(2,*)"norb"     ,     norb
            call flush(2)
            call Dh_gen(orb%nsub,orb%act,orb%total,&
                       nact,norb,T,mat2%D,mat2%Dh,orb%grouptable)
            ! call d_h_gen_c(orb%nsub,orb%act,orb%total,nact,norb,T,mat2%D,mat2%Dh,orb%grouptable)
            write(2,*)"Dh_gen finish"
            call flush(2)
            call Y_gen(orb%nsub,orb%act,orb%total,&
                       nact,norb,U,mat2%P,mat2%Y,orb%grouptable)
            write(2,*)"Y_gen finish"
            call flush(2)
            call Hessian3(norb,mat2%Dh,mat2%A,mat2%Y,mat2%H,mat2%Hdiag)
          end if
        !end
        !拉伸Hessian
          ij=0
          do i=1,ndim1
            do j=1,ndim1
              ij=ij+1
              kl=0
              do k=1,ndim1
                do l=1,ndim1
                  kl=kl+1
                  H(ij,kl)=mat2%H(i,j,k,l) 
                end do
              end do
            end do
          end do
        !end
        ! write(*,*)"H(:,1)"
        ! write(*,*),H(:200,1)
        ! write(*,*)"A(:,1)"
        ! write(*,*),mat2%A(:200,1)
        !This is the orbiatl gradients 
          redundant=1.0d-9
          
                !  ij=0
                !  nij=0
                !  do i=1,ndim1
                !    do j=1,ndim1
                !       ij=ij+1
                !       Y(ij)=-1.0d0*(mat2%A(i,j)-mat2%A(j,i))
                !       if(dabs(Y(ij)).lt.redundant)then
                !       else
                !         nij=nij+1
                !       end if
                !       write(6,"(f9.5)",advance='no')Y(ij)
                !    end do
                !    write(6,"(f9.5)")Y(ij)
                !  end do
                !  write(*,*) "Y"
                ! write(*,*),Y(:100)
                ! write(*,*),Y(4)
          ! do j=i+1,ndim1 only half is needed 
            grad_digit=0.0d0
            ij=0
            nij=0
            do i=1,ndim1
              do j=i+1,ndim1 ! only half is needed
              !do j=i,ndim1 ! only half is needed
                 ij=ij+1
                 ! Residual, only A=UB , same as NR, AH if U=1
                 Y(ij)= 1.0d0*(mat2%A(i,j)-mat2%A(j,i))
                 G(ij)= 1.0d0*(mat2%A(i,j)-mat2%A(j,i))
                 if(dabs(Y(ij)).lt.redundant)then
                 else
              !                 write(1,*)"in else"
                   IF(LMCSCF)then
                     ! Check for degeneracy irreps
                     if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                       write(1,*)icycle,"-th macro-loop, ini-grad",i,j,Y(ij)
                       nij=nij+1
                       if(dabs(G(ij)).gt.grad_digit)then
                         grad_digit=dabs(G(ij))
                       end if
                     end if
                     ! end Check for degeneracy irreps
                   else
                     ! no closed shell 
                     if(IABS(redu%orbirp(i)).lt.10.and.IABS(redu%orbirp(j)).lt.10)then
                       if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
              !                     write(1,*)"in else : casscf type"                
              !                     call flush(1)
                       ! Check for degeneracy irreps
                         if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                           write(1,*)icycle,"-th macro-loop, ini-grad",i,j,Y(ij)
                           nij=nij+1
                           if(dabs(G(ij)).gt.grad_digit)then
                             grad_digit=dabs(G(ij))
                           end if
                         end if
                       ! end Check for degeneracy irreps
                       end if
                     ! i is closed shell
                     else if(IABS(redu%orbirp(i)).gt.10.and.IABS(redu%orbirp(j)).lt.10)then
                       !if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                         if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                           write(1,*)icycle,"-th macro-loop, ini-grad",i,j,Y(ij)
                           nij=nij+1
                           if(dabs(G(ij)).gt.grad_digit)then
                             grad_digit=dabs(G(ij))
                           end if
                         end if
                       !end if
                     ! j is closed shell
                     else if(IABS(redu%orbirp(j)).gt.10.and.IABS(redu%orbirp(i)).lt.10)then
                       !if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                         if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                           write(1,*)icycle,"-th macro-loop, ini-grad",i,j,Y(ij)
                           nij=nij+1
                           if(dabs(G(ij)).gt.grad_digit)then
                             grad_digit=dabs(G(ij))
                           end if
                         end if
                       !end if  
                     end if ! close shell
                   end if ! MCSCF or CASSCF
                 end if
              end do
            end do
          ! end do j=i+1,ndim1
          
          
          ! do j=1,i only half is needed
            do i=1,ndim1
              do j=1,i ! only half is needed
              !do j=i,ndim1 ! only half is needed
                 ij=ij+1
                 ! Residual, only A=UB , same as NR, AH if U=1
                 Y(ij)= 1.0d0*(mat2%A(i,j)-mat2%A(j,i))
                 G(ij)= 1.0d0*(mat2%A(i,j)-mat2%A(j,i))
                 if(dabs(Y(ij)).lt.redundant)then
                 else
                   IF(LMCSCF)then
                     ! Check for degeneracy irreps
                     if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                       write(1,*)icycle,"-th macro-loop, ini-grad",i,j,Y(ij)
                       nij=nij+1
                       if(dabs(G(ij)).gt.grad_digit)then
                         grad_digit=dabs(G(ij))
                       end if
                     end if
                     ! end Check for degeneracy irreps
                   else
                     ! no closed shell 
                     if(IABS(redu%orbirp(i)).lt.10.and.IABS(redu%orbirp(j)).lt.10)then
                       if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                         ! Check for degeneracy irreps
                         if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                           write(1,*)icycle,"-th macro-loop, ini-grad",i,j,Y(ij)
                           nij=nij+1
                           if(dabs(G(ij)).gt.grad_digit)then
                             grad_digit=dabs(G(ij))
                           end if
                         end if
                         ! end Check for degeneracy irreps
                       end if
                     ! i is closed shell
                     else if(IABS(redu%orbirp(i)).gt.10.and.IABS(redu%orbirp(j)).lt.10)then
                       !if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                         if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                           write(1,*)icycle,"-th macro-loop, ini-grad",i,j,Y(ij)
                           nij=nij+1
                           if(dabs(G(ij)).gt.grad_digit)then
                             grad_digit=dabs(G(ij))
                           end if
                         end if
                       !end if
                     ! j is closed shell 
                     else if(IABS(redu%orbirp(j)).gt.10.and.IABS(redu%orbirp(i)).lt.10)then
                       !if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                         if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                           write(1,*)icycle,"-th macro-loop, ini-grad",i,j,Y(ij)
                           nij=nij+1
                           if(dabs(G(ij)).gt.grad_digit)then
                             grad_digit=dabs(G(ij))
                           end if
                         end if
                       !end if
                     end if ! check closed shell
                   end if ! MCSCF or CASSCF
                 end if
              end do
            end do
          !end do j=1,i only half is needed
          
        
  
  
            nij0=nij
          ! re-run it to get the valid rotations index
            allocate(valid(nij,2)); valid=0
            ij=0
            kl=0
              !do j=i,ndim1  ! only half is needed
            do i=1,ndim1
              do j=i+1,ndim1  !half
                 ij=ij+1
                 if(dabs(Y(ij)).lt.redundant)then
                 else
                   IF(LMCSCF)then
                     if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                       kl=kl+1
                       valid(kl,1)=i
                       valid(kl,2)=j
                     end if
                   else
                     ! no closed shell 
                     if(IABS(redu%orbirp(i)).lt.10.and.IABS(redu%orbirp(j)).lt.10)then
                       if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                         if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                           kl=kl+1
                           valid(kl,1)=i
                           valid(kl,2)=j
                         end if
                       end if
                     ! i is closed shell | i-a,i-v
                     else if(IABS(redu%orbirp(i)).gt.10.and.IABS(redu%orbirp(j)).lt.10)then
                       !if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                       if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                         kl=kl+1
                         valid(kl,1)=i
                         valid(kl,2)=j
                       end if
                       !end if
                     ! j is closed shell | j-a,j-v 
                     else if(IABS(redu%orbirp(j)).gt.10.and.IABS(redu%orbirp(i)).lt.10)then
                       !if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                       if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                         kl=kl+1
                         valid(kl,1)=i
                         valid(kl,2)=j
                       end if
                       !end if
                     end if
                   end if
                 end if
              end do
            end do
            do i=1,ndim1
              do j=1,i  !half
                 ij=ij+1
                 if(dabs(Y(ij)).lt.redundant)then
                 else
                   IF(LMCSCF)then
                     if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                       kl=kl+1
                       valid(kl,1)=i
                       valid(kl,2)=j
                     end if
                   else
                     ! no closed shell 
                     if(IABS(redu%orbirp(i)).lt.10.and.IABS(redu%orbirp(j)).lt.10)then
                       if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                         if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                           kl=kl+1
                           valid(kl,1)=i
                           valid(kl,2)=j
                         end if
                       end if
                     ! i is closed shell | i-a,i-v
                     else if(IABS(redu%orbirp(i)).gt.10.and.IABS(redu%orbirp(j)).lt.10)then
                       !if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                       if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                         kl=kl+1
                         valid(kl,1)=i
                         valid(kl,2)=j
                       end if
                       !end if
                     ! j is closed shell | j-a,j-v 
                     else if(IABS(redu%orbirp(j)).gt.10.and.IABS(redu%orbirp(i)).lt.10)then
                       !if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                       if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                         kl=kl+1
                         valid(kl,1)=i
                         valid(kl,2)=j
                       end if
                       !end if
                     end if ! 
                   end if
                 end if
              end do
            end do
  
            write(1,*)"valid rotations"
            do i=1,nij
              write(1,"(I3)",advance='no')valid(i,1)
            end do
            write(1,*)" "
            do i=1,nij
              write(1,"(I3)",advance='no')valid(i,2)
            end do
          ! end  re-run it to get the valid rotations index
        !end

            nij=nij0/2

        allocate(porb(nij));   porb=0
        allocate(YP(nij));     YP=0.0d0
        allocate(HP(nij,nij)); HP=0.0d0
        allocate(HP2(nij,nij)); HP2=0.0d0
        allocate(SHP(nij,nij)); SHP=0.0
        allocate(xS(nij));     xS=1.0

        !        ij=0
        !        kl=0
        !        do i=1,ndim1
        !          do j=1,ndim1
        !             ij=ij+1
        !             if(dabs(Y(ij)).lt.redundant)then
        !             else 
        !               kl=kl+1
        !               porb(kl)=ij             
        !             end if
        !          end do
        !        end do

        !        do i=1,nij
        !          YP(i)=Y(porb(i))
        !          do j=1,nij
        !            HP(i,j)=H(porb(i),porb(j))
        !          end do
        !        end do

          ij=0
          kl=0
          do i=1,ndim1
            do j=1,ndim1 ! only half is needed
            !do j=1,ndim1 ! all
               ij=ij+1
               G(ij)=(mat2%A(i,j)-mat2%A(j,i)) ! gradients
               !Check the valid rotations  
               Lvalid=.false.
               do k=1,nij
                 !write(*,*)"k",k                 
                 if(i.eq.valid(k,1).and.j.eq.valid(k,2))then
                   Lvalid=.true.
                   exit
                 end if
               end do
               if(Lvalid)then  ! can be deleted
                 kl=kl+1
                 porb(kl)=ij
                 !write(*,*)"ij",ij,"kl",kl
               end if
            end do
          end do

          grad_digit=0.0d0
          do i=1,nij
            YP(i)= G(porb(i))
            if(dabs(YP(i)).gt.grad_digit)then  ! ??
              grad_digit=dabs(YP(i))
            end if
            do j=1,nij    ! needed Hessian
              HP(i,j)=H(porb(i),porb(j))
            end do
          end do
          YP=-YP

          allocate(x_estimate(nij));x_estimate=0.0d0
          allocate(xxx(nij));xxx=0.0d0
          allocate(xxxs(nij));xxxs=0.0d0
          allocate(XP(nij));XP=0.0d0
          allocate(xD(nij));xD=1.0d0
          itr_max=2
          mr=nij
          tol_abs = 1.0D-08
          tol_rel = 1.0D-08
          write(*,*) nij
          write(*,*)
          !$OMP PARALLEL DEFAULT(shared) PRIVATE(i,j)
          !$OMP DO
          do i=1,nij
            !$OMP PARALLEL DO
            do j=1,nij
              if(dabs(HP(j,i)).lt.redundant)then
                HP(j,i)=0.0
              end if
            end do
          end do
          !$OMP END PARALLEL
          HP2=HP
          ! write(*,*) HP

          do i=1,nij
            do j=1,nij
              sum_h=sum_h+HP(j,i)
            end do
          end do

          sum_h=sum_h/(nij*nij)
          write(*,*) sum_h
          ! write(*,*)HP(:,1)
          do i=1,nij
            do j=1,nij
              if(dabs(HP(j,i)).lt.dabs(sum_h))then
                SHP(j,i)=HP(j,i)
                HP(j,i)=0.0d0
              else
           
              end if
            end do
          end do
          ! ! write(*,*)SHP(:,1)
          
          ! write(*,*)"sum_h",sum_h

          ! call Smatrix_to_csr(SHP,nij,nij)

          ! call ax_cr (nij,SHPnonzeros,SHProwoffset,SHPcolindex,SHPvalues, xS,xxxs)
          ! write(*,*)xxxs(:100)

          ! call matrix_to_csr(HP,nij,nij)

          ! call ax_crDouble(nij,nonzeros,rowoffset,colindex,values,xD,xxx)
          ! write(*,*)xxx(:100)

          ! deallocate(rowoffset)
          ! deallocate(colindex)
          ! deallocate(values)

          ! write(*,*)SHPnonzeros
          ! write(*,*)SHProwoffset
          ! write(*,*)SHPcolindex(:170)
          ! write(*,*)SHPvalues(:170)
          ! xS=YP
          
          ! call entrywise_csr(nij,nonzeros,rowoffset,colindex,values,YP,x_estimate,SHPnonzeros,SHProwoffset,SHPcolindex,SHPvalues,xS)
          ! write(*,*) x_estimate
          ! x_estimate=0.0d0

          call matrix_to_csr(HP2,nij,nij)

          ! xxx=1.0d0
          ! call ax_crDouble(nij,nonzeros,rowoffset,colindex,values,xD,xxx)
          ! write(*,*)xxx(:100)
          call pmgmres_crd( nij, nonzeros, rowoffset, colindex, values, x_estimate, YP, itr_max, mr, tol_abs, tol_rel)

          ! call matrix_to_coo(HP,nij,nij)
          !call mgmres_st ( nij, nonzeros, coo_rows, coo_cols, coo_values, x_estimate, YP, itr_max, mr, tol_abs, tol_rel )
          nonzeros=0
          SHPnonzeros=0
          deallocate(rowoffset)
          deallocate(colindex)
          deallocate(values)

          ! deallocate(SHProwoffset)
          ! deallocate(SHPcolindex)
          ! deallocate(SHPvalues)

          ! deallocate(coo_rows)
          ! deallocate(coo_cols)
          ! deallocate(coo_values)
          XP=x_estimate
        !   write(*,*),XP(:200)

        !   if(XP(1).lt.0)then
        !     XP=-1.0d0*XP
        !   end if

          write(1,*)"AH XP"
          !call print_mat(nij,1,XP,1) 
          !          stop

          do i=1,nij
            X(porb(i))=XP(i)
          end do 
          ! write(*,*) X(:200)

  !更新mat2%R d0_sum d2_sum
          ij=0;
          d0=0.0d0; d2=0.0d0;
          d0_sum=0.0d0; d2_sum=0.0d0
          do i=1,nij
            i1=valid(i,1)
            j1=valid(i,2)
            ! recover the anti-symmetric rotations
            mat2%R(i1,j1)=  XP(i)/1.0d0
            mat2%R(j1,i1)= -XP(i)/1.0d0
            write(1,*)"rotation R in",iloop,"micro-iter",&
                      i1,j1,mat2%R(i1,j1)
            if(dabs(mat2%R(i1,j1)).gt.0.1)then
              write(*,*)"notice : large rotation between",i1,j1,&
                        " : ",mat2%R(i1,j1)
              ! reduce the redundant affect 
            !   mat2%R=mat2%R/5.0d0
            end if
            d0=d0+dabs(XP(i))
            d2=d2+(XP(i))**2
          end do
          d0_sum = d0_sum + d0
          d2_sum = d2_sum + d2
  !end
  ! recover the degenrate rotations
          is0=0
          do isym=1,orb%nsub
            js0=0
            do jsym=1,orb%nsub
              if(redu%irreps(isym).eq.redu%irreps(jsym).and.isym.lt.jsym)then
                do i=1,orb%total(isym)
                  do j=1,orb%total(isym)
                    ii=i+is0; jj=j+is0
                    kk=i+js0; ll=j+js0
            ! avoid possible over-shooting by 1/2 (?)
            mat2%R(ii,jj)=mat2%R(ii,jj)/1.0d0
            mat2%R(kk,ll)=mat2%R(ii,jj)*redu%sign_mat(kk,ll)
           write(1,*)ii,jj,mat2%R(ii,jj),kk,ll,mat2%R(kk,ll)
                  end do
                end do
              end if
              js0=js0+orb%total(jsym)
            end do
            is0=is0+orb%total(isym)
          end do
      Rabs(icycle)=d0

    !   deallocate(HPA)
      deallocate(XP)
    !   deallocate(D)
    if(.true.)then
      call CAL_T(norb,mat2%R,d3,150,mat2%T)
      call T_TO_U(norb,mat2%T,mat2%U)
    else  ! debugging, delete 
      allocate(T1(norb,norb));T1=0.0d0
      allocate(T2(norb,norb));T2=0.0d0
      call CAL_T(norb,mat2%R,d3,150,T1)
      call MXM(norb,mat2%U,T1,T2)
      mat2%T=mat2%T+T2
      call T_TO_U(norb,mat2%T,mat2%U)
      call print_mat(norb,norb,mat2%U,6)
      deallocate(T1) 
      deallocate(T2) 
    end if
  !end
    !call print_mat(norb,norb,mat2%R,6)
    write(*,2435,advance='no')"    ",Rabs(icycle)
    write(*,2436)" gmers  "
2435    format(A4,f8.4)
2436    format(A3,A7)

    deallocate(X)
    deallocate(Y,YP)
    deallocate(H,HP)

end Subroutine Solver_gmres