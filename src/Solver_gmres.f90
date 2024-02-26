Subroutine Solver_gmres(icycle)
    use global_control
    use matrix

    integer::icycle
    integer::itr_max,mr
    double precision::tol_abs,tol_rel

    integer,allocatable::porb(:)
    double precision d0,dv,redundant
    double precision,allocatable::D(:),X(:),Y(:),XP(:),YP(:)
    double precision,allocatable::H(:,:),HP(:,:),HPA(:,:)

    double precision,allocatable::T1(:,:),T2(:,:),T3(:,:)
    double precision,allocatable::x_estimate(:)

! Used for the augmented Hessian part
    double precision lamda   ! damping parameter that control the level shift

    integer,allocatable::valid(:,:)
    double precision,allocatable::G(:)
    logical Lvalid

! Parameter initialization
    lamda=1.0d0 ! "1" means the standard AH 
    ! The larger lamda, the shorter step

! initialize of part of Hessian matrix that will be used in this solver

    ndim1=norb
    

    allocate(X(ndim1**2)); X=0.0d0
    allocate(Y(ndim1**2)); Y=0.0d0
    allocate(G(ndim1**2)); G=0.0d0
    allocate(H(ndim1**2,ndim1**2)); H=0.0d0
    allocate(x_estimate(ndim1**2));x_estimate=0.0d0

    !计算梯度
    ij=0
    nij=0
    do i=1,ndim1
        do j=1,ndim1
        ij=ij+1
        Y(ij)=-1.0d0*(mat2%A(i,j)-mat2%A(j,i))
        if(dabs(Y(ij)).lt.redundant)then
        else
            nij=nij+1
        end if
        end do
    end do

    !计算Hessian 此时是4维的
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
        write(2,*)"Dh_gen finish"
        call flush(2)
        call Y_gen(orb%nsub,orb%act,orb%total,&
                    nact,norb,U,mat2%P,mat2%Y,orb%grouptable)
        write(2,*)"Y_gen finish"
        call flush(2)
        call Hessian3(norb,mat2%Dh,mat2%A,mat2%Y,mat2%H,mat2%Hdiag)
    end if
    !拉伸Hessian 4维>2维
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

    !matrix_to_csr
    ! call matrix_to_csr(H,ndim1**2,ndim1**2)
    call matrix_to_coo(H,ndim1**2,ndim1**2)
    itr_max=2
    mr=10
    tol_abs = 1.0D-08
    tol_rel = 1.0D-08
    ! call pmgmres_ilu_cr ( ndim1**2,nonzeros, rowoffset, colindex, values, x_estimate, Y, itr_max, &
    !   mr, tol_abs, tol_rel )
    call mgmres_st ( ndim1**2, nonzeros, coo_rows, coo_cols, coo_values, x_estimate, Y, itr_max, mr, &
      tol_abs, tol_rel )
    ! write(*,*)"x_estimate",x_estimate(1:ndim1**2)
    
    do i=1,ndim1
        do j=1,ndim1
            mat2%R(i,j)= x_estimate(i*ndim1+j)/1.0d0
            mat2%R(j,i)= -x_estimate(i*ndim1+j)/1.0d0
        end do
    end do
    ! mat2%R(i1,j1)=  x_estimate(i+1)/1.0d0
    ! mat2%R(j1,i1)= -x_estimate(i+1)/1.0d0
    deallocate(coo_rows)
    deallocate(coo_cols)
    deallocate(coo_values)
end Subroutine Solver_gmres

