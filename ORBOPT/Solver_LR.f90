      Subroutine Solver_LR(ndim,LHS,RHS,VEC)

        integer::ndim
        double precision::LHS(ndim,ndim) 
        double precision::RHS(ndim) 
        double precision::VEC(ndim) 

        double precision TM1(ndim,ndim),TM2(ndim,ndim)
        double precision TL1(ndim),TL2(ndim)

!        call print_mat(ndim,ndim,LHS,21)
!        call print_mat(   1,ndim,RHS,21)
!        call flush(21)
!        stop

        if(.true.)then 
          TM1=LHS
          TL1=RHS
          call inv(ndim,TM1)
          write(1,*)" inv of LHS "
          call print_mat(ndim,ndim,TM1,1)
          call MXV(ndim,ndim,TM1,TL1,VEC)
          write(1,*)" Lagrange parameters "
          call print_mat( 1,ndim,VEC,1)            
        else
          TM1=LHS
          call eigtql2(ndim,TM1,TL1) 
          write(1,*)" eig-value "
          call print_mat(   1,ndim,TL1,1)
          write(1,*)" eig-vector "
          call print_mat(ndim,ndim,TM1,1)
          call flush(1)
        end if

        stop

      End Subroutine Solver_LR

! Simple standard CG
      Subroutine Solver_CG(ndim,LHS,RHS,VEC)

        integer::ndim
        double precision::LHS(ndim,ndim)
        double precision::RHS(ndim)
        double precision::VEC(ndim)

        double precision A(ndim,ndim)
        double precision X(ndim),R(ndim),B(ndim),AP(ndim)
        double precision TL(ndim),P(ndim)
        double precision alpha,dv,rsold,rsnew

        ! initial trail vector
        X=VEC 
        ! Ax=b form
        A=LHS
        b=RHS

        call MXV(ndim,ndim,A,X,TL)
        R=B-TL
        P=R

        call inner_product(ndim,R,R,dv)
        rsold=dv

        do i=1,ndim        
          call MXV(ndim,ndim,A,P,AP)
          call inner_product(ndim,P,AP,dv)          
          alpha=rsold/dv
          X=X+alpha*P
          write(1,*)i,"th trail vector",X            
          R=R-alpha*AP
          call inner_product(ndim,R,R,dv)   
          rsnew=dv       
          if(dsqrt(rsnew).lt.1.0d-10)then
            exit
          end if
          P=R+(rsnew/rsold)*P
          rsold=rsnew
        end do

        VEC=X
        write(1,*),"th vector   X ",X 
        write(1,*),"th vector VEC ",VEC

      End Subroutine Solver_CG


! Simple standard PCG
      Subroutine Solver_PCG(ndim,LHS,RHS,PM,VEC)

        integer::ndim
        double precision::LHS(ndim,ndim)
        double precision:: PM(ndim,ndim)
        double precision::RHS(ndim)
        double precision::VEC(ndim)

        double precision A(ndim,ndim)
        double precision X(ndim),R(ndim),B(ndim),Z(ndim),AP(ndim)
        double precision R0(ndim),Z0(ndim) 
        double precision TL(ndim),P(ndim)
        double precision alpha,beta,dv,rsold,rsnew

        ! For debugging/checking
        double precision TM1(ndim,ndim),TL1(ndim)

        if(.true.)then
!          write(1,*)" The initial residual vector : "
!          call print_mat(1,ndim,RHS,1) 
!          write(1,*)" The initial trial vector : "
!          call print_mat(1,ndim,VEC,1) 
          TM1=LHS
!          write(1,*)" The o-o diagional Hessian : "
!          call print_mat1(ndim,ndim,TM1,1)
          call eigtql2(ndim,TM1,TL1)
!          write(1,*)" The  eigenvalues of SA Hessian : "
!          call print_mat(1,ndim,TL1,1)
!          write(1,*)" The eigenvectors of SA Hessian : "
!          call print_mat(ndim,ndim,TM1,1)
        end if

        ! initial trail vector
        X=VEC
        ! Ax=b form
        A=LHS
        b=RHS

        ! initial residual
        call MXV(ndim,ndim,A,X,TL)
!        write(1,*)" The TL : ", TL
        R=B-TL

        write(123,*)"initial  preconditioner"   
        call print_mat(ndim,ndim,PM,123) 
        ! inverse the preconditioner
        call inv(ndim,PM)      
        write(124,*)"inversed preconditioner"   
        call print_mat(ndim,ndim,PM,124) 

!        stop

        ! initial Z,P vector
        call MXV(ndim,ndim,PM,R,Z) 
        P=Z

        do i=1,200
          call inner_product(ndim,R,Z,dv)
          rsold=dv
          call MXV(ndim,ndim,A,P,AP)
          call inner_product(ndim,P,AP,dv)
          alpha=rsold/dv

          Z0=Z ! save Zk
          R0=R ! save Rk
          X=X+alpha*P
          R=R-alpha*AP

!          write(1,*)"X-vectors i-th ",i
!          call print_mat(1,ndim,X,1)

!          write(1,*)" residual i-th ",i
!          write(1,*)i,"th residual",R
          !call print_mat(1,ndim,R,1)

          call inner_product(ndim,R,R,dv)
          rsnew=dv
          if(dsqrt(rsnew).lt.1.0d-10)then
            exit
          end if

          call MXV(ndim,ndim,PM,R,Z) 
          call inner_product(ndim,Z,R,dv)
          rsnew=dv
          call inner_product(ndim,Z0,R0,dv)
          beta=rsnew/dv 

          P=Z+beta*P
!          rsold=rsnew
        end do

        VEC=X
        write(1,*)"X vector, only orbital case",vec
        
      End Subroutine Solver_PCG


! PCG with orthogonal to SA space
      Subroutine Solver_PCG_SA(ndim,LHS,RHS,PM,VEC,SAvec)

        integer::ndim
        double precision::LHS(ndim,ndim)
        double precision:: PM(ndim,ndim)
        double precision::  RHS(ndim)
        double precision::  VEC(ndim)
        double precision::SAvec(ndim)

        double precision A(ndim,ndim)
        double precision X(ndim),R(ndim),B(ndim),Z(ndim),AP(ndim)
        double precision R0(ndim),Z0(ndim) 
        double precision TL(ndim),P(ndim)
        double precision alpha,beta,dv,rsold,rsnew

        double precision VS(ndim,2)

        ! For debugging/checking
        double precision TM1(ndim,ndim),TL1(ndim)

        if(.true.)then
          TM1=LHS
          call eigtql2(ndim,TM1,TL1)
!          write(1,*)" The  eigenvalues of SA Hessian : "
!          call print_mat(1,ndim,TL1,1)
!          write(1,*)" The eigenvectors of SA Hessian : "
!          call print_mat(ndim,ndim,TM1,1)
        end if

        ! state-averaged space
        VS(:,1)=SAvec

        write(1,*)"The SA vec-space" 
        call print_mat(1,ndim,VS(:,1),1)

        ! initial trail vector
        X=VEC
        ! Ax=b form
        A=LHS
        b=RHS

        ! initial residual
        call MXV(ndim,ndim,A,X,TL)
        R=B-TL

        ! inverse the preconditioner
        call inv(ndim,PM)       
        ! initial Z,P vector
        call MXV(ndim,ndim,PM,R,Z) 
        P=Z

        do i=1,200
          call inner_product(ndim,R,Z,dv)
          rsold=dv
          call MXV(ndim,ndim,A,P,AP)
          call inner_product(ndim,P,AP,dv)
          alpha=rsold/dv

          Z0=Z ! save Zk
          R0=R ! save Rk
          X=X+alpha*P
          R=R-alpha*AP

!          write(1,*)i,"-th X-vectors"
!          call print_mat(1,ndim,X,1)
!          write(1,*)i,"-th  residual"
!          call print_mat(1,ndim,R,1)

          call inner_product(ndim,R,R,dv)
          rsnew=dv
          if(dsqrt(rsnew).lt.1.0d-10)then
            write(1,*)"========================================="
            write(1,*)"     Convergence for PCG is reached"
            write(1,*)"========================================="
            exit
          end if

          call MXV(ndim,ndim,PM,R,Z) 
          call inner_product(ndim,Z,R,dv)
          rsnew=dv
          call inner_product(ndim,Z0,R0,dv)
          beta=rsnew/dv 

          P=Z+beta*P

          ! direction vector is chosen orthogonal to SA space
          VS(:,2)=P
          call MODIFIED_ORTHOGONALIZE(ndim,2,VS(:,1:2))
          P=VS(:,2)          

          call inner_product(ndim,P,SAvec,dv)      
          write(1,*)"check for orthogonal",dv 

!          rsold=rsnew
        end do

        VEC=X
        write(1,*)"X vector",vec

      End Subroutine Solver_PCG_SA

! PCG only SA space



! CG with projectting out the SA space
      Subroutine Solver_CG_SA(nij,nC0,LHS,RHS,VEC,nS0,conf,SAvec)

        integer::nij,nC0,nS0
        integer::conf(nS0)
        double precision::LHS(nij+nC0,nij+nC0)
        double precision::RHS(nij+nC0)
        double precision::VEC(nij+nC0)
        double precision::SAvec(nC0)

        double precision A(nij+nC0,nij+nC0)
        double precision X(nij+nC0), R(nij+nC0)
        double precision B(nij+nC0),AP(nij+nC0)
        double precision TL(nij+nC0),P(nij+nC0)
        double precision alpha,dv,rsold,rsnew

        ! For debugging/checking
        double precision TM1(nij+nC0,nij+nC0),TL1(nij+nC0)

        ! For projection
        double precision,allocatable::PAvec(:),SSvec(:)
        double precision,allocatable::VS(:,:)

        write(1,*)"Entering Solver_CG_SA solver"
        call flush(1) 
         
        ndim=nij+nC0

        if(.true.)then
          TM1=LHS
          write(1,*)" The  SA Hessian : "
          call print_mat(ndim,ndim,TM1,1)
          call eigtql2(ndim,TM1,TL1)
!          write(1,*)" The  eigenvalues of SA Hessian : "
!          call print_mat(1,ndim,TL1,1)
!          call flush(1)
!          write(1,*)" The eigenvectors of SA Hessian : "
!          call print_mat(ndim,ndim,TM1,1)
!          call flush(1)
        end if

!        VS=0.0d0

        ! initial trail vector
        X=VEC
        ! Ax=b form
        A=LHS
        b=RHS

        call MXV(ndim,ndim,A,X,TL)
        R=B-TL
        P=R

        write(1,*)"The state-averaged vectors"
        write(1,*) SAvec 
        write(1,*)"The  x       "
        call print_mat(1,ndim, x,1)
        write(1,*)"The  b       "
        call print_mat(1,ndim, b,1)
        write(1,*)"The Ax       "
        call print_mat(1,ndim,TL,1)
        write(1,*)"The  R (b-Ax)"
        call print_mat(1,ndim, R,1)

        !VS(:,1)=SAvec
!        write(1,*)"th residual vector"
!        call print_mat(1,ndim,R,1) 
!        write(1,*)"                    "
        call inner_product(ndim,R,R,dv)
        rsold=dv

        do i=1,50

          i0=0 ! offset for MPS vectors
          do ist=1,nS0 
            nC=conf(ist)
            allocate(PAvec(nC)); PAvec=0.0d0
            allocate(SSvec(nC)); SSvec=0.0d0
            allocate(VS(nC,2));     VS=0.0d0

            PAvec=P(nij+i0+1:nij+i0+nC)
            SSvec=SAvec(i0+1:i0+nC)
             
            call inner_product(nC,PAvec,SSvec,dv)
!            write(1,*)"check initial ORTHOGONALIZE", dv   

            VS(:,1)=SSvec
            VS(:,2)=PAvec
            write(1,*)"VS matrix ",ist 
            call print_mat(nC,2,VS,1)

!           call MODIFIED_GRAM_SCHMIDT(ndim,i+1,VS(:,1:i+1))
            call MODIFIED_ORTHOGONALIZE(nC,2,VS(:,1:2))

            if(nC.gt.1)then
              PAvec=VS(:,2)
            else
              PAvec=0.0d0
            end if 

            write(1,*)"PAvec for state ",ist 
            call print_mat(1,nC,PAvec,1)
         
            P(nij+i0+1:nij+i0+nC)=PAvec

            call inner_product(nC,PAvec,SSvec,dv)
            write(1,*)"check after ORTHOGONALIZE", dv   

            deallocate(PAvec)
            deallocate(SSvec)
            deallocate(VS)

            i0=i0+conf(ist)
          end do

          call MXV(ndim,ndim,A,P,AP)
          call inner_product(ndim,P,AP,dv)
          alpha=rsold/dv
          write(1,*)" ============================== "  
          write(1,*) "rsold,dv,alpha",rsold,dv,alpha
          call print_mat1(1,ndim, P,1)
          write(1,*)" ============ AP ============ "  
          call print_mat1(1,ndim,AP,1)
          write(1,*)" ============  X ============ "  
          call print_mat1(1,ndim,X,1)
          write(1,*)" ============  R ============ "  
          call print_mat1(1,ndim,R,1)
          write(1,*)" ============================== "  

          X=X+alpha*P
          R=R-alpha*AP

          call inner_product(ndim,R,R,dv)
          rsnew=dv
          if(dsqrt(rsnew).lt.1.0d-10)then
            exit
          end if

!!          do i=1,nS0  
!!            write(1,*)"SA vector space" 
!!            write(1,*) VS(:,1)

          P=R+(rsnew/rsold)*P

          write(1,*)" ============ new X ============ "  
          call print_mat1(1,ndim,X,1)
          write(1,*)" ============ new R ============ "  
          call print_mat1(1,ndim,R,1)
          write(1,*)" ============ new P ============ "  
          call print_mat1(1,ndim,P,1)
          write(1,*) "rsnew",rsnew
          write(1,*)" =============================== "  

          rsold=rsnew
!          write(1,*)i,"th    trail vector"
!          call print_mat(1,ndim,P,1) 
!          write(1,*)i,"th vector"
!          call print_mat(1,ndim,X,1) 
!          write(1,*)i,"th residual vector"
!          call print_mat(1,ndim,R,1) 
!          write(1,*)"                    "

        end do

        VEC=X
        write(1,*)i,"th    final vector"
        call print_mat(1,ndim,VEC,1) 

      End Subroutine Solver_CG_SA

