! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      SUBROUTINE transform_subspaces(nsub,act,tot,norb,TRANS,T,U,group)
        ! Yingjin in around 2013
        ! use omp_lib
        use date_time
        use omp_lib          

        INTEGER::norb
        DOUBLE PRECISION::TRANS(norb,norb)
        DOUBLE PRECISION::T(norb,norb)
        DOUBLE PRECISION::U(norb,norb,norb,norb)

        integer::act(nsub) 
        integer::tot(nsub) 
        integer::group(8,8) 

        type::transform
          double precision,allocatable::U(:,:) 
        end type transform
        type(transform),allocatable::UM(:) 

        double precision,allocatable::TM1(:,:),TM2(:,:),TM3(:,:),TM(:,:)
        double precision,allocatable::TG1(:,:,:,:),TG2(:,:,:,:)
 
        allocate(UM(nsub))
        i0=0
        do i=1,nsub
          if(tot(i).ne.0)then ! avoid the void array 
            allocate(UM(i)%U(tot(i),tot(i)))
            UM(i)%U=0.0d0
            UM(i)%U=TRANS(i0+1:i0+tot(i),i0+1:i0+tot(i))
            !            call print_mat(tot(i),tot(i),UM(i)%U)
            i0=i0+tot(i)
          end if
        end do
       
        allocate(TM(norb,norb));TM=0.0d0 
        i0=0
        do i=1,nsub 
          if(tot(i).ne.0)then ! avoid the void array
            allocate(TM1(tot(i),tot(i)));TM1=0.0d0
            allocate(TM2(tot(i),tot(i)));TM2=0.0d0
            TM1=T(i0+1:i0+tot(i),i0+1:i0+tot(i))
            call MdXM(tot(i),UM(i)%U,TM1,TM2) 
            call  MXM(tot(i),TM2,UM(i)%U,TM1) 
            TM(i0+1:i0+tot(i),i0+1:i0+tot(i))=TM1
            deallocate(TM1)
            deallocate(TM2)
            i0=i0+tot(i)
          end if 
        end do     
        ! The transformed 1-body integrals
        T=TM
        ! Two Electron Integrals
        ! Transform
        ALLOCATE(TG1(norb,norb,norb,norb));TG1=0.0d0
        ALLOCATE(TG2(norb,norb,norb,norb));TG2=0.0d0
        walltime(80) = wtime()
        i0=0
        DO i=1,nsub
          j0=0
          DO j=1,nsub
            k0=0   
            DO k=1,nsub
              l0=0
              DO l=1,nsub
                if(group(i,group(j,group(k,l))).eq.1)then
                  allocate(TM1(tot(k),tot(l)))
                  allocate(TM2(tot(k),tot(l)))
                  allocate(TM3(tot(k),tot(l)))
                  if(tot(k).ne.0.and.tot(l).ne.0)then
                    !$OMP PARALLEL DEFAULT(shared) PRIVATE(ii, jj, TM1, TM2,TM3,dtmp)
                    !$OMP DO
                    do ii=1,tot(i)
                      do jj=1,tot(j) 
                          TM1=0.0d0
                          TM2=0.0d0
                          TM3=0.0d0
                          TM1=U(I0+ii,J0+jj,k0+1:k0+tot(k),l0+1:l0+tot(l))
                          call MXMG(tot(k),tot(l),tot(k),UM(k)%U,TM1,TM2,"TN")
                          call MXMG(tot(k),tot(l),tot(l),TM2,UM(l)%U,TM3,"NN")
                          TG1(I0+ii,J0+jj,k0+1:k0+tot(k),l0+1:l0+tot(l))=TM3
                      end do
                    end do
                    !$OMP END PARALLEL
                  end if
                  deallocate(TM1)
                  deallocate(TM2)
                  deallocate(TM3)
                end if
                l0=l0+tot(l)
              END DO
              k0=k0+tot(k)
            END DO
            j0=j0+tot(j)
          END DO
          i0=i0+tot(i)
        END DO
        walltime(81) = wtime()

        !        call print_gat(norb,norb,norb,norb,MIDU1,11)
        !        stop

        i0=0
        DO i=1,nsub
          j0=0
          DO j=1,nsub
            k0=0
            if(tot(i).ne.0.and.tot(j).ne.0)then
              allocate(TM1(tot(i),tot(j)))
              allocate(TM2(tot(i),tot(j)))
              allocate(TM3(tot(i),tot(j)))
              DO k=1,nsub
                l0=0
                DO l=1,nsub
                  if(group(i,group(j,group(k,l))).eq.1)then
                    !$OMP PARALLEL DEFAULT(shared) PRIVATE(ii, jj, TM1, TM2,TM3,dtmp)
                    !$OMP DO
                    do kk=1,tot(k)
                      do ll=1,tot(l) 
                          TM1=0.0d0
                          TM2=0.0d0
                          TM3=0.0d0
                          TM1=TG1(i0+1:i0+tot(i),j0+1:j0+tot(j),k0+kk,l0+ll)
                          call MXMG(tot(i),tot(j),tot(i),UM(i)%U,TM1,TM2,"TN")
                          call MXMG(tot(i),tot(j),tot(j),TM2,UM(j)%U,TM3,"NN")
                          TG2(i0+1:i0+tot(i),j0+1:j0+tot(j),k0+kk,l0+ll)=TM3
                      end do
                    end do
                    !$OMP END PARALLEL
                  end if
                  l0=l0+tot(l)
                END DO
                k0=k0+tot(k)
              END DO
              deallocate(TM1)
              deallocate(TM2)
              deallocate(TM3)
            end if
            j0=j0+tot(j)
          END DO
          i0=i0+tot(i)
        END DO

        walltime(82) = wtime()
        write(*,*)"*******************T_1*****",walltime(81)-walltime(80)
        write(*,*)"*******************T_2*****",walltime(82)-walltime(81)


        !        call print_gat(norb,norb,norb,norb,TG2,15)
        !        stop  !!!!=======================!!!!
        ! The transformerd 2-body matrix
        U=TG2 
        deallocate(TM)
        deallocate(UM)
        deallocate(TG1)
        deallocate(TG2)

      END SUBROUTINE transform_subspaces

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! Print the matrix elegantly
      Subroutine PRINT_MAT(m,n,M1,ifort)

        implicit none 

!        character*10 name_mat
        integer::m,n,ifort
        integer i,j
        double precision::M1(m,n)

!         write(ifort,*)"== Print Mat =="
!        if(n.lt.9)then
          do i=1,m
            do j=1,n
              write(ifort,"(f7.4,1X)",advance='no')M1(i,j)
              !write(ifort,"(f9.6,1X)",advance='no')M1(i,j)
              !write(ifort,"(f12.8,1X)",advance='no')M1(i,j)
            end do
            write(ifort,*)
          end do    
          call flush(ifort) 
!        else

      end Subroutine PRINT_MAT
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      Subroutine PRINT_MAT1(m,n,M1,ifort)

        implicit none

!        character*10 name_mat
        integer::m,n,ifort
        integer i,j
        double precision::M1(m,n)

!         write(ifort,*)"== Print Mat =="
!        if(n.lt.9)then
          do i=1,m
            do j=1,n
              !write(ifort,"(f7.3,1X)",advance='no')M1(i,j)
              !write(ifort,"(f4.1,1X)",advance='no')M1(i,j)
              write(ifort,"(f12.8,1X)",advance='no')M1(i,j)
            end do
            write(ifort,*)
          end do
!        else

      end Subroutine PRINT_MAT1


! Print the matrix elegantly
      Subroutine PRINT_MAT2(m,n,M1,ifort,if1,if2)

        implicit none

!        character*10 name_mat
        integer::m,n,ifort
        integer i,j,if1,if2
        character*20 ctmp
        double precision::M1(m,n)

        ctmp=""
         

!        write(ifort,*)"== Print Mat =="
!        if(n.lt.9)then
          do i=1,m
            do j=1,n
              !write(ifort,"(f7.3,1X)",advance='no')M1(i,j)
              write(ifort,trim(ctmp),advance='no')M1(i,j)
            end do
            write(ifort,*)
          end do
!        else

      end Subroutine PRINT_MAT2


! -------------------------------------------------------
! Print the 3-dims matrix elegantly
      Subroutine PRINT_TAT(ic,m,n,M1,ifort)

        implicit none

!        character*10 name_mat
        integer::ic,m,n,ifort
        integer i0,i,j
        double precision::M1(ic,m,n)

!        if(n.lt.9)then
        do i0=1,ic
          write(ifort,"(A3,I3,A8)")"DIM ",i0," is shown"
          do i=1,m
            do j=1,n
              write(ifort,"(f7.3,1X)",advance='no')M1(i0,i,j)
            end do
            write(ifort,*)
          end do
          write(ifort,*)"======"
        end do
!        else

      end Subroutine PRINT_TAT

! --------------------------------------------------------
! Print the 4-tensor elegantly
      Subroutine PRINT_GAT(m,n,p,q,M1,ifort)

        implicit none 

!        character*10 name_mat
        integer::m,n,p,q,ifort
        integer i,j,k,l
        double precision::M1(m,n,p,q)

!        if(n.lt.9)then
          do i=1,m
            do j=1,n
              write(ifort,"(A3,I3,I3,A8)")"DIM ",i,j," is shown"
              do k=1,p    
                do l=1,q
                  !write(ifort,"(f7.3,1X)",advance='no')M1(i,k,l,j)
                  ! write(ifort,"(f12.8,1X)",advance='no')M1(i,k,l,j) ! origin 
                  write(ifort,"(f12.8,1X)",advance='no')M1(i,j,k,l) 
                  !write(ifort,"(f7.4,1X)",advance='no')M1(i,j,k,l)
                  !write(ifort,"(f12.8,1X)",advance='no')M1(k,i,l,j)
                end do
                write(ifort,*)
              end do
              write(ifort,*)"======"
            end do
          end do  
          call flush(ifort)
     
!        else

      end Subroutine PRINT_GAT
! --------------------------------------------------------
! Trace if a matrix
      Subroutine trace(n,M1,dv)

        implicit none

        integer::n
        double precision::dv
        double precision::M1(n,n)

        integer i
       
        dv=0.0d0 
        do i=1,n
          dv=dv+M1(i,i) 
        end do 
 
!        write(*,*)dv

      End subroutine trace 
! ---------------------------------------------------------
      Subroutine traceG(n1,n2,M1,dv)

        implicit none

        integer::n1,n2
        double precision::dv
        double precision::M1(n1,n1)

        integer i,n
       
        dv=0.0d0 
        if(n1.gt.n2)then
          n=n2
        else
          n=n1
        end if 
 
        do i=1,n
          dv=dv+M1(i,i) 
        end do 
 
      End subroutine traceG 
! ---------------------------------------------------------
!  Matrix multi Matrix
      Subroutine MXM(n,M1,M2,M3)

        implicit none

        integer::n
        double precision::M1(n,n)
        double precision::M2(n,n)
        double precision::M3(n,n)

        M3=0.0d0

        call dgemm("N","N",n,n,n,1.0d0,M1,n,M2,n,0.0d0,M3,n)

      end subroutine
! ---------------------------------------------------------
!  Matrix multi Vector
      Subroutine MXV(m,n,M1,V1,V2)

        integer::m,n
        double precision::M1(m,n)
        double precision::V1(n)
        double precision V2(m)
       
        V2=0.0d0 
        call dgemv('N',m,n,1.0d0,M1,m,V1,1,0.0d0,V2,1)

      end subroutine 

! ---------------------------------------------------------
!  Matrix'dagger multi Matrix
      Subroutine MdXM(n,M1,M2,M3)

        implicit none

        integer::n
        double precision::M1(n,n)
        double precision::M2(n,n)
        double precision::M3(n,n)

        M3=0.0d0

        call dgemm("T","N",n,n,n,1.0d0,M1,n,M2,n,0.0d0,M3,n)

      end subroutine
! ---------------------------------------------------------

! ---------------------------------------------------------
!  Matrix multi Matrix'dagger
      Subroutine MXMd(n,M1,M2,M3)

        implicit none

        integer::n
        double precision::M1(n,n)
        double precision::M2(n,n)
        double precision::M3(n,n)

        M3=0.0d0

        call dgemm("N","T",n,n,n,1.0d0,M1,n,M2,n,0.0d0,M3,n)

      end subroutine
! ---------------------------------------------------------
!  General MXM
      Subroutine MXMG(m,n,k,M1,M2,M3,NT)

! M1, M2, M3 should be different 

! M1: m*k
! M2: k*n
! M3: m*n
 
        implicit none
        integer::m,n,k
        character*2::NT 
        double precision::M1(m,k)
        double precision::M2(k,n)
        double precision  M3(m,n)

        M3=0.0d0

!        do i=1,n
!          write(*,*)M2(i,:)
!        end do
       
         if(NT.eq."NN")then
          call dgemm("N","N",m,n,k,1.0d0,M1,m,M2,k,0.0d0,M3,m)
         else if(NT.eq."NT")then
!          call dgemm("N","T",m,n,k,1.0d0,M1,m,M2,k,0.0d0,M3,m)
          call dgemm("N","T",m,n,k,1.0d0,M1,m,M2,n,0.0d0,M3,m)
         else if(NT.eq."TN")then
          call dgemm("T","N",m,n,k,1.0d0,M1,k,M2,k,0.0d0,M3,m)
         else if(NT.eq."TT")then
          call dgemm("T","T",m,n,k,1.0d0,M1,k,M2,n,0.0d0,M3,m)
         end if 

!        do i=1,n
!          write(*,*)M3(i,:)
!        end do

      end Subroutine MXMG
! ---------------------------------------------------------
!   Ul*M*UR'
      Subroutine trans2(n,ul,T,ur)

        integer::n
        double precision::ul(n,n),ur(n,n)
        double precision::T(n,n)
        double precision Ts(n,n)

        Ts=0.0d0

        call dgemm("N","N",n,n,n,1.0d0,ul,n,T,n,0.0d0,Ts,n)
        call dgemm("N","T",n,n,n,1.0d0,Ts,n,ur,n,0.0d0,T,n)

      end Subroutine trans2

!============================================================
!  U1*M*UR
      Subroutine trans(n,ul,T,ur)

        integer::n
        double precision::ul(n,n),ur(n,n)
        double precision::T(n,n)
        double precision Ts(n,n)

        Ts=0.0d0

        call dgemm("N","N",n,n,n,1.0d0,ul,n,T,n,0.0d0,Ts,n)
        call dgemm("N","N",n,n,n,1.0d0,Ts,n,ur,n,0.0d0,T,n)

      end Subroutine trans

!================================================================
!  diagonalization
      SUBROUTINE EIGTQL2(DIM,EIGVEC,EIGVAL)

        INTEGER::DIM
        DOUBLE PRECISION::EIGVEC(DIM,DIM)
        DOUBLE PRECISION::EIGVAL(DIM)

        character*1 uplo
        integer n,lda,lwork
        double precision,allocatable::a(:,:),tmatrix(:,:)
        double precision,allocatable::d(:),e(:),tau(:)
        double precision,allocatable::work(:),work2(:)

        character*1 compz
        integer ldz
        double precision,allocatable::z(:,:)

        allocate(a(dim,dim),tmatrix(dim,dim))
        allocate(d(dim),e(dim),tau(dim))
        allocate(work(dim))

        !c       write(*,*)eigvec

        uplo ='U'
        n = DIM
        a = EIGVEC
        lda = DIM
        lwork = DIM

        d=0.0d0;e=0.0d0;tau=0.0d0

        CALL DSYTRD(uplo,n,a,lda,d,e,tau,work,lwork,info)

        lda = dim; lwork = dim
        CALL DORGTR(uplo,n,a,lda,tau,work,lwork,info)

        allocate(z(dim,dim))
        allocate(work2(2*n-2))
        compz = 'V'
        ldz = Dim
        z = a

        CALL DSTEQR(compz,n,d,e,z,ldz,work2,info)

        EIGVEC = z
        EIGVAL = d

        deallocate(a,tmatrix,d,e,z,tau,work,work2)

      END SUBROUTINE EIGTQL2

!==========================================================
!     matrxi inversion
      Subroutine inv(n,TM)

        integer::n
        double precision::TM(n,n)
        integer info
        !integer ipiv(n)
        integer,allocatable::ipiv(:)
        
        double precision,allocatable::T2(:,:)
        double precision,allocatable::OM(:)

        write(6,*)"entring the inv ",n
        call flush(6)

        allocate(T2(n,n));T2=0.0d0 
        allocate(OM(n*n));OM=0.0d0
        allocate(ipiv(n));ipiv=0
        T2=TM
        OM=0.0d0
        call dgetrf(n,n,TM,n,ipiv,info)
        call dgetri(n,TM,n,ipiv,OM,n*n,info)
        deallocate(T2)
        deallocate(OM)
        deallocate(ipiv)

      end Subroutine inv

!=========================================================
!   T -> U in MCSCF 
! to Stefan : Eq.12
      SUBROUTINE T_TO_U(iORB,TMATRIX,UMATRIX)

        double precision UMATRIX(iorb,iorb)
        double precision TMATRIX(iorb,iorb)
        INTEGER::iORB

        UMATRIX=0.0d0

        DO i=1,iORB
           DO j=1,iORB
             UMATRIX(i,j)=TMATRIX(i,j)
           end do
             UMATRIX(i,i)=UMATRIX(i,i)+1.0d0
        end do

      end SUBROUTINE T_TO_U

! =========================================================
!  R-> T in MCSCF
! to Stefan : Eq.13
      SUBROUTINE CAL_T(iORB,RM,DVA,INTP,TMATRIX)

        INTEGER::iORB,INTP
        DOUBLE PRECISION::DVA
        DOUBLE PRECISION::RM(iORB,iORB)
        double precision::TMATRIX(iorb,iorb)

        DOUBLE PRECISION,ALLOCATABLE::TMAT(:,:)
        DOUBLE PRECISION,ALLOCATABLE::TIMG(:,:)
        DOUBLE PRECISION,ALLOCATABLE::RT1(:,:)
        DOUBLE PRECISION,ALLOCATABLE::RT2(:,:)
        DOUBLE PRECISION,ALLOCATABLE::MTEMP1(:,:)
        DOUBLE PRECISION,ALLOCATABLE::MTEMP2(:,:)
        DOUBLE PRECISION PRODN


        ALLOCATE(TMAT(iorb,iorb))
        ALLOCATE(TIMG(iorb,iorb))
        ALLOCATE(RT1(iorb,iorb))
        ALLOCATE(RT2(iorb,iorb))
        ALLOCATE(MTEMP1(iorb,iorb))
        ALLOCATE(MTEMP2(iorb,iorb))

        ipod=1;
        RT1=0.0d0
        RT2=0.0d0

        if(DVA.GT.10.0d0)then
          ipod=DVA/10.0d0
        end if

        !write(*,*)ipod
!        stop
!        write(711,*)"================================="
        do i=1,iorb
          do j=1,iorb
            RT1(i,j)=RM(i,j)/ipod
!            write(711,*)i,j,RT1(i,j)
          end do
        end do

        MTEMP1=RT1
        MTEMP2=0.0d0
        TMAT=0.0d0
        TIMG=0.0d0
        TMATRIX=0.0d0

        DO i=1,INTP
          ii=i
          PRODN=i*1.0d0
          IF(i.GT.1)THEN
            CALL DGEMM('N','N',iorb,iorb,iorb,1.0d0,MTEMP1,iorb,RT1,&
                       iorb,0.0d0,MTEMP2,iorb)
          ELSE
            MTEMP2=MTEMP1
          END IF
          !         WRITE(*,*)MTEMP2
          !         STOP
          IF(ii.GT.1)THEN
            DO
              IF(ii.EQ.1)EXIT
              PRODN=PRODN*(ii-1)
              ii=ii-1
            END DO
          END IF
          !         WRITE(*,*)PRODN
          !         STOP
          DO i2=1,iORB
            DO j2=1,iORB
              TMAT(i2,j2)=TMAT(i2,j2)+1.0d0/PRODN*MTEMP2(i2,j2)
          !             WRITE(*,*)i2,j2,TMAT(i2,j2)
              MTEMP1(i2,j2)=MTEMP2(i2,j2)
            END DO
          END DO
        END DO

        do i2=1,iorb
           TMAT(i2,i2)=TMAT(i2,i2)+1
        end do

        do i=1,iorb
          do j=1,iorb
            TIMG(i,j)=TMAT(i,j)
          end do
        end do
        if(ipod.GT.1)then
          do i=1,ipod-1
            CALL DGEMM('N','N',iORB,iORB,iORB,1.0d0,TMAT,iORB,TIMG,&
                      iORB,0.0d0,RT2,iORB)
            TMAT=RT2
          end do
        end if

!        OPEN(UNIT=2009,FILE='TMAT.txt')
          DO i=1,iORB
            TMAT(i,i)=TMAT(i,i)-1
            DO j=1,iORB
!              WRITE(2009,1010)i,j,TMAT(i,j)
              TMATRIX(i,j)=TMAT(i,j)
            END DO
          END DO
!1010      FORMAT(1X,'TMAT(',I2,',',I2,')=',D24.16,';')

        DEALLOCATE(TMAT,TIMG,RT1,RT2,MTEMP1,MTEMP2)

      END SUBROUTINE CAL_T

! =============================================================
!  Integrals transform
      SUBROUTINE transform(iORB,TRANSFMATRIX,T,U)

        INTEGER::iORB
        DOUBLE PRECISION::TRANSFMATRIX(iorb,iorb)
        DOUBLE PRECISION::T(iORB,iORB)
        DOUBLE PRECISION::U(iORB,iORB,iORB,iORB)
        DOUBLE PRECISION,allocatable::NEWT(:,:),NEWU(:,:,:,:)
        DOUBLE PRECISION,ALLOCATABLE::MTR1(:,:),MTR2(:,:),MTR3(:,:)
        DOUBLE PRECISION,ALLOCATABLE::MIDU1(:,:,:,:)
        
        double precision,allocatable::tempu(:,:,:,:)
        double precision,allocatable::tempnewu(:,:,:,:)

        allocate(NEWT(iorb,iorb))
        allocate(NEWU(iorb,iorb,iorb,iorb))

        NEWT=0.0d0;NEWU=0.0d0

! One Electron Integrals
! Transform and output the One electron Integrals
        ALLOCATE(MTR1(iorb,iorb))

        CALL DGEMM('T','N',iORB,iORB,iORB,1.0D0,TRANSFMATRIX,&
                   iORB,T,iORB,0.0D0,MTR1,iORB)

        CALL DGEMM('N','N',iORB,iORB,iORB,1.0D0,MTR1,iORB,&
                    TRANSFMATRIX,iORB,0.0D0,NEWT,iORB)
        DEALLOCATE(MTR1)
        

! Two Electron Integrals

        allocate(tempu(iorb,iorb,iorb,iorb))
        do i=1,iorb; do j=1,iorb; do k=1,iorb; do l=1,iorb
        tempu(i,k,l,j)=u(i,j,k,l)
        end do; end do; end do; end do

        allocate(tempnewu(iorb,iorb,iorb,iorb))
! Transform
        ALLOCATE(MIDU1(iORB,iORB,iORB,iORB))
        ALLOCATE(MTR1(iORB,iORB))
        ALLOCATE(MTR2(iORB,iORB))
        ALLOCATE(MTR3(iORB,iORB))
        DO L=1,iORB
          DO K=1,iORB

                DO J=1,iORB; DO I=1,iORB
                MTR1(I,J)=tempU(I,J,K,L)
                END DO; END DO

        CALL DGEMM('T','N',iORB,iORB,iORB,1.0D0,TRANSFMATRIX,&
                   iORB,MTR1,iORB,0.0D0,MTR2,iORB)

        CALL DGEMM('T','T',iORB,iORB,iORB,1.0D0,TRANSFMATRIX,&
                   iORB,MTR2,iORB,0.0D0,MTR3,iORB)

                DO JJ=1,iORB; DO II=1,iORB
                MIDU1(II,JJ,K,L)=MTR3(II,JJ)
                END DO; END DO

          END DO
        END DO

        DO II=1,iORB
          DO JJ=1,iORB
                DO L=1,iORB; DO K=1,iORB
                MTR1(K,L)=MIDU1(II,JJ,K,L)
                END DO; END DO

        CALL DGEMM('N','N',iORB,iORB,iORB,1.0D0,MTR1,&
                   iORB,TRANSFMATRIX,iORB,0.0D0,MTR2,iORB)

        CALL DGEMM('T','N',iORB,iORB,iORB,1.0D0,MTR2,iORB,&
                    TRANSFMATRIX,iORB,0.0D0,MTR3,iORB)

                DO LL=1,iORB; DO KK=1,iORB
                tempNEWU(II,JJ,KK,LL)=MTR3(KK,LL)
                END DO; END DO

          END DO
        END DO

        do i=1,iorb; do j=1,iorb; do k=1,iorb; do l=1,iorb
        newu(i,j,k,l)=tempnewu(i,k,l,j)
        end do; end do; end do; end do
        deallocate(tempu,tempnewu)
        DEALLOCATE(MTR1,MTR2,MTR3)
        DEALLOCATE(MIDU1)

        OPEN(UNIT=105,FILE='SEINTEGRALS_NEW',FORM='UNFORMATTED',&
             ACCESS='SEQUENTIAL')
        DO I=1,iORB; DO J=1,iORB
          WRITE(105)NEWT(i,j)
        end do;end do
        close(105)

        OPEN(UNIT=106,FILE='DEINTEGRALS_NEW',FORM='UNFORMATTED',&
             ACCESS='SEQUENTIAL')
        DO I=1,iORB; DO J=1,iORB; DO K=1,iORB; DO L=1,iORB
          WRITE(106)NEWU(I,J,K,L)
        END DO; END DO; END DO; END DO
        CLOSE(106)

        deallocate(newt,newu)

      END SUBROUTINE transform

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      SUBROUTINE transform_RAM(iORB,TRANSFMATRIX,T,U)

        INTEGER::iORB
        DOUBLE PRECISION::TRANSFMATRIX(iorb,iorb)
        DOUBLE PRECISION::T(iORB,iORB)
        DOUBLE PRECISION::U(iORB,iORB,iORB,iORB)
        DOUBLE PRECISION,allocatable::NEWT(:,:),NEWU(:,:,:,:)
        DOUBLE PRECISION,ALLOCATABLE::MTR1(:,:),MTR2(:,:),MTR3(:,:)
        DOUBLE PRECISION,ALLOCATABLE::MIDU1(:,:,:,:)
        
        double precision,allocatable::tempu(:,:,:,:)
        double precision,allocatable::tempnewu(:,:,:,:)

        allocate(NEWT(iorb,iorb))
        allocate(NEWU(iorb,iorb,iorb,iorb))

        NEWT=0.0d0;NEWU=0.0d0

! One Electron Integrals
! Transform and output the One electron Integrals

        ALLOCATE(MTR1(iorb,iorb))
        CALL DGEMM('T','N',iORB,iORB,iORB,1.0D0,TRANSFMATRIX,&
                   iORB,T,iORB,0.0D0,MTR1,iORB)

        CALL DGEMM('N','N',iORB,iORB,iORB,1.0D0,MTR1,iORB,&
                   TRANSFMATRIX,iORB,0.0D0,NEWT,iORB)
        DEALLOCATE(MTR1)
        
! Two Electron Integrals

        allocate(tempu(iorb,iorb,iorb,iorb))
        do i=1,iorb; do j=1,iorb; do k=1,iorb; do l=1,iorb
        tempu(i,k,l,j)=u(i,j,k,l)
        end do; end do; end do; end do
!        call print_gat(iorb,iorb,iorb,iorb,tempU,8)

        allocate(tempnewu(iorb,iorb,iorb,iorb))
! Transform
        ALLOCATE(MIDU1(iORB,iORB,iORB,iORB))
        ALLOCATE(MTR1(iORB,iORB))
        ALLOCATE(MTR2(iORB,iORB))
        ALLOCATE(MTR3(iORB,iORB))
        DO L=1,iORB
          DO K=1,iORB

                DO J=1,iORB; DO I=1,iORB
                MTR1(I,J)=tempU(I,J,K,L)
                END DO; END DO

!                call print_mat(iorb,iorb,MTR1) 
!                stop

        CALL DGEMM('T','N',iORB,iORB,iORB,1.0D0,TRANSFMATRIX,&
                   iORB,MTR1,iORB,0.0D0,MTR2,iORB)

        CALL DGEMM('T','T',iORB,iORB,iORB,1.0D0,TRANSFMATRIX,&
                   iORB,MTR2,iORB,0.0D0,MTR3,iORB)

                DO JJ=1,iORB; DO II=1,iORB
                MIDU1(II,JJ,K,L)=MTR3(II,JJ)
                END DO; END DO

          END DO
        END DO

!        call print_gat(iorb,iorb,iorb,iorb,MIDU1,12)
!        stop

        DO II=1,iORB
          DO JJ=1,iORB
                DO L=1,iORB; DO K=1,iORB
                MTR1(K,L)=MIDU1(II,JJ,K,L)
                END DO; END DO

        CALL DGEMM('N','N',iORB,iORB,iORB,1.0D0,MTR1,&
                   iORB,TRANSFMATRIX,iORB,0.0D0,MTR2,iORB)

        CALL DGEMM('T','N',iORB,iORB,iORB,1.0D0,MTR2,iORB,&
                    TRANSFMATRIX,iORB,0.0D0,MTR3,iORB)

                DO LL=1,iORB; DO KK=1,iORB
                tempNEWU(II,JJ,KK,LL)=MTR3(KK,LL)
                END DO; END DO

          END DO
        END DO

!        call print_gat(iorb,iorb,iorb,iorb,tempNEWU,14)
!        stop

        do i=1,iorb; do j=1,iorb; do k=1,iorb; do l=1,iorb
        newu(i,j,k,l)=tempnewu(i,k,l,j)
        end do; end do; end do; end do
        deallocate(tempu,tempnewu)
        DEALLOCATE(MTR1,MTR2,MTR3)
        DEALLOCATE(MIDU1)

        t=newt
        u=newu
!        call print_gat(iorb,iorb,iorb,iorb,U,16)
!        stop

        deallocate(newt,newu)

      END SUBROUTINE transform_RAM

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      SUBROUTINE GRAM_SCHMIDT(iORB,UMAT)
 
        INTEGER::iorb
        DOUBLE PRECISION::UMAT(iORB,iORB)
        DOUBLE PRECISION,ALLOCATABLE::OMAT(:,:)
        DOUBLE PRECISION,ALLOCATABLE::EV(:)
        DOUBLE PRECISION,ALLOCATABLE::TEMP_E(:)
        DOUBLE PRECISION,ALLOCATABLE::OMATRIX(:,:)
        DOUBLE PRECISION L

        ALLOCATE(OMAT(IORB,IORB),EV(IORB),TEMP_E(IORB))
        ALLOCATE(OMATRIX(IORB,IORB))


        OMATRIX=0.0
        EV=0.0
        OMAT=0.0

        DO i=1,iORB
          DO j=1,iORB
             OMATRIX(i,j)=UMAT(j,i)
             OMAT(i,j)=UMAT(j,i)
          END DO
        END DO

        CALL LENGTH(OMATRIX(1,:),iORB,L)

        DO i=1,iORB
          OMATRIX(1,i)=OMATRIX(1,i)/L
        END DO

        DO i=2,iORB
          TEMP_E=0.0
          DO j=1,i-1
            CALL GSM(OMAT(i,:),OMATRIX(j,:),iORB,EV)
            TEMP_E=TEMP_E+EV
          END DO
          OMATRIX(i,:)=OMAT(i,:)-TEMP_E
          CALL LENGTH(OMATRIX(i,:),iORB,L)
          OMATRIX(i,:)=OMATRIX(i,:)/L
        END DO

        OPEN(UNIT=84,file='EU.txt')
        DO i=1,iORB
          DO j=1,iORB
            UMAT(i,j)=OMATRIX(j,i)
            write(84,*)i,j,UMAT(i,j)
          END DO
        END DO

        DEALLOCATE(OMAT,EV,TEMP_E,OMATRIX)

      END SUBROUTINE GRAM_SCHMIDT

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      SUBROUTINE MODIFIED_GRAM_SCHMIDT(iORB,norb,UMAT)

        INTEGER::iorb,norb
        DOUBLE PRECISION::UMAT(iORB,nORB)

        DOUBLE PRECISION,ALLOCATABLE::OMAT(:,:)
        DOUBLE PRECISION,ALLOCATABLE::EV(:)
        DOUBLE PRECISION,ALLOCATABLE::TEMP_E(:)
        DOUBLE PRECISION,ALLOCATABLE::OMATRIX(:,:)
        DOUBLE PRECISION L

        ALLOCATE(OMAT(nORB,IORB),EV(IORB),TEMP_E(IORB))
        ALLOCATE(OMATRIX(nORB,IORB))

        OMATRIX=0.0
        EV=0.0
        OMAT=0.0

        DO i=1,nORB
          DO j=1,iORB
             OMATRIX(i,j)=UMAT(j,i)
             OMAT(i,j)=UMAT(j,i)
          END DO
        END DO

        CALL LENGTH(OMATRIX(1,:),iORB,L)

        DO i=1,iORB
          OMATRIX(1,i)=OMATRIX(1,i)/L
        END DO

        DO i=2,nORB
          TEMP_E=0.0
          DO j=1,i-1
            CALL GSM(OMAT(i,:),OMATRIX(j,:),iORB,EV)
            TEMP_E=TEMP_E+EV
          END DO
          OMATRIX(i,:)=OMAT(i,:)-TEMP_E
          CALL LENGTH(OMATRIX(i,:),iORB,L)
          OMATRIX(i,:)=OMATRIX(i,:)/L
        END DO

        DO i=1,iORB
          DO j=1,nORB
            UMAT(i,j)=OMATRIX(j,i)
          END DO
        END DO
 
        DEALLOCATE(OMAT,EV,TEMP_E,OMATRIX)

      END SUBROUTINE MODIFIED_GRAM_SCHMIDT

! only do ORTHOGONALIZATION
      SUBROUTINE MODIFIED_ORTHOGONALIZE(iORB,norb,UMAT)

        INTEGER::iorb,norb
        DOUBLE PRECISION::UMAT(iORB,nORB)
        DOUBLE PRECISION,ALLOCATABLE::OMAT(:,:)
        DOUBLE PRECISION,ALLOCATABLE::EV(:)
        DOUBLE PRECISION,ALLOCATABLE::TEMP_E(:)
        DOUBLE PRECISION,ALLOCATABLE::OMATRIX(:,:)
        DOUBLE PRECISION,ALLOCATABLE::UMATRIX(:,:)
        DOUBLE PRECISION L

        ALLOCATE(OMAT(nORB,IORB),EV(IORB),TEMP_E(IORB))
        ALLOCATE(OMATRIX(nORB,IORB))
        ALLOCATE(UMATRIX(nORB,IORB))

        OMATRIX=0.0
        EV=0.0
        OMAT=0.0

        DO i=1,nORB
          DO j=1,iORB
             OMATRIX(i,j)=UMAT(j,i)
             OMAT(i,j)=UMAT(j,i)
          END DO
        END DO

        CALL LENGTH(OMATRIX(1,:),iORB,L)

        DO i=1,iORB
          OMATRIX(1,i)=OMATRIX(1,i)/L
          UMATRIX(1,i)=OMATRIX(1,i)
        END DO

        DO i=2,nORB
          TEMP_E=0.0
          DO j=1,i-1
            CALL GSM(OMAT(i,:),OMATRIX(j,:),iORB,EV)
            TEMP_E=TEMP_E+EV
          END DO
          OMATRIX(i,:)=OMAT(i,:)-TEMP_E
!          CALL LENGTH(OMATRIX(i,:),iORB,L)
!          OMATRIX(i,:)=OMATRIX(i,:)/L
          UMATRIX(i,:)=OMATRIX(i,:)
        END DO

        DO i=1,iORB
          DO j=1,nORB
!            UMAT(i,j)=OMATRIX(j,i)
            UMAT(i,j)=UMATRIX(j,i)
          END DO
        END DO

        DEALLOCATE(OMAT,EV,TEMP_E,OMATRIX,UMATRIX)

      END SUBROUTINE MODIFIED_ORTHOGONALIZE 


      SUBROUTINE GSM(UM,RM,iORB,EV)

        DOUBLE PRECISION::UM(iORB)
        DOUBLE PRECISION::RM(iORB)
        DOUBLE PRECISION CE
        DOUBLE PRECISION::EV(IORB)

        CE=0.0
        EV=0.0

        DO i=1,iORB
           CE=CE+UM(i)*RM(i)
        END DO

        EV=CE*RM

      END SUBROUTINE GSM

! -------------------------------

      SUBROUTINE LENGTH(OM,iORB,LEN)

        DOUBLE PRECISION::OM(iORB)
        DOUBLE PRECISION LEN

        LEN=0.0
        DO i=1,iORB
          LEN=LEN+OM(i)**2
        END DO

        LEN=sqrt(LEN)

      END

! -------------------------------
    
      Subroutine matrix_copy(ndim1,ndim2,M1,M2)

        integer::ndim1,ndim2
        double precision::M1(ndim1,ndim2)
        double precision::M2(ndim1,ndim2)
     
        M2=M1 

      end subroutine matrix_copy

! -------------------------------

      Subroutine tensor_copy(ndim1,ndim2,ndim3,ndim4,M1,M2)

        integer::ndim1,ndim2,ndim3,ndim4
        double precision::M1(ndim1,ndim2,ndim3,ndim4)
        double precision::M2(ndim1,ndim2,ndim3,ndim4)
     
        M2=M1 

      end subroutine tensor_copy

! -------------------------------

      Subroutine inner_product(ndim,v1,v2,dv)

        integer::ndim
        double precision::v1(ndim)
        double precision::v2(ndim) 
        double precision dv 
        
        dv=0.0d0
        do i=1,ndim
          dv=dv+v1(i)*v2(i)
        end do

      End subroutine inner_product

! --- record the date and time

      Subroutine record_time(ifort,iadvance)
  
        integer::ifort,iadvance 
        character* 8 date
        character*10 time
        character* 5 zone
        integer      values(8)
 
        call date_and_time(date, time, zone, values)

        if(iadvance.eq.0)then
          write(ifort,"(80X, I4,A,I4,A,I4,A,I4,A,I4,A,I4,A,I4,A,I4,A)")&
          VALUES(1),"Y",VALUES(2),"M",VALUES(3),"D",VALUES(4),"U",&
          VALUES(5),"h",VALUES(6),"m",VALUES(7),"s",VALUES(8),"ms"
        else
          write(ifort,"(80X, I4,A,I4,A,I4,A,I4,A,I4,A,I4,A,I4,A,I4,A)",advance='no')&
          VALUES(1),"Y",VALUES(2),"M",VALUES(3),"D",VALUES(4),"U",&
          VALUES(5),"h",VALUES(6),"m",VALUES(7),"s",VALUES(8),"ms"
        end if
        call flush(ifort)

      End Subroutine 

      Subroutine directproductV(n1,n2,TM,v1,v2)

        integer::n1,n2
        double precision::TM(n1,n2)
        double precision::V1(n1)  
        double precision::V2(n2)  

        TM=0.0d0
        do i=1,n1
          do j=1,n2
            TM(i,j)=V1(i)*V2(j)
          end do
        end do

      end Subroutine directproductV  

      Subroutine get_position(ip,n,idx,pos,ipos,iflag)

        integer::ip,n,iflag      
        integer::idx(n),pos(n)

!        write(6,*)"idx",idx 
!        write(6,*)"pos",pos
        ipos=0 
        do i=1,n
          if(ip.eq.pos(i))then
            if(iflag.eq.0)then
              ipos=i
              exit                
            else
              ipos=idx(i)
              exit
            end if   
          end if
        end do  

      end Subroutine get_position

