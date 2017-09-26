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


        ALLOCATE(MTR1(iorb,iorb))
        MTR1=MATMUL(TRANSPOSE(TRANSFMATRIX),T)

        NEWT=MATMUL(MTR1,TRANSFMATRIX)
        DEALLOCATE(MTR1)


        allocate(tempu(iorb,iorb,iorb,iorb))
        do i=1,iorb; do j=1,iorb; do k=1,iorb; do l=1,iorb
        tempu(i,k,l,j)=u(i,j,k,l)
        end do; end do; end do; end do

        allocate(tempnewu(iorb,iorb,iorb,iorb))
        ALLOCATE(MIDU1(iORB,iORB,iORB,iORB))
        ALLOCATE(MTR1(iORB,iORB))
        ALLOCATE(MTR2(iORB,iORB))
        ALLOCATE(MTR3(iORB,iORB))
        DO L=1,iORB
          DO K=1,iORB

                DO J=1,iORB; DO I=1,iORB
                MTR1(I,J)=tempU(I,J,K,L)
                END DO; END DO



         MTR2=MATMUL(TRANSPOSE(TRANSFMATRIX),MTR1)


         MTR3=MATMUL(TRANSPOSE(TRANSFMATRIX),TRANSPOSE(MTR2))

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


         MTR2=MATMUL(MTR1,TRANSFMATRIX)


         MTR3=MATMUL(TRANSPOSE(MTR2),TRANSFMATRIX)

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

        t=newt
        u=newu

        deallocate(newt,newu)

      END SUBROUTINE transform_RAM

