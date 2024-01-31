      Subroutine SA_prec_mps()

        use global_control
        use matrix

        double precision dv
        double precision,allocatable::S(:,:,:)

        nroots=dmrg_nstates 

        allocate(S(nroots,nroots,nroots)); S=0.0d0

        do i0=1,nroots          
          Do i=1,nroots
            Do j=1,nroots
!              Do j=1,max(MPS(i)%nC,)

!              end do
            End Do
          End Do              
        end do


        deallocate(S)

      ! The preconditioner for mps/(Davidson vector) part

!        S=0.0d0
!        Do i=1,nroots
!          Do j=1,nroots
!            Do k=1,ncsf(State_Sym)
!              dnum=rdia(k)-Ene
!              dnum=Sign(Max(Abs(dnum),1.0d-16),dnum)
!              S(i+1,j+1)=S(i+1,j+1)+ &
!              CI(i*ncsf(State_Sym)+k)*CI(j*ncsf(State_Sym)+k)/dnum
!
!          write(6,*)i,j,k,CI(i*ncsf(State_Sym)+k),
!     &                    CI(j*ncsf(State_Sym)+k),dnum
!
!            End Do
!          End Do
!        End Do

      End Subroutine  

      Subroutine modified_redundant_residual(nij,nC,nredu,Oredu,YP) 

        integer::nij,nredu,nC
        integer::Oredu(nredu,2)
        double precision:: YP(nij+nC)
        double precision,allocatable::L1(:)

        nall=nij+nC
        neff=nij+nC-nredu

        allocate(L1(nall)); L1=0.0d0

        ii=0
        do i=1,nall
          iflag=0
          do i2=1,nredu
            if(i.eq.Oredu(i2,1))then
              iflag=1
            end if
          end do
          if(iflag.ne.1)then
            ii=ii+1
            ! Rotation (R)
            L1(i)=YP(ii)
          end if
        end do

        do i=1,nredu
          L1(Oredu(i,1))=-1.0d0*YP(Oredu(i,2))
        end do

        write(1,*)" residual (with redundant)-- "
!        call print_mat(1,nall,L1,1)
        call flush(1)

        ! Updated residual
        YP=L1

        deallocate(L1)

      end Subroutine modified_redundant_residual

      Subroutine rm_redundant(nij,nC,nredu,Oredu,YP,YP2)

        integer::nij,nredu,nC
        integer::Oredu(nredu,2)
        double precision:: YP(nij+nC),YP2(nij+nC-nredu)
        double precision,allocatable::L1(:)

        nall=nij+nC
        neff=nij+nC-nredu

        ! get the non-redundant 
        allocate(L1(neff)); L1=0.0d0

        ii=0
        do i=1,nall
          iflag=0
          do i2=1,nredu
            if(i.eq.Oredu(i2,1))then
              iflag=1
            end if
          end do
          if(iflag.ne.1)then
            ii=ii+1
            ! with redundant (YP)
            L1(ii)=YP(i)
          end if
        end do
    
        YP2=0.0d0 
        write(1,*)"The non-redundant, as YP2"
!        call print_mat(1,neff,L1,1)
        YP2=L1

        deallocate(L1)
 
      End Subroutine rm_redundant

      Subroutine RminusAx(nij,nC,LP,nredu,Oredu,YP,rvec,ralpha)

        integer::nij,nredu,nC
        integer::Oredu(nredu,2)
        double precision:: LP(nij+nC)
        double precision:: YP(nij+nC)

        double precision,allocatable::L1(:),L2(:)
        double precision::rvec(nij+nC-nredu)
        double precision::ralpha  

        nall=nij+nC
        neff=nij+nC-nredu       

        !write(1,*)"Ax" 
        !call print_mat(1,nall,LP,1)
        !write(1,*)"G/R"
        !call print_mat(1,nall,YP,1)
 
        allocate(L1(nall)); L1=0.0d0

        ii=0
        do i=1,nall
          iflag=0
          do i2=1,nredu
            if(i.eq.Oredu(i2,1))then
              iflag=1
            end if
          end do
          if(iflag.ne.1)then
            ii=ii+1
            ! Rotation (R)
            L1(i)=LP(ii)
          end if          
        end do  

        do i=1,nredu
          L1(Oredu(i,1))=-1.0d0*LP(Oredu(i,2)) 
        end do

        write(1,*)" AX (with redundant)-- " 
!        call print_mat(1,nall,L1,1) 
        call flush(1)

        ! Ax  
        LP=L1
        ! Updated residual
        YP=YP-LP*ralpha

        ! Save the non-redundant residuals
        allocate(L2(neff)); L2=0.0d0 

        ii=0
        do i=1,nall
          iflag=0
          do i2=1,nredu
            if(i.eq.Oredu(i2,1))then
              iflag=1
            end if
          end do
          if(iflag.ne.1)then
            ii=ii+1
            ! residual (Y)
            L2(ii)=YP(i)
          end if
        end do

        write(1,*)"The non-redundant residual, as rvec"
!        call print_mat(1,neff,L2,1)
        !r for orb 
        rvec=L2

        deallocate(L1,L2)

      End Subroutine RminusAx

      Subroutine Hall_Rall(nij,nC,HP2,LP,YP,nredu,Oredu,sigma,xvec)

        integer::nij,nredu,nC
        integer::Oredu(nredu,2)
        double precision:: LP(nij+nC)
        double precision:: YP(nij+nC)
        double precision::HP2(nij+nC,nij+nC)
 
!        double precision ralpha,res,dv
        double precision,allocatable::L1(:),L2(:) 
        double precision,allocatable::TM1(:,:) 

        double precision::sigma(nij+nC-nredu)
        double precision:: xvec(nij+nC-nredu)

        nall=nij+nC
        neff=nij+nC-nredu       

        write(1,*)"nall",nall 
        write(1,*)"neff",neff
        write(1,*)"nredu, Oredu ",nredu, " :: ",Oredu

        allocate( L1(neff)) 
        allocate( L2(neff)) 
        allocate(TM1(neff,neff)) 

         L1=0.0d0
         L2=0.0d0
        TM1=0.0d0

        ii=0
        do i=1,nall

          iflag=0
          do i2=1,nredu
            if(i.eq.Oredu(i2,1))then 
              iflag=1
            end if
          end do
          if(iflag.ne.1)then
            ii=ii+1
            ! Rotation (R)
            L1(ii)=LP(i)
          end if

          jj=0 
          do j=1,nall

            jflag=0
            do j2=1,nredu
              if(j.eq.Oredu(j2,1))then
                jflag=1
              end if
            end do
            if(jflag.ne.1)jj=jj+1

            if(iflag.ne.1.and.jflag.ne.1)then
!              write(1,*)i,j,"",ii,jj
              ! Hessian (H)
              TM1(ii,jj)=HP2(i,j)
            end if

          end do
        end do

        write(1,*)"The non-redundant rotations, as lagrange (xvec)"
!        call print_mat(1,neff,L1,1)
        !X for orb 
        xvec=L1
        !call inner_product(neff,L1,L1,res)          

        write(1,*)" ----> The effective Hessian <----"
!        call print_mat(neff,neff,TM1,1)

        call MXV(neff,neff,TM1,L1,L2)
        write(1,*)"The updated Ax (non-redundant, as sigma)"
!        call print_mat(1,neff,L2,1)
       
        ! Ax as sigma
        sigma=L2  
        ! <x|A|x> for orb 
        !call inner_product(neff,L1,L1,ralpha)          

        Yp=0.0d0
        YP(1:neff)=L2(1:neff)

        deallocate(L1,L2,TM1) 

      End Subroutine Hall_Rall


