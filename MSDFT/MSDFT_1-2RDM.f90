!C*MODULE BLWCAS  *DECK OneAndTwoRDM 
      SUBROUTINE OneAndTwoRDM(iMETHOD,NATOM,NOrb,NAct,NConf,RDM1,RDM2, &
                 T,U,NuENG,IDFT,iSTATE,WEIGHTS,iTOLER,iROOT,GEOM,ATOMCHG)
      IMPLICIT NONE
      integer Norb,NAct,Nconf,IDFT,iSTATE,iMETHOD,iTOLER,iROOT
      integer NATOM
      integer IConfA(NConf,NAct),IConfB(NConf,NAct)
      INTEGER ATOMCHG(NATOM)
      REAL*8  GEOM(NATOM,3)
      real*8 WEIGHTS(iSTATE),tolerance
      real*8 Coeff(NConf,NConf),T(NOrb,NOrb),U(NOrb,NOrb,NOrb,NOrb)
      real*8 RDM1(NAct,NAct),RDM2(NAct,NAct,NAct,NAct)
      integer ITCon1(NAct),ITCon2(NAct)
      integer ITCon3(NAct),ITCon4(NAct)
      integer i,j,k,l,m,n,icount,i1,i2,i3,i4,NA,NB
      integer DETs(NConf,2,NAct),conf_index(iSTATE)
      real*8 Dsum(NConf,NConf)
      real*8 RDM1MN(NAct,NAct,NConf,NConf)
      real*8 RDM2MN(NAct,NAct,NAct,NAct,NConf,NConf)
      real*8 HH(NConf,NConf),NuENG,HD(NConf,NConf),XHF(NConf,NConf)
      real*8 EXC(NConf,NConf),DSIGN,DSIGN1,DSIGN2
      real*8,allocatable :: TMP1(:,:,:),TMP2(:,:,:,:,:)
      real*8 XHFJ(NConf,NConf),XHFK(NConf,NConf)
      integer IHFJ(NConf,NConf),IHFK(NConf,NConf)
      real*8 DFTFCA(NORB,NORB,NConf),DFTFCB(NORB,NORB,NConf)
      real*8 Energy_XC(Nconf),IFAC(Nconf,Nconf)

      tolerance=(10.0D0)**(-iTOLER)
      WRITE(8406,*)'tolerance:',tolerance
!C READ DETs
      OPEN(23,file='RDMDETs')
      READ(23,*)NA,NB
      Do i=1,NConf
        READ(23,*)(DETs(i,1,j),j=1,NA)
        READ(23,*)(DETs(i,2,j),j=1,NB)
      End Do
      CLOSE(23)
      RDM1MN=0.0D0
      RDM2MN=0.0D0
      
      IConfA=0
      IConfB=0
      Do i=1,NConf
        Do j=1,NA
          k=DETs(i,1,j)
          IConfA(i,k)=1
        end do
        Do j=1,NB
          k=DETs(i,2,j)
          IConfB(i,k)=1
        end do
      End do

!C Calculate 1-RDM
      do i=1,NAct
        do j=1,NAct
          Dsum=0.0D0
          do m=1,NConf
            do n=1,NConf
!C <Phi(m)|ai+aj-|Phi(n)>
             icount=0
             do i1=1,NAct
               icount=icount+ABS(IConfB(m,i1)-IConfB(n,i1))
             end do
             IF(icount==0) THEN
              ITCon1(:)=IConfA(m,:)
              ITCon2(:)=IConfA(n,:)
              IF(ITCon2(j)==1) then
                ITCon2(j)=0
              else
                goto 100
              end if
              IF(ITCon2(i)==0) then
                ITCon2(i)=1
              else
                goto 100
              end if
              icount=0
              do i1=1,NAct
                icount=icount+ABS(ITCon1(i1)-ITCon2(i1))
              end do
              IF(icount==0) then
                CALL RDM1SIGN(IConfA(n,:),NACT,I,J,DSIGN)
                Dsum(m,n)=Dsum(m,n)+1.0D0*DSIGN
              end if
             END IF
100   continue
!C  <Phi(m)|bi+bj-|Phi(n)>
             icount=0
             do i1=1,NAct
               icount=icount+ABS(IConfA(m,i1)-IConfA(n,i1))
             end do
             IF(icount==0) THEN
              ITCon1(:)=IConfB(m,:)
              ITCon2(:)=IConfB(n,:)
              IF(ITCon2(j)==1) then
                ITCon2(j)=0
              else
                goto 120
              end if
              IF(ITCon2(i)==0) then
                ITCon2(i)=1
              else
                goto 120
              end if
              icount=0
              do i1=1,NAct
                icount=icount+ABS(ITCon1(i1)-ITCon2(i1))
              end do
              IF(icount==0) then
                CALL RDM1SIGN(IConfB(n,:),NACT,I,J,DSIGN)
                Dsum(m,n)=Dsum(m,n)+1.0D0*DSIGN
              end if
             END IF
120   continue
            RDM1MN(i,j,m,n)=Dsum(m,n)
            end do
          end do
        end do
      end do

!C Calculate 2-RDM
      do i=1,NAct
        do j=1,NAct
          do k=1,NAct
            do l=1,NAct
              Dsum=0.0D0
              do m=1,NConf
                do n=1,NConf
                  Dsum=0.0D0
                  DSIGN=1.0D0
!C <Phi(m)|ai+ak+al-aj-|Phi(n)>
                  ITCon1(:)=IConfA(m,:)
                  ITCon2(:)=IConfA(n,:)
                  icount=0
                  do i1=1,NAct
                    icount=icount+ABS(IConfB(m,i1)-IConfB(n,i1))
                  end do
                  IF(icount==0) THEN
                    if(ITCon2(J)==1) then
                      ITCon2(J)=0
                    else
                      goto 200
                    end if
                    if(ITCon2(L)==1) then
                      ITCon2(L)=0
                    else
                      goto 200
                    end if
                    if(ITCon2(K)==0) then
                      ITCon2(K)=1
                    else
                      goto 200
                    end if
                    if(ITCon2(I)==0) then
                      ITCon2(I)=1
                    else
                      goto 200
                    end if
                    icount=0
                    do i1=1,NAct
                      icount=icount+ABS(ITCon1(i1)-ITCon2(i1))
                    end do
                    if(icount==0) then
                      CALL RDM2SIGN(IConfA(n,:),NACT,I,J,K,L,DSIGN)
                      Dsum(m,n)=Dsum(m,n)+1.0D0*DSIGN
                    end if
                  END IF
200   Continue

!C <Phi(m)|bi+bk+bl-bj-|Phi(n)>
                  ITCon1(:)=IConfB(m,:)
                  ITCon2(:)=IConfB(n,:)
                  icount=0
                  do i1=1,NAct
                    icount=icount+ABS(IConfA(m,i1)-IConfA(n,i1))
                  end do
                  IF(icount==0) THEN
                    if(ITCon2(J)==1) then
                      ITCon2(J)=0
                    else
                      goto 220
                    end if
                    if(ITCon2(L)==1) then
                      ITCon2(L)=0
                    else
                      goto 220
                    end if
                    if(ITCon2(K)==0) then
                      ITCon2(K)=1
                    else
                      goto 220
                    end if
                    if(ITCon2(I)==0) then
                      ITCon2(I)=1
                    else
                      goto 220
                    end if
                    icount=0
                    do i1=1,NAct
                      icount=icount+ABS(ITCon1(i1)-ITCon2(i1))
                    end do
                    if(icount==0) then
                      CALL RDM2SIGN(IConfB(n,:),NACT,I,J,K,L,DSIGN)
                      Dsum(m,n)=Dsum(m,n)+1.0D0*DSIGN
                    end if
                  END IF
220   Continue

!C <Phi(m)|ai+bk+bl-aj-|Phi(n)>
                  ITCon1(:)=IConfA(m,:)
                  ITCon2(:)=IConfA(n,:)
                  ITCon3(:)=IConfB(m,:)
                  ITCon4(:)=IConfB(n,:)
                  if(ITCon2(J)==1) then
                    ITCon2(J)=0
                  else
                    goto 240
                  end if
                  if(ITCon4(L)==1) then
                      ITCon4(L)=0
                  else
                    goto 240
                  end if
                  if(ITCon4(K)==0) then
                    ITCon4(K)=1
                  else
                    goto 240
                  end if
                  if(ITCon2(I)==0) then
                    ITCon2(I)=1
                  else
                    goto 240
                  end if
                  icount=0
                  do i1=1,NAct
                    icount=icount+ABS(ITCon1(i1)-ITCon2(i1))
                    icount=icount+ABS(ITCon3(i1)-ITCon4(i1))
                  end do
                  if(icount==0) then
                    CALL RDM1SIGN(IConfA(n,:),NACT,I,J,DSIGN1)
                    CALL RDM1SIGN(IConfB(n,:),NACT,K,L,DSIGN2)
                    DSIGN=DSIGN1*DSIGN2
                    Dsum(m,n)=Dsum(m,n)+1.0D0*DSIGN
                  end if
240   Continue

!C <Phi(m)|bi+ak+al-bj-|Phi(n)>
                  ITCon1(:)=IConfA(m,:)
                  ITCon2(:)=IConfA(n,:)
                  ITCon3(:)=IConfB(m,:)
                  ITCon4(:)=IConfB(n,:)
                  if(ITCon4(J)==1) then
                    ITCon4(J)=0
                  else
                    goto 260
                  end if
                  if(ITCon2(L)==1) then
                      ITCon2(L)=0
                  else
                    goto 260
                  end if
                  if(ITCon2(K)==0) then
                    ITCon2(K)=1
                  else
                    goto 260
                  end if
                  if(ITCon4(I)==0) then
                    ITCon4(I)=1
                  else
                    goto 260
                  end if
                  icount=0
                  do i1=1,NAct
                    icount=icount+ABS(ITCon1(i1)-ITCon2(i1))
                    icount=icount+ABS(ITCon3(i1)-ITCon4(i1))
                  end do
                  if(icount==0) then
                    CALL RDM1SIGN(IConfA(n,:),NACT,K,L,DSIGN1)
                    CALL RDM1SIGN(IConfB(n,:),NACT,I,J,DSIGN2)
                    DSIGN=DSIGN1*DSIGN2
                    Dsum(m,n)=Dsum(m,n)+1.0D0*DSIGN
                  end if
260   Continue
                RDM2MN(i,j,k,l,m,n)=Dsum(m,n)/2.0D0
                end do
              end do
            end do
          end do
        end do
      end do

! Calculate the CI Matrix
      HD=0.0D0
      XHF=0.0D0
      HH=0.0D0
      XHFJ=0.0D0
      XHFK=0.0D0
      IHFJ=0
      IHFK=0
      Do m=1,NConf
        do n=1,NConf
          do i=1,NAct
            do j=1,NAct 
              HD(m,n)=HD(m,n)+RDM1MN(i,j,m,n)*T(i,j)
            end do
          end do
          do i=1,NAct
            do j=1,NAct
              do k=1,NAct
                do l=1,NAct
                  XHF(m,n)=XHF(m,n)+RDM2MN(i,j,k,l,m,n)*U(i,j,k,l)
                  IF(iMETHOD>0) THEN
                  if(RDM2MN(i,j,k,l,m,n)>0.0D0) THEN
                    XHFJ(m,n)=XHFJ(m,n)+RDM2MN(i,j,k,l,m,n)*U(i,j,k,l)
                    IHFJ(m,n)=IHFJ(m,n)+RDM2MN(i,j,k,l,m,n)
                  else if(RDM2MN(i,j,k,l,m,n)<0.0D0) THEN
                    XHFK(m,n)=XHFK(m,n)+RDM2MN(i,j,k,l,m,n)*U(i,j,k,l)
                     IHFK(m,n)=IHFK(m,n)+RDM2MN(i,j,k,l,m,n)
                  endif 
                  ENDIF
                end do
              end do
            end do
          end do
          if(DABS(XHF(m,n)) < tolerance) THEN
            XHF(m,n)=0.0D0
          endif
          if(DABS(HD(m,n )) < tolerance) THEN
            HD(m,n) =0.0D0
          endif
          HH(m,n)=HD(m,n)+XHF(m,n)
          write(*,'(I4,I4,3X,F18.8,3X,F18.8)'),m,n,HD(m,n),XHF(m,n)
        end do
      end do
      Do i=1,NConf
        HH(i,i)=HH(i,i)+NuENG
      end do

! calculate DFT dEc
      EXC=0.0D0
      IF(IDFT==1) THEN
        WRITE(8406,*)'iMETHOD=',iMETHOD
        CALL DFT_Fc(NATOM,NORB,NCONF,DFTFCA,DFTFCB,Energy_XC, &
                    GEOM,ATOMCHG,0,Coeff,iROOT,iSTATE,WEIGHTS)
        CALL CALDEXCR(HH,HD,XHFJ,XHFK,NConf,EXC,iMETHOD,NORB,NA, &
                      NB,Energy_XC,IFAC)
        HH=HH+EXC
      ENDIF

      CALL MSDFTENG(HH,Coeff,NConf)
      ALLOCATE(TMP1(NAct,NAct,iROOT))
      ALLOCATE(TMP2(NAct,NAct,NAct,NAct,iROOT))

      IF(IDFT==1) THEN
        CALL DFT_Fc(NATOM,NORB,NCONF,DFTFCA,DFTFCB,Energy_XC, &
                    GEOM,ATOMCHG,1,Coeff,iROOT,iSTATE,WEIGHTS)
        CALL MSDFT_Fc(NORB,NCONF,Coeff,DFTFCA,DFTFCB,IFAC,EXC, &
                      Energy_XC(1),1,iROOT,iSTATE,WEIGHTS)
        
      ENDIF

      DO I1=1,iROOT
      RDM1=0.0D0
      RDM2=0.0D0
      do i=1,NAct
        do j=1,NAct
          do m=1,NConf
            do n=1,NConf
              RDM1(i,j)=RDM1(i,j)+Coeff(m,I1)*Coeff(n,I1)*RDM1MN(i,j,m,n)
            end do
          end do
          do k=1,NAct
            do l=1,NAct
              do m=1,NConf
              do n=1,NConf
              RDM2(i,j,k,l)=RDM2(i,j,k,l)+Coeff(m,I1)*Coeff(n,I1)*RDM2MN(i,j,k,l,m,n)
              end do
              end do      
            end do
          end do
        end do
      end do
      TMP1(:,:,I1)=RDM1(:,:)
      TMP2(:,:,:,:,I1)=RDM2(:,:,:,:)

      ENDDO
      write(8406,*)'WEIGHTS:'
      write(8406,*) WEIGHTS
      write(8406,*)'ROOTs AND STATES:'
      write(8406,*) iROOT,iSTATE
      RDM1=0.0D0
      RDM2=0.0D0
      conf_index=0
      IF(ABS(iROOT-iSTATE)>0) THEN
        call CONF_PURE(Coeff(:,1:iROOT),NConf,iROOT,iSTATE,conf_index)
      ENDIF
      DO I2=1,iSTATE
        IF(ABS(iROOT-iSTATE)>0) THEN
          I1=conf_index(I2)
        ELSE
          I1=I2 
        ENDIF
        RDM1(:,:)=RDM1(:,:)+TMP1(:,:,I1)*WEIGHTS(I2)
        RDM2(:,:,:,:)=RDM2(:,:,:,:)+TMP2(:,:,:,:,I1)*WEIGHTS(I2)
      ENDDO
      write(8406,*)'The target states are:'
      write(8406,*) conf_index
      DEALLOCATE(TMP1)
      DEALLOCATE(TMP2)
      END


      SUBROUTINE RDM1SIGN(IVEC,DIM,I1,J1,DSIGN)
      INTEGER I1,J1,DIM
      INTEGER IVEC(DIM)
      REAL*8 DSIGN
      INTEGER Ni,Nj,I,iCOUNT

      DSIGN=1.0D0      
      iCOUNT=0
      DO I=1,I1-1
        IF(IVEC(I)==1) THEN
          iCOUNT=iCOUNT+1
        ENDIF
      ENDDO
      Ni=iCOUNT

      iCOUNT=0
      DO I=1,J1-1
        IF(IVEC(I)==1) THEN
          iCOUNT=iCOUNT+1
        ENDIF
      ENDDO
      Nj=iCOUNT

      iCOUNT=0
      IF(I1<J1) THEN
        iCOUNT=Ni+Nj
      ELSEIF(I1>J1) THEN
        iCOUNT=Ni+Nj-1
      ELSEIF(I1==J1) THEN
        iCOUNT=Ni+Nj
      ENDIF
      
      DSIGN=(-1.0D0)**(iCOUNT)
      END




      SUBROUTINE RDM2SIGN(IVEC,DIM,I1,J1,K1,L1,DSIGN)
      INTEGER I1,J1,K1,L1,DIM
      INTEGER IVEC(DIM)
      REAL*8 DSIGN
      INTEGER Ni,Nj,Nk,Nl,I,iCOUNT

      DSIGN=1.0D0
      iCOUNT=0
      DO I=1,I1-1
        IF(IVEC(I)==1) THEN
          iCOUNT=iCOUNT+1
        ENDIF
      ENDDO
      Ni=iCOUNT

      iCOUNT=0
      DO I=1,J1-1
        IF(IVEC(I)==1) THEN
          iCOUNT=iCOUNT+1
        ENDIF
      ENDDO
      Nj=iCOUNT

      iCOUNT=0
      DO I=1,K1-1
        IF(IVEC(I)==1) THEN
          iCOUNT=iCOUNT+1
        ENDIF
      ENDDO
      Nk=iCOUNT

      iCOUNT=0
      DO I=1,L1-1
        IF(IVEC(I)==1) THEN
          iCOUNT=iCOUNT+1
        ENDIF
      ENDDO
      Nl=iCOUNT

      IF(L1>J1) Nl=Nl-1

      IF    (K1>J1 .AND. K1<=L1) THEN
        Nk=Nk-1
      ELSEIF(K1>J1 .AND. K1>L1) THEN
        Nk=Nk-2
      ELSEIF(K1<=J1 .AND. K1>L1) THEN
        Nk=Nk-1
      ENDIF

      IF    (I1>J1 .AND. I1>K1 .AND. I1>L1) THEN
        Ni=Ni-1
      ELSEIF(I1>J1 .AND. I1>K1 .AND. I1<=L1) THEN
        Ni=Ni
      ELSEIF(I1>J1 .AND. I1<=K1 .AND. I1>L1) THEN
        Ni=Ni-2
      ELSEIF(I1<=J1 .AND. I1>K1 .AND. I1>L1) THEN
        Ni=Ni
      ELSEIF(I1<=J1 .AND. I1<=K1 .AND. I1>L1) THEN
        Ni=Ni-1
      ELSEIF(I1>J1 .AND. I1<=K1 .AND. I1<=L1) THEN
        Ni=Ni-1
      ELSEIF(I1<=J1 .AND. I1>K1 .AND. I1<=L1) THEN
        Ni=Ni+1
      ENDIF
      iCOUNT=Ni+Nj+Nk+Nl

      DSIGN=(-1.0D0)**(iCOUNT)
      END

      SUBROUTINE CONF_PURE(COEFF,NDET,iROOT,iSTATE,iOUT)
      INTEGER NDET,iROOT,iSTATE
      INTEGER iOUT(iSTATE)
      INTEGER SS(iROOT),TT(iROOT)
      REAL*8 COEFF(NDET,iROOT),DSUM1,DSUM2
      INTEGER I,J,K1,K2,N1(NDET),N2(NDET),iCOUNT,NUM
      INTEGER iS,iT
      N1=0
      N2=0
      iCOUNT=0
      OPEN(23,FILE='couple_index.tmp')
      DO WHILE(.TRUE.)
        READ(23,*,END=601)I,J
        iCOUNT=iCOUNT+1
        N1(iCOUNT)=I
        N2(iCOUNT)=J
      ENDDO
601   CLOSE(23)
      IF(iCOUNT==0) THEN
        WRITE(*,*)'ERROR: CAN NOT FIND couple_index.tmp FILE'
        GOTO 999
      ENDIF
      NUM=iCOUNT

      iS=0
      iT=0
      DO I=1,iROOT
        DSUM1=0.0D0
        DO J=1,NUM
          K1=N1(J)
          K2=N2(J)
          IF(DABS(COEFF(K2,I))>1.0D-10) THEN
            DSUM1=DSUM1+COEFF(K1,I)/COEFF(K2,I)
          ENDIF
        ENDDO
        IF(DSUM1>0.5) THEN
          iS=iS+1
          SS(iS)=I
        ELSEIF(DSUM1<-0.5) THEN
          iT=iT+1
          TT(iT)=I
        ELSE
          WRITE(*,*)'ERROR IN FIND iROOT'
        ENDIF
      ENDDO

      IF(iS<iSTATE) THEN
        WRITE(*,*)'NOT ENOUGH ROOTS FOR iSTATES'
        GOTO 999
      ENDIF

      iOUT(:)=SS(1:iSTATE)
999   CONTINUE
      END







