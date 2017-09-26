!C*MODULE BLWCAS  *DECK GEN_DETMO
      SUBROUTINE GEN_DETMO(Norb,NDet)
!C Generate CASSCFMO
      Real*8 COEFF(Norb,Norb)
      Integer DETs(NDet,2,Norb)     
      Integer i,j,k,l,NA,NB

!C READ DETs
      OPEN(23,file='DETs')
      READ(23,*)NA,NB
      Do i=1,NDet
        READ(23,*)(DETs(i,1,j),j=1,NA)
        READ(23,*)(DETs(i,2,j),j=1,NB)
      End Do
      CLOSE(23)

!C READ Initial MOs
      OPEN(23,file='INT_MO')
      Do i=1,Norb
        READ(23,*)(COEFF(j,i),j=1,Norb)
      End Do
      CLOSE(23)

!C WRITE CASSCFMO
      OPEN(24,file='CASSCFMO')
      Do i=1,NDet
        Do j=1,NA
          k=DEts(i,1,j)
          WRITE(24,*)COEFF(:,k)
        End do
        Do j=1,NB
          k=DEts(i,2,j)
          WRITE(24,*)COEFF(:,k)
        End do
      End do
      CLOSE(24)
      END

