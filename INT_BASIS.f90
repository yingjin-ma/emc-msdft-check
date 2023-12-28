subroutine BASIS_INT(ncenters,nconts,coord,atomchg,Hcore,TWOEIout,baselable)
implicit none
! input data zone
integer       :: ncenters,nconts
integer       :: atomchg(ncenters) ! natoms, multiplicity, charge
character     :: baselable*30
real(8)       :: coord(ncenters,3) ! natoms, coordinate
! output data zone
integer       :: iconv                   ! converge
integer,external :: INDEX_2E
integer       :: nRec
real(8) :: S(nconts,nconts),Hcore(nconts,nconts)
real(8),allocatable :: TWOEI(:)
integer       :: i,j,k,l
real(8)       :: TWOEIout(nconts,nconts,nconts,nconts)

nRec = index_2e(nconts-1,nconts-1,nconts-1,nconts-1)+1
print*,nconts,nRec
allocate(TWOEI(nRec))
call basisint(ncenters,coord,atomchg,baselable,nconts,S,Hcore,nRec,TWOEI)

do i = 1,nconts
   do j=1,nconts
      do k=1,nconts
         do l=1,nconts
             TWOEIout(i,j,k,l)  = TWOEI(INDEX_2E(i-1,j-1,k-1,l-1 )+1)
        enddo
      enddo
   enddo
enddo
deallocate(TWOEI)

print*,"mayj : done the BASIS_INT"

end subroutine

subroutine INTXC_ALLOCATE(N)
use MOL_info
implicit none
integer N
allocate(atoms(N))
end subroutine

subroutine INTXC_DEALLOCATE()
use MOL_info
implicit none
deallocate(atoms)
end subroutine


!===============================
subroutine basisint(ncenters,coord,atomchg,baselable,n1,Sout,Hcoreout,n2,TWOEIout)
! read basis set and integral
! build S matrix and 2-e integral
!===============================

use MOL_info
    implicit none
    integer       :: ncenters,atomchg(ncenters)  ! natoms, multiplicity, charge
    character     :: baselable*30      ! basis lable, functional
    real(8)       :: coord(ncenters,3) ! coordinate, atomcharge
    character     :: filename*50
    integer    :: i,j,k,l,n,m,start,info
    integer    :: ij,kl,ijkl,kl_max
    integer    :: n1,n2
    integer,allocatable  ::  ishls(:),jshls(:)
    real(8)    :: Sout(n1,n1),Hcoreout(n1,n1),TWOEIout(n2)
    integer    :: numShell,numGauss,nBases
    character  :: base_temp*30,basis_char*30
    integer    :: nblock1,nblock2,nblock1_bas,nblock2_bas
    real(8)    :: sum_tmp
    integer,allocatable :: atm(:,:),atm2(:,:),atm3(:,:)
    integer,allocatable :: bas(:,:)
    real(8),allocatable :: env(:)
    integer :: shls(4)
    real(8),allocatable :: buf1e(:,:), buf2e(:,:,:,:),buf1eV(:,:),buf1eT(:,:)
    integer :: offpoint 
    integer :: di, dj, dk, dl, nRec
    real(8),external :: CINTgto_norm
    integer,external :: CINTcgto_cart,CINTcgto_spheric
    integer,external :: INDEX_2E
    real(8),allocatable    ::   P1(:),P2(:)
    real(8) :: a,b
    real(8) :: dis
    integer   :: LWORK,LIWORK
    real(8),allocatable   ::  WORK(:),IWORK(:)
    real(8),allocatable   ::  S_e(:),S_temp(:,:),S_X(:,:)
    real(8),parameter   :: ans2bohr = 1.889725989
    real(8),parameter   :: PI = 3.14159265359
    real(8),parameter   :: covrad(26) = [0.38,0.32,1.34,0.9,0.82,0.77,0.75,&
                                     0.73,0.71,0.69,1.54,1.30,1.18,1.11,&
                                     1.06,1.02,0.99,0.97,1.96,1.74,1.44,& 
                                     1.36,1.25,1.27,1.39,1.25 ]

print*,'******'
print*,baselable
    filename=trim(baselable)//'.bas'
print*,filename
    open(unit=201,file=filename)

nconts = 0
nContssph = 0
Natoms = ncenters
!allocate(atoms(Natoms))
nRec=n2

! read coordinate

do i = 1,Natoms
  atoms(i)%coor = coord(i,:)
  atoms(i)%charge = atomchg(i) 
  write (atoms(i)%base,"(I0.2,A1,A27)")  atomchg(i),'-',baselable
enddo
!print *,atoms(1)%base



!  1  Read base informaion in
    nBases = 0 
    do i = 1,Natoms
       base_temp =""
       atoms(i)%nsShell = 0
       atoms(i)%npShell = 0
       atoms(i)%ndShell = 0
       atoms(i)%nfShell = 0
       atoms(i)%ngShell = 0
       do while(atoms(i)%base /= base_temp)
            read(201,*), base_temp
       enddo
print*,base_temp
       do while(.True.)
           read(201,*), basis_char,basis_char
           if(basis_char == "END") exit
           if(basis_char == "s") atoms(i)%nsShell = atoms(i)%nsShell+1
           if(basis_char == "p") atoms(i)%npShell = atoms(i)%npShell+1
           if(basis_char == "d") atoms(i)%ndShell = atoms(i)%ndShell+1
           if(basis_char == "f") atoms(i)%nfShell = atoms(i)%nfShell+1
           if(basis_char == "g") atoms(i)%ngShell = atoms(i)%ngShell+1
       enddo
       rewind(201)
       base_temp =""

       do while(atoms(i)%base /= base_temp)
            read(201,*), base_temp
       enddo

       numShell = atoms(i)%nsShell + atoms(i)%npShell + atoms(i)%ndShell + atoms(i)%nfShell + atoms(i)%ngShell
       nBases = nBases + numShell
       ! counts ncontractions Number
       atoms(i)%nconts = atoms(i)%nsShell + 3*atoms(i)%npShell + 6*atoms(i)%ndShell + 10*atoms(i)%nfShell + 15*atoms(i)%ngShell
       atoms(i)%nshell = numShell
   !!!!! 1.2 Read exponents and contrCoeff in Shell by shell
       allocate(atoms(i)%shell(numShell))
       ! 1.2.2 Read exponents and contrCoeff
       do j = 1,numShell
          Read(201,*),numGauss,basis_char
          atoms(i)%shell(j)%nGauss = numGauss
          if (basis_char == "s") atoms(i)%shell(j)%angMoment = 0
          if (basis_char == "p") atoms(i)%shell(j)%angMoment = 1
          if (basis_char == "d") atoms(i)%shell(j)%angMoment = 2
          if (basis_char == "f") atoms(i)%shell(j)%angMoment = 3
          if (basis_char == "g") atoms(i)%shell(j)%angMoment = 4
          allocate(atoms(i)%shell(j)%exponents(numGauss))
          allocate(atoms(i)%shell(j)%contrCoeff(numGauss))
          do k = 1,numGauss
             Read(201,*),atoms(i)%shell(j)%exponents(k),atoms(i)%shell(j)%contrCoeff(k)
          enddo
       enddo
       rewind(201)
    enddo 
    close(201)



!  2  Prepare the data 

    allocate(atm(6,nAtoms))
    allocate(atm2(6,nAtoms))
    allocate(atm3(6,nAtoms))
    allocate(bas(8,nBases))
    allocate(env(50*nBases+12*nAtoms))
    atm =0
    bas =0
    env =0
    offpoint = 1  ! start offpoint set to 1
    ! 2.1 Prepare Atm data.
    do i = 1,nAtoms
       atm(1,i) = atoms(i)%charge
       atm(2,i) = offpoint
       env(offpoint+1) = atoms(i)%coor(1)*ans2bohr
       env(offpoint+2) = atoms(i)%coor(2)*ans2bohr
       env(offpoint+3) = atoms(i)%coor(3)*ans2bohr
       offpoint = offpoint + 3
    enddo
    offpoint = offpoint +1
    ! 2.2 Prepare Basis data.
    start = 0
    do i =1,nAtoms
        numShell = atoms(i)%nsShell + atoms(i)%npShell + atoms(i)%ndShell + atoms(i)%nfShell + atoms(i)%ngShell
        do j = 1,numShell
           bas(1,j+start) = i-1
           bas(2,j+start) = atoms(i)%shell(j)%angMoment
           bas(3,j+start) = atoms(i)%shell(j)%nGauss
           bas(4,j+start) = 1
           ! 2.2.1 prepare exponents 
           bas(6,j+start) = offpoint
           do k = 1,atoms(i)%shell(j)%nGauss
               env(offpoint+k) = atoms(i)%shell(j)%exponents(k)
           enddo
           offpoint = offpoint + atoms(i)%shell(j)%nGauss

           ! 2.2.2 prepare contractcoeff
           bas(7,j+start) = offpoint
           do k = 1,atoms(i)%shell(j)%nGauss        
!!                 env(offpoint+k) = atoms(i)%shell(j)%contrCoeff(k)
              env(offpoint+k) = atoms(i)%shell(j)%contrCoeff(k)*CINTgto_norm( &
                                atoms(i)%shell(j)%angMoment, atoms(i)%shell(j)%exponents(k))
!              atoms(i)%shell(j)%contrCoeff(k) = atoms(i)%shell(j)%contrCoeff(k)*CINTgto_norm( &
!                                atoms(i)%shell(j)%angMoment, atoms(i)%shell(j)%exponents(k))
           enddo
           offpoint = offpoint + atoms(i)%shell(j)%nGauss

        enddo
        start  = start +numShell
    enddo


 do i = 1,nAtoms
    nConts = atoms(i)%nconts+nConts
 enddo


 do i = 1,nAtoms
    nContssph = nContssph + atoms(i)%nsShell + 3*atoms(i)%npShell &
              + 5*atoms(i)%ndShell + 7*atoms(i)%nfShell + 9*atoms(i)%ngShell
 enddo

nRec = index_2e(nconts-1,nconts-1,nconts-1,nconts-1)+1
print *,"nconts:   nRec:  "
print *,nConts, nRec, n2


allocate(NorVEC(nConts))
allocate(NorVECsph(nBases))


NorVEC = 0
!    Normalization 
do i = 0,nBases-1
  shls(1)=i
  shls(2)=i
  di = CINTcgto_cart(i, bas)
  ! Cartesian Normalization
  allocate (buf1e(di,di))
  call cint1e_ovlp_cart(buf1e, shls , atm, nAtoms,bas,nBases,env,0_8)
  call Normal(i,di,bas,nBases,buf1e,NorVEC,nConts)
  deallocate (buf1e)
  ! Spherical Normalization
  di = CINTcgto_spheric(i, bas)
  allocate (buf1e(di,di))
  call cint1e_ovlp_sph(buf1e,shls,atm,nAtoms,bas,nBases,env,0_8)
  NorVECsph(1+i) = dsqrt(1.0/buf1e(1,1)) 
  deallocate (buf1e)
enddo

! Absorb Shell Normalization factor into contcoeff
k = 1
do i =1,natoms
   do j = 1,atoms(i)%nshell
     ! print *,i,j
      atoms(i)%shell(j)%contrCoeff(:) =atoms(i)%shell(j)%contrCoeff(:)&
                                       *NorVECsph(k)
      k = k +1
   enddo
enddo
! 3.


!########## 1-electron integral    #######

allocate( S (nConts,nConts) ) 
allocate( Hcore (nConts,nConts))
S = 0
Hcore = 0
do i = 0,nBases-1 
   do j = i,nBases-1
      shls(1)=i
      shls(2)=j
      di = CINTcgto_cart(i, bas)
      dj = CINTcgto_cart(j, bas)
      allocate (buf1e(di,dj))
      allocate (buf1eT(di,dj))
      allocate (buf1eV(di,dj))
! 3.1 Computing Overlap matrix
      call cint1e_ovlp_cart(buf1e, shls , atm, nAtoms,bas,nBases,env,0_8)
      call store1e(shls,di,dj,bas,nBases,buf1e,NorVEC,nConts,S)
! 3.2 Computing G=T+V core Hamiltonian
      call cint1e_nuc_cart(buf1eV, shls , atm, nAtoms,bas,nBases,env,0_8)
      call cint1e_kin_cart(buf1eT, shls , atm, nAtoms,bas,nBases,env,0_8)
      call store1e(shls,di,dj,bas,nBases,buf1eT+buf1eV,NorVEC,nConts,Hcore)      
      deallocate (buf1e)
      deallocate (buf1eT)
      deallocate (buf1eV)
   enddo
enddo

do i = 1,nConts 
  do j = 1,i 
      write(6,*)"S and Hcore (",i,",",j,") = ", S(i,j),Hcore(i,j)
  end do 
end do 



!4.  
!##########   2-electron integral  ########
allocate(TWOEI(nRec))
allocate(ishls(nBases*nBases))
allocate(jshls(nBases*nBases))

ij = 0

do i = 0,nBases-1
    do j=0,i
       ishls(ij) = i
       jshls(ij) = j
       ij = ij +1
    enddo
enddo

do ij = 0, nBases*(nBases+1)/2-1
    i = ishls(ij)
    j = jshls(ij)
    di = CINTcgto_cart(i, bas)
    dj = CINTcgto_cart(j, bas)
    do kl = 0, ij
             k = ishls(kl)
             l = jshls(kl)
             dk = CINTcgto_cart(k, bas)
             dl = CINTcgto_cart(l, bas)
             shls(1) = i
             shls(2) = j
             shls(3) = k
             shls(4) = l
             allocate (buf2e(di,dj,dk,dl))
             ! integral based on shell
             call cint2e_cart(buf2e, shls, atm, nAtoms,bas,nBases,env,0_8)
             ! Store 2-e integral into TWOEI array 
             call store2e(shls,bas,buf2e,TWOEI,di,dj,dk,dl,nConts,nRec,nBases,NorVEC)
             deallocate(buf2e)
    enddo
enddo

do i = 1,nRec
  write(6,*)"U",TWOEI(i) 
end do 

deallocate(ishls)
deallocate(jshls)

Hcoreout =Hcore
Sout = S
TWOEIout = TWOEI

deallocate(S,Hcore,TWOEI)
deallocate(NorVECsph,NorVEC)
!deallocate(atoms)
deallocate(atm,atm2,atm3,bas,env)

write(6,*)" mayj : done basisint"
call flush(6)

end subroutine basisint
!=======================
subroutine Normal(iL,di,bas,nBasis,buf1e,Norfac,nConts)
!  Build Normalazation array Norfac
!=======================
    implicit none
    integer :: iL,di,nConts,nBasis
    real(8) :: buf1e(di,di)
    integer :: bas(8,nbasis)
    real(8) :: Norfac(nConts)

    real(8) :: Frac
    integer :: startnum
    integer :: i,j,k
    integer,external :: CINTcgto_cart
    ! Get StartNumber 
    startnum = 0 
    do i = 0,iL-1
       startnum = startnum + CINTcgto_cart(i, bas)
    enddo
    
    do i = 1,di
       Frac = dsqrt(1.0/buf1e(i,i))
       Norfac(startnum+i) = Frac 
    enddo
end subroutine


    integer function INDEX_2E(i,j,k,l)
    implicit none
    integer  ::   i,j,k,l
    integer  ::   ij,kl
    
    if (i>j) then
        ij = i*(i+1)/2 + j
    else
        ij = j*(j+1)/2 + i
    endif

    if (k>l) then
        kl = k*(k+1)/2 + l
    else
        kl = l*(l+1)/2 + k
    endif
  
    if (ij > kl) then 
        INDEX_2E = (ij+1)*ij/2 + kl
    else
        INDEX_2E = (kl+1)*kl/2 + ij
    endif
end function 


!=================
subroutine  store2e(shls,bas,buf2e,TWOEI,di,dj,dk,dl,nConts,nRec,nBases,Norfac)
!  Store 2-e integral into TWOEI
!================
implicit none
integer :: di, dj, dk, dl, nConts, nRec, nBases
integer :: shls(4)
integer :: bas(8,nBases)
real(8) :: buf2e(di,dj,dk,dl)
real(8) :: TWOEI(nRec)
real(8) :: Norfac(nConts)
real(8) :: intVal
integer,external :: CINTcgto_cart
integer,external :: INDEX_2E

integer :: is,js,ks,ls
integer :: i,j,k,l,num
integer :: e1,e2,e3,e4
integer :: ii


! 1. Get start Number of this shell
num = 0 
do ii = 0,shls(1)-1
   num = num + CINTcgto_cart(ii, bas)
enddo
is = num

num = 0 
do ii = 0,shls(2)-1
   num = num + CINTcgto_cart(ii, bas)
enddo
js = num

num = 0 
do ii = 0,shls(3)-1
   num = num + CINTcgto_cart(ii, bas)
enddo
ks = num

num = 0 
do ii = 0,shls(4)-1
   num = num + CINTcgto_cart(ii, bas)
enddo
ls = num


! 2. Store integral into TWOEI 

do i = 1,di
   do j = 1,dj
      do k = 1,dk
         do l = 1,dl
              e1 = is+i-1
              e2 = js+j-1
              e3 = ks+k-1
              e4 = ls+l-1
        !      if (e1.ge.e2 .and. e3.ge.e4 .and. &
        !          (e1+1)*e1/2+e2 .ge. (e3+1)*e3/2+e4)  then
              intVal = buf2e(i,j,k,l)*Norfac(e1+1)*Norfac(e2+1)*Norfac(e3+1)*Norfac(e4+1)
!              intVal = buf2e(i,j,k,l)
              TWOEI(INDEX_2E(e1,e2,e3,e4)+1) = intVal
!              print *,e1,e2,e3,e4,intVal
        !      endif
         enddo
      enddo
   enddo
enddo
end subroutine

!=========================
subroutine  store1e(shls,di,dj,bas,nBases,buf1e,Norfac,nConts,S)
! Store 1-e integral in to low-matrix 
!
!========================
implicit none
integer :: di, dj, nConts,  nBases
integer :: shls(4)
integer :: bas(8,nBases)
real(8) :: S(nConts,nConts)
real(8) :: Norfac(nConts)
real(8) :: buf1e(di,di)
integer,external :: CINTcgto_cart

real(8) :: intVal
integer :: is,js
integer :: i,j,num
integer :: e1,e2
integer :: ii

num = 0 
do ii = 0,shls(1)-1
   num = num + CINTcgto_cart(ii, bas)
enddo
is = num

num = 0 
do ii = 0,shls(2)-1
   num = num + CINTcgto_cart(ii, bas)
enddo
js = num

do i = 1,di
   do j = 1,dj
        e1 = is+i
        e2 = js+j
        if (e2 .ge. e1) then
          intVal = buf1e(i,j)*Norfac(e1)*Norfac(e2)
!          PRINT *,e1,e2,e2*(e2-1)/2+e1,intVal
          S(e1,e2) = intVal
          S(e2,e1) = intVal
        endif
   enddo
enddo
end subroutine

