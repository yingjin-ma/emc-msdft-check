      Subroutine Closed_shell_update(dm1,dm2,rdm1,rdm2,dCv)

        use global_control
        use matrix

        double precision::dCv

        double precision::dm1(nact,nact)
        double precision::dm2(nact,nact,nact,nact)

        double precision::rdm1(nocc,nocc)
        double precision::rdm2(nocc,nocc,nocc,nocc)

!        write(6,*)"nocc,nact",nocc,nact
!        call flush(6)

! CS means Closed Shell
        call CS_RDMs_update&
             (orb%nsub,nocc,nact,orb%closed,orb%act,orb%total,&
              dm1,dm2,rdm1,rdm2,orb%grouptable,dCv)

        call CS_dim_cls2fro()

      End Subroutine Closed_shell_update

      Subroutine CS_dim_cls2fro()

        use global_control
        use matrix

!        integer::n1
!        double precision::D1(n1,n1)
!        double precision::P1(n1,n1,n1,n1)

        nact=nocc

        orb%frozen=0 !cls -> fro treatment
        orb%closed=0
        orb%   occ=orbbk%   occ
        orb%   act=orbbk%   occ
        orb%  act2=orbbk%  act2
        orb% extnl=orbbk% extnl
        orb% total=orbbk% total
        orb%   ele=orbbk%   ele + orb%closed*2
        ncore=0
        ifcore=.false.
        ifcore2=ifcore

! The change of RDMs dimension (Separately)
!        deallocate(mat2%d)
!        deallocate(mat2%P)
!        allocate(mat2%d(nocc,nocc))
!        allocate(mat2%p(nocc,nocc,nocc,nocc))
!        mat2%d=D1
!        mat2%p=P1

! Nuclues  Energy
        FRONRE=FRONRE0

      End Subroutine CS_dim_cls2fro

      Subroutine CS_RDMs_update&
        (nsub,nocc,nact,cls,act,tot,D,P,D1,P1,group,dCv)

        integer::nsub,nocc,nact
        integer::cls(nsub),act(nsub),tot(nsub)
        integer::group(8,8)
        double precision::dCv
        double precision::D (nact,nact),P (nact,nact,nact,nact)
        double precision::D1(nocc,nocc),P1(nocc,nocc,nocc,nocc)

        double precision,allocatable::GM1(:,:,:,:)

        double precision dij,dil,dik,djl,dkl,djk
        double precision d12,d13,d23,d34,d14
        double precision dv

!        write(6,*) nsub,nocc,nact
!        call flush(6)

        D1=0.0d0
        P1=0.0d0

!        call print_mat(nocc,nocc,D1,6)
!        call flush(6)

        i0=0
        iC=0
        do i=1,nsub
          iC=iC+cls(i)
          D1(iC+1:iC+act(i),iC+1:iC+act(i))=&
          D (i0+1:i0+act(i),i0+1:i0+act(i))
          do ii=1,cls(i)
            do jj=1,cls(i)+act(i)
              if(ii.eq.jj)then
                D1(iC-cls(i)+ii,iC-cls(i)+ii)=2.0d0*(dCv)
              end if
            end do
          end do
          i0=i0+act(i)
          iC=iC+act(i)
        end do

!        call print_mat(nocc,nocc,D1,6)
!        print *,group
!        call flush(6)
!        stop

        allocate(GM1(nocc,nocc,nocc,nocc))
        GM1=0.0d0

        i0=0
        iC=0
        DO i=1,nsub
          iC=iC+cls(i)
          j0=0
          jC=0
          DO j=1,nsub
            jC=jC+cls(j)
            k0=0
            kC=0
            DO k=1,nsub
              kC=kC+cls(k)
              l0=0
              lC=0
              DO l=1,nsub
                lC=lC+cls(l)

                if(group(i,group(j,group(k,l))).eq.1)then

                  do ii=1,cls(i)+act(i)
                    do jj=1,cls(j)+act(j)
                      do kk=1,cls(k)+act(k)
                        do ll=1,cls(l)+act(l)

                          ia=ii+iC-cls(i)
                          ja=jj+jC-cls(j)
                          ka=kk+kC-cls(k)
                          la=ll+lC-cls(l)

                          dij=0.0d0
                          djk=0.0d0
                          dil=0.0d0
                          dkl=0.0d0
                          if(ii.eq.jj)dij=1.0d0
                          if(jj.eq.kk)djk=1.0d0
                          if(ii.eq.ll)dil=1.0d0
                          if(kk.eq.ll)dkl=1.0d0

                          ! closed shell approximation
                          ! p482 with additional anti-commutation
                          ! additonal deduced is needed for this part
                          ! talked with xiaoyu
                          if(ii.le.cls(i))then
                            ! i closed
                             GM1(ia,ja,ka,la)&
                            =2.0d0*dij*D1(ka,la)-dil*D1(ja,ka)
                          end if
                          if(jj.le.cls(j))then
                            ! j closed
                             GM1(ia,ja,ka,la)&
                            =2.0d0*dij*D1(ka,la)-djk*D1(ia,la)
                          end if
                          if(kk.le.cls(k))then
                          ! k closed  ! == i
                             GM1(ia,ja,ka,la)&
                            =2.0d0*dkl*D1(ia,ja)-djk*D1(ia,la)
                          end if
                          if(ll.le.cls(l))then
                          ! l closed  ! == j
                             GM1(ia,ja,ka,la)&
                            =2.0d0*dkl*D1(ia,ja)-dil*D1(ja,ka)
                          end if

                        end do
                      end do
                    end do
                  end do

                end if

                l0=l0+act(l)
                lc=lc+act(l)
              END DO
              k0=k0+act(k)
              kc=kc+act(k)
            END DO
            j0=j0+act(j)
            jc=jc+act(j)
          END DO
          i0=i0+act(i)
          ic=ic+act(i)
        END DO

        do i=1,nocc
          do j=1,nocc
            do k=1,nocc
              do l=1,nocc
                !P1(i,j,k,l)=GM1(i,l,j,k)
                P1(i,k,l,j)=GM1(i,j,k,l)
              end do
            end do
          end do
        end do

!        do i=1,nocc
!          do j=1,nocc
!            do k=1,nocc
!              do l=1,nocc
!                write(6,*)i,j,k,l,P1(i,j,k,l)
!              end do
!            end do
!          end do
!        end do

!        call print_gat(nocc,nocc,nocc,nocc,P1,6)
!        stop


! Fill with the all active part
        i0=0
        iC=0
        DO i=1,nsub
          iC=iC+cls(i)
          j0=0
          jC=0
          DO j=1,nsub
            jC=jC+cls(j)
            k0=0
            kC=0
            DO k=1,nsub
              kC=kC+cls(k)
              l0=0
              lC=0
              DO l=1,nsub
                lC=lC+cls(l)

                if(group(i,group(j,group(k,l))).eq.1)then

                  P1(iC+1:iC+act(i), &
                     jC+1:jC+act(j), &
                     kC+1:kC+act(k), &
                     lC+1:lC+act(l))=&
                  P (i0+1:i0+act(i), &
                     j0+1:j0+act(j), &
                     k0+1:k0+act(k), &
                     l0+1:l0+act(l))
                end if

                l0=l0+act(l)
                lc=lc+act(l)
              END DO
              k0=k0+act(k)
              kc=kc+act(k)
            END DO
            j0=j0+act(j)
            jc=jc+act(j)
          END DO
          i0=i0+act(i)
          ic=ic+act(i)
        END DO

!        call print_gat(nact,nact,nact,nact, P,6)
!        call print_gat(nocc,nocc,nocc,nocc,P1,6)

      End Subroutine CS_RDMs_update


      Subroutine CS_dim_fro2cls()

        use global_control
        use matrix

        orb%frozen=orbbk%frozen !fro -> cls treatment
        orb%closed=orbbk%closed
        orb%   occ=orbbk%   occ
        orb%   act=orbbk%   act
        orb%  act2=orbbk%  act2
        orb% extnl=orbbk% extnl
        orb% total=orbbk% total
        orb%   ele=orbbk%   ele

        ncore=0
        nact=0
        nocc=0
        do i=1,orb%nsub
          ncore=ncore+orb%closed(i)+orb%frozen(i)
          nact = nact+orb%act(i)
          nocc = nocc+orb%occ(i)
        end do
        ifcore=.true.
        ifcore2=ifcore

      End Subroutine CS_dim_fro2cls

