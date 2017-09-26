          subroutine casconf(n_orbital,n_alpha,n_beta,n_froz,n_closed,icount,ierro)
          implicit none
          integer n_orbital,n_alpha,n_beta,n_froz,n_closed,ierro
          integer a(n_orbital,99999)
          integer b(n_orbital,99999)
          integer ica,icb,it,i,j,k,icount
          integer vec(n_orbital),vocc(n_orbital+n_froz+n_closed)

          ierro=0
          call comb(n_orbital,n_alpha,a,ica)
          call comb(n_orbital,n_beta,b,icb)
          do i=1,n_froz+n_closed
            vocc(i)=i
          enddo
          open(23,file='DETs')
! need to improve
          write(23,*)n_alpha+n_froz+n_closed,n_beta+n_froz+n_closed
!          write(23,*)n_alpha,n_beta

          it=0
          do i=1,ica
          do j=1,icb
          it=it+1
          icount=0
          do k=1,n_orbital
           if(a(k,i).eq.1) then
             icount=icount+1
             vec(icount)=k+n_froz+n_closed
           endif
          enddo
          if(icount/=n_alpha) then
            ierro=-1
            write(8406,*)'ERROR IN CONF: ALPHA'
            goto 201
          endif
          vocc(n_froz+n_closed+1:n_froz+n_closed+n_alpha)=vec(1:n_alpha)
          write(23,*)(vocc(k),k=1,n_froz+n_closed+n_alpha)
          icount=0
          do k=1,n_orbital
           if(b(k,j).eq.1) then
             icount=icount+1
             vec(icount)=k+n_froz+n_closed
           endif
          enddo
          if(icount/=n_beta) then
            ierro=-1
            write(8406,*)'ERROR IN CONF: BETA'
            goto 201
          endif
          vocc(n_froz+n_closed+1:n_froz+n_closed+n_beta)=vec(1:n_beta)
          write(23,*)(vocc(k),k=1,n_froz+n_closed+n_beta)

          end do
          end do
          close(23)
          icount=ica*icb
201       continue
          end

          subroutine comb(from,choose,combination,ic)
          ! Ref:
          ! https://liam0205.me/2016/01/31/binomial-in-cpp/
          ! From `from' items select `choose' ones
          ! All the possible results will be store in `combination'
          ! the number of combination is `ic'
          implicit none
          integer :: from, choose
          integer :: work(from)
          integer ::combination(from,99999)
          integer :: i
          logical,external :: found10
          integer :: pos
          integer ::num0, num1,j,ic
          ! initialization
          ic=1
          work(1:choose)=1
          work(choose+1:from)=0
          !write(*,*)work
          combination(:,ic)=work
          do while(found10(work,from,pos))
          ic=ic+1
          work(pos)=0
          work(pos+1)=1
          if(pos.gt.1)then
              num1=0
              num0=0
              do j=1,pos-1
              if(work(j).eq.1) num1=num1+1
              !if(work(j).eq.0) num0=num0+1
              end do
              work(1:num1)=1
              work(num1+1:pos-1)=0
          end if

          !write(*,*)work
          combination(:,ic)=work
          end do
          return
          end subroutine

          logical function found10(work,from,pos)
          implicit none
          integer:: from, pos,i
          integer:: work(from)
          found10=.false.
          do i=1,from-1
          if(work(i).eq.1.and.work(i+1).eq.0) then
              pos=i
              found10=.true.
              exit
          end if
          end do
          return
          end function

          subroutine find_couple(icount,n_orbital,ierro)
          integer na,nb,icount,n_orbital,ierro
          integer veca(n_orbital,icount)
          integer vecb(n_orbital,icount)
          integer t1a(n_orbital),t1b(n_orbital)
          integer t2a(n_orbital),t2b(n_orbital)
          integer i,j,k,suma,sumb

          ierro=0
          open(23,file='DETs')
          veca=0
          vecb=0
          read(23,*)na,nb
          do i=1,icount
            read(23,*)(veca(k,i),k=1,na)
            read(23,*)(vecb(k,i),k=1,nb)
          enddo
          close(23)
          open(24,file='couple_index.tmp')
          do i=1,icount
            t1a=0
            t1b=0
            t1a(:)=veca(1:na,i)
            t1b(:)=vecb(1:nb,i)
            do j=i+1,icount
              t2a=0
              t2b=0
              t2a(:)=veca(1:na,j)
              t2b(:)=vecb(1:nb,j)
              suma=0
              sumb=0
              do k=1,n_orbital
                suma=suma+ABS(t2a(k)-t1b(k))
                sumb=sumb+ABS(t2b(k)-t1a(k))
              enddo
              if(suma==0 .and. sumb==0) then
                write(24,*),i,j
              endif
            enddo
          enddo
          close(24)
          end






