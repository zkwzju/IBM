      subroutine post_ellip
      use ellipsoid_m
      use common_m
      implicit none
      
      INTEGER n,l,l1,l2,i,j,k,i1,i2,j1,j2,k1,k2
      integer i11,i12,i21,i22,j11,j12,j21,j22,k11,k12,k21,k22
      INTEGER i_ddf,j_ddf,k_ddf,i_t,j_t
      integer i_ddf1,i_ddf2,j_ddf1,j_ddf2,k_ddf1,k_ddf2
      real*8  area,vol(nz)
!      real*8  aa(n_ll,n_ll),bb(n_ll)
      integer info

 
      if(ldebug) write(*,*) 'starting: post_ellip'
      
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(n,l,i1,i2,j1,j2,k1,k2,
C$OMP&    i,i_ddf,j,j_ddf,k)  
      do n=1,num_p
         do l=1,n_l(n)
            call del_fn(l,n)
            i1=p_iw(l,n)
            i2=p_ie(l,n)
            j1=p_js(l,n)
            j2=p_jn(l,n)
            k1=p_kb(l,n)
            k2=p_kt(l,n)

            do i=i1,i2
               i_ddf=i-i1+1
               if(i .lt. 1) then
                  i_t=nx+i
               elseif (i .gt. nx) then
                  i_t=i-nx
               else
                  i_t=i
               endif
               do j=j1,j2
                  j_ddf=j-j1+1
                  if(j .lt.1) then
                     j_t=ny+j
                  elseif(j .gt. ny) then
                     j_t=j-ny
                  else
                     j_t=j
                  endif
                  do k=k1,k2
                     ddf(l,n,i_ddf,j_ddf,k-k1+1)=ddf_dum(i_t,j_t,k)
                  enddo
               enddo
            enddo
         enddo
      enddo
C$OMP END PARALLEL DO

c     to calculate the Lagrangian volume

      area = deltax*deltay
      vol(1) = area*(zets(1)-zets(2))/2.d0
      vol(nz) = area*(zets(nz-1)-zets(nz))/2.d0
      do k=2,nz-1
         vol(k) = area*(zets(k-1)-zets(k+1))/2.d0
      end do

      
      if(lident) then
         do n=1,num_p
            select case (iweight)
            case(3)
               dv_l = 1.0d0
            case(4)
               call read_weights(dv_l,n_ll,num_p)
            end select
         enddo
      else
         select case (iweight)
         case(3)
            dv_l = 1.0d0
         case(4)
            call read_weights(dv_l,n_ll,num_p)
         case default
            write(*,*) 'Only allow reading weight from file'
            write(*,*) 'or const for non-identical particle'
            stop
         end select
      endif
      dv_l = dv_l * scale_dv      


 
!      call correlation(n_ll,vol)
!      call single_point_test

      if(ldebug) write(*,*) 'end: post_ellip'
      end subroutine post_ellip


! -------------------------------------------------- !
      subroutine read_weights(dv_l,n_ll,np)
      use ellipsoid_m,only : sname_w
      implicit none

      real*8 dv_l(n_ll,np)
      integer  n_ll,np
      integer istatus,l,n

      open(unit=170, file=sname_w, STATUS='OLD',iostat=istatus)
      if(istatus .ne. 0) then
         write(*,*) 'Problem to open file',sname_w
         stop
      end if

      do n=1,np
         do l=1,n_ll
            read(170,*,IOSTAT=istatus)  dv_l(l,n)
            if (istatus > 0)  then
               write(*,*) sname_w, 'cannot be opened'
               stop
            else if (istatus < 0) then ! end of file reached
               write(*,*) 'Too few weights given in file'
               stop
            end if
         end do
      end do
      close(170)
      end subroutine read_weights

 
