      SUBROUTINE ibm_ellip(cfnl2,cfdif)
      use ellipsoid_m, only: for_px,for_py,for_pz,
     &     torq_x,torq_y,torq_z,lcheck,lconvect
      use p_dyn_m
      use common_m
      use ellip_common_m, only : lunformatted,lcheckmarker,nmarker,
     &     marker_ind
      implicit none
c
c This subroutine creates the forcing function (forcing_xyz)
c that will be added to the rhs of the helmholtz equation 
c to account for the immersed body (based on H. Lee)
c
      INTEGER n,l,i,j,k,i1,i2,j1,j2,k1,k2,nn_l,n_l_refu,
     &           n_l_refv,n_l_refw
      REAL*8 u_tilde(nx0,ny0,nz),v_tilde(nx0,ny0,nz)
      REAL*8 w_tilde(nx0,ny0,nz)
      REAL*8 error_u,error_v,error_w
      REAL*8 max_error_u,max_error_v,max_error_w
      REAL*8 up_tilde(n_ll,num_p),vp_tilde(n_ll,num_p)
      REAL*8 wp_tilde(n_ll,num_p)
      REAL*8 u_forcing(n_ll,num_p)
      REAL*8 v_forcing(n_ll,num_p),w_forcing(n_ll,num_p)
      REAL*8 cfnl2,cfdif
      REAL*8 c_d(num_p),c_d_dum(num_p),c_d_dumm(num_p),c_dd(num_p),
     &       c_lx(num_p),c_lz(num_p),c_lx_dum(num_p),c_lz_dum(num_p)
     &       ,c_d_dum_test(num_p)
c      REAL*8 u_dum(nx0,ny0,nz),v_dum(nx0,ny0,nz),w_dum(nx0,ny0,nz)
c      REAL*8 h1_dum(nx0,ny0,nz),h2_dum(nx0,ny0,nz),h3_dum(nx0,ny0,nz)
c      REAL*8 ru_dum(nx0,ny0,nz),rv_dum(nx0,ny0,nz),rw_dum(nx0,ny0,nz)
      REAL*8 testfft,testfft2
      REAL*8 u_p_n(n_ll,num_p),v_p_n(n_ll,num_p)
      REAL*8 w_p_n(n_ll,num_p)
      REAL*8 u_c_n(num_p),v_c_n(num_p),w_c_n(num_p)
      REAL*8 x_c_n(num_p),y_c_n(num_p),z_c_n(num_p)
      REAL*8 dum_p
      REAL*8 u_rot,v_rot,w_rot
      REAL*8 c_d_cluster,c_lx_cluster,c_lz_cluster
      INTEGER i_ddf,j_ddf,i_t,j_t,i_f,j_f
c      REAL*8 p_dum(nz,ny,nx),dp_dum(nz,ny,nx)
      REAL*8 axis1(3),axis2(3),axis3(3)
      real*8 temp_var
c--------------------------------------------------------------------
c--------------------------------------------------------------------
c initializing the variables to zero upon every entry for 
c a new time step
c
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(l,n)
      do n=1,num_p
         do l=1,n_ll
            up_tilde(l,n)=0.d0
            vp_tilde(l,n)=0.d0
            wp_tilde(l,n)=0.d0
            u_forcing(l,n)=0.d0
            v_forcing(l,n)=0.d0
            w_forcing(l,n)=0.d0
         enddo
      enddo
C$OMP END PARALLEL DO
c
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do k=1,nz
         do j=1,ny0
            do i=1,nx0
               forcing_x(i,j,k)=0.d0
               forcing_y(i,j,k)=0.d0
               forcing_z(i,j,k)=0.d0
            enddo
         enddo
      enddo


C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do k=1,nz
         do j=1,ny0
            do i=1,nx0
               u_tilde(i,j,k)=u(i,j,k)
     &                        +cfnl2*h1(i,j,k)+2.d0*cfdif*ru(i,j,k)
     &                                       -2.d0*cfdif*dpx(i,j,k)
               v_tilde(i,j,k)=v(i,j,k)
     &                        +cfnl2*h2(i,j,k)+2.d0*cfdif*rv(i,j,k)
     &                                       -2.d0*cfdif*dpy(i,j,k)
               w_tilde(i,j,k)=w(i,j,k)
     &                        +cfnl2*h3(i,j,k)+2.d0*cfdif*rw(i,j,k)
     &                                       -2.d0*cfdif*dpz(i,j,k)
!>>>
               if(lcheck) then
                  if(lconvect) then
                     forcing_x03(i,j,k) = (
     &                    +cfnl2*h1(i,j,k)
     &                    +2.d0*cfdif*ru(i,j,k)
     &                    -2.d0*cfdif*dpx(i,j,k))/dt
                     forcing_y03(i,j,k) = (
     &                    +cfnl2*h2(i,j,k)
     &                    +2.d0*cfdif*rv(i,j,k)
     &                    -2.d0*cfdif*dpy(i,j,k))/dt
                     forcing_z03(i,j,k) = (
     &                    +cfnl2*h3(i,j,k)
     &                    +2.d0*cfdif*rw(i,j,k)
     &                    -2.d0*cfdif*dpz(i,j,k))/dt

                  else
                     forcing_x03(i,j,k) = (
!     &                    +cfnl2*h1(i,j,k)
     &                    +2.d0*cfdif*ru(i,j,k)
     &                    -2.d0*cfdif*dpx(i,j,k))/dt
                     forcing_y03(i,j,k) = (
!     &                    +cfnl2*h2(i,j,k)
     &                    +2.d0*cfdif*rv(i,j,k)
     &                    -2.d0*cfdif*dpy(i,j,k))/dt
                     forcing_z03(i,j,k) = (
!     &                    +cfnl2*h3(i,j,k)
     &                    +2.d0*cfdif*rw(i,j,k)
     &                    -2.d0*cfdif*dpz(i,j,k))/dt
                  endif
               endif
            enddo
         enddo
      enddo
C$OMP END PARALLEL DO

!>>>  
      if(lcheck) then
         call fft_f2r_2d_new(forcing_x03)
         call fft_f2r_2d_new(forcing_y03)
         call fft_f2r_2d_new(forcing_z03)
         if(irkk .eq.1) then
            forcing_x01=forcing_x03
            forcing_y01=forcing_y03
            forcing_z01=forcing_z03
         elseif(irkk .eq. 2) then
            forcing_x02=forcing_x03
            forcing_y02=forcing_y03
            forcing_z02=forcing_z03
         endif
      endif
!<<<         
      
      
c      This section follows the flow solver of M. Uhlmann (2005) to
c     define the forcing function.
c
cc  Transforming the variables to real space
            call fft_f2r_2d_new(u_tilde)
            call fft_f2r_2d_new(v_tilde)
            call fft_f2r_2d_new(w_tilde)
          
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(l,n,i1,i2,j1,j2,k1,k2
C$OMP& ,i,j,k,i_t,j_t,i_ddf,j_ddf)
      do n=1,num_p
         do l=1,n_l(n)
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
                     up_tilde(l,n)=up_tilde(l,n)+u_tilde(i_t,j_t,k)
     &                    *ddf(l,n,i_ddf,j_ddf,k-k1+1)
                     vp_tilde(l,n)=vp_tilde(l,n)+v_tilde(i_t,j_t,k)
     &                    *ddf(l,n,i_ddf,j_ddf,k-k1+1)
                     wp_tilde(l,n)=wp_tilde(l,n)+w_tilde(i_t,j_t,k)
     &                    *ddf(l,n,i_ddf,j_ddf,k-k1+1)
                  enddo
               enddo
            enddo

            u_forcing(l,n)=(u_p(l,n)-up_tilde(l,n))/dt
            v_forcing(l,n)=(v_p(l,n)-vp_tilde(l,n))/dt
            w_forcing(l,n)=(w_p(l,n)-wp_tilde(l,n))/dt
         enddo
      enddo
C$OMP END PARALLEL DO

      do n=1,num_p
         for_px(n)=0.d0
         for_py(n)=0.d0
         for_pz(n)=0.d0
      enddo
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(l,n,temp_var)
      do n=1,num_p
         do l=1,n_l(n)
            temp_var = marker_v(l,n)*dv_l(l,n)
            for_px(n)=for_px(n)+u_forcing(l,n)*temp_var   
            for_py(n)=for_py(n)+v_forcing(l,n)*temp_var
            for_pz(n)=for_pz(n)+w_forcing(l,n)*temp_var
         enddo

      enddo
C$OMP END PARALLEL DO
!>>>
      if(lcheckmarker) then
!     particle velocity difference temperart output
         do i=1,nmarker
            n = marker_ind(i,1)
            l = marker_ind(i,2)
            write(355,'(I10,I6,3E15.7)') itime,i,up_tilde(n,l),
     &           vp_tilde(n,l),wp_tilde(n,l)
         enddo
      endif
!<<<
      do n=1,num_p
         torq_x(n)=0.d0
         torq_y(n)=0.d0
         torq_z(n)=0.d0
      enddo

C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(l,n,temp_var)
      do n=1,num_p
         do l=1,n_l(n)
            temp_var = marker_v(l,n)*dv_l(l,n)
            torq_x(n)=torq_x(n)+(ry_l(l,n)*w_forcing(l,n)
     &           -rz_l(l,n)*v_forcing(l,n))*temp_var
            torq_y(n)=torq_y(n)+(rz_l(l,n)*u_forcing(l,n)
     &           -rx_l(l,n)*w_forcing(l,n))*temp_var
            torq_z(n)=torq_z(n)+(rx_l(l,n)*v_forcing(l,n)
     &           -ry_l(l,n)*u_forcing(l,n))*temp_var
         enddo
      enddo
C$OMP END PARALLEL DO

c Caculating the forcing term on the Eulerian grid
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(l,n,i1,i2,j1,j2,k1,k2
C$OMP& ,i,j,k,i_f,j_f,i_ddf,j_ddf)
      do n=1,num_p
         do l=1,n_l(n)
            i1=p_iw(l,n)
            i2=p_ie(l,n)
            j1=p_js(l,n)
            j2=p_jn(l,n)
            k1=p_kb(l,n)
            k2=p_kt(l,n)

            do i=i1,i2
             i_ddf=i-i1+1
             if(i .lt. 1) then
                i_f=nx+i
             elseif (i .gt. nx) then
                i_f=i-nx
             else
                i_f=i
             endif
             do j=j1,j2
              j_ddf=j-j1+1
              if(j .lt.1) then
                 j_f=ny+j
              elseif(j .gt. ny) then
                 j_f=j-ny
              else
                 j_f=j
              endif
              do k=k1,k2
               forcing_x(i_f,j_f,k)=forcing_x(i_f,j_f,k)+u_forcing(l,n)
     &                           *ddf(l,n,i_ddf,j_ddf,k-k1+1)*dv_l(l,n)
               forcing_y(i_f,j_f,k)=forcing_y(i_f,j_f,k)+v_forcing(l,n)
     &                           *ddf(l,n,i_ddf,j_ddf,k-k1+1)*dv_l(l,n)
               forcing_z(i_f,j_f,k)=forcing_z(i_f,j_f,k)+w_forcing(l,n)
     &              *ddf(l,n,i_ddf,j_ddf,k-k1+1)*dv_l(l,n)

              enddo
             enddo
            enddo
         enddo
      enddo
C$OMP END PARALLEL DO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
 100  format(2(1x,I8),6(3x,g15.8))
c Lagrangian based approximation for drag
      do n=1,num_p
         c_d(n)=0.0d0
         c_d_dum(n)=0.0d0

         c_lx(n)=0.0d0
         c_lx_dum(n)=0.0d0

         c_lz(n)=0.0d0
         c_lz_dum(n)=0.0d0
      enddo

C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(l,n)
      do n=1,num_p
         c_d_dum(n)  = -for_py(n)
         c_lx_dum(n) = -for_px(n)
         c_lz_dum(n) = -for_pz(n)

         if(irkk .eq. 1) then 
            c_d_total (n)=0.d0
            c_lx_total(n)=0.d0
            c_lz_total(n)=0.d0
            moment_x(n)=0.d0
            moment_y(n)=0.d0
            moment_z(n)=0.d0
         endif
         c_d(n)=c_d_dum(n)      !*2.d0/pi/r_c**2
         c_d_total(n)=c_d_total(n)+c_d(n)

         c_lx(n)=c_lx_dum(n)    !*2.d0/pi/r_c**2
         c_lx_total(n)=c_lx_total(n)+c_lx(n)
c
         c_lz(n)=c_lz_dum(n)    !*2.d0/pi/r_c**2
         c_lz_total(n)=c_lz_total(n)+c_lz(n)
c
         moment_x(n)=moment_x(n)+torq_x(n) !*2.d0/pi/r_c**3
         moment_y(n)=moment_y(n)+torq_y(n) !*2.d0/pi/r_c**3
         moment_z(n)=moment_z(n)+torq_z(n) !*2.d0/pi/r_c**3
      enddo
C$OMP END PARALLEL DO
c
      if(irkk .eq. 3) then
         c_d_cluster =0.d0
         c_lx_cluster=0.d0
         c_lz_cluster=0.d0
         do n=1,num_p
!            if(mod((itime-istart),oi_chan).eq.0) then
!               write(175,100) itime,n,c_d_total(n),
!     &              c_lx_total(n),c_lz_total(n),
!     &              moment_x(n),moment_y(n),moment_z(n)
!            endif
            c_d_cluster =c_d_cluster+c_d_total(n)
            c_lx_cluster=c_lx_cluster+c_lx_total(n)
            c_lz_cluster=c_lz_cluster+c_lz_total(n)
         enddo
         write(114,*) itime,c_lx_cluster
         write(115,*) itime,c_d_cluster
         write(116,*) itime,c_lz_cluster

         if(lunformatted) then
            write(223,rec=itime) itime, c_lx_total
            write(224,rec=itime) itime, c_d_total
            write(225,rec=itime) itime, c_lz_total
            
            write(123,rec=itime) itime, moment_x
            write(124,rec=itime) itime, moment_y
            write(125,rec=itime) itime, moment_z
         else
            
            write(223,310) itime,(c_lx_total(n),n=1,num_p)
            write(224,310) itime,(c_d_total(n) ,n=1,num_p)
            write(225,310) itime,(c_lz_total(n),n=1,num_p)
            
            write(123,310) itime,(moment_x(n),n=1,num_p)
            write(124,310) itime,(moment_y(n),n=1,num_p)
            write(125,310) itime,(moment_z(n),n=1,num_p)
         endif
      endif

 310  format (I7,<num_p>E15.7)
      
!>>>
      if(lcheck) then
         if(irkk .eq. 1) then
            forcing_x1 = forcing_x
            forcing_y1 = forcing_y
            forcing_z1 = forcing_z
         elseif(irkk .eq. 2) then
            forcing_x2 = forcing_x
            forcing_y2 = forcing_y
            forcing_z2 = forcing_z
         else
            forcing_x3 = forcing_x
            forcing_y3 = forcing_y
            forcing_z3 = forcing_z
         endif
      endif
!<<<      
         
c     Transforming the forcing function to fourier space

      call fft_r2f_2d_new(forcing_x)
      call fft_r2f_2d_new(forcing_y)
      call fft_r2f_2d_new(forcing_z)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      END SUBROUTINE ibm_ellip

      ! -------------------------------------------------- !
      subroutine int_f_fluid
      use common_m
      use ellip_common_m,only:x_0,x_1,y_0,y_1,z_0,z_1
      implicit none

      INTEGER n,i,j,k,id,jd
      
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(n,i,id,j,jd,k)
      do n=1,num_p
         int_fx(n) = 0.d0
         int_fy(n) = 0.d0
         int_fz(n) = 0.d0
         
         do i=x_0(n),x_1(n)
            if(i<1) then
               id = i+nx
            elseif(i>nx) then
               id = i-nx
            else
               id = i
            end if
          
            do j=y_0(n),y_1(n)
               if(j<1) then
                  jd = j+ny
               elseif(i>nx) then
                  jd = j-ny
               else
                  jd = j
               end if
               do k=z_0(n),z_1(n)

                  int_fx(n) = int_fx(n)+forcing_x(id,jd,k)*cell_v(k)
                  int_fy(n) = int_fy(n)+forcing_y(id,jd,k)*cell_v(k)
                  int_fz(n) = int_fz(n)+forcing_z(id,jd,k)*cell_v(k)

               end do 
            end do 
         end do 
      end do 
!$OMP END PARALLEL DO
      end subroutine int_f_fluid


