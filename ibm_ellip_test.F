      SUBROUTINE ibm_ellip(cfnl2,cfdif)
      use ellipsoid_m, only: for_px,for_py,for_pz,
     &     torq_x,torq_y,torq_z,lcheck,single_point_test
      use p_dyn_m
      use common_m
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

!>>>
      real*8 f_in(3,3,6),f_out(3,3,6)
!<<<
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
               u_tilde(i,j,k)=
     &                        (cfnl2*h1(i,j,k)+2.d0*cfdif*ru(i,j,k))/dt
               v_tilde(i,j,k)=
     &                        (cfnl2*h2(i,j,k)+2.d0*cfdif*rv(i,j,k))/dt
               w_tilde(i,j,k)=
     &                        (cfnl2*h3(i,j,k)+2.d0*cfdif*rw(i,j,k))/dt

            enddo
         enddo
      enddo
C$OMP END PARALLEL DO
c This section follows the flow solver of M. Uhlmann (2005) to
c define the forcing function.
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

            f_in = u_tilde(i1:i2,j1:j2,k1:k2)
            call single_point_test(f_in,f_out)
            forcing_x(i1:i2,j1:j2,k1:k2) = f_out(1:3,1:3,1:6)

c$$$            write(*,*) i1,i2,j1,j2,k1,k2,'---'
c$$$            write(*,*) 'f_in'
c$$$            write(*,*) f_in
c$$$            write(*,*) 'f_out'
c$$$            write(*,*) f_out
c$$$            
            f_in = v_tilde(i1:i2,j1:j2,k1:k2)
            call single_point_test(f_in,f_out)
            forcing_y(i1:i2,j1:j2,k1:k2) = f_out(1:3,1:3,1:6)

            f_in = w_tilde(i1:i2,j1:j2,k1:k2)
            call single_point_test(f_in,f_out)
            forcing_z(i1:i2,j1:j2,k1:k2) = f_out(1:3,1:3,1:6)

         enddo
      enddo
C$OMP END PARALLEL DO

      
c$$$      if(lcheck) then
c$$$         call int_f_fluid
c$$$         
c$$$         write(401,*) itime,int_fx,for_px
c$$$         write(402,*) itime,int_fy,for_py
c$$$         write(403,*) itime,int_fz,for_pz
c$$$      endif


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
      
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(n)
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
C$OMP END PARALLEL DO
      end subroutine int_f_fluid


