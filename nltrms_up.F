      SUBROUTINE nltrms_f
      use common_m

c     This subroutine computes the de-aliased nonlinear terms for NS eq.
c     in Fourier space
c     It uses Arakawa method, ie
c         odd times steps computes FFT[ u_i*ddx_i(u_j) ]
c         even time steps computes FFT[ ddx_i(u_i*u_j) ] with div(u) correction
c     The computed terms are returned in ru, rv and rw
c     After the subroutine finish the values of ul, vl and wl 
c     (that are in real space)should be available to be used 
c     by subroutine tt_nltrms and should not be used
c     until tt_nltrms has been called (in tt_comp_h)
c 
      INTEGER i,j,k
c
      if (FLTIMER) call begin_timer(56)
c
c     enlarge in Fourier space
      call enlarge_f_new(u,ul)
      call enlarge_f_new(v,vl)
      call enlarge_f_new(w,wl)
c     IFFT large variables
      call fftl_f2r_2d_new(ul)
      call fftl_f2r_2d_new(vl)
      call fftl_f2r_2d_new(wl)
c
      if (mod(itime,2).eq.1) then !odd time steps
c     on odd time-steps, use the convective form of the nonlinear term
        if (FLTIMER) call begin_timer(51)
c
c     u*ddx(u)+v*ddy(u)+w*ddz(u) --> ru
        call ddx_f_new(u,ru)
        call enlarge_f_new(ru,ddxul)
        call fftl_f2r_2d_new(ddxul)
        call ddy_f_new(u,ru)
        call enlarge_f_new(ru,ddyul)
        call fftl_f2r_2d_new(ddyul)
        call ddz(u,ru)!ddz computes in real and Fourier spaces
        call enlarge_f_new(ru,ddzul)
        call fftl_f2r_2d_new(ddzul)
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
        do k=1,nz
          do j=1,nyl
            do i=1,nxl
              ddxul(i,j,k)=ul(i,j,k)*ddxul(i,j,k)
     &                    +vl(i,j,k)*ddyul(i,j,k)
     &                    +wl(i,j,k)*ddzul(i,j,k)
            enddo
          enddo
        enddo
C$OMP END PARALLEL DO
        call fftl_r2f_2d_new(ddxul)
        call reduce_f_new(ddxul,ru)
c
c     u*ddx(v)+v*ddy(v)+w*ddz(v) --> rv
        call ddx_f_new(v,rv)
        call enlarge_f_new(rv,ddxul)
        call fftl_f2r_2d_new(ddxul)
        call ddy_f_new(v,rv)
        call enlarge_f_new(rv,ddyul)
        call fftl_f2r_2d_new(ddyul)
        call ddz(v,rv)!ddz computes in real and Fourier spaces
        call enlarge_f_new(rv,ddzul)
        call fftl_f2r_2d_new(ddzul)
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
        do k=1,nz
          do j=1,nyl
            do i=1,nxl
              ddxul(i,j,k)=ul(i,j,k)*ddxul(i,j,k)
     &                    +vl(i,j,k)*ddyul(i,j,k)
     &                    +wl(i,j,k)*ddzul(i,j,k)
            enddo
          enddo
        enddo
C$OMP END PARALLEL DO
        call fftl_r2f_2d_new(ddxul)
        call reduce_f_new(ddxul,rv)
c
c     u*ddx(w)+v*ddy(w)+w*ddz(w) --> rw
        call ddx_f_new(w,rw)
        call enlarge_f_new(rw,ddxul)
        call fftl_f2r_2d_new(ddxul)
        call ddy_f_new(w,rw)
        call enlarge_f_new(rw,ddyul)
        call fftl_f2r_2d_new(ddyul)
        call ddz(w,rw)!ddz computes in real and Fourier spaces
        call enlarge_f_new(rw,ddzul)
        call fftl_f2r_2d_new(ddzul)
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
        do k=1,nz
          do j=1,nyl
            do i=1,nxl
              ddxul(i,j,k)=ul(i,j,k)*ddxul(i,j,k)
     &                    +vl(i,j,k)*ddyul(i,j,k)
     &                    +wl(i,j,k)*ddzul(i,j,k)
            enddo
          enddo
        enddo
C$OMP END PARALLEL DO
        call fftl_r2f_2d_new(ddxul)
        call reduce_f_new(ddxul,rw)
c
        if (FLTIMER) call end_timer(51)
c
      else !even time steps
c     on even-time steps, use the divergence form of the nonlinear term
c
        if (FLTIMER) call begin_timer(52)
c
c     u*u, u*v, u*w
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
        do k=1,nz
          do j=1,nyl
            do i=1,nxl
              ddxul(i,j,k)=ul(i,j,k)*ul(i,j,k)
              ddyul(i,j,k)=vl(i,j,k)*ul(i,j,k)
              ddzul(i,j,k)=wl(i,j,k)*ul(i,j,k)
            enddo
          enddo
        enddo
C$OMP END PARALLEL DO
c
c     ddx(u*u)+ddy(u*v)+ddz(u*w) --> ru
        call fftl_r2f_2d_new(ddxul)
        call reduce_f_new(ddxul,rw)
        call ddx_f_new(rw,ru)
        call fftl_r2f_2d_new(ddyul)
        call reduce_f_new(ddyul,rp)
        call ddy_f_new(rp,rv)
        call fftl_r2f_2d_new(ddzul)
        call reduce_f_new(ddzul,rtt)
        call ddz(rtt,rw)!ddz computes in real and Fourier spaces
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
        do k=1,nz
          do j=1,ny
            do i=1,nx
              ru(i,j,k)=ru(i,j,k)+rv(i,j,k)+rw(i,j,k)
            enddo
          enddo
        enddo
C$OMP END PARALLEL DO
c
c     v*v, v*w, w*w
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
        do k=1,nz
          do j=1,nyl
            do i=1,nxl
              ddxul(i,j,k)=vl(i,j,k)*vl(i,j,k)
              ddyul(i,j,k)=vl(i,j,k)*wl(i,j,k)
              ddzul(i,j,k)=wl(i,j,k)*wl(i,j,k)
            enddo
          enddo
        enddo
C$OMP END PARALLEL DO
c
c     ddx(v*u)+ddy(v*v)+ddz(v*w) --> rv
        call ddx_f_new(rp,rv)
        call fftl_r2f_2d_new(ddxul)
        call reduce_f_new(ddxul,rp)
        call ddy_f_new(rp,rw)
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
        do k=1,nz
          do j=1,ny
            do i=1,nx
              rv(i,j,k)=rv(i,j,k)+rw(i,j,k)
            enddo
          enddo
        enddo
C$OMP END PARALLEL DO
        call fftl_r2f_2d_new(ddyul)
        call reduce_f_new(ddyul,rp)   
        call ddz(rp,rw)!ddz computes in real and Fourier spaces
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
        do k=1,nz
          do j=1,ny
            do i=1,nx
              rv(i,j,k)=rv(i,j,k)+rw(i,j,k)
            enddo
          enddo
        enddo
C$OMP END PARALLEL DO
c
c     ddx(w*u)+ddy(w*v)+ddz(w*w) --> rw
        call ddx_f_new(rtt,rw)
        call ddy_f_new(rp,rtt)
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
        do k=1,nz
          do j=1,ny
            do i=1,nx
              rw(i,j,k)=rw(i,j,k)+rtt(i,j,k)
            enddo
          enddo
        enddo
C$OMP END PARALLEL DO
        call fftl_r2f_2d_new(ddzul)
        call reduce_f_new(ddzul,rp)
        call ddz(rp,rtt)!ddz computes in real and Fourier spaces
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
        do k=1,nz
          do j=1,ny
            do i=1,nx
              rw(i,j,k)=rw(i,j,k)+rtt(i,j,k)
            enddo
          enddo
        enddo
C$OMP END PARALLEL DO
c
        if (FLTIMER) call end_timer(52)
c
      endif
c
      if (FLTIMER) call end_timer(56)
      return
      END
c
c
      SUBROUTINE tt_nltrms_f
      use common_m

c     This subroutine computes the de-aliased nonlinear terms
c     for thermal equation in Fourier space, ie
c         FFT[ ddx_i(up_i*tt) ] (sum implyed en repeated i)
c     where up is fluid or particles velocity 
c     (see FLSTLE below)
c     If fluid velocity is used, the values of ul, vl and wl (in real space) 
c     computed in subroutine nltrs_f are used. Otherwise, these
c     variables are overwritten by the enlarged particle velocities 
c     The computed non-linear term is returned in rtt in Fourier space    
c
c     We might want to separate this in tt*div(up)+up_k*dtt/dx_k
c     if it behaves bad, or try some version of Arakawa's method
c
c NOTE: these large variables are defined as a common block
c      REAL*8 ul(nxl0,nyl0,nz),vl(nxl0,nyl0,nz),wl(nxl0,nyl0,nz)
c      REAL*8 ddxul(nxl0,nyl0,nz),ddyul(nxl0,nyl0,nz),ddzul(nxl0,nyl0,nz)
c END NOTE
c      REAL*8 igmod
      INTEGER i,j,k
c
c      if (FLTIMER) call begin_timer()
c
c     enlarge tt
      call enlarge_f_new(tt,ddxul)
      call fftl_f2r_2d_new(ddxul)
c
c     ddx(u*tt)+ddy(v*tt)+ddz(w*tt) --> rtt
      if(FLSTLE)then
C recall ws is magnitude (>0)
c         igmod=sqrt(gy(1,1)**2+gz(1,1)**2)
c         igmod=1.d0/igmod
C NOTE by mic: may want to compute settling velocity components
C              in initial.F to speed up
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
         do k=1,nz
            do j=1,nyl
               do i=1,nxl
                  ul(i,j,k)=ul(i,j,k)*ddxul(i,j,k)
                  vl(i,j,k)=(vl(i,j,k)
     &             +ws*gy(1,k)/sqrt(gy(1,k)**2+gz(1,k)**2))
     &                 *ddxul(i,j,k)
                  wl(i,j,k)=(wl(i,j,k)
     &             +ws*gz(1,k)/sqrt(gy(1,k)**2+gz(1,k)**2))
     &                 *ddxul(i,j,k)
               enddo
            enddo
         enddo
C$OMP END PARALLEL DO
      else !advection by fluid velocity
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
         do k=1,nz
            do j=1,nyl
               do i=1,nxl
                  ul(i,j,k)=ul(i,j,k)*ddxul(i,j,k)
                  vl(i,j,k)=vl(i,j,k)*ddxul(i,j,k)
                  wl(i,j,k)=wl(i,j,k)*ddxul(i,j,k)
               enddo
            enddo
         enddo
C$OMP END PARALLEL DO
      endif!FLSTLE
c
c     reduce variables
      call fftl_r2f_2d_new(ul)
      call reduce_f_new(ul,rp)
      call ddx_f_new(rp,ru)
      call fftl_r2f_2d_new(vl)
      call reduce_f_new(vl,rp)
      call ddy_f_new(rp,rv)
      call fftl_r2f_2d_new(wl)
      call reduce_f_new(wl,rp)
      call ddz(rp,rw)!ddz computes in real and Fourier spaces
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
      do k=1,nz
         do j=1,ny
            do i=1,nx
               rtt(i,j,k)=ru(i,j,k)+rv(i,j,k)+rw(i,j,k)
            enddo
         enddo
      enddo
C$OMP END PARALLEL DO
c
c      if (FLTIMER) call end_timer()
      return
      END
