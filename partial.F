      SUBROUTINE ddx(val,dval)
      use common_m

      REAL*8 val(nx0,ny0,nz),dval(nx0,ny0,nz)
c
      INTEGER i,j,k
      REAL*8     rtmp(tnxhp0)
      COMPLEX*16 ctmp( nxhp0)
c
#ifdef SCSL
c for SCSL
      REAL*8 fft_work(tnxhp0)
#endif
c
#ifdef ESSL
c for ESSL
c fft_work must be 8-byte real array
      REAL*8 fft_work(nfft_aux)
#endif
c
c      if (FLTIMER) call begin_timer(80)
c
#ifdef FFTW3
c for FFTW3
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,rtmp,ctmp)
      do k=1,nz
         do j = 1,ny
            call dcopy(nx,val(1,j,k),1,rtmp,1)
c scale input data for normalization
            call dscal(nx,1.d0/dble(nx),rtmp,1)
            call dfftw_execute_dft_r2c(plan_1_r2c_x,rtmp,ctmp)
            do i=1,nxh
               ctmp(i)=dcmplx((0.d0,1.d0)*wavex(i))*ctmp(i)
            enddo
            ctmp(nxhp)=dcmplx(0.d0)
            call dfftw_execute_dft_c2r(plan_1_c2r_x,ctmp,rtmp)
            call dcopy(nx,rtmp,1,dval(1,j,k),1)
         enddo
      enddo
C$OMP END PARALLEL DO
#endif
c
#ifdef SCSL
c for SCSL
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,rtmp,ctmp,fft_work)
      do k=1,nz
         do j = 1,ny
            call dcopy(nx,val(1,j,k),1,rtmp,1)
            call dzfftm(-1,nx,1,1.d0/dble(nx),rtmp,tnxhp0,ctmp,nxhp0,
     &           rx_tb_f,fft_work,1)
            do i=1,nxh
               ctmp(i)=dcmplx((0.d0,1.d0)*wavex(i))*ctmp(i)
            enddo
            ctmp(nxhp)=dcmplx(0.d0)
            call zdfftm(1,nx,1,1.d0,ctmp,nxhp0,rtmp,tnxhp0,
     &           rx_tb_b,fft_work,1)
            call dcopy(nx,rtmp,1,dval(1,j,k),1)
         enddo
      enddo
C$OMP END PARALLEL DO
#endif
c
#ifdef ESSL
c for ESSL
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,rtmp,ctmp,fft_work)
      do k=1,nz
         do j = 1,ny
            call dcopy(nx,val(1,j,k),1,rtmp,1)
            call drcft(0,rtmp,tnxhp0,ctmp,nxhp0,nx,1,1,1.d0/dble(nx),
     &           rx_tb_f,nfft_aux,fft_work,nfft_aux)
            do i=1,nxh
               ctmp(i)=dcmplx(((0.d0,1.d0)*wavex(i)))*ctmp(i)
            enddo
            ctmp(nxhp)=dcmplx(0.d0)
            call dcrft(0,ctmp,nxhp0,rtmp,tnxhp0,nx,1,-1,1.d0,
     &           rx_tb_b,nfft_aux,fft_work,nfft_aux)
            call dcopy(nx,rtmp,1,dval(1,j,k),1)
         enddo
      enddo
C$OMP END PARALLEL DO
#endif
c
c      if (FLTIMER) call end_timer(80)
      return
      END
c
c
      SUBROUTINE ddy(val,dval)
      use common_m

      REAL*8 val(nx0,ny0,nz),dval(nx0,ny0,nz)
c
      INTEGER i,j,k
      REAL*8     rtmp(tnyhp0)
      COMPLEX*16 ctmp( nyhp0)
c
#ifdef SCSL
c for SCSL
      REAL*8 fft_work(tnyhp0)
#endif
c
#ifdef ESSL
c for ESSL
c fft_work must be 8-byte real array
      REAL*8 fft_work(nfft_aux)
#endif
c
c      if (FLTIMER) call begin_timer(82)
c
#ifdef FFTW3
c for FFTW3
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,rtmp,ctmp)
      do k=1,nz
         do i=1,nx
            call dcopy(ny,val(i,1,k),nx0,rtmp,1)
c scale input data for normalization
            call dscal(ny,1.d0/dble(ny),rtmp,1)
            call dfftw_execute_dft_r2c(plan_1_r2c_y,rtmp,ctmp)
            do j=1,nyh
               ctmp(j)=dcmplx((0.d0,1.d0)*wavey(j))*ctmp(j)
            enddo
            ctmp(nyhp)=dcmplx(0.d0)
            call dfftw_execute_dft_c2r(plan_1_c2r_y,ctmp,rtmp)
            call dcopy(ny,rtmp,1,dval(i,1,k),nx0)
         enddo
      enddo
C$OMP END PARALLEL DO
#endif
c
#ifdef SCSL
c for SCSL
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,rtmp,ctmp,fft_work)
      do k=1,nz
         do i=1,nx
            call dcopy(ny,val(i,1,k),nx0,rtmp,1)
            call dzfftm(-1,ny,1,1.d0/dble(ny),rtmp,tnyhp0,ctmp,nyhp0,
     &           ry_tb_f,fft_work,1)
            do j=1,nyh
               ctmp(j)=dcmplx(((0.d0,1.d0)*wavey(j)))*ctmp(j)
            enddo
            ctmp(nyhp)=dcmplx(0.d0)
            call zdfftm(1,ny,1,1.d0,ctmp,nyhp0,rtmp,tnyhp0,
     &           ry_tb_b,fft_work,1)
            call dcopy(ny,rtmp,1,dval(i,1,k),nx0)
         enddo
      enddo
C$OMP END PARALLEL DO
#endif
c
#ifdef ESSL
c for ESSL
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,rtmp,ctmp,fft_work)
      do k=1,nz
         do i=1,nx
            call dcopy(ny,val(i,1,k),nx0,rtmp,1)
            call drcft(0,rtmp,tnyhp0,ctmp,nyhp0,ny,1,1,1.d0/dble(ny),
     &           ry_tb_f,nfft_aux,fft_work,nfft_aux)
            do j=1,nyh
               ctmp(j)=dcmplx(((0.d0,1.d0)*wavey(j)))*ctmp(j)
            enddo
            ctmp(nyhp)=dcmplx(0.d0)
            call dcrft(0,ctmp,nyhp0,rtmp,tnyhp0,ny,1,-1,1.d0,
     &           ry_tb_b,nfft_aux,fft_work,nfft_aux)
            call dcopy(ny,rtmp,1,dval(i,1,k),nx0)
         enddo
      enddo
C$OMP END PARALLEL DO
#endif
c
c      if (FLTIMER) call end_timer(82)
      return
      END
c
c
      SUBROUTINE ddz(val,dval)!ddz_f_new(val,dval)
      use common_m

      REAL*8 val(nx0,ny0,nz),dval(nx0,ny0,nz)
c
c     This subroutine computes ddz in real and Fourier spaces
c     There is not difference in its implementation
c
      INTEGER i,j,jcomp,k
c
c      if (FLTIMER) call begin_timer()
c
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,jcomp,k)
      do j = 1,nyh
c     First part. This is done this way to reuse auxt1/auxt2(nz0,nx0,nyh)
c     traspose data
        do k = 1,nz
          call dcopy(nx,val(1,j,k),1,auxt1(k,1,j),nz0)
        enddo
c     cumpute
       call dgemm('N','N',nz,nx,nz,1.d0,chbd1,nz0,auxt1(1,1,j),nz0,
     &                                       0.d0,auxt2(1,1,j),nz0)
c     transpose back data
        do i=1,nx
          call dcopy(nz,auxt2(1,i,j),1,dval(i,j,1),nx0y0)
        enddo
c     Second part.
        jcomp=nyh+j
c     traspose data
        do k = 1,nz
          call dcopy(nx,val(1,jcomp,k),1,auxt1(k,1,j),nz0)
        enddo
c     compute
       call dgemm('N','N',nz,nx,nz,1.d0,chbd1,nz0,auxt1(1,1,j),nz0,
     &                                       0.d0,auxt2(1,1,j),nz0)
c     transpose back data
        do i=1,nx
          call dcopy(nz,auxt2(1,i,j),1,dval(i,jcomp,1),nx0y0)
        enddo
      enddo
C$OMP END PARALLEL DO
c
c      if (FLTIMER) call end_timer()
      return
      END
c
c
      SUBROUTINE ddx_f_new(var,dvar)
      use common_m

      REAL*8 var(nx0,ny0,nz)
      REAL*8 dvar(nx0,ny0,nz)
c     The input var is preserved
      INTEGER i,j,k,jimag
c
c      if (FLTIMER) call begin_timer()
c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,jimag,k)
        do k=1,nz
         do j=1,ny
           call dcopy(nx,var(1,j,k),1,dvar(1,j,k),1)
         enddo
         do i=1,nx
           call dscal(nyh, wavex(i),dvar(i,1   ,k),nx0)
           call dscal(nyh,-wavex(i),dvar(i,nyhp,k),nx0)
         enddo
         do j=1,nyh
           jimag=nyh+j
           call dswap(nx,dvar(1,j,k),1,dvar(1,jimag,k),1)
         enddo
        enddo
C$OMP END PARALLEL DO
c
c      if (FLTIMER) call end_timer()
      return
      END
c
c
      SUBROUTINE ddy_f_new(var,dvar)
      use common_m

      REAL*8 var(nx0,ny0,nz)
      REAL*8 dvar(nx0,ny0,nz)
c     The input var is preserved
      INTEGER j,k,jimag
c
c      if (FLTIMER) call begin_timer()
c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(j,jimag,k)
        do k=1,nz
         do j=1,nyh
           jimag=nyh+j
           call dcopy(nx,var(1,j    ,k),1,dvar(1,j    ,k),1)
           call dcopy(nx,var(1,jimag,k),1,dvar(1,jimag,k),1)
c
           call dscal(nx, wavey(j),dvar(1,j    ,k),1)
           call dscal(nx,-wavey(j),dvar(1,jimag,k),1)
c
           call dswap(nx,dvar(1,j,k),1,dvar(1,jimag,k),1)
         enddo
        enddo
C$OMP END PARALLEL DO
c
c      if (FLTIMER) call end_timer()
      return
      END
c
c
      SUBROUTINE ddt(irk,val0,val1,val2,val3,dval)
      use common_m

      INTEGER irk
      REAL*8 val0(nx0,ny0,nz),val1(nx0,ny0,nz),val2(nx0,ny0,nz)
     &       ,val3(nx0,ny0,nz),dval(nx0,ny0,nz)
c
      REAL*8 hh0,hh1,hh2,hh3,a0,a1,a2,a3
      INTEGER i,j,k
c
      if (FLTIMER) call begin_timer(86)
c initialize values 
      if(irk.eq.1)then
        hh0=0.d0
      elseif(irk.eq.2)then
        hh0=dt/3.d0
      elseif(irk.eq.3)then
        hh0=3.d0*dt/4.d0
      endif
      hh1=dt+hh0
      hh2=dt+hh1
      hh3=dt+hh2
c for restart with different dt
      if (istart.gt.0.and.itime.eq.istart+1)then
        hh1=dt_old+hh0
        hh2=dt_old+hh1
        hh3=dt_old+hh2
      elseif(istart.gt.0.and.itime.eq.istart+2)then
        hh2=dt_old+hh1
        hh3=dt_old+hh2
      elseif(istart.gt.0.and.itime.eq.istart+3)then
        hh3=dt_old+hh2
      endif     
c 3rd order
      a0=-(hh2*hh3+hh1*(hh2+hh3))/(hh0-hh1)/(hh0-hh2)/(hh0-hh3)
      a1=-(hh2*hh3+hh0*(hh2+hh3))/(hh1-hh0)/(hh1-hh2)/(hh1-hh3)
      a2=-(hh1*hh3+hh0*(hh1+hh3))/(hh2-hh0)/(hh2-hh1)/(hh2-hh3)
      a3=-(hh1*hh2+hh0*(hh1+hh2))/(hh3-hh0)/(hh3-hh1)/(hh3-hh2)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
       do k=1,nz
        do j=1,ny
         do i=1,nx
c 3rd order 
          dval(i,j,k)=a0*val0(i,j,k)+a1*val1(i,j,k)
     &               +a2*val2(i,j,k)+a3*val3(i,j,k)
         enddo
        enddo
       enddo
C$OMP END PARALLEL DO
      if (FLTIMER) call end_timer(86)
      return
      END
