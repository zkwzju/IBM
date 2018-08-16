      subroutine initfft
      use common_m
c
#ifdef FFTW3	
      REAL*8 raux2(tnyhp0)
      COMPLEX*16 caux(nxhp0),caux2(nyhp0)
c     
c for fft_r2f_2d
      call dfftw_plan_many_dft_r2c(plan_nx_r2c_y,1,ny,nx,
     &     u,ny0,nx0,1,
     &     tmpf,nyhp0,1,nyhp0,
     &     FFTW_UNALIGNED+FFTW_MEASURE)
      call dfftw_plan_many_dft(plan_nyh_c2cf_x,1,nx,nyh,
     &     tmpf,nx0,nyhp0,1,
     &     tmpf,nx0,nyhp0,1,
     &     FFTW_FORWARD,FFTW_UNALIGNED+FFTW_MEASURE)
c for fft_f2r_2d
      call dfftw_plan_many_dft(plan_nyh_c2cb_x,1,nx,nyh,
     &     tmpf,nx0,nyhp0,1,
     &     tmpf,nx0,nyhp0,1,
     &     FFTW_BACKWARD,FFTW_UNALIGNED+FFTW_MEASURE)
      call dfftw_plan_many_dft_c2r(plan_nx_c2r_y,1,ny,nx,
     &     tmpf,ny0,1,nyhp0,    !nyhp0,1,nyhp0,!
     &     u,ny0,nx0,1,
     &     FFTW_BACKWARD,FFTW_UNALIGNED+FFTW_MEASURE)
c for fftl_r2f_2d
      call dfftw_plan_many_dft_r2c(plan_nxl_r2c_y,1,nyl,nxl,
     &     ul,nyl0,nxl0,1,
     &     tmplf,nylhp0,1,nylhp0,
     &     FFTW_UNALIGNED+FFTW_MEASURE)
      call dfftw_plan_many_dft(plan_nylh_c2cf_x,1,nxl,nylh,
     &     tmplf,nxl0,nylhp0,1,
     &     tmplf,nxl0,nylhp0,1,
     &     FFTW_FORWARD,FFTW_UNALIGNED+FFTW_MEASURE)
c for fftl_f2r_2d
      call dfftw_plan_many_dft(plan_nylh_c2cb_x,1,nxl,nylh,
     &     tmplf,nxl0,nylhp0,1,
     &     tmplf,nxl0,nylhp0,1,
     &     FFTW_BACKWARD,FFTW_UNALIGNED+FFTW_MEASURE)
      call dfftw_plan_many_dft_c2r(plan_nxl_c2r_y,1,nyl,nxl,
     &     tmplf,nyl0,1,nylhp0, !nylhp0,1,nylhp0,!
     &     ul,nyl0,nxl0,1,
     &     FFTW_BACKWARD,FFTW_UNALIGNED+FFTW_MEASURE)
c for ddx
      call dfftw_plan_dft_r2c(plan_1_r2c_x,1,nx,
     &     u,
     &     caux,
     &     FFTW_UNALIGNED+FFTW_MEASURE)
      call dfftw_plan_dft_c2r(plan_1_c2r_x,1,nx,
     &     caux,
     &     u,
     &     FFTW_UNALIGNED+FFTW_MEASURE)
c for ddy
      call dfftw_plan_dft_r2c(plan_1_r2c_y,1,ny,
     &     raux2,
     &     caux2,
     &     FFTW_UNALIGNED+FFTW_MEASURE)
      call dfftw_plan_dft_c2r(plan_1_c2r_y,1,ny,
     &     caux2,
     &     raux2,
     &     FFTW_UNALIGNED+FFTW_MEASURE)
#endif
c
#ifdef SCSL
c for SCSL library
      call dzfftm(0,ny,1,1.d0,ru,1,ru,1,ry_tb_f,ru,1)
      call zzfftm(0,nx,1,1.d0,ru,1,ru,1,cx_tb,ru,1)
      call zdfftm(0,ny,1,1.d0,ru,1,ru,1,ry_tb_b,ru,1)
c     
      call dzfftm(0,nyl,1,1.d0,ru,1,ru,1,ry_ltb_f,ru,1)
      call zzfftm(0,nxl,1,1.d0,ru,1,ru,1,cx_ltb,ru,1)
      call zdfftm(0,nyl,1,1.d0,ru,1,ru,1,ry_ltb_b,ru,1)
c     
      call dzfftm(0,nx,1,1.d0,ru,1,ru,1,rx_tb_f,ru,1)
      call zdfftm(0,nx,1,1.d0,ru,1,ru,1,rx_tb_b,ru,1)
#endif
c
#ifdef ESSL
      REAL*8 fft_r(1)
      COMPLEX*16 fft_c(1)
c  fft_work must be 8-byte real array
      REAL*8 fft_work(nfft_aux)
c
      call drcft(1,fft_r,tnyhp0,fft_c,nyhp0,ny,1,1,1.d0/dble(ny),
     &     ry_tb_f,nfft_aux,fft_work,nfft_aux)
      call dcft(1,fft_c,nyhp0,1,fft_c,1,nx0,nx,1,1,1.d0/dble(nx),
     &     cx_tb_f,nfft_aux,fft_work,nfft_aux)
      call dcft(1,fft_c,1,nx0,fft_c,nyhp0,1,nx,1,-1,1.d0,
     &     cx_tb_b,nfft_aux,fft_work,nfft_aux)
      call dcrft(1,fft_c,nyhp0,fft_r,tnyhp0,ny,1,-1,1.d0,
     &     ry_tb_b,nfft_aux,fft_work,nfft_aux)
      call drcft(1,fft_r,tnxhp0,fft_c,nxhp0,nx,1,1,1.d0/dble(nx),
     &     rx_tb_f,nfft_aux,fft_work,nfft_aux)
      call dcrft(1,fft_c,nxhp0,fft_r,tnxhp0,nx,1,-1,1.d0,
     &     rx_tb_b,nfft_aux,fft_work,nfft_aux)
      call drcft(1,fft_r,tnylhp0,fft_c,nylhp0,nyl,1,1,1.d0/dble(nyl),
     &     ry_ltb_f,nfft_aux,fft_work,nfft_aux)
      call dcft(1,fft_c,nylhp0,1,fft_c,1,nxl0,nxl,1,1,
     &     1.d0/dble(nxl),cx_ltb_f,nfft_aux,fft_work,nfft_aux)
      call dcft(1,fft_c,1,nxl0,fft_c,nylhp0,1,nxl,1,-1,1.d0,
     &     cx_ltb_b,nfft_aux,fft_work,nfft_aux)
      call dcrft(1,fft_c,nylhp0,fft_r,tnylhp0,nyl,1,-1,1.d0,
     &     ry_ltb_b,nfft_aux,fft_work,nfft_aux)
#endif
c
      return
      end
c
c
      subroutine cleanfft
      use common_m
c
#ifdef FFTW3
      call dfftw_destroy_plan(plan_nx_r2c_y)
      call dfftw_destroy_plan(plan_nyh_c2cf_x)
c
      call dfftw_destroy_plan(plan_nyh_c2cb_x)
      call dfftw_destroy_plan(plan_nx_c2r_y)
c
      call dfftw_destroy_plan(plan_nxl_r2c_y)
      call dfftw_destroy_plan(plan_nylh_c2cf_x)
c
      call dfftw_destroy_plan(plan_nylh_c2cb_x)
      call dfftw_destroy_plan(plan_nxl_c2r_y)
c
      call dfftw_destroy_plan(plan_1_r2c_x)
      call dfftw_destroy_plan(plan_1_c2r_x)
c
      call dfftw_destroy_plan(plan_1_r2c_y)
      call dfftw_destroy_plan(plan_1_c2r_y)
c
#endif
      return
      end
c
c
      subroutine fft_r2f_2d_new(ur)
      use common_m
      REAL*8 ur(nx0,ny0,nz)
c
c     2-d real to fourier transform in x and y directions
c     On input:  ur contains data on 1:nx,1:ny,1:nz
c     On output: ur contains data on (1:nx),1:nyh,1:nz   (real component)
c     On output: ur contains data on (1:nx),nyhp:ny,1:nz (imag component)
c
      INTEGER i,j,jimag,k
c      COMPLEX*16 tmp(nyhp0,nx0,nz)!this is a common block now ->tmpf
c
#ifdef SCSL
c for SCSL
      REAL*8 tp(tnyhp0)
      COMPLEX*16 tf(nx0)
      COMPLEX*16 tx(nx0)
      REAL*8 ry_wk(tnyhp0),cx_wk(tnx)
#endif
c
#ifdef ESSL
c for ESSL
c fft_work must be 8-byte real array
      REAL*8 tp(tnyhp0)
      COMPLEX*16 tf(nx0)
      REAL*8 fft_work(nfft_aux)
#endif
c
      if (FLTIMER) call begin_timer(70)
c
#ifdef FFTW3
c FFTW3
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,jimag,k)
      do k=1,nz
c scale input data for normalization
         call dscal(nx0y,1.d0/dble(nx)/dble(ny),ur(1,1,k),1)
c compute fft in y-direction
         call dfftw_execute_dft_r2c(plan_nx_r2c_y
     &        ,ur(1,1,k),tmpf(1,1,k))
c compute fft in x-direction
         call dfftw_execute_dft(plan_nyh_c2cf_x
     &        ,tmpf(1,1,k),tmpf(1,1,k))         
c copy to output variable
         do j=1,nyh 
            jimag=nyh+j
            do i=1,nx
               ur(i,j    ,k) = dreal(tmpf(j,i,k))
               ur(i,jimag,k) = dimag(tmpf(j,i,k))
            enddo
         enddo
      enddo
C$OMP END PARALLEL DO
#endif
c
#ifdef SCSL
c for SCSL
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP&            PRIVATE(i,j,jimag,k,tp,tf,tx,cx_wk,ry_wk)
      do k=1,nz
c compute fft in y-direction
         do i=1,nx
            call dcopy(ny,ur(i,1,k),nx0,tp,1)
            call dzfftm(-1,ny,1,1.d0/dble(ny),tp,tnyhp0,tmpf(1,i,k),
     &           nyhp0,ry_tb_f,ry_wk,1)
         enddo
c compute fft in x-direction
         do j=1,nyh 
            call zcopy(nx,tmpf(j,1,k),nyhp0,tx,1)
            call zzfftm(-1,nx,1,1.d0/dble(nx),tx,nx0,tf,nx0,
     &           cx_tb,cx_wk,1)
            jimag=nyh+j
            do i=1,nx
               ur(i,j    ,k) = dreal(tf(i))
               ur(i,jimag,k) = dimag(tf(i))
            enddo
         enddo
      enddo
C$OMP END PARALLEL DO
#endif
c
#ifdef ESSL
c for ESSL
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,jimag,k,tp,tf,fft_work)
      do k=1,nz
c compute fft in y-direction
         do i=1,nx
            call dcopy(ny,ur(i,1,k),nx0,tp,1)
            call drcft(0,tp,tnyhp0,tmpf(1,i,k),nyhp0,ny,1,1,
     &           1.d0/dble(ny),ry_tb_f,nfft_aux,fft_work,nfft_aux)
         enddo
c compute fft in x-direction
         do j=1,nyh 
            call dcft(0,tmpf(j,1,k),nyhp0,1,tf,1,nx0,nx,1,1,
     &           1.d0/dble(nx),cx_tb_f,nfft_aux,fft_work,nfft_aux)
            jimag=nyh+j
            do i=1,nx
               ur(i,j    ,k) = dreal(tf(i))
               ur(i,jimag,k) = dimag(tf(i))
            enddo
         enddo
      enddo
C$OMP END PARALLEL DO
#endif
c
      if (FLTIMER) call end_timer(70)
      return
      end
c
c
      subroutine fft_f2r_2d_new(ur)
      use common_m
      REAL*8 ur(nx0,ny0,nz)
c
c     2-d fourier to real transform in x and y directions
c     On input:  ur contains data on 1:nx,   1:nyh,1:nz (real component)
c     On input:  ur contains data on 1:nx,nyhp:ny ,1:nz (imag component)
c     On output: ur contains data on 1:nx,1:ny,1:nz
c
      INTEGER i,j,k,jimag
c      COMPLEX*16 tmp(nyhp0,nx0,nz)!this is a common block now ->tmpf
c
#ifdef SCSL
c for SCSL
      REAL*8 tp(tnyhp0)
      COMPLEX*16 tf(nx0)
      COMPLEX*16 tx(nx0)
      REAL*8 ry_wk(tnyhp0),cx_wk(tnx)
#endif
c
#ifdef ESSL
c for ESSL
c fft_work must be 8-byte real array
      REAL*8 tp(tnyhp0)
      COMPLEX*16 tf(nx0)
      REAL*8 fft_work(nfft_aux)
#endif
c
      if (FLTIMER) call begin_timer(75)
c
#ifdef FFTW3
c FFTW3
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,jimag,k)
      do k=1,nz
         do j=1,nyh
            jimag=nyh+j
            do i=1,nx
               tmpf(j,i,k) = dcmplx(ur(i,j,k),ur(i,jimag,k))
            enddo
         enddo
c cumpute fft in x-direction
         call dfftw_execute_dft(plan_nyh_c2cb_x,
     &        tmpf(1,1,k),tmpf(1,1,k))
c cumpute fft in y-direction
         do i=1,nx
            tmpf(nyhp,i,k)=dcmplx(0.d0)
         enddo
         call dfftw_execute_dft_c2r(plan_nx_c2r_y
     &        ,tmpf(1,1,k),ur(1,1,k))         
c         do i=1,nx0
c            ur(i,ny0,k)=0.d0
c         enddo
c         do j=1,ny
c            ur(nx0,j,k)=0.d0
c         enddo
      enddo
C$OMP END PARALLEL DO
#endif
c
#ifdef SCSL
c for SCSL
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP&            PRIVATE(i,j,jimag,k,tf,tp,tx,ry_wk,cx_wk)
      do k=1,nz
c cumpute fft in x-direction
         do j=1,nyh
            jimag=nyh+j
            do i=1,nx
               tf(i) = dcmplx(ur(i,j,k),ur(i,jimag,k))
            enddo
c            tf(nxhp)=(0.d0,0.d0) !this to zero is an error
            call zzfftm(1,nx,1,1.d0,tf,nx0,tx,nx0,cx_tb,cx_wk,1)
            call zcopy(nx,tx,1,tmpf(j,1,k),nyhp0)
         enddo
c cumpute fft in y-direction
         do i=1,nx
            tmpf(nyhp,i,k)=(0.d0,0.d0)
            call zdfftm(1,ny,1,1.d0,tmpf(1,i,k),nyhp0,tp,tnyhp0,
     &           ry_tb_b,ry_wk,1)
            call dcopy(ny,tp,1,ur(i,1,k),nx0)
         enddo
      enddo
C$OMP END PARALLEL DO
#endif
c
#ifdef ESSL
c for ESSL
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,jimag,k,tf,tp,fft_work)
      do k=1,nz
c cumpute fft in x-direction
         do j=1,nyh
            jimag=nyh+j
            do i=1,nx
               tf(i) = dcmplx(ur(i,j,k),ur(i,jimag,k))
            enddo
c            tf(nxhp)=(0.d0,0.d0) !this to zero is an error
            call dcft(0,tf,1,nx0,tmpf(j,1,k),nyhp0,1,nx,1,-1,1.d0,
     &           cx_tb_b,nfft_aux,fft_work,nfft_aux)
         enddo
c cumpute fft in y-direction
         do i=1,nx
            tmpf(nyhp,i,k)=(0.d0,0.d0)
            call dcrft(0,tmpf(1,i,k),nyhp0,tp,tnyhp0,ny,1,-1,1.d0,
     &           ry_tb_b,nfft_aux,fft_work,nfft_aux)
            call dcopy(ny,tp,1,ur(i,1,k),nx0)
         enddo
      enddo
C$OMP END PARALLEL DO
#endif
c
      if (FLTIMER) call end_timer(75)
      return
      end
c
c
      subroutine fftl_r2f_2d_new(ur)
      use common_m
      REAL*8 ur(nxl0,nyl0,nz)
c
c     2-d real to fourier transform in x and y directions
c     On input:  ur contains data on 1:nxl,1:nyl,1:nz       (real component)
c     On output: ur contains data on (1:nxl),nylhp:nyl,1:nz (imag component)
c
      INTEGER i,j,k,jimag
c      COMPLEX*16 tmp(nylhp0,nxl0,nz)!this is a common block now ->tmplf
c
#ifdef SCSL
c for SCSL
      REAL*8 tp(tnylhp0)
      COMPLEX*16 tf(nxl0)
      COMPLEX*16 tx(nxl0)
      REAL*8 ry_wk(tnylhp0),cx_wk(tnxl)
#endif
c
#ifdef ESSL
c for ESSL
c fft_work must be 8-byte real array
      REAL*8 tp(tnylhp0)
      COMPLEX*16 tf(nxl0)
      REAL*8 fft_work(nfft_aux)
#endif
c
      if (FLTIMER) call begin_timer(71)
c
#ifdef FFTW3
c FFTW3
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,jimag,k)
      do k=1,nz
c scale input data for normalization
         call dscal(nxl0yl,1.d0/dble(nxl)/dble(nyl),ur(1,1,k),1)
c compute fft in y-direction
         call dfftw_execute_dft_r2c(plan_nxl_r2c_y
     &        ,ur(1,1,k),tmplf(1,1,k))
c compute fft in x-direction
         call dfftw_execute_dft(plan_nylh_c2cf_x
     &        ,tmplf(1,1,k),tmplf(1,1,k))         
c copy to output variable
         do j=1,nylh 
            jimag=nylh+j
            do i=1,nxl
               ur(i,j    ,k) = dreal(tmplf(j,i,k))
               ur(i,jimag,k) = dimag(tmplf(j,i,k))
            enddo
         enddo
      enddo
C$OMP END PARALLEL DO
#endif
c
#ifdef SCSL
c for SCSL
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP&            PRIVATE(i,j,jimag,k,tp,tf,tx,cx_wk,ry_wk)
      do k=1,nz
c compute fft in y-direction
         do i=1,nxl
            call dcopy(nyl,ur(i,1,k),nxl0,tp,1)
            call dzfftm(-1,nyl,1,1.d0/dble(nyl),tp,tnylhp0,tmplf(1,i,k),
     &           nylhp0,ry_ltb_f,ry_wk,1)
         enddo
c compute fft in x-direction
         do j=1,nylh 
            call zcopy(nxl,tmplf(j,1,k),nylhp0,tx,1)
            call zzfftm(-1,nxl,1,1.d0/dble(nxl),tx,nxl0,tf,nxl0,
     &           cx_ltb,cx_wk,1)
            jimag=nylh+j
            do i=1,nxl
               ur(i,j    ,k) = dreal(tf(i))
               ur(i,jimag,k) = dimag(tf(i))   
            enddo
         enddo
      enddo
C$OMP END PARALLEL DO
#endif
c
#ifdef ESSL
c for ESSL
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,jimag,k,tp,tf,fft_work)
      do k=1,nz
c compute fft in y-direction
         do i=1,nxl
            call dcopy(nyl,ur(i,1,k),nxl0,tp,1)
            call drcft(0,tp,tnylhp0,tmplf(1,i,k),nylhp0,nyl,1,1,
     &           1.d0/dble(nyl),ry_ltb_f,nfft_aux,fft_work,nfft_aux)
         enddo
c compute fft in x-direction
         do j=1,nylh 
            call dcft(0,tmplf(j,1,k),nylhp0,1,tf,1,nxl0,nxl,1,1,
     &           1.d0/dble(nxl),cx_ltb_f,nfft_aux,fft_work,nfft_aux)
            jimag=nylh+j
            do i=1,nxl
               ur(i,j    ,k) = dreal(tf(i))
               ur(i,jimag,k) = dimag(tf(i))   
            enddo
         enddo
      enddo
C$OMP END PARALLEL DO
#endif
c
      if (FLTIMER) call end_timer(71)
      return
      end
c
c
      subroutine fftl_f2r_2d_new(ur)
      use common_m
      REAL*8 ur(nxl0,nyl0,nz)
c
c     2-d fourier to real transform in x and y directions
c     On input:  ur contains data on 1:nxl,    1:nylh,1:nz (real component)
c     On input:  ur contains data on 1:nxl,nylhp:nyl ,1:nz (imag component)
c     On output: ur contains data on 1:nxl,1:nyl ,1:nz
c
      INTEGER i,j,k,jimag
c      COMPLEX*16 tmp(nylhp0,nxl0,nz)!this is a common block now ->tmplf
c
#ifdef SCSL
c for SCSL
      REAL*8 tp(tnylhp0)
      COMPLEX*16 tf(nxl0)
      COMPLEX*16 tx(nxl0)
      REAL*8 ry_wk(tnylhp0),cx_wk(tnxl)
#endif
c
#ifdef ESSL
c for ESSL
c fft_work must be 8-byte real array
      REAL*8 tp(tnylhp0)
      COMPLEX*16 tf(nxl0)
      REAL*8 fft_work(nfft_aux)
#endif
c
      if (FLTIMER) call begin_timer(76)
c
#ifdef FFTW3
c FFTW3
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,jimag,k)
      do k=1,nz
         do j=1,nylh
            jimag=nylh+j
            do i=1,nxl
               tmplf(j,i,k) = dcmplx(ur(i,j,k),ur(i,jimag,k))
            enddo
         enddo
c cumpute fft in x-direction
         call dfftw_execute_dft(plan_nylh_c2cb_x,
     &        tmplf(1,1,k),tmplf(1,1,k))
c cumpute fft in y-direction
         do i=1,nxl
            tmplf(nylhp,i,k)=dcmplx(0.d0)
         enddo
         call dfftw_execute_dft_c2r(plan_nxl_c2r_y
     &        ,tmplf(1,1,k),ur(1,1,k))         
c         do i=1,nxl0
c            ur(i,nyl0,k)=0.d0
c         enddo
c         do j=1,nyl
c            ur(nxl0,j,k)=0.d0
c         enddo
      enddo
C$OMP END PARALLEL DO
#endif
c
#ifdef SCSL
c for SCSL
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP&            PRIVATE(i,j,jimag,k,tf,tp,tx,ry_wk,cx_wk)
      do k=1,nz
c cumpute fft in x-direction
         do j=1,nylh
            jimag=nylh+j
            do i=1,nxl
               tf(i) = dcmplx(ur(i,j,k),ur(i,jimag,k))
            enddo
c            tf(nxlhp)=(0.d0,0.d0) !this to zero is an error
            call zzfftm(1,nxl,1,1.d0,tf,nxl0,tx,nxl0,cx_ltb,cx_wk,1)
            call zcopy(nxl,tx,1,tmplf(j,1,k),nylhp0)
         enddo
c cumpute fft in y-direction
         do i=1,nxl
            tmplf(nylhp,i,k)=(0.d0,0.d0)
            call zdfftm(1,nyl,1,1.d0,tmplf(1,i,k),nylhp0,tp,tnylhp0,
     &           ry_ltb_b,ry_wk,1)
            call dcopy(nyl,tp,1,ur(i,1,k),nxl0)
         enddo
      enddo
C$OMP END PARALLEL DO
#endif
c
#ifdef ESSL
c for ESSL
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,jimag,k,tp,tf,fft_work)
      do k=1,nz
c compute fft in y-direction
         do i=1,nxl
            call dcopy(nyl,ur(i,1,k),nxl0,tp,1)
            call drcft(0,tp,tnylhp0,tmplf(1,i,k),nylhp0,nyl,1,1,
     &           1.d0/dble(nyl),ry_ltb_f,nfft_aux,fft_work,nfft_aux)
         enddo
c compute fft in x-direction
         do j=1,nylh 
            call dcft(0,tmplf(j,1,k),nylhp0,1,tf,1,nxl0,nxl,1,1,
     &           1.d0/dble(nxl),cx_ltb_f,nfft_aux,fft_work,nfft_aux)
            jimag=nylh+j
            do i=1,nxl
               ur(i,j    ,k) = dreal(tf(i))
               ur(i,jimag,k) = dimag(tf(i))   
            enddo
         enddo
      enddo
C$OMP END PARALLEL DO
#endif
c
      if (FLTIMER) call end_timer(76)
      return
      end
c
c
      subroutine fft_r2f_2d_new_slice(ur)
      use common_m
      REAL*8 ur(nx0,ny0)
c
c     2-d real to fourier transform in x and y directions
c     On input:  ur contains data on 1:nx,1:ny
c     On output: ur contains data on (1:nx),1:nyh   (real component)
c     On output: ur contains data on (1:nx),nyhp:ny (imag component)
c
      INTEGER i,j,jimag
c      COMPLEX*16 tmp(nyhp0,nx0,nz)!this is a common block now ->tmpf
c
#ifdef SCSL
c for SCSL
      REAL*8 tp(tnyhp0)
      COMPLEX*16 tf(nx0)
      COMPLEX*16 tx(nx0)
      REAL*8 ry_wk(tnyhp0),cx_wk(tnx)
#endif
c
#ifdef ESSL
c for ESSL
c fft_work must be 8-byte real array
      REAL*8 tp(tnyhp0)
      COMPLEX*16 tf(nx0)
      REAL*8 fft_work(nfft_aux)
#endif
c
c      if (FLTIMER) call begin_timer()
c
#ifdef FFTW3
c FFTW3
c scale input data for normalization
      call dscal(nx0y,1.d0/dble(nx)/dble(ny),ur,1)
c compute fft in y-direction
      call dfftw_execute_dft_r2c(plan_nx_r2c_y
     &     ,ur,tmpf(1,1,1))
c compute fft in x-direction
      call dfftw_execute_dft(plan_nyh_c2cf_x
     &     ,tmpf(1,1,1),tmpf(1,1,1))         
c copy to output variable
      do j=1,nyh 
         jimag=nyh+j
         do i=1,nx
            ur(i,j    ) = dreal(tmpf(j,i,1))
            ur(i,jimag) = dimag(tmpf(j,i,1))
         enddo
      enddo
#endif
c
#ifdef SCSL
c SCSL
c compute fft in y-direction
      do i=1,nx
         call dcopy(ny,ur(i,1),nx0,tp,1)
         call dzfftm(-1,ny,1,1.d0/dble(ny),tp,tnyhp0,tmpf(1,i,1),nyhp0,
     &        ry_tb_f,ry_wk,1)
      enddo
c compute fft in x-direction
      do j=1,nyh 
         call zcopy(nx,tmpf(j,1,1),nyhp0,tx,1)
         call zzfftm(-1,nx,1,1.d0/dble(nx),tx,nx0,tf,nx0,
     &        cx_tb,cx_wk,1)
         jimag=nyh+j
         do i=1,nx
            ur(i,j    ) = dreal(tf(i))
            ur(i,jimag) = dimag(tf(i))
         enddo
      enddo
#endif
c
#ifdef ESSL
c for ESSL
c compute fft in y-direction
      do i=1,nx
         call dcopy(ny,ur(i,1),nx0,tp,1)
         call drcft(0,tp,tnyhp0,tmpf(1,i,1),nyhp0,ny,1,1,1.d0/dble(ny),
     &        ry_tb_f,nfft_aux,fft_work,nfft_aux)
      enddo
c compute fft in x-direction
      do j=1,nyh 
         call dcft(0,tmpf(j,1,1),nyhp0,1,tf,1,nx0,nx,1,1,1.d0/dble(nx),
     &        cx_tb_f,nfft_aux,fft_work,nfft_aux)
         jimag=nyh+j
         do i=1,nx
            ur(i,j    ) = dreal(tf(i))
            ur(i,jimag) = dimag(tf(i))
         enddo
      enddo
#endif
c      if (FLTIMER) call end_timer()
      return
      end
c
c
      subroutine fft_f2r_2d_new_slice(ur)
      use common_m
      REAL*8 ur(nx0,ny0)
c
c     2-d fourier to real transform in x and y directions
c     On input:  ur contains data on 1:nx,   1:nyh,1:nz (real component)
c     On input:  ur contains data on 1:nx,nyhp:ny ,1:nz (imag component)
c     On output: ur contains data on 1:nx,1:ny,1:nz
c
      INTEGER i,j,jimag
c      COMPLEX*16 tmp(nyhp0,nx0,nz)!this is a common block now ->tmpf
c
#ifdef SCSL
c for SCSL
      REAL*8 tp(tnyhp0)
      COMPLEX*16 tf(nx0)
      COMPLEX*16 tx(nx0)
      REAL*8 ry_wk(tnyhp0),cx_wk(tnx)
#endif
c
#ifdef ESSL
c for ESSL
c fft_work must be 8-byte real array
      REAL*8 tp(tnyhp0)
      COMPLEX*16 tf(nx0)
      REAL*8 fft_work(nfft_aux)
#endif
c
c      if (FLTIMER) call begin_timer()
c
#ifdef FFTW3
c FFTW3
      do j=1,nyh
         jimag=nyh+j
         do i=1,nx
            tmpf(j,i,1) = dcmplx(ur(i,j),ur(i,jimag))
         enddo
      enddo
c cumpute fft in x-direction
      call dfftw_execute_dft(plan_nyh_c2cb_x,
     &     tmpf(1,1,1),tmpf(1,1,1))
c cumpute fft in y-direction
      do i=1,nx
         tmpf(nyhp,i,1)=dcmplx(0.d0)
      enddo
      call dfftw_execute_dft_c2r(plan_nx_c2r_y
     &     ,tmpf(1,1,1),ur(1,1))         
c      do i=1,nx0
c         ur(i,ny0,k)=0.d0
c      enddo
c      do j=1,ny
c         ur(nx0,j,k)=0.d0
c      enddo
#endif
c
#ifdef SCSL
c cumpute fft in x-direction
      do j=1,nyh
         jimag=nyh+j
         do i=1,nx
            tf(i) = dcmplx(ur(i,j),ur(i,jimag))
         enddo
c         tf(nxhp)=(0.d0,0.d0) !this to zero is an error
         call zzfftm(1,nx,1,1.d0,tf,nx0,tx,nx0,cx_tb,cx_wk,1)
         call zcopy(nx,tx,1,tmpf(j,1,1),nyhp0)
      enddo
c cumpute fft in y-direction
      do i=1,nx
         tmpf(nyhp,i,1)=(0.d0,0.d0)
         call zdfftm(1,ny,1,1.d0,tmpf(1,i,1),nyhp0,tp,tnyhp0,
     &        ry_tb_b,ry_wk,1)
         call dcopy(ny,tp,1,ur(i,1),nx0)
      enddo
#endif
c
#ifdef ESSL
c for ESSL
c cumpute fft in x-direction
      do j=1,nyh
         jimag=nyh+j
         do i=1,nx
            tf(i) = dcmplx(ur(i,j),ur(i,jimag))
         enddo
c         tf(nxhp)=(0.d0,0.d0) !this to zero is an error
         call dcft(0,tf,1,nx0,tmpf(j,1,1),nyhp0,1,nx,1,-1,1.d0,
     &        cx_tb_b,nfft_aux,fft_work,nfft_aux)
      enddo
c cumpute fft in y-direction
      do i=1,nx
         tmpf(nyhp,i,1)=(0.d0,0.d0)
         call dcrft(0,tmpf(1,i,1),nyhp0,tp,tnyhp0,ny,1,-1,1.d0,
     &        ry_tb_b,nfft_aux,fft_work,nfft_aux)
         call dcopy(ny,tp,1,ur(i,1),nx0)
      enddo
#endif
c
c      if (FLTIMER) call end_timer()
      return
      end
