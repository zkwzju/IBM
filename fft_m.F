      module fft_m
      use global_m
c  Common blocks for temporary storage in ffts
      COMPLEX*16, allocatable, dimension(:,:,:) :: tmpf,tmplf

c
#ifdef FFTW3
c  Support for FFT routines from FFTW library (fftw-3.1.2)
c
#include "fftw3.f"
c
      INTEGER*8 plan_nx_r2c_y,plan_nyh_c2cf_x
     &         ,plan_nx_c2r_y,plan_nyh_c2cb_x
     &         ,plan_nxl_r2c_y,plan_nylh_c2cf_x
     &         ,plan_nxl_c2r_y,plan_nylh_c2cb_x
     &         ,plan_1_r2c_x,plan_1_c2r_x
     &         ,plan_1_r2c_y,plan_1_c2r_y
c
      COMMON /fftw3p/ plan_nx_r2c_y,plan_nyh_c2cf_x
     &               ,plan_nx_c2r_y,plan_nyh_c2cb_x
     &               ,plan_nxl_r2c_y,plan_nylh_c2cf_x
     &               ,plan_nxl_c2r_y,plan_nylh_c2cb_x
     &               ,plan_1_r2c_x,plan_1_c2r_x
     &               ,plan_1_r2c_y,plan_1_c2r_y
#endif
c
#ifdef SCSL
c  Support for FFT routines from SCSL library
c
      INTEGER nxpnfr,nypnfr,tnxpfr
      PARAMETER (nxpnfr=nx+256, nypnfr=ny+256)
      PARAMETER (tnxpfr=2*nx+256)
c
      INTEGER nylpnfr,tnxlpfr
      PARAMETER (nylpnfr=nyl+256)
      PARAMETER (tnxlpfr=2*nxl+256)
c
      INTEGER tnx,tnxl
      PARAMETER(tnx=2*nx,tnxl=2*nxl)
c
c  all of these arrays must be 8-byte reals
      REAL*8 ry_tb_f(nypnfr),ry_tb_b(nypnfr)
      REAL*8 cx_tb(tnxpfr)
      REAL*8 rx_tb_f(nxpnfr),rx_tb_b(nxpnfr)
c
      REAL*8 ry_ltb_f(nylpnfr),ry_ltb_b(nylpnfr)
      REAL*8 cx_ltb(tnxlpfr)
c
      COMMON /fftscsl/  ry_tb_b,ry_tb_f,cx_tb,rx_tb_b,rx_tb_f,
     &                              ry_ltb_b,ry_ltb_f,cx_ltb
#endif
c
#ifdef ESSL
c  Support for FFT routines from ESSL library
c
      INTEGER nfft_aux
      PARAMETER (nfft_aux = 20000)
c
c  all of these arrays must be 8-byte reals
      REAL*8 cx_tb_f(nfft_aux),cx_tb_b(nfft_aux)
      REAL*8 rx_tb_f(nfft_aux),rx_tb_b(nfft_aux)
      REAL*8 ry_tb_f(nfft_aux),ry_tb_b(nfft_aux)
c
      REAL*8 cx_ltb_f(nfft_aux),cx_ltb_b(nfft_aux)
      REAL*8 ry_ltb_f(nfft_aux),ry_ltb_b(nfft_aux)
c
      COMMON /fft/ cx_tb_f,cx_tb_b,rx_tb_f,rx_tb_b,ry_tb_f,ry_tb_b,
     &	           cx_ltb_f,cx_ltb_b,ry_ltb_f,ry_ltb_b
#endif

      contains
      subroutine fft_init
        allocate(tmpf(nyhp0,nx0,nz),tmplf(nylhp0,nxl0,nz))
      end subroutine
      end module 
