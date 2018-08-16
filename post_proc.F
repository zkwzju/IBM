      SUBROUTINE comp_vort(varx,vary,varz)
      use common_m

      REAL*8 varx(nx0,ny0,nz),vary(nx0,ny0,nz),varz(nx0,ny0,nz)
c
      INTEGER i,j,k
      if (FLTIMER) call begin_timer(87)
c
      call ddy(w,varx)
      call ddz(v,ru)
c
      call ddz(u,vary)
      call ddx(w,rv)
c
      call ddx(v,varz)
      call ddy(u,rw)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
       do k=1,nz
        do j=1,ny
          do i=1,nx
          varx(i,j,k)=varx(i,j,k)-ru(i,j,k)
          vary(i,j,k)=vary(i,j,k)-rv(i,j,k)
          varz(i,j,k)=varz(i,j,k)-rw(i,j,k)
          enddo
        enddo
       enddo
C$OMP END PARALLEL DO
c
      if (FLTIMER) call end_timer(87)
      return
      END
c
c
      SUBROUTINE comp_bott_shear(bott_shear_zx,bott_shear_zy)
      use common_m

      REAL*8 bott_shear_zx(nx0,ny0),bott_shear_zy(nx0,ny0)
c
      INTEGER i,j
c
      REAL*8 ddot
      EXTERNAL ddot
c
      if (FLTIMER) call begin_timer(93)
c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
      do j = 1,ny
        do i = 1,nx
         bott_shear_zx(i,j)= ddot(nz,chbd1(nz,1),nz0,u(i,j,1),nx0*ny0)
c
         bott_shear_zy(i,j)= ddot(nz,chbd1(nz,1),nz0,v(i,j,1),nx0*ny0)
        enddo
      enddo
C$OMP END PARALLEL DO
c
      if (FLTIMER) call end_timer(93)
      return
      END
c
c     
      SUBROUTINE comp_erosion(auxeros)
      use common_m

      REAL*8 auxeros(nx0,ny0)
c     
c     This subroutine computes local erosion assuming variables 
c     in Fourier space
c     Local erosion is returned in physical space 
c     The entrainment function used is Garcia and Parker (1993), JGR-Oceans
c     To get to use this subroutine (ws.ne.0)
c     Recall ws is magnitude (ws>0)
c
      REAL*8 bott_shear_zx(nx0,ny0),bott_shear_zy(nx0,ny0)
      REAL*8 ushear,auxzu5
      INTEGER i,j
c
c      if (FLTIMER) call begin_timer()
c
c     compute bottom shear stress
c     here it is assumed that u,v, w are in fourier space
      call comp_bott_shear(bott_shear_zx,bott_shear_zy)
c     transform bottom shear to physical space
      call fft_f2r_2d_new_slice(bott_shear_zx)
      call fft_f2r_2d_new_slice(bott_shear_zy)
c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,ushear,auxzu5)
      do j = 1,ny
        do i = 1,nx
c     compute ushear
           ushear = sqrt(1.d0/re*sqrt(bott_shear_zx(i,j)**2
     &                               +bott_shear_zy(i,j)**2))
c     compute zu
           auxzu5 = cferos1*(ushear/ws)*(rey_p**cferos2)
c     compute A*zu**5
           auxzu5 = 1.3d-7*(auxzu5**5)
c     compute Es
           auxeros(i,j)=0.3d0*auxzu5/(0.3d0+auxzu5)/Tmax
        enddo
      enddo
C$OMP END PARALLEL DO
c
c      if (FLTIMER) call end_timer()
      return
      END
c
      SUBROUTINE make_output_1d(LU)
      use common_m

      INTEGER LU
c
      REAL*8 bott_shear_zx,bott_shear_zy
      INTEGER iout,jout,kout
c
      REAL*8 ddot
      EXTERNAL ddot
c
      if (FLTIMER) call begin_timer(96)
c
      iout=nxh+1
      jout=nyh+1+nyh/2
      kout=nzmm
c compute bottom shear stress
      bott_shear_zx=ddot(nz,chbd1(nz,1),nz0,u(iout,jout,1),nx0*ny0)
      bott_shear_zy=ddot(nz,chbd1(nz,1),nz0,v(iout,jout,1),nx0*ny0)
c
      write(LU,300)rtime,u(iout,jout,kout),v(iout,jout,kout),
     &bott_shear_zx,bott_shear_zy
c
      if (FLTIMER) call end_timer(96)
 300  FORMAT(1p5e18.10)
      return
      END
