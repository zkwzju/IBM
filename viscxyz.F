c
cjpf  Collapsed all subroutines to fix error in compiling with -O2
c
      SUBROUTINE viscxyz(val,dval,reyn)
      use common_m

      REAL*8 val(nx0,ny0,nz),dval(nx0,ny0,nz),reyn
c
c     compute dval = laplacian(val)/reyn
c
      INTEGER i,j,jcomp,k
      REAL*8 visc
      REAL*8     rxtmp(tnxhp0),rytmp(tnyhp0)
      COMPLEX*16 cxtmp( nxhp0),cytmp( nyhp0)
c
#ifdef SCSL
c for SCSL
      REAL*8 rx_wk(nx+2),ry_wk(ny+2)
#endif
c
#ifdef ESSL
c for ESSL
c  fft_work must be an 8-byte real array
      REAL*8 fft_work(nfft_aux)
#endif
c
c      if (FLTIMER) call begin_timer()
c
      visc = 1.d0/reyn
c
C$OMP PARALLEL DEFAULT(SHARED), SHARED(visc)
C$OMP DO PRIVATE(i,j,jcomp,k)
      do j = 1,nyh              !ny
cmic  First part, we do this in two part to reuse auxt1/auxt2(nz0,nx0,nyh)
c     traspose data
         do k = 1,nz
            call dcopy(nx,val(1,j,k),1,auxt1(k,1,j),nz0)
         enddo
c     compute
         call dgemm('N','N',nz,nx,nz,visc,chbd2,nz0,auxt1(1,1,j),nz0
     &        ,0.d0,auxt2(1,1,j),nz0)
c     transpose back data
         do i=1,nx
            call dcopy(nz,auxt2(1,i,j),1,dval(i,j,1),nx0y0)
         enddo
c     Second part
         jcomp=nyh+j
c     traspose data
         do k = 1,nz
            call dcopy(nx,val(1,jcomp,k),1,auxt1(k,1,j),nz0)
         enddo
c     compute
         call dgemm('N','N',nz,nx,nz,visc,chbd2,nz0,auxt1(1,1,j),nz0
     &        ,0.d0,auxt2(1,1,j),nz0)
c     transpose back data
         do i=1,nx
            call dcopy(nz,auxt2(1,i,j),1,dval(i,jcomp,1),nx0y0)
         enddo
c
      enddo
C$OMP END DO
c
#ifdef FFTW3
c for FFTW3
C$OMP DO PRIVATE(i,j,k,rxtmp,cxtmp,rytmp,cytmp)
      do k=1,nz
c     d2dx
         do j = 1,ny
            call dcopy(nx,val(1,j,k),1,rxtmp,1)
c scale input data for normalization
            call dscal(nx,1.d0/dble(nx),rxtmp,1)
            call dfftw_execute_dft_r2c(plan_1_r2c_x,rxtmp,cxtmp)
            do i=1,nxh
               cxtmp(i)=dcmplx(-wavexs(i))*cxtmp(i)
            enddo
            cxtmp(nxhp)=dcmplx(0.d0)
            call dfftw_execute_dft_c2r(plan_1_c2r_x,cxtmp,rxtmp)
c     
c     add to laplacian
            call daxpy(nx,visc,rxtmp,1,dval(1,j,k),1)
         enddo
c     end d2dx
c
c     d2dy
         do i = 1,nx
            call dcopy(ny,val(i,1,k),nx0,rytmp,1)
c scale input data for normalization
            call dscal(ny,1.d0/dble(ny),rytmp,1)
            call dfftw_execute_dft_r2c(plan_1_r2c_y,rytmp,cytmp)
            do j=1,nyh
               cytmp(j)=dcmplx(-waveys(j))*cytmp(j)
            enddo
            cytmp(nyhp)=dcmplx(0.d0)
            call dfftw_execute_dft_c2r(plan_1_c2r_y,cytmp,rytmp)
c     
c     add to laplacian
            call daxpy(ny,visc,rytmp,1,dval(i,1,k),nx0)
         enddo
c     end d2dy
c     rezero above nx
         do i=nx+1,nx0
            do j=1,ny
               dval(i,j,k) = 0.d0
            enddo
         enddo
      enddo
C$OMP END DO
C$OMP END PARALLEL
#endif
c
#ifdef SCSL
c for SCSL
C$OMP DO PRIVATE(i,j,k,rxtmp,cxtmp,rytmp,cytmp,rx_wk,ry_wk)
      do k=1,nz
c     d2dx
         do j = 1,ny
            call dcopy(nx,val(1,j,k),1,rxtmp,1)
            call dzfftm(-1,nx,1,1.d0/dble(nx),rxtmp,tnxhp0,cxtmp,nxhp0,
     &           rx_tb_f,rx_wk,1)
            do i=1,nxh
               cxtmp(i)=dcmplx(-wavexs(i))*cxtmp(i)
            enddo
            cxtmp(nxhp)=dcmplx(0.d0)
            call zdfftm(1,nx,1,1.d0,cxtmp,nxhp0,rxtmp,tnxhp0,
     &           rx_tb_b,rx_wk,1)
c
c     add to laplacian
            call daxpy(nx,visc,rxtmp,1,dval(1,j,k),1)
         enddo
c     end d2dx
c
c     d2dy
         do i = 1,nx
            call dcopy(ny,val(i,1,k),nx0,rytmp,1)
            call dzfftm(-1,ny,1,1.d0/dble(ny),rytmp,tnyhp0,cytmp,nyhp0,
     &           ry_tb_f,ry_wk,1)
            do j=1,nyh
               cytmp(j)=dcmplx(-waveys(j))*cytmp(j)
            enddo
            cytmp(nyhp)=dcmplx(0.d0)
            call zdfftm(1,ny,1,1.d0,cytmp,nyhp0,rytmp,tnyhp0,
     &           ry_tb_b,ry_wk,1)
c
c     add to laplacian
            call daxpy(ny,visc,rytmp,1,dval(i,1,k),nx0)
         enddo
c     end d2dy
c     rezero above nx
         do i=nx+1,nx0
            do j=1,ny
               dval(i,j,k) = 0.d0
            enddo
         enddo
      enddo
C$OMP END DO
C$OMP END PARALLEL
#endif
c
#ifdef ESSL
c for ESSL
C$OMP DO PRIVATE(i,j,k,rxtmp,cxtmp,rytmp,cytmp,fft_work)
      do k=1,nz
c     d2dx
         do j = 1,ny
            call dcopy(nx,val(1,j,k),1,rxtmp,1)
            call drcft(0,rxtmp,tnxhp0,cxtmp,nxhp0,nx,1,1,1.d0/dble(nx),
     &           rx_tb_f,nfft_aux,fft_work,nfft_aux)
            do i=1,nxh
               cxtmp(i)=dcmplx(-wavexs(i))*cxtmp(i)
            enddo
            cxtmp(nxhp)=dcmplx(0.d0)
            call dcrft(0,cxtmp,nxhp0,rxtmp,tnxhp0,nx,1,-1,1.d0,
     &           rx_tb_b,nfft_aux,fft_work,nfft_aux)
c
c     add to laplacian
            call daxpy(nx,visc,rxtmp,1,dval(1,j,k),1)
         enddo
c     end d2dx
c
c     d2dy
         do i = 1,nx
            call dcopy(ny,val(i,1,k),nx0,rytmp,1)
            call drcft(0,rytmp,tnyhp0,cytmp,nyhp0,ny,1,1,1.d0/dble(ny),
     &           ry_tb_f,nfft_aux,fft_work,nfft_aux)
            do j=1,nyh
               cytmp(j)=dcmplx(-waveys(j))*cytmp(j)
            enddo
            cytmp(nyhp)=dcmplx(0.d0)
            call dcrft(0,cytmp,nyhp0,rytmp,tnyhp0,ny,1,-1,1.d0,
     &           ry_tb_b,nfft_aux,fft_work,nfft_aux)
c
c     add to laplacian
            call daxpy(ny,visc,rytmp,1,dval(i,1,k),nx0)
         enddo
c     end d2dy
c     rezero above nx
         do i=nx+1,nx0
            do j=1,ny
               dval(i,j,k) = 0.d0
            enddo
         enddo
      enddo
C$OMP END DO
C$OMP END PARALLEL
#endif
c
c      if (FLTIMER) call end_timer()
      return
      END
c
c
      SUBROUTINE viscxyz_f(val,dval,reyn,flgsv)
      use common_m

      REAL*8 val(nx0,ny0,nz),dval(nx0,ny0,nz),reyn
      LOGICAL flgsv !this allows to use the same subroutine for vel and tt
c
c     compute dval = laplacian(val)/reyn in Fourier space
c
      INTEGER i,j,jimag,k
      REAL*8 visc,mltx,mlty
c
      if (FLTIMER) call begin_timer(90)
c
      visc = 1.d0/reyn
C$OMP PARALLEL DEFAULT(SHARED), SHARED(visc)
C$OMP DO PRIVATE(i,j,jimag,k)
      do j = 1,nyh              !ny
cmic  First part, we do this in two part to reuse auxt1/auxt2(nz0,nx0,nyh)
c     traspose data
         do k = 1,nz
            call dcopy(nx,val(1,j,k),1,auxt1(k,1,j),nz0)
         enddo
c     compute
         call dgemm('N','N',nz,nx,nz,visc,chbd2,nz0,auxt1(1,1,j),nz0
     &        ,0.d0,auxt2(1,1,j),nz0)
c     transpose back data
         do i=1,nx
            call dcopy(nz,auxt2(1,i,j),1,dval(i,j,1),nx0y0)
         enddo
c     Second part
         jimag=nyh+j
c     traspose data
         do k = 1,nz
            call dcopy(nx,val(1,jimag,k),1,auxt1(k,1,j),nz0)
         enddo
c     compute
         call dgemm('N','N',nz,nx,nz,visc,chbd2,nz0,auxt1(1,1,j),nz0
     &        ,0.d0,auxt2(1,1,j),nz0)
c     transpose back data
         do i=1,nx
            call dcopy(nz,auxt2(1,i,j),1,dval(i,jimag,1),nx0y0)
         enddo
c
      enddo
C$OMP END DO
c
      if(flgsv)then             !use anisotropic spectral viscosity
C$OMP DO PRIVATE(i,j,jimag,k,mltx,mlty)
         do k=1,nz
            do j=1,nyh
               jimag=nyh+j
               mlty=1.d0+svy(j)
               do i=1,nx
                  mltx=1.d0+svx(i)
                  dval(i,j    ,k)=dval(i,j    ,k)
     &                 -visc*(mltx*wavexs(i)+mlty*waveys(j))
     &                 *val(i,j    ,k)
                  dval(i,jimag,k)=dval(i,jimag,k)
     &                 -visc*(mltx*wavexs(i)+mlty*waveys(j))
     &                 *val(i,jimag,k)
               enddo
            enddo
         enddo
C$OMP END DO
      else
C$OMP DO PRIVATE(i,j,jimag,k)
         do k=1,nz
            do j=1,nyh
               jimag=nyh+j
               do i=1,nx
                  dval(i,j    ,k)=dval(i,j    ,k)
     &                 -visc*(wavexs(i)+waveys(j))
     &                 *val(i,j    ,k)
                  dval(i,jimag,k)=dval(i,jimag,k)
     &                 -visc*(wavexs(i)+waveys(j))
     &                 *val(i,jimag,k)
               enddo
            enddo
         enddo
C$OMP END DO
      endif
C$OMP END PARALLEL
c
      if (FLTIMER) call end_timer(90)
      return
      END
