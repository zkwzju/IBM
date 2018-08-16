      SUBROUTINE comp_h(cfnl1,irk)
      use common_m

      REAL*8 cfnl1
      INTEGER irk
c
c     This subroutine computes h for RK-3 for the three components 
c     of velocity in real space
c     Uses the non-linear terms computed by nltrms
c     with convective and divergence forms alternately,
c     these are stored in ru, rv and rw
c     Use a variable mean pressure gradient in the streamwise direction (y), 
c     set gradpbar = 1.d0 in main (or shell script) to keep force balance
c
c
      REAL*8 thrmfac
      INTEGER i,j,k
c
      if (FLTIMER) call begin_timer(50)
c
      thrmfac = 0.d0
      if (FLTHRM) thrmfac = ra
c
c     add non-linear terms and mean pressure gradient
c     compute de-aliased non-linear terms in Fourier space, 
c     results in ru, rv and rw
      call nltrms_f
c
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do k=1,nz
        do j=1,ny
          do i=1,nx
           h1(i,j,k)=cfnl1*h1(i,j,k)-dt*ru(i,j,k)
           h2(i,j,k)=cfnl1*h2(i,j,k)-dt*rv(i,j,k)
           h3(i,j,k)=cfnl1*h3(i,j,k)-dt*rw(i,j,k)
          enddo
        enddo
c       add gradpbar to zero mode only
        h2(1,1,k)=h2(1,1,k)+dt*gradpbar
      enddo
C$OMP END PARALLEL DO
c     
c     ru, rv, rw, rtt, rp can be re-used now
c     ul, vl, and wl cannot be reused until tt_nltrms has been called
c
c     add body force terms computed in de-aliased way if gravity varies
      if(FLGVARY)then
         write(26,*)'Check implementation of varying gravity'
         write(26,*)'which is not completely suited' 
         write(26,*)'for sloping channels'
         stop
c     enlarge tt --> ddxul
         call enlarge_f_new(tt,ddxul)
         call fftl_f2r_2d_new(ddxul)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
         do k=1,nz
            do j=1,nyl
               do i=1,nxl
                  ddyul(i,j,k)=ddxul(i,j,k)*gz(j,k)
                  ddxul(i,j,k)=ddxul(i,j,k)*gy(j,k)
               enddo
            enddo
         enddo
C$OMP END PARALLEL DO
         call fftl_r2f_2d_new(ddxul)
         call reduce_f_new(ddxul,rv)
         call fftl_r2f_2d_new(ddyul)
         call reduce_f_new(ddyul,rw)
c
c     add to h
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
         do k=1,nz
            do j=1,ny 
               do i=1,nx
                  h2(i,j,k)=h2(i,j,k)+dt*thrmfac*rv(i,j,k)
                  h3(i,j,k)=h3(i,j,k)+dt*thrmfac*rw(i,j,k)
               enddo
            enddo
C NOTE by mic: here i substract the mean of the tt*gz 
C              to avoid press BC problem - channel flow only
C NOTE by mic: gz cannot be a function of y to substract mean here
            h3(1,1,k)=h3(1,1,k)+dt*thrmfac*(-rw(1,1,k)) 
         enddo
C$OMP END PARALLEL DO
      else!if gravity does not vary with y not need de-aliasing
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
         do k=1,nz
            do j=1,ny
               do i=1,nx
                  h2(i,j,k)=h2(i,j,k)+dt*thrmfac*tt(i,j,k)*gy(1,k)!grav fixed y
                  h3(i,j,k)=h3(i,j,k)+dt*thrmfac*tt(i,j,k)*gz(1,k)!grav fixed y
               enddo
            enddo
C NOTE by mic: here i substract the mean of the tt*gz 
C              to avoid press BC problem - channel flow only 
            h3(1,1,k)=h3(1,1,k)+dt*thrmfac*(-tt(1,1,k))*gz(1,k) 
         enddo
C$OMP END PARALLEL DO         
      endif
c
c     ru, rv, rw, rtt, rp can be re-used now
c     ul, vl, and wl cannot be reused until tt_nltrms has been called
c
      if (FLTIMER) call end_timer(50)
      return
      END
c
c
      SUBROUTINE rhs(cfnl2,cfdif)
      use ellipsoid_m
      use common_m
      REAL*8 cfnl2,cfdif
c
c     find right hand side of helmholtz eqn for the three components of 
c     velocity, in real space, storing results in u, v and w
c     uses h (advection + forces) computed before in comp_h
c     h has not aliasing arrors here
c
      REAL*8 rhsfac
      INTEGER i,j,k
c
      rhsfac  = -re/cfdif
c
      if (FLTIMER) call begin_timer(53)
c
c     compute 1/re*div grad u for the i-momentum equation and
c     sum the right hand side terms for the i-momentum equation, i=x,y,z
c
c      call viscxyz(u,ru,re)
c      call viscxyz(v,rv,re)
c      call viscxyz(w,rw,re)
c       do k=1,nz
c        do j=1,ny
c          do i=1,nx
c             if (abs(v(i,j,k)) .gt. 0.d0) then
c                write(118,*) i,j,k,v(i,j,k)
c             endif
c          enddo
c        enddo
c       enddo
      call viscxyz_f(u,ru,re,FLSV_NO)
      call viscxyz_f(v,rv,re,FLSV_NO)
      call viscxyz_f(w,rw,re,FLSV_NO)
cccccccccccccccccccccccccccccccccccccccccccccccc GA cccccccccccccccccc
      if (flagibm .eq. 1) then
         if (ibm_moving .eq. 1 .and. icount .eq. 1) then
            call ellipsoid(lgeneration=.false.)
            call post_ellip

         endif
         icount=1
         if (ibm_moving .eq. 2 .and. irkk .eq. 3) then
            call ellipsoid(lgeneration=.false.)
            call post_ellip

         endif
            
         if (FLTIMER) call begin_timer(54)

         call ibm_ellip(cfnl2,cfdif)

         if (FLTIMER) call end_timer(54)

C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
         do k=1,nz
            do j=1,ny
               do i=1,nx
                  u(i,j,k)=(u(i,j,k)+cfnl2*h1(i,j,k)+cfdif*ru(i,j,k)
     &                                  +dt*forcing_x(i,j,k))*rhsfac
c     &             +dt*forcing_x(i,j,k)-2.d0*cfdif*dpx(i,j,k))*rhsfac
                  v(i,j,k)=(v(i,j,k)+cfnl2*h2(i,j,k)+cfdif*rv(i,j,k)
     &                                  +dt*forcing_y(i,j,k))*rhsfac
c     &             +dt*forcing_y(i,j,k)-2.d0*cfdif*dpy(i,j,k))*rhsfac
                  w(i,j,k)=(w(i,j,k)+cfnl2*h3(i,j,k)+cfdif*rw(i,j,k)
     &                                  +dt*forcing_z(i,j,k))*rhsfac
c     &             +dt*forcing_z(i,j,k)-2.d0*cfdif*dpz(i,j,k))*rhsfac
!>>>               
c               forcing_x0(i,j,k) = (cfnl2*h1(i,j,k)+cfdif*ru(i,j,k))/dt
c     &              -2.d0*cfdif*dpx(i,j,k))/dt
c               forcing_y0(i,j,k) = (cfnl2*h2(i,j,k)+cfdif*rv(i,j,k))/dt
c     &              -2.d0*cfdif*dpy(i,j,k))/dt
c               forcing_z0(i,j,k) = (cfnl2*h3(i,j,k)+cfdif*rw(i,j,k))/dt
c     &              -2.d0*cfdif*dpz(i,j,k))/dt
!<<<
               enddo
            enddo
         enddo
C$OMP END PARALLEL DO
      elseif (flagibm .eq.0) then
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
         do k=1,nz
            do j=1,ny
               do i=1,nx
              u(i,j,k)=(u(i,j,k)+cfnl2*h1(i,j,k)+cfdif*ru(i,j,k))*rhsfac
              v(i,j,k)=(v(i,j,k)+cfnl2*h2(i,j,k)+cfdif*rv(i,j,k))*rhsfac
              w(i,j,k)=(w(i,j,k)+cfnl2*h3(i,j,k)+cfdif*rw(i,j,k))*rhsfac
               enddo
            enddo
         enddo
C$OMP END PARALLEL DO
      endif
c
      if (FLTIMER) call end_timer(53)
      return
      END
