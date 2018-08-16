c by mic
      SUBROUTINE tt_comp_h(cfnl1)
      use common_m

      REAL*8 cfnl1
c
c     This subroutine finds h for RK-3 for temperature,
c     in Furier space
c     Uses the non-linear term computed in tt_nltrms_f (in Fourier space)
c     with divergence form, this is stored in rtt
c
      INTEGER i,j,k
c
      if (FLTIMER) call begin_timer(54)
c
c     compute non-linear term in Fourier space, result is stored in rtt
      call tt_nltrms_f
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            htt(i,j,k) = cfnl1*htt(i,j,k) + dt*(-rtt(i,j,k))
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
c
c     ru, rv, rw, rtt, rp and g can be re-used now
c
      if (FLTIMER) call end_timer(54)
      return
      END
c
c
      SUBROUTINE tt_rhs(cfnl2,cfdif)
      use common_m

      REAL*8 cfnl2,cfdif
c
c     find right hand side of helmholtz eqn for temperature,
c     in Fourier space, result stored in tt
c
      REAL*8 rhsfac
      INTEGER i,j,k
c
      if (FLTIMER) call begin_timer(55)
c
      rhsfac = -re*pr/cfdif
c
c     compute 1/(re pr)*div grad tt for the temperature equation
c     and sum the right hand side terms for the temperature equation
c
c      call viscxyz(tt,rtt,re*pr)
      call viscxyz_f(tt,rtt,re*pr,FLSV_YES)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do k=1,nz
       do j=1,ny
        do i=1,nx
         tt(i,j,k)=(tt(i,j,k)+cfnl2*htt(i,j,k)+cfdif*rtt(i,j,k))*rhsfac
        enddo
       enddo
      enddo
C$OMP END PARALLEL DO
c
      if (FLTIMER) call end_timer(55)
      return
      END
