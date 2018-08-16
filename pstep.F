      SUBROUTINE pstep
      use common_m

c
c     original by ted bagwell (4-91)
c     solve poisson eqn for pressure
c     the right hand side of eqn for pressure is stored in rw
c
      REAL*8 tmp1(nz),tmp2(nz),tmp3(nz),tmp4(nz)
      REAL*8 tmp5(nzmm),tmp6(nzmm),tmp7(nzmm),tmp8(nzmm)
      REAL*8 p1(nz),p2(nz),dp1dz,dp2dz
      INTEGER i,j,jimag,k,l
      REAL*8 adiagr,fac
      REAL*8 alpha,beta
c
      REAL*8 ddot
      EXTERNAL ddot
c
      if (FLTIMER) call begin_timer(65)
c     setting BC for pressure at wall = 0
      rw(1,1,1)  = 0.d0
      rw(1,nyhp,1)  = 0.d0
      rw(1,1,nz) = 0.d0
      rw(1,nyhp,nz) = 0.d0
c
      do k=2,nzm
        tmp1(k) = rw(1,1,k)
        tmp2(k) = rw(1,1,k) - 1.d0
      enddo
      tmp1(1)  = 0.d0
      tmp2(1)  = 0.d0
      tmp1(nz) = 0.d0
      tmp2(nz) = 0.d0
c
      call dgemv('N',nz,nz,1.d0,ci_nodr,nz0,tmp1,1,0.d0,tmp3,1)
      call dgemv('N',nz,nz,1.d0,ci_nodr,nz0,tmp2,1,0.d0,tmp4,1)
c
      do l = 1,nz
        tmp3(l) = tmp3(l)/e_nodr(l)
        tmp4(l) = tmp4(l)/e_nodr(l)
      enddo
c
      call dgemv('N',nz,nz,1.d0,cn_nodr,nz0,tmp3,1,0.d0,p1,1)
      call dgemv('N',nz,nz,1.d0,cn_nodr,nz0,tmp4,1,0.d0,p2,1)
c
c     compute dp1dz at z=-1 and dp2dz at z=-1
c
      dp1dz = ddot(nz,p1,1,chbd1(nz,1),nz0)
      dp2dz = ddot(nz,p2,1,chbd1(nz,1),nz0)
c
      alpha = dp2dz/(dp2dz-dp1dz)
      beta  = dp1dz/(dp1dz-dp2dz)
      if (mod(itime-istart,oi_cfl).eq.0) write(26,900) alpha,beta
c
      do k=1,nz
        p(1,1,k) = alpha*p1(k)+beta*p2(k)
        p(1,nyhp,k) = 0.d0
      enddo
c
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP&  PRIVATE(i,j,jimag,k,adiagr,fac,tmp5,tmp6,tmp7,tmp8)
      do j = 1,nyh
        jimag=nyh+j
        do i = 1,nx
          if (i.ne.1.or.j.ne.1) then
c
            call dcopy(nzmm,rw(i,j,2),nx0y0,tmp5,1)
            call dcopy(nzmm,rw(i,jimag,2),nx0y0,tmp6,1)
c
            call dgemv('N',nzmm,nzmm,1.d0,ci_nono,nzmm0,tmp5,1,0.d0,
     &        tmp7,1)
            call dgemv('N',nzmm,nzmm,1.d0,ci_nono,nzmm0,tmp6,1,0.d0,
     &        tmp8,1)
c
            adiagr = wavexs(i) + waveys(j)
            do k=1,nzmm
              fac = 1.d0/(e_nono(k)-adiagr)
              tmp5(k)=tmp7(k)*fac
              tmp6(k)=tmp8(k)*fac
            enddo
c
            call dgemv('N',nzmm,nzmm,1.d0,cn_nono,nzmm0,tmp5,1,0.d0,
     &        tmp7,1)
            call dgemv('N',nzmm,nzmm,1.d0,cn_nono,nzmm0,tmp6,1,0.d0,
     &        tmp8,1)
c
            call dcopy(nzmm,tmp7,1,p(i,j,2),nx0y0)
            call dcopy(nzmm,tmp8,1,p(i,jimag,2),nx0y0)
c
            p(i,j,1)  = ddot(nzmm,prbc1,1,tmp7,1)
            p(i,jimag,1)  = ddot(nzmm,prbc1,1,tmp8,1)
            p(i,j,nz) = ddot(nzmm,prbcn,1,tmp7,1)
            p(i,jimag,nz) = ddot(nzmm,prbcn,1,tmp8,1)
          endif
        enddo
      enddo
C$OMP END PARALLEL DO
c
      if (FLTIMER) call end_timer(65)
      return
 900  FORMAT('pstep: alpha = ',f12.8,' beta = ',f12.8)
      END
