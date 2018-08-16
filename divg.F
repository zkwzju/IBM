      SUBROUTINE divg
      use common_m
c
c     original by bala
c
c     compute divergence of the velocity field
c     find max div(u), max div(u) on boundary, max slip (u,v,w)
c     also max slip tt if thermal field present
c
      INTEGER i,j,k,i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4,i5,j5,k5
      INTEGER ind(nz),ind1,ind2,ind3,ind3p,ind4,ind4p,ind5,ind5p
      INTEGER i6,j6,k6,ind6,ind6p
      REAL*8 val(nz),bmax,dmax,slipu,slipv,slipw,sliptt
      REAL*8 u_dum(nx0,ny0,nz),rp_dum(nx0,ny0,nz)
c
      INTEGER idamax
      EXTERNAL idamax
c
      if (FLTIMER) call begin_timer(25)
c
c     compute divergence
c
      call ddx(u,ru)
      call ddy(v,rp_dum)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(k)
      do k=1,nz
        call daxpy(nx0y,1.d0,rp_dum(1,1,k),1,ru(1,1,k),1)
      enddo
C$OMP END PARALLEL DO
      call ddz(w,rp_dum)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(k)
      do k=1,nz
        ind(k) = 0
        call daxpy(nx0y,1.d0,rp_dum(1,1,k),1,ru(1,1,k),1)
c     check that zeros have been maintained beyond data region
        do i1=nx+1,nx0
          do j1=1,ny
            if (ru(i1,j1,k).ne.0.d0) then
              ind(k) = j1
              ru(i1,j1,k) = 0.d0
            endif
          enddo
        enddo
        do j1=ny+1,ny0
          do i1=1,nx0
            if (ru(i1,j1,k).ne.0.d0) then
              ind(k) = i1
              ru(i1,j1,k) = 0.d0
            endif
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
c
      do k=1,nz
        if (ind(k).ne.0) then
          write(26,*) 'Warning: data leak for k, ind(k) = ', k, ind(k)
          go to 100
        endif
      enddo
 100  continue
c
C$OMP PARALLEL DEFAULT(SHARED),
C$OMP&  SHARED(ind3,ind3p,ind4,ind4p,ind5,ind5p,ind6,ind6p)
C$OMP DO PRIVATE(i,j,k)
      do k=1,nz
        ind(k) = idamax(nx0y,ru(1,1,k),1)
        i = mod(ind(k)-1,nx0) + 1
        j = (ind(k)-1)/nx0 + 1
        val(k) = ru(i,j,k)
      enddo
C$OMP END DO NOWAIT
C$OMP SECTIONS
C$OMP SECTION
      ind3  = idamax(nx0y,u(1,1, 1),1)
C$OMP SECTION
      ind3p = idamax(nx0y,u(1,1,nz),1)
C$OMP SECTION
      ind4  = idamax(nx0y,v(1,1, 1),1)
C$OMP SECTION
      ind4p = idamax(nx0y,v(1,1,nz),1)
C$OMP SECTION
      ind5  = idamax(nx0y,w(1,1, 1),1)
C$OMP SECTION
      ind5p = idamax(nx0y,w(1,1,nz),1)
C$OMP SECTION
      if (FLTHRM) ind6  = idamax(nx0y,tt(1,1, 1),1)
C$OMP SECTION
      if (FLTHRM) ind6p = idamax(nx0y,tt(1,1,nz),1)
C$OMP END SECTIONS
C$OMP END PARALLEL
c
c     find maximum divergence over the interior field
c
      k1 = 1 + idamax(nzmm,val(2),1)
      ind1 = ind(k1)
      bmax = val(k1)
      i1 = mod(ind1-1,nx0) + 1
      j1 = (ind1-1)/nx0 + 1
c
c     find maximum divergence on the boundary at z=+/-1
c
      k2 = 1
      ind2 = ind(1)
      if (abs(val(nz)).gt.abs(val(1))) then
        k2 = nz
        ind2 = ind(nz)
      endif
      i2 = mod(ind2-1,nx0) + 1
      j2 = (ind2-1)/nx0 + 1
      dmax = ru(i2,j2,k2)
c
c     find maximum u on boundary at z=+/-1 
c
      i3 = mod(ind3-1,nx0) + 1
      j3 = (ind3-1)/nx0 + 1
      k3 = 1
      i = mod(ind3p-1,nx0) + 1
      j = (ind3p-1)/nx0 + 1
      k = nz
c
      if (abs(u(i,j,k)).gt.abs(u(i3,j3,k3))) then
        i3 = i
        j3 = j
        k3 = k
        ind3 = ind3p
      endif
      slipu = u(i3,j3,k3)
c
c     find maximum v on boundary at z=+/-1 
c
      i4 = mod(ind4-1,nx0) + 1
      j4 = (ind4-1)/nx0 + 1
      k4 = 1
      i = mod(ind4p-1,nx0) + 1
      j = (ind4p-1)/nx0 + 1
      k = nz
c
      if (abs(v(i,j,k)).gt.abs(v(i4,j4,k4))) then
        i4 = i
        j4 = j
        k4 = k
        ind4 = ind4p
      endif
      slipv = v(i4,j4,k4)
c
c     find maximum w on boundary at z=+/-1 
c
      i5 = mod(ind5-1,nx0) + 1
      j5 = (ind5-1)/nx0 + 1
      k5 = 1
      i = mod(ind5p-1,nx0) + 1
      j = (ind5p-1)/nx0 + 1
      k = nz
c
      if (abs(w(i,j,k)).gt.abs(w(i5,j5,k5))) then
        i5 = i
        j5 = j
        k5 = k
        ind5 = ind5p
      endif
      slipw = w(i5,j5,k5)
c
c     find maximum tt on boundary at z=+/-1 
c
      sliptt = 0.d0
      if (FLTHRM) then
        i6 = mod(ind6-1,nx0) + 1
        j6 = (ind6-1)/nx0 + 1
        k6 = 1
        i = mod(ind6p-1,nx0) + 1
        j = (ind6p-1)/nx0 + 1
        k = nz
c
        if (abs(tt(i,j,k)).gt.abs(tt(i6,j6,k6))) then
          i6 = i
          j6 = j
          k6 = k
          ind6 = ind6p
        endif
        sliptt = tt(i6,j6,k6)
      endif
c
      write(26,101) bmax, i1,k1,j1
      write(26,102) dmax, i2,k2,j2
      write(26,103) slipu,i3,k3,j3
      write(26,104) slipw,i5,k5,j5
      write(26,105) slipv,i4,k4,j4
      if (FLTHRM) write(26,106) sliptt,i6,k6,j6
c
c      write(41,110) itime,abs(bmax)
c      write(42,110) itime,abs(dmax)
c      write(43,110) itime,abs(slipu)
c      write(44,110) itime,abs(slipw)
c      write(45,110) itime,abs(slipv)
c      if (FLTHRM) write(46,110) itime,abs(sliptt)
c
      if (FLTIMER) call end_timer(25)
      return
 101  FORMAT(' divg: max divergence =',1pe15.7,' at cell ',
     &       'i=',i4,' j=',i4,' k=',i4)
 102  FORMAT(' divg: max divergence on boundary=',1pe15.7,' at cell ',
     &       'i=',i4,' j=',i4,' k=',i4)
 103  FORMAT(' divg: max slip u=',1pe15.7,' at cell ',
     &       'i=',i4,' j=',i4,' k=',i4)
 104  FORMAT(' divg: max slip v=',1pe15.7,' at cell ',
     &       'i=',i4,' j=',i4,' k=',i4)
 105  FORMAT(' divg: max slip w=',1pe15.7,' at cell ',
     &       'i=',i4,' j=',i4,' k=',i4)
 106  FORMAT(' divg: max slip t=',1pe15.7,' at cell ',
     &       'i=',i4,' j=',i4,' k=',i4)
 110  FORMAT(1x,i6,5x,1pe14.7)
      END
