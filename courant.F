      SUBROUTINE courant(dt0,ux,uy,uz)
      use common_m
      REAL*8 dt0
      REAL*8 ux(nx0,ny0,nz),uy(nx0,ny0,nz),uz(nx0,ny0,nz)
c
      REAL*8 bmax,idzm,nidzp,chk,cfl
      REAL*8 workx(nx0),worky(ny0),val(2:nzm)
      INTEGER i,j,k,i1,j1,k1,ind1,ind(2:nzm),imax(ny0),jmax
c
c      INTEGER idmax
c      EXTERNAL idmax
      INTEGER idamax
      EXTERNAL idamax
c
      if (FLTIMER) call begin_timer(20)
c
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP&            PRIVATE(i,j,k,idzm,nidzp,imax,jmax,workx,worky)
      do k=2,nzm
        idzm =  idz(k)
        nidzp = -idz(k-1)
        do j=1,ny
          do i=1,nx
c     uses  w dt / dz- for w > 0, -w dt / dz+ for w < 0,
c     (bear in mind the 1 -> -1 storage in z)
            workx(i) = abs(ux(i,j,k)*idx) +
     &        abs(uy(i,j,k)*idy) + max(uz(i,j,k)*idzm,uz(i,j,k)*nidzp)
          enddo
          do i=nx+1,nx0
            workx(i) = 0.d0
          enddo
c          imax(j)=idmax(nx0,workx,1)
          imax(j)=idamax(nx0,workx,1)
          worky(j)=workx(imax(j))
        enddo
        do j=ny+1,ny0
          worky(j) = 0.d0
        enddo
c        jmax   = idmax(ny0,worky,1)
        jmax   = idamax(ny0,worky,1)
        val(k) = worky(jmax)
        ind(k) = imax(jmax)+(jmax-1)*nx0
      enddo
C$OMP END PARALLEL DO
c      k1 = 1 + idmax(nzmm,val,1)
      k1 = 1 + idamax(nzmm,val,1)
      ind1 = ind(k1)
      bmax = val(k1)
      i1 = mod(ind1-1,nx0) + 1
      j1 = (ind1-1)/nx0 + 1
      chk = abs(ux(i1,j1,k1)*idx) + abs(uy(i1,j1,k1)*idy) +
     &      max(uz(i1,j1,k1)*(1.d0/(zets(k1)-zets(k1+1))),
     &          uz(i1,j1,k1)*(1.d0/(zets(k1)-zets(k1-1))))
      if (i1.gt.nx.or.abs(chk-bmax).gt.1.d-12*abs(bmax)) then
        write(26,*) 'Error in courant ', i1,nx,chk,bmax
        write(*,*) 'Error in courant ', i1,nx,abs(chk-bmax)
        if (i1.gt.nx.or.abs(chk-bmax).gt.1.d-6*abs(bmax)) stop
      endif
c
      if (dt.lt.1.d-8) write(26,*)' sum(abs(ui/dxi)) = ',bmax
      cfl = bmax*dt
      dt0 = cflmax/bmax
      write(26,101) cfl,i1,j1,k1
      write( *,101) cfl,i1,j1,k1
      write(*,*)
      if (cfl.gt.cflmax) then
        write(26,*)'courant: Warning, cfl > cflmax'
        write(26,*)'courant: dt = ',dt,' max allowed = dt0 = ',dt0
      endif
c
      if (FLTIMER) call end_timer(20)
      return
 101  FORMAT('   CFL number: ',1p1e14.5,3i5)
c 101  FORMAT(/2x,' CFL number: ',1p1e14.5,3i5)
      END
