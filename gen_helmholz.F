c  General Helmholz solver for mixed-mixed boundary conditions
c  by mic 
      SUBROUTINE gen_helmholz(cfre,var,bct,bcb,
     &     ea_mimi,cia_mimi,cna_mimi,
     &     mv1,mbct1,mbcb1,mvn,mbctn,mbcbn)
      use common_m
      REAL*8 cfre,var(nx0,ny0,nz),bct(nx0,ny0),bcb(nx0,ny0)
      REAL*8 ea_mimi(nzmm),cia_mimi(nzmm0,nzmm),cna_mimi(nzmm0,nzmm)
      REAL*8 mv1(nzmm),mvn(nzmm)
      REAL*8 mbct1,mbcb1,mbctn,mbcbn
c 
      REAL*8 fac,fac2,fac3
      INTEGER i,j,k,jimag
c
      REAL*8 ddot
      EXTERNAL ddot
c
c      if (FLTIMER) call begin_timer()
c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,jimag,k,fac,fac2,fac3)
      do j=1,nyh
         jimag=nyh+j
c
c     Modify rhs with BC
         do k=2,nzm
            do i=1,nx
               var(i,j,k)     = var(i,j    ,k)
     &            +chbd2(k,1 )*(-mbct1*bct(i,j    )+mbcb1*bcb(i,j    ))
     &            +chbd2(k,nz)*( mbctn*bct(i,j    )-mbcbn*bcb(i,j    ))
c
               var(i,jimag,k) = var(i,jimag,k)
     &            +chbd2(k,1 )*(-mbct1*bct(i,jimag)+mbcb1*bcb(i,jimag))
     &            +chbd2(k,nz)*( mbctn*bct(i,jimag)-mbcbn*bcb(i,jimag))
            enddo
         enddo
c
c     Compute multiplicative factor
         fac3 = waveys(j) + cfre
c i=1
         fac2 = wavexs(1) + fac3
         do k=1,nzmm
            fac = 1.d0/(ea_mimi(k) - fac2)
            mlt(k,1,j)      = fac
         enddo
c all the others
         do i=2,nxh
            fac2 = wavexs(i) + fac3
            do k=1,nzmm
               fac = 1.d0/(ea_mimi(k) - fac2)
               mlt(k,i,j)      = fac
               mlt(k,nx+2-i,j) = fac
            enddo
         enddo
c
c     Real component
c     transpose data
         do k = 1,nzmm
            call dcopy(nx,var(1,j,k+1),1,auxt1(k,1,j),nz0)
         enddo
c     solve
        call dgemm('N','N',nzmm,nx,nzmm,1.d0,cia_mimi,nzmm0,
     &    auxt1(1,1,j),nz0,0.d0,auxt2(1,1,j),nz0)
        do i = 1,nx
          do k = 1,nzmm
            auxt2(k,i,j) = auxt2(k,i,j)*mlt(k,i,j)
          enddo
        enddo
        call dgemm('N','N',nzmm,nx,nzmm,1.d0,cna_mimi,nzmm0,
     &    auxt2(1,1,j),nz0,0.d0,auxt1(1,1,j),nz0)
c     transpose back data
         do i=1,nx
            call dcopy(nzmm,auxt1(1,i,j),1,var(i,j,2),nx0y0)
         enddo
c     compute top and bottom
         do i = 1,nx
            var(i,j,1)  = ddot(nzmm,mv1,1,auxt1(1,i,j),1)
     &                   +mbct1*bct(i,j)-mbcb1*bcb(i,j)
            var(i,j,nz) = ddot(nzmm,mvn,1,auxt1(1,i,j),1)
     &                   -mbctn*bct(i,j)+mbcbn*bcb(i,j)
         enddo
c
c     Imaginary component
c     transpose data
         do k = 1,nzmm
            call dcopy(nx,var(1,jimag,k+1),1,auxt1(k,1,j),nz0)
         enddo
c     solve
        call dgemm('N','N',nzmm,nx,nzmm,1.d0,cia_mimi,nzmm0,
     &    auxt1(1,1,j),nz0,0.d0,auxt2(1,1,j),nz0)
        do i = 1,nx
          do k = 1,nzmm
            auxt2(k,i,j) = auxt2(k,i,j)*mlt(k,i,j)
          enddo
        enddo
        call dgemm('N','N',nzmm,nx,nzmm,1.d0,cna_mimi,nzmm0,
     &    auxt2(1,1,j),nz0,0.d0,auxt1(1,1,j),nz0)
c     transpose back data
         do i=1,nx
            call dcopy(nzmm,auxt1(1,i,j),1,var(i,jimag,2),nx0y0)
         enddo
c     compute top and bottom
         do i = 1,nx
            var(i,jimag,1)  = ddot(nzmm,mv1,1,auxt1(1,i,j),1)
     &                       +mbct1*bct(i,jimag)-mbcb1*bcb(i,jimag)
            var(i,jimag,nz) = ddot(nzmm,mvn,1,auxt1(1,i,j),1)
     &                       -mbctn*bct(i,jimag)+mbcbn*bcb(i,jimag)
         enddo
c     
      enddo
C$OMP END PARALLEL DO
c
c      if (FLTIMER) call end_timer()
      return
      END
