      SUBROUTINE reduce_f_new(varl,var)
      use common_m
      REAL*8 varl(nxl0,nyl0,nz)
      REAL*8 var(nx0,ny0,nz)
c     This subroutine copy in Fourier space a large variable 
c     to a small variable
c     The input varl is preserved
c
      INTEGER j,k,jimag,jimagl
c
c      if (FLTIMER) call begin_timer()
c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(j,jimag,jimagl,k)
      do k=1,nz
        do j=1,nyh !nyh=nylch
c real component
          call dcopy(nxh ,varl(     1,j,k),1,var(   1,j,k),1)
          call dcopy(nxhm,varl(cnxlch,j,k),1,var(cnxh,j,k),1)
c imaginary component
          jimag =nyh +j
          jimagl=nylh+j
          call dcopy(nxh ,varl(     1,jimagl,k),1,var(   1,jimag,k),1)
          call dcopy(nxhm,varl(cnxlch,jimagl,k),1,var(cnxh,jimag,k),1)
        enddo
      enddo
C$OMP END PARALLEL DO
c
c      if (FLTIMER) call end_timer()
      return
      END
c
c
      SUBROUTINE enlarge_f_new(var,varl)
      use common_m
      REAL*8 var(nx0,ny0,nz)
      REAL*8 varl(nxl0,nyl0,nz)
c     This subroutine copy and zero-pad in Fourier space a small 
c     variable to a large variable
c     The input var is preserved
c
      INTEGER i,j,k,jimag,jimagl
c
c      if (FLTIMER) call begin_timer()
c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,jimag,jimagl,k)
      do k=1,nz
        do j=1,nyh! nylch=nyh
c real component
          call dcopy(nxh ,var(   1,j,k),1,varl(     1,j,k),1)
          call dcopy(nxhm,var(cnxh,j,k),1,varl(cnxlch,j,k),1)
c imaginary component
          jimag =nyh +j
          jimagl=nylh+j
          call dcopy(nxh ,var(   1,jimag,k),1,varl(     1,jimagl,k),1)
          call dcopy(nxhm,var(cnxh,jimag,k),1,varl(cnxlch,jimagl,k),1)
c zero padding 
          do i=nxlchp,cnxlchm
            varl(i,j     ,k) = 0.d0
            varl(i,jimagl,k) = 0.d0
          enddo
        enddo
        do j=nylchp,nylh
          jimagl=nylh+j
          do i=1,nxl
            varl(i,j     ,k) = 0.d0
            varl(i,jimagl,k) = 0.d0
          enddo            
        enddo
      enddo
C$OMP END PARALLEL DO
c
c      if (FLTIMER) call end_timer()
      return
      END
