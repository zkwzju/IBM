      SUBROUTINE del_fn(l1,n1)
      use ellipsoid_m
      use common_m
      implicit none
c
c Dirac Delta Function based on Roma et al. (1999) and Hyungoo Lee
c
      INTEGER i,i1,i2,j,j1,j2,k,k1,k2,zref,l1,n1,i_t,j_t
      REAL*8 dum1(nx),dum2(ny),dum3(nz),dum33(nz),dum333(nz)
      REAL*8 R,dzz,zminimum
      REAL*8 zarray(3),zarray1(2),hzp,hzm,hz
      REAL*8 ddf_control,dumc1,dumc2,dumc3
c
c Setting the limits of the Lagragian points domain to reduce
c computation
      i1=p_iw(l1,n1)
      i2=p_ie(l1,n1)
      j1=p_js(l1,n1)
      j2=p_jn(l1,n1)
      k1=p_kb(l1,n1)
      k2=p_kt(l1,n1)
c
c Calculating the dirac delta function in the x and y directions
      do i=i1,i2
         R=dabs((xets4(i)-x_o(l1,n1)))/deltax
         call FI_3point(R)
         if(i .lt. 1) then
            i_t=nx+i
         elseif(i .gt. nx) then
            i_t=i-nx
         else
            i_t=i
         endif
         dum1(i_t)=FI
      enddo
c
      do j=j1,j2
         R=dabs((yets4(j)-y_o(l1,n1)))/deltay
         call FI_3point(R)
         if(j .lt. 1) then
            j_t=ny+j
         elseif(j .gt. ny) then
            j_t=j-ny
         else
            j_t=j
         endif
         dum2(j_t)=FI
      enddo
c
      !>>> markers can only stay in the middle along z direction
      do k=k1,k2
         hz = abs(zets(k-1)-zets(k+1))/2.d0    
         R=dabs(zets(k)-z_o(l1,n1))/hz
         call FI_3point(R)
         dum3(k)=FI
      enddo

cC$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)

      do i=i1,i2
         if(i .lt. 1) then
            i_t=nx+i
         elseif(i .gt. nx) then
            i_t=i-nx
         else
            i_t=i
         endif
         do j=j1,j2
            if(j .lt. 1) then
               j_t=ny+j
            elseif(j .gt. ny) then
               j_t=j-ny
            else
               j_t=j
            endif
            do k=k1,k2
               ddf_dum(i_t,j_t,k)=dum1(i_t)*dum2(j_t)*dum3(k)
            enddo
         enddo
      enddo
c     C$OMP END PARALLEL DO
         

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      END SUBROUTINE del_fn
c
c
c
      SUBROUTINE FI_3point(R)
      use common_m
c 3 point delta function by Peskin (2002)
      REAL*8 R
c
      if (R .le. 0.5d0) then
         FI=1.d0/3.d0*(1.d0+sqrt(-3.d0*R**2.d0+1.d0))
      elseif (R .le. 1.5d0) then
         FI=1.d0/6.0d0*(5.d0-3.d0*R-sqrt(-3.d0*(1.d0-R)**2.d0+1.d0))
      else
         FI=0.d0
      endif
      return
      END SUBROUTINE FI_3point
c
c
c
      SUBROUTINE FI_2point(z_0,k,reference)
      use common_m
      REAL*8 z_0,sum
      INTEGER k,reference
c
      if (reference .eq. 1) then
         sum=dabs(z_0-zets(1))+dabs(z_0-zets(2))
         if (k .eq. 1) then
            FI=(zets(2)-z_0)/sum/(zets(2)-zets(1))*2.d0
         elseif (k .eq. 2) then
            FI=(z_0-zets(1))/sum/(zets(3)-zets(1))*2.d0
         else
            FI=0.d0
         endif
      endif
      if (reference .eq. nz) then
         sum=dabs(z_0-zets(nz))+dabs(z_0-zets(nzm)) 
         if (k .eq. nz) then
            FI=(zets(nzm)-z_0)/sum/(zets(nzm)-zets(nz))*2.d0
         elseif (k .eq. nzm) then
            FI=(z_0-zets(nz))/sum/(zets(nzmm)-zets(nz))*2.d0
         else
            FI=0.d0
         endif
      endif
c
      return
      END SUBROUTINE FI_2point
