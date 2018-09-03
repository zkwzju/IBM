      SUBROUTINE modified_windows_fn(z_0,kk1,kk2,dum3)
      use common_m
      implicit none
c
c Calculates the coefficients to modify the ddf for nonuniform mesh
c Based on A. Pinelli 2010 "Immersed_boundary Methods for General
c Finite-difference and Finite-volume Navier-Stokes Solvers"
c
      INTEGER k,kk,kk1,kk2,lwork,i22,i33,ipiv1(3),i,j
      PARAMETER (lwork=3*(3+4))
      REAL*8 z_0,ms(5),Mi(3,3),work(lwork)
      real*8 dum3(nz)
c
      do k=1,5
         ms(k)=0
      enddo
      do i=1,3
         ibmb(i)=0.d0
      enddo
c
      do k=1,5
         do kk=kk1,kk2
            if (kk .eq. 1) then 
               ms(k)=ms(k)+(zets(1)-z_0)**(k-1)*dum3(kk)
     &                 *(zets(1)-zets(2))/2.d0
            elseif (kk .eq. nz) then
               ms(k)=ms(k)+(zets(nz)-z_0)**(k-1)*dum3(kk)
     &                 *(zets(nzm)-zets(nz))/2.d0
            else
               ms(k)=ms(k)+(zets(kk)-z_0)**(k-1)*dum3(kk)
     &                 *(zets(kk-1)-zets(kk+1))/2.d0
            endif
         enddo
      enddo
c
      do i=1,3
         do j=1,3
            Mi(i,j)=0
         enddo
      enddo
c
      Mi(1,1)=ms(1)
      Mi(1,2)=ms(2)
      Mi(1,3)=ms(3)
      Mi(2,1)=ms(2)
      Mi(2,2)=ms(3)
      Mi(2,3)=ms(4)
      Mi(3,1)=ms(3)
      Mi(3,2)=ms(4)
      Mi(3,3)=ms(5)
c
cc Inverting Mi (result stored in Mi)
      call dgetrf(3,3,Mi,3,ipiv1,i22)
      call dgetri(3,Mi,3,ipiv1,work,lwork,i33)
c
cc Multiplying Mi(inverse) and ibme
      do i=1,3
         do j=1,3
            ibmb(i)=ibmb(i)+Mi(i,j)*ibme(j)
         enddo
      enddo
c      call dgemm('N','N',3,1,3,1.d0,Mi,3,
c     &                   ibme,3,0.d0,ibmb,3)
      write(117,*) ibmb(1),ibmb(2),ibmb(3)
      return
      END SUBROUTINE modified_windows_fn
