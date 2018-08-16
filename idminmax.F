C NOTE: idmax and idmin are not part of blas library and we needed to 
C       implement our own to compile with standard blas.
C       These two subroutines have been downloaded form the web 
      integer function auxidmax(n, x, incx)
c
c     PURPOSE
c        finds the index of element having maximum value
c        without taking into account the nan(s)...
c
c     NOTES
c        - original version modified by Bruno (for the nan problem...)
c          (01/01/2003)
c        - this function return 1 if x has only nan components : may be
c          this is not a good behavior
c        - this function doesn't test if n<1 or incx<1 : this is done
c          by the scilab interface
c
      implicit none
      integer n, incx
      real*8 x(incx,*)

      real*8 xmax
      integer i, j
c      external isnan
c      logical  isnan

      auxidmax = 1

c     initialize the max with the first component being not a nan
      j = 1      
#ifdef GF41
#else
      do while ( isnan(x(1,j)) )
         j = j + 1
         if ( j .gt. n ) return
      enddo
#endif
      xmax = x(1,j)
      auxidmax = j

c     the usual loop
      do i = j+1, n
         if ( x(1,i) .gt. xmax) then!a test with a nan must always return false
            xmax = x(1,i)
            auxidmax = i
         endif
      enddo
 
      end
c
      integer function auxidmin(n, x, incx)
c
c     PURPOSE
c        finds the index of the first element having minimum value
c        without taking into account the nan(s)
c
c     NOTE
c        - original version modified by Bruno (for the nan problem...)
c          (01/01/2003)
c        - this function return 1 if x has only nan components : may be
c          this is not a good behavior
c        - this function doesn't test if n<1 or incx<1 : this is done
c          by the scilab interface
c
      implicit none
      integer n, incx
      real*8 x(incx,*)

      real*8 xmin
      integer i, j
c      external isnan
c      logical  isnan

      auxidmin = 1

c     initialize the min with the first component being not a nan
      j = 1      
#ifdef GF41
#else
      do while ( isnan(x(1,j)) )
         j = j + 1
         if ( j .gt. n ) return
      enddo
#endif
      xmin = x(1,j)
      auxidmin = j

c     the usual loop
      do i = j+1, n
         if ( x(1,i) .lt. xmin) then!a test with a nan must always return false
            xmin = x(1,i)
            auxidmin = i
         endif
      enddo
 
      end
