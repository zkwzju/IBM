module interp_m
  use precision_m
  implicit none

contains
  subroutine Chebyshev_uniform(xa,ya,n,x,y,dy_out)
    real(wp), intent(in)  :: xa(n),ya(n)
    integer , intent(in)  :: n
    real(wp), intent(in)  :: x
    real(wp), intent(out) :: y
    real(wp), optional    :: dy_out
    !integer, parameter    :: NMAX=10
    INTEGER :: i,m,ns
    REAL(wp) :: den,dif,dift,ho,hp,w,dy
    !real(wp) :: c(NMAX),d(NMAX)
    real(wp),allocatable  :: c(:),d(:)

    allocate(c(n),d(n))
    
    ns=1
    dif=abs(x-xa(1))

    do i=1,n
       dift=abs(x-xa(i))
       if (dift.lt.dif) then
          ns=i
          dif=dift
       endif
       c(i)=ya(i) 
       d(i)=ya(i)
    enddo
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
       do i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.0_wp) pause 'failure in Chebyshev_uniform'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
       enddo
       if (2*ns.lt.n-m) then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       endif
       y=y+dy
    enddo
    if(present(dy_out)) dy_out = dy

    deallocate(c,d)
    return
  end subroutine Chebyshev_uniform
  ! -------------------------------------------------- !
  subroutine spline(x,y,n,yp1,ypn,y2)
    real(wp), intent(in)  :: x(n), y(n)
    integer,  intent(in)  :: n
    real(wp), intent(in)  :: yp1,ypn
    real(wp), intent(out) :: y2(n)
    !integer, parameter :: NMAX=500
    integer  :: i,k
    real(wp) :: p,qn,sig,un
    !real(wp) :: u(NMAX)
    real(wp), allocatable :: u(:)

    allocate(u(n))
    
    if (yp1.gt..99e30) then
       y2(1)=0.0_wp
       u(1)=0.0_wp
    else
       y2(1)=-0.5_wp
       u(1)=(3.0_wp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif
    do i=2,n-1
       sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
       p=sig*y2(i-1)+2.0_wp
       y2(i)=(sig-1.0_wp)/p
       u(i)=(6.0_wp*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo
    if (ypn.gt..99e30) then
       qn=0.0_wp
       un=0.0_wp
    else
       qn=0.50_wp
       un=(3.0_wp/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif
    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0_wp)
    do k=n-1,1,-1
       y2(k)=y2(k)*y2(k+1)+u(k)
    enddo
    deallocate(u)
    return
  end subroutine spline
  ! -------------------------------------------------- !
  subroutine uniform_Chebyshev(xa,ya,y2a,n,x,y)
    integer,  intent(in)  ::  n
    real(wp), intent(in)  :: xa(n),ya(n),y2a(n),x
    real(wp), intent(out) :: y
    integer  :: k,khi,klo
    real(wp) :: a,b,h
    klo=1
    khi=n
1    if (khi-klo.gt.1) then
       k=(khi+klo)/2
       if(xa(k).lt.x) then          ! xa in decreasing order
          khi=k
       else
          klo=k
       endif
       goto 1
    endif
    h=xa(khi)-xa(klo)
    if (h.eq.0.) pause 'Distinct input xa is required in uniform_Chebyshev'
    a=(xa(khi)-x)/h 
    b=(x-xa(klo))/h
    y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_wp

    return
  end subroutine uniform_Chebyshev
end module interp_m

!!$
!!$program test
!!$  use interp_m
!!$  
!!$  real(wp) :: xa(9)
!!$  real(wp) :: ya(9)
!!$  real(wp) :: x,y,dy,y2(9)
!!$  real(wp) :: yp1,ypn
!!$  integer  :: n
!!$
!!$  n=9;
!!$  
!!$  ya=[0.000000000000000E+000,0.000000000000000E+000 , 0.000000000000000E+000, -87500.0000000000 ,      -350000.000000000    ,   -87500.0000000000, 0.000000000000000E+000 , 0.000000000000000E+000 , 0.000000000000000E+000]
!!$  xa=[1.00000000000000  ,     0.750000000000000    ,   0.500000000000000,0.250000000000000   ,    0.000000000000000E+000, -0.250000000000000,-0.500000000000000   , -0.750000000000000,       -1.00000000000000]
!!$
!!$  x=0.000000000000000E+000
!!$  x= 0.382683432365090
!!$  yp1 = 0.0_wp
!!$  ypn = 0.0_wp
!!$  call spline(xa,ya,n,0.0_wp,0.0_wp,y2)
!!$
!!$  call uniform_Chebyshev(xa,ya,y2,n,x,y)
!!$  !call Chebyshev_uniform(xa,ya,n,x,y,dy)
!!$
!!$  write(*,*) ya
!!$  write(*,*) '---',y2
!!$  write(*,*)'===',y
!!$
!!$end program test

