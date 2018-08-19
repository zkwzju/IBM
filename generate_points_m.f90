module generate_points_m
  use ellip_common_m

contains
  subroutine generate_points
    implicit none
    integer  :: nu,nv,i,j,l,n,k
    integer  :: i1,i2,j1,j2,k1,k2
    integer  :: np

    write(*,*) 'generate_points begin ------'

    do n=1,num_p
       select case(imethod)
       case(1)
          np = n_l(n)
          call rand_ellipsoid_points(xa,xb,xc,np,xlp,ylp,zlp)
       case(2)
          call layer_ellipsoid_points(xa,xb,xc,np,xlp,ylp,zlp,mesh_length)
          n_l(n) = np

          if(np > n_ll) then
             write(*,*) 'Fatal: too many Lagrangian points'
             stop
          end if
       case(3)
          call read_ellipsoid_points(xlp,ylp,zlp,np)
          n_l(n) = np
       end select
    end do

    xlp = scale_p * xlp     ! scale particle size
    ylp = scale_p * ylp
    zlp = scale_p * zlp
    

  end subroutine generate_points

  ! -------------------------------------------------- !
  subroutine rand_ellipsoid_points(a,b,c,np,x,y,z)
    real(wp), intent(in)  :: a,b,c
    real(wp), intent(out) :: x(:),y(:),z(:)
    integer,  intent(in)  :: np
    real(wp) :: r1_v(np),r2_v(np),phi_v(np),theta_v(np),S_v(np), Sr_v(np)
    real(wp) :: S_max, r1,r2,phi,theta,S,Sr,xp,yp,zp
    integer  :: i
    logical  :: redo_flag

    call rand_num(r1_v)
    call rand_num(r2_v)

    phi_v=2*pi*r1_v-pi;
    theta_v=acos(2*r2_v-1);

    ! uniform distribution on unit sphere
    x=sin(theta_v)*cos(phi_v)
    y=sin(theta_v)*sin(phi_v)
    z=cos(theta_v)

    ! 
    S_max = max(a*b,b*c,a*c)     ! Maximum 
    S_v = sqrt((a*c*y)**2+(a*b*z)**2+(b*c*x)**2)   ! ellipsoid surface element area
    Sr_v = S_v/S_max

    do i = 1,np 
       call rand_num(r1)
       if(r1<Sr_v(i)) then    ! map the point on sphere to ellisoid 
          x(i)=a*x(i)
          y(i)=b*y(i)
          z(i)=c*z(i)
       else                   ! regenerate point on sphere and try rejection again
          redo_flag = .true.
          do while(redo_flag)
             call rand_num(r1)
             call rand_num(r2)
             phi=2*pi*r1-pi;
             theta=acos(2*r2-1);
             xp=sin(theta)*cos(phi)
             yp=sin(theta)*sin(phi)
             zp=cos(theta)

             S = sqrt((a*c*yp)**2+(a*b*zp)**2+(b*c*xp)**2)   ! ellipsoid surface element area
             Sr = S/S_max
             call rand_num(r1)
             if(r1<Sr) then
                x(i)=a*xp
                y(i)=b*yp
                z(i)=c*zp
                redo_flag = .false.
             end if
          end do
       end if
    end do

    return
  end subroutine rand_ellipsoid_points
  ! -------------------------------------------------- !
  subroutine read_ellipsoid_points(x,y,z,np)
    real(wp), intent(inout) :: x(:),y(:),z(:)
    integer,  intent(out)   :: np
    integer :: istatus

    np = 1
    open(unit=170, file=sname, STATUS='OLD',iostat=istatus)
    if(istatus .ne. 0) then
       write(*,*) 'Problem to open file',sname
       stop
    end if

    read(170,*)              ! heading line
    do
       read(170,*,IOSTAT=istatus)  x(np),y(np),z(np)
       if (istatus > 0)  then
          write(*,*) sname, 'cannot be opened'
          stop
       else if (istatus < 0) then  ! end of file reached 
          exit
       else
          np = np + 1
          if(np > n_ll+1) then
             write(*,*) 'Number of Lagrangian markers in file larger than allocated',n_ll
             stop
          end if
       end if
    end do
    np = np - 1     ! adjust
    close(170)
  end subroutine read_ellipsoid_points
  ! -------------------------------------------------- !

  ! -------------------------------------------------- !


  subroutine swap_array(a,b,s)
    real(wp), intent(inout) :: a(:),b(:)
    integer :: s
    real(wp) :: temp(s)

    temp = a
    a = b
    b = temp
  end subroutine swap_array
  ! -------------------------------------------------- !
  subroutine swap(a,b)
    real(wp), intent(inout) :: a,b
    real(wp) :: temp

    temp = a
    a = b
    b = temp
  end subroutine swap
  ! -------------------------------------------------- !
  ! to generate points along longest axis with approximate interval dh
  subroutine layer_ellipsoid_points(a_in,b_in,c_in,np,x,y,z,dh)
    real(wp), intent(in)  :: a_in,b_in,c_in,dh
    real(wp), intent(out) :: x(:),y(:),z(:)
    integer,  intent(out) :: np

    integer  :: n,i,j,nh
    real(wp) :: Circ    ! circumference
    real(wp) :: h, la, lb, lc, r, lx, dtheta, theta
    real(wp) :: a,b,c, temp   ! a >=b >=c
    logical  :: swap_ac,swap_ab,swap_bc
    a = a_in
    b = b_in
    c = c_in

    ! sort a, b, c to make a >=b >=c
    swap_ac = .false.
    swap_ab = .false.
    swap_bc = .false.

    if(a < c) then
       call swap(a,c)
       swap_ac = .true.
    end if
    if(a < b) then
       call swap(a,b)
       swap_ab = .true.
    end if
    if(b < c) then
       call swap(b,c)
       swap_bc = .true.
    end if


    np = 0
    n = int(a/dh)

    if( n*int(b/dh) > n_ll) then
       write(*,*) 'Fatal: too many Lagrangian makers due to small mesh length', n*int(b/dh),'>',n_ll
       stop
    end if

    ! to generate points between one tip and the middle plane 
    do i = 1,n-1
       lx = i*dh
       r = sqrt(1.0_wp - lx**2/a**2)
       lb = r*b
       lc = r*c
       h = (lb-lc)**2/(lb+lc)**2
       Circ = pi*(lb+lc) * (1 + 3*h/(10+sqrt(4-3*h)))
       nh = int(Circ/dh)
       dtheta = 2*pi/nh
       do j = 1,nh
          theta = j*dtheta
          np = np + 1
          x(np) = lx
          y(np) = lb * cos(theta)
          z(np) = lc * sin(theta)
       end do
    end do
    ! to generate the tip point
    np = np + 1
    x(np) = a
    y(np) = 0.0_wp
    z(np) = 0.0_wp

    ! copy to the symmetry plane
    x(np+1:2*np) = -x(1:np)
    y(np+1:2*np) = y(1:np)
    z(np+1:2*np) = z(1:np)
    np = 2*np

    ! last layer in the middle
    h = (b-c)**2/(b+c)**2
    Circ = pi*(b+c) * (1 + 3*h/(10+sqrt(4-3*h)))
    nh = int(Circ/dh)
    dtheta = 2*pi/nh

    do j = 1,nh
       theta = j*dtheta
       np = np + 1
       x(np) = 0.0_wp
       y(np) = b * cos(theta)
       z(np) = c * sin(theta)
    end do

    if(swap_bc) call swap_array(y,z,np)    
    if(swap_ab) call swap_array(x,y,np)
    if(swap_ac) call swap_array(x,z,np)

    return
  end subroutine layer_ellipsoid_points



  ! -------------------------------------------------- !

end module generate_points_m
