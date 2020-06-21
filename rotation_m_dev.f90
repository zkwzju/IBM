module precision_m
  implicit none

  integer,parameter :: wp=selected_real_kind(15,307)
  private
  public :: wp
end module precision_m
! -------------------------------------------------- !
module rotation_m
  use precision_m
  implicit none

  type qtrn
     real(wp) :: w
     real(wp) :: x
     real(wp) :: y
     real(wp) :: z
  end TYPE qtrn
  
  private
  public :: rotation_int,rotation_int2,rotation_int3,rotation_leapfrog
  
contains
    subroutine rotation_int3(om,i,j,k,torq,Ixyz,dt,om_b)
    real(wp), intent(inout)  :: om(3),i(3),j(3),k(3)
    real(wp), intent(in)     :: torq(3)
    real(wp), intent(in)     :: Ixyz(3)
    real(wp), intent(in)     :: dt
    real(wp), intent(out),optional :: om_b(3)
    
    real(wp) :: mom(3)                   
    real(wp) :: om2(3),om3(3),om4(3)            ! angular velocity
    real(wp) :: p1(3)                           ! rotation angle
    real(wp) :: phi2
    real(wp) :: fac,fac_33(3,3),fac_31(3)       ! temporary coef
    real(wp) :: i1(3),j1(3),k1(3)               ! triad of unit vectors
    real(wp) :: A(3,3),AT(3,3)                  ! rotation matrix
    real(wp) :: A1(3,3),A1T(3,3)                ! rotation matrix
    real(wp) :: I_m(3,3), Iinv_m(3,3)
    real(wp) :: sinf,cosf,phi
    I_m    = reshape((/ Ixyz(1),0.0_wp,0.0_wp, 0.0_wp,Ixyz(2),0.0_wp, 0.0_wp,0.0_wp,Ixyz(3) /), shape(I_m))
    Iinv_m = reshape((/ 1.0_wp/Ixyz(1),0.0_wp,0.0_wp, 0.0_wp,1.0_wp/Ixyz(2),0.0_wp, 0.0_wp,0.0_wp,1.0_wp/Ixyz(3) /), shape(Iinv_m)) ! inverse of I_m
    A(:,1)  = i
    A(:,2)  = j
    A(:,3)  = k
    AT = transpose(A)
    mom = torq * dt      ! angular momentum

    ! step 1
    p1   = om*dt         ! time step dt
    phi  = sqrt(sum(p1**2)) 
    phi2 = phi**2
    sinf = sin(phi)
    cosf = cos(phi)
    
    i1 = (1-cosf)*(sum(p1*i)/phi2) * p1 + i*cosf + cross(p1,i)*sinf/phi
    j1 = (1-cosf)*(sum(p1*j)/phi2) * p1 + j*cosf + cross(p1,j)*sinf/phi
    k1 = (1-cosf)*(sum(p1*k)/phi2) * p1 + k*cosf + cross(p1,k)*sinf/phi

    
    ! step 2

    A1(:,1) = i1       ! re-use intemediate variables 
    A1(:,2) = j1
    A1(:,3) = k1
    A1T = transpose(A1)
    fac_33 = matmul(matmul(A1,Iinv_m),A1T)
    fac_31 = matmul(matmul(matmul(A,I_m),AT),om) + mom

    om = matmul(fac_33,fac_31)

    i = i1           ! update body triad of unit vectors
    j = j1
    k = k1
    
    if(present(om_b)) then
       A1(:,1) = i1
       A1(:,2) = j1
       A1(:,3) = k1
       A1T = transpose(A1)
       om_b = matmul(A1T,om)
    end if
    
    return
  end subroutine rotation_int3
  ! __________________________________________________ !

  subroutine rotation_int2(om,i,j,k,torq,Ixyz,dt,om_b)
    real(wp), intent(inout)  :: om(3),i(3),j(3),k(3)
    real(wp), intent(in)     :: torq(3)
    real(wp), intent(in)     :: Ixyz(3)
    real(wp), intent(in)     :: dt
    real(wp), intent(out), optional :: om_b(3)
    
    real(wp) :: mom(3)                   
    real(wp) :: om2(3),om3(3),om4(3)            ! angular velocity
    real(wp) :: p1(3)                           ! rotation angle
    real(wp) :: phi2
    real(wp) :: fac,fac_33(3,3),fac_31(3)       ! temporary coef
    real(wp) :: i1(3),j1(3),k1(3)               ! triad of unit vectors
    real(wp) :: A(3,3),AT(3,3)                  ! rotation matrix
    real(wp) :: A1(3,3),A1T(3,3)                ! rotation matrix
    real(wp) :: I_m(3,3), Iinv_m(3,3)

    I_m    = reshape((/ Ixyz(1),0.0_wp,0.0_wp, 0.0_wp,Ixyz(2),0.0_wp, 0.0_wp,0.0_wp,Ixyz(3) /), shape(I_m))
    Iinv_m = reshape((/ 1.0_wp/Ixyz(1),0.0_wp,0.0_wp, 0.0_wp,1.0_wp/Ixyz(2),0.0_wp, 0.0_wp,0.0_wp,1.0_wp/Ixyz(3) /), shape(Iinv_m)) ! inverse of I_m
    A(:,1)  = i
    A(:,2)  = j
    A(:,3)  = k
    AT = transpose(A)
    mom = torq * dt      ! angular momentum

    ! step 1
    p1   = om*dt         ! time step dt
    phi2 = sum(p1**2)  
    fac  = 1/sqrt(1+phi2)

    i1 = (1-fac)*(sum(p1*i)/phi2) * p1 + fac*(i + cross(p1,i))
    j1 = (1-fac)*(sum(p1*j)/phi2) * p1 + fac*(j + cross(p1,j))
    k1 = (1-fac)*(sum(p1*k)/phi2) * p1 + fac*(k + cross(p1,k))

    
    ! step 2

    A1(:,1) = i1       ! re-use intemediate variables 
    A1(:,2) = j1
    A1(:,3) = k1
    A1T = transpose(A1)
    fac_33 = matmul(matmul(A1,Iinv_m),A1T)
    fac_31 = matmul(matmul(matmul(A,I_m),AT),om) + mom

    om = matmul(fac_33,fac_31)

    i = i1           ! update body triad of unit vectors
    j = j1
    k = k1
    if(present(om_b)) then
       A1(:,1) = i1
       A1(:,2) = j1
       A1(:,3) = k1
       A1T = transpose(A1)
       om_b = matmul(A1T,om)
    end if
    
    return
  end subroutine rotation_int2
  ! __________________________________________________ !
  subroutine rotation_int(om,i,j,k,torq,Ixyz,dt,om_b)
    real(wp), intent(inout)  :: om(3),i(3),j(3),k(3)
    real(wp), intent(in)     :: torq(3)
    real(wp), intent(in)     :: Ixyz(3)
    real(wp), intent(in)     :: dt
    real(wp), intent(out),optional :: om_b(3)
    
    real(wp) :: mom(3)                   
    real(wp) :: om2(3),om3(3),om4(3)            ! angular velocity
    real(wp) :: p1(3)                           ! rotation angle
    real(wp) :: phi2
    real(wp) :: fac,fac_33(3,3),fac_31(3)       ! temporary coef
    real(wp) :: i1(3),j1(3),k1(3)               ! triad of unit vectors
    real(wp) :: A(3,3),AT(3,3)                  ! rotation matrix
    real(wp) :: A1(3,3),A1T(3,3)                ! rotation matrix
    real(wp) :: I_m(3,3), Iinv_m(3,3)

    I_m    = reshape((/ Ixyz(1),0.0_wp,0.0_wp, 0.0_wp,Ixyz(2),0.0_wp, 0.0_wp,0.0_wp,Ixyz(3) /), shape(I_m))
    Iinv_m = reshape((/ 1.0_wp/Ixyz(1),0.0_wp,0.0_wp, 0.0_wp,1.0_wp/Ixyz(2),0.0_wp, 0.0_wp,0.0_wp,1.0_wp/Ixyz(3) /), shape(Iinv_m)) ! inverse of I_m
    A(:,1)  = i
    A(:,2)  = j
    A(:,3)  = k
    AT = transpose(A)
    mom = torq * dt      ! angular momentum

    ! step 1
    p1   = om*dt/2        ! time step dt/2
    phi2 = sum(p1**2)
    fac  = 1/sqrt(1+phi2)

    i1 = (1-fac)*(sum(p1*i)/phi2) * p1 + fac*(i + cross(p1,i))
    j1 = (1-fac)*(sum(p1*j)/phi2) * p1 + fac*(j + cross(p1,j))
    k1 = (1-fac)*(sum(p1*k)/phi2) * p1 + fac*(k + cross(p1,k))

    ! step 2

    A1(:,1) = i1       ! re-use intemediate variables 
    A1(:,2) = j1
    A1(:,3) = k1
    A1T = transpose(A1)
    fac_33 = matmul(matmul(A1,Iinv_m),A1T)
    fac_31 = matmul(matmul(matmul(A,I_m),AT),om) + mom

    om2 = matmul(fac_33,fac_31)

    p1   = om2*dt/2      ! time step dt/2
    phi2 = sum(p1**2)
    fac  = 1/sqrt(1+phi2)

    i1 = (1-fac)*(sum(p1*i)/phi2) * p1 + fac*(i + cross(p1,i))
    j1 = (1-fac)*(sum(p1*j)/phi2) * p1 + fac*(j + cross(p1,j))
    k1 = (1-fac)*(sum(p1*k)/phi2) * p1 + fac*(k + cross(p1,k))

    ! step 3

    A1(:,1) = i1
    A1(:,2) = j1
    A1(:,3) = k1
    A1T = transpose(A1)
    fac_33 = matmul(matmul(A1,Iinv_m),A1T)
    fac_31 = matmul(matmul(matmul(A,I_m),AT),om) + mom
    
    om3 = matmul(fac_33,fac_31)

    p1   = om3*dt        ! time step dt
    phi2 = sum(p1**2)
    fac  = 1/sqrt(1+phi2)

    i1 = (1-fac)*(sum(p1*i)/phi2) * p1 + fac*(i + cross(p1,i))
    j1 = (1-fac)*(sum(p1*j)/phi2) * p1 + fac*(j + cross(p1,j))
    k1 = (1-fac)*(sum(p1*k)/phi2) * p1 + fac*(k + cross(p1,k))

    ! step 4

    A1(:,1) = i1
    A1(:,2) = j1
    A1(:,3) = k1
    A1T = transpose(A1)
    fac_33 = matmul(matmul(A1,Iinv_m),A1T)
    fac_31 = matmul(matmul(matmul(A,I_m),AT),om) + mom

    om4 = matmul(fac_33,fac_31)

    ! step 5

    om = (om + 2*om2 + 2*om3 + om4)/6    ! update angular velocity

    ! step 6
    p1   = om*dt
    phi2 = sum(p1**2)
    fac  = 1/sqrt(1+phi2)

    i1 = (1-fac)*(sum(p1*i)/phi2) * p1 + fac*(i + cross(p1,i))
    j1 = (1-fac)*(sum(p1*j)/phi2) * p1 + fac*(j + cross(p1,j))
    k1 = (1-fac)*(sum(p1*k)/phi2) * p1 + fac*(k + cross(p1,k))

    ! step 7

    i = i1           ! update body triad of unit vectors
    j = j1
    k = k1

    if(present(om_b)) then
       A1(:,1) = i1
       A1(:,2) = j1
       A1(:,3) = k1
       A1T = transpose(A1)
       om_b = matmul(A1T,om)
    end if
    
    return
  end subroutine rotation_int

  ! -------------------------------------------------- !
   subroutine rotation_leapfrog(om,om_b,i,j,k,torq,Ixyz,dt,time)
    real(wp), intent(inout)  :: om(3),om_b(3),i(3),j(3),k(3)
    real(wp), intent(in)     :: torq(3)
    real(wp), intent(in)     :: Ixyz(3)
    real(wp), intent(in)     :: dt,time

    real(wp) :: mom(3)                   
    real(wp) :: omb(3),omb0(3),omb1(3)          ! angular velocity
    real(wp) :: torqb(3)
    real(wp) :: A(3,3),AT(3,3)                  ! rotation matrix
    real(wp) :: Ap(3,3),ApT(3,3)                ! rotation matrix at next time step
    real(wp) :: W(3,3),P(3,3),D(3,3)            ! eq.(8) in ref. PRE v58,p1169,1998, Omelyan
    real(wp), parameter,dimension(3,3) :: Iunit = reshape((/1.0_wp,0.0_wp,0.0_wp, 0.0_wp,1.0_wp,0.0_wp, 0.0_wp,0.0_wp,1.0_wp/),(/3,3/))
    integer  :: nstep, step

    if(dt>0) then
       nstep = int(time/dt)
    else
       write(*,*) 'Fatal: time step is not positive'
       stop
    end if
    A(1,:)  = i
    A(2,:)  = j
    A(3,:)  = k
    AT = transpose(A)
    
    omb = matmul(A,om)        ! from inertial to body coordinates
    torqb = matmul(A,torq)
    
    ! initialize 
    call rotation_leapfrog_init(omb,omb0,Ixyz,torqb,dt)
    !write(*,*) 'omb',omb
    !write(*,*) 'omb0',omb0

    do step=1, nstep
       call rotation_leapfrog_step(omb,omb0,omb1,A,Ap,Ixyz,torqb,dt)
       A = Ap
       om_b = 0.5_wp * (omb0+omb1)     ! interpolate to time step
       omb0 = omb1

       write(*,*) sum(A(1,:)**2),sum(A(2,:)**2),sum(A(3,:)**2)
       write(200,'(4ES15.5)') step*dt,omb1
    end do

    om = matmul(AT,om_b)      ! from body to inertial coordinates
  end subroutine rotation_leapfrog
    

  ! -------------------------------------------------- !
  ! to prepare initial value for leapfrog by explicit Euler integration with half time step
  subroutine rotation_leapfrog_init(omb,omb0,Ixyz,torqb,dt)
    real(wp), intent(in)  :: omb(3),Ixyz(3),torqb(3),dt
    real(wp), intent(out) :: omb0(3)

    omb0(1) = omb(1) + 0.5_wp*dt/Ixyz(1) * (torqb(1) + (Ixyz(2)-Ixyz(3))*omb(2)*omb(3))
    omb0(2) = omb(2) + 0.5_wp*dt/Ixyz(2) * (torqb(2) + (Ixyz(3)-Ixyz(1))*omb(1)*omb(3))
    omb0(3) = omb(3) + 0.5_wp*dt/Ixyz(3) * (torqb(3) + (Ixyz(1)-Ixyz(2))*omb(1)*omb(2))
    return
  end subroutine rotation_leapfrog_init
  ! -------------------------------------------------- !
  subroutine rotation_leapfrog_step(omb,omb0,omb1,A,Ap,Ixyz,torqb,dt)
    real(wp), intent(in)  :: omb(3), omb0(3), A(3,3), Ixyz(3), torqb(3), dt
    real(wp), intent(out) :: omb1(3), Ap(3,3)
    real(wp) :: W(3,3),P(3,3),D(3,3)            ! eq.(8) in ref. PRE v58,p1169,1998, Omelyan
    real(wp), parameter,dimension(3,3) :: Iunit = reshape((/1.0_wp,0.0_wp,0.0_wp, 0.0_wp,1.0_wp,0.0_wp, 0.0_wp,0.0_wp,1.0_wp/),(/3,3/))
    real(wp)  :: fac, ombs2
    integer :: i,j
    
    W = reshape((/0.0_wp,-omb0(3),omb0(2), omb0(3),0.0_wp,-omb0(1), -omb0(2),omb0(1),0.0_wp/), shape(W))
    do i=1,3
       do j=1,3
          P(i,j) = omb0(i)*omb0(j)
       end do
    end do

    ombs2 = sum(omb0**2)

    fac = (dt**2/4)*ombs2
    D = ((1-fac)*Iunit + dt*W + (dt**2/2)*P) / (1+fac)
    
    ! time marching for ratotaion matrix
    Ap = matmul(D,A)

    call sol_ite(omb0,omb1,Ixyz,torqb,dt)
  end subroutine rotation_leapfrog_step
  ! -------------------------------------------------- !
!!$  subroutine rotation_leapfrog_step_q(omb,omb0,omb1,A,Ap,Ixyz,torqb,dt)
!!$    real(wp), intent(in)  :: omb(3), omb0(3), A(3,3), Ixyz(3), torqb(3), dt
!!$    real(wp), intent(out) :: omb1(3), Ap(3,3)
!!$    real(wp) :: G(4,4),Q(4,4)            ! eq.(8) in ref. PRE v58,p1169,1998, Omelyan
!!$    real(wp), parameter,dimension(4,4) :: I4 = reshape((/1.0_wp,0.0_wp,0.0_wp,0.0_wp, &
!!$         0.0_wp,1.0_wp,0.0_wp,0.0_wp, &
!!$         0.0_wp,0.0_wp,1.0_wp,0.0_wp, &
!!$         0.0_wp,0.0_wp,1.0_wp,0.0_wp/),(/4,4/))
!!$    real(wp)  :: fac, ombs2
!!$    
!!$    ombs2 = sum(omb0**2)
!!$    
!!$    fac = -(dt**2/4)*ombs2
!!$    D = ((1-fac)*Iunit + dt*W + (dt**2/2)*P) / (1+fac)
!!$    
!!$    ! time marching for ratotaion matrix
!!$    Ap = matmul(D,A)
!!$
!!$    call sol_ite(omb0,omb1,Ixyz,torqb,dt)
!!$  end subroutine rotation_leapfrog_step_q
  ! -------------------------------------------------- !
  subroutine sol_ite(omb0,omb1,Ixyz,torqb,dt)
    real(wp), intent(in)  :: omb0(3), Ixyz(3), torqb(3)
    real(wp), intent(out) :: omb1(3)
    real(wp), intent(in)  :: dt

    real(wp), parameter :: rtol=1.0e-6_wp, atol=1.0e-8_wp
    integer,  parameter :: ntrial=100
    real(wp) :: omb_old,omb_new
    real(wp) :: rerr, term1,term2,term3
    integer  :: i
    logical  :: ldone
    
    omb1 = omb0
    ldone = .false.

    do i=1,ntrial
       omb_old = sqrt(sum(omb1**2))
       term1 = 0.5_wp*(omb0(2)*omb0(3) + omb1(2)*omb1(3))
       term2 = 0.5_wp*(omb0(1)*omb0(3) + omb1(1)*omb1(3))
       term3 = 0.5_wp*(omb0(1)*omb0(2) + omb1(1)*omb1(2))

       omb1(1) = omb0(1) + dt/Ixyz(1) * (torqb(1) + (Ixyz(2)-Ixyz(3))*term1)
       omb1(2) = omb0(2) + dt/Ixyz(2) * (torqb(2) + (Ixyz(3)-Ixyz(1))*term2)
       omb1(3) = omb0(3) + dt/Ixyz(3) * (torqb(3) + (Ixyz(1)-Ixyz(2))*term3)
       omb_new = sqrt(sum(omb1**2))
       if(abs(omb_new - omb_old) /(omb_old + atol) < rtol) then
          ldone = .true.
          exit
       end if
    end do
    
    if(ldone .eq. .false.) then
       write(*,*) 'Fatal: failed to find solution for rotational velocity in allowed number of iterations!'
       stop
    end if
  end subroutine sol_ite
  ! -------------------------------------------------- !
  function cross(a, b)
    real(wp), dimension(3) :: cross
    real(wp), dimension(3), intent(in) :: a, b

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  end function cross
  ! -------------------------------------------------- !
  subroutine m2q(m,q)
    type(qtrn), intent(out) :: q
    real(wp), intent(in)    :: m(3,3)
    real(wp) :: tr, s
    integer  :: i
    tr = m(1,1) + m(2,2) + m(3,3)
    if(tr .ge. 0) then
       s = sqrt(tr+1)
       q%w = 0.5_wp * s
       s = 0.5_wp / s
       q%x = (m(3,2) - m(2,3)) * s
       q%y = (m(1,3) - m(3,1)) * s
       q%z = (m(2,1) - m(1,2)) * s
    else
       i = 0
       if(m(2,2) > m(1,1)) i = 1
       if(m(3,3) > m(2,2)) i = 2
       select case (i)
       case (0)
          s = sqrt((m(1,1)- (m(2,2) + m(3,3))) + 1)
          q%x = 0.5_wp * s
          s = 0.5_wp / s
          q%y = (m(1,2) - m(2,1)) * s
          q%z = (m(3,1) - m(1,3)) * s
          q%w = (m(3,2) - m(2,3)) * s
       case (1)
          s = sqrt((m(2,2)- (m(3,3) + m(1,1))) + 1)
          q%y = 0.5_wp * s
          s = 0.5_wp / s
          q%z = (m(2,3) - m(3,2)) * s
          q%x = (m(1,2) - m(2,1)) * s
          q%w = (m(1,3) - m(3,1)) * s
       case (2)
          s = sqrt((m(3,3)- (m(1,1) + m(2,2))) + 1)
          q%z = 0.5_wp * s
          s = 0.5_wp / s
          q%x = (m(1,2) - m(2,1)) * s
          q%y = (m(3,1) - m(1,3)) * s
          q%w = (m(3,2) - m(2,3)) * s
       end select
    end if
          
  end subroutine m2q
  ! -------------------------------------------------- !
  subroutine q2m(q,m)
    type(qtrn), intent(in) :: q
    real(wp),   intent(out)    :: m(3,3)
    real(wp) :: w,x,y,z
    w=q%w
    x=q%x
    y=q%y
    z=q%z
    m(1,1) = 1.d0-2.d0*(y**2+z**2)
    m(1,2) = 2.d0*(x*y-w*z)
    m(1,3) = 2.d0*(x*z+w*y)
    m(2,1) = 2.d0*(x*y+w*z)
    m(2,2) = 1.d0-2.d0*(x**2+z**2)
    m(2,3) = 2.d0*(y*z-w*x)
    m(3,1) = 2.d0*(x*z-w*y)
    m(3,2) = 2.d0*(y*z+w*x)
    m(3,3) = 1.d0-2.d0*(x**2+y**2)

  end subroutine q2m
  
end module rotation_m


program test
  use precision_m
  use rotation_m
  implicit none
  real(wp) :: Ixyz(3)
  real(wp) :: om(3), torq(3), i(3),j(3),k(3),om_b(3)
  real(wp) :: dt, time
  integer  :: step, nt
  Ixyz = (/1,1,2/)
  om   = (/0,1,100/)
  torq = (/0,0,0/)
  i = (/1,0,0/)
  j = (/0,1,0/)
  k = (/0,0,1/)
  dt = 0.01_wp
  time = 0.25_wp
  nt = int(time/dt)
  write(100,'(4ES15.5)') 0,om
  do step=1,nt
     
     call rotation_int2(om,i,j,k,torq,Ixyz,dt,om_b)
     !write(*,*) step*dt,'-----------------'
     write(100,'(4ES15.5)') step*dt,om_b
  end do
  
  om   = (/0,1,100/)
  i = (/1,0,0/)
  j = (/0,1,0/)
  k = (/0,0,1/)
  call rotation_leapfrog(om,om_b,i,j,k,torq,Ixyz,dt,time)
end program test
