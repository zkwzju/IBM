module rotation_m
  use precision_m
  implicit none

  private
  public :: rotation_leapfrog
  
contains
  ! -------------------------------------------------- !
  subroutine rotation_leapfrog(om,om_b,i,j,k,torq,Ixyz,dt,time)
    real(wp), intent(inout) :: om(3)   ! rotation vector in Eulerian framework 
    real(wp), intent(inout) :: om_b(3) ! rotation vector in body framework
    real(wp), intent(inout) :: i(3),j(3),k(3) ! particle principal axis vectors
    real(wp), intent(in)    :: torq(3) ! torque in Eulerian
    real(wp), intent(in)    :: Ixyz(3) ! inertial in body pricipal-axis
    real(wp), intent(in)    :: dt    ! integration time step
    real(wp), intent(in)    :: time  ! total integration time range

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

    do step=1, nstep
       call rotation_leapfrog_step(omb,omb0,omb1,A,Ap,Ixyz,torqb,dt)
       A = Ap
       om_b = 0.5_wp * (omb0+omb1)     ! interpolate to time step
       omb0 = omb1
    end do
    
    i = A(1,:)
    j = A(2,:)
    k = A(3,:)
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
  
end module rotation_m


!!$program test
!!$  use precision_m
!!$  use rotation_m
!!$  implicit none
!!$  real(wp) :: Ixyz(3)
!!$  real(wp) :: om(3), torq(3), i(3),j(3),k(3),om_b(3)
!!$  real(wp) :: dt, time
!!$  integer  :: step, nt
!!$  Ixyz = (/1,1,2/)
!!$  om   = (/0,1,100/)
!!$  torq = (/0,0,0/)
!!$  i = (/1,0,0/)
!!$  j = (/0,1,0/)
!!$  k = (/0,0,1/)
!!$  dt = 0.01_wp
!!$  time = 0.25_wp
!!$  nt = int(time/dt)
!!$  write(100,'(4ES15.5)') 0,om
!!$  do step=1,nt
!!$     call rotation_leapfrog(om,om_b,i,j,k,torq,Ixyz,dt,dt)
!!$     om = om_b
!!$     write(100,'(4ES15.5)') step*dt,om_b
!!$  end do
!!$
!!$
!!$end program test
