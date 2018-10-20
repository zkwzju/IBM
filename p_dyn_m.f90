module p_dyn_m
  use precision_m
  
  implicit none
contains
  subroutine part_translation(x_p,u_p,F,grav,vol,dt)
    use ellip_common_m, only: rho_p,rho_f
    real(wp),intent(inout) :: x_p(3),u_p(3)
    real(wp),intent(in) :: F(3),grav(3)
    real(wp),intent(in) :: vol,dt
    real(wp) :: u_p0(3)

    u_p0 = u_p
    u_p = u_p + dt*(rho_f*F/vol + grav*(rho_p-rho_f))/rho_p
    x_p = x_p + (u_p+u_p0)*dt/2.0_wp 
    
  end subroutine part_translation
  ! -------------------------------------------------- !
  subroutine part_translation2(x_p,u_p,F,grav,vol,dt)
    use ellip_common_m, only: rho_p,rho_f
    real(wp),intent(inout) :: x_p(3),u_p(3)
    real(wp),intent(in) :: F(3),grav(3)
    real(wp),intent(in) :: vol,dt
    real(wp) :: u_p0(3)

    u_p0 = u_p
    u_p = u_p + dt*(rho_f*F/(vol*(rho_p-rho_f)) + grav)
    x_p = x_p + (u_p+u_p0)*dt/2.0_wp 
    
  end subroutine part_translation2
  ! -------------------------------------------------- !
  subroutine part_rotation(ax1,ax2,ax3,om,om_b,torq,I_ellip,dt)
    use ellip_common_m, only: rho_p,rho_f
    use rotation_m
    
    real(wp),intent(inout) :: ax1(3),ax2(3),ax3(3) ! orientation axes
    real(wp),intent(inout) :: om(3) ! rotation in Eulerian frame
    real(wp),intent(inout) :: om_b(3) ! rotation in body frame
    real(wp),intent(in)    :: torq(3) ! Eulerian frame
    real(wp),intent(in)    :: I_ellip(3)
    real(wp),intent(in)    :: dt

    call rotation_leapfrog(om,om_b,ax1,ax2,ax3,torq,I_ellip,dt,dt)
    
  end subroutine part_rotation
  ! -------------------------------------------------- !
  subroutine part_rotation2(om,torq,I_ellip,dt)
    use ellip_common_m, only: rho_p,rho_f
    use rotation_m

    real(wp),intent(inout) :: om(3) ! rotation in Eulerian frame
    real(wp),intent(in)    :: torq(3) ! Eulerian frame
    real(wp),intent(in)    :: I_ellip(3)
    real(wp),intent(in)    :: dt

    om=om+dt*rho_p*rho_f/(I_ellip*(rho_p-rho_f))*torq
  end subroutine part_rotation2
  ! -------------------------------------------------- !
  subroutine force_on_particle(rho_f,rho_p,m_p,u,u_old,FcrossV,g,F_all)
    real(wp), intent(in)  :: rho_f,rho_p,m_p
    real(wp), intent(in)  :: u(3), u_old(3), FcrossV(3), g(3)
    real(wp), intent(out) :: F_all(3)

    F_all = -rho_f*rho_p/(m_p*(rho_p-rho_f)) * (FcrossV+g)

  end subroutine force_on_particle
    

  ! -------------------------------------------------- !

  subroutine torq_on_particle(rho_f,rho_p,m_p,u,u_old,FcrossV,g,torq_all)
    real(wp), intent(in)  :: rho_f,rho_p,m_p
    real(wp), intent(in)  :: u(3), u_old(3), FcrossV(3), g(3)
    real(wp), intent(out) :: torq_all(3)

    torq_all = -rho_f*rho_p/(m_p*(rho_p-rho_f)) * FcrossV



  end subroutine torq_on_particle
  ! -------------------------------------------------- !
  subroutine part_inner_integral(u,v,w,int_u,int_v,int_w,int_r_u,int_r_v,int_r_w,n)
    use common_m,only:xets,yets,zets,deltax,deltay,deltaz,x_c,y_c,z_c,num_p,cell_v,nx0,ny0,nz,nx,ny,num_p
    use ellip_common_m,only:xa,xb,xc,x_0,x_1,y_0,y_1,z_0,z_1
    real(wp), intent(in)  :: u(:,:,:),v(:,:,:),w(:,:,:)
    real(wp), intent(out) :: int_u,int_v,int_w,int_r_u,int_r_v,int_r_w
    integer, intent(in)   :: n
    integer :: i,j,k,id,jd
    real(wp) :: alpha,alpha_v,rx,ry,rz


    int_u = 0.0_wp
    int_v = 0.0_wp
    int_w = 0.0_wp
    int_r_u = 0.0_wp
    int_r_v = 0.0_wp
    int_r_w = 0.0_wp

    
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do i=x_0(n),x_1(n)
       if(i<1) then
          id = i+nx
       elseif(i>nx) then
          id = i-nx
       else
          id = i
       end if

       do j=y_0(n),y_1(n)
          if(j<1) then
             jd = j+ny
          elseif(i>nx) then
             jd = j-ny
          else
             jd = j
          end if

          do k=z_0(n),z_1(n)
             call volume_fraction(xets(id),yets(jd),zets(k),deltax,deltay,deltaz(k),x_c(n),y_c(n),z_c(n),xa,xb,xc,alpha)

             alpha_v = alpha*cell_v(k)

             int_u = int_u + u(id,jd,k)*alpha_v
             int_v = int_v + v(id,jd,k)*alpha_v
             int_w = int_w + w(id,jd,k)*alpha_v
             rx = xets(id) - x_c(n)
             ry = yets(jd) - y_c(n)
             rz = zets(k ) - z_c(n)
             int_r_u = int_r_u + (rx*w(id,jd,k) - rz*v(id,jd,k))*alpha_v
             int_r_v = int_r_v + (rz*u(id,jd,k) - rx*w(id,jd,k))*alpha_v
             int_r_w = int_r_w + (rx*v(id,jd,k) - ry*u(id,jd,k))*alpha_v
          end do
       end do
    end do
    !$OMP END PARALLEL DO         
  end subroutine part_inner_integral
  
  ! -------------------------------------------------- !
  ! refer to Fig.9 in Kempe & Frohlich, JCP,(2012),p3663
  subroutine volume_fraction(x,y,z,dx,dy,dz,xc,yc,zc,a,b,c,alpha)
    real(wp), intent(in)  :: x,y,z,dx,dy,dz,xc,yc,zc,a,b,c
    real(wp), intent(out) :: alpha

    real(wp) :: phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8
    real(wp) :: factor,a2,b2,c2

    a2 = a**2
    b2 = b**2
    c2 = c**2
    factor = (x-xc)**2/a2 + (y-yc)**2/b2 + (z-zc)**2/c2
    phi1 = sqrt(factor) - 1.0_wp
    phi2 = sqrt(factor + (2*dx*(x-xc)+dx**2)/a2) - 1.0_wp
    phi3 = sqrt(factor + (2*dx*(x-xc)+dx**2)/a2 + (2*dy*(y-yc)+dy**2)/b2) - 1.0_wp
    phi4 = sqrt(factor + (2*dy*(y-yc)+dy**2)/b2) - 1.0_wp

    factor = factor + (2*dz*(z-zc)+dz**2)/c2   ! lift to upper layer
    phi5 = sqrt(factor) - 1.0_wp
    phi6 = sqrt(factor + (2*dx*(x-xc)+dx**2)/a2) - 1.0_wp
    phi7 = sqrt(factor + (2*dx*(x-xc)+dx**2)/a2 + (2*dy*(y-yc)+dy**2)/b2) - 1.0_wp
    phi8 = sqrt(factor + (2*dy*(y-yc)+dy**2)/b2) - 1.0_wp


    factor = abs(phi1)+abs(phi2)+abs(phi3)+abs(phi4)+abs(phi5)+abs(phi6)+abs(phi7)+abs(phi8)
    alpha = (-phi1*H(-phi1)-phi2*H(-phi2)-phi3*H(-phi3)-phi4*H(-phi4)-phi5*H(-phi5)-phi6*H(-phi6)-phi7*H(-phi7)-phi8*H(-phi8)) / factor

  end subroutine volume_fraction
  ! -------------------------------------------------- !
  ! Heaviside function
  function H(x)
    real(wp) :: H, x
    H = sign(0.5_wp,x) + 0.5_wp
  end function H
  ! -------------------------------------------------- !


end module p_dyn_m
