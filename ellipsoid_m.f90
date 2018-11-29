module ellipsoid_m
  use ellip_common_m
  use common_m, only : irkk
  implicit none
contains
  ! -------------------------------------------------- !
  subroutine read_settings
    use parser_m
    real(wp) :: g_y    

    call parser_read('IBM on flag',iflag_ibm,1)
    call parser_read('Moving wall velocity (top)',vt_l,0.0_wp)
    call parser_read('Moving wall velocity (bottom)',vb_l,0.0_wp)
    call parser_read('Initial v field',v0_l,1.0_wp)
    call parser_read('Linear v field flag',llinear_v,.false.)

    ! a,b,c are along x,y,z, respectively
    call parser_read('Ellipsoid semi-axis a',xa)
    call parser_read('Ellipsoid semi-axis b',xb)
    call parser_read('Ellipsoid semi-axis c',xc)
    call parser_read('Particle mass density',rho_p,2.0_wp)
    call parser_read('Fluid mass density',rho_f,1.0_wp)
    call parser_read('Gravity (y)',g_y,0.0_wp)
    grav=(/0.0_wp,g_y,0.0_wp/) ! gravity
    call parser_read('Method for weight dv',iweight,1)
    if(iweight .ne.1) then
       call parser_read('Update weight dv',lupdate_dv,.false.)
    else
       lupdate_dv = .false.
    end if
    call parser_read('Scale factor for particle',scale_p,1.0_wp)
    call parser_read('Scale factor for dv',scale_dv,1.0_wp)
    call parser_read('Rigid body assumption',lrigid,.true.)
    call parser_read('Translational motion',ltranslation,.true.)
    call parser_read('Rotational motion',lrotation,.false.)
    call parser_read('Prestep fraction',frac_step,0.0_wp)
    if(frac_step < 0.0_wp .or. frac_step > 1.0_wp) then
       write(*,*) 'Prestep fraction wrong value:',frac_step
       stop
    endif
    call parser_read('Non-uniform z-grid',lnonuniform,.false.)
    call parser_read('Z-grid stretching parameter',aalpha,0.99_wp)
    call parser_read('Additional checking flag',lcheck,.false.)
    call parser_read('Sphere particle',lsphere,.true.)
    call parser_read('Convection',lconvect,.false.)

    call parser_read('clip x translation',lclip_x,.false.)
    call parser_read('clip y translation',lclip_y,.false.)
    call parser_read('clip z translation',lclip_z,.false.)
    call parser_read('clip x rotation',lclip_ox,.false.)
    call parser_read('clip y rotation',lclip_oy,.false.)
    call parser_read('clip z rotation',lclip_oz,.false.)
  end subroutine read_settings

  ! -------------------------------------------------- !
  subroutine init_location
    use parser_m
    real(wp) :: x0,y0,z0
    integer  :: n

    call parser_read('Particle initial location (x)',x0,rlenx/2.0_wp);
    call parser_read('Particle initial location (y)',y0,rleny/2.0_wp);
    call parser_read('Particle initial location (z)',z0,0.0_wp);

    do n=1,num_p
       if(x0 <0 .or. y0<0 .or. z0 > rlenz/2.0_wp .or. z0 < -rlenz/2.0_wp) then
          write(*,*) "Particle is initially outside of the channel!"
          stop
       endif
       x_c(n) = mod(x0,rlenx)
       y_c(n) = mod(y0,rleny)
       z_c(n) = z0
    end do
  end subroutine init_location
  ! -------------------------------------------------- !

  subroutine initialize_ellip(ibm_moving)
    use parser_m
    implicit none
    integer, intent(inout)  :: ibm_moving
    character(len=120)      :: input,command_buffer
    real(wp) :: p,mass 
    integer  :: imove
    integer  :: i,Np
    logical  :: lflag,laxis1,laxis2,laxis3,linertia
    real(wp) :: axis1(3),axis2(3),axis3(3)

    ! particle orientation vetors
    allocate(axis_1(3,num_p),axis_2(3,num_p),axis_3(3,num_p))
    ! particle Lagrangian markers coordinates wrt particle center
    allocate(xlp(n_ll),ylp(n_ll),zlp(n_ll))
    allocate(x_0(num_p),x_1(num_p),y_0(num_p),y_1(num_p),z_0(num_p),z_1(num_p))
    allocate(for_px(num_p),for_py(num_p),for_pz(num_p))
    allocate(torq_x(num_p),torq_y(num_p),torq_z(num_p))

    ! setting control perameters

    max_len_ellip = max(max(xa,xb),xc)
    vol_ellip = 4.0_wp/3.0_wp*pi*xa*xb*xc     ! ellipsoid volume
    mass = vol_ellip * rho_p

    ! inertia of ellipsoid
    I_ellip(1) = mass/5.0_wp * (xb**2 + xc**2)
    I_ellip(2) = mass/5.0_wp * (xa**2 + xc**2)
    I_ellip(3) = mass/5.0_wp * (xa**2 + xb**2)

    ! to read rotational inertial
    call parser_is_defined('Rotational inertia',linertia)
    if(linertia) then
       call parser_read('Rotational inertia',I_ellip)
    end if

    ! initial orientation vectors
    ! axis_1
    call parser_is_defined('Principal axis-1 direction',laxis1)
    if(laxis1) then
       call parser_read('Principal axis-1 direction',axis1)
       ! normalize
       axis1 = axis1/sqrt(axis1(1)**2+axis1(2)**2+axis1(3)**2)
    else
       axis1=(/1.0_wp,0.0_wp,0.0_wp/)
    end if
    ! axis_2
    call parser_is_defined('Principal axis-2 direction',laxis2)
    if(laxis2) then
       call parser_read('Principal axis-2 direction',axis2)
       ! normalize
       axis2 = axis2/sqrt(axis2(1)**2+axis2(2)**2+axis2(3)**2)
    else
       axis2=(/0.0_wp,1.0_wp,0.0_wp/)
    end if
    ! axis3 = axis1 x axis2
    axis3(1) =  axis1(2)*axis2(3) - axis2(2)*axis1(3)
    axis3(2) = -axis1(1)*axis2(3) + axis2(1)*axis1(3)
    axis3(3) =  axis1(1)*axis2(2) - axis2(1)*axis1(2)

    do i=1,num_p      ! all particles have same initial orientation
       axis_1(:,i) = axis1
       axis_2(:,i) = axis2
       axis_3(:,i) = axis3
    end do


    ! initial rotational velocity
    om_ellip   = 0.0_wp
    om_ellip_b = 0.0_wp

    call parser_read('Movement type',imove,0)
    select case(imove)
    case(0)     ! still
       ibm_moving = 0
    case(1)     ! free movement 
       ibm_moving = 1
    case(2)     ! specific movement
       ibm_moving = 2
    end select



    call parser_read('Lagrangian markers generation method',imethod,1)
    ! approximate surface area of an ellipsoid
    p=1.6075_wp
    S_ellp=4*pi*((xa**p * xb**p + xa**p * xc**p +xb**p * xc**p)/3.0_wp)**(1.0_wp/p)

    select case(imethod)
    case(1)
       call random_seed()       ! random function seeding
       call parser_is_defined('Number of Lagrangian points',lflag)
       if(lflag) then
          call parser_read('Number of Lagrangian points',Np)
          do i=1,num_p   
             n_l(i) = Np
          end do
       else              ! default random distribution determined by Eulerian grid

          do i=1,num_p
             n_l(i) = int(S_ellp/(deltax*deltay))
          end do
       end if
       ! ------------------------------ !
       do i=1,num_p
          if(n_l(i) .eq. 0) then
             write(*,*) 'Fatal: Number of Lagranitan points is zero'
             stop
          end if
          if(n_l(i)>n_ll) then
             write(*,*) 'Fatal: Number of Lagrangian points is bigger than allocated', n_l(i),'>', n_ll
             stop
          end if
       end do
       ! -------------------------------------------------- !
    case(2)        ! non-random distribution
       call parser_is_defined('Surface mesh length',lflag)
       if(lflag) then
          call parser_read('Surface mesh length',mesh_length)
       else
          ! default surface mesh length
          mesh_length = deltax
       end if
    case(3)
       call parser_read('Lagrangian makers file',sname,'lagrangian_pts.plt')
    case default
       write(*,*) 'Method can only be 1, 2, or 3, for random, length specified, and reading from file, respectively'
       stop
    end select

    call parser_read('Lagrangian weights file',sname_w,'lagrangian_weights.dat')

!    write(*,*) 'initialize_ellip done'
  end subroutine initialize_ellip
  ! ------------------------------------------------------------ !

  subroutine ellipsoid(lgeneration)
    use generate_points_m
    implicit none
    logical,intent(in) :: lgeneration

    if(lgeneration) then    ! generate Lagrangian markers on still ellipsoid
       call generate_points
    end if
    call Lag_marker_moving  ! moving markers according to ellipsoid body movement
    call LP_support         ! finding support index for delta function
    call Lag_vel            ! Lagrangian markers velocity

    call particle_Euler_index

  end subroutine ellipsoid

  ! ------------------------------------------------------------ !  
  ! -------------------------------------------------- !
  ! Identifying domain of each lagrangian point to reduce computing time
  subroutine LP_support
    integer  :: n,l,i,j,k

    !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(n,l)
    do n=1,num_p
       do l=1,n_l(n)
          do i=1,nx+1
             if (xets(i) .ge. x_o(l,n)) then
                if(xets(i) - deltax/2.0_wp .ge. x_o(l,n)) then
                   p_iw(l,n)=i-2
                   p_ie(l,n)=i
                else
                   p_iw(l,n)=i-1
                   p_ie(l,n)=i+1
                end if
                goto 200
             endif
          enddo
200       continue
          do j=1,ny+1
             if (yets(j) .ge. y_o(l,n)) then
                if(yets(j) - deltay/2.0_wp .ge. y_o(l,n)) then
                   p_js(l,n)=j-2
                   p_jn(l,n)=j
                else
                   p_js(l,n)=j-1
                   p_jn(l,n)=j+1
                end if
                goto 300
             endif
          enddo
300       continue          
          do k=2,nz-2               ! k is not allowed to reach the boundary
             if (zets(k) .le. z_o(l,n)) then       ! zets decreasing along k
                p_kb(l,n) = k-2
                p_kt(l,n) = k+1
                goto 400
             endif
          enddo
400       continue
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine LP_support

  ! -------------------------------------------------- !

  subroutine Lag_marker_moving
    integer  :: np,n,i,j,l
    real(wp) :: A(3,3), xyzc

    !$OMP PARALLEL DO DEFAULT(SHARED),private(n)
    do n=1,num_p
       A(:,1) = axis_1(:,n)
       A(:,2) = axis_2(:,n)
       A(:,3) = axis_3(:,n)

       np = n_l(n)

       ! X'= A*X        ! new coordinates = rotation matrix * old coordinates
       rx_l(1:np,n)  = A(1,1)*xlp(1:np) + A(1,2)*ylp(1:np) + A(1,3)*zlp(1:np)
       ry_l(1:np,n)  = A(2,1)*xlp(1:np) + A(2,2)*ylp(1:np) + A(2,3)*zlp(1:np)
       rz_l(1:np,n)  = A(3,1)*xlp(1:np) + A(3,2)*ylp(1:np) + A(3,3)*zlp(1:np)

       ! shift particle center to the channel domain

       call s_shift2range(x_c(n),rlenx)
       call s_shift2range(y_c(n),rleny)
       call s_shift2range_z(z_c(n))

       x_o(1:np,n) = rx_l(1:np,n) + x_c(n)
       y_o(1:np,n) = ry_l(1:np,n) + y_c(n)
       z_o(1:np,n) = rz_l(1:np,n) + z_c(n)
       ! shift points x->[0,rlenx],y->[0,rleny]
       call v_shift2range(x_o(1:np,n),rlenx)
       call v_shift2range(y_o(1:np,n),rleny)
       ! z->[-rlenz/2,rlenz/2]
       call v_shift2range_z(z_o(1:np,n))

    end do
    !$OMP END PARALLEL DO

    if(irkk .eq. 3) then
       write(88,'(I7,9f10.6)') itime, ((A(i,j),i=1,3),j=1,3)
       write(1110, 300) itime,(x_c(i),i=1,num_p)
       write(1111, 300) itime,(y_c(i),i=1,num_p)
       write(1112, 300) itime,(z_c(i),i=1,num_p)
       write(1120, 300) itime,(u_c(i),i=1,num_p)
       write(1121, 300) itime,(v_c(i),i=1,num_p)
       write(1122, 300) itime,(w_c(i),i=1,num_p)
       write(120,  300) itime,(om_x(i),i=1,num_p)
       write(121,  300) itime,(om_y(i),i=1,num_p)
       write(122,  300) itime,(om_z(i),i=1,num_p)
    endif
300 format (I7,<num_p>E15.7)

    
  end subroutine Lag_marker_moving


  ! -------------------------------------------------- !
  subroutine Lag_vel
    implicit none
    integer  :: l,n
    real(wp) :: u_rot,v_rot,w_rot

    !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(l)
    do n=1,num_p
       do l=1,n_l(n)

          u_rot=om_y(n)*rz_l(l,n)-om_z(n)*ry_l(l,n)
          v_rot=om_z(n)*rx_l(l,n)-om_x(n)*rz_l(l,n)
          w_rot=om_x(n)*ry_l(l,n)-om_y(n)*rx_l(l,n)

          u_p(l,n)=u_c(n)+u_rot
          v_p(l,n)=v_c(n)+v_rot
          w_p(l,n)=w_c(n)+w_rot
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine Lag_vel

  ! -------------------------------------------------- !
  ! -------------------------------------------------- !
  ! -------------------------------------------------- !
  subroutine finalize_ellip
    ! 
  end subroutine finalize_ellip
  ! ------------------------------------------------------------ !

  ! ------------------------------------------------------------ !
  subroutine particle_Euler_index
    integer :: n,l


    !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(n,l)
    do n=1,num_p
       x_0(n) = p_iw(1,n)
       x_1(n) = p_ie(1,n)
       y_0(n) = p_js(1,n)
       y_1(n) = p_jn(1,n)
       z_0(n) = p_kb(1,n)
       z_1(n) = p_kt(1,n)
       
       do l=2,n_l(n)
          x_0(n) = min(p_iw(l,n),x_0(n))
          x_1(n) = max(p_ie(l,n),x_1(n))
          y_0(n) = min(p_js(l,n),y_0(n))
          y_1(n) = max(p_jn(l,n),y_1(n))
          z_0(n) = min(p_kb(l,n),z_0(n))
          z_1(n) = max(p_kt(l,n),z_1(n))
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine particle_Euler_index
  ! -------------------------------------------------- !

  ! -------------------------------------------------- !
  subroutine v_shift2range(lx,x_range)
    real(wp), intent(inout) :: lx(:)
    real(wp), intent(in)    :: x_range

    lx=mod(lx,x_range) + (sign(0.5_wp,-lx) + 0.5)*x_range
    return
  end subroutine v_shift2range
  ! -------------------------------------------------- !

  ! -------------------------------------------------- !
  subroutine v_shift2range_z(lz)
    real(wp), intent(inout) :: lz(:)

    lz=mod(lz-hrlenz,rlenz) + (sign(0.5_wp,-(lz-hrlenz)) + 0.5)*rlenz - hrlenz
    return
  end subroutine v_shift2range_z
  ! -------------------------------------------------- !

  ! -------------------------------------------------- !
  subroutine s_shift2range(lx,x_range)
    real(wp), intent(inout) :: lx
    real(wp), intent(in)    :: x_range

    lx=mod(lx,x_range) + (sign(0.5_wp,-lx) + 0.5)*x_range
    return
  end subroutine s_shift2range
  ! -------------------------------------------------- !

  ! -------------------------------------------------- !
  subroutine s_shift2range_z(lz)
    real(wp), intent(inout) :: lz

    lz=mod(lz-hrlenz,rlenz) + (sign(0.5_wp,-(lz-hrlenz)) + 0.5)*rlenz - hrlenz
    return
  end subroutine s_shift2range_z
  ! -------------------------------------------------- !


  ! -------------------------------------------------- !
  subroutine single_point_test(f_in,f_out)
    use common_m
    implicit none
    real(wp),intent(in) :: f_in(3,3,6)
    real(wp),intent(out) :: f_out(3,3,6)
    
    real(wp),dimension(54,1) :: vec_ddf,vec_f_in,vec_f_out
    real(wp),dimension(1,54) :: vecT_ddf
    real(wp),dimension(1,1)  :: mat_ddf,matInv_ddf
    real(wp),dimension(54,54) :: mat_proj
    integer :: i,j

    vec_f_in = reshape(f_in,(/54,1/))
    
    vec_ddf = reshape(ddf,(/54,1/))
    vecT_ddf = reshape(vec_ddf,(/1,54/))
    mat_ddf = matmul(vecT_ddf,vec_ddf)
    matInv_ddf = 1.0_wp/mat_ddf
    mat_proj = matmul(matmul(vec_ddf,matInv_ddf),vecT_ddf)

    vec_f_out = matmul(mat_proj,vec_f_in)-2.0_wp*matmul(matmul(vec_ddf,vecT_ddf),vec_f_in)

    f_out = reshape(vec_f_out,(/3,3,6/))
    

  end subroutine single_point_test



  ! -------------------------------------------------- !
  subroutine correlation(nlm)
    use common_m

    implicit none

    integer, intent(in) :: nlm
    integer :: l1,i11,i21,j11,j21,k11,k21
    integer :: l2,i12,i22,j12,j22,k12,k22,i1,i2,j1,j2,k1,k2
    integer :: i_ddf1,i_ddf2,j_ddf1,j_ddf2,k_ddf1,k_ddf2
    integer :: info,lwork
    integer :: n,nl,i,j,k
    integer :: n_eig
    real(wp) :: inv_lambda
    real(wp) :: eig_A(nlm)
    real(wp) :: work(3*nlm-1)
    real(wp) :: corr_A(nlm,nlm)
    real(wp) :: mat_G(nlm,nlm),mat_Gsum(nlm,nlm)
    real(wp) :: vec_V(nlm),vec_Vt(1,nlm)
    real(wp) :: inv_mass_non


    lwork = 3*n_ll-1
    corr_A = 0.0_wp
    mat_Gsum = 0.0_wp
    
    if(ibm_moving .eq. 0) then
       inv_mass_non = 0.0_wp
    else
       inv_mass_non = (rho_f*deltax*deltay*maxval(deltaz))/(vol_ellip * rho_p)
    endif
   
    do n=1,num_p
       nl = n_l(n)
       do l1=1,nl
          i11=p_iw(l1,n)
          i21=p_ie(l1,n)
          j11=p_js(l1,n)
          j21=p_jn(l1,n)
          k11=p_kb(l1,n)
          k21=p_kt(l1,n)
          do l2=l1,nl
             i12=p_iw(l2,n)
             i22=p_ie(l2,n)
             j12=p_js(l2,n)
             j22=p_jn(l2,n)
             k12=p_kb(l2,n)
             k22=p_kt(l2,n)

             i1=max(i11,i12)
             i2=min(i21,i22)
             j1=max(j11,j12)
             j2=min(j21,j22)
             k1=max(k11,k12)
             k2=min(k21,k22)

             do i=i1,i2
                i_ddf1=i-i11+1
                i_ddf2=i-i12+1
                do j=j1,j2
                   j_ddf1=j-j11+1
                   j_ddf2=j-j12+1
                   do k=k1,k2
                      k_ddf1=k-k11+1
                      k_ddf2=k-k12+1
                      corr_A(l1,l2)=corr_A(l1,l2)+ddf(l1,n,i_ddf1,j_ddf1,k_ddf1)*ddf(l2,n,i_ddf2,j_ddf2,k_ddf2)
                   end do
                enddo
             enddo
             corr_A(l2,l1)=corr_A(l1,l2)
          end do
       end do

       corr_A = corr_A + inv_mass_non       ! motion correction 
       
       call DSYEV('V','U',nlm,corr_A,nlm,eig_A, work, lwork, info)
       
       if(info .eq. 3) then

          n_eig = count(eig_A .gt. 1.0e-3_wp,dim=1)

          do i=0,n_eig-1
             vec_V = corr_A(:,nlm-i)
             inv_lambda = 1.0_wp/eig_A(nlm-i)
             call dgemm('N','T',nlm,nlm,1,inv_lambda,vec_V,nlm,vec_V,nlm,0.0_wp,mat_G,nlm)   ! vec_V*vec_V^T / lambda

             mat_Gsum = mat_Gsum + mat_G
             write(*,*) 'eigem decom',i,inv_lambda
             
             write(*,*) 'eigen vector'
             write(*,*) vec_V
             write(*,*) '-----'
             write(*,*) mat_G

          end do
       end if

       if(scale_dv>2.51_wp/maxval(eig_A)) then
          write(*,*) 'The maximum eigenvalue in the correlation matrix is:', maxval(eig_A)
          write(*,*) 'Warning: violation of stability condition, too big weight!'
          write(*,*) scale_dv,'>',2.51_wp/maxval(eig_A)
          write(*,*) '!!!!!!!!!!!!!!!!!!'
       endif
    end do


  end subroutine correlation

  ! -------------------------------------------------- !
  function finv(A) result(Ainv)
    implicit none
    real(wp),intent(in) :: A(:,:)
    real(wp)            :: Ainv(size(A,1),size(A,2))
    real(wp)            :: work(size(A,1)) ! work array for LAPACK
    integer         :: n,info,ipiv(size(A,1)) ! pivot indices
    
    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)
    ! SGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call SGETRF(n,n,Ainv,n,ipiv,info)
    if (info.ne.0) stop 'Matrix is numerically singular!'
    ! SGETRI computes the inverse of a matrix using the LU factorization
    ! computed by SGETRF.
    call SGETRI(n,Ainv,n,ipiv,work,n,info)
    if (info.ne.0) stop 'Matrix inversion failed!'
  end function finv
  ! -------------------------------------------------- !
end module ellipsoid_m
