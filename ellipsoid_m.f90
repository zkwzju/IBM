module ellipsoid_m
  use ellip_common_m
  use common_m, only : irkk,marker_v,lident,ldebug,ibm_moving,nzm,zevn,u_fft,v_fft,w_fft,nx0,ny0,nz,forcing_x,forcing_x1,forcing_x2,forcing_x3,forcing_y,forcing_y1,forcing_y2,forcing_y3,forcing_z,forcing_z1,forcing_z2,forcing_z3,forcing_x01,forcing_x02,forcing_x03,forcing_y01,forcing_y02,forcing_y03,forcing_z01,forcing_z02,forcing_z03,istart,ch_fin
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
    call parser_read('Movement type',ibm_moving,0)
    
    ! a,b,c are along x,y,z, respectively
    call parser_read('Ellipsoid semi-axis a',xa,1.0_wp)
    call parser_read('Ellipsoid semi-axis b',xb,1.0_wp)
    call parser_read('Ellipsoid semi-axis c',xc,1.0_wp)
    call parser_read('Particle mass density',rho_p,2.0_wp)
    call parser_read('Fluid mass density',rho_f,1.0_wp)
    call parser_read('Gravity (y)',g_y,0.0_wp)
    grav=(/0.0_wp,g_y,0.0_wp/) ! gravity
    call parser_read('Method for weight dv',iweight,3)

    call parser_read('Scale factor for particle',scale_p,1.0_wp)
    call parser_read('Scale factor for dv',scale_dv,1.0_wp)
    call parser_read('Rigid body assumption',lrigid,.true.)
    call parser_read('Translational motion',ltranslation,.true.)
    call parser_read('Rotational motion',lrotation,.true.)
    call parser_read('Prestep fraction',frac_step,0.0_wp)
    if(frac_step < 0.0_wp .or. frac_step > 1.0_wp) then
       write(*,*) 'Prestep fraction wrong value:',frac_step
       stop
    endif
    call parser_read('Non-uniform z-grid',lnonuniform,.false.)
    call parser_read('Z-grid stretching parameter',aalpha,0.0_wp)
    call parser_read('Additional checking flag',lcheck,.false.)
    call parser_read('Sphere particle',lsphere,.true.)
    call parser_read('Convection',lconvect,.false.)

    call parser_read('clip x translation',lclip_x,.false.)
    call parser_read('clip y translation',lclip_y,.false.)
    call parser_read('clip z translation',lclip_z,.false.)
    call parser_read('clip x rotation',lclip_ox,.false.)
    call parser_read('clip y rotation',lclip_oy,.false.)
    call parser_read('clip z rotation',lclip_oz,.false.)

    call parser_read('unformatted output',lunformatted,.false.)

    if(lcheck) then
       allocate(   forcing_x1(nx0,ny0,nz),forcing_y1(nx0,ny0,nz),forcing_z1(nx0,ny0,nz) )
       allocate(   forcing_x2(nx0,ny0,nz),forcing_y2(nx0,ny0,nz),forcing_z2(nx0,ny0,nz) )
       allocate(   forcing_x3(nx0,ny0,nz),forcing_y3(nx0,ny0,nz),forcing_z3(nx0,ny0,nz) )
       allocate(   forcing_x01(nx0,ny0,nz),forcing_y01(nx0,ny0,nz),forcing_z01(nx0,ny0,nz) )
       allocate(   forcing_x02(nx0,ny0,nz),forcing_y02(nx0,ny0,nz),forcing_z02(nx0,ny0,nz) )
       allocate(   forcing_x03(nx0,ny0,nz),forcing_y03(nx0,ny0,nz),forcing_z03(nx0,ny0,nz) )
    endif
  end subroutine read_settings

  ! -------------------------------------------------- !
  subroutine init_particle_properties
    use parser_m
    implicit none
    character(len=120)      :: strFile

    logical  :: lrho,lvol,lvolume
    integer  :: i,istatus
    allocate(V_particle(num_p))      ! particles volume
    allocate(rho_particle(num_p))    ! paritcles density

    vol_ellip = 4.0_wp/3.0_wp*pi*xa*xb*xc     ! default ellipsoid volume

    rho_particle = rho_p
    V_particle   = vol_ellip
    
  ! read each particle density if given
    call parser_is_defined('Particles density (filename)',lrho)
    if(lrho) then
       call parser_read('Particles density (filename)',strFile)
       open(unit=170, file=strFile, STATUS='OLD',iostat=istatus)
       if(istatus .ne. 0) then
          write(*,*) 'Problem to open file: ',strFile
          stop
       end if

       do i=1,num_p
          read(170,*,IOSTAT=istatus) rho_p
          if (istatus > 0)  then
             write(*,*) strFile, 'problematic'
             stop
          else if (istatus < 0) then  ! end of file reached
             write(*,*) "particles rotational inertia in file doesn't match particles number",i,num_p
             stop
          else
             rho_particle(i) = rho_p
          end if
       end do
       close(170)
    endif


    ! to read particle volume if given (all the same volume)
    call parser_is_defined('Particle volume',lvol)
    if(lvol) then
       call parser_read('Particle volume',vol_ellip)
    end if

    ! read each particle volume if given in file (different)
    call parser_is_defined('Particles volume (filename)',lvolume)
    if(lvolume) then
       call parser_read('Particles volume (filename)',strFile)
       open(unit=170, file=strFile, STATUS='OLD',iostat=istatus)
       if(istatus .ne. 0) then
          write(*,*) 'Problem to open file: ',strFile
          stop
       end if

       do i=1,num_p
          read(170,*,IOSTAT=istatus) vol_ellip
          if (istatus > 0)  then
             write(*,*) strFile, 'problematic'
             stop
          else if (istatus < 0) then  ! end of file reached
             write(*,*) "particles rotational inertia in file doesn't match particles number",i,num_p
             stop
          else
             V_particle(i) = vol_ellip
          end if
       end do
       close(170)
    endif
  
    
  end subroutine init_particle_properties
  
  ! -------------------------------------------------- !

  subroutine init_location
    use parser_m
    real(wp) :: x0,y0,z0
    character(len=120) :: strFile
    integer  :: n,istatus
    logical  :: llocation
    call parser_read('Particle initial location (x)',x0,rlenx/2.0_wp);
    call parser_read('Particle initial location (y)',y0,rleny/2.0_wp);
    call parser_read('Particle initial location (z)',z0,0.0_wp);

    call parser_is_defined('Particles initial location (filename)',llocation)
    
    if(llocation) then       ! read particles location from file
       call parser_read('Particles initial location (filename)',strFile)

       open(unit=170, file=strFile, STATUS='OLD',iostat=istatus)
       if(istatus .ne. 0) then
          write(*,*) 'Problem to open file: ',strFile
          stop
       end if

       do n=1,num_p
          read(170,*,IOSTAT=istatus)  x0,y0,z0
          if (istatus > 0)  then
             write(*,*) strFile, 'problematic'
             stop
          else if (istatus < 0) then  ! end of file reached
             write(*,*) "Particles location in file doesn't match particles number",n,num_p
             stop
          else
             if(x0 <0 .or. y0<0 .or. z0 > rlenz/2.0_wp .or. z0 < -rlenz/2.0_wp) then
                write(*,*) "Particle is initially outside of the channel!"
                stop
             endif

             x_c(n) = mod(x0,rlenx)
             y_c(n) = mod(y0,rleny)
             z_c(n) = z0
          end if
       end do

       close(170)
    else                   ! set particles location from single input
       do n=1,num_p
          if(x0 <0 .or. y0<0 .or. z0 > rlenz/2.0_wp .or. z0 < -rlenz/2.0_wp) then
             write(*,*) "Particle is initially outside of the channel!"
             stop
          endif
          x_c(n) = mod(x0,rlenx)
          y_c(n) = mod(y0,rleny)
          z_c(n) = z0
       end do
    endif

    
  end subroutine init_location
  ! -------------------------------------------------- !
  subroutine init_orientation
    use parser_m
    implicit none
    character(len=120) :: strFile
    logical  :: linertia,laxis1,laxis2,laxis3,ldirection
    real(wp) :: axis1(3),axis2(3),axis3(3),mass
    integer  :: i,istatus
    
    ! particle orientation vetors
    allocate(axis_1(3,num_p),axis_2(3,num_p),axis_3(3,num_p))
    allocate(I_particle(3,num_p))

    mass = vol_ellip * rho_p      ! default for identical ellipsoid
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
    axis3 = cross(axis1,axis2)

    ! All particles have same mass and initial orientation, if not specified
    do i=1,num_p      
       axis_1(:,i) = axis1
       axis_2(:,i) = axis2
       axis_3(:,i) = axis3
       I_particle(:,i) = I_ellip
    end do

    ! read initial orientation if defined
    call parser_is_defined('Particles 1 & 2 axes direction (filename)',ldirection)
    if(ldirection) then
       call parser_read('Particles 1 & 2 axes direction (filename)',strFile)
       open(unit=170, file=strFile, STATUS='OLD',iostat=istatus)
       if(istatus .ne. 0) then
          write(*,*) 'Problem to open file: ',strFile
          stop
       end if

       do i=1,num_p
          read(170,*,IOSTAT=istatus) axis1,axis2
          ! normalize
          axis1 = axis1/sqrt(axis1(1)**2+axis1(2)**2+axis1(3)**2)
          axis2 = axis2/sqrt(axis2(1)**2+axis2(2)**2+axis2(3)**2)
          if (istatus > 0)  then
             write(*,*) strFile, 'problematic'
             stop
          else if (istatus < 0) then  ! end of file reached
             write(*,*) "particles direction in file doesn't match particles number",i,num_p
             stop
          else
             axis3 = cross(axis1,axis2)
             axis_1(:,i) = axis1
             axis_2(:,i) = axis2
             axis_3(:,i) = axis3
          end if
       end do
       close(170)
    endif

    ! read particles inertia if defined
    call parser_is_defined('Particles rotational inertia (filename)',linertia)
    if(linertia) then
       call parser_read('Particles rotational inertia (filename)',strFile)
       open(unit=170, file=strFile, STATUS='OLD',iostat=istatus)
       if(istatus .ne. 0) then
          write(*,*) 'Problem to open file: ',strFile
          stop
       end if

       do i=1,num_p
          read(170,*,IOSTAT=istatus) I_ellip
          if (istatus > 0)  then
             write(*,*) strFile, 'problematic'
             stop
          else if (istatus < 0) then  ! end of file reached
             write(*,*) "particles rotational inertia in file doesn't match particles number",i,num_p
             stop
          else
             I_particle(:,i) = I_ellip
          end if
       end do
       close(170)
    endif

    
  end subroutine init_orientation
  ! -------------------------------------------------- !

  subroutine init_velocity
    use parser_m
    implicit none
    character(len=120), allocatable :: buffer(:)
    logical :: lvelocity
    integer :: i, nvelocity
        ! setting marker moving velocity
    call parser_is_defined('Particles initial translational velocity (x,y,z)',lvelocity)
    if(lvelocity) then
       call parser_getsize('Particles initial translational velocity (x,y,z)',nvelocity) 
       if(nvelocity/3 .ne. num_p) then
          write(*,*) "Problem: number of values for Particles initial velocity!"
          stop
       endif
       allocate(buffer(3*num_p))
       call parser_read('Particles initial translational velocity (x,y,z)',buffer)
       do i=1,num_p
          read(buffer(i*3-2),*) u_c(i)
          read(buffer(i*3-1),*) v_c(i)
          read(buffer(i*3  ),*) w_c(i)
       enddo
       deallocate(buffer)
    else
       u_c = 0.0_wp                ! zero initial velocity
       v_c = 0.0_wp
       w_c = 0.0_wp
    endif
    
    ! read initial rotational speed
    call parser_is_defined('Particles initial angular velocity (x,y,z)',lvelocity)
    if(lvelocity) then
       call parser_getsize('Particles initial angular velocity (x,y,z)',nvelocity)   ! nmarker re-used
       if(nvelocity/3 .ne. num_p) then
          write(*,*) "Problem: number of values for Particles initial angular velocity !"
          stop
       endif
       allocate(buffer(3*num_p))
       call parser_read('Particles initial angular velocity (x,y,z)',buffer)
       do i=1,num_p
          read(buffer(i*3-2),*) om_x(i)
          read(buffer(i*3-1),*) om_y(i)
          read(buffer(i*3  ),*) om_z(i)
       enddo
       deallocate(buffer)
    else
       om_x = 0.0_wp               ! zero angular velocity
       om_y = 0.0_wp
       om_z = 0.0_wp
    endif

    ! for forced motion only 
    if(ibm_moving .eq. 2) then
       allocate(n_wave(num_p))
       allocate(phase0(num_p))
       n_wave = 0.0_wp
       phase0 = 0.0_wp
       call parser_is_defined('Particles oscilation wave number and initial phase',lvelocity)
       if(lvelocity) then
          allocate(buffer(2*num_p))
          call parser_read('Particles oscilation wave number and initial phase',buffer)
          do i=1,num_p
             read(buffer(i*2-1),*) n_wave(i)
             read(buffer(i*2)  ,*) phase0(i)
          enddo
          deallocate(buffer)
       endif
    endif

       
  end subroutine init_velocity
  
  ! -------------------------------------------------- !
  subroutine init_markers
    use parser_m
    use generate_points_m
    implicit none
    logical  :: lflag
    character(len=120), allocatable :: buffer(:)
    integer  :: i
    
    ! particle Lagrangian markers coordinates wrt particle center
    if(lident) then       ! identical particles
       allocate(xlp(n_ll,1),ylp(n_ll,1),zlp(n_ll,1))
    else
       allocate(xlp(n_ll,num_p),ylp(n_ll,num_p),zlp(n_ll,num_p))
    endif
    
    call parser_read('Lagrangian markers generation method',imethod,2)
    ! approximate surface area of an ellipsoid

    select case(imethod)
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
       write(*,*) 'Method can only be 2 or 3, for length specified, and reading from file, respectively'
       stop
    end select

    call parser_read('Lagrangian weights file',sname_w,'lagrangian_weights.dat')

    ! read checking points
    call parser_is_defined('Markers to check fluid velocity',lcheckmarker) 
    if(lcheckmarker) then
       open(unit=355, file='fluid_vel_at_markers')
       call parser_getsize('Markers to check fluid velocity',nmarker)
       nmarker = nmarker/2
       allocate(u_interp(nmarker),v_interp(nmarker),w_interp(nmarker))
       allocate(u_fft(nx0,ny0,nz),v_fft(nx0,ny0,nz),w_fft(nx0,ny0,nz))
       
       allocate(buffer(nmarker*2))
       allocate(marker_ind(nmarker,2))
       call parser_read('Markers to check fluid velocity',buffer)
       do i=1,nmarker
          read(buffer(2*(i-1)+1),*) marker_ind(i,1)
          read(buffer(2*(i-1)+2),*) marker_ind(i,2)
       enddo
       deallocate(buffer)
    end if

    call generate_points
  end subroutine init_markers

  ! -------------------------------------------------- !
  
  subroutine initialize_ellip
    use parser_m
    implicit none

    allocate(x_0(num_p),x_1(num_p),y_0(num_p),y_1(num_p),z_0(num_p),z_1(num_p))
    allocate(for_px(num_p),for_py(num_p),for_pz(num_p))
    allocate(torq_x(num_p),torq_y(num_p),torq_z(num_p))
    for_px = 0.0_wp
    for_py = 0.0_wp
    for_pz = 0.0_wp
    torq_x = 0.0_wp
    torq_y = 0.0_wp
    torq_z = 0.0_wp


    call init_particle_properties
    call init_location
    call init_orientation
    call init_velocity
    call init_markers

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

    if(ldebug) write(*,*) 'start: LP_suport'

    !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(n,l,i,j,k)
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
          if(z_o(l,n) > zets(1) .or. z_o(l,n) < zets(nz)) then
             write(*,*) "z_o out of bound",z_o(l,n)
             stop
          endif
          if(aalpha .eq. -1) then
             if(z_o(l,n) > 0) then
                if(z_o(l,n) > zevn(1)-0.5_wp*deltaz(1)) then
                   p_kb(l,n) = 1
                   p_kt(l,n) = 2
                   marker_v(l,n) = deltax*deltay*0.5_wp*deltaz(1) !deltaz uniform
                else  ! interial points
                   do k=2,nzm/2
                      if(z_o(l,n) .ge. zevn(k)) then
                         if(z_o(l,n) .ge. zevn(k)+0.5_wp*deltaz(1)) then
                            p_kb(l,n) = k
                            p_kt(l,n) = k+2
                         else
                            p_kb(l,n) = k-1
                            p_kt(l,n) = k+1
                         endif
                         marker_v(l,n) = deltax*deltay*deltaz(1)
                         goto 400
                      endif
                   enddo
                endif
             else     ! z_o(l,n) <=0
                if(z_o(l,n) < zevn(nz)+0.5_wp*deltaz(1)) then
                   p_kb(l,n) = nz-1
                   p_kt(l,n) = nz
                   marker_v(l,n) = deltax*deltay*0.5_wp*deltaz(1) !deltaz uniform
                else
                   do k=nzm/2,nzm
                      if(z_o(l,n) .ge. zevn(k)) then
                         if(z_o(l,n) .ge. zevn(k)+0.5_wp*deltaz(1)) then
                            p_kb(l,n) = k
                            p_kt(l,n) = k+2
                         else
                            p_kb(l,n) = k-1
                            p_kt(l,n) = k+1
                         endif
                         marker_v(l,n) = deltax*deltay*deltaz(1)
                         goto 400
                      endif
                   enddo
                end if
             endif
          else  ! aalpha neq -1
             if(z_o(l,n) < zets(1) .and. z_o(l,n) > zets(2)) then    ! zets decreasing along k
                p_kb(l,n) = 1
                p_kt(l,n) = 2
                marker_v(l,n) = deltax*deltay*dabs(zets(1)-zets(2))/2.0_wp
             elseif(z_o(l,n) > zets(nz) .and. z_o(l,n) < zets(nz-1)) then
                p_kb(l,n) = nz-1
                p_kt(l,n) = nz
                marker_v(l,n) = deltax*deltay*dabs(zets(nz)-zets(nz-1))/2.0_wp
             else
                if(z_o(l,n) .ge. 0.0_wp) then
                   do k=3,floor((nz+1)/2.0_wp)
                      if(zets(k) .le. z_o(l,n)) then
                         p_kb(l,n) = k-2
                         p_kt(l,n) = k+1
                         marker_v(l,n) = deltax*deltay*dabs(zets(k+1)-zets(k-1))/2.0_wp
                         goto 400
                      endif
                   enddo
                else                               ! anti-symmetric treatment 
                   do k = nz-2, floor((nz+1)/2.0_wp),-1
                      if(zets(k) .ge. z_o(l,n)) then
                         p_kb(l,n) = k-1
                         p_kt(l,n) = k+2
                         marker_v(l,n) = deltax*deltay*dabs(zets(k+1)-zets(k-1))/2.0_wp
                         goto 400
                      endif
                   enddo
                endif
             endif
          endif       ! endif aalpha
400       continue
       enddo
    enddo
    !$OMP END PARALLEL DO

    if(ldebug) write(*,*) 'end: LP_suport'
  end subroutine LP_support

  ! -------------------------------------------------- !

  subroutine Lag_marker_moving
    integer  :: np,n,i,j,l
    real(wp) :: A(3,3), xyzc

    !$OMP PARALLEL DO DEFAULT(SHARED),private(n,A,np)
    do n=1,num_p
       A(:,1) = axis_1(:,n)
       A(:,2) = axis_2(:,n)
       A(:,3) = axis_3(:,n)

       np = n_l(n)

       ! X'= A*X        ! new coordinates = rotation matrix * original coordinates
       if(lident) then
          rx_l(1:np,n)  = A(1,1)*xlp(1:np,1) + A(1,2)*ylp(1:np,1) + A(1,3)*zlp(1:np,1)
          ry_l(1:np,n)  = A(2,1)*xlp(1:np,1) + A(2,2)*ylp(1:np,1) + A(2,3)*zlp(1:np,1)
          rz_l(1:np,n)  = A(3,1)*xlp(1:np,1) + A(3,2)*ylp(1:np,1) + A(3,3)*zlp(1:np,1)
       else
          rx_l(1:np,n)  = A(1,1)*xlp(1:np,n) + A(1,2)*ylp(1:np,n) + A(1,3)*zlp(1:np,n)
          ry_l(1:np,n)  = A(2,1)*xlp(1:np,n) + A(2,2)*ylp(1:np,n) + A(2,3)*zlp(1:np,n)
          rz_l(1:np,n)  = A(3,1)*xlp(1:np,n) + A(3,2)*ylp(1:np,n) + A(3,3)*zlp(1:np,n)
       endif

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
       if((ch_fin.eq.1) .and. (itime.eq.istart+1)) goto 5000
       
       if(lunformatted) then
          write(88  ),itime,axis_1,axis_2
          write(1110),itime,x_c
          write(1111),itime,y_c
          write(1112),itime,z_c
          write(1120),itime,u_c
          write(1121),itime,v_c
          write(1122),itime,w_c
          write(120 ),itime,om_x
          write(121 ),itime,om_y
          write(122 ),itime,om_z
       else
          write(88,   312) itime,axis_1,axis_2
          write(1110, 311) itime,x_c
          write(1111, 311) itime,y_c
          write(1112, 311) itime,z_c
          write(1120, 311) itime,u_c
          write(1121, 311) itime,v_c
          write(1122, 311) itime,w_c
          write(120,  311) itime,om_x
          write(121,  311) itime,om_y
          write(122,  311) itime,om_z
       endif
5000   continue ! skip first step output when restarting
    endif
311 format (I7,<num_p>E15.7)
312 format (I7,<num_p*6>f10.6)   

    
  end subroutine Lag_marker_moving


  ! -------------------------------------------------- !
  subroutine Lag_vel
    implicit none
    integer  :: l,n
    real(wp) :: u_rot,v_rot,w_rot

    !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(n,l,u_rot,v_rot,w_rot)
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

    lx = mod(lx,x_range*(1.0_wp+1.0e-15_wp))
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
    if(lx .ne. x_range) then
       lx = mod(lx,x_range)
    endif
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
    use, intrinsic :: iso_c_binding, only: c_size_t
    use IFPORT
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
    integer (C_SIZE_T) array_len, array_size

    lwork = 3*n_ll-1
    corr_A = 0.0_wp
    mat_Gsum = 0.0_wp
    
    if(ibm_moving .eq. 0) then
       inv_mass_non = 0.0_wp
    else
       inv_mass_non = (rho_f*deltax*deltay*maxval(abs(deltaz)))/(vol_ellip * rho_p)
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

       if(ibm_moving .eq. 1) then
          do i=1,nlm
             corr_A(i,i) = corr_A(i,i) + inv_mass_non       ! motion correction
          enddo
       endif
       
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


       ! sort the eigenvalues
       array_len = nlm
       array_size = 8
       call qsort(eig_A,array_len, array_size,com_function)
       write(*,*) '=============================='
       write(*,*) 'The maximum/minimum eigenvalue in the correlation matrix is:', eig_A(1),'/',eig_A(nlm)
       if(n_ll<5) then          ! write largest 5 eigenvalues
          write(*,*) 'Correlation matrix largest eigenvalues:',eig_A(1:n_ll)
       else
          write(*,*) 'Correlation matrix largest eigenvalues:',eig_A(1:5)
       endif
 
       if(scale_dv>2.51_wp/eig_A(1)) then
          write(*,*) 'Warning: violation of stability condition, too big weight!'
          write(*,*) scale_dv,'>',2.51_wp/eig_A(1)
          write(*,*) '!!!!!!!!!!!!!!!!!!'
       endif
       write(*,*) '=============================='
    end do


  end subroutine correlation
  ! -------------------------------------------------- !
  integer(2) function com_function(par1,par2)
    real(wp) :: par1, par2
    com_function = par1-par2
  end function com_function
  
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

  FUNCTION cross(a, b)
    real(wp), DIMENSION(3) :: cross
    real(wp), DIMENSION(3), INTENT(IN) :: a, b

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  END FUNCTION cross
  ! -------------------------------------------------- !
end module ellipsoid_m
